# Phase 3: Decoupling the Downstream Electrostatics Packages — Design

Date: 2026-09-15
Status: in progress — §5.1 ParticleMeshEwald complete and reviewed; §5.4 EwaldSummations descoped by the user

## 1. Context

`~/project/q2dmd` holds ExTinyMD plus five electrostatics packages that plug into it.
Phases 1 and 2 gave ExTinyMD its own electrostatics standard library — `Ewald3D`, `Ewald2D`,
`ICMEwald2D`, `ICMEwald3D`, `PME3D`, `ICMPME3D` — built on a two-layer contract: a
framework-free core taking a plan object plus array-of-structs positions and charges, and a
thin MD adapter implementing `ExTinyMD.energy` / `ExTinyMD.update_acceleration!`.

Phase 3 applies that same contract to the five downstream packages. Two requirements, both
from the user:

1. **Make them work with current ExTinyMD.** All five pin `ExTinyMD = "0.2"`; ExTinyMD is at
   0.3 and moved to CellListMap 0.10, so none of them resolve.
2. **Decouple them.** "One should never need to init an ExTinyMD type thing when using the
   standalone force packages" — they should take positions as AoS, charges, and their own
   parameters, nothing more.

### Measured state

| package | src LOC | test LOC | dominant ExTinyMD coupling | exported names |
|---|---|---|---|---|
| QuasiEwald | 1789 | 510 | **79× `Point`**, `SimulationInfo`, `MDSys` | 50 |
| FastSpecSoG | 1465 | 478 | `CellList3D`, `SimulationInfo`, `Atom` | 48 |
| EwaldSummations | 1164 | 320 | **100× `Point`**, `CellList3D`, `position_check3D` | 43 |
| SoEwald2D | 823 | 332 | `SimulationInfo`, `MDSys` | 25 |
| ParticleMeshEwald | ~200 | 104 | **none in `src/`** | 4 |

Three facts that shape the work:

- **`Point` is the coupling.** Roughly 180 uses across QuasiEwald and EwaldSummations. Since
  `SVector{3,T}` supports everything `Point` does — arithmetic, iteration, indexing — this is
  a large but mechanical substitution. The exception is ExTinyMD's `dist2`, which becomes
  `sum(abs2, a - b)`.
- **No package declares `[extensions]`.** ParticleMeshEwald's `[weakdeps]` on ExTinyMD and
  EwaldSummations have therefore been inert since they were written: without an
  `[extensions]` entry and an `ext/` directory, Julia never loads anything for them. So PME
  is already free of ExTinyMD in `src/` but gains nothing from it either.
- **Two packages cannot drive an MD run at all.** EwaldSummations and FastSpecSoG define
  `energy` and a free-function `force`, but never `ExTinyMD.update_acceleration!` — the method
  `simulate!` actually calls.

## 2. Goals

1. All five packages resolve against current ExTinyMD and their test suites pass.
2. Each exposes a framework-free core: plan object + AoS positions + charges, no ExTinyMD type
   constructed.
3. ExTinyMD moves to `[weakdeps]` in all five, with the MD glue in `ext/`.
4. EwaldSummations and FastSpecSoG gain the `update_acceleration!` they have never had.

## 3. Non-goals

- Rewriting any numerical algorithm. This phase is restructuring; the physics stays as it is.
- Adding features. No new methods, no GPU, no threading that does not already exist.
- Performance work beyond what decoupling implies.
- Deciding whether any package should be retired. See §8.

## 4. The contract, restated for these packages

Identical in shape to what Phases 1–2 established inside ExTinyMD.

### 4.1 Core layer — no ExTinyMD

```julia
plan = QuasiEwald.QuasiEwaldPlan(n_atoms, L; α, s, γ, ...)

E = QuasiEwald.energy(plan, poses, charges)
F = QuasiEwald.force(plan, poses, charges)
QuasiEwald.force!(F, plan, poses, charges)
```

- `poses` is AoS. `Vector{SVector{3,T}}` canonical; kernels index `p[1]`, `p[2]`, `p[3]` only,
  so `NTuple{3,T}` and ExTinyMD's `Point` also work without a conversion layer.
- The plan owns all scratch — cell lists, k-sets, quadrature nodes, FFT plans, force buffers —
  so queries allocate nothing after construction.
- An optional `neighbor_list =` keyword lets a caller pass a list it already maintains.

### 4.2 Naming: `energy`/`force` are defined but NOT exported

Each package defines `energy`, `force` and `force!` in its own namespace and **does not export
them**. Callers write `QuasiEwald.energy(plan, poses, charges)`.

The reason: if five packages each exported `energy`, then `using EwaldSummations, QuasiEwald`
would make the name ambiguous and error on first use. The alternatives are worse — a shared
owner for the generic would reintroduce a hard common dependency, which is the coupling this
phase exists to remove, and package-prefixed names like `quasi_ewald_energy` are noise. Plan
**constructors** are exported, since those are what a user needs to discover.

Note ExTinyMD's own stdlib uses `coulomb_energy`/`coulomb_force` instead, because
`ExTinyMD.energy` already exists with a four-argument MD signature. The downstream packages
have no such collision and use the plain names.

### 4.3 Adapter layer — in `ext/`

```julia
# ext/QuasiEwaldExTinyMDExt.jl
ExTinyMD.energy(interaction, finder, sys, info)
ExTinyMD.update_acceleration!(interaction, finder, sys, info)
```

Each adapter extracts charges and positions, honouring ExTinyMD's id/slot indirection —
`sys.atoms` is indexed by particle **id**, `info.particle_info` by storage **slot** — then
calls the core and, for forces, divides by mass and accumulates. Phase 1's
`src/interactions/electrostatics/adapter.jl` is the reference implementation; the gather
pattern and its `_finder_list` fallback for `NoNeighborFinder` transfer directly.

### 4.3a The interaction type must live in the extension, not `src/`

**Discovered while decoupling ParticleMeshEwald, and it governs the remaining four packages.**

`MDSys` requires its interactions to subtype `ExTinyMD.AbstractInteraction`. A supertype is
fixed at struct definition and **no extension can retrofit one**. So a type defined in `src/`,
where ExTinyMD does not exist, can never be placed in `sys.interactions` — it is not a matter
of writing the adapter differently.

ParticleMeshEwald absorbed this without loss: it has no forces, so it could never drive
`simulate!` regardless, and its extension supplies `ExTinyMD.energy` alone. QuasiEwald,
SoEwald2D and FastSpecSoG all currently *do* put their interactions in `sys.interactions`, so
for them the naive reading of this design would trade away MD integration entirely — which is
half of what the user asked for.

The resolution follows the two-layer contract already in force, taken one step further:

- **`src/` defines the plan** — parameters plus scratch, no ExTinyMD, constructed and queried
  from plain arrays. This is what a standalone user touches.
- **`ext/` defines a thin interaction wrapper** subtyping `ExTinyMD.AbstractInteraction` and
  holding a plan, plus the `ExTinyMD.energy` / `ExTinyMD.update_acceleration!` methods on it.
  This is what goes in `sys.interactions`.

```julia
# src/  — no ExTinyMD anywhere
struct QuasiEwaldPlan{T, ...}
    ...
end
QuasiEwald.energy(plan, poses, charges)

# ext/QuasiEwaldExTinyMDExt.jl
struct QuasiEwaldInteraction{P} <: ExTinyMD.AbstractInteraction
    plan::P
end
ExTinyMD.energy(i::QuasiEwaldInteraction, finder, sys, info) = ...
ExTinyMD.update_acceleration!(i::QuasiEwaldInteraction, finder, sys, info) = ...
```

The wrapper is a few lines and carries no physics. Standalone users never see it; MD users
construct it from a plan. Both requirements are met without ExTinyMD becoming a hard
dependency.

A consequence worth stating: the interaction *type name* changes for these packages, since the
old name is the plan now. That is a breaking change, appropriate at 0.x, and the migration is
one line at each construction site.

#### 4.3a-bis Preserving the old name: the dispatcher-function pattern

**Discovered while decoupling QuasiEwald, which is the first package to actually need a wrapper
*type* rather than just a method.** It supersedes the closing paragraph above: the interaction
type name does *not* have to change, though preserving it costs something.

§4.3a is right that a supertype is fixed at struct definition. What it does not say is the
stronger constraint underneath: **a struct definition cannot be dot-qualified at all.**
`struct QuasiEwald.Foo <: ExTinyMD.AbstractInteraction ... end` is not legal Julia, so an
extension cannot define a type *into* its parent package's namespace the way it defines a
method into a parent's generic function. The wrapper type genuinely lives in the extension
module's own namespace.

To keep `using QuasiEwald, ExTinyMD; QuasiEwaldShortInteraction(...)` working unchanged, `src/`
declares the name as a plain **dispatcher function** that forwards through
`Base.get_extension`:

```julia
# src/ — no ExTinyMD, and no type of this name either
for name in (:QuasiEwaldShortInteraction, :QuasiEwaldLongInteraction, :SortingFinder)
    @eval function $name(args...; kwargs...)
        ext = Base.get_extension(QuasiEwald, :QuasiEwaldExTinyMDExt)
        ext === nothing && error($(string(name)) * " requires ExTinyMD to be loaded")
        return getfield(ext, $(QuoteNode(name)))(args...; kwargs...)
    end
end
```

**The cost, and it must be checked per package before choosing this:** the preserved name is a
*function*, not a type. Every construction site keeps working; every **type-position** use
breaks — `::QuasiEwaldShortInteraction`, `Vector{QuasiEwaldShortInteraction}`,
`isa QuasiEwaldShortInteraction`, and any method dispatching on it. Before adopting the pattern
for a package, grep the package *and its dependents* for type-position uses of the name. For
QuasiEwald this was zero (construction only, 20 sites), so the pattern was free and the
migration was one convenience constructor instead of edits to five test files and two examples.
If a package has type-position uses, rename instead, per §4.3a's original advice, and let the
plan take the old name.

A second, unrelated Julia constraint surfaced at the same time: a package cannot define a bare
`function energy(...)` (per §4.2) while its module still does a blanket `using ExTinyMD`, since
ExTinyMD exports `energy` and Julia refuses to shadow a `using`-imported binding with a new
local definition — even for a disjoint signature. **Switch to `import ExTinyMD` plus a narrow
`using ExTinyMD: <the names actually used unqualified>`.** This bites during the transition,
while decoupled and not-yet-decoupled code coexist in one module, so expect it in SoEwald2D and
FastSpecSoG too.

#### 4.3b Version numbers: a decoupling is a breaking change

Moving `ExTinyMD` to `[weakdeps]` removes exported names (the old `Pkg_Es`/`Pkg_Fs!`-style
adapter entry points become `ExTinyMD.energy`/`update_acceleration!` methods), raises the julia
floor, and changes dependency bounds. For any package **registered in General**, that requires a
minor bump under 0.x semver. Check registration with
`Pkg.Registry.reachable_registries()` rather than assuming — of the five, QuasiEwald, SoEwald2D,
FastSpecSoG and ExTinyMD are registered; ParticleMeshEwald is not.

Applied: QuasiEwald 0.2.1 → **0.3.0**.

### 4.4 `Project.toml` shape

```toml
[deps]
# ExTinyMD removed

[weakdeps]
ExTinyMD = "fec76197-d59f-46dd-a0ed-76a83c21f7aa"

[extensions]
<Pkg>ExTinyMDExt = "ExTinyMD"

[compat]
ExTinyMD = "0.3"
CellListMap = "0.10"     # where the package uses it

[extras]
ExTinyMD = "fec76197-d59f-46dd-a0ed-76a83c21f7aa"

[targets]
test = ["Test", "ExTinyMD", ...]
```

ExTinyMD in `[extras]`/`test` is required so the extension is exercised by `Pkg.test()` —
the same arrangement Phase 2 used for FINUFFT. It does not make ExTinyMD a runtime
dependency.

## 5. Per-package work

Ordered by increasing coupling, so the pattern is validated cheaply before it is applied to
the packages where a mistake is expensive.

### 5.1 ParticleMeshEwald — validate the pattern

Already free of ExTinyMD in `src/`. Work: convert the SoA query
(`energy(pme, x, y, z, q)`) to AoS, write the `ext/` that its `[weakdeps]` has always
implied, add `[extensions]`, bump compat. Also fix two defects Phase 1's survey found:
`energy_long` scales the caller's coordinate arrays in place and divides them back, and
`energy_short_kernel!` indexes a thread-local accumulator by `Threads.threadid() - 1`, which
is element 0 when `threadid()` is 1.

**This package is now largely superseded by ExTinyMD's `PME3D`** — see §8. This phase makes it
work as asked and does not act on that.

**Outcome, as completed.** CellListMap bumped to 0.10 (the co-resolution blocker), AoS query
API with no caller mutation, `energy`/`energy_short`/`energy_long` un-exported, the
`threadid()-1` accumulator and the whole KernelAbstractions kernel removed (no GPU backend was
ever used anywhere in the repository, and the kernel was the bug's only source), the
`examples/utils.jl` include removed after confirming by grep that what it defined was
referenced nowhere, and an `ext/` supplying `ExTinyMD.energy` only. Baseline energies
unchanged across eight configurations. 16 tests, up from 9.

Two things it exposed that the plan did not anticipate. First, capturing baseline energies is
harder than it looks: ExTinyMD's `SimulationInfo` consumes `rand()` internally, so a
comparison script that constructs one shifts the RNG stream and the "before" numbers are not
comparable to the "after" ones. The implementer caught this itself and re-derived the baseline
from unmodified source via `git stash`. Any per-package before/after check in this phase must
either avoid `SimulationInfo` or seed immediately before each measurement. Second, the first
failure after the compat bump was not the expected `MethodError` but an unsatisfiable
`Pkg.test()` resolve, because the registered EwaldSummations and ExTinyMD both still pin
CellListMap 0.9 — so part of the `[weakdeps]` cleanup had to be pulled forward.

### 5.2 SoEwald2D — lightest of the four coupled packages

823 src LOC, coupling concentrated in `SimulationInfo`/`MDSys` rather than `Point`. Already
defines `ExTinyMD.update_acceleration!`, so the adapter mostly moves rather than being written.

### 5.3 FastSpecSoG

1465 src LOC. Coupling via `CellList3D`, `SimulationInfo`, `Atom`. Uses the **older** interface
generation — `energy(interaction, neighbor, info, atoms)` — and has **no**
`update_acceleration!`, so its adapter is new work rather than a move.

### 5.4 EwaldSummations — OUT OF SCOPE

**The user removed this package from the phase on 2026-09-15:** "no need to include
EwaldSummations, the others are good."

It therefore stays on `ExTinyMD = "0.2"` / `CellListMap = "0.9"`, is not decoupled, is not
thinned, and its k-space implementations are not deleted. Two consequences:

- **The migration check is no longer a gate.** Phase 1's spec §8.2 deferred to this phase a
  script proving ExTinyMD's stdlib reproduces EwaldSummations before the latter's code was
  removed. Nothing is being removed, so nothing needs gating. The check retains independent
  value as cross-validation against an implementation nobody has touched, but it is optional
  and unscheduled.
- **FastSpecSoG cannot keep EwaldSummations as a test dependency.** Its tests use
  `Ewald2DInteraction`, `Ewald2D_short_energy_N` and `Ewald2D_long_energy_N` as accuracy
  references. Since EwaldSummations stays on CellListMap 0.9, FastSpecSoG's test target will
  not resolve once it moves to 0.10. The substitution is ExTinyMD's own `Ewald2D`, which is a
  strictly better reference: same physics, and 1128 tests behind it rather than an
  unmaintained package. See §5.3.

*Original plan, retained for the record:*

1164 src LOC, **100 `Point` uses**. The user chose to thin this package: ExTinyMD's stdlib now
owns the k-space Ewald methods, and EwaldSummations keeps what is genuinely its own — direct
lattice summation, `ICM_reflect`, and the error estimates from arXiv:2503.18126. Its Ewald3D,
Ewald2D and ICM+Ewald implementations are superseded.

Before deleting any of it, run the migration check Phase 1's spec §8.2 deferred to this phase:
a script that compares ExTinyMD's stdlib against EwaldSummations' implementations on shared
configurations, proving the replacement reproduces the original before the original goes away.
That script is the gate on the deletion, not an afterthought.

Also drop ForwardDiff — it exists only for the per-pair `ForwardDiff.derivative` in the force
kernels that Phase 1 replaced with an analytic expression — and drop the
`@reexport using ExTinyMD`, which is the single largest source of the coupling.

### 5.5 QuasiEwald — heaviest

1789 src LOC, **79 `Point` uses**, plus `SimulationInfo` and `MDSys` in the adapter path.
Already on the current interface generation. Its `SortingFinder` is a custom neighbour finder
subtyping `ExTinyMD.AbstractNeighborFinder`, which cannot live in the core once ExTinyMD is a
weak dependency — it moves to the extension, or the core grows its own z-sorting helper and
the extension wraps it. Decide when the code is in front of you; note it here so it is not a
surprise.

## 6. Testing

Each package keeps its existing suite, which is the regression net for the restructuring: if
the numbers move, the decoupling broke something. Additionally, per package:

- **A standalone smoke test** that constructs a plan and evaluates energy and force from plain
  arrays, with `ExTinyMD` **not loaded**. This is the requirement the whole phase exists for,
  and it is the one test that would fail if the coupling crept back.
- **An adapter test** driving the interaction through `simulate!`, which is what catches the
  id/slot, mass-division and accumulation faults that per-call tests cannot.
- **A before/after numerical check** on at least one configuration per package, comparing the
  decoupled implementation against the pre-Phase-3 result. Record the reference values in the
  test so a future change has something to violate.

Phase 1 and 2 saw ten defective test specifications, every one a test that passed for the
wrong reason rather than one that failed. The smoke test above is the one most at risk of
being vacuous: it must genuinely run in a session where ExTinyMD was never loaded, not merely
avoid mentioning it.

## 6a. Lessons from ParticleMeshEwald, to be applied to the other four

Phase 3a was deliberately sequenced first to find these cheaply. All were confirmed by its
reviewer independently rather than taken on report.

**1. Task 1 will fail to resolve before it fails to compile.** The plan predicted a
`MethodError` from CellListMap's renamed keyword as the first failure after the compat bump.
The actual first failure was an unsatisfiable `Pkg.test()` resolve: the *published* versions
of ExTinyMD and EwaldSummations still pin CellListMap 0.9, so any package listing them as
test dependencies cannot resolve once it moves to 0.10. None of the four remaining packages'
published versions are on 0.10 either, so **each plan must schedule dropping or re-pinning
stale registry test dependencies as part of the compat task**, not discover it mid-flight.

**2. The wrapper pattern of §4.3a is validated, but PME did not exercise it.** PME has no
forces, so its extension supplies `ExTinyMD.energy` only and never needed a wrapper at all.
The four remaining packages do need one. Their adapter tests must therefore **drive
`simulate!` with the wrapper placed in `sys.interactions`**, not merely call
`ExTinyMD.energy` directly — that is the only thing that tests the wrapper under load, and
it is what catches the id/slot, mass-division and accumulation faults that per-call tests
cannot.

**3. Baseline captures are contaminated by `SimulationInfo`.** It consumes `rand()`
internally, so a before/after comparison script that constructs one shifts the RNG stream and
the two sets of numbers are not comparable. Seed immediately before each measurement, or
avoid `SimulationInfo` in the comparison entirely. PME's implementer caught its own
contaminated capture and re-derived from unmodified source via `git stash`.

**4. `@inbounds arr[0] += x` does not throw.** The `threadid() - 1` accumulator bug was worse
than "would raise a `BoundsError`" — under `@inbounds` it is an undefined-behaviour write that
produces a wrong energy silently. Worth knowing wherever this idiom appears in the remaining
packages; grep for `threadid` in each.

**5. Verify every "this test would fail" claim by making it fail.** PME's implementer
hand-verified the standalone test by adding `using ExTinyMD` to its subprocess script, and the
threading bug by observing the silent corruption. Across Phases 1–2, eleven specifications
proved defective and every one was a test that passed for the wrong reason. This is the
single discipline that has caught the most.

## 6b. Release checklist — `[sources]` must be removed across all five together

ExTinyMD 0.3 is **not in the General registry** (only 0.2.7 is), so each decoupled package
needs a `[sources]` override to resolve during this phase.

**Use the git URL, not a sibling path.** A `{path = "../ExTinyMD.jl"}` pin resolves locally but
cannot work in CI, where only the one repository is checked out and `actions/checkout` will not
write outside the workspace. Verified on ParticleMeshEwald: with a path pin, `Pkg.test()` fails
with `expected package ExTinyMD [fec76197] to exist at path .../ExTinyMD.jl`. Use
`{url = "https://github.com/HPMolSim/ExTinyMD.jl", rev = "main"}`.

A package pinned to an **unmerged branch** of a sibling (SoEwald2D needs QuasiEwald, whose
decoupling is on `decouple-extinymd`) must pin that branch by name and be updated to
`rev = "main"` once it lands. That is a second, per-package gate on top of the registration
gate below.

That is not merely inconvenient for outside users. **General's automerge rejects any package
whose `Project.toml` carries a `[sources]` section**, so no package can be tagged or
registered while the override is present. Five inline comments will not reliably be
remembered, so it is recorded here as one coordinated item:

- [x] **Register ExTinyMD 0.3 — DONE.** `v0.3.0` is tagged at `0c2c5ff1` (the PR #13 merge) and
      resolves from General. Verified by execution, not by the tag's existence: a scratch
      environment doing `Pkg.add(PackageSpec(name="ExTinyMD", version=v"0.3"))` installs it, and
      the installed copy has `Ewald2D`, `coulomb_energy`, `PME3D` and `ICMEwald2D` defined.

      **That verification was necessary, not ceremonial.** This package's `Project.toml` has read
      `0.3.0` since 2023 (PR #7), long before any of the electrostatics work. Had a `0.3.0` been
      registered back then, every downstream `ExTinyMD = "0.3"` bound would have silently
      resolved to a version with none of the Phase 1/2 API — a failure that looks like a
      dependency resolving fine and then `UndefVarError` at first use. It did not happen, but
      the check is the only thing that distinguishes the two cases.
- [ ] Repoint SoEwald2D's QuasiEwald pin from `rev = "decouple-extinymd"` to `rev = "main"`
      once QuasiEwald PR #5 lands
- [x] Remove `[sources]` from ParticleMeshEwald — done on `main`; resolves from the registry
      (PR #8 merged 2026-09-15; repo is `flatironinstitute/ParticleMeshEwald.jl`, **not**
      HPMolSim, and it is the one package of the five **not** registered in General)
- [x] Remove `[sources]` from QuasiEwald — done on `decouple-extinymd` (PR #5; 0.2.1 → 0.3.0).
      It pinned only ExTinyMD, so the whole block went.
- [ ] Remove `[sources]` from SoEwald2D (0.1.5 → 0.2.0) — **still needed**, and it is now the
      only one. Not for ExTinyMD but for **QuasiEwald**, whose 0.3.0 is unregistered while PR #5
      is open. So SoEwald2D's gate is no longer ExTinyMD's registration but QuasiEwald's: pin
      `rev = "decouple-extinymd"`, repoint to `rev = "main"` when #5 lands, then remove the block
      once QuasiEwald 0.3.0 is tagged and registered.
- [x] FastSpecSoG needs **no** `[sources]` at all (0.1.0 → 0.2.0) — ExTinyMD was its only
      unregistered dependency, so the registration removed the need before it was ever added
- [ ] ~~Remove `[sources]` from EwaldSummations~~ — out of scope, never decoupled, stays on
      CellListMap 0.9 and ExTinyMD 0.2
- [ ] Confirm each resolves from the registry with no local path and no git URL

### Landed so far

| | state |
|---|---|
| ExTinyMD Phase 1 (electrostatics stdlib) | merged, PRs #11/#12 |
| ExTinyMD Phase 2 (PME3D/ICMPME3D via FINUFFT) | merged, PR #13 |
| ParticleMeshEwald (Phase 3a) | merged, PR #8 |
| QuasiEwald (Phase 3b) | PR #5 open |

### A CI trap worth one line

**Check each package's CI matrix against the `julia` compat floor you just raised.** QuasiEwald's
matrix tested `'1.9'` while the decoupling raised its floor to `1.10`, so a job failed for no
reason but the mismatch. Prefer `'lts'` + `'1'` + `'nightly'`, which tracks the floor
automatically, as ExTinyMD's own CI does. SoEwald2D (`'1'`) and FastSpecSoG
(`'1.10'`/`'1.11'`/`'nightly'`) were already fine.

**As of ExTinyMD 0.3.0's registration this is nearly cleared.** ParticleMeshEwald, QuasiEwald
and FastSpecSoG now resolve entirely from the registry. Only SoEwald2D still carries a pin, and
only because it test-depends on QuasiEwald — so the chain of gates is now
`QuasiEwald PR #5 → tag → register → SoEwald2D's pin can go`.

## 7. Sequencing

**Order: QuasiEwald → SoEwald2D → FastSpecSoG.** ParticleMeshEwald is done; EwaldSummations
is out of scope.

This is the reverse of ordering by coupling weight, and the reason is the inter-package test
dependency graph, checked rather than assumed:

| package | test target depends on | consequence |
|---|---|---|
| QuasiEwald | `Test` only | nothing blocks it — goes first |
| SoEwald2D | **QuasiEwald** (`IcmSys`, `IcmSysInit`, used as an ICM reference in 6 places) | needs a decoupled QuasiEwald on CellListMap 0.10, so it goes second |
| FastSpecSoG | **EwaldSummations** (`Ewald2DInteraction` and friends, as accuracy references) | out of scope now, so this must be swapped for ExTinyMD's `Ewald2D` — goes last |

The order is forced, but **not** by CellListMap, which is what an earlier draft of this section
claimed. Decoupling QuasiEwald showed that it never called CellListMap at all — it only ever
consumed neighbor lists built by ExTinyMD's own `CellListQ2D`/`CellListDirQ2D` — so the
dependency was dead and has been dropped outright (its Task-1 "bump to 0.10" was cosmetic).

The real constraint is the **ExTinyMD version**. A decoupled package requires ExTinyMD 0.3, and
a not-yet-decoupled package caps it at 0.2. So SoEwald2D cannot pick up QuasiEwald 0.3 in its
test environment while its own `[compat]` still says `ExTinyMD = "0.2"`: the two bounds have no
common solution. SoEwald2D must therefore move to ExTinyMD 0.3 in the same change that picks up
the decoupled QuasiEwald. Same ordering, different mechanism — and the mechanism matters,
because a package that never depended on CellListMap would otherwise look unblocked when it is
not.

Worth checking per package for the same reason: **a CellListMap bump may be dead work.** Grep
for `InPlaceNeighborList`/`neighborlist!`/`update!` before assuming the 0.9 → 0.10 migration has
any call sites to fix.

The cost is that the §4.3a wrapper pattern gets its first `simulate!` exercise on the
heaviest package rather than the lightest. That is acceptable because the pattern's
*soundness* is already established — ParticleMeshEwald's reviewer verified by execution that
a struct defined inside an extension module can subtype `ExTinyMD.AbstractInteraction`, and
that `PME <: AbstractInteraction` is false even with the extension loaded. What remains
untested is ergonomics under load, not whether it works.

QuasiEwald also has a second relocation the others do not: its `SortingFinder` subtypes
`ExTinyMD.AbstractNeighborFinder`, so it faces the same problem as the interaction types and
moves to the extension alongside them. Learning that on the first package rather than the
last is a small compensation for the reordering.

### 7a. Phase 3b outcome, and what it hands Phase 3c

QuasiEwald is done (branch `decouple-extinymd`, 6 commits, suite 3874/3874 verified by the
controller independently of the implementer's report). Three things it learned that change the
work for SoEwald2D and FastSpecSoG:

1. **§4.3a-bis** — the dispatcher-function pattern, which preserved all three interaction/finder
   names. Check for type-position uses before reusing it.
2. **§4.3b** — registered packages need a minor version bump. SoEwald2D and FastSpecSoG are both
   registered, so both need one.
3. **`import ExTinyMD`, not `using`** — otherwise a bare `function energy(...)` will not compile
   while coupled and decoupled code coexist in one module.

**SoEwald2D needs a `[sources]` override for QuasiEwald, not only for ExTinyMD.** Verified: its
test target lists QuasiEwald with no `[compat]` bound, so it resolves the *registry* QuasiEwald
0.2.x, which is pinned to CellListMap 0.9. The moment SoEwald2D moves to 0.10 that resolve
fails. This adds QuasiEwald to the set of git-URL pins its `Project.toml` carries, and hence to
the §6b removal checklist.

A second, smaller trap in the same file: **SoEwald2D declares its test environment twice** —
once as `[extras]`/`[targets]` in the top-level `Project.toml` and once as a separate
`test/Project.toml`. Julia honours `test/Project.toml` when present and ignores the
`[extras]`/`[targets]` pair, so edits made only to the latter will appear to do nothing. Decide
which one survives before starting, and put the `[sources]` block where it is actually read.

Each package is independently mergeable and leaves its repository working, so the phase can
stop cleanly after any of them.

## 8. Open question for the user: ParticleMeshEwald's future

Phase 2 gave ExTinyMD a `PME3D` that computes the same sum as `Ewald3D` to 4.2e-16, with
forces, ICM integration and an independent finite-difference check. `ParticleMeshEwald.jl` is
energy-only, has no forces at all, and its `[weakdeps]` have never been active.

So after this phase it would be a decoupled, AoS-interfaced package whose entire function is
available in ExTinyMD, better tested and with forces. The same question the user already
answered for EwaldSummations — thin it to what is genuinely its own — has an awkward answer
here, because once the PME is removed there is nothing left.

**This phase does not decide that.** The user asked for these packages to work with ExTinyMD
and to have a simpler interface, and that is what Phase 3 delivers for all five. The question
of whether ParticleMeshEwald should continue to exist is theirs, and is better asked with the
decoupled version in front of them than in the abstract.

## 9. Success criteria

1. All five packages resolve against ExTinyMD 0.3 / CellListMap 0.10 and their suites pass.
2. Each has a passing smoke test that evaluates energy and force from plain arrays in a
   session where ExTinyMD was never loaded.
3. ExTinyMD appears in `[weakdeps]` — not `[deps]` — in all five, with a working `ext/`.
4. EwaldSummations and FastSpecSoG can drive `simulate!`, which neither could before.
5. EwaldSummations' superseded k-space code is removed only after the migration check proves
   ExTinyMD's stdlib reproduces it.
6. No numerical result changes except where a Phase 1/2 fix deliberately corrected one, and
   every such case is recorded.
