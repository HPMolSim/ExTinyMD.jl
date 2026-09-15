# Phase 3: Decoupling the Downstream Electrostatics Packages — Design

Date: 2026-09-15
Status: draft, awaiting review

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

### 5.2 SoEwald2D — lightest of the four coupled packages

823 src LOC, coupling concentrated in `SimulationInfo`/`MDSys` rather than `Point`. Already
defines `ExTinyMD.update_acceleration!`, so the adapter mostly moves rather than being written.

### 5.3 FastSpecSoG

1465 src LOC. Coupling via `CellList3D`, `SimulationInfo`, `Atom`. Uses the **older** interface
generation — `energy(interaction, neighbor, info, atoms)` — and has **no**
`update_acceleration!`, so its adapter is new work rather than a move.

### 5.4 EwaldSummations — thin it, per the user's decision

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

## 7. Sequencing

One plan per package, executed in the §5 order. Each is independently mergeable and leaves the
repository working, so the phase can stop cleanly after any package.

ParticleMeshEwald first is deliberate: it is the cheapest place to discover that something
about the `[weakdeps]`/`[extensions]`/`ext/` arrangement does not work as expected.

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
