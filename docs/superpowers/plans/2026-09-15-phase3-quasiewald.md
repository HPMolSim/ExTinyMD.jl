# Phase 3b: QuasiEwald — Decoupling Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `QuasiEwald.jl` usable from plain arrays with no ExTinyMD type constructed, while keeping it fully usable *inside* ExTinyMD — including driving `simulate!`, which ParticleMeshEwald could not exercise.

**Architecture:** The pattern established in the phase spec §4.3a: `src/` holds framework-free *plans*, `ext/` holds thin *wrappers* subtyping `ExTinyMD.AbstractInteraction` / `AbstractNeighborFinder` and holding a plan, plus the adapter methods.

**Tech Stack:** Julia 1.10+, CellListMap 0.10, StaticArrays, GaussQuadrature, SpecialFunctions.

**Spec:** `docs/superpowers/specs/2026-09-15-downstream-decoupling-design.md` — this is §5.5, now sequenced **first** among the remaining packages (see §7).

**Repository:** `/mnt/home/xgao1/project/q2dmd/QuasiEwald.jl` — **a different repository from the one this plan lives in.** ExTinyMD is at `../ExTinyMD.jl`.

## Why this one goes first

Ordering by coupling weight would put QuasiEwald last — it is the heaviest at 1789 src LOC. The test-dependency graph overrides that:

| package | test target needs | consequence |
|---|---|---|
| **QuasiEwald** | `Test` only | nothing blocks it |
| SoEwald2D | QuasiEwald's `IcmSys`/`IcmSysInit`, in 6 places | needs QuasiEwald on CellListMap 0.10 first |
| FastSpecSoG | EwaldSummations (now out of scope) | must swap to ExTinyMD's `Ewald2D` |

Anything left on CellListMap 0.9 breaks the resolve of anything depending on it that has moved to 0.10.

This also means QuasiEwald is where the §4.3a wrapper pattern gets its first real exercise. Its *soundness* is settled — ParticleMeshEwald's reviewer verified by execution that a struct defined inside an extension module can subtype `ExTinyMD.AbstractInteraction`, and that retrofitting one onto a `src/`-defined type is impossible. What is untested is how it feels in a running MD loop. That is what this package must demonstrate.

## Measured starting state

```
src: 13 files, 1789 LOC       test: 510 LOC       exports: 50
```

ExTinyMD API surface, counted:

| symbol | uses | disposition |
|---|---|---|
| `Point` | ~70 lines across 5 files | → `SVector{3,T}` |
| `SimulationInfo` | 11 | → extension only |
| `MDSys` | 8 | → extension only |
| `ExTinyMD.AbstractNeighborFinder` | 3 | `SortingFinder` → extension |
| `CellListQ2D` | 3 | → extension only |
| `position_checkQ2D` | 2 | → local minimum-image helper |
| `ExTinyMD.update_acceleration!` | 2 | → extension |
| `ExTinyMD.energy` | 2 | → extension |
| `ExTinyMD.AbstractInteraction` | 2 | both interaction types → extension wrappers |
| `dist2` | 2 | → `sum(abs2, a - b)` |
| `ExTinyMD.update_finder!` | 1 | → extension |

`Point` concentration: `force/force_long.jl` 33, `tools/Icm.jl` 15, `energy/energy_long.jl` 11, `types.jl` 4, `force/force_short.jl` 4.

`StaticArrays` is **already a dependency**, so `SVector` needs no new dep.

### The three types, and what each actually needs

Read `src/types.jl:142-230` before starting. The situation is better than the LOC count suggests:

- **`QuasiEwaldShortInteraction`** (types.jl:142) — fields are *pure parameters*: `γ_1, γ_2, ϵ_0, L, rbe, accuracy, α, n_atoms, r_c, n_t, gauss_para`. **Nothing ExTinyMD-typed in it.** Only its supertype couples it. It is already a plan in all but name.
- **`QuasiEwaldLongInteraction`** (types.jl:161) — parameters plus scratch, and the scratch is the interesting part: `q`, `mass`, `coords::Vector{Point{3,T}}`, `acceleration::Vector{Point{3,T}}`. `Point` becomes `SVector`. But `mass` and `acceleration` are **MD concerns leaking into what should be a plan** — see Task 3.
- **`SortingFinder`** (types.jl:208) — fields `z_coords`, `z_list` are clean data. Only its *constructor* takes `SimulationInfo`, and `ExTinyMD.update_finder!` is defined on it at types.jl:219.

## Global Constraints

- `julia = "1.10"` minimum, no newer syntax.
- **CellListMap 0.10 keyword API**: `InPlaceNeighborList(xpositions = ...)`, `update!(cl, xpositions = ...)`, `neighborlist(xpositions = ..., ...)`. The 0.9 spellings raise `MethodError`. `../ExTinyMD.jl/src/MD_core/neighbor_finder/cell_list.jl` is the reference.
- **`r_c` strictly less than half the smallest *periodic* box side.** Under the quasi-2D geometry only `L[1]` and `L[2]` count — z is not periodic. CellListMap 0.10 raises `ArgumentError: UNIT CELL CHECK FAILED` otherwise. **Annotate every configuration you write with its computed `r_c` inline.**
- **ExTinyMD in `[weakdeps]` + `[extensions]`, never `[deps]`**, plus `[extras]`/`test` target so `Pkg.test()` exercises the extension.
- **`[sources]` must use the git URL**, not a sibling path: `ExTinyMD = {url = "https://github.com/HPMolSim/ExTinyMD.jl", rev = "main"}`. A path pin resolves locally but fails in CI, where only this repository is checked out — verified on ParticleMeshEwald, where the path form gave `expected package ExTinyMD [fec76197] to exist at path .../ExTinyMD.jl`. Copy ParticleMeshEwald's comment block explaining the registration gate.
- **`energy`, `force`, `force!` defined but NOT exported.** Five packages each exporting `energy` would collide. Plan constructors stay exported.
- **Never mutate the caller's `poses` or `charges`.** Scatter into plan-owned scratch.
- Positions accessed only via `p[1]`, `p[2]`, `p[3]`.
- No `Threads.threadid()` accumulator indexing. **Grep for `threadid` first** — under `@inbounds`, index 0 is an undefined-behaviour write producing a silently wrong answer, not a `BoundsError`. That was ParticleMeshEwald's bug.
- **Preserve the physics.** This phase is restructuring. A number that moves is a finding to report, not an expected value to update.
- Work on a branch in the QuasiEwald repository. Commit after every task.

## Baseline capture — read this before Task 1

**`ExTinyMD.SimulationInfo` consumes `rand()` internally.** A before/after comparison script that constructs one shifts the RNG stream, so the "before" numbers are not comparable to the "after" ones. ParticleMeshEwald's implementer hit exactly this and had to re-derive from unmodified source via `git stash`.

So: seed immediately before each measurement, or build your comparison configurations without `SimulationInfo` at all. Record the baseline **before touching anything** and state in your report how you avoided contamination.

---

### Task 1: Compat — CellListMap 0.10, and expect a resolve failure first

- [ ] **Step 1:** Record baseline energies and forces, per the warning above. Several configurations, spanning the short-range, long-range and ICM paths. Put the numbers in your report.

- [ ] **Step 2:** In `Project.toml`: `CellListMap = "0.9"` → `"0.10"`, `ExTinyMD = "0.2"` → `"0.3"`, and check every stdlib compat entry — ParticleMeshEwald had `LinearAlgebra = "1.12.0"` pinning a stdlib to a version requiring Julia 1.12 while claiming 1.10 support. Widen `julia` to `"1.10"`.

- [ ] **Step 3: Expect the first failure to be a resolve, not a compile.** ParticleMeshEwald's plan predicted a `MethodError` from the renamed keyword; the actual first failure was an unsatisfiable `Pkg.test()`, because published packages in the test target still pin CellListMap 0.9. QuasiEwald's test target is only `["Test"]`, so it may be spared — but check, and if a resolve fails, deal with the pin before the API.

- [ ] **Step 4:** Fix the CellListMap call sites. Grep for `InPlaceNeighborList`, `update!`, `neighborlist`.

- [ ] **Step 5:** Suite passes **and** baseline numbers unchanged. Report both sets.

- [ ] **Step 6:** Commit.

---

### Task 2: `Point` → `SVector{3,T}`

~70 lines across five files. Mechanical, but the largest single change in the package, so do it as its own commit with the suite green on both sides.

`SVector{3,T}` supports everything `Point` does — arithmetic, iteration, indexing. Two things do not carry over:

- **`dist2`** is ExTinyMD's. `dist2(a, b)` → `sum(abs2, a - b)`; `dist2(a)` → `sum(abs2, a)`. Two call sites.
- **`position_checkQ2D`** is ExTinyMD's, and it has a trap worth knowing: it scans `m ∈ -1:1` and returns the *first* image inside the cutoff together with a `dist_sq` of exactly zero when none is, so every caller carries an `iszero` guard. Phase 1 replaced it inside ExTinyMD with a true nearest-image helper, `dx - L*round(dx/L)` per periodic axis, which is correct for any cutoff rather than only `r_c < L/2`. Write a local equivalent — quasi-2D wraps x and y only, leaving z as a plain difference — and drop the sentinel-zero dance. Two call sites.

- [ ] Substitute, run the suite, confirm the baseline numbers are unchanged, commit.

**Report any site where the substitution was not mechanical.** Those are where a latent assumption about `Point` lived, and they are worth knowing about.

---

### Task 3: Split the plans from the MD concerns

**Files:** `src/types.jl`, `src/energy/`, `src/force/`

This is the task with design judgement in it.

`QuasiEwaldShortInteraction` is already pure parameters — dropping `<: ExTinyMD.AbstractInteraction` makes it a plan, nothing more.

`QuasiEwaldLongInteraction` is not: it carries `mass::Vector{T}` and `acceleration::Vector{Point{3,T}}`. A framework-free core should compute **forces** and leave mass to the caller — mass is an MD property, and dividing by it inside the plan is what makes this type MD-specific rather than a solver. The core's output should be a force buffer; the extension's `update_acceleration!` divides by mass and accumulates, exactly as ExTinyMD's own adapter does (`../ExTinyMD.jl/src/interactions/electrostatics/adapter.jl`).

**Decide and justify in your report:**

1. **Naming.** The cleanest option is that `src/` gets `QuasiEwaldShortPlan` / `QuasiEwaldLongPlan` and `ext/` keeps the *existing* names `QuasiEwaldShortInteraction` / `QuasiEwaldLongInteraction` as the wrappers — so existing MD-side code keeps compiling unchanged, and only the standalone path is new. Alternatively `src/` keeps the old names and `ext/` invents new ones. **Pick one and say why**; whichever you choose, the same convention should hold for SoEwald2D and FastSpecSoG afterwards.
2. **Whether `mass`/`acceleration` leave the long plan.** I believe they should, per the argument above, but you are looking at the code. If removing them costs more than it buys — for instance if `force_long_total!` is written to accumulate accelerations in a way that is awkward to restructure — say so with specifics rather than doing a half-conversion.
3. **`SortingFinder`.** Its fields are clean; its constructor takes `SimulationInfo` and `ExTinyMD.update_finder!` is defined on it. Split it: a plan-side sorter constructed from `poses` (or from a z-coordinate vector) in `src/`, and the `AbstractNeighborFinder` wrapper plus `update_finder!` in `ext/`.

- [ ] Write the core queries: `QuasiEwald.energy(plan, poses, charges)`, `QuasiEwald.force(plan, poses, charges)`, `QuasiEwald.force!(F, plan, poses, charges)`. Short and long are separate plans, as they are separate interactions today.
- [ ] Suite green, baseline unchanged, commit.

---

### Task 4: The extension

**Files:** `ext/QuasiEwaldExTinyMDExt.jl` (new), `Project.toml`

- [ ] Define the wrappers subtyping `ExTinyMD.AbstractInteraction` and `ExTinyMD.AbstractNeighborFinder`, holding plans.
- [ ] Move `ExTinyMD.energy` (×2), `ExTinyMD.update_acceleration!` (×2) and `ExTinyMD.update_finder!` (×1) onto the wrappers. Honour the id/slot indirection: `sys.atoms` is indexed by particle **id**, `info.particle_info` by storage **slot**.
- [ ] `[weakdeps]`, `[extensions]`, `[extras]`, `test` target, `[sources]` git URL.

- [ ] **The test that matters: drive `simulate!` with the wrappers in `sys.interactions`.**

ParticleMeshEwald has no forces, so its adapter test only called `ExTinyMD.energy` directly — the wrapper pattern has never run inside an MD loop. This test is the first one that exercises it, and it is what catches the id/slot, mass-division and accumulation faults that per-call tests cannot. Run a short microcanonical trajectory and assert bounded energy drift.

Model it on `../ExTinyMD.jl/test/electrostatics/test_adapter.jl`, which has both the `simulate!` drift test and the permuted slot/id test — the latter caught a real bug in ExTinyMD's own `SubLennardJones`, where a storage slot was used to index the id-keyed `sys.atoms`. Consider porting that permutation test too; it is cheap and it is the only thing that catches that class of fault.

- [ ] Commit.

---

### Task 5: The standalone test, and documentation

- [ ] **Step 1:** The test the whole phase exists for — a subprocess asserting ExTinyMD is absent from `Base.loaded_modules`, then evaluating an energy and a force from plain arrays. ParticleMeshEwald's `test/standalone.jl` is the working model; copy its shape, including the UUID check.

**Verify it fails when you add `using ExTinyMD` to the subprocess script.** Without that check the test proves nothing, and across these phases eleven specifications have turned out to be tests that passed for the wrong reason.

- [ ] **Step 2:** README — the AoS API, the un-exported query names, the `r_c < min(Lx,Ly)/2` rule, the wrapper types for MD use, and the `[sources]`/registration note.

- [ ] **Step 3:** Full suite, commit.

---

## Self-Review

**Spec coverage.** §4.1 core → Task 3. §4.2 naming → Global Constraints, Task 5. §4.3 adapter → Task 4. §4.3a wrapper → Task 3 decision 1 and Task 4. §4.4 `Project.toml` → Tasks 1 and 4. §6 testing: standalone → Task 5; adapter under `simulate!` → Task 4; before/after numerical → Tasks 1–3, where recording the baseline is an explicit step in each.

**What is genuinely uncertain.** Task 3's three decisions, especially whether `mass`/`acceleration` leave the long plan. The rest is mechanical substitution plus relocation, on a package whose suite already passes.

**Delegated with justification required:** plan/wrapper naming (sets the convention for two more packages), whether `mass`/`acceleration` leave the long plan, and how `SortingFinder` splits.
