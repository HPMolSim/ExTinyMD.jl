# Phase 3d: Decoupling FastSpecSoG Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task.

**Goal:** Make `FastSpecSoG.jl` usable with no ExTinyMD dependency, and replace its unmaintained `EwaldSummations` test reference with ExTinyMD's own `Ewald2D`.

**Architecture:** Two-layer contract. `src/` holds framework-free plans queried from AoS positions + charges; `ext/FastSpecSoGExTinyMDExt.jl` supplies `ExTinyMD.energy`. `ExTinyMD` moves to `[weakdeps]`.

**Tech Stack:** Julia 1.10+, ChebParticleMesh, FastChebInterp, FFTW, DoubleFloats, LoopVectorization.

**Spec:** `docs/superpowers/specs/2026-09-15-downstream-decoupling-design.md` — §4.1, §4.2, §4.3, §4.3a, §4.3a-bis, §4.3b, §5.3, §5.4, §6a, §6b, §7, §7a.

## Global Constraints

- `ExTinyMD` in `[weakdeps]`, NEVER `[deps]`. `src/` must contain **zero ExTinyMD references in code** (comments/docstrings fine). Grep to verify after every task.
- Positions are AoS and this package is **already** AoS-native: `Vector{NTuple{3,T}}`. Keep read-only position parameters **untyped** so `NTuple{3,T}`, `SVector{3,T}` and ExTinyMD's `Point{3,T}` all work by `p[1]`/`p[2]`/`p[3]` indexing alone. Do NOT introduce a conversion layer, and do NOT churn `NTuple` → `SVector` where nothing requires it.
- Queries never mutate caller-supplied positions/charges.
- `r_c` strictly less than `min(Lx, Ly)/2`, annotated inline at every test site.
- **Version bump required** (§4.3b): FastSpecSoG is registered in General. 0.1.0 → **0.2.0**.
- Do not change the physics. Any change to a computed value must be reported, not absorbed.

---

## What makes this package different from the previous three

**Good news first.** Grep finds **no `Point` usage anywhere in `src/`** — positions are already
`Vector{NTuple{3,T}}`. The single largest task in the other three packages does not exist here.

**It is energy-only.** There is no `update_acceleration!` and no force code at all. So, exactly
as with ParticleMeshEwald (§5.1), it can never drive `simulate!` regardless of this work, and
the extension supplies `ExTinyMD.energy` alone. You do **not** need the §4.3a-bis
dispatcher-function pattern unless you conclude otherwise — see Decision 1.

**Its `ExTinyMD.energy` methods have the wrong signature and are already dead.** They are
declared as `ExTinyMD.energy(interaction, neighbor, info::SimulationInfo, atoms::Vector{Atom})`
(`src/energy/energy.jl:12,40`). ExTinyMD 0.3's MD loop calls
`energy(interaction, neighborfinder, sys::MDSys, info::SimulationInfo)` — a different arity and
order. So these methods can never be dispatched by `simulate!`; only this package's own tests
call them, directly. **Verify this yourself** before relying on it (read
`../ExTinyMD.jl/src/MD_core/recorder/energy_logger.jl:51` and grep for other call sites), then
say so in your report: it changes what "keep MD integration working" means here, because there
is no working MD integration to preserve.

**The `boundary::Boundary{T}` struct field is redundant.** It appears in `FSSoGInteraction` and
`FSSoGThinInteraction` (`src/types.jl:30,76`) but is only ever constructed internally as
`Q2dBoundary(L...)` (`src/FSSoGInteraction.jl:54`) from an `L` the struct already carries. So
it can be dropped outright rather than abstracted — no information is lost. Confirm by grep
that nothing reads it except the `position_check3D` calls you are replacing.

---

## THE MIN-IMAGE REPLACEMENT — read before Task 2

`position_check3D` is called with a `Q2dBoundary`, which is `Boundary((Lx,Ly,Lz), (1,1,0))`.
The z period is 0, so its `mz` loop runs only at `mz = 0`. Effective behaviour:

- wrap **x and y** to the nearest periodic image under `Lx`/`Ly`;
- leave **z** alone (slab axis, not periodic);
- return the **FULL 3-D squared distance** `dx² + dy² + dz²`.

Use exactly this helper (identical to the one Phase 3c adds to SoEwald2D — the geometry is the
same slab). Do **not** substitute an in-plane-only distance: this package's short-range sum is
a function of the true separation `r`, and its cutoff test is on `r`.

```julia
@inline _wrap(dx::T, L::T) where {T} = dx - L * round(dx / L)

"""
    _min_image_slab(pos_i, pos_j, L) -> (coord_i, coord_j, r_sq)

Slab-geometry nearest image: `x`/`y` wrap under `L[1]`/`L[2]`, `z` is a plain difference
(the slab axis is not periodic). Returns `pos_i` shifted to its nearest in-plane image of
`pos_j`, `pos_j` unchanged, and the **full three-dimensional** squared distance.
Inputs need only support `p[1]`/`p[2]`/`p[3]` indexing.
"""
@inline function _min_image_slab(pos_i, pos_j, L::NTuple{3, T}) where {T}
    dx = _wrap(T(pos_i[1]) - T(pos_j[1]), L[1])
    dy = _wrap(T(pos_i[2]) - T(pos_j[2]), L[2])
    dz = T(pos_i[3]) - T(pos_j[3])
    r_sq = dx^2 + dy^2 + dz^2
    coord_i = SVector{3, T}(T(pos_j[1]) + dx, T(pos_j[2]) + dy, T(pos_i[3]))
    coord_j = SVector{3, T}(T(pos_j[1]), T(pos_j[2]), T(pos_j[3]))
    return coord_i, coord_j, r_sq
end
```

(If you prefer to return `NTuple{3,T}` here to match this package's existing convention and
avoid adding a StaticArrays dependency, that is fine and arguably better — decide, and say
which you did. The returned coords are only ever indexed.)

**The sentinel-to-explicit-test trap.** `position_check3D` returns an all-zero sentinel when no
image is in range, so callers guard with `iszero(r_sq)`. That guard was *also* incidentally
skipping `r == 0` pairs. Replacing it with an explicit `r_sq ≥ r_c^2` test keeps them, and in
QuasiEwald that shipped a `0/0` NaN reaching production, caught only in final review. Guard
`r == 0` explicitly here and add a test asserting finiteness. Since this is energy-only there
may be no division by `r` at all — **check** rather than assume, and report which.

---

## File Structure

| File | Responsibility |
|---|---|
| `Project.toml` | `ExTinyMD` → `[weakdeps]`+`[extensions]`; `[sources]` for ExTinyMD; drop `EwaldSummations` from test extras; version 0.2.0 |
| `src/FastSpecSoG.jl` | drop `using ExTinyMD`; update exports |
| `src/types.jl` | drop the three `<: ExTinyMD.AbstractInteraction` supertypes; drop the `boundary` fields; add `_min_image_slab` |
| `src/FSSoGInteraction.jl` | constructors stop building a `Boundary` |
| `src/energy/energy_short.jl`, `energy_short_naive.jl` | `position_check3D` → `_min_image_slab`; drop `Boundary`/`CellList3D` annotations |
| `src/energy/energy.jl` | ExTinyMD-coupled entry points move to `ext/` |
| `ext/FastSpecSoGExTinyMDExt.jl` | NEW — `ExTinyMD.energy` with ExTinyMD 0.3's real signature |
| `test/runtests.jl`, `test/energy.jl` | EwaldSummations → ExTinyMD `Ewald2D` |
| `test/standalone.jl` | NEW — proves ExTinyMD is never loaded |

---

## Task 1: Baseline, compat, `[sources]`

- [ ] **Step 1: capture a numerical baseline first.** Record to full precision, for at least
three parameter sets: `energy_naive`, `energy_short`, `energy_mid`, `energy_long`,
`energy_long_naive`, `short_energy_naive`, `short_energy_Cheb`, `energy_per_atom`, and the
`FSSoGThinInteraction` path. Seed once with `Random.seed!` and **build positions/charges
directly** — do not construct a `SimulationInfo`, which consumes `rand()` internally and makes
a baseline non-reproducible across changes that alter RNG draw counts. Re-run after every task
and report bit-identity or the exact deviation.

- [ ] **Step 2:** `ExTinyMD` → `[weakdeps]` + `[extensions]`:

```toml
[weakdeps]
ExTinyMD = "fec76197-d59f-46dd-a0ed-76a83c21f7aa"

[extensions]
FastSpecSoGExTinyMDExt = "ExTinyMD"
```

Version `0.1.0` → `0.2.0`. `[compat]`: `ExTinyMD = "0.3"` (was `0.2.6`).

- [ ] **Step 3:** `[sources]` git-URL pin for ExTinyMD (0.3 unregistered; General has 0.2.7).
Use the git URL, not a sibling path — a path pin resolves locally but breaks CI, where only
this repo is checked out. Copy ParticleMeshEwald's comment block about General's automerge
rejecting `[sources]`.

```toml
[sources]
ExTinyMD = {url = "https://github.com/HPMolSim/ExTinyMD.jl", rev = "main"}
```

- [ ] **Step 4:** remove `EwaldSummations` from `[extras]`/`[targets]`, add `ExTinyMD` to the
test target. EwaldSummations is out of scope for this project and stays on CellListMap 0.9, so
it cannot co-resolve with ExTinyMD 0.3.

- [ ] **Step 5:** `Pkg.test()`, record the assertion count, re-run the baseline, commit.
Expect the suite to FAIL at this point because the tests still reference EwaldSummations —
that is fine and expected; Task 4 fixes it. Record the failure rather than papering over it.

---

## Task 2: Drop the ExTinyMD types from `src/`

**Files:** `src/types.jl`, `src/FSSoGInteraction.jl`, `src/energy/energy_short.jl`, `src/energy/energy_short_naive.jl`, `src/FastSpecSoG.jl`

- [ ] **Step 1:** add `_min_image_slab` + `_wrap` (verbatim from the section above) to
`src/types.jl`.

- [ ] **Step 2:** drop `<: ExTinyMD.AbstractInteraction` from `FSSoG_naive`,
`FSSoGInteraction`, `FSSoGThinInteraction`. They become plain structs — the plan layer.

**Keep all three names unchanged.** Unlike the other packages, nothing here needs the
§4.3a-bis dispatcher pattern, because these types do not need an ExTinyMD supertype for
anything that works today (see Decision 1). Keeping the names means every existing test and
downstream call site compiles untouched.

- [ ] **Step 3:** delete the `boundary::Boundary{T}` field from `FSSoGInteraction` and
`FSSoGThinInteraction`, and stop constructing it in `src/FSSoGInteraction.jl`. Every reader of
that field is a `position_check3D` call being replaced in Step 4; confirm by grep first.

- [ ] **Step 4:** replace each `position_check3D(...)` with `_min_image_slab(pos_i, pos_j, L)`
and each `iszero(r_sq)` sentinel guard with an explicit `r_sq ≥ r_c^2` test. Handle `r == 0`
per the trap section above.

- [ ] **Step 5:** remove `Boundary{T}` and `CellList3D{T}` type annotations from every
signature. The neighbor-list parameter becomes untyped (or a keyword `neighbor_list = nothing`
per Task 3) so any `(i, j, ...)`-yielding iterable works. **Always recompute the true distance
from the positions** — never trust a supplied list's reported distance. That is the Phase 1
`Ewald2D` bug: `short_energy` trusted a `CellListQ2D` list's `r`, which was in-plane, and
returned +0.0238 against a true −0.1539.

- [ ] **Step 6:** drop `using ExTinyMD` from `src/FastSpecSoG.jl`. Suite + baseline + commit.
The `Point`-free nature of this package means Step 2-5 should be bit-identical; if not, report
the deviation.

---

## Task 3: The AoS query API

**Files:** `src/energy/*.jl`, `test/plan.jl` (new)

- [ ] **Step 1:** give the plans a uniform, documented query surface taking positions and
charges as **arguments** rather than reading them from a struct field:

```julia
FastSpecSoG.energy(plan, poses, charges; neighbor_list = nothing)
```

`FSSoGInteraction`/`FSSoGThinInteraction` currently store `position` and `charge` as fields
and the energy functions read them. Keep those as plan-owned **scratch** and have the query
scatter the caller's arrays into them (the ParticleMeshEwald `_scatter_long!` pattern), so the
caller's data is never mutated and nothing is allocated per call. Preserve the existing
field-reading functions as internal helpers if that keeps the diff small — say what you chose.

`energy` is defined but **NOT exported** (§4.2): with five sibling packages each exporting an
`energy`, `using FastSpecSoG, SoEwald2D` would make the bare name ambiguous. Keep the existing
descriptive exports (`energy_short`, `energy_mid`, `energy_long`, `energy_naive`, …) exactly
as they are — those are not ambiguous and renaming them is out of scope.

- [ ] **Step 2:** `test/plan.jl` — construct the plans and query them from plain arrays with
no ExTinyMD type anywhere; assert against `energy_naive`/`long_energy_naive` as the internal
reference; assert positions/charges are not mutated; assert finiteness at `r = 0`.

- [ ] **Step 3:** suite + baseline + commit.

---

## Task 4: Swap EwaldSummations for ExTinyMD's `Ewald2D`

**Files:** `test/runtests.jl`, `test/energy.jl`

This is a strict improvement, not just a dependency removal: EwaldSummations is unmaintained
and pinned to CellListMap 0.9, while ExTinyMD's `Ewald2D` is the Phase 1 stdlib implementation
with 1128 tests behind it, validated to ~5e-9 α-independence and against ICM to 4.7e-8.

- [ ] **Step 1:** map the old API onto the new one. Old (`test/energy.jl:20-27, 86-93`):

```julia
Ewald2D_interaction = Ewald2DInteraction(n_atoms, 5.0, 0.25, (L, L, L), ϵ = 1.0)
Ewald2D_neighbor = CellList3D(info, Ewald2D_interaction.r_c, boundary, 1)
energy_ewald = energy(Ewald2D_interaction, Ewald2D_neighbor, info, atoms)
energy_ewald_s = Ewald2D_short_energy_N(10, Ewald2D_interaction, p, q)
energy_ewald_l = Ewald2D_long_energy_N(10, Ewald2D_interaction, p, q)
```

New — ExTinyMD 0.3 exports a **framework-free** API, so the reference needs no
`SimulationInfo`, no `MDSys` and no neighbor finder at all:

```julia
using ExTinyMD: Ewald2D, Ewald2DLong, EwaldShort, coulomb_energy, long_energy, short_energy

ewald = Ewald2D(n_atoms, (L, L, L); α = 0.25, s = 5.0, ϵ = 1.0)
energy_ewald = coulomb_energy(ewald, poses, charges)
```

**Note the argument-order change**: the old positional call was `(n_atoms, s, α, L; ϵ)` — `s`
before `α`. The new one is `(n_atoms, L; α, s, ϵ)` with `α` and `s` as **keywords**. Getting
these backwards silently changes the Ewald splitting rather than erroring, and the total energy
is α-independent, so **the total will still be right while the short/long split is wrong** —
which is exactly what the `_short_`/`_long_` assertions below test. Read
`../ExTinyMD.jl/src/interactions/electrostatics/ewald.jl:107` and confirm the signature before
writing the call.

For the split references, use `short_energy(ewald.short, poses, charges)` and
`long_energy(ewald.long, poses, charges)` — check the actual field names on `EwaldInteraction`
in `ewald.jl` rather than guessing.

- [ ] **Step 2:** the old `_N` variants took a truncation count (`10`). ExTinyMD's `Ewald2DLong`
takes `α`/`s` and derives its own k-cutoff (`k_c = 2αs`). So there is no direct `N`
equivalent — construct the reference at an accuracy tight enough that the comparison is
meaningful and **state the tolerance you chose and why**. Do not loosen an existing assertion
to make it pass; if a tolerance must change, report the old and new values with the measured
error.

- [ ] **Step 3:** drop `EwaldSummations` from `test/runtests.jl`'s `using` line.

- [ ] **Step 4:** suite green. Report the assertion count and every tolerance you touched,
with before/after numbers.

---

## Task 5: Extension, standalone proof, docs

**Files:** `ext/FastSpecSoGExTinyMDExt.jl` (new), `src/energy/energy.jl`, `test/standalone.jl` (new), `README.md`

- [ ] **Step 1:** create `ext/FastSpecSoGExTinyMDExt.jl` and move the ExTinyMD-coupled entry
points from `src/energy/energy.jl` into it, **with ExTinyMD 0.3's real signature**:

```julia
ExTinyMD.energy(interaction, neighborfinder, sys::ExTinyMD.MDSys, info::ExTinyMD.SimulationInfo)
```

not the current `(interaction, neighbor, info, atoms)`, which cannot be dispatched by the MD
loop. Gather positions in **slot** order and charges via
`sys.atoms[info.particle_info[slot].id]` — `info.particle_info` is indexed by storage slot,
`sys.atoms` by particle id. Model on
`../ExTinyMD.jl/src/interactions/electrostatics/adapter.jl` and
`../ParticleMeshEwald.jl/ext/`.

Add a test that the adapter is **correct when slot order differs from id order**: reverse
`info.particle_info`, give every id a distinct charge, and check `ExTinyMD.energy` against a
direct plan query on the same positions/charges. An id/slot bug is invisible when every
particle is identical.

- [ ] **Step 2:** `test/standalone.jl` — a subprocess builds the plans from plain arrays,
computes energies, and asserts ExTinyMD is absent from `Base.loaded_modules` throughout. Then
a **second** subprocess loads ExTinyMD first and re-runs the identical assertion, and the test
asserts that process **fails** — capturing its stderr and asserting it failed for the intended
reason (the specific `AssertionError`), not merely that it exited non-zero. `!success(proc)`
alone passes on a missing dependency or a precompile error, proving nothing; that was a review
finding on QuasiEwald.

- [ ] **Step 3:** README — "Standalone usage (no ExTinyMD)" and "MD usage via ExTinyMD"
sections: the AoS API, `energy` unexported while the descriptive names stay exported, the
`r_c < min(Lx,Ly)/2` rule, the no-mutation guarantee, that this package is **energy-only** and
therefore cannot drive `simulate!`, and a `[sources]`/registration note pointing at the design
doc §6b. No Documenter `(@ref)` links.

- [ ] **Step 4:** confirm `src/` has zero ExTinyMD code references, `ExTinyMD` is in
`[weakdeps]` only, and every `r_c` in `test/` carries its inline `min(Lx,Ly)/2` comparison.
Full suite. Commit.

---

## Decisions delegated to you — justify each in writing

1. **Whether to add §4.3a-bis dispatcher wrappers so the three interaction types can sit in
   `sys.interactions`.** My reading: no. The package is energy-only, so it has no
   `update_acceleration!` and `simulate!` would fail on it regardless; and its `ExTinyMD.energy`
   methods already have a signature the MD loop cannot call, so there is no working MD
   integration to preserve. That makes it ParticleMeshEwald's case (§5.1): plans in `src/`, and
   the extension supplies `ExTinyMD.energy` alone. **Verify the signature claim yourself** and
   say whether you agree.
2. **Whether `_min_image_slab` returns `SVector{3,T}` or `NTuple{3,T}`.** This package is
   `NTuple`-native and has no StaticArrays dependency; adding one for two return values may not
   be worth it. Decide and justify.
3. **The `r == 0` behaviour**, and whether any division by `r` even exists on the energy-only
   path.

---

## Self-review checklist

1. Does `src/` mention ExTinyMD anywhere outside a comment or docstring?
2. Is `_min_image_slab` returning the **3-D** distance, with a test that would fail if someone
   swapped in an in-plane version?
3. Was the `Ewald2D` argument order (`α`/`s` as keywords, not positional `s, α`) confirmed
   against the source rather than assumed? Is there an assertion on the short/long **split**,
   not just the α-independent total, so a swapped pair would be caught?
4. Every tolerance you changed reported with before/after measured error?
5. Does any test pass regardless of the code under test?
6. Baseline bit-identical after every task, or the deviation explained?
