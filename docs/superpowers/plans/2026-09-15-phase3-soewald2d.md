# Phase 3c: Decoupling SoEwald2D Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task.

**Goal:** Make `SoEwald2D.jl` usable with no ExTinyMD dependency — a framework-free plan layer queried from AoS position/charge arrays — while keeping `simulate!` integration through a package extension.

**Architecture:** The two-layer contract from the spec. `src/` holds `SoEwald2DShortPlan`/`SoEwald2DLongPlan`, queried via `SoEwald2D.energy`/`force`/`force!`. `ext/SoEwald2DExTinyMDExt.jl` holds `ExTinyMD.AbstractInteraction` wrappers plus the `ExTinyMD.energy`/`update_acceleration!` methods. `ExTinyMD` moves to `[weakdeps]`.

**Tech Stack:** Julia 1.11+, StaticArrays, Enzyme (reverse-mode autodiff, unchanged), GaussQuadrature, SpecialFunctions.

**Spec:** `docs/superpowers/specs/2026-09-15-downstream-decoupling-design.md` — read §4.1, §4.2, §4.3, §4.3a, §4.3a-bis, §4.3b, §6a, §6b, §7, §7a.

## Global Constraints

- `ExTinyMD` in `[weakdeps]`, NEVER `[deps]`. `src/` must contain **zero ExTinyMD references in code** (comments/docstrings are fine). Verify by grep at the end of every task.
- Interchange type: `SVector{3,T}`. Kernels index only `p[1]`/`p[2]`/`p[3]` so `Point{3,T}` and `NTuple{3,T}` also work with no conversion layer. Leave read-only `poses`/`charges` parameters **untyped**; only concretely type freshly-built output/accumulator buffers.
- Queries must never mutate caller-supplied `poses`/`charges`. Scatter into plan-owned scratch instead.
- `r_c` must be strictly less than `min(Lx, Ly)/2`. Annotate every `r_c` in tests with the comparison inline.
- `update_acceleration!` ACCUMULATES (`+=`) into `info.particle_info[i].acceleration`; never overwrites.
- id/slot indirection: `info.particle_info` is indexed by storage **slot**, `sys.atoms` by particle **id**. Gather positions in slot order and look up charge/mass via `sys.atoms[info.particle_info[slot].id]`.
- `energy`/`force`/`force!` are defined but **NOT exported** (§4.2) — always called qualified.
- **Version bump required** (§4.3b): SoEwald2D is registered in General. 0.1.5 → **0.2.0**.
- Do not change the physics. Any change to a computed value must be reported, not absorbed.

---

## THE ONE THING MOST LIKELY TO GO WRONG — read before Task 2

SoEwald2D's short-range code calls ExTinyMD's **`position_check3D`**, and its boundary is
**`Q2dBoundary`**, which is `Boundary((Lx,Ly,Lz), (1,1,0))`. Because the z period is 0,
`position_check3D`'s `mz` loop runs only at `mz = 0`. So in this package it behaves as:

- wrap **x and y** to the nearest periodic image under `Lx`/`Ly`;
- leave **z** alone (the slab axis is not periodic);
- but return the **FULL 3-D squared distance**, `dx² + dy² + dz²`.

That last point is where this will go wrong. QuasiEwald's `_min_image_q2d`, which you may be
tempted to copy, returns the **in-plane** `ρ_sq` because quasi-2D Ewald genuinely needs the
in-plane distance. **SoEwald2D needs the 3-D distance.** Copying QuasiEwald's helper silently
changes every short-range pair's cutoff test and its `erfc(α·r)/r` argument.

This exact mistake has already been made once in this project: in Phase 1, `Ewald2D` +
`CellListQ2D` returned a wrong-signed energy (+0.0238 against a true −0.1539) because the
short-range sum trusted a neighbor list's reported `r`, which was in-plane for a 2-D finder.
Write a **new** helper for this package, do not copy, and name it so the distinction is
unmissable (suggested: `_min_image_slab`, docstring stating "returns the full 3-D squared
distance; x/y wrapped, z free").

Also: `position_check3D` returns an **all-zero sentinel triple** when no image is inside the
cutoff, so callers currently guard with `iszero(r_sq)`. Replace that with an explicit
`r_sq ≥ r_c^2` test — and then heed the §Task 2 warning about ρ = 0 below, which is the bug
that guard was accidentally hiding in QuasiEwald.

---

## File Structure

| File | Responsibility |
|---|---|
| `Project.toml` | `ExTinyMD` → `[weakdeps]` + `[extensions]`; `[sources]` for **both** ExTinyMD and QuasiEwald; version 0.2.0 |
| `src/SoEwald2D.jl` | `import ExTinyMD` removed entirely; exports updated; dispatcher functions for the two wrapper names |
| `src/types.jl` | `SoEwald2DShortPlan`, `SoEwald2DLongPlan`, `_min_image_slab`; `revise_interaction!` and `SoEwald2D_init` removed |
| `src/energy/energy_short.jl` | AoS `energy` query for the short plan |
| `src/energy/energy_long.jl` | AoS `energy` query for the long plan (kernels below `SoEwald2D_El` are already framework-free) |
| `src/force/force_short.jl` | AoS `force`/`force!` for the short plan; `Point`→`SVector` |
| `src/force/force_long.jl` | AoS `force`/`force!` for the long plan; delete the `Base.real(::Point)` piracy |
| `src/tools/direct_sum.jl`, `diff_direct_sum.jl` | AoS signatures (these are validation references) |
| `ext/SoEwald2DExTinyMDExt.jl` | NEW — wrapper structs + `ExTinyMD.energy`/`update_acceleration!` |
| `test/plan.jl` | NEW — core queries, FD self-consistency, no-mutation, ρ=0 |
| `test/adapter.jl` | NEW — wrappers driven through `simulate!` |
| `test/standalone.jl` | NEW — proves ExTinyMD is never loaded |

---

## Task 1: Baseline capture, then compat and `[sources]`

**Files:** `Project.toml`, `test/Project.toml`

- [ ] **Step 1: capture a numerical baseline BEFORE touching anything.**

Write a script to `/tmp/.../soe_baseline.jl` that records, to full precision, for at least
three parameter sets: `SoEwald2D_El`, `SoEwald2D_Es`, the force arrays from `SoEwald2D_Fl!`
and `SoEwald2D_Fs!`, and `direct_sum`/`soe_direct_sum`/`diff_direct_sum`.

**Critical:** `SimulationInfo` consumes `rand()` internally, so a baseline that constructs one
is not reproducible across code changes that alter RNG draw counts. Seed once with
`Random.seed!(...)` and **build positions/charges directly** (`rand()*L` etc.), constructing
only what you must. This is how the QuasiEwald baseline was made trustworthy; skipping it
produced hours of false "the physics changed" alarms in an earlier phase.

Also note `SoEwald2DLongInteraction`'s own default `rng = MersenneTwister(123)` and the `rbm`
(random batch) path — capture with `rbm = false` so the baseline is deterministic, and
separately confirm the `rbm = true` path still runs.

Re-run this baseline after EVERY task and report bit-identity or the exact deviation.

- [ ] **Step 2: move ExTinyMD to `[weakdeps]`, add `[extensions]`, bump versions.**

```toml
[weakdeps]
ExTinyMD = "fec76197-d59f-46dd-a0ed-76a83c21f7aa"

[extensions]
SoEwald2DExTinyMDExt = "ExTinyMD"
```

Version `0.1.5` → `0.2.0`. `[compat]`: `ExTinyMD = "0.3"`, add `QuasiEwald = "0.3"`,
`StaticArrays = "1.6"`.

- [ ] **Step 3: `[sources]` for BOTH ExTinyMD and QuasiEwald.**

Neither ExTinyMD 0.3 nor QuasiEwald 0.3 is registered yet. Use git URLs, not sibling paths — a
path pin resolves locally but breaks CI, where only this repo is checked out (verified on
ParticleMeshEwald: `expected package ExTinyMD [fec76197] to exist at path ...`).

```toml
[sources]
ExTinyMD = {url = "https://github.com/HPMolSim/ExTinyMD.jl", rev = "main"}
QuasiEwald = {url = "https://github.com/HPMolSim/QuasiEwald.jl", rev = "decouple-extinymd"}
```

`QuasiEwald` must point at the `decouple-extinymd` branch, because that work is not yet merged
to `main`. Add a comment saying so and that it must become `rev = "main"` once PR #5 lands.
Copy ParticleMeshEwald's comment block about General's automerge rejecting `[sources]`.

- [ ] **Step 4: resolve the duplicate test-environment declaration.**

This package declares its test env **twice**: `[extras]`/`[targets]` in the top-level
`Project.toml` AND a separate `test/Project.toml`. **Julia honours `test/Project.toml` when it
is present and ignores the `[extras]`/`[targets]` pair entirely** — so a `[sources]` block or
compat edit made only in the top-level file will silently do nothing.

Decide which survives and say why in the commit message. Recommended: keep `test/Project.toml`
(it is what Julia actually reads), delete the now-dead `[extras]`/`[targets]` from the
top-level file, and put the test-side `[sources]`/`[compat]` in `test/Project.toml`.
**Verify by execution** that `Pkg.test()` resolves — do not infer it.

- [ ] **Step 5:** run `Pkg.test()`, record the assertion count, re-run the baseline, commit.

---

## Task 2: `Point` → `SVector`, and the min-image replacement

**Files:** `src/types.jl`, `src/force/force_short.jl`, `src/force/force_long.jl`,
`src/energy/energy_short.jl`, `src/tools/diff_direct_sum.jl`, `src/tools/direct_sum.jl`

- [ ] **Step 1: add `_min_image_slab` to `src/types.jl`.**

Read the "THE ONE THING MOST LIKELY TO GO WRONG" section above first. The helper wraps x and y
and returns the **full 3-D** squared distance:

```julia
@inline _wrap(dx::T, L::T) where {T} = dx - L * round(dx / L)

"""
    _min_image_slab(pos_i, pos_j, L) -> (coord_i, coord_j, r_sq)

Slab-geometry nearest image: `x` and `y` wrap under `L[1]`/`L[2]`, `z` is a plain
difference (the slab axis is not periodic). Returns `pos_i` shifted to its nearest
in-plane image of `pos_j`, `pos_j` unchanged, and the **full three-dimensional**
squared distance `dx^2 + dy^2 + dz^2` between them.

The 3-D distance is the point of this helper and the reason it is not QuasiEwald's
`_min_image_q2d`, which returns the in-plane distance instead: SoEwald2D's real-space
sum is `erfc(α·r)/r` in the true separation `r`, and its cutoff test is on `r`, not on
the in-plane `ρ`. `pos_i`/`pos_j` need only support `p[1]`/`p[2]`/`p[3]` indexing.
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

- [ ] **Step 2: replace the `iszero(r_sq)` sentinel guards with `r_sq ≥ r_c^2`** at both
short-range call sites (`src/energy/energy_short.jl:18`, `src/force/force_short.jl:15`).

**A regression this WILL introduce if you are not careful.** The old sentinel had a side
effect nobody intended: it also skipped pairs at `r == 0`. An explicit `r_sq ≥ r_c^2` test
keeps them, and `SoEwald2D_Fs_pair` builds the force from a unit vector — so a `0/0` NaN
appears for any coincident pair. In QuasiEwald this exact substitution shipped a NaN that
fired for any two particles sharing a column **and** for any pair whose in-plane separation
was an exact multiple of `Lx`/`Ly` (so the minimum image wraps to exactly zero) — i.e. any
lattice initialisation.

For SoEwald2D `r == 0` means truly coincident particles, which is unphysical, but it is still
reachable by a wrapped separation. **Guard it explicitly** and add a test asserting finiteness.
Decide and justify whether the correct behaviour is a zero contribution or an error; a zero
contribution matching the `r → 0` limit is the safer default. Report which you chose.

- [ ] **Step 3: delete the `Base.real(::Point)` type piracy** (`src/force/force_long.jl:1`).

`Base.real(x::Point) = Point(real.(x.coo))` is a method on ExTinyMD's type owned by neither
package. Grep suggests every `real(...)` call in this package is on a complex scalar, not a
`Point`, so the method is likely dead — **verify that by deleting it and running the suite**
before concluding. If something does need it, the `SVector` equivalent is `real.(v)`
broadcast, which needs no new method at all.

- [ ] **Step 4:** `Point` → `SVector{3,T}` for `acceleration` and every accumulator. Leave
read-only position inputs untyped. `dist2(a,b)` → `sum(abs2, a .- b)` or the `r_sq` returned
by `_min_image_slab`.

- [ ] **Step 5:** run the suite and the baseline. `Point`→`SVector` must be bit-identical.
If it is not, **operation order** is the likely cause: QuasiEwald hit exactly this, where
`F * (d * (1/ρ))` and `(F * d) / ρ` differed in the last 1–2 ULPs. Match the original order
rather than accepting the drift. Commit.

---

## Task 3: The framework-free plans

**Files:** `src/types.jl`, `src/energy/*.jl`, `src/force/*.jl`, `src/tools/*.jl`, `test/plan.jl`

- [ ] **Step 1: define the plans.**

`SoEwald2DShortPlan` is easy — the existing `SoEwald2DShortInteraction` is already pure
parameters (`ϵ_0, L, s, α, n_atoms, r_c`). Only its supertype couples it. Copy the fields.

`SoEwald2DLongPlan` keeps the solver state (`k_set`, `soepara`, `rbm`, `rbm_p`, `P`, `prob`,
`indice`, `iterpara`, `adpara`, `parallel`, `rng`) and keeps `q`, `x`, `y`, `z` as
**plan-owned scratch** — the kernels below `SoEwald2D_El` are structure-of-arrays and already
framework-free, so a query scatters AoS `poses`/`charges` into those buffers rather than
mutating the caller's data (the ParticleMeshEwald `_scatter_long!` pattern).

**`mass` and `acceleration` leave the plan.** A framework-free solver returns forces and lets
the caller divide by mass. In QuasiEwald this was provably a bitwise no-op because every
`./ mass` sat inside a plain sum over a fixed divisor; **verify the same holds here** for
`SoEwald2D_Fl!`'s `acceleration -= Point(Fx,Fy,Fz)/mass[i]` and report the check.

- [ ] **Step 2: sign convention — state it explicitly.**

`SoEwald2D_Fl!` currently does `acceleration -= Point(Fx,Fy,Fz)/mass`, so `force_sum` returns
the energy **gradient**, not the force. The new `SoEwald2D.force`/`force!` must return the
**force** (i.e. `-gradient`), so the extension accumulates with `+=` like every other adapter
in this family. Document the flip in the docstring and make a test assert the sign against a
finite difference of the energy — a sign error here is invisible to any magnitude-only test.

- [ ] **Step 3: the AoS queries.**

```julia
SoEwald2D.energy(plan::SoEwald2DShortPlan, poses, charges; neighbor_list = nothing)
SoEwald2D.energy(plan::SoEwald2DLongPlan,  poses, charges)
SoEwald2D.force!(F, plan, poses, charges; ...)   # fills F, returns F
SoEwald2D.force(plan, poses, charges; ...)        # allocates
```

Defined but NOT exported. With no `neighbor_list`, the short plan tests every pair `O(n²)` —
it owns no cell list. Always recompute the true distance from `poses`; never trust a supplied
list's reported distance (this is the Phase 1 `Ewald2D` bug).

- [ ] **Step 4: port `direct_sum`, `soe_direct_sum`, `diff_direct_sum`** to `(plan, poses, charges)`.
These are the accuracy references the tests compare against, so they must work standalone or
the test suite stays coupled.

- [ ] **Step 5: delete `revise_interaction!` and `SoEwald2D_init`.**

`revise_interaction!` is the sys/info gather — it becomes the extension's job. `SoEwald2D_init`
(`src/types.jl:165`) references `n_atoms`, which is **not one of its parameters** — every call
raises `UndefVarError` unconditionally. Grep confirms nothing calls it (its only other mention
is the export list). Same class of dead-broken code as `QuasiEwaldRbeInit`. **Verify the grep
yourself**, then delete it and drop it from the exports.

- [ ] **Step 6: `test/plan.jl`.** Core queries against the `direct_sum` reference; `force`
vs `force!` agreement; no-mutation of `poses`/`charges`; a finite-difference check of
**every** component of both plans (F = -dE/dr) asserting the **sign**; the `r = 0` finiteness
test from Task 2.

On the FD check: if a component disagrees, **do not exclude it and move on.** In QuasiEwald an
apparent `F_z ≠ -dE/dz` was reported as a formula bug and was actually two coupled convergence
knobs — tightening the truncation parameter alone made it worse. Sweep the accuracy parameters
in two dimensions before concluding anything, and record the table.

- [ ] **Step 7:** suite + baseline + commit.

---

## Task 4: The extension

**Files:** `ext/SoEwald2DExTinyMDExt.jl` (new), `src/SoEwald2D.jl`, `test/adapter.jl` (new)

- [ ] **Step 1: wrapper structs in the extension.**

Per §4.3a a struct's supertype is fixed at definition, and per §4.3a-bis a struct definition
**cannot be dot-qualified at all**, so the wrapper types must live in the extension's
namespace. Use the dispatcher-function pattern to preserve the names:

```julia
# src/SoEwald2D.jl
export SoEwald2DShortInteraction, SoEwald2DLongInteraction
for name in (:SoEwald2DShortInteraction, :SoEwald2DLongInteraction)
    @eval function $name(args...; kwargs...)
        ext = Base.get_extension(SoEwald2D, :SoEwald2DExTinyMDExt)
        ext === nothing && error(...)   # see below
        return getfield(ext, $(QuoteNode(name)))(args...; kwargs...)
    end
end
```

**First check for type-position uses** of both names (`::Name`, `Vector{Name}`, `isa Name`,
`Name{T}` in a signature) across this repo AND `FastSpecSoG.jl`. The pattern makes the name a
**function**, so construction keeps working but every type-position use breaks. QuasiEwald had
zero such uses; this package has `SoEwald2DLongInteraction{T}` in many internal signatures, so
**expect to find some** — those internal ones become the plan type, which is fine, but any that
must stay on the wrapper mean you should rename instead, per §4.3a. Report the grep.

Make the error message distinguish "ExTinyMD is not loaded" from "ExTinyMD is loaded but the
extension failed to precompile" — the latter is what you will actually hit while developing,
and a message claiming the former wastes real debugging time.

- [ ] **Step 2: the adapter methods.** `ExTinyMD.energy` and `ExTinyMD.update_acceleration!`
for both wrappers. Gather positions in **slot** order and charge/mass via
`sys.atoms[info.particle_info[slot].id]`; accumulate with `+=`; divide by mass here.
Model on `../ExTinyMD.jl/src/interactions/electrostatics/adapter.jl`.

Restrict the neighbor-finder argument to the finder types this package actually supports
(`CellList3D`/`CellListDir3D` — check what it really uses, it is 3-D here, unlike QuasiEwald)
plus `NoNeighborFinder`, with an informative error otherwise. An untyped finder parameter
silently accepts a mismatched finder and under-counts; in QuasiEwald that was a review finding.

- [ ] **Step 3: `test/adapter.jl` — two tests, and make them able to fail.**

1. **Permutation test.** Reverse `info.particle_info` so slot order ≠ id order, give every id a
   **distinct mass and charge**, and check `ExTinyMD.energy` and `ExTinyMD.update_acceleration!`
   against direct plan queries on the same poses/charges. This is the test that catches id/slot
   and mass-division faults — in QuasiEwald it caught both (36/39 assertions failed) while the
   trajectory test caught neither.
2. **Conservation test under `simulate!`.** Put both wrappers in `sys.interactions` and run a
   trajectory. **Assert on `KE + E_elec`, not on `E_elec` alone.** A bound on electrostatic
   energy drift cannot catch a force error: in QuasiEwald, an exact 2× force fault moved it
   from 9.46e-4 to 9.41e-4 — identical to two digits — because `E_elec` is ~1% of the total and
   is set by thermal motion. `VerletProcess` with `NoThermoStat` conserves `KE + PE`, and a 2×
   force makes it conserve `KE + 2·PE`, so `KE + E_elec` drifts at first order.

**Verify both by sabotage:** introduce (a) an id/slot fault, (b) a dropped mass division,
(c) an exact 2× force error; confirm the tests fail; revert; confirm they pass. Report the
numbers for each. Set every bound from a measured value with a stated margin, not by guess.

- [ ] **Step 4:** suite + baseline + commit. The baseline may need a scratch environment with
`[sources]` path pins, since ExTinyMD is no longer resolvable from this package's main env.

---

## Task 5: Standalone proof, docs, and cleanup

**Files:** `test/standalone.jl` (new), `test/runtests.jl`, `README.md`

- [ ] **Step 1: `test/standalone.jl`.** A subprocess builds both plans from plain `SVector`
arrays, evaluates energy and force, and asserts ExTinyMD is absent from
`Base.loaded_modules` throughout.

Then make the check load-bearing: a **second** subprocess does `using ExTinyMD, SoEwald2D`
first and re-runs the identical assertion, and the test asserts that process **fails**.
Capture its stderr and assert it failed **for the intended reason** (the specific
`AssertionError`), not merely that it exited non-zero — `!success(proc)` alone passes on a
missing dependency or a precompile error, which is a test that proves nothing. This was a
review finding on QuasiEwald.

- [ ] **Step 2: README.** Add "Standalone usage (no ExTinyMD)" and "MD usage via ExTinyMD"
sections: the AoS API, `energy`/`force`/`force!` deliberately unexported, the
`r_c < min(Lx,Ly)/2` rule, the no-mutation guarantee, the force **sign** convention, the
wrapper types and that they are functions rather than types (with the
`Base.get_extension(...)` workaround), and a `[sources]`/registration note pointing at the
design doc §6b. No Documenter `(@ref)` links — the README is not rendered by Documenter.

- [ ] **Step 3:** confirm `src/` has zero ExTinyMD code references, `ExTinyMD` is in
`[weakdeps]` only, `energy`/`force`/`force!` are unexported, and every `r_c` in `test/` carries
its inline `min(Lx,Ly)/2` comparison. Full suite. Commit.

---

## Self-review checklist

1. Does `src/` mention ExTinyMD anywhere outside a comment or docstring?
2. Is `_min_image_slab` returning the **3-D** distance, and is that asserted by a test that
   would fail if someone swapped in an in-plane version?
3. Does any test pass regardless of the code under test? Every bound traceable to a measured
   number with a stated margin?
4. Was each of the three sabotage checks actually executed, with numbers?
5. Is the force **sign** asserted, not just its magnitude?
6. Baseline bit-identical after every task, or the deviation explained?
