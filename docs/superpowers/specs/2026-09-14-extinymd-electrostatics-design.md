# ExTinyMD Electrostatics Standard Library — Design

Date: 2026-09-14
Status: approved, not yet implemented

## 1. Context

`~/project/q2dmd` holds six Julia packages: `ExTinyMD.jl` (a small MD engine) and five
electrostatics packages that plug into it. A survey on 2026-09-14 found the following.

ExTinyMD is at 0.3.0 and moved to CellListMap 0.10. Every downstream package still pins
`ExTinyMD = "0.2"` and `CellListMap = "0.9"`, so none of them resolve against the current
engine.

Two incompatible interaction interface generations are in use:

| Package | pins ExTinyMD | CellListMap | entry points |
|---|---|---|---|
| ExTinyMD 0.3.0 | — | 0.10 | `energy(i, finder, sys, info)`, `update_acceleration!(i, finder, sys, info)` |
| EwaldSummations 0.2.0 | 0.2 | 0.9 | `energy(i, nb, info, atoms)`, free `force(...)` |
| FastSpecSoG 0.1.0 | 0.2.6 | — | `energy(i, nb, info, atoms)` |
| QuasiEwald 0.2.1 | 0.2 | 0.9 | current generation |
| SoEwald2D 0.1.5 | 0.2 | — | current generation |
| ParticleMeshEwald 0.1.0 | weakdep | 0.9.14 | `energy(pme, x, y, z, q)` only |

Specific defects confirmed by reading the source:

- EwaldSummations and FastSpecSoG define `energy` and a free-function `force`, but never
  `ExTinyMD.update_acceleration!` — the method `simulate!` actually calls. Neither package
  can drive an MD run.
- ParticleMeshEwald declares `[weakdeps] ExTinyMD, EwaldSummations` but ships no `ext/`
  directory. Those extensions do not exist. It is also energy-only.
- ExTinyMD 0.3 is internally inconsistent: `AllNeighborFinder` and `NoNeighborFinder` name
  their field `neighborlist` (`src/types.jl:86,97`), while every interaction reads
  `.neighbor_list` (`src/interactions/lennard_jones.jl:15,36`). `LennardJones` +
  `AllNeighborFinder` throws `FieldError`. The `CellList*` finders use `neighbor_list`, which
  is why the shipped tests (25/25 passing, Julia 1.13) do not catch it.
- Using any of these packages standalone requires constructing ExTinyMD types
  (`SimulationInfo`, `Atom`, `MDSys`, `Boundary`, a cell-list finder) even when the caller
  only wants an energy from arrays.

## 2. Goals

1. Implement Ewald3D, Ewald2D, ICM+Ewald2D, ICM+Ewald3D+ELC, PME3D, and ICM+PME3D+ELC in
   ExTinyMD's standard library, with both energies and forces.
2. Give every method a framework-free array API: a plan object plus AoS positions and
   charges. No ExTinyMD type is ever constructed to evaluate an energy or a force.
3. Make the five downstream packages resolve and run against current ExTinyMD, and decouple
   their numerical cores from it.
4. Document ExTinyMD with a Documenter site, including a method-selection guide.

## 3. Non-goals

- Mesh (NUFFT) Ewald2D and ICM + mesh-Ewald2D. The 2D long-range kernel carries
  `exp(±kz) erfc(k/2α ± αz)`, which does not factor in k-space; a mesh version needs a
  z-quadrature or Chebyshev expansion on top of a 2D NUFFT — essentially the FastSpecSoG /
  SoEwald2D approach. Explicitly deferred to its own project.
- Migrating ExTinyMD's internal `Point{N,T}` to `SVector`. Large refactor, not required.
- GPU execution. `ParticleMeshEwald`'s KernelAbstractions usage is CPU-only in practice and
  is not being extended.
- Re-deriving ICM force conventions. See §7.

## 4. Architecture: the two-layer contract

Every method in every package is split into a framework-free numerical core and a thin MD
adapter.

### 4.1 Core layer

```julia
plan = Ewald3D(n_atoms, L; α = 0.2, s = 4.0, ϵ = 1.0, ϵ_inf = Inf)

E = energy(plan, poses, charges)             # ::T
F = force(plan, poses, charges)              # ::Vector{SVector{3,T}}
force!(F, plan, poses, charges)              # in-place
```

- `poses::AbstractVector` of any 3-component indexable value. `Vector{SVector{3,T}}` is
  canonical and documented; kernels access `p[1], p[2], p[3]` only, so `Vector{Point{3,T}}`
  and `Vector{NTuple{3,T}}` work with no conversion layer and no method explosion.
- `charges::AbstractVector{T}`.
- `L::NTuple{3,T}` is passed to the constructor, not the query. Geometry is part of the plan.
- The plan owns **all** scratch: k-set, cell list, NUFFT plan, per-thread force
  accumulators, and (for PME) scaled-coordinate buffers. Query functions allocate nothing
  after construction, except `force` which allocates its return value; `force!` does not.
- Both queries accept `neighbor_list = nothing` as a keyword. When a list is supplied it is
  used as-is; when `nothing`, the plan rebuilds from its own cell list. This is what lets the
  MD adapter hand over the list the simulation loop already maintains instead of
  duplicating neighbor search every step.

`SVector{3,T}` is the interchange type across package boundaries — StaticArrays is a neutral
dependency already used by ExTinyMD, whereas `Point` is ExTinyMD-specific and would
reintroduce the coupling being removed.

### 4.2 Adapter layer

```julia
ExTinyMD.energy(interaction, finder, sys, info)
ExTinyMD.update_acceleration!(interaction, finder, sys, info)
```

The adapter extracts charges and positions from `sys.atoms` / `info.particle_info`, forwards
`finder.neighbor_list`, calls the core, and for forces divides by mass and accumulates into
`info.particle_info[i].acceleration`.

Charge and position extraction must respect ExTinyMD's id indirection. The established
pattern, which this design keeps, is:

```julia
charge   = [atoms[info.particle_info[i].id].charge for i in 1:n_atoms]
position = [info.particle_info[i].position         for i in 1:n_atoms]
```

The `atoms` vector is indexed by particle **id**, while `particle_info` is indexed by current
storage slot; `info.id_dict` maps id to slot. Positions are read in slot order and charges
gathered to match, so core-layer index `i` means "slot `i`" consistently. These two
allocations per call are hoisted into the plan as reusable buffers.

`ExTinyMD.energy` for the current generation takes `(interaction, finder, sys, info)`. That
signature is canonical; the older `(interaction, neighbor, info, atoms)` form used by
EwaldSummations and FastSpecSoG is dropped.

### 4.3 Placement

| | core | adapter | ExTinyMD dep |
|---|---|---|---|
| ExTinyMD stdlib | `src/interactions/electrostatics/` | `.../adapter.jl` | n/a |
| Other five packages | `src/` | `ext/<Pkg>ExTinyMDExt.jl` | `[weakdeps]` |

ExTinyMD's own stdlib necessarily lives inside ExTinyMD, so `using ExTinyMD` is required to
reach it. The constraint being satisfied is that no ExTinyMD *type* need be *constructed*:
`energy(Ewald3D(n, L; α, s), poses, charges)` is the whole standalone path.

## 5. Stdlib structure

```
src/interactions/electrostatics/
  common.jl        # position/charge buffers, neutrality check, k-set generation, prefactors
  short.jl         # EwaldShort
  long_ewald3d.jl  # Ewald3DLong
  long_ewald2d.jl  # Ewald2DLong
  icm.jl           # ICM{Inner}, ICM_reflect, ELC slab term
  ewald.jl         # EwaldInteraction{S,L} composite, convenience constructors
  adapter.jl       # ExTinyMD.energy, ExTinyMD.update_acceleration!
ext/ExTinyMDFINUFFTExt.jl   # PME3DLong
```

### 5.1 `common.jl`

- `k_set_3D(k_c, L) -> Vector{NTuple{4,T}}` — `(kx, ky, kz, k)` for `0 < k ≤ k_c` over a
  sphere. Ported from `EwaldSummations/src/Ewald3D/Ewald3D.jl`.
- `k_set_2D(k_c, L) -> Vector{NTuple{3,T}}` — `(kx, ky, k)`, spherical cutoff in the xy
  plane. Ported from `EwaldSummations/src/Ewald2D/Ewald2D.jl`.
- `check_neutrality(charges; atol)` — warns once when `|Σq| > atol`. Ewald3D's conditionally
  convergent sum assumes neutrality; a non-neutral system silently produces a
  volume-dependent offset. None of the current packages check this.
- Parameter conversion `r_c = s/α`, `k_c = 2αs` — the existing convention across all
  packages, retained so parameters transfer between them.

### 5.2 `short.jl` — `EwaldShort`

One real-space kernel, used by every method.

```
E_s = 1/(4πϵ) [ Σ_{i<j, r<r_c} q_i q_j erfc(α r_ij)/r_ij  −  (α/√π) Σ_i q_i² ]
```

Fields: `α`, `r_c`, `ϵ`, `n_atoms`, a boundary convention tag, and a cell list. The
convention tag selects minimum-image handling: `Periodic3D` uses `position_check3D`
semantics, `PeriodicQ2D` uses `position_checkQ2D` (periodic in x,y only). This is the single
axis on which Ewald3D and Ewald2D differ in their short-range part; everything else is
shared.

Force, analytic:

```
dE/dr = −q_i q_j [ erfc(αr)/r² + (2α/√π) e^{−α²r²}/r ]
F_ij  = −dE/dr · (r_i − r_j)/r
```

This replaces the per-pair `ForwardDiff.derivative` call in
`EwaldSummations/src/Ewald3D/Ewald3D_short.jl:Ewald3D_Fs_pair` and its Ewald2D twin. One
closed-form expression, no dual numbers in the inner loop, and ForwardDiff leaves the
dependency list.

Threading: pairs partitioned across tasks with per-task force buffers reduced at the end.
Not `Threads.threadid()` indexing — that is unsound under task migration and is the cause of
the `output[Threads.threadid() - 1]` indexing bug in
`ParticleMeshEwald/src/energy.jl:energy_short_kernel!`, where `threadid() == 1` indexes
element 0.

### 5.3 `long_ewald3d.jl` — `Ewald3DLong`

```
E_l  = 1/(2Vϵ) Σ_{0<k≤k_c} |ρ_k|² e^{−k²/4α²}/k²         ρ_k = Σ_j q_j e^{i k·r_j}
E_k0 = |P|² / (2Vϵ(2ϵ_inf+1))                            P   = Σ_j q_j r_j
```

with `V = Lx Ly Lz`, k-set covering both `±k`. `ϵ_inf = Inf` (conducting boundary) zeroes the
dipole term; finite `ϵ_inf` gives the surface correction. Ported from
`EwaldSummations/src/Ewald3D/Ewald3D_long.jl`, which uses the equivalent form
`2π/V · |ρ_k|² e^{−k²/4α²}/k²` later divided by `4πϵ`.

Forces are analytic and already correct in the reference
(`Ewald3D_long_force_k!`, `Ewald3D_long_force_k0!`); they are ported with the structure
factor hoisted out of the per-particle loop.

One correctness note carried over: the reference accumulates `ρ_k` in `ComplexF64`
unconditionally (`zero(ComplexF64)`, `1.0im`), which silently demotes `Float32` and blocks
higher-precision `T`. The port uses `Complex{T}`.

### 5.4 `long_ewald2d.jl` — `Ewald2DLong`

Periodic in x and y, free in z. Per
`EwaldSummations/src/Ewald2D/Ewald2D_long.jl`, for each `k` in the 2D k-set:

```
E_k  = 1/ϵ Σ_i Σ_j q_i q_j cos(k·ρ_ij)
         [ e^{k z_ij} erfc(k/2α + α z_ij) + e^{−k z_ij} erfc(k/2α − α z_ij) ] / (8 Lx Ly k)

E_k0 = −1/ϵ Σ_i Σ_j q_i q_j [ e^{−(α z_ij)²}/(α√π) + z_ij erf(α z_ij) ] / (4 Lx Ly)
```

Note the long-range prefactor is `1/ϵ`, not `1/(4πϵ)` — the `4π` is absorbed into the
expressions above. The short-range part of the same method uses `1/(4πϵ)`. This asymmetry is
inherited from the reference and is deliberate; it is called out in the docstrings because it
looks like a bug and is not.

This is `O(N²K)`, inherently — it is the exact 2D Ewald sum, kept as the accuracy reference
for quasi-2D systems, not as a production method for large N.

Two structural changes from the reference:

- `Ewald2D_long_force_k0!` and `Ewald2D_long_force_k!` open a `@threads` region *inside* the
  per-`i` loop, so a thread region is launched once per particle (and in `force_k!`, once per
  particle per wavevector) with `Atomic{T}` accumulation inside. The port threads the outer
  loop over `i` with plain per-task accumulators.
- Guard `exp(±k z_ij)` against overflow. `ICM_Ewald2D_long_force_k!` already does this with
  an `abs(k*z_ij) > 650` test, but the non-ICM `Ewald2D_long_force_k!` does not, and quasi-2D
  systems with large `L_z` can reach it. The paired `erfc` factor underflows to zero at the
  same time, so the product is zero; compute it as zero rather than `Inf * 0 = NaN`.

### 5.5 `icm.jl` — `ICM{Inner}` and ELC

ICM is a wrapper over any inner interaction, not a per-method variant. This is what makes
ICM+Ewald2D, ICM+Ewald3D+ELC and ICM+PME3D+ELC one implementation instead of three.

`ICM_reflect(γ, L, N_image, poses, charges) -> (ref_poses, ref_charges)` builds the image
series, real particles first, then images interleaved up/down. Recurrence from
`EwaldSummations/src/direct_sum/direct_sum_Q2D.jl`:

```
z_up[1]   = 2L_z − z          q_up[1]   = γ_up   q
z_down[1] = −z                q_down[1] = γ_down q
z_up[m]   = 2L_z − z_down[m−1]    q_up[m]   = γ_up   q_down[m−1]
z_down[m] = −z_up[m−1]            q_down[m] = γ_down q_up[m−1]
```

`γ = (γ_up, γ_down)` are the dielectric contrast ratios at the two confining walls,
conventionally `γ = (ϵ_mid − ϵ_out)/(ϵ_mid + ϵ_out)` for each wall, with `γ = 0` for a
dielectrically matched wall and `γ → ±1` for the perfect-conductor / perfect-insulator
limits. The API takes `γ` directly as a number pair rather than deriving it from
permittivities, matching the existing packages. Image count `N_image` per side.
The port preallocates into plan-owned buffers of length `n_atoms(2 N_image + 1)` rather than
`push!`-ing into fresh vectors on every call.

Semantics preserved verbatim from `ICM_short.jl` / `ICM_long.jl`:

- Pairs are partitioned into real–real and real–image. `CellListICM` builds both lists with a
  unitcell of `(Lx, Ly, (2 N_image + 1) L_z + 2 r_c)`.
- Real–image **energies** carry a factor ½; real–real do not.
- Self-energy is summed over real particles only.
- Forces accumulate on the real index only: real–real pairs get equal and opposite
  contributions, real–image pairs contribute to the real particle alone.
- Long-range `ρ_k` is summed over **all** reflected charges, while the force loop runs over
  real indices only.

ELC slab term, for the 3D variant, from `ICM_Ewald3D_long_energy_slab`:

```
E_slab = −1/(4πϵ) · π/(Lx Ly (2 N_pad + 1) L_z) Σ_i q_i Σ_j q_j (z_i − z_j)²
```

with the matching analytic z-force. `N_pad` sets the z-padding factor; the inner Ewald3D's
k-set is generated for `(Lx, Ly, (2 N_pad + 1) L_z)`.

### 5.6 `ewald.jl` — composition

`EwaldInteraction{S,L} <: AbstractInteraction` holds a short and a long part. `energy` is the
sum of both; `force!` calls both accumulating into one buffer.

A single composite interaction rather than two `sys.interactions` entries (the pattern
QuasiEwald and SoEwald2D use). It receives one finder — the cell list, which the short part
needs and the long part ignores. It also avoids widening `MDSys`'s `T_INTERACTION`
parameter, which is a concrete tuple type; heterogeneous interaction vectors already degrade
to abstract element types there and adding more entries makes it worse.

Convenience constructors, one per requested method:

```julia
Ewald3D(n_atoms, L; α, s, ϵ=1, ϵ_inf=Inf)
Ewald2D(n_atoms, L; α, s, ϵ=1)
ICMEwald2D(n_atoms, L; α, s, γ, N_image, ϵ=1)
ICMEwald3D(n_atoms, L; α, s, γ, N_image, N_pad, ϵ=1)   # ICM + Ewald3D + ELC
PME3D(n_atoms, L; α, s, ϵ=1, ϵ_inf=Inf)                # requires `using FINUFFT`
ICMPME3D(n_atoms, L; α, s, γ, N_image, N_pad, ϵ=1)     # requires `using FINUFFT`
```

### 5.7 `ext/ExTinyMDFINUFFTExt.jl` — `PME3DLong`

Loaded only when the user does `using FINUFFT`, keeping the binary artifact off
LJ-only users. Replaces `Ewald3DLong`'s direct k-sum inside the same `EwaldInteraction`
composite, so the short-range part, the adapter, and ICM are all shared unchanged.

Energy: type-1 NUFFT (nonuniform → uniform) maps charges at particle positions to `ρ_k` on
the `(2n_k+1)³` grid, then `E_l = 1/(2Vϵ) Σ |ρ_k|² D_k` with
`D_k = e^{−k²/4α²}/k²`, `D_0 = 0`. This matches `ParticleMeshEwald/src/energy.jl` and is
algebraically identical to §5.3's normalization.

Forces: type-2 NUFFT (uniform → nonuniform) evaluates the k-space field back at particle
positions. Three transforms give the three force components, or one vectorized transform with
`ik D_k ρ_k` as input. **This does not exist in ParticleMeshEwald today** and is new work.

Three defects in the existing PME code that the port fixes:

1. `energy_long` scales the caller's coordinate arrays in place (`x .*= 2π/L[1]`) and divides
   them back afterwards. This mutates caller data and leaves it corrupted if the NUFFT
   throws. The port scatters into plan-owned scaled-coordinate buffers.
2. `energy_short_kernel!` accumulates into `output[Threads.threadid() - 1]`, which indexes
   element 0 when `threadid() == 1`, and relies on `threadid()` for correctness under a
   migrating task scheduler. Replaced by the task-partitioned reduction in §5.2.
3. `PME`'s k-space truncation is a **rectangular box** `|m_i| ≤ n_k[i]`, whereas
   `Ewald3DLong`'s is a **sphere** `|k| ≤ k_c`. These do not agree at finite cutoff. See
   §8.2 for how the cross-validation test handles this.

## 6. Data flow

Standalone:

```
poses, charges ──► plan (owns k-set, cell list, NUFFT plan, buffers)
                     ├─ short: neighbor list ──► erfc pair sum + self term
                     └─ long:  ρ_k (direct sum or NUFFT) ──► k-space sum
                                                    └──► E or F::Vector{SVector{3,T}}
```

Inside MD, per `Verlet.info_update!` step:

```
simulate! ──► info_update! ──► update_acceleration!(interaction, finder, sys, info)
                                    │
                                    ├─ update_finder!(finder, info)       # existing cell list
                                    ├─ gather charges/positions into plan buffers
                                    ├─ core force!(F, plan, poses, charges;
                                    │              neighbor_list = finder.neighbor_list)
                                    └─ info.particle_info[i].acceleration += F[i]/mass[i]
```

With ICM inserted, the gather step is followed by `ICM_reflect` into the extended buffers,
the inner method runs on those, and only the first `n_atoms` entries of the force buffer are
written back.

## 7. Risk: ICM force conventions

Image positions depend on the real particles' `z` coordinates, so differentiating the ICM
energy involves a chain rule through the reflection map. EwaldSummations' force routines
encode a specific convention — forces accumulate on the real index only, real–image pairs
contribute once — which is the "frozen image" treatment.

This design **ports that convention exactly and pins it with tests**. It does not re-derive
it. The reasoning: the convention is what the user's existing published results were computed
with, so a silent change of physics is worse than a preserved quirk, and the
finite-difference test in §8.3 will report disagreement between force and energy gradient if
the convention is internally inconsistent — turning an unknown into a measured number rather
than a guess.

If the finite-difference check does fail for ICM specifically, that is a finding to report,
not a licence to change the formula. Stop and surface it.

## 8. Testing

Test files under `test/electrostatics/`, included from `test/runtests.jl`.
`EwaldSummations` is added to `[extras]` and the `test` target as the independent oracle.

### 8.1 Unit

- `k_set_3D` / `k_set_2D`: count, spherical cutoff respected, `±k` symmetry, no `k=0`.
- `ICM_reflect`: image count `n(2N+1)`, recurrence values against hand-computed small cases,
  charge signs for `γ = (0,0)`, `(1,1)`, `(−1,−1)`.
- `EwaldShort`: against a brute-force `O(N²)` real-space sum over periodic images, both
  boundary conventions.
- `check_neutrality`: warns for `Σq ≠ 0`, silent otherwise.

### 8.2 Cross-validation

- Ewald3D vs direct lattice summation (`EwaldSummations`' `Energy_3D`) at small N.
- Ewald2D vs quasi-2D direct summation (`Energy_Q2D`).
- ICM+Ewald2D vs ICM + direct summation.
- ICM+Ewald3D+ELC vs ICM+Ewald2D — different algorithms, same physics; agreement here is the
  strongest available check.
- Ewald3D vs PME3D. **Both cutoffs must be pushed until truncation error is below the
  comparison tolerance**, because the two use different k-space truncation shapes (sphere vs
  box, §5.7). Comparing them at matched `s` and loose tolerance would pass while hiding a
  real error; comparing at converged cutoffs against the direct sum is the meaningful test.
  The test asserts both converge to the direct-sum value, not that they agree with each other
  at finite cutoff.

### 8.3 Energy–force consistency

For every method: `F_i ≈ −∂E/∂r_i` by central finite differences on a small random
configuration, componentwise, with step size chosen for the `Float64` sweet spot and a
relative tolerance stated per method.

**No package in the repository currently tests this.** It is the only test that catches a
wrong force, and every force in this design is either newly written (PME3D), newly
analytic (short-range), or restructured for threading (Ewald2D long). This tier is not
optional.

### 8.4 Integration

Short microcanonical Verlet run through `simulate!` for each method, asserting bounded total
energy drift. Catches adapter bugs — mass division, id/slot mix-ups, accumulation into the
wrong particle, stale neighbor lists — that per-call unit tests cannot.

### 8.5 Regression

A test for the `neighborlist` / `neighbor_list` field inconsistency of §1:
`LennardJones` + `AllNeighborFinder` must run. This is the bug the current suite misses.

## 9. Documentation

Documenter.jl site in `ExTinyMD.jl/docs/`, deployed by CI.

- Getting started: install, a complete LJ fluid example, a complete charged-system example.
- MD core reference: `Boundary`, `MDSys`, `SimulationInfo`, simulators, thermostats, loggers,
  neighbor finders.
- Interactions: Lennard-Jones, substrate LJ, external field.
- **Electrostatics guide**, the substantive new prose:
  - Method selection table: boundary condition (3D periodic / quasi-2D) × dielectric setup
    (none / confined) × system size → recommended method.
  - Parameter selection: the `r_c = s/α`, `k_c = 2αs` convention, how `s` sets accuracy, cost
    scaling of each method.
  - Error estimates, from `EwaldSummations/src/error_estimate.jl`
    (`icm_energy_error`, `elc_energy_error`, `icm_elc_energy_error`), citing arXiv:2503.18126.
  - The `1/ϵ` vs `1/(4πϵ)` prefactor asymmetry of §5.4, stated explicitly.
- Two-layer API reference: standalone and in-MD usage for each method.
- Docstrings on all public types and functions. ExTinyMD currently has essentially none;
  EwaldSummations has a handful.

## 10. Phasing

Each phase gets its own implementation plan. **The plan following this spec covers Phase 1
only** — the three phases are too large for one plan, and Phases 2 and 3 both depend on
interfaces that Phase 1 establishes, so planning them in detail now would be planning
against assumptions rather than code.

**Phase 1 — core and direct methods.** `neighbor_list` field fix and its regression test;
`common.jl`; `EwaldShort`; `Ewald3DLong`; `Ewald2DLong`; `EwaldInteraction`; `ICM` + ELC;
adapter; test tiers 8.1–8.4 for these; docs skeleton with the site building in CI.

**Phase 2 — particle mesh.** `ext/ExTinyMDFINUFFTExt.jl` with PME3D energy and forces;
`ICMPME3D`; the §8.2 convergence test; docs pages for both.

**Phase 3 — downstream decoupling.** Per package: add the AoS core API, move ExTinyMD to
`[weakdeps]` with the adapter in `ext/`, bump `ExTinyMD` and `CellListMap` compat, get the
existing suite green.

- `ParticleMeshEwald`: AoS API over the current SoA one; write the missing `ext/`.
- `EwaldSummations`: thin to direct summation, `ICM_reflect`, and error estimates — its
  k-space Ewald implementations are superseded by the stdlib. Drop ForwardDiff and the
  `@reexport using ExTinyMD`.
- `QuasiEwald`, `SoEwald2D`: already on the current adapter generation; needs the AoS core
  extracted and compat bumps.
- `FastSpecSoG`: also needs `update_acceleration!` written, which it has never had.

Phase 3 is the largest and least certain phase — it is effectively a public-API change to
four packages. It is sequenced last deliberately, so it targets a finished ExTinyMD
interface rather than a moving one.

## 11. Decisions on record

| Question | Decision |
|---|---|
| Sequencing | Stdlib first, then downstream compat |
| FINUFFT dependency | Package extension (`weakdeps`), not a hard dep |
| EwaldSummations' fate | Thin to reference/validation; stdlib owns the k-space methods |
| Particle-mesh scope | PME3D and ICM+PME3D+ELC only |
| Mesh Ewald2D | Out of scope |
| Canonical signature | `(interaction, finder, sys, info)` |
| Interchange type | `SVector{3,T}` |
| ICM force convention | Ported verbatim, pinned by test, not re-derived |

## 12. Success criteria

1. All six methods callable standalone from arrays with no ExTinyMD type constructed.
2. All six callable inside `simulate!` via `update_acceleration!`.
3. Every method agrees with direct summation to a stated tolerance at small N.
4. Every method passes the finite-difference force check.
5. ExTinyMD's suite passes, including the `AllNeighborFinder` regression test.
6. All five downstream packages resolve against current ExTinyMD and their suites pass.
7. Documenter site builds in CI with the electrostatics guide and full docstring coverage of
   the new public API.
