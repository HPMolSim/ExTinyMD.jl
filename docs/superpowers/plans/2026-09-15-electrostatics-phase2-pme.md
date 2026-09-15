# ExTinyMD Electrostatics — Phase 2: Particle Mesh Ewald

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `PME3D` and `ICMPME3D` to ExTinyMD as a FINUFFT-backed package extension, reducing the `O(N·K)` reciprocal-space sum to `O(N log N)` while remaining bit-comparable to the direct `Ewald3DLong`.

**Architecture:** A new long-range solver, `PME3DLong`, that plugs into the existing `EwaldInteraction` composite. The real-space kernel (`EwaldShort`), the MD adapter, and `ICM` are all reused unchanged. FINUFFT enters as a weak dependency so the binary artifact stays off users who never ask for it.

**Tech Stack:** Julia 1.10+, FINUFFT.jl 3.x (weak dep), StaticArrays, SpecialFunctions.

**Spec:** `docs/superpowers/specs/2026-09-14-extinymd-electrostatics-design.md` §5.7, as amended by commit `93135b8`.

**Baseline:** Phase 1 merged to `main` (`e7aa93b`), 816 tests passing, docs live.

## Global Constraints

- `julia = "1.10"`. No newer syntax.
- FINUFFT is a **weak** dependency. It goes in `[weakdeps]` and `[extensions]`, **never** `[deps]`. `Pkg.test()` without FINUFFT installed must still pass; the extension's tests are skipped or run in the test environment where FINUFFT *is* available (see Task 1).
- **The k-space mask is spherical.** `D_k = exp(−k²/4α²)/k²` for `0 < |k| ≤ k_c`, and **exactly zero** everywhere else on the rectangular FINUFFT grid, including `k = 0`. This is what makes `PME3DLong` and `Ewald3DLong` the same sum, and every correctness test in this plan depends on it. Do not substitute a box cutoff.
- **Zero every k-space buffer before filling it.** Only masked-in grid points are written, so `similar()` leaves NaN garbage in the rest, which a type-2 transform then spreads to every particle. This is not hypothetical — it bit the controller while validating this design, and the symptom was two force components silently becoming `NaN` while the third was exact.
- Never mutate the caller's position arrays. Scale into plan-owned buffers. (`ParticleMeshEwald`'s `energy_long` does `x .*= 2π/L[1]` and divides back, which corrupts caller data if the transform throws.)
- Never index a thread-local accumulator by `Threads.threadid()`. Add no threading this plan does not ask for.
- Component `*_force!` accumulate into their buffer; the composite `coulomb_force!` zeroes first. `PME3DLong.long_force!` therefore **accumulates**.
- `n_target` semantics are identical to `Ewald3DLong`: targets are `1:n_target`, sources are all `1:n_atoms`, default makes them coincide. `ICMPME3D` is the only caller passing a smaller value.
- `r_c = s/α` strictly less than half the smallest periodic box side.
- Work on branch `pme-finufft`. Commit after every task.
- Suite: `julia --project=. -e 'using Pkg; Pkg.test()'` from the package root.

## Measured reference values

The controller validated this entire design in a standalone script before writing this plan. Quote these when your tests need a target:

| quantity | measured agreement with `Ewald3DLong` |
|---|---|
| long-range energy, spherical mask, `n=40, L=12, α=0.8, s=4` | **6.4e-16** (identical 7688-vector k-sets) |
| energy with `n_target` = 30, 20, 10, 1 | 2.8e-16, 7.9e-16, 1.9e-16, 4.2e-16 |
| long-range **force**, via three type-2 transforms | **5.9e-15** relative, 2.6e-17 absolute |
| reusing one type-2 plan for three `exec` calls | exact, 0.0 difference |

So the tests here assert agreement at `rtol` around `1e-12`, not a convergence study.

## File Structure

| File | Responsibility |
|---|---|
| `src/interactions/electrostatics/long_pme3d.jl` | `PME3DLong` **struct** (no FINUFFT reference), plus `PME3D`/`ICMPME3D` stubs and their "load FINUFFT" fallbacks |
| `ext/ExTinyMDFINUFFTExt.jl` | the constructor and the `long_energy`/`long_force!` methods — everything that touches FINUFFT |
| `test/electrostatics/test_pme3d.jl` | all Phase 2 tests |
| `Project.toml` | `[weakdeps]`, `[extensions]`, `[compat]`, and FINUFFT in the test target |
| `docs/src/electrostatics.md` | method table, extension loading, PME guidance |

**Why the struct lives in the parent and its methods in the extension:** a Julia extension cannot add exported names to its parent module. Putting the struct in `src/` with free type parameters for the two FINUFFT plans lets `ExTinyMD` export `PME3DLong`, `PME3D` and `ICMPME3D` normally, and lets a user who has not loaded FINUFFT get a clear error instead of `UndefVarError`.

---

### Task 1: `PME3DLong` struct, weak dependency, and energy

**Files:**
- Create: `src/interactions/electrostatics/long_pme3d.jl`
- Create: `ext/ExTinyMDFINUFFTExt.jl`
- Create: `test/electrostatics/test_pme3d.jl`
- Modify: `Project.toml`, `src/ExTinyMD.jl`, `test/runtests.jl`

**Interfaces:**
- Consumes: `k_set_3D` is *not* used here — the grid is rectangular and masked. Consumes `ewald_cutoffs` from `common.jl`.
- Produces:
  - `PME3DLong{T,P1,P2}` struct, fields listed below
  - `PME3DLong(n_atoms, L; α, s, ϵ=one(T), ϵ_inf=T(Inf))` (in the extension)
  - `ExTinyMD.long_energy(long::PME3DLong{T}, poses, charges; n_target = long.n_atoms)::T`
  - `PME3D`, `ICMPME3D` stub functions, exported, with informative fallbacks

- [ ] **Step 1: Declare the weak dependency**

In `Project.toml`, add (uuid is FINUFFT.jl's, verify with `julia -e 'using Pkg; Pkg.add("FINUFFT")'` then `Pkg.status`):

```toml
[weakdeps]
FINUFFT = "d8beea63-0952-562e-9c6a-8e8ef7364055"

[extensions]
ExTinyMDFINUFFTExt = "FINUFFT"
```

and to `[compat]`:

```toml
FINUFFT = "3"
```

and add FINUFFT to `[extras]` and the `test` target, so the extension is actually exercised by `Pkg.test()`:

```toml
[extras]
Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
FINUFFT = "d8beea63-0952-562e-9c6a-8e8ef7364055"

[targets]
test = ["Test", "FINUFFT"]
```

Then `julia --project=. -e 'using Pkg; Pkg.resolve()'`.

Note this is the one place FINUFFT is allowed outside `[weakdeps]` — a test-target entry does not make it a runtime dependency of the package.

- [ ] **Step 2: Write the failing test**

Create `test/electrostatics/test_pme3d.jl`:

```julia
using FINUFFT

@testset "PME3DLong energy equals Ewald3DLong exactly" begin
    # With the spherical D_k mask the two solvers sum over an identical k-set, so
    # this is an exact comparison rather than a convergence study. Controller
    # measured 6.4e-16 on this configuration.
    Random.seed!(4242)
    n = 40
    L = (12.0, 12.0, 12.0)
    α, s = 0.8, 4.0          # r_c = 5.0 < L/2 = 6
    poses = [SVector(rand() * L[1], rand() * L[2], rand() * L[3]) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]

    ewald = Ewald3DLong(n, L; α = α, s = s)
    pme   = PME3DLong(n, L; α = α, s = s)

    @test isapprox(long_energy(pme, poses, charges),
                   long_energy(ewald, poses, charges), rtol = 1e-12)
end

@testset "PME3DLong honours n_target" begin
    # ICMPME3D needs targets over real particles and sources over all reflected
    # charges. Controller measured 2e-16..8e-16 across these values.
    Random.seed!(99)
    n = 30
    L = (12.0, 12.0, 12.0)
    α, s = 0.8, 3.5
    poses = [SVector(rand() * L[1], rand() * L[2], rand() * L[3]) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]

    ewald = Ewald3DLong(n, L; α = α, s = s)
    pme   = PME3DLong(n, L; α = α, s = s)

    for nt in (n, 20, 10, 1)
        @test isapprox(long_energy(pme, poses, charges; n_target = nt),
                       long_energy(ewald, poses, charges; n_target = nt),
                       rtol = 1e-12)
    end
end

@testset "PME3DLong k-space buffers are zeroed, not merely allocated" begin
    # Only masked-in grid points are written each call, so a buffer left
    # uninitialised (or dirty from a previous call with different positions)
    # leaks into the transform. Calling twice with different configurations must
    # not let the first contaminate the second.
    Random.seed!(31337)
    n = 20
    L = (12.0, 12.0, 12.0)
    α, s = 0.8, 3.5
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]
    A = [SVector(rand() * L[1], rand() * L[2], rand() * L[3]) for _ in 1:n]
    B = [SVector(rand() * L[1], rand() * L[2], rand() * L[3]) for _ in 1:n]

    pme = PME3DLong(n, L; α = α, s = s)
    E_B_first = long_energy(pme, B, charges)      # fresh plan
    long_energy(pme, A, charges)                  # dirty it
    E_B_after = long_energy(pme, B, charges)      # must be unchanged
    @test isapprox(E_B_first, E_B_after, rtol = 1e-14)
    @test isfinite(E_B_after)
end

@testset "PME3D constructors error helpfully without FINUFFT" begin
    # Not testable in this file — FINUFFT is loaded. See Task 1 Step 6, which
    # checks the fallback message by other means.
    @test isdefined(ExTinyMD, :PME3D)
    @test isdefined(ExTinyMD, :ICMPME3D)
end
```

Add to `test/runtests.jl`, inside the existing `@testset "electrostatics"` block:

```julia
    include("electrostatics/test_pme3d.jl")
```

- [ ] **Step 3: Run to verify it fails**

```
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: FAIL with `UndefVarError: PME3DLong not defined`.

- [ ] **Step 4: Write the struct and stubs in the parent**

Create `src/interactions/electrostatics/long_pme3d.jl`:

```julia
"""
    PME3DLong(n_atoms, L; α, s, ϵ = 1.0, ϵ_inf = Inf)

Particle-mesh reciprocal-space solver for the triply periodic Ewald sum, computing
the same quantity as [`Ewald3DLong`](@ref) in `O(N log N)` instead of `O(N·K)`.

The structure factor is evaluated with a type-1 NUFFT onto a rectangular grid, and
the spherical cutoff `0 < |k| ≤ k_c` is recovered by zeroing the Green's function
`D_k` outside it. The two solvers therefore sum over an **identical** k-set and
agree to machine precision, rather than merely converging to the same limit.

Requires `using FINUFFT`; the implementation lives in a package extension so the
binary artifact stays off users who never ask for it.
"""
struct PME3DLong{T, P1, P2}
    α::T
    r_c::T
    k_c::T
    ϵ::T
    ϵ_inf::T
    L::NTuple{3,T}
    n_atoms::Int
    n_k::NTuple{3,Int}
    # D is the spherical-masked Green's function: exp(-k²/4α²)/k² inside the
    # cutoff, exactly zero outside it and at k = 0.
    D::Array{T,3}
    ρ_src::Array{Complex{T},3}
    ρ_tgt::Array{Complex{T},3}
    hx::Array{Complex{T},3}
    hy::Array{Complex{T},3}
    hz::Array{Complex{T},3}
    xs::Vector{T}
    ys::Vector{T}
    zs::Vector{T}
    ox::Vector{Complex{T}}
    oy::Vector{Complex{T}}
    oz::Vector{Complex{T}}
    plan1::P1
    plan2::P2
end

Base.show(io::IO, l::PME3DLong) =
    print(io, "PME3DLong(α = $(l.α), k_c = $(l.k_c), ϵ = $(l.ϵ), ϵ_inf = $(l.ϵ_inf), " *
              "grid = $(size(l.D)))")

"""
    PME3D(n_atoms, L; α, s, ϵ = 1.0, ϵ_inf = Inf)

Particle-mesh Ewald summation for a triply periodic system: [`EwaldShort`](@ref)
paired with [`PME3DLong`](@ref). Computes the same energy and forces as
[`Ewald3D`](@ref) at `O(N log N)` instead of `O(N·K)`.

Requires `using FINUFFT`.
"""
function PME3D end

"""
    ICMPME3D(n_atoms, L; α, s, γ, N_image, N_pad, ϵ = 1.0)

The image-charge method combined with particle-mesh Ewald and the electrostatic
layer correction — the `O(N log N)` counterpart of [`ICMEwald3D`](@ref).

Requires `using FINUFFT`.
"""
function ICMPME3D end

# Fallbacks so a user who forgot the extension gets a sentence rather than a
# MethodError listing zero candidates.
const _FINUFFT_HINT = "requires FINUFFT. Run `using FINUFFT` (and add it to your " *
                      "project) to load ExTinyMD's particle-mesh extension."

PME3D(args...; kwargs...)    = error("PME3D " * _FINUFFT_HINT)
ICMPME3D(args...; kwargs...) = error("ICMPME3D " * _FINUFFT_HINT)
PME3DLong(args...; kwargs...) = error("PME3DLong " * _FINUFFT_HINT)
```

In `src/ExTinyMD.jl`, add the include **before** `adapter.jl` (which must stay last):

```julia
include("interactions/electrostatics/long_pme3d.jl")
```

and export:

```julia
export PME3DLong, PME3D, ICMPME3D
```

- [ ] **Step 5: Write the extension's constructor and energy**

Create `ext/ExTinyMDFINUFFTExt.jl`:

```julia
module ExTinyMDFINUFFTExt

using ExTinyMD
using ExTinyMD: PME3DLong, EwaldShort, Periodic3D, ewald_cutoffs, EwaldInteraction,
                ICM, ICMShort, long_energy, long_force!
using FINUFFT
using StaticArrays

# FINUFFT tolerance. 1e-14 keeps the NUFFT's own error far below the k-space
# truncation error, so PME3DLong is limited by the same cutoff as Ewald3DLong
# rather than by the transform.
const NUFFT_TOL = 1e-14

function ExTinyMD.PME3DLong(n_atoms::Int, L::NTuple{3,T}; α::T, s::T, ϵ::T = one(T),
                            ϵ_inf::T = T(Inf)) where {T}
    r_c, k_c = ewald_cutoffs(s, α)
    n_k = (ceil(Int, k_c / (2π / L[1])),
           ceil(Int, k_c / (2π / L[2])),
           ceil(Int, k_c / (2π / L[3])))
    dims = (2n_k[1] + 1, 2n_k[2] + 1, 2n_k[3] + 1)

    # Spherical mask: zero outside the cutoff and at k = 0, so the k-set matches
    # Ewald3DLong's exactly.
    D = zeros(T, dims)
    @inbounds for i in 1:dims[1], j in 1:dims[2], m in 1:dims[3]
        k_x = (i - n_k[1] - 1) * 2π / L[1]
        k_y = (j - n_k[2] - 1) * 2π / L[2]
        k_z = (m - n_k[3] - 1) * 2π / L[3]
        k = sqrt(k_x^2 + k_y^2 + k_z^2)
        if 0 < k <= k_c
            D[i, j, m] = exp(-k^2 / (4 * α^2)) / k^2
        end
    end

    plan1 = finufft_makeplan(1, [dims...], +1, 1, NUFFT_TOL, dtype = T)
    plan2 = finufft_makeplan(2, [dims...], -1, 1, NUFFT_TOL, dtype = T)

    return PME3DLong{T, typeof(plan1), typeof(plan2)}(
        α, r_c, k_c, ϵ, ϵ_inf, L, n_atoms, n_k,
        D,
        zeros(Complex{T}, dims), zeros(Complex{T}, dims),
        zeros(Complex{T}, dims), zeros(Complex{T}, dims), zeros(Complex{T}, dims),
        zeros(T, n_atoms), zeros(T, n_atoms), zeros(T, n_atoms),
        zeros(Complex{T}, n_atoms), zeros(Complex{T}, n_atoms), zeros(Complex{T}, n_atoms),
        plan1, plan2)
end

# Scale the first `m` positions into the plan's own buffers. Never touch `poses`.
function _scale!(long::PME3DLong{T}, poses, m::Int) where {T}
    @inbounds for j in 1:m
        p = poses[j]
        long.xs[j] = 2π * T(p[1]) / long.L[1]
        long.ys[j] = 2π * T(p[2]) / long.L[2]
        long.zs[j] = 2π * T(p[3]) / long.L[3]
    end
    return nothing
end

# Type-1 transform of the first `m` charges into `out`. `out` is zeroed by FINUFFT.
function _structure_factor!(out, long::PME3DLong{T}, poses, charges, m::Int,
                            plan) where {T}
    _scale!(long, poses, m)
    # Always re-set the points: positions move every timestep, so there is no
    # count-based shortcut worth taking.
    finufft_setpts!(plan, view(long.xs, 1:m), view(long.ys, 1:m), view(long.zs, 1:m))
    q = Complex{T}[charges[j] for j in 1:m]
    finufft_exec!(plan, q, out)
    return out
end

function ExTinyMD.long_energy(long::PME3DLong{T}, poses, charges;
                              n_target::Int = long.n_atoms) where {T}
    n = long.n_atoms
    _structure_factor!(long.ρ_src, long, poses, charges, n, long.plan1)

    V = long.L[1] * long.L[2] * long.L[3]
    E = zero(T)

    if n_target == n
        @inbounds for idx in eachindex(long.D)
            E += abs2(long.ρ_src[idx]) * long.D[idx]
        end
    else
        _structure_factor!(long.ρ_tgt, long, poses, charges, n_target, long.plan1)
        @inbounds for idx in eachindex(long.D)
            E += real(conj(long.ρ_src[idx]) * long.ρ_tgt[idx]) * long.D[idx]
        end
        # ρ_src was overwritten by the second transform's setpts? No — the output
        # array differs, but the plan's points changed. Recompute nothing here;
        # long_force! re-derives ρ_src itself.
    end
    E /= (2 * V * long.ϵ)

    # Surface term, identical in form to Ewald3DLong's.
    P_src = _dipole(poses, charges, n, T)
    P_tgt = n_target == n ? P_src : _dipole(poses, charges, n_target, T)
    E += dot(P_tgt, P_src) / (2 * V * long.ϵ * (2 * long.ϵ_inf + one(T)))

    return E
end

_dipole(poses, charges, m::Int, ::Type{T}) where {T} =
    sum(k -> charges[k] * SVector{3,T}(T(poses[k][1]), T(poses[k][2]), T(poses[k][3])),
        1:m; init = zero(SVector{3,T}))

end # module
```

**Two things to get right here**, both of which the controller hit while validating:

1. `_structure_factor!` calls `finufft_setpts!` on every invocation, because the positions move every timestep. There is no count-based shortcut worth taking, and no Ref tracking the last point count — an earlier draft of this plan carried `plan1_npts`/`plan2_npts` fields for that purpose and they were deleted during Task 1 as write-only.

2. `long.ρ_src` and `long.ρ_tgt` must be distinct arrays, and the second transform must not clobber the first. FINUFFT writes its whole output array, so no manual zeroing is needed for `ρ_*` — but `hx/hy/hz` in Task 2 **are** only partially written and must be zeroed explicitly.

- [ ] **Step 6: Check the no-FINUFFT fallback**

In a separate environment without FINUFFT:

```
julia --startup-file=no -e 'using Pkg; Pkg.activate(mktempdir()); Pkg.develop(path="."); using ExTinyMD; try; PME3D(2, (10.0,10.0,10.0); α=1.0, s=4.0); catch e; println(sprint(showerror, e)); end'
```

Expected: the message from `_FINUFFT_HINT`, not an `UndefVarError` or a bare `MethodError`. Paste the actual output into your report.

- [ ] **Step 7: Run tests**

```
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS, including the three new PME testsets. Report the measured agreement figures.

- [ ] **Step 8: Commit**

```bash
git add Project.toml src/ ext/ test/
git commit -m "feat: add PME3DLong, a FINUFFT-backed reciprocal-space solver

Type-1 NUFFT for the structure factor with a spherical D_k mask, so the k-set
is identical to Ewald3DLong's and the two agree to machine precision rather
than merely converging to the same limit. FINUFFT is a weak dependency."
```

---

### Task 2: `PME3DLong` forces

**Files:**
- Modify: `ext/ExTinyMDFINUFFTExt.jl`, `test/electrostatics/test_pme3d.jl`

**Interfaces:**
- Consumes: Task 1's struct and `ρ_src`.
- Produces: `ExTinyMD.long_force!(F, long::PME3DLong{T}, poses, charges; n_target = long.n_atoms)` — **accumulates**.

**The formulation**, validated by the controller against `Ewald3DLong` to 5.9e-15:

```
F_i[d] = −(q_i /(V ϵ)) · Im( Σ_k k_d · D_k · ρ_src[k] · e^{−i k·r_i} )
```

The inner sum is a type-2 transform whose input is `k_d · D_k · ρ_src` — three transforms, one per Cartesian component, all sharing one plan. Reusing a single plan across the three `exec` calls is safe (measured: exactly 0.0 difference against three fresh plans).

- [ ] **Step 1: Write the failing test**

Append to `test/electrostatics/test_pme3d.jl`:

```julia
@testset "PME3DLong force equals Ewald3DLong force" begin
    # Controller measured 5.9e-15 relative / 2.6e-17 absolute on this shape.
    Random.seed!(7)
    n = 25
    L = (12.0, 12.0, 12.0)
    α, s = 0.8, 3.5
    poses = [SVector(rand() * L[1], rand() * L[2], rand() * L[3]) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]

    ewald = Ewald3DLong(n, L; α = α, s = s)
    pme   = PME3DLong(n, L; α = α, s = s)

    for nt in (n, 12)
        Fe = [zero(SVector{3,Float64}) for _ in 1:n]
        Fp = [zero(SVector{3,Float64}) for _ in 1:n]
        long_force!(Fe, ewald, poses, charges; n_target = nt)
        long_force!(Fp, pme,   poses, charges; n_target = nt)
        for i in 1:nt, d in 1:3
            @test isapprox(Fp[i][d], Fe[i][d], rtol = 1e-10, atol = 1e-16)
        end
    end
end

@testset "PME3DLong force accumulates, does not zero" begin
    # The composite coulomb_force! zeroes once and lets short and long accumulate.
    Random.seed!(8)
    n = 10
    L = (12.0, 12.0, 12.0)
    poses = [SVector(rand() * L[1], rand() * L[2], rand() * L[3]) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]
    pme = PME3DLong(n, L; α = 0.8, s = 3.5)

    F1 = [zero(SVector{3,Float64}) for _ in 1:n]
    long_force!(F1, pme, poses, charges)
    F2 = [zero(SVector{3,Float64}) for _ in 1:n]
    long_force!(F2, pme, poses, charges)
    long_force!(F2, pme, poses, charges)
    for i in 1:n
        @test isapprox(F2[i], 2 .* F1[i], rtol = 1e-12)
    end
end
```

- [ ] **Step 2: Run to verify it fails**

Expected: FAIL — `long_force!` has no `PME3DLong` method, so the buffers stay zero and the comparison against `Ewald3DLong` fails (or a `MethodError`, depending on dispatch).

- [ ] **Step 3: Implement**

Add to `ext/ExTinyMDFINUFFTExt.jl`:

```julia
function ExTinyMD.long_force!(F::Vector{SVector{3,T}}, long::PME3DLong{T}, poses,
                              charges; n_target::Int = long.n_atoms) where {T}
    n = long.n_atoms
    _structure_factor!(long.ρ_src, long, poses, charges, n, long.plan1)

    n_k = long.n_k
    dims = size(long.D)

    # Only masked-in grid points are written below, so these MUST be zeroed first.
    # Leaving them uninitialised puts NaN on the grid, and the type-2 transform
    # then spreads it to every particle.
    fill!(long.hx, zero(Complex{T}))
    fill!(long.hy, zero(Complex{T}))
    fill!(long.hz, zero(Complex{T}))

    @inbounds for i in 1:dims[1], j in 1:dims[2], m in 1:dims[3]
        d = long.D[i, j, m]
        iszero(d) && continue
        k_x = (i - n_k[1] - 1) * 2π / long.L[1]
        k_y = (j - n_k[2] - 1) * 2π / long.L[2]
        k_z = (m - n_k[3] - 1) * 2π / long.L[3]
        g = long.ρ_src[i, j, m] * d
        long.hx[i, j, m] = k_x * g
        long.hy[i, j, m] = k_y * g
        long.hz[i, j, m] = k_z * g
    end

    _scale!(long, poses, n_target)
    finufft_setpts!(long.plan2, view(long.xs, 1:n_target),
                    view(long.ys, 1:n_target), view(long.zs, 1:n_target))

    ox = view(long.ox, 1:n_target)
    oy = view(long.oy, 1:n_target)
    oz = view(long.oz, 1:n_target)
    finufft_exec!(long.plan2, long.hx, ox)
    finufft_exec!(long.plan2, long.hy, oy)
    finufft_exec!(long.plan2, long.hz, oz)

    V = long.L[1] * long.L[2] * long.L[3]
    pref = one(T) / (V * long.ϵ)
    @inbounds for i in 1:n_target
        F[i] -= pref * charges[i] * SVector{3,T}(imag(ox[i]), imag(oy[i]), imag(oz[i]))
    end

    # Surface term, matching Ewald3DLong: F_i = -q_i P_src /(V ϵ (2ϵ_inf + 1))
    P_src = _dipole(poses, charges, n, T)
    surf = one(T) / (V * long.ϵ * (2 * long.ϵ_inf + one(T)))
    @inbounds for i in 1:n_target
        F[i] -= (surf * charges[i]) * P_src
    end

    return F
end
```

If `finufft_exec!` refuses a `view` for its output, allocate `ox/oy/oz` at full length and index `1:n_target` when accumulating; say which you did in your report.

- [ ] **Step 4: Run tests**

Expected: PASS. Report the worst force difference you measure.

- [ ] **Step 5: Commit**

```bash
git add ext/ test/
git commit -m "feat: add PME3D forces via type-2 NUFFT

Three type-2 transforms of k_d * D_k * rho, sharing one plan. ParticleMeshEwald
has no force implementation at all, so this is new work; it matches
Ewald3DLong's analytic force to 5.9e-15."
```

---

### Task 3: `PME3D` constructor and MD integration

**Files:**
- Modify: `ext/ExTinyMDFINUFFTExt.jl`, `test/electrostatics/test_pme3d.jl`

**Interfaces:**
- Produces: `ExTinyMD.PME3D(n_atoms, L; α, s, ϵ=one(T), ϵ_inf=T(Inf))::EwaldInteraction`

- [ ] **Step 1: Write the failing test**

```julia
@testset "PME3D composite matches Ewald3D" begin
    Random.seed!(1234)
    n = 30
    L = (12.0, 12.0, 12.0)
    α, s = 0.8, 3.5
    poses = [SVector(rand() * L[1], rand() * L[2], rand() * L[3]) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]

    e = Ewald3D(n, L; α = α, s = s)
    p = PME3D(n, L; α = α, s = s)
    @test p isa ExTinyMD.AbstractInteraction

    @test isapprox(coulomb_energy(p, poses, charges),
                   coulomb_energy(e, poses, charges), rtol = 1e-12)

    Fe = coulomb_force(e, poses, charges)
    Fp = coulomb_force(p, poses, charges)
    for i in 1:n, d in 1:3
        @test isapprox(Fp[i][d], Fe[i][d], rtol = 1e-10, atol = 1e-16)
    end
end

@testset "PME3D drives simulate! through the adapter" begin
    # The adapter and EwaldShort are reused unchanged; this confirms a PME3D
    # composite is a drop-in for Ewald3D in an actual MD run.
    Random.seed!(20260915)
    n, L = 30, 12.0
    boundary = Boundary((L, L, L), (1, 1, 1))
    atoms = Vector{Atom{Float64}}()
    for i in 1:(n ÷ 2); push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0)); end
    for i in (n ÷ 2 + 1):n; push!(atoms, Atom(type = 2, mass = 1.0, charge = -1.0)); end
    info = SimulationInfo(n, atoms, (0.0, L, 0.0, L, 0.0, L), boundary;
                          min_r = 1.0, temp = 1.0)
    info.running_step = 1

    inter = PME3D(n, (L, L, L); α = 0.8, s = 3.5)   # r_c = 4.375 < 6
    finder = CellList3D(info, inter.short.r_c, boundary, 1)
    sys = MDSys(n_atoms = n, atoms = atoms, boundary = boundary,
                interactions = [(inter, finder)],
                loggers = [TemperatureLogger(1000; output = false)],
                simulator = VerletProcess(dt = 1e-4))

    E0 = energy(inter, finder, sys, info)
    simulate!(sys.simulator, sys, info, 200)
    E1 = energy(inter, finder, sys, info)
    @test isfinite(E1)
    @test abs(E1 - E0) < 0.05 * max(abs(E0), 1.0)
end
```

- [ ] **Step 2: Run to verify it fails**

Expected: the `error("PME3D requires FINUFFT...")` fallback fires, because the extension has not yet added a real method.

- [ ] **Step 3: Implement**

```julia
function ExTinyMD.PME3D(n_atoms::Int, L::NTuple{3,T}; α::T, s::T, ϵ::T = one(T),
                        ϵ_inf::T = T(Inf)) where {T}
    short = EwaldShort(n_atoms, L; α = α, s = s, ϵ = ϵ, convention = Periodic3D())
    long  = PME3DLong(n_atoms, L; α = α, s = s, ϵ = ϵ, ϵ_inf = ϵ_inf)
    return EwaldInteraction(short, long, n_atoms)
end
```

Note `EwaldInteraction`'s constructor infers `T` from `short.α` and allocates the adapter's gather buffers, so nothing further is needed for MD integration.

- [ ] **Step 4: Run tests, then commit**

```bash
git add ext/ test/
git commit -m "feat: add the PME3D constructor and confirm MD integration"
```

---

### Task 4: `ICMPME3D`

**Files:**
- Modify: `ext/ExTinyMDFINUFFTExt.jl`, `test/electrostatics/test_pme3d.jl`

**Interfaces:**
- Produces: `ExTinyMD.ICMPME3D(n_atoms, L; α, s, γ, N_image, N_pad, ϵ=one(T))::ICM`

This is where `n_target` earns its keep: `ICM` calls `long_energy`/`long_force!` with `n_target = icm.n_atoms` on the reflected configuration, and Task 1 and 2 already implement that for `PME3DLong`. So `ICMPME3D` is the same shape as `ICMEwald3D` with the solver swapped.

- [ ] **Step 1: Write the failing test**

```julia
@testset "ICMPME3D matches ICMEwald3D" begin
    # Same physics, same k-set, different reciprocal-space engine.
    Random.seed!(20260928)
    n = 8
    L = (5.0, 5.0, 10.0)
    γ = (0.3, 0.3)
    poses = [SVector(rand() * L[1], rand() * L[2], 2.0 + 6.0 * rand()) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]

    # r_c = 2.35 < min(Lx,Ly)/2 = 2.5; N_image = 3 needs N_pad >= 2
    a = ICMEwald3D(n, L; α = 1.7, s = 4.0, γ = γ, N_image = 3, N_pad = 2)
    b = ICMPME3D(n, L;  α = 1.7, s = 4.0, γ = γ, N_image = 3, N_pad = 2)

    @test isapprox(coulomb_energy(b, poses, charges),
                   coulomb_energy(a, poses, charges), rtol = 1e-10)

    Fa = coulomb_force(a, poses, charges)
    Fb = coulomb_force(b, poses, charges)
    for i in 1:n, d in 1:3
        @test isapprox(Fb[i][d], Fa[i][d], rtol = 1e-8, atol = 1e-14)
    end
end

@testset "ICMPME3D force matches -grad(energy)" begin
    # Independent of the ICMEwald3D comparison above: a finite-difference check
    # constrains the PME force path on its own.
    Random.seed!(20260929)
    n = 6
    L = (6.0, 6.0, 5.0)
    γ = (0.8, 0.8)
    poses = [SVector(1.0, 1.0, 0.5), SVector(3.0, 1.5, 4.6), SVector(1.5, 3.5, 2.5),
             SVector(4.0, 4.0, 0.7), SVector(2.0, 4.5, 4.3), SVector(4.5, 2.0, 2.0)]
    charges = [1.0, -1.0, 1.0, -1.0, 1.0, -1.0]

    inter = ICMPME3D(n, L; α = 1.3, s = 3.5, γ = γ, N_image = 3, N_pad = 2)
    F = coulomb_force(inter, poses, charges)
    f = p -> coulomb_energy(inter, p, charges)
    for i in 1:n, d in 1:3
        @test isapprox(F[i][d], -fd_gradient(f, poses, i, d; h = 1e-5),
                       rtol = 1e-4, atol = 1e-8)
    end
end
```

- [ ] **Step 2: Run to verify it fails**, then implement:

```julia
function ExTinyMD.ICMPME3D(n_atoms::Int, L::NTuple{3,T}; α::T, s::T, γ::Tuple{T,T},
                           N_image::Int, N_pad::Int, ϵ::T = one(T)) where {T}
    n_ref = n_atoms * (1 + 2 * N_image)
    L_pad = (L[1], L[2], (2 * N_pad + 1) * L[3])
    long  = PME3DLong(n_ref, L_pad; α = α, s = s, ϵ = ϵ, ϵ_inf = T(Inf))
    short = ICMShort(n_atoms, L; α = α, s = s, ϵ = ϵ, N_image = N_image)
    return ICM(long, short, γ, N_image, n_atoms, L; elc = true, N_pad = N_pad)
end
```

Note the long solver is built with `n_ref` and the **z-padded** box, exactly as `ICMEwald3D` does — its `n_atoms` field is what its source loops and buffers are sized by.

- [ ] **Step 3: Run tests, report the measured agreement, then commit.**

---

### Task 5: Documentation

**Files:**
- Modify: `docs/src/electrostatics.md`, `README.md`

- [ ] **Step 1: Update the method-selection table**

`PME3D` and `ICMPME3D` are currently described as forthcoming. Replace that with real rows:

| Boundary condition | Dielectric walls | Method | Cost | Notes |
|---|---|---|---|---|
| Triply periodic | none | `Ewald3D` | `O(N·K)` | reference implementation |
| Triply periodic | none | `PME3D` | `O(N log N)` | needs `using FINUFFT` |
| Slab, periodic in x,y | none | `Ewald2D` | `O(N²K)` | exact; accuracy reference only |
| Slab, periodic in x,y | confined | `ICMEwald2D` | `O(N²K)` | exact for confined slabs |
| Slab, periodic in x,y | confined | `ICMEwald3D` | `O(N·K)` | Ewald3D + ELC |
| Slab, periodic in x,y | confined | `ICMPME3D` | `O(N log N)` | PME3D + ELC; needs `using FINUFFT` |

- [ ] **Step 2: Add a "Particle mesh" section** covering:
  - `using FINUFFT` is required, and what the error says if you forget
  - that `PME3D` computes **the same sum** as `Ewald3D`, not an approximation of it — the spherical `D_k` mask gives both an identical k-set, and they agree to ~1e-15. Say this plainly; it is the most useful fact a user can have when choosing between them.
  - the practical consequence: pick `Ewald3D` for small systems and as a reference, `PME3D` when `N·K` starts to hurt. Do not invent a crossover particle count — if you want to quote one, measure it and say how.
  - that `α` and `s` mean exactly what they do for `Ewald3D`, including the `r_c < min(L)/2` rule.

- [ ] **Step 3: Update the README's Electrostatics section** to list the two new methods and note the FINUFFT requirement.

- [ ] **Step 4: Rebuild docs, confirm clean, commit.**

```
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=".")); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

The docs environment will need FINUFFT for any `@docs` block referencing the extension's methods; add it to `docs/Project.toml` if the build complains, and say so in your report.

---

## Self-Review

**Spec coverage.** §5.7 is implemented by Tasks 1–4; its three noted defects in `ParticleMeshEwald` are addressed as follows: in-place mutation of caller coordinates → `_scale!` into plan-owned buffers (Task 1); `Threads.threadid()` indexing → no threading at all, and the real-space part is `EwaldShort`, already fixed in Phase 1; box-vs-sphere truncation → the spherical mask, which turns the problem into an exact test (spec amended in `93135b8`). §8.2's PME row is implemented by Task 1's first testset.

**Placeholders.** None. Every code step carries the code; the one judgement call left open (a crossover particle count in the docs) is explicitly marked as "measure it or omit it".

**Type consistency.** `PME3DLong` is constructed in the extension and its struct declared in the parent with free plan parameters. `long_energy`/`long_force!` are methods on ExTinyMD's existing generics, matching `Ewald3DLong`'s signatures including `n_target`. `PME3D` returns an `EwaldInteraction` and `ICMPME3D` returns an `ICM`, so both reuse the Phase 1 `coulomb_*` methods and the adapter unchanged.

**Resolved during Task 1.** The `plan1_npts` / `plan2_npts` Refs this plan originally carried
were deleted as write-only: `setpts!` must be called every timestep regardless of point count.
A later measurement showed skipping it would in fact have been harmless — FINUFFT.jl's guru
binding keeps a live reference to the position buffers rather than snapshotting them at
`setpts!` time, so `exec` re-reads whatever is in them. Calling it unconditionally avoids
depending on that undocumented detail, which is the reason to prefer it.

**On measuring force agreement.** Two error metrics are in play and they differ by three
orders of magnitude on the same data. Per *particle* — dividing by that particle's largest
component — gives ~6e-15. Per *component* — dividing by each component individually, some of
which are as small as 2.5e-5 against a largest of 0.034 — gives ~1e-12. The absolute error is
1.25e-16 either way. Quote which metric you mean; a number without it invites a false alarm.
