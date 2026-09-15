module ExTinyMDFINUFFTExt

using ExTinyMD
using ExTinyMD: PME3DLong, EwaldShort, Periodic3D, ewald_cutoffs, EwaldInteraction,
                ICM, ICMShort, long_energy, long_force!, _dipole
using FINUFFT
using StaticArrays
using LinearAlgebra: dot

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

# Scale the first `m` positions into the plan's own buffers. Never touch `poses`:
# ParticleMeshEwald's energy_long scales the caller's array in place and divides
# back afterwards, which corrupts caller data if the transform throws.
function _scale!(long::PME3DLong{T}, poses, m::Int) where {T}
    @inbounds for j in 1:m
        p = poses[j]
        long.xs[j] = 2π * T(p[1]) / long.L[1]
        long.ys[j] = 2π * T(p[2]) / long.L[2]
        long.zs[j] = 2π * T(p[3]) / long.L[3]
    end
    return nothing
end

# Type-1 transform of the first `m` charges into `out`. `out` is zeroed by FINUFFT
# (it writes its whole output array), so no manual zeroing is needed here.
#
# `finufft_setpts!` is called unconditionally: positions move every timestep, not
# only when `m` changes, so there is no valid point at which to skip it. (An
# earlier draft carried a `Ref{Int}` per plan to skip `setpts!` when the point
# count was unchanged; it was deleted — see Task 1 report.)
function _structure_factor!(out, long::PME3DLong{T}, poses, charges, m::Int,
                            plan) where {T}
    _scale!(long, poses, m)
    finufft_setpts!(plan, view(long.xs, 1:m), view(long.ys, 1:m), view(long.zs, 1:m))
    q = Complex{T}[charges[j] for j in 1:m]
    finufft_exec!(plan, q, out)
    return out
end

"""
    long_energy(long::PME3DLong, poses, charges; n_target = long.n_atoms) -> T

Particle-mesh counterpart of `Ewald3DLong`'s `long_energy`: identical semantics
(`n_target` targets, all `long.n_atoms` sources), computed via a type-1 NUFFT
structure factor summed against the spherically masked Green's function `D`
instead of a direct sum over `k_set`.
"""
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
        # ρ_src and ρ_tgt are distinct arrays: this second transform must not
        # clobber the first, since both are needed below.
        _structure_factor!(long.ρ_tgt, long, poses, charges, n_target, long.plan1)
        @inbounds for idx in eachindex(long.D)
            E += real(conj(long.ρ_src[idx]) * long.ρ_tgt[idx]) * long.D[idx]
        end
    end
    E /= (2 * V * long.ϵ)

    # Surface term, identical in form to Ewald3DLong's.
    P_src = _dipole(poses, charges, n, T)
    P_tgt = n_target == n ? P_src : _dipole(poses, charges, n_target, T)
    E += dot(P_tgt, P_src) / (2 * V * long.ϵ * (2 * long.ϵ_inf + one(T)))

    return E
end

"""
    long_force!(F, long::PME3DLong, poses, charges; n_target = long.n_atoms)

Particle-mesh counterpart of `Ewald3DLong`'s `long_force!`: identical semantics
(accumulates into `F`, does not zero it; `n_target` targets, all `long.n_atoms`
sources), computed via three type-2 NUFFTs — one per Cartesian component — of
`k_d · D_k · ρ_src`, sharing a single plan, instead of a direct sum over `k_set`.
"""
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

function ExTinyMD.PME3D(n_atoms::Int, L::NTuple{3,T}; α::T, s::T, ϵ::T = one(T),
                        ϵ_inf::T = T(Inf)) where {T}
    short = EwaldShort(n_atoms, L; α = α, s = s, ϵ = ϵ, convention = Periodic3D())
    long  = PME3DLong(n_atoms, L; α = α, s = s, ϵ = ϵ, ϵ_inf = ϵ_inf)
    return EwaldInteraction(short, long, n_atoms)
end

end # module
