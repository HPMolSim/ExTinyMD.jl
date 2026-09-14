"""
    Ewald3DLong(n_atoms, L; α, s, ϵ = 1.0, ϵ_inf = Inf)

Reciprocal-space part of the triply periodic Ewald sum, evaluated as a direct sum
over the k-set:

    E_l = 1/(2Vϵ) Σ_{0<|k|≤k_c} |ρ_k|² exp(−k²/4α²)/k²,   ρ_k = Σ_j q_j exp(i k·r_j)

plus the surface (dipole) term `|P|²/(2Vϵ(2ϵ_inf+1))` with `P = Σ_j q_j r_j`.
`ϵ_inf = Inf` is the conducting (tin-foil) boundary and drops the surface term;
finite `ϵ_inf` applies the correction for a medium of that permittivity at infinity.

Cost is `O(N·K)`. For large systems use the FINUFFT-backed `PME3D` instead, which
computes the same quantity with the same normalisation.
"""
struct Ewald3DLong{T}
    α::T
    r_c::T
    k_c::T
    ϵ::T
    ϵ_inf::T
    L::NTuple{3,T}
    n_atoms::Int
    k_set::Vector{NTuple{4,T}}
end

function Ewald3DLong(n_atoms::Int, L::NTuple{3,T}; α::T, s::T, ϵ::T = one(T),
                     ϵ_inf::T = T(Inf)) where {T}
    r_c, k_c = ewald_cutoffs(s, α)
    return Ewald3DLong{T}(α, r_c, k_c, ϵ, ϵ_inf, L, n_atoms, k_set_3D(k_c, L))
end

Base.show(io::IO, l::Ewald3DLong) =
    print(io, "Ewald3DLong(α = $(l.α), k_c = $(l.k_c), ϵ = $(l.ϵ), ϵ_inf = $(l.ϵ_inf), " *
              "$(length(l.k_set)) k-vectors)")

# Σ_j q_j exp(i k·r_j), accumulated in Complex{T} so Float32 and higher precisions survive.
@inline function _structure_factor(k_x::T, k_y::T, k_z::T, poses, charges,
                                   n::Int) where {T}
    ρ = zero(Complex{T})
    @inbounds for j in 1:n
        p = poses[j]
        ρ += charges[j] * cis(k_x * T(p[1]) + k_y * T(p[2]) + k_z * T(p[3]))
    end
    return ρ
end

@inline function _dipole(poses, charges, n::Int, ::Type{T}) where {T}
    P = zero(SVector{3,T})
    @inbounds for j in 1:n
        p = poses[j]
        P += charges[j] * SVector{3,T}(T(p[1]), T(p[2]), T(p[3]))
    end
    return P
end

"""
    long_energy(long, poses, charges; n_target = long.n_atoms) -> T

Reciprocal-space Ewald energy plus the dipole surface term. `n_target` lets the
targets (the first `n_target` entries of `poses`/`charges`) differ from the
sources (all `long.n_atoms` of them) — used by the image-charge method (Task 8),
where real particles are targets and all reflected charges are sources. When
`n_target == long.n_atoms` (the default) the two coincide and this reduces to
the ordinary Ewald energy.
"""
function long_energy(long::Ewald3DLong{T}, poses, charges;
                     n_target::Int = long.n_atoms) where {T}
    α, n = long.α, long.n_atoms          # n = number of SOURCE charges
    V = long.L[1] * long.L[2] * long.L[3]

    E = zero(T)
    @inbounds for (k_x, k_y, k_z, k) in long.k_set
        ρ_src = _structure_factor(k_x, k_y, k_z, poses, charges, n)
        # real(conj(ρ_src) * ρ_tgt) collapses to abs2(ρ) when the sets coincide,
        # so there is one code path, not two.
        ρ_tgt = n_target == n ? ρ_src :
                _structure_factor(k_x, k_y, k_z, poses, charges, n_target)
        E += real(conj(ρ_src) * ρ_tgt) * exp(-k^2 / (4 * α^2)) / k^2
    end
    E /= (2 * V * long.ϵ)

    # Surface term. 1/(2*Inf+1) evaluates to zero, so the conducting case needs no branch.
    P_src = _dipole(poses, charges, n, T)
    P_tgt = n_target == n ? P_src : _dipole(poses, charges, n_target, T)
    E += dot(P_tgt, P_src) / (2 * V * long.ϵ * (2 * long.ϵ_inf + one(T)))

    return E
end

"""
    long_force!(F, long, poses, charges; n_target = long.n_atoms)

Accumulate the reciprocal-space force into `F`. **Does not zero `F` first.**
Sources span all `long.n_atoms` charges; the force is written only for the
first `n_target` of them, using the full coefficient (no factor of ½) — the
image-charge convention, self-consistent because an image moves at twice the
rate of its source.
"""
function long_force!(F::Vector{SVector{3,T}}, long::Ewald3DLong{T}, poses,
                     charges; n_target::Int = long.n_atoms) where {T}
    α, n = long.α, long.n_atoms          # n = number of SOURCE charges
    V = long.L[1] * long.L[2] * long.L[3]
    pref = one(T) / (V * long.ϵ)

    @inbounds for (k_x, k_y, k_z, k) in long.k_set
        ρ_src = _structure_factor(k_x, k_y, k_z, poses, charges, n)
        D = exp(-k^2 / (4 * α^2)) / k^2
        kvec = SVector{3,T}(k_x, k_y, k_z)
        for i in 1:n_target
            p = poses[i]
            phase = cis(-(k_x * T(p[1]) + k_y * T(p[2]) + k_z * T(p[3])))
            F[i] -= (pref * charges[i] * D * imag(ρ_src * phase)) * kvec
        end
    end

    # Surface term force: F_i = -q_i P_src / (V ϵ (2ϵ_inf + 1))
    P_src = _dipole(poses, charges, n, T)
    surf = one(T) / (V * long.ϵ * (2 * long.ϵ_inf + one(T)))
    @inbounds for i in 1:n_target
        F[i] -= (surf * charges[i]) * P_src
    end

    return F
end
