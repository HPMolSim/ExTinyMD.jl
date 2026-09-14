"""
    Ewald2DLong(n_atoms, L; α, s, ϵ = 1.0)

Reciprocal-space part of the exact Ewald sum for a slab periodic in x and y and
free in z. For each in-plane wavevector `k`:

    E_k  = 1/ϵ Σ_i Σ_j q_i q_j cos(k·ρ_ij)
             [ e^{k z_ij} erfc(k/2α + α z_ij) + e^{−k z_ij} erfc(k/2α − α z_ij) ]
             / (8 L_x L_y k)

    E_k0 = −1/ϵ Σ_i Σ_j q_i q_j [ e^{−(α z_ij)²}/(α√π) + z_ij erf(α z_ij) ]
             / (4 L_x L_y)

!!! note "Prefactor"
    The long-range prefactor is `1/ϵ`, not `1/(4πϵ)` — the `4π` is already folded
    into the expressions above. The short-range part of the same method does carry
    `1/(4πϵ)`. This asymmetry is deliberate; α-independence of the total energy is
    the test that pins it.

Cost is `O(N²K)`: this is the *exact* 2D sum, intended as the accuracy reference
for quasi-2D systems rather than as a production method for large `N`.
"""
struct Ewald2DLong{T}
    α::T
    r_c::T
    k_c::T
    ϵ::T
    L::NTuple{3,T}
    n_atoms::Int
    k_set::Vector{NTuple{3,T}}
end

function Ewald2DLong(n_atoms::Int, L::NTuple{3,T}; α::T, s::T, ϵ::T = one(T)) where {T}
    r_c, k_c = ewald_cutoffs(s, α)
    return Ewald2DLong{T}(α, r_c, k_c, ϵ, L, n_atoms, k_set_2D(k_c, L))
end

Base.show(io::IO, l::Ewald2DLong) =
    print(io, "Ewald2DLong(α = $(l.α), k_c = $(l.k_c), ϵ = $(l.ϵ), " *
              "$(length(l.k_set)) k-vectors)")

# exp(±k z) erfc(k/2α ± α z), guarded. For large positive argument the exp overflows
# while the erfc underflows; the product tends to zero, so return zero rather than
# letting Inf * 0 produce NaN.
@inline function _exp_erfc(kz::T, arg::T) where {T}
    kz > T(600) && return zero(T)
    return exp(kz) * erfc(arg)
end

"""
    long_energy(long::Ewald2DLong, poses, charges; n_target = long.n_atoms) -> T

Reciprocal-space Ewald energy for the quasi-2D slab. `n_target` lets the targets
(the first `n_target` entries of `poses`/`charges`) differ from the sources (all
`long.n_atoms` of them) — used by the image-charge method (Task 8), where real
particles are targets and all reflected charges are sources. When
`n_target == long.n_atoms` (the default) the two coincide and this reduces to the
ordinary Ewald2D energy.
"""
function long_energy(long::Ewald2DLong{T}, poses, charges;
                     n_target::Int = long.n_atoms) where {T}
    α, n = long.α, long.n_atoms          # n = number of SOURCE charges
    A = long.L[1] * long.L[2]

    E = zero(T)

    # k = 0 term. Targets i run to n_target, sources j over all n.
    @inbounds for i in 1:n_target, j in 1:n
        z = T(poses[i][3]) - T(poses[j][3])
        E -= charges[i] * charges[j] *
             (exp(-(α * z)^2) / (α * sqrt(T(π))) + z * erf(α * z)) / (4 * A)
    end

    # k != 0 terms
    @inbounds for (k_x, k_y, k) in long.k_set
        acc = zero(T)
        for i in 1:n_target
            p_i = poses[i]
            for j in 1:n
                p_j = poses[j]
                x = T(p_i[1]) - T(p_j[1])
                y = T(p_i[2]) - T(p_j[2])
                z = T(p_i[3]) - T(p_j[3])
                g = _exp_erfc(k * z, k / (2α) + α * z) +
                    _exp_erfc(-k * z, k / (2α) - α * z)
                acc += charges[i] * charges[j] * cos(k_x * x + k_y * y) * g
            end
        end
        E += acc / (8 * A * k)
    end

    return E / long.ϵ
end

"""
    long_force!(F, long::Ewald2DLong, poses, charges; n_target = long.n_atoms)

Accumulate the reciprocal-space force into `F`. **Does not zero `F` first.**
Sources span all `long.n_atoms` charges; the force is written only for the
first `n_target` of them, using the full coefficient (no factor of ½) — the
image-charge convention, self-consistent because an image moves at twice the
rate of its source.
"""
function long_force!(F::Vector{SVector{3,T}}, long::Ewald2DLong{T}, poses,
                     charges; n_target::Int = long.n_atoms) where {T}
    α, n = long.α, long.n_atoms          # n = number of SOURCE charges
    A = long.L[1] * long.L[2]
    ϵ = long.ϵ

    # k = 0 term: F_i,z = q_i Σ_j q_j erf(α z_ij) / (2 L_x L_y ϵ)
    @inbounds for i in 1:n_target
        fz = zero(T)
        for j in 1:n
            z = T(poses[i][3]) - T(poses[j][3])
            fz += charges[j] * erf(α * z)
        end
        F[i] += SVector{3,T}(zero(T), zero(T), charges[i] * fz / (2 * A * ϵ))
    end

    # k != 0 terms
    @inbounds for (k_x, k_y, k) in long.k_set
        for i in 1:n_target
            p_i = poses[i]
            sx = zero(T); sy = zero(T); sz = zero(T)
            for j in 1:n
                p_j = poses[j]
                x = T(p_i[1]) - T(p_j[1])
                y = T(p_i[2]) - T(p_j[2])
                z = T(p_i[3]) - T(p_j[3])
                qq = charges[i] * charges[j]
                phase = k_x * x + k_y * y

                ee_p = _exp_erfc(k * z, k / (2α) + α * z)    # e^{kz}  erfc(k/2α + αz)
                ee_m = _exp_erfc(-k * z, k / (2α) - α * z)   # e^{-kz} erfc(k/2α - αz)

                sum_xy = -qq * sin(phase) * (ee_p + ee_m)
                sx += k_x * sum_xy / k
                sy += k_y * sum_xy / k

                # d/dz of (e^{kz} erfc(k/2α+αz) + e^{-kz} erfc(k/2α-αz)), same
                # overflow guard as _exp_erfc: exp overflows while the Gaussian
                # underflows, and the product tends to zero.
                gauss_p = k * z > T(600) ? zero(T) :
                          exp(k * z) * exp(-(k / (2α) + α * z)^2)
                gauss_m = -k * z > T(600) ? zero(T) :
                          exp(-k * z) * exp(-(k / (2α) - α * z)^2)
                sz += qq * cos(phase) * (k * ee_p - k * ee_m -
                                         2α / sqrt(T(π)) * gauss_p +
                                         2α / sqrt(T(π)) * gauss_m) / k
            end
            F[i] -= SVector{3,T}(sx, sy, sz) / (4 * A * ϵ)
        end
    end

    return F
end
