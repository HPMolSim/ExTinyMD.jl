# Boundary conventions. These select which axes wrap under the minimum-image
# convention, and nothing else.
abstract type AbstractBoundaryConvention end

"""
    Periodic3D()

Boundary convention for a triply periodic box: the minimum-image displacement
wraps all three axes. Used by [`Ewald3D`](@ref).
"""
struct Periodic3D  <: AbstractBoundaryConvention end

"""
    PeriodicQ2D()

Boundary convention for a quasi-2D slab: the minimum-image displacement wraps
x and y only, z is unwrapped. Used by [`Ewald2D`](@ref) and the ICM methods.
"""
struct PeriodicQ2D <: AbstractBoundaryConvention end

"""
    ewald_cutoffs(s, α) -> (r_c, k_c)

Real- and reciprocal-space cutoffs from the splitting parameter `α` and the
dimensionless accuracy parameter `s`: `r_c = s/α`, `k_c = 2αs`. Larger `s` means
more accuracy and more work in both spaces. This convention is shared with
EwaldSummations, QuasiEwald, SoEwald2D and ParticleMeshEwald, so parameters
transfer between them unchanged.
"""
ewald_cutoffs(s::T, α::T) where {T} = (s / α, 2 * α * s)

"""
    k_set_3D(k_c, L) -> Vector{NTuple{4,T}}

Reciprocal lattice vectors `(k_x, k_y, k_z, |k|)` of the box `L` with
`0 < |k| ≤ k_c`. The set is closed under `k -> -k`; the long-range sums assume
that and therefore carry no factor-of-two correction.
"""
function k_set_3D(k_c::T, L::NTuple{3,T}) where {T}
    mx_max = ceil(Int, k_c * L[1] / 2π) + 1
    my_max = ceil(Int, k_c * L[2] / 2π) + 1
    mz_max = ceil(Int, k_c * L[3] / 2π) + 1

    k_set = Vector{NTuple{4,T}}()
    for m_x in -mx_max:mx_max, m_y in -my_max:my_max, m_z in -mz_max:mz_max
        k_x = m_x * 2π / L[1]
        k_y = m_y * 2π / L[2]
        k_z = m_z * 2π / L[3]
        k = sqrt(k_x^2 + k_y^2 + k_z^2)
        if 0 < k <= k_c
            push!(k_set, (k_x, k_y, k_z, k))
        end
    end
    return k_set
end

"""
    k_set_2D(k_c, L) -> Vector{NTuple{3,T}}

In-plane reciprocal lattice vectors `(k_x, k_y, |k|)` with `0 < |k| ≤ k_c`.
`L[3]` is ignored: the z axis is not periodic in the quasi-2D geometry.
"""
function k_set_2D(k_c::T, L::NTuple{3,T}) where {T}
    mx_max = ceil(Int, k_c * L[1] / 2π) + 1
    my_max = ceil(Int, k_c * L[2] / 2π) + 1

    k_set = Vector{NTuple{3,T}}()
    for m_x in -mx_max:mx_max, m_y in -my_max:my_max
        k_x = m_x * 2π / L[1]
        k_y = m_y * 2π / L[2]
        k = sqrt(k_x^2 + k_y^2)
        if 0 < k <= k_c
            push!(k_set, (k_x, k_y, k))
        end
    end
    return k_set
end

"""
    check_neutrality(charges; atol) -> Σq

Return the net charge, warning once when it is non-zero. Ewald summation of a
non-neutral system is conditionally convergent and picks up a box-volume-dependent
offset, which shows up as an energy that drifts with `L` rather than as an error.
"""
function check_neutrality(charges::AbstractVector{T}; atol::T = sqrt(eps(T))) where {T}
    net = sum(charges)
    if abs(net) > atol
        @warn "System is not charge neutral (Σq = $net); Ewald energies carry a " *
              "volume-dependent offset." maxlog = 1
    end
    return net
end

# Nearest-image displacement. `dx - L*round(dx/L)` is the true minimum image, unlike
# ExTinyMD's position_check3D, which returns the first image inside the cutoff and is
# equivalent only while r_c < L/2.
@inline _wrap(dx::T, L::T) where {T} = dx - L * round(dx / L)

@inline function min_image_disp(p_i, p_j, L::NTuple{3,T}, ::Periodic3D) where {T}
    return SVector{3,T}(_wrap(T(p_i[1]) - T(p_j[1]), L[1]),
                        _wrap(T(p_i[2]) - T(p_j[2]), L[2]),
                        _wrap(T(p_i[3]) - T(p_j[3]), L[3]))
end

@inline function min_image_disp(p_i, p_j, L::NTuple{3,T}, ::PeriodicQ2D) where {T}
    return SVector{3,T}(_wrap(T(p_i[1]) - T(p_j[1]), L[1]),
                        _wrap(T(p_i[2]) - T(p_j[2]), L[2]),
                        T(p_i[3]) - T(p_j[3]))
end
