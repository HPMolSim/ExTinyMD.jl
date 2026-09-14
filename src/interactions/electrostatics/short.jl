"""
    EwaldShort(n_atoms, L; α, s, ϵ = 1.0, convention = Periodic3D())

Real-space part of an Ewald split, shared by every method in this library:

    E_s = 1/(4πϵ) [ Σ_{i<j, r<r_c} q_i q_j erfc(α r_ij)/r_ij − (α/√π) Σ_i q_i² ]

`convention` selects which axes wrap — [`Periodic3D`](@ref) for a triply periodic
box, [`PeriodicQ2D`](@ref) for a slab periodic in x and y only.

The struct owns a `CellListMap` neighbour list so that [`short_energy`](@ref) and
[`short_force!`](@ref) allocate nothing after construction. Pass
`neighbor_list = ...` to reuse a list maintained elsewhere, as the MD adapter does.
"""
mutable struct EwaldShort{T, C <: AbstractBoundaryConvention, TC}
    α::T
    r_c::T
    k_c::T
    ϵ::T
    L::NTuple{3,T}
    n_atoms::Int
    convention::C
    cell_list::TC
    pos_buffer::Vector{SVector{3,T}}
end

function EwaldShort(n_atoms::Int, L::NTuple{3,T}; α::T, s::T, ϵ::T = one(T),
                    convention::C = Periodic3D()) where {T, C <: AbstractBoundaryConvention}
    r_c, k_c = ewald_cutoffs(s, α)
    pos_buffer = [zero(SVector{3,T}) for _ in 1:n_atoms]
    cell_list = InPlaceNeighborList(xpositions = pos_buffer, cutoff = r_c,
                                    unitcell = _cell_unitcell(L, r_c, convention),
                                    parallel = true)
    return EwaldShort{T, C, typeof(cell_list)}(α, r_c, k_c, ϵ, L, n_atoms, convention,
                                               cell_list, pos_buffer)
end

Base.show(io::IO, s::EwaldShort) =
    print(io, "EwaldShort(α = $(s.α), r_c = $(s.r_c), ϵ = $(s.ϵ), $(s.convention))")

# For a non-periodic axis, inflate the unitcell so CellListMap finds no images along it.
_cell_unitcell(L::NTuple{3,T}, r_c::T, ::Periodic3D) where {T} =
    SVector{3,T}(L[1], L[2], L[3])
_cell_unitcell(L::NTuple{3,T}, r_c::T, ::PeriodicQ2D) where {T} =
    SVector{3,T}(L[1], L[2], max(L[3] + 2 * r_c, T(2) * r_c))

function _refresh_neighbors!(short::EwaldShort{T}, poses) where {T}
    @inbounds for i in 1:short.n_atoms
        p = poses[i]
        short.pos_buffer[i] = SVector{3,T}(T(p[1]), T(p[2]), T(p[3]))
    end
    update!(short.cell_list, xpositions = short.pos_buffer)
    return neighborlist!(short.cell_list)
end

"""
    short_energy(short, poses, charges; neighbor_list = nothing) -> T

Real-space Ewald energy. `poses` is AoS — any vector whose elements support
`p[1]`, `p[2]`, `p[3]`.
"""
function short_energy(short::EwaldShort{T}, poses, charges;
                      neighbor_list = nothing) where {T}
    nb = neighbor_list === nothing ? _refresh_neighbors!(short, poses) : neighbor_list
    α, r_c = short.α, short.r_c

    E = zero(T)
    @inbounds for (i, j, r) in nb
        (r < r_c && r > zero(T)) || continue
        E += charges[i] * charges[j] * erfc(α * r) / r
    end

    @inbounds for i in 1:short.n_atoms
        E -= charges[i]^2 * α / sqrt(T(π))
    end

    return E / (4π * short.ϵ)
end

# -dE/dr for E(r) = q_i q_j erfc(α r)/r
@inline function _short_pair_dEdr(q_i::T, q_j::T, α::T, r::T) where {T}
    return q_i * q_j * (erfc(α * r) / r^2 + 2α / sqrt(T(π)) * exp(-(α * r)^2) / r)
end

"""
    short_force!(F, short, poses, charges; neighbor_list = nothing)

Accumulate the real-space force into `F`. **Does not zero `F` first** — the
composite interaction zeroes once and lets short and long parts accumulate.
"""
function short_force!(F::Vector{SVector{3,T}}, short::EwaldShort{T}, poses, charges;
                      neighbor_list = nothing) where {T}
    nb = neighbor_list === nothing ? _refresh_neighbors!(short, poses) : neighbor_list
    α, r_c, ϵ = short.α, short.r_c, short.ϵ
    conv, L = short.convention, short.L
    pref = one(T) / (4π * ϵ)

    @inbounds for (i, j, r) in nb
        (r < r_c && r > zero(T)) || continue
        d = min_image_disp(poses[i], poses[j], L, conv)
        F_ij = _short_pair_dEdr(charges[i], charges[j], α, r) * d / r * pref
        F[i] += F_ij
        F[j] -= F_ij
    end
    return F
end
