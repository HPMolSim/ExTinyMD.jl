"""
    EwaldInteraction(short, long, n_atoms)

An Ewald-split electrostatic interaction: a real-space kernel plus a
reciprocal-space solver. Composing them here rather than registering two separate
entries in `sys.interactions` keeps the pair together and hands the MD adapter a
single object.

Query it with [`coulomb_energy`](@ref) and [`coulomb_force`](@ref) for standalone
use, or add it to an `MDSys` and let the adapter drive it.
"""
struct EwaldInteraction{T, S, L} <: AbstractInteraction
    short::S
    long::L
    n_atoms::Int
    force_buffer::Vector{SVector{3,T}}
    # Gather buffers for the MD adapter (Task 9). Declared here, with the struct,
    # rather than bolted on later: the adapter is called every timestep and must
    # not allocate, and adding fields in a later task would force every consumer
    # of this struct to be re-tested for no benefit.
    pos_scratch::Vector{SVector{3,T}}
    charge_scratch::Vector{T}
end

function EwaldInteraction(short::S, long::L, n_atoms::Int) where {S, L}
    T = typeof(short.α)
    return EwaldInteraction{T, S, L}(short, long, n_atoms,
                                     [zero(SVector{3,T}) for _ in 1:n_atoms],
                                     [zero(SVector{3,T}) for _ in 1:n_atoms],
                                     zeros(T, n_atoms))
end

Base.show(io::IO, i::EwaldInteraction) =
    print(io, "EwaldInteraction($(i.n_atoms) atoms)\n  short: $(i.short)\n  long:  $(i.long)")

"""
    coulomb_energy(interaction, poses, charges; neighbor_list = nothing) -> T

Total electrostatic energy. `poses` is AoS; no ExTinyMD type is required.
"""
function coulomb_energy(inter::EwaldInteraction{T}, poses, charges;
                        neighbor_list = nothing) where {T}
    return short_energy(inter.short, poses, charges; neighbor_list = neighbor_list) +
           long_energy(inter.long, poses, charges)
end

"""
    coulomb_force!(F, interaction, poses, charges; neighbor_list = nothing) -> F

Total electrostatic force, written into `F`. `F` is zeroed first, then the short-
and long-range parts accumulate into it.
"""
function coulomb_force!(F::Vector{SVector{3,T}}, inter::EwaldInteraction{T}, poses,
                        charges; neighbor_list = nothing) where {T}
    fill!(F, zero(SVector{3,T}))
    short_force!(F, inter.short, poses, charges; neighbor_list = neighbor_list)
    long_force!(F, inter.long, poses, charges)
    return F
end

"""
    coulomb_force(interaction, poses, charges; neighbor_list = nothing) -> Vector{SVector{3,T}}

Allocating form of [`coulomb_force!`](@ref). In an MD loop prefer the in-place
version, or let the adapter reuse `interaction.force_buffer`.
"""
function coulomb_force(inter::EwaldInteraction{T}, poses, charges;
                       neighbor_list = nothing) where {T}
    F = [zero(SVector{3,T}) for _ in 1:inter.n_atoms]
    return coulomb_force!(F, inter, poses, charges; neighbor_list = neighbor_list)
end

"""
    Ewald3D(n_atoms, L; α, s, ϵ = 1.0, ϵ_inf = Inf)

Standard Ewald summation for a triply periodic system. `α` splits real and
reciprocal space, `s` sets the accuracy (`r_c = s/α`, `k_c = 2αs`).

```jldoctest
julia> using StaticArrays

julia> inter = Ewald3D(2, (10.0, 10.0, 10.0); α = 1.0, s = 4.0);

julia> poses = [SVector(0.0, 0.0, 0.0), SVector(5.0, 0.0, 0.0)];

julia> round(coulomb_energy(inter, poses, [1.0, -1.0]); digits = 6)
-0.021815
```
"""
function Ewald3D(n_atoms::Int, L::NTuple{3,T}; α::T, s::T, ϵ::T = one(T),
                 ϵ_inf::T = T(Inf)) where {T}
    short = EwaldShort(n_atoms, L; α = α, s = s, ϵ = ϵ, convention = Periodic3D())
    long  = Ewald3DLong(n_atoms, L; α = α, s = s, ϵ = ϵ, ϵ_inf = ϵ_inf)
    return EwaldInteraction(short, long, n_atoms)
end
