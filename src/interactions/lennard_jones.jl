"""
    LennardJones(; ϵ = 1.0, cutoff = 3.5, σ = 1.0)

Pairwise Lennard-Jones interaction, `E(r) = 4ϵ[(σ/r)^12 - (σ/r)^6]` for
`r < cutoff`, evaluated over the pairs in the associated neighbour finder.
"""
struct LennardJones{T} <: AbstractInteraction
    ϵ::T
    cutoff::T
    σ::T
end

Base.show(io::IO, interaction::LennardJones) = print(io, "LennardJones with ϵ = $(interaction.ϵ), cutoff = $(interaction.cutoff), σ = $(interaction.σ)")

LennardJones(;ϵ::T = 1.0, cutoff::T = 3.5, σ::T = 1.0) where T = LennardJones(ϵ, cutoff, σ)

function update_acceleration!(interaction::LennardJones{T}, neighborfinder::T_NIEGHBER, sys::MDSys{T}, info::SimulationInfo{T}) where {T<:Number, T_NIEGHBER<:AbstractNeighborFinder}
    atoms = sys.atoms
    boundary = sys.boundary
    update_finder!(neighborfinder, info)
    for (i, j, r) in neighborfinder.neighbor_list
        coord_1, coord_2, dist_sq = position_check3D(info.particle_info[i].position, info.particle_info[j].position, boundary, interaction.cutoff)
        if iszero(dist_sq)
            nothing
        else
            temp = (interaction.σ)^2 / dist_sq
            dist = sqrt(dist_sq)
            lj_force = T(24) * interaction.ϵ * (T(2) * temp^T(6) - temp^T(3)) / dist
            direction = (coord_1 - coord_2) / dist
            lj_force_vec = lj_force * direction
            info.particle_info[i].acceleration += lj_force_vec / atoms[i].mass
            info.particle_info[j].acceleration -= lj_force_vec / atoms[j].mass
        end
    end
    return nothing
end

"""
    energy(interaction, neighborfinder, sys, info) -> T

The energy contributed by `interaction`. This is the interface every
energy-reporting `AbstractInteraction` implements (used by
[`EnergyLogger`](@ref)); not every interaction defines a method — one that
only ever appears inside a logger-free `MDSys` need not.
"""
function energy(interaction::LennardJones{T}, neighborfinder::T_NIEGHBER, sys::MDSys{T}, info::SimulationInfo{T}) where {T<:Number, T_NIEGHBER<:AbstractNeighborFinder}

    boundary = sys.boundary
    lj_energy = zero(T)
    for (i, j, r) in neighborfinder.neighbor_list
        coord_1, coord_2, dist_sq = position_check3D(info.particle_info[i].position, info.particle_info[j].position, boundary, interaction.cutoff)
        if iszero(dist_sq)
            nothing
        else
            temp = (interaction.σ)^2 / dist_sq
            lj_energy += T(4) * interaction.ϵ * (temp^T(6) - temp^T(3))
        end
    end
    return lj_energy
end