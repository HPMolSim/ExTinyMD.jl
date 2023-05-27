export LennardJones, acceleration

struct LennardJones{T} <: AbstractInteraction
    ϵ::T
    cutoff::T
    σ::T
end

LennardJones(;ϵ::T = 1.0, cutoff::T = 3.5, σ::T = 1.0) where T = LennardJones(ϵ, cutoff, σ)

function update_acceleration!(interaction::LennardJones{T}, neighborfinder::T_NIEGHBER, sys::MDSys{T}, info::SimulationInfo{T}) where {T<:Number, T_NIEGHBER<:AbstractNeighborFinder}
    atoms = sys.atoms
    boundary = sys.boundary
    update_finder!(neighborfinder, info)
    for (i, j, r) in neighborfinder.neighbor_list
        coord_1, coord_2, dist_sq = position_check3D(info.coords[i], info.coords[j], boundary, interaction.cutoff)
        if iszero(dist_sq)
            nothing
        else
            temp = (interaction.σ)^2 / dist_sq
            dist = sqrt(dist_sq)
            lj_force = T(24) * interaction.ϵ * (T(2) * temp^T(6) - temp^T(3)) / dist
            direction = (coord_1 - coord_2) / dist
            lj_force_vec = lj_force * direction
            info.acceleration[i] += lj_force_vec / atoms[i].mass
            info.acceleration[j] -= lj_force_vec / atoms[j].mass
        end
    end
    return nothing
end

function energy(interaction::LennardJones{T}, neighborfinder::T_NIEGHBER, sys::MDSys{T}, info::SimulationInfo{T}) where {T<:Number, T_NIEGHBER<:AbstractNeighborFinder}

    boundary = sys.boundary
    lj_energy = zero(T)
    for (i, j, r) in neighborfinder.neighbor_list
        coord_1, coord_2, dist_sq = position_check3D(info.coords[i], info.coords[j], boundary, interaction.cutoff)
        if iszero(dist_sq)
            nothing
        else
            temp = (interaction.σ)^2 / dist_sq
            lj_energy += T(4) * interaction.ϵ * (temp^T(6) - temp^T(3))
        end
    end
    return lj_energy
end