export SubLennardJones, update_acceleration!

struct SubLennardJones{T} <: AbstractInteraction
    ϵ::T
    cutoff::T
    σ::T
    sub_down::T
    sub_up::T
end

SubLennardJones(sub_down::T, sub_up::T;ϵ::T = 1.0, cutoff::T = 3.5, σ::T = 1.0) where T<:Number = SubLennardJones(ϵ, cutoff, σ, sub_down, sub_up)

function update_acceleration!(interaction::SubLennardJones{T}, neighborfinder::SubNeighborFinder{T, TI}, atoms::Vector{Atom{T}}, boundary::Boundary{T}, info::SimulationInfo{T}) where {T<:Number, TI<:Integer}
    for i in neighborfinder.down_neighbor
        dz = info.coords[i][3] - interaction.sub_down
        if zero(T) < dz < interaction.cutoff
            temp = (interaction.σ)^2 / dz^2
            lj_force = T(24) * interaction.ϵ * (T(2) * temp^T(6) - temp^T(3)) / dz
            mass = atoms[i].mass
            info.acceleration[i] += Point(zero(T), zero(T), lj_force / mass)
        end
    end

    for i in neighborfinder.up_neighbor
        dz = interaction.sub_up - info.coords[i][3]
        if zero(T) < dz < interaction.cutoff
            temp = (interaction.σ)^2 / dz^2
            lj_force = - T(24) * interaction.ϵ * (T(2) * temp^T(6) - temp^T(3)) / dz
            mass = atoms[i].mass
            info.acceleration[i] += Point(zero(T), zero(T), lj_force / mass)
        end
    end
    return nothing
end

