"""
    ExternalField(; E = Point((0.0, 0.0, 0.0)))

A uniform external electric field `E`; each particle feels a force `charge * E`.
Pair with a [`NoNeighborFinder`](@ref), since no neighbour list is needed.
"""
struct ExternalField{T} <: AbstractInteraction
    E::Point{3, T}
end

Base.show(io::IO, interaction::ExternalField) = print(io, "ExternalField with E = $(interaction.E)")

function ExternalField(;E::Point{3, T} = Point((0.0, 0.0, 0.0))) where {T} 
    return ExternalField{T}(E)
end

function update_acceleration!(interaction::ExternalField{T}, neighborfinder::NoNeighborFinder{T}, sys::MDSys{T}, info::SimulationInfo{T}) where {T<:Number}
    atoms = sys.atoms
    for i in 1:length(info.particle_info)
        info.particle_info[i].acceleration += atoms[i].charge * interaction.E / atoms[i].mass
    end
    return nothing
end