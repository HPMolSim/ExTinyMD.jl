# The `cell_list` field holds a `CellListMap.InPlaceNeighborList`, whose concrete type
# depends on CellListMap internals, so it is kept as a free type parameter `NL`.
"""
    CellList3D{T, TI, NL} <: AbstractNeighborFinder

A 3D neighbour finder backed by `CellListMap.InPlaceNeighborList`, rebuilt in
place every `update_steps` steps.

    CellList3D(info, cutoff, boundary, update_steps)

Build one from the current positions in `info` with the given `cutoff`. A
non-periodic axis of `boundary` has its cell-list unit-cell length inflated so
no periodic images are found along it.
"""
mutable struct CellList3D{T, TI, NL} <: AbstractNeighborFinder
    cell_list::NL
    neighbor_list::Vector{Tuple{Int64, Int64, T}}
    update_steps::TI
end

Base.show(io::IO, cell_list::CellList3D) = print(io, "CellList3D")

function CellList3D(info::SimulationInfo{T}, cutoff::T, boundary::Boundary{T}, update_steps::TI) where {T<:Number, TI<:Integer}
    coords = [SVector{3, T}(p_info.position[1], p_info.position[2], p_info.position[3]) for p_info in info.particle_info]

    # if the system is non-periodic in some direction, set the unitcell length at that direction as 2 L_max so that no periodic images will be counted
    unitcell = [isone(boundary.period[i]) ? boundary.length[i] : T(1.5) * maximum(boundary.length) for i in 1:3]
    cell_list = InPlaceNeighborList(xpositions = coords, cutoff = cutoff, unitcell = unitcell, parallel=true)
    update!(cell_list, xpositions = coords)
    neighbor_list = neighborlist!(cell_list)

    return CellList3D{T, TI, typeof(cell_list)}(cell_list, neighbor_list, update_steps)
end

function update_finder!(neighborfinder::T_NIEGHBOR, info::SimulationInfo{T}) where {T<:Number, T_NIEGHBOR <: CellList3D}
    if iszero(info.running_step % neighborfinder.update_steps)
        coords = [SVector{3, T}(p_info.position[1], p_info.position[2], p_info.position[3]) for p_info in info.particle_info]
        update!(neighborfinder.cell_list, xpositions = coords)
        neighborfinder.neighbor_list = neighborlist!(neighborfinder.cell_list)
    end
    return nothing
end

"""
    CellListQ2D{T, TI, NL} <: AbstractNeighborFinder

Like [`CellList3D`](@ref), but the neighbour search is only in x and y (z is
ignored), for quasi-2D slabs.

    CellListQ2D(info, cutoff, boundary, update_steps)
"""
mutable struct CellListQ2D{T, TI, NL} <: AbstractNeighborFinder
    cell_list::NL
    neighbor_list::Vector{Tuple{Int64, Int64, T}}
    update_steps::TI
end

Base.show(io::IO, cell_list::CellListQ2D) = print(io, "CellListQ2D")

function CellListQ2D(info::SimulationInfo{T}, cutoff::T, boundary::Boundary{T}, update_steps::TI) where {T<:Number, TI<:Integer}
    coords = [SVector{2, T}(p_info.position[1], p_info.position[2]) for p_info in info.particle_info]

    # if the system is non-periodic in some direction, set the unitcell length at that direction as 2 L_max so that no periodic images will be counted
    unitcell = SVector{2, T}([isone(boundary.period[i]) ? boundary.length[i] : T(1.5) * maximum(boundary.length) for i in 1:2])
    cell_list = InPlaceNeighborList(xpositions = coords, cutoff = cutoff, unitcell = unitcell, parallel=true)
    neighbor_list = neighborlist!(cell_list)

    return CellListQ2D{T, TI, typeof(cell_list)}(cell_list, neighbor_list, update_steps)
end

function update_finder!(neighborfinder::T_NIEGHBOR, info::SimulationInfo{T}) where {T<:Number, T_NIEGHBOR <: CellListQ2D}
    
    if iszero(info.running_step % neighborfinder.update_steps)
        coords = [SVector{2, T}(p_info.position[1], p_info.position[2]) for p_info in info.particle_info]
        update!(neighborfinder.cell_list, xpositions = coords)
        neighborfinder.neighbor_list = neighborlist!(neighborfinder.cell_list)
    end
    return nothing
end

"""
    CellListDir3D{T, TI} <: AbstractNeighborFinder

A 3D neighbour finder that recomputes the full neighbour list directly (via
`CellListMap.neighborlist`, not an in-place cell list) every `update_steps`
steps. Simpler and more predictable than [`CellList3D`](@ref) at the cost of
reallocating each rebuild.

    CellListDir3D(info, cutoff, boundary, update_steps)
"""
mutable struct CellListDir3D{T, TI} <: AbstractNeighborFinder
    unitcell::SVector{3, T}
    cutoff::T
    neighbor_list::Vector{Tuple{Int64, Int64, T}}
    update_steps::TI
end

function CellListDir3D(info::SimulationInfo{T}, cutoff::T, boundary::Boundary{T}, update_steps::TI) where {T<:Number, TI<:Integer}
    coords = [SVector{3, T}(p_info.position[1], p_info.position[2], p_info.position[3]) for p_info in info.particle_info]

    # if the system is non-periodic in some direction, set the unitcell length at that direction as 2 L_max so that no periodic images will be counted
    unitcell = SVector{3, T}([isone(boundary.period[i]) ? boundary.length[i] : T(1.5) * maximum(boundary.length) for i in 1:3])
    neighbor_list = neighborlist(xpositions = coords, cutoff = cutoff, unitcell = unitcell, parallel = false)

    return CellListDir3D{T, TI}(unitcell, cutoff, neighbor_list, update_steps)
end

function update_finder!(neighborfinder::CellListDir3D{T, TI}, info::SimulationInfo{T}) where {T<:Number, TI<:Integer}
    if iszero(info.running_step % neighborfinder.update_steps)
        coords = [SVector{3, T}(p_info.position[1], p_info.position[2], p_info.position[3]) for p_info in info.particle_info]
        neighborfinder.neighbor_list = neighborlist(xpositions = coords, cutoff = neighborfinder.cutoff, unitcell = neighborfinder.unitcell, parallel = false)
    end
    return nothing
end

"""
    CellListDirQ2D{T, TI} <: AbstractNeighborFinder

Like [`CellListDir3D`](@ref), but the neighbour search is only in x and y, for
quasi-2D slabs.

    CellListDirQ2D(info, cutoff, boundary, update_steps)
"""
mutable struct CellListDirQ2D{T, TI} <: AbstractNeighborFinder
    unitcell::SVector{2, T}
    cutoff::T
    neighbor_list::Vector{Tuple{Int64, Int64, T}}
    update_steps::TI
end

function CellListDirQ2D(info::SimulationInfo{T}, cutoff::T, boundary::Boundary{T}, update_steps::TI) where {T<:Number, TI<:Integer}
    coords = [SVector{2, T}(p_info.position[1], p_info.position[2]) for p_info in info.particle_info]

    # if the system is non-periodic in some direction, set the unitcell length at that direction as 2 L_max so that no periodic images will be counted
    unitcell = SVector{2, T}([isone(boundary.period[i]) ? boundary.length[i] : T(1.5) * maximum(boundary.length) for i in 1:2])
    neighbor_list = neighborlist(xpositions = coords, cutoff = cutoff, unitcell = unitcell, parallel = true)

    return CellListDirQ2D{T, TI}(unitcell, cutoff, neighbor_list, update_steps)
end

function update_finder!(neighborfinder::CellListDirQ2D{T, TI}, info::SimulationInfo{T}) where {T<:Number, TI<:Integer}
    if iszero(info.running_step % neighborfinder.update_steps)
        coords = [SVector{2, T}(p_info.position[1], p_info.position[2]) for p_info in info.particle_info]
        neighborfinder.neighbor_list = neighborlist(xpositions = coords, cutoff = neighborfinder.cutoff, unitcell = neighborfinder.unitcell, parallel = true)
    end
    return nothing
end
