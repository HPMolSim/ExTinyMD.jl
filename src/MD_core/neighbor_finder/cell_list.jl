mutable struct CellList3D{T, TI} <: AbstractNeighborFinder
    cell_list::InPlaceNeighborList{Box{OrthorhombicCell, 3, T, T, 9, T}, CellList{3, T}, CellListMap.AuxThreaded{3, T}, CellListMap.NeighborList{T}}
    neighbor_list::Vector{Tuple{Int64, Int64, T}}
    update_steps::TI
end

Base.show(io::IO, cell_list::CellList3D) = print(io, "CellList3D")

function CellList3D(info::SimulationInfo{T}, cutoff::T, boundary::Boundary{T}, update_steps::TI) where {T<:Number, TI<:Integer}
    coords = [SVector{3, T}(p_info.position[1], p_info.position[2], p_info.position[3]) for p_info in info.particle_info]

    # if the system is non-periodic in some direction, set the unitcell length at that direction as 2 L_max so that no periodic images will be counted
    unitcell = [isone(boundary.period[i]) ? boundary.length[i] : T(1.5) * maximum(boundary.length) for i in 1:3]
    cell_list = InPlaceNeighborList(x = coords, cutoff = cutoff, unitcell = unitcell, parallel=true)
    update!(cell_list, coords)
    neighbor_list = neighborlist!(cell_list)

    return CellList3D{T, TI}(cell_list, neighbor_list, update_steps)
end

function update_finder!(neighborfinder::T_NIEGHBOR, info::SimulationInfo{T}) where {T<:Number, T_NIEGHBOR <: CellList3D}
    if iszero(info.running_step % neighborfinder.update_steps)
        coords = [SVector{3, T}(p_info.position[1], p_info.position[2], p_info.position[3]) for p_info in info.particle_info]
        update!(neighborfinder.cell_list, coords)
        neighborfinder.neighbor_list = neighborlist!(neighborfinder.cell_list)
    end
    return nothing
end

mutable struct CellListQ2D{T, TI} <: AbstractNeighborFinder
    cell_list::InPlaceNeighborList{Box{OrthorhombicCell, 2, T, T, 4, T}, CellList{2, T}, CellListMap.AuxThreaded{2, T}, CellListMap.NeighborList{T}}
    neighbor_list::Vector{Tuple{Int64, Int64, T}}
    update_steps::TI
end

Base.show(io::IO, cell_list::CellListQ2D) = print(io, "CellListQ2D")

function CellListQ2D(info::SimulationInfo{T}, cutoff::T, boundary::Boundary{T}, update_steps::TI) where {T<:Number, TI<:Integer}
    coords = [SVector{2, T}(p_info.position[1], p_info.position[2]) for p_info in info.particle_info]

    # if the system is non-periodic in some direction, set the unitcell length at that direction as 2 L_max so that no periodic images will be counted
    unitcell = SVector{2, T}([isone(boundary.period[i]) ? boundary.length[i] : T(1.5) * maximum(boundary.length) for i in 1:2])
    cell_list = InPlaceNeighborList(x = coords, cutoff = cutoff, unitcell = unitcell, parallel=true)
    neighbor_list = neighborlist!(cell_list)

    return CellListQ2D{T, TI}(cell_list, neighbor_list, update_steps)
end

function update_finder!(neighborfinder::T_NIEGHBOR, info::SimulationInfo{T}) where {T<:Number, T_NIEGHBOR <: CellListQ2D}
    
    if iszero(info.running_step % neighborfinder.update_steps)
        coords = [SVector{2, T}(p_info.position[1], p_info.position[2]) for p_info in info.particle_info]
        update!(neighborfinder.cell_list, coords)
        neighborfinder.neighbor_list = neighborlist!(neighborfinder.cell_list)
    end
    return nothing
end

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
    neighbor_list = neighborlist(coords, cutoff; unitcell = unitcell, parallel = false)

    return CellListDir3D{T, TI}(unitcell, cutoff, neighbor_list, update_steps)
end

function update_finder!(neighborfinder::CellListDir3D{T, TI}, info::SimulationInfo{T}) where {T<:Number, TI<:Integer}
    if iszero(info.running_step % neighborfinder.update_steps)
        coords = [SVector{3, T}(p_info.position[1], p_info.position[2], p_info.position[3]) for p_info in info.particle_info]
        neighborfinder.neighbor_list = neighborlist(coords, neighborfinder.cutoff; unitcell = neighborfinder.unitcell, parallel = false)
    end
    return nothing
end

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
    neighbor_list = neighborlist(coords, cutoff; unitcell = unitcell, parallel = true)

    return CellListDirQ2D{T, TI}(unitcell, cutoff, neighbor_list, update_steps)
end

function update_finder!(neighborfinder::CellListDirQ2D{T, TI}, info::SimulationInfo{T}) where {T<:Number, TI<:Integer}
    if iszero(info.running_step % neighborfinder.update_steps)
        coords = [SVector{2, T}(p_info.position[1], p_info.position[2]) for p_info in info.particle_info]
        neighborfinder.neighbor_list = neighborlist(coords, neighborfinder.cutoff; unitcell = neighborfinder.unitcell, parallel = true)
    end
    return nothing
end