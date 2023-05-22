export CellList3D, CellList2D, PointToStaticArrays3D!, PointToStaticArrays2D!

mutable struct CellList3D{T, TI} <: AbstractNeighborFinder
    cell_list::InPlaceNeighborList{Box{OrthorhombicCell, 3, T, T, 9, T}, CellList{3, T}, CellListMap.AuxThreaded{3, T}, CellListMap.NeighborList{T}}
    neighbor_list::Vector{Tuple{Int64, Int64, T}}
    coords::Vector{SVector{3, T}}
    update_steps::TI
end

function CellList3D(info::SimulationInfo{T}, cutoff::T, boundary::Boundary{T}, update_steps::TI) where {T<:Number, TI<:Integer}
    coords = [SVector{3, T}(xi[1], xi[2], xi[3]) for xi in info.coords]

    # if the system is non-periodic in some direction, set the unitcell length at that direction as 2 L_max so that no periodic images will be counted
    unitcell = [isone(boundary.period[i]) ? boundary.length[i] : T(1.5) * maximum(boundary.length) for i in 1:3]
    cell_list = InPlaceNeighborList(x = coords, cutoff = cutoff, unitcell = unitcell, parallel=false)
    update!(cell_list, coords)
    neighbor_list = neighborlist!(cell_list)

    return CellList3D{T, TI}(cell_list, neighbor_list, coords, update_steps)
end

function update_finder!(neighborfinder::T_NIEGHBOR, info::SimulationInfo{T}) where {T<:Number, T_NIEGHBOR <: CellList3D}
    if isone(info.running_step % neighborfinder.update_steps)
        PointToStaticArrays3D!(info.coords, neighborfinder.coords)
        update!(neighborfinder.cell_list, neighborfinder.coords)
        neighborfinder.neighbor_list = neighborlist!(neighborfinder.cell_list)
    end
    return nothing
end


function PointToStaticArrays3D!(x_new::Vector{Point{3, T}}, x_old::Vector{SVector{3, T}}) where T
    @inbounds for i in 1:length(x_new)
        x_old[i] = SVector{3, T}(x_new[i][1], x_new[i][2], x_new[i][3])
    end
    return nothing
end

mutable struct CellListQ2D{T, TI} <: AbstractNeighborFinder
    cell_list::InPlaceNeighborList{Box{OrthorhombicCell, 2, T, T, 4, T}, CellList{2, T}, CellListMap.AuxThreaded{2, T}, CellListMap.NeighborList{T}}
    neighbor_list::Vector{Tuple{Int64, Int64, T}}
    coords::Vector{SVector{2, T}}
    update_steps::TI
end

function CellListQ2D(info::SimulationInfo{T}, cutoff::T, boundary::Boundary{T}, update_steps::TI) where {T<:Number, TI<:Integer}
    coords = [SVector{2, T}(xi[1], xi[2]) for xi in info.coords]

    # if the system is non-periodic in some direction, set the unitcell length at that direction as 2 L_max so that no periodic images will be counted
    unitcell = [boundary.length[1], boundary.length[2]]
    cell_list = InPlaceNeighborList(x = coords, cutoff = cutoff, unitcell = unitcell, parallel=false)
    update!(cell_list, coords)
    neighbor_list = neighborlist!(cell_list)

    return CellListQ2D{T, TI}(cell_list, neighbor_list, coords, update_steps)
end

function update_finder!(neighborfinder::T_NIEGHBOR, info::SimulationInfo{T}) where {T<:Number, T_NIEGHBOR <: CellListQ2D}

    PointToStaticArrays2D!(info.coords, neighborfinder.coords)
    update!(neighborfinder.cell_list, neighborfinder.coords)
    neighborfinder.neighbor_list = neighborlist!(neighborfinder.cell_list)
    return nothing
end

function PointToStaticArrays2D!(x_new::Vector{Point{3, T}}, x_old::Vector{SVector{2, T}}) where T
    @inbounds for i in 1:length(x_new)
        x_old[i] = SVector{2, T}(x_new[i][1], x_new[i][2])
    end
    return nothing
end