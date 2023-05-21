export PointToStaticArrays3D, PointToStaticArrays2D

mutable struct CellList3D{T} <: AbstractNeighborFinder
    celllist::InPlaceNeighborList
end



function PointToStaticArrays3D(x::Vector{Point{3, T}}) where T
    return [SVector{3, T}(xi[1], xi[2], xi[3]) for xi in x]
end

function PointToStaticArrays2D(x::Vector{Point{3, T}}) where T
    return [SVector{2, T}(xi[1], xi[2]) for xi in x]
end