struct Point{N,T}
    coo::NTuple{N,T}
end
Point(arg::T, args::T...) where T<:Number = Point((arg, args...))
Base.:(+)(x::Point{N,T}, y::Point{N,T}) where {N, T} = Point(x.coo .+ y.coo)
Base.:(-)(x::Point{N,T}, y::Point{N,T}) where {N, T} = Point(x.coo .- y.coo)
Base.:(-)(x::Point{N,T}) where {N, T} = Point(Base.:(-).(x.coo))
Base.adjoint(x::Point) = x
Base.:(*)(x::Number, y::Point) = Point(y.coo .* x)
Base.:(*)(y::Point, x::Number) = Point(y.coo .* x)
Base.:(/)(y::Point, x::Number) = Point(y.coo ./ x)
Base.iterate(x::Point, args...) = Base.iterate(x.coo, args...)
Base.getindex(x::Point, i::Int) = x.coo[i]

dist2(x::Number, y::Number) = abs2(x - y)
dist2(x::Point, y::Point) = sum(abs2, x - y)
dist2(x::Point) = sum(abs2, x)

struct Atom{T}
    type::Int
    mass::T
    charge::T
end

Atom(;type::Int = 1, mass::T = 1.0, charge::T = 0.0) where T<:Number = Atom{T}(type, mass, charge)

struct Boundary{T}
    length::NTuple{3, T}
    period::NTuple{3, Int}
end

abstract type AbstractLogger end
abstract type AbstractNeighborFinder end
abstract type AbstractInteraction end
abstract type AbstractThermoStat end
abstract type AbstractSimulator end

struct MDSys{T_NUM, T_INTERACTION, T_LOGGER, T_SIMULATOR}
    n_atoms::Int64
    atoms::Vector{Atom{T_NUM}}
    boundary::Boundary{T_NUM}
    interactions::Vector{T_INTERACTION}
    loggers::Vector{T_LOGGER}
    simulator::T_SIMULATOR
end

Base.show(io::IO, sys::MDSys) = print(io, " MDSys with $(sys.n_atoms) atoms \n boundary: $(sys.boundary) \n simulator: $(sys.simulator) \n interactions: $(sys.interactions) \n loggers: $(sys.loggers)")

function MDSys(;
    n_atoms::Int64,
    atoms::Vector{Atom{T_NUM}},
    boundary::Boundary{T_NUM},
    interactions::Vector{T_INTERACTION},
    loggers::Vector{T_LOGGER},
    simulator::T_SIMULATOR,
) where {T_NUM <: Number, T_INTERACTION <: Tuple{AbstractInteraction, AbstractNeighborFinder}, T_LOGGER <: AbstractLogger, T_SIMULATOR <: AbstractSimulator}
    return MDSys{T_NUM, T_INTERACTION, T_LOGGER, T_SIMULATOR}(n_atoms, atoms, boundary, interactions, loggers, simulator)
end

mutable struct PatricleInfo{T}
    id::Int
    position::Point{3, T}
    velocity::Point{3, T}
    acceleration::Point{3, T}
end

mutable struct SimulationInfo{T}
    running_step::Int64
    particle_info::Vector{PatricleInfo{T}}
    id_dict::Dict{Int, Int}
end

Base.show(io::IO, info::SimulationInfo) = print(io, "SimulationInfo: $(info.running_step) steps, $(length(info.particle_info)) particles")

struct NoThermoStat <: AbstractThermoStat
    nostat::Bool
end
NoThermoStat() = NoThermoStat(true)

function thermostat_update!(thermostat::NoThermoStat, sys::MDSys{T}, info::SimulationInfo{T}) where T <: Number
    return nothing
end


struct AllNeighborFinder{T} <: AbstractNeighborFinder
    neighborlist::Vector{Tuple{Int64, Int64, T}}
end
AllNeighborFinder(n_atoms::TI, T::Type = Float64) where {TI <: Integer} = AllNeighborFinder{T}([(i, j, zero(T)) for i in 1:n_atoms - 1 for j in i+1:n_atoms])

Base.show(io::IO, neighborfinder::AllNeighborFinder) = print(io, "AllNeighborFinder with $(length(neighborfinder.neighborlist)) pairs")

function update_finder!(neighborfinder::T_NIEGHBOR, info::SimulationInfo{T}) where {T<:Number, T_NIEGHBOR <: AllNeighborFinder}
    return nothing
end

struct NoNeighborFinder{T} <: AbstractNeighborFinder
    neighborlist::Vector{Tuple{Int64, Int64, T}}
end
NoNeighborFinder(T::Type = Float64) = NoNeighborFinder{T}([(0, 0, zero(T))])

Base.show(io::IO, neighborfinder::NoNeighborFinder) = print(io, "NoNeighborFinder")

function update_finder!(neighborfinder::T_NIEGHBOR, info::SimulationInfo{T}) where {T<:Number, T_NIEGHBOR <: NoNeighborFinder}
    return nothing
end

struct NoInteraction <: AbstractInteraction
    nointeaction::Bool
end
NoInteraction() = NoInteraction(true)

Base.show(io::IO, interaction::NoInteraction) = print(io, "NoInteraction")

function update_acceleration!(interaction::NoInteraction, neighborfinder::T_NEIGHBOR, sys::MDSys{T}, info::SimulationInfo{T}) where {T<:Number, T_NEIGHBOR<:AbstractNeighborFinder}
    return nothing
end

