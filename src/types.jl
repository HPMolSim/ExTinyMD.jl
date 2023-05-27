export Point, Atom, Boundary, Q2dBoudary, CubicBoundary, MDSys, position_check3D, position_checkQ2D, BoundaryCheck!, SimulationInfo, thermostat_update!, update_acceleration!, update_finder!, NoInteraction, NoNeighborFinder, NoThermoStat, dist2

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


struct Atom{T}
    mass::T
    charge::T
end

Atom(;mass::T, charge::T) where T<:Number = Atom{T}(mass, charge)

struct Boundary{T}
    length::NTuple{3, T}
    period::NTuple{3, Int}
end

function Boundary(L::NTuple{3, T}, period_set::NTuple{3, Char}) where T
    if period_set == ('p', 'p', 'p')
        return Boundary(L, (1, 1, 1))
    elseif period_set == ('p', 'p', 'f')
        return Boundary(L, (1, 1, 0))
    elseif period_set == ('p', 'f', 'f')
        return Boundary(L, (1, 0, 0))
    end
    error("Illegale Input!")
end

function Q2dBoudary(Lx::T, Ly::T, Lz::T) where T
    return Boundary((Lx, Ly, Lz), (1, 1, 0))
end

function CubicBoundary(L::T) where T
    return Boundary((L, L, L), (1, 1, 1))
end


function BoundaryCheck!(coords::Vector{Point{3, T}}, boundary::Boundary{T}) where T <: Number
    Lx, Ly, Lz = boundary.length
    px, py, pz = boundary.period
    for i in 1:length(coords)
        coords[i] -= Point(
            ((coords[i][1]>zero(T) && coords[i][1]<Lx) || iszero(px)) ? zero(T) : Lx * div(coords[i][1], Lx, RoundDown),
            ((coords[i][2]>zero(T) && coords[i][2]<Ly) || iszero(py)) ? zero(T) : Ly * div(coords[i][2], Ly, RoundDown),
            ((coords[i][3]>zero(T) && coords[i][3]<Lz) || iszero(pz)) ? zero(T) : Lz * div(coords[i][3], Lz, RoundDown)
        )
    end
    return nothing
end


@inbounds function position_check3D(coord_1::Point{3, T}, coord_2::Point{3, T}, boundary::Boundary{T}, cutoff::T) where T

    for mx = -boundary.period[1]:boundary.period[1]
        for my = -boundary.period[2]:boundary.period[2]
            for mz = -boundary.period[3]:boundary.period[3]
                dist_sq = dist2(coord_1 + Point(mx * boundary.length[1], my * boundary.length[2], mz * boundary.length[3]), coord_2)
                if dist_sq < abs2(cutoff)
                    return (coord_1 + Point(mx * boundary.length[1], my * boundary.length[2], mz * boundary.length[3]), coord_2, dist_sq)
                end
            end
        end
    end
    return (Point(zero(T), zero(T), zero(T)), Point(zero(T), zero(T), zero(T)), zero(T))
end

@inbounds function position_checkQ2D(coord_1::Point{3, T}, coord_2::Point{3, T}, boundary::Boundary{T}, cutoff::T) where T

    for mx = -boundary.period[1]:boundary.period[1]
        for my = -boundary.period[2]:boundary.period[2]
            dist_sq = abs2(coord_1[1] + mx * boundary.length[1] - coord_2[1]) + abs2(coord_1[2] + mx * boundary.length[2] - coord_2[2])
            if dist_sq < abs2(cutoff)
                return (coord_1 + Point(mx * boundary.length[1], my * boundary.length[2], zero(T)), coord_2, dist_sq)
            end
        end
    end
    return (Point(zero(T), zero(T), zero(T)), Point(zero(T), zero(T), zero(T)), zero(T))
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

mutable struct SimulationInfo{T}
    running_step::Int64
    coords::Vector{Point{3, T}}
    velcoity::Vector{Point{3, T}}
    acceleration::Vector{Point{3, T}}
end

function SimulationInfo(n_atoms::TI, atoms::Vector{Atom{T}}, place::NTuple{6, T}, boundary::Boundary{T}; min_r=zero(T), max_attempts::TI=100, rng=Random.GLOBAL_RNG, temp::T = 1.0) where {T, TI<:Integer}
    atoms_coords = random_position(n_atoms, place, boundary; min_r = min_r, max_attempts = max_attempts)
    atoms_velocity = random_velocity(;temp = temp, atoms = atoms, rng=rng)
    acceleration = [Point(zero(T), zero(T), zero(T)) for _=1:n_atoms]
    return SimulationInfo{T}(zero(Int64), atoms_coords, atoms_velocity, acceleration)
end

struct NoThermoStat <: AbstractThermoStat
    nostat::Bool
end

NoThermoStat() = NoThermoStat(true)

function thermostat_update!(thermostat::NoThermoStat, sys::MDSys{T}, info::SimulationInfo{T}) where T <: Number
    return nothing
end


struct NoNeighborFinder{T} <: AbstractNeighborFinder
    neighborlist::Vector{Tuple{Int64, Int64, T}}
end

NoNeighborFinder(n_atoms::TI, T::Type = Float64) where {TI <: Integer} = NoNeighborFinder{T}([(i, j, zero(T)) for i in 1:n_atoms - 1 for j in i+1:n_atoms])

function update_finder!(neighborfinder::T_NIEGHBOR, info::SimulationInfo{T}) where {T<:Number, T_NIEGHBOR <: NoNeighborFinder}
    return nothing
end



struct NoInteraction <: AbstractInteraction
    nointeaction::Bool
end

NoInteraction() = NoInteraction(true)

function update_acceleration!(interaction::NoInteraction, neighborfinder::T_NEIGHBOR, atoms::Vector{Atom{T}}, boundary::Boundary{T}, info::SimulationInfo{T}) where {T<:Number, T_NEIGHBOR<:AbstractNeighborFinder}
    return nothing
end

