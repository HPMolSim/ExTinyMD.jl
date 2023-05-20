export Point, Atom, Boundary, Q2dBoudary, MDsys

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

function BoundaryCheck(coo::Point{3, T}, boundary::Boundary{T}) where T <: Number
    return Point(
        coo[1] - boundary.period[1] * boundary.length[1] * div(coo[1], boundary.length[1], RoundDown), 
        coo[2] - boundary.period[2] * boundary.length[2] * div(coo[2], boundary.length[2], RoundDown), 
        coo[3] - boundary.period[3] * boundary.length[3] * div(coo[3], boundary.length[3], RoundDown))
end

abstract type AbstractLogger end
abstract type AbstractNeighborFinder end
abstract type AbstractInteraction end
abstract type AbstractThermoStat end
abstract type AbstractSimulator end

struct MDSys{T}
    n_atoms::Integer
    atoms::Vector{Atom{T}}
    boundary::Boundary{T}
    interactions::Tuple{AbstractInteraction}
    logger::Tuple{AbstractLogger}
    simulator::AbstractSimulator
end

function MDSys(;
    n_atoms::TI,
    atoms::Vector{Atom{T}},
    boundary::Boundary{T},
    interactions::Tuple{AbstractInteraction},
    thermostat::Tuple{AbstractThermoStat},
    logger::Tuple{AbstractLogger},
    simulator::AbstractSimulator,
) where {TI <: Integer, T}
    return System{T}(n_atoms, atoms, boundary, interactions, thermostat, logger, simulator)
end

mutable struct SimulationInfo{T}
    running_step::Int64
    coords::Vector{Point{3, T}}
    velcoity::Vector{Point{3, T}}
end

function SimulationInfo(mdsys::MDSys, place::NTuple{6, T}; min_r=zero(T), max_attempts::TI=100, rng=Random.GLOBAL_RNG, temp::T = 1.0) where {T, TI<:Integer}
    atoms_coords = random_position(mdsys.n_atoms, place, mdsys.boundary; min_r = min_r, max_attempts = max_attempts)
    atoms_velocity = random_velocity(;temp = temp, atoms = mdsys.atoms, rng=rng)

    return SimulationInfo{T}(zero(Int64), atoms_coords, atoms_velocity)
end

struct NoThermoStat <: AbstractThermoStat
    nostat::Bool
end

NoThermoStat() = NoThermoStat(true)