"""
    Point{N,T}(coo::NTuple{N,T})

A small fixed-size vector of `N` coordinates of type `T`, used throughout
ExTinyMD for positions, velocities, accelerations, and the electric field.
Supports `+`, `-`, scalar `*`/`/`, iteration, and indexing (`p[i]`).
"""
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

"""
    dist2(x, y) -> Number
    dist2(x) -> Number

Squared Euclidean distance between `x` and `y` (numbers or [`Point`](@ref)s), or
the squared norm of `x` alone.
"""
dist2(x::Number, y::Number) = abs2(x - y)
dist2(x::Point, y::Point) = sum(abs2, x - y)
dist2(x::Point) = sum(abs2, x)

"""
    Atom{T}(type, mass, charge)

A particle species: an integer `type` tag, a `mass`, and a `charge`. Construct
with keywords: `Atom(type = 1, mass = 1.0, charge = 0.0)`.
"""
struct Atom{T}
    type::Int
    mass::T
    charge::T
end

Atom(;type::Int = 1, mass::T = 1.0, charge::T = 0.0) where T<:Number = Atom{T}(type, mass, charge)

"""
    Boundary{T}(length, period)

The simulation box: `length` is the `(Lx, Ly, Lz)` edge lengths and `period` is
`(px, py, pz)` with `1` for a periodic axis and `0` for a free (non-periodic)
axis. See also [`CubicBoundary`](@ref) and [`Q2dBoundary`](@ref) for common
cases, and the `Boundary(L, ('p','p','f'))`-style constructor in
`boundary.jl`.
"""
struct Boundary{T}
    length::NTuple{3, T}
    period::NTuple{3, Int}
end

abstract type AbstractLogger end
abstract type AbstractNeighborFinder end
abstract type AbstractInteraction end
abstract type AbstractThermoStat end
abstract type AbstractSimulator end

"""
    MDSys(; n_atoms, atoms, boundary, interactions, loggers, simulator)

The top-level MD system: the atoms, the box, a list of
`(interaction, neighbor_finder)` pairs, the loggers to run each step, and the
simulator driving the integration. Passed to [`simulate!`](@ref) together with
a [`SimulationInfo`](@ref).
"""
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

"""
    SimulationInfo{T}

Mutable simulation state: the current `running_step`, a `Vector` of per-particle
`(id, position, velocity, acceleration)` records in storage-slot order, and
`id_dict` mapping a particle `id` to its current slot. Also constructible as
`SimulationInfo(n_atoms, atoms, place, boundary; min_r, max_attempts, rng, temp)`,
which randomly places `n_atoms` atoms inside the box `place =
(xlo, xhi, ylo, yhi, zlo, zhi)` and draws velocities from a Maxwell-Boltzmann
distribution at temperature `temp`.
"""
mutable struct SimulationInfo{T}
    running_step::Int64
    particle_info::Vector{PatricleInfo{T}}
    id_dict::Dict{Int, Int}
end

Base.show(io::IO, info::SimulationInfo) = print(io, "SimulationInfo: $(info.running_step) steps, $(length(info.particle_info)) particles")

"""
    NoThermoStat()

A thermostat that does nothing, i.e. NVE (microcanonical) dynamics.
"""
struct NoThermoStat <: AbstractThermoStat
    nostat::Bool
end
NoThermoStat() = NoThermoStat(true)

"""
    thermostat_update!(thermostat, sys, info)

Apply one thermostat step, typically rescaling or resampling velocities in
`info.particle_info`. Every concrete `AbstractThermoStat` implements a method;
[`NoThermoStat`](@ref) is a no-op.
"""
function thermostat_update!(thermostat::NoThermoStat, sys::MDSys{T}, info::SimulationInfo{T}) where T <: Number
    return nothing
end


"""
    AllNeighborFinder(n_atoms, T = Float64)

A neighbor finder whose list is every distinct pair `(i, j)`, `i < j`, with no
cutoff and no update logic. Useful for small systems or as a reference when a
cutoff-based finder is not needed.
"""
struct AllNeighborFinder{T} <: AbstractNeighborFinder
    neighbor_list::Vector{Tuple{Int64, Int64, T}}
end
AllNeighborFinder(n_atoms::TI, T::Type = Float64) where {TI <: Integer} = AllNeighborFinder{T}([(i, j, zero(T)) for i in 1:n_atoms - 1 for j in i+1:n_atoms])

Base.show(io::IO, neighborfinder::AllNeighborFinder) = print(io, "AllNeighborFinder with $(length(neighborfinder.neighbor_list)) pairs")

function update_finder!(neighborfinder::T_NIEGHBOR, info::SimulationInfo{T}) where {T<:Number, T_NIEGHBOR <: AllNeighborFinder}
    return nothing
end

"""
    NoNeighborFinder(T = Float64)

A neighbor finder with an always-empty list, for interactions such as
[`ExternalField`](@ref) that need no pairwise neighbours at all.
"""
struct NoNeighborFinder{T} <: AbstractNeighborFinder
    neighbor_list::Vector{Tuple{Int64, Int64, T}}
end
NoNeighborFinder(T::Type = Float64) = NoNeighborFinder{T}(Tuple{Int64, Int64, T}[])

Base.show(io::IO, neighborfinder::NoNeighborFinder) = print(io, "NoNeighborFinder")

"""
    update_finder!(neighborfinder, info)

Refresh `neighborfinder`'s neighbour list from the current positions in
`info`, typically every `update_steps` steps. Every concrete
`AbstractNeighborFinder` implements a method; [`NoNeighborFinder`](@ref) and
[`AllNeighborFinder`](@ref) are no-ops.
"""
function update_finder!(neighborfinder::T_NIEGHBOR, info::SimulationInfo{T}) where {T<:Number, T_NIEGHBOR <: NoNeighborFinder}
    return nothing
end

"""
    NoInteraction()

A placeholder interaction that contributes no force and no energy.
"""
struct NoInteraction <: AbstractInteraction
    nointeaction::Bool
end
NoInteraction() = NoInteraction(true)

Base.show(io::IO, interaction::NoInteraction) = print(io, "NoInteraction")

"""
    update_acceleration!(interaction, neighborfinder, sys, info)

Accumulate the acceleration contributed by `interaction` into
`info.particle_info[i].acceleration` for every particle `i`. This is the
interface every `AbstractInteraction` implements to plug into `simulate!`;
[`NoInteraction`](@ref) is a no-op. See `README.md` for how to add a new
interaction.
"""
function update_acceleration!(interaction::NoInteraction, neighborfinder::T_NEIGHBOR, sys::MDSys{T}, info::SimulationInfo{T}) where {T<:Number, T_NEIGHBOR<:AbstractNeighborFinder}
    return nothing
end