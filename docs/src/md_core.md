```@meta
CurrentModule = ExTinyMD
```

# MD Core

This page covers the pieces used to build and drive an ordinary (non-electrostatic)
molecular dynamics simulation: the box, the atoms, the simulation state, the
integrators and thermostats, the loggers, and the neighbour finders. See
[Interactions](@ref) for how forces get added to the system, and
[Electrostatics](@ref) for the electrostatics standard library.

## Box and atoms

A simulation lives in a [`Boundary`](@ref) box populated with [`Atom`](@ref)s.
[`CubicBoundary`](@ref) and [`Q2dBoundary`](@ref) build the two shapes used
throughout this package: a fully periodic cube, and a slab periodic in x and y
only.

```@docs
Boundary
CubicBoundary
Q2dBoundary
Atom
create_atoms
```

## Simulation state

[`SimulationInfo`](@ref) holds the mutable state — positions, velocities,
accelerations, and the running step count — while [`MDSys`](@ref) bundles the
(immutable, for the run) description of the system: the atoms, the box, the
interactions, the loggers, and the simulator.

```@docs
SimulationInfo
random_position
random_velocity
MDSys
BoundaryCheck!
position_check3D
position_checkQ2D
dist2
Point
```

## Running a simulation

```@docs
simulate!
VerletProcess
NHVerletProcess
```

## Thermostats

```@docs
NoThermoStat
AndersenThermoStat
BerendsenThermoStat
thermostat_update!
```

## Loggers

Loggers are attached to `MDSys` via the `loggers` keyword and run once per step
inside [`simulate!`](@ref).

```@docs
TemperatureLogger
TrajectoryLogger
EnergyLogger
```

## Neighbour finders

A neighbour finder is paired with each interaction as `(interaction, finder)`
in `MDSys`'s `interactions` list. [`update_finder!`](@ref) refreshes it; most
interactions call this themselves before using the list.

```@docs
update_finder!
NoNeighborFinder
AllNeighborFinder
CellList3D
CellListDir3D
CellListQ2D
CellListDirQ2D
SubNeighborFinder
```

## Data loading and analysis

```@docs
load_trajectory
load_lammpstrj
data2info
z_hist
hist_init
distance_hist!
MSD
```
