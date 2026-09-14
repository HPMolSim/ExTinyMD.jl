```@meta
CurrentModule = ExTinyMD
```

# ExTinyMD

`ExTinyMD.jl` ("Extremely Tiny Molecular Dynamics") is a small, hackable
molecular dynamics package written in Julia. It provides the pieces of an MD
code — a box and atoms, an integrator, thermostats, neighbour finders,
loggers — plus a growing set of interactions, including a full electrostatics
standard library (see [Electrostatics](@ref)).

It is intentionally minimal: adding a new interaction is a couple of methods
(see [Interactions](@ref)), and every piece is a plain Julia struct with no
hidden global state.

## Installation

```julia
julia> ]
pkg> add ExTinyMD
```

## A first simulation: a 3D Lennard-Jones fluid

The following walkthrough (adapted from `example/LJ_fluid_rdf.jl`) sets up a
periodic cube of Lennard-Jones particles, runs it under a Verlet integrator
with an Andersen thermostat, and samples a radial distribution function (rdf).

First, build the box and the atoms:

```julia
using ExTinyMD

n_atoms = 1000
L = 100.0
boundary = CubicBoundary(L)
atoms = create_atoms([(n_atoms, Atom(type = 1, mass = 1.0, charge = 0.0))])
```

Randomly place them and give them thermal velocities:

```julia
info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, L), boundary;
                      min_r = 0.1, temp = 1.0)
```

Attach a Lennard-Jones interaction, using a cell list to find neighbour pairs
within its cutoff, plus a temperature logger:

```julia
interactions = [(LennardJones(), CellList3D(info, 4.5, boundary, 100))]
loggers = [TemperatureLogger(100, output = false)]
simulator = VerletProcess(dt = 0.001, thermostat = AndersenThermoStat(1.0, 0.05))

sys = MDSys(
    n_atoms = n_atoms, atoms = atoms, boundary = boundary,
    interactions = interactions, loggers = loggers, simulator = simulator,
)
```

Equilibrate, then sample the rdf by accumulating a pair-distance histogram
every 100 steps:

```julia
simulate!(simulator, sys, info, 1_000_000)   # equilibrate

N, bin_num = 20_000, 100
hist, volume, r, dr = hist_init(N, bin_num, 4.6)

for i in 1:N
    simulate!(simulator, sys, info, 100)
    distance_hist!(hist, sys.interactions[1][2].neighbor_list, dr)
end

rdf = hist ./ (N .* volume)
```

`rdf[k]` is the pair correlation function at radius `r[k]`; a bulk liquid
shows the familiar first-neighbour peak near `r ≈ 2^(1/6) σ` and decays to 1
at large `r`. Plotting it (with `Plots.jl` or any plotting package of your
choice — not a dependency of ExTinyMD itself) is exactly the `plot(r, 2 .*
rdf, ...)` call in `example/LJ_fluid_rdf.jl`.

## Where to go next

- [MD Core](@ref) — the box, atoms, simulators, thermostats, loggers, and
  neighbour finders used above.
- [Interactions](@ref) — the built-in interactions, and how to add your own.
- [Electrostatics](@ref) — Ewald summation and the image-charge method for
  dielectrically confined slabs, usable standalone or inside an `MDSys` like
  the Lennard-Jones interaction above.
