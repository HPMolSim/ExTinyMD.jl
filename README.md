# ExTinyMD

[![Build Status](https://github.com/HPMolSim/ExTinyMD.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/HPMolSim/ExTinyMD.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/HPMolSim/ExTinyMD.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/HPMolSim/ExTinyMD.jl)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://HPMolSim.github.io/ExTinyMD.jl/dev)


`ExTinyMD.jl` for Extremely Tiny Molecular Dynamic is a simple package written in `Julia`, which provide a simple software for MD simulations.

## Getting Started

You can simple type `]` to enter the package in Julia REPL and type
```julia
pkg> add ExTinyMD
```
to install the package.

Here is a minimal simulation — a Lennard-Jones fluid, thermostatted, run for a
few thousand steps:

```julia
using ExTinyMD

n_atoms  = 1000
L        = 100.0
boundary = CubicBoundary(L)
atoms    = create_atoms([(n_atoms, Atom(type = 1, mass = 1.0, charge = 0.0))])

info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, L), boundary;
                      min_r = 0.1, temp = 1.0)

sys = MDSys(
    n_atoms      = n_atoms,
    atoms        = atoms,
    boundary     = boundary,
    interactions = [(LennardJones(), CellList3D(info, 4.5, boundary, 100))],
    loggers      = [TemperatureLogger(100)],
    simulator    = VerletProcess(dt = 0.001,
                                 thermostat = AndersenThermoStat(1.0, 0.05)),
)

simulate!(sys.simulator, sys, info, 10_000)
```

The [documentation](https://HPMolSim.github.io/ExTinyMD.jl/dev) carries the full
version of this walkthrough, including how to accumulate a radial distribution
function from the run, along with reference pages for the simulators,
thermostats, loggers and neighbour finders.

## Electrostatics

ExTinyMD also includes a small electrostatics standard library: Ewald
summation and the image-charge method for dielectrically confined slabs,
usable either standalone (a plan object plus plain arrays) or as an ordinary
`MDSys` interaction. Six methods are available:

- `Ewald3D` — triply periodic Ewald summation.
- `PME3D` — particle-mesh Ewald for the same triply periodic system, computing
  the identical sum as `Ewald3D` in `O(N log N)` instead of `O(N·K)`. Requires
  `using FINUFFT` (a weak dependency, loaded via a package extension).
- `Ewald2D` — exact Ewald summation for a slab periodic in x, y and free in z.
- `ICMEwald2D` — `Ewald2D` plus the image-charge method, for a slab confined
  between two dielectric walls.
- `ICMEwald3D` — the image-charge method combined with `Ewald3D` and an
  electrostatic layer correction, for the same confined slab at lower cost.
- `ICMPME3D` — the image-charge method combined with `PME3D` and the same
  layer correction, the `O(N log N)` counterpart of `ICMEwald3D`. Also
  requires `using FINUFFT`.

See the [Electrostatics](https://HPMolSim.github.io/ExTinyMD.jl/dev/electrostatics/)
page of the documentation for a full guide, including how to choose `α`/`s`
and the image-charge parameters.

## How to Contribute

If you find any bug or have any suggestion, please open an [issue](https://github.com/HPMolSim/ExTinyMD.jl/issues).

If you want to add some new features, such as force or loggers, you can simply define something in your own package as
```julia
function ExTinyMD.update_acceleration!(
    interaction::YourForce, 
    neighborfinder::YourFinder, 
    sys::MDSys{T}, 
    info::SimulationInfo{T}) where {T<:Number}

    update_finder!(neighborfinder, info)
    YourForce!(interaction, neighborfinder, sys.atoms, sys.boundary, info)

    return nothing
end
```
and run them together — ExTinyMD will call it like any built-in interaction.
The electrostatics library under `src/interactions/electrostatics/` is a worked
example of the same contract, and the
[Interactions](https://HPMolSim.github.io/ExTinyMD.jl/dev/interactions/) page
documents it.

(The companion packages [QuasiEwald.jl](https://github.com/HPMolSim/QuasiEwald.jl),
[SoEwald2D.jl](https://github.com/HPMolSim/SoEwald2D.jl) and
[FastSpecSoG.jl](https://github.com/HPMolSim/FastSpecSoG.jl) implement this
interface too, but currently pin `ExTinyMD = "0.2"` and do not yet resolve
against the current release.)
