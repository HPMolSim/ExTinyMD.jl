# ExTinyMD

[![Build Status](https://github.com/ArrogantGao/ExTinyMD.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ArrogantGao/ExTinyMD.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/ArrogantGao/ExTinyMD.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ArrogantGao/ExTinyMD.jl)


`ExTinyMD.jl` for Extremely Tiny Molecular Dynamic is a simple package written in `Julia`, which provide a simple software for MD simulations.

## Getting Started

You can simple type `]` to enter the package in Julia REPL and type
```julia
pkg> add ExTinyMD
```
to install the package.

Here is an example, which simulate a 3D LJ fluid and plot its rdf:
```julia
using ExTinyMD, Plots

begin
    # create the atoms and the box
    n_atoms = 1000
    n_atoms = Int64(round(n_atoms))
    L = 100.0
    boundary = CubicBoundary(L)
    atoms = Vector{Atom{Float64}}()

    for i in 1:n_atoms
        push!(atoms, Atom(type = 1, mass = 1.0, charge = 0.0))
    end

    # random init, position and velocity
    info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, L), boundary; min_r = 0.1, temp = 1.0)

    # set up the interaction needed, here only LJ, neighbors by cell_list
    interactions = [(LennardJones(), CellList3D(info, 4.5, boundary, 100))]

    # loggers, will store the data during simulations
    loggers = [TemperatureLogger(100, output = false), TrajectionLogger(step = 1000, output = false)]

    # simulator and thermostat
    simulator = VerletProcess(dt = 0.001, thermostat = AndersenThermoStat(1.0, 0.05))

    # create MDSys
    sys = MDSys(
        n_atoms = n_atoms,
        atoms = atoms,
        boundary = boundary,
        interactions = interactions,
        loggers = loggers,
        simulator = simulator
    )

    # run 1e6 steps
    simulate!(simulator, sys, info, 1000000)

    # sample 2e6 steps to get rdf
    N = 20000
    bin_num = 100

    hist, volume, r, dr = hist_init(N, bin_num, 4.6)

    for i in 1:N
        simulate!(simulator, sys, info, 100)
        distance_hist!(hist, sys.interactions[1][2].neighbor_list, dr)
    end

    rdf = hist ./ (N .* volume)
    plot(r, 2 .* rdf, xlim = (0.0, 4.0), ylim = (0.0, 3.0))
    savefig("rdf_LJ_fluid.png")
end
```

## How to Contribute

If you find any bug or have any suggestion, please open an [issue](https://github.com/ArrogantGao/ExTinyMD.jl/issues).

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
and run them together. The package [QuasiEwald.jl](https://github.com/ArrogantGao/QuasiEwald.jl) can be an example.