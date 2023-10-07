# ExTinyMD

[![Build Status](https://github.com/ArrogantGao/ExTinyMD.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ArrogantGao/ExTinyMD.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Build Status](https://travis-ci.com/ArrogantGao/ExTinyMD.jl.svg?branch=main)](https://travis-ci.com/ArrogantGao/ExTinyMD.jl)
[![Coverage](https://codecov.io/gh/ArrogantGao/ExTinyMD.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ArrogantGao/ExTinyMD.jl)


`ExTinyMD.jl` for Extremely Tiny Molecular Dynamic is a simple package for MD simulations based on `julia`.

Currectly, to use this package, you have to clone this repo onto you machine and manually add it.
```julia
pkg> add ExTinyMD
```

Here is an simple example, which will simulate a 3D LJ fluid and plot it rdf:
```julia
using ExTinyMD, Plots

begin
    # init the particles and the boundary condition
    n_atoms = 1000
    L = 100.0
    boundary = CubicBoundary(L)
    atoms = Vector{Atom{Float64}}()

    for i in 1:n_atoms
        push!(atoms, Atom(mass = 1.0, charge = 0.0))
    end

    info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, L), boundary; min_r = 0.1, temp = 1.0)

    # init the interactions and the loggers, using CellList as neighbor finder
    interactions = [(LennardJones(), CellListDir3D(info, 4.5, boundary, 100))]
    loggers = [TempartureLogger(100, output = false), TrajectionLogger(info, 1000, output = false)]
    simulator = VerletProcess(dt = 0.001, thermostat = AndersenThermoStat(1.0, 0.05))

    sys = MDSys(
        n_atoms = n_atoms,
        atoms = atoms,
        boundary = boundary,
        interactions = interactions,
        loggers = loggers,
        simulator = simulator
    )

    # run for 1e6 steps
    simulate!(simulator, sys, info, 1000000)

    # sample for another 1e6 steps
    N = 10000
    bin_num = 100

    hist, volume, r, dr = hist_init(N, bin_num, 4.6)

    for i in 1:N
        simulate!(simulator, sys, info, 100)
        distance_hist!(hist, sys.interactions[1][2].neighbor_list, dr)
    end

    # plot the rdf
    rdf = hist ./ (N .* volume)
    plot(r, 2 .* rdf, xlim = (0.0, 4.0), ylim = (0.0, 3.0))
    savefig("rdf_LJ_fluid.png")
end
```
