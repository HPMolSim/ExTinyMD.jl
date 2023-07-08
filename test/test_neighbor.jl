using Random, Distributions, CellListMap, StaticArrays, BenchmarkTools, StatProfilerHTML, Plots, DelimitedFiles

include("../src/types.jl")
include("../src/MD_core/system_init.jl")
include("../src/MD_core/Verlet.jl")
include("../src/MD_core/simulator.jl")
include("../src/MD_core/Andersen.jl")
include("../src/MD_core/loggers.jl")
include("../src/MD_core/cell_list.jl")
include("../src/interactions/lennard_jones.jl")


# for n_atoms in 10.0.^[2.0, 3.0]
n_atoms = 1000
n_atoms = Int64(round(n_atoms))
# ρ = 1000 / (100.0)^3
L = 100.0
boundary = CubicBoundary(L)
atoms = Vector{Atom{Float64}}()

for i in 1:n_atoms
    push!(atoms, Atom(mass = 1.0, charge = 0.0))
end

info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, L), boundary; min_r = 0.1, temp = 1.0)

interactions = [(LennardJones(), CellListDir3D(info, 4.5, boundary, 100))]
loggers = [TempartureLogger(100)]
simulator = VerletProcess(dt = 0.001, thermostat = AndersenThermoStat(1.0, 0.05))

sys = MDSys(
    n_atoms = n_atoms,
    atoms = atoms,
    boundary = boundary,
    interactions = interactions,
    loggers = loggers,
    simulator = simulator
)
@show (n_atoms, L)
simulate!(simulator, sys, info, 1)
@time simulate!(simulator, sys, info, 1000000)
# end

@test 1 == 1



N = 100
bin_num = 200

hist, volume, r, dr = hist_init(N, bin_num, 4.6)

for i in 1:N
    simulate!(simulator, sys, info, 100)
    distance_hist!(hist, sys.interactions[1][2].neighbor_list, dr)
end

rdf = hist ./ (N .* volume)
plot(r, 2 .* rdf, xlim = (0.0, 4.0), ylim = (0.0, 2.5))