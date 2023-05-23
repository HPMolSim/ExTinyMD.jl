using Random, Distributions, CellListMap, StaticArrays, BenchmarkTools, StatProfilerHTML

include("../src/types.jl")
include("../src/MD_core/system_init.jl")
include("../src/MD_core/Verlet.jl")
include("../src/MD_core/simulator.jl")
include("../src/MD_core/Andersen.jl")
include("../src/MD_core/loggers.jl")
include("../src/MD_core/cell_list.jl")
include("../src/interactions/lennard_jones.jl")


# for n_atoms in 10.0.^[3.0]
n_atoms = 1000
n_atoms = Int64(round(n_atoms))
ρ = 1000 / (100.0)^3
L = (n_atoms / ρ)^(1/3)
boundary = CubicBoundary(L)
atoms = Vector{Atom{Float64}}()

for i in 1:n_atoms
    push!(atoms, Atom(mass = 1.0, charge = 1.0))
end

info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, L), boundary; min_r = 2.0, temp = 1.0)

interactions = [(LennardJones(), CellList3D(info, 4.5, boundary, 100))]
loggers = [TempartureLogger(100)]
simulator = VerletProcess(dt = 0.001, thermostat = AndersenThermoStat(2.5, 0.01))

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

@time info_update_once!(simulator, sys, info, 10000)
# @time (for i in 1:10000 info_update!(simulator, sys, info) end)
# end

# @profilehtml info_update_once!(simulator, sys, info, 100000)
