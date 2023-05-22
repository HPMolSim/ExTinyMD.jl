include("../src/types.jl")
include("../src/MD_core/system_init.jl")
include("../src/MD_core/simulator.jl")
include("../src/MD_core/Verlet.jl")
include("../src/MD_core/Andersen.jl")
include("../src/MD_core/loggers.jl")

using Random, Distributions, CellListMap, StaticArrays, BenchmarkTools

n_atoms = 1000
boundary = CubicBoundary(100.0)
atoms = Vector{Atom{Float64}}()

for i in 1:n_atoms/2
    push!(atoms, Atom(mass = 1.0, charge = 1.0))
    push!(atoms, Atom(mass = 1.0, charge = - 1.0))
end

interactions = [(NoInteraction(), NoNeighborFinder(n_atoms))]
logger = [TempartureLogger(100)]
simulator = VerletProcess(dt = 0.01, thermostat = AndersenThermoStat(1.0, 0.01))

sys = MDSys(
    n_atoms = n_atoms,
    atoms = atoms,
    boundary = boundary,
    interactions = interactions,
    logger = logger,
    simulator = simulator
)

info = SimulationInfo(n_atoms, atoms, (0.0, 100.0, 0.0, 100.0, 0.0, 100.0), boundary; min_r = 0.1)

cell3d = CellList3D(info, 1.0, boundary, 10)

info.running_step += 7

@btime update_finder!(cell3d, info)

boundary_q2d = Q2dBoudary(10.0, 10.0, 10.0)
cellq2d = CellListQ2D(info, 1.0, boundary_q2d, 10)
@btime update_finder!(cellq2d, info)