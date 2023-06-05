# for n_atoms in 10.0.^[2.0, 3.0]
n_atoms = 1000
n_atoms = Int64(round(n_atoms))
ρ = 1000 / (100.0)^3
L = (n_atoms / ρ)^(1/3)
boundary = Q2dBoundary(L, L, L)
atoms = Vector{Atom{Float64}}()

for i in 1:n_atoms
    push!(atoms, Atom(mass = 1.0, charge = 1.0))
end

info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.5, L - 0.5), boundary; min_r = 2.0, temp = 1.0)

interactions = [(LennardJones(), CellListDir3D(info, 4.5, boundary, 100)), (SubLennardJones(0.0, L; cutoff = 1.0, σ = 0.5), SubNeighborFinder(1.5, info.coords, 0.0, L))]

# interactions = [(SubLennardJones(0.0, L; cutoff = 1.0, σ = 0.5), SubNeighborFinder(1.5, info.coords, 0.0, L))]

loggers = [TempartureLogger(100), TrajectionLogger(info, 1000)]
# loggers = [TempartureLogger(100)]
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
@time simulate!(simulator, sys, info, 10000)
# end