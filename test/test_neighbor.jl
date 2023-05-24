using Random, Distributions, CellListMap, StaticArrays, BenchmarkTools, StatProfilerHTML, Plots

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
ρ = 1000 / (100.0)^3
L = (n_atoms / ρ)^(1/3)
boundary = CubicBoundary(L)
atoms = Vector{Atom{Float64}}()

for i in 1:n_atoms
    push!(atoms, Atom(mass = 1.0, charge = 1.0))
end

info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, L), boundary; min_r = 2.0, temp = 1.0)

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

function distance_hist!(hist::Vector{T}, n_list::Vector{Tuple{Int64, Int64, T}}, d::T) where T
    for (i, j, r) in n_list
        bin = Int64(round(r / d)) + 1
        hist[bin] += 1
    end    
    return nothing
end

N = 100000
bin_num = 200
bin_dr = 4.5 / (bin_num - 1)
hist = zeros(Float64, bin_num);

r = [bin_dr * i for i in 1:bin_num]
volume = 4π .* r.^2 * bin_dr

for i in 1:N
    simulate!(simulator, sys, info, 100)
    distance_hist!(hist, sys.interactions[1][2].neighbor_list, bin_dr)
end

rdf = hist ./ (N .* volume)
plot(r, 2 .* rdf, xlim = (0.0, 4.0), ylim = (0.0, 2.5))