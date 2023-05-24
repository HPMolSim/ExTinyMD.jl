module ExTinyMD

using LinearAlgebra, Random, Distributions, CellListMap, StaticArrays, BenchmarkTools, StatProfilerHTML, Plots


include("types.jl")


# this part will be about the MD processes
include("MD_core/system_init.jl")
include("MD_core/Andersen.jl")
include("MD_core/loggers.jl")
include("MD_core/Verlet.jl")
include("MD_core/simulator.jl")
include("MD_core/cell_list.jl")

# this part will be about the interactions
include("interactions/lennard_jones.jl")

end
