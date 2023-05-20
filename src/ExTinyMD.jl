module ExTinyMD

using LinearAlgebra, CellListMap, Random, Distributions

include("types.jl")
include("system_init.jl")

# this part will be about the MD processes
include("MD_core/Anderson.jl")
include("MD_core/loggers.jl")
include("MD_core/Verlet.jl")
include("MD_core/simulator.jl")

# this part will be about the interactions
include("interactions/lennard_jones.jl")

end
