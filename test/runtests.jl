using ExTinyMD
using Test

@testset "ExTinyMD.jl" begin
    include("simulation.jl")
    include("cell_list.jl")
    include("regression_finder.jl")
end
