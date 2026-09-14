using ExTinyMD
using Test
using Random
using SpecialFunctions

include("electrostatics/reference.jl")

@testset "ExTinyMD.jl" begin
    include("simulation.jl")
    include("cell_list.jl")
    include("regression_finder.jl")
end

@testset "electrostatics" begin
    include("electrostatics/test_reference.jl")
    include("electrostatics/test_common.jl")
    include("electrostatics/test_short.jl")
    include("electrostatics/test_long_ewald3d.jl")
    include("electrostatics/test_ewald.jl")
end
