using ExTinyMD
using Test

@testset "ExTinyMD.jl" begin
    include("simulation.jl")
    include("cell_list.jl")
    include("rbe.jl")
end
