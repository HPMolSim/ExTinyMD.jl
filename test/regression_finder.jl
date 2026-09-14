@testset "neighbor finders expose neighbor_list" begin
    # Regression: AllNeighborFinder/NoNeighborFinder named the field `neighborlist`
    # while every interaction reads `.neighbor_list`, so LJ + AllNeighborFinder threw.
    @test hasfield(ExTinyMD.AllNeighborFinder{Float64}, :neighbor_list)
    @test hasfield(ExTinyMD.NoNeighborFinder{Float64}, :neighbor_list)

    n_atoms = 20
    L = 20.0
    boundary = CubicBoundary(L)
    atoms = create_atoms([(n_atoms, Atom(type = 1, mass = 1.0))])
    info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, L), boundary;
                          min_r = 1.0, temp = 1.0)
    info.running_step = 1

    interaction = LennardJones(ϵ = 1.0, σ = 1.0, cutoff = 4.0)
    all_finder = AllNeighborFinder(n_atoms, Float64)

    # must not throw, and must produce a finite energy
    E = energy(interaction, all_finder, MDSys(
        n_atoms = n_atoms, atoms = atoms, boundary = boundary,
        interactions = [(interaction, all_finder)],
        loggers = [TemperatureLogger(100; output = false)],
        simulator = VerletProcess(dt = 0.001),
    ), info)
    @test isfinite(E)

    # NoNeighborFinder yields exactly zero pair energy for a cutoff-respecting interaction
    no_finder = NoNeighborFinder(Float64)
    @test energy(interaction, no_finder, MDSys(
        n_atoms = n_atoms, atoms = atoms, boundary = boundary,
        interactions = [(interaction, no_finder)],
        loggers = [TemperatureLogger(100; output = false)],
        simulator = VerletProcess(dt = 0.001),
    ), info) == 0.0
end
