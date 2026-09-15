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

@testset "SubLennardJones looks up mass by id, not slot" begin
    # Regression: update_acceleration! converted the neighbour id to a storage slot
    # via id_dict, then indexed the id-keyed sys.atoms array with that slot. In stock
    # ExTinyMD slot == id, so the bug was invisible — and test/simulation.jl exercises
    # SubLennardJones only with uniform masses, which masks it further. It becomes
    # observable with distinct per-id masses and a permuted slot order.
    n = 4
    L = 10.0
    boundary = Q2dBoundary(L, L, L)
    atoms = [Atom(type = 1, mass = Float64(i), charge = 0.0) for i in 1:n]
    info = SimulationInfo(n, atoms, (0.0, L, 0.0, L, 0.4, 0.6), boundary; temp = 1.0)

    # place every particle at the same height, inside the substrate cutoff, so the
    # force is identical and only the mass division distinguishes the particles
    for p in info.particle_info
        p.position = Point(p.position[1], p.position[2], 0.5)
    end

    interaction = SubLennardJones(0.0, L; cutoff = 1.0, σ = 0.5)
    finder = SubNeighborFinder(1.5, info, 0.0, L)
    sys = MDSys(n_atoms = n, atoms = atoms, boundary = boundary,
                interactions = [(interaction, finder)],
                loggers = [TemperatureLogger(100; output = false)],
                simulator = VerletProcess(dt = 0.001))

    # baseline, with slot == id
    for p in info.particle_info
        p.acceleration = Point(0.0, 0.0, 0.0)
    end
    ExTinyMD.update_acceleration!(interaction, finder, sys, info)
    a_by_id = Dict(p.id => p.acceleration[3] for p in info.particle_info)

    # permute the slot order; the acceleration each *id* receives must be unchanged
    reverse!(info.particle_info)
    for i in eachindex(info.particle_info)
        info.id_dict[info.particle_info[i].id] = i
    end
    @test info.particle_info[1].id != 1
    for p in info.particle_info
        p.acceleration = Point(0.0, 0.0, 0.0)
    end
    ExTinyMD.update_acceleration!(interaction, finder, sys, info)

    for p in info.particle_info
        @test isapprox(p.acceleration[3], a_by_id[p.id]; rtol = 1e-12)
    end
end
