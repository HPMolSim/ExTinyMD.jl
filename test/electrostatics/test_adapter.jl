function _charged_system(n, L; charge = 1.0)
    boundary = Boundary((L, L, L), (1, 1, 1))
    atoms = Vector{Atom{Float64}}()
    for i in 1:(n ÷ 2)
        push!(atoms, Atom(type = 1, mass = 1.0, charge = charge))
    end
    for i in (n ÷ 2 + 1):n
        push!(atoms, Atom(type = 2, mass = 1.0, charge = -charge))
    end
    info = SimulationInfo(n, atoms, (0.0, L, 0.0, L, 0.0, L), boundary;
                          min_r = 1.0, temp = 1.0)
    info.running_step = 1
    return boundary, atoms, info
end

@testset "adapter energy matches the core API" begin
    Random.seed!(20260930)
    n, L = 20, 10.0
    boundary, atoms, info = _charged_system(n, L)
    inter = Ewald3D(n, (L, L, L); α = 1.0, s = 4.0)   # r_c = 4.0 < L/2 = 5
    finder = CellList3D(info, inter.short.r_c, boundary, 1)

    sys = MDSys(n_atoms = n, atoms = atoms, boundary = boundary,
                interactions = [(inter, finder)],
                loggers = [TemperatureLogger(100; output = false)],
                simulator = VerletProcess(dt = 0.001))

    poses = [SVector(p.position[1], p.position[2], p.position[3])
             for p in info.particle_info]
    charges = [atoms[p.id].charge for p in info.particle_info]

    @test isapprox(energy(inter, finder, sys, info),
                   coulomb_energy(inter, poses, charges), rtol = 1e-10)
end

@testset "adapter writes acceleration = force/mass" begin
    Random.seed!(20260931)
    n, L = 16, 10.0
    boundary, atoms, info = _charged_system(n, L)
    # give the two species different masses so a missing division shows up
    atoms = [Atom(type = a.type, mass = a.type == 1 ? 1.0 : 4.0, charge = a.charge)
             for a in atoms]
    inter = Ewald3D(n, (L, L, L); α = 1.0, s = 4.0)   # r_c = 4.0 < L/2 = 5
    finder = CellList3D(info, inter.short.r_c, boundary, 1)
    sys = MDSys(n_atoms = n, atoms = atoms, boundary = boundary,
                interactions = [(inter, finder)],
                loggers = [TemperatureLogger(100; output = false)],
                simulator = VerletProcess(dt = 0.001))

    poses = [SVector(p.position[1], p.position[2], p.position[3])
             for p in info.particle_info]
    charges = [atoms[p.id].charge for p in info.particle_info]
    F = coulomb_force(inter, poses, charges)

    for p in info.particle_info
        p.acceleration = Point(0.0, 0.0, 0.0)
    end
    ExTinyMD.update_acceleration!(inter, finder, sys, info)

    for (slot, p) in enumerate(info.particle_info)
        m = atoms[p.id].mass
        for d in 1:3
            @test isapprox(p.acceleration[d], F[slot][d] / m, rtol = 1e-10)
        end
    end
end

@testset "adapter is correct when slot order differs from id order" begin
    # In stock ExTinyMD `particle_info[i].id == i`, so a gather that confuses slot
    # with id looks correct forever. Permute the mapping so the two differ and the
    # confusion becomes observable: charges must follow ids, positions must follow
    # slots. `substrate_lennard_jones.jl` already relies on this indirection via
    # `info.id_dict`, so it is a real invariant, not a hypothetical one.
    Random.seed!(20260934)
    n, L = 12, 10.0
    boundary, atoms, info = _charged_system(n, L)

    # give every id a distinct charge so a slot/id mix-up cannot cancel out
    atoms = [Atom(type = a.type, mass = 1.0 + 0.1 * i, charge = (isodd(i) ? 1.0 : -1.0) * (1 + 0.01 * i))
             for (i, a) in enumerate(atoms)]

    # reverse the slot order, keeping ids attached to their particles
    reverse!(info.particle_info)
    for i in eachindex(info.particle_info)
        info.id_dict[info.particle_info[i].id] = i
    end
    @test info.particle_info[1].id != 1        # the mapping really is permuted

    inter = Ewald3D(n, (L, L, L); α = 1.0, s = 4.0)
    finder = CellList3D(info, inter.short.r_c, boundary, 1)
    sys = MDSys(n_atoms = n, atoms = atoms, boundary = boundary,
                interactions = [(inter, finder)],
                loggers = [TemperatureLogger(100; output = false)],
                simulator = VerletProcess(dt = 0.001))

    poses   = [SVector(p.position[1], p.position[2], p.position[3])
               for p in info.particle_info]
    charges = [atoms[p.id].charge for p in info.particle_info]
    @test isapprox(energy(inter, finder, sys, info),
                   coulomb_energy(inter, poses, charges), rtol = 1e-10)

    # and the acceleration must land on the right particle
    F = coulomb_force(inter, poses, charges)
    for p in info.particle_info
        p.acceleration = Point(0.0, 0.0, 0.0)
    end
    ExTinyMD.update_acceleration!(inter, finder, sys, info)
    for (slot, p) in enumerate(info.particle_info)
        m = atoms[p.id].mass
        for d in 1:3
            @test isapprox(p.acceleration[d], F[slot][d] / m, rtol = 1e-10)
        end
    end
end

@testset "adapter runs inside simulate! with bounded energy drift" begin
    Random.seed!(20260932)
    n, L = 30, 12.0
    boundary, atoms, info = _charged_system(n, L)
    inter = Ewald3D(n, (L, L, L); α = 0.8, s = 3.5)   # r_c = 4.375 < L/2 = 6
    lj = LennardJones(ϵ = 1.0, σ = 1.0, cutoff = 3.0)
    finder = CellList3D(info, max(inter.short.r_c, 3.0), boundary, 1)

    sys = MDSys(n_atoms = n, atoms = atoms, boundary = boundary,
                interactions = [(lj, finder), (inter, finder)],
                loggers = [TemperatureLogger(1000; output = false)],
                simulator = VerletProcess(dt = 1e-4))

    E0 = energy(inter, finder, sys, info)
    simulate!(sys.simulator, sys, info, 200)
    E1 = energy(inter, finder, sys, info)

    @test isfinite(E1)
    # microcanonical Verlet at this dt should not let the electrostatic energy run
    # away; a sign error or a mass bug shows up as an unbounded value
    @test abs(E1 - E0) < 0.5 * max(abs(E0), 1.0)
end

@testset "adapter works for ICM with NoNeighborFinder" begin
    Random.seed!(20260933)
    n, L = 12, 8.0
    boundary, atoms, info = _charged_system(n, L)
    inter = ICMEwald2D(n, (L, L, L); α = 1.0, s = 3.0, γ = (0.3, 0.3), N_image = 3)   # r_c = 3.0 < 4
    finder = NoNeighborFinder(Float64)   # ICM keeps its own cell list
    sys = MDSys(n_atoms = n, atoms = atoms, boundary = boundary,
                interactions = [(inter, finder)],
                loggers = [TemperatureLogger(100; output = false)],
                simulator = VerletProcess(dt = 0.001))

    poses = [SVector(p.position[1], p.position[2], p.position[3])
             for p in info.particle_info]
    charges = [atoms[p.id].charge for p in info.particle_info]
    @test isapprox(energy(inter, finder, sys, info),
                   coulomb_energy(inter, poses, charges), rtol = 1e-10)
end
