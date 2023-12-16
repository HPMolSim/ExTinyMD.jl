# for n_atoms in 10.0.^[2.0, 3.0]
@testset "a run for Q2D LJ fluid" begin
    n_atoms = 1000
    n_atoms = Int64(round(n_atoms))
    ρ = 1000 / (100.0)^3
    L = (n_atoms / ρ)^(1/3)
    boundary = Q2dBoundary(L, L, L)
    atoms = create_atoms([(n_atoms ÷ 2, Atom(type = 1, mass = 1.0, charge = 1.0)), (n_atoms ÷ 2, Atom(type = 2, mass = 1.0, charge = - 1.0))])

    info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.5, L - 0.5), boundary; min_r = 2.0, temp = 1.0)

    interactions = [(LennardJones(), CellList3D(info, 4.5, boundary, 100)), (SubLennardJones(0.0, L; cutoff = 1.0, σ = 0.5), SubNeighborFinder(1.5, info, 0.0, L)), (ExternalField(E = Point((0.0, 0.0, 1.0))), NoNeighborFinder())]

    loggers = [TemperatureLogger(100, output = false), TrajectoryLogger(step = 100, output = false)]
    simulator = VerletProcess(dt = 0.001, thermostat = AndersenThermoStat(1.0, 0.05))

    sys = MDSys(
        n_atoms = n_atoms,
        atoms = atoms,
        boundary = boundary,
        interactions = interactions,
        loggers = loggers,
        simulator = simulator
    )
    simulate!(simulator, sys, info, 10000)

    @test info.running_step == 10000
end

@testset "a run for Q2D LJ fluid, BerendsenThermoStat" begin
    n_atoms = 1000
    n_atoms = Int64(round(n_atoms))
    ρ = 1000 / (100.0)^3
    L = (n_atoms / ρ)^(1/3)
    boundary = Q2dBoundary(L, L, L)
    atoms = create_atoms([(n_atoms ÷ 2, Atom(type = 1, mass = 1.0, charge = 1.0)), (n_atoms ÷ 2, Atom(type = 2, mass = 1.0, charge = - 1.0))])

    info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.5, L - 0.5), boundary; min_r = 2.0, temp = 1.0)

    interactions = [(LennardJones(), CellList3D(info, 4.5, boundary, 100)), (SubLennardJones(0.0, L; cutoff = 1.0, σ = 0.5), SubNeighborFinder(1.5, info, 0.0, L)), (ExternalField(E = Point((0.0, 0.0, 1.0))), NoNeighborFinder())]

    loggers = [TemperatureLogger(100, output = false), TrajectoryLogger(step = 100, output = false)]
    simulator = VerletProcess(dt = 0.001, thermostat = BerendsenThermoStat(1.0, 100.0))

    sys = MDSys(
        n_atoms = n_atoms,
        atoms = atoms,
        boundary = boundary,
        interactions = interactions,
        loggers = loggers,
        simulator = simulator
    )
    simulate!(simulator, sys, info, 10000)

    @test info.running_step == 10000
end

@testset "a run for Q2D LJ fluid, NHThermoStat" begin
    n_atoms = 1000
    n_atoms = Int64(round(n_atoms))
    ρ = 1000 / (100.0)^3
    L = (n_atoms / ρ)^(1/3)
    boundary = Q2dBoundary(L, L, L)
    atoms = create_atoms([(n_atoms ÷ 2, Atom(type = 1, mass = 1.0, charge = 1.0)), (n_atoms ÷ 2, Atom(type = 2, mass = 1.0, charge = - 1.0))])

    info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.5, L - 0.5), boundary; min_r = 2.0, temp = 1.0)

    interactions = [(LennardJones(), CellList3D(info, 4.5, boundary, 100)), (SubLennardJones(0.0, L; cutoff = 1.0, σ = 0.5), SubNeighborFinder(1.5, info, 0.0, L)), (ExternalField(E = Point((0.0, 0.0, 1.0))), NoNeighborFinder())]

    simulator = NHVerletProcess(dt = 0.001, Q = 1.0, temperature = 1.0)

    loggers_1 = [TemperatureLogger(100, output = false), TrajectoryLogger(step = 1, output = false)]
    
    sys_1 = MDSys(
        n_atoms = n_atoms,
        atoms = atoms,
        boundary = boundary,
        interactions = interactions,
        loggers = loggers_1,
        simulator = simulator
    )

    simulate!(simulator, sys_1, info, 1000)

    loggers_2 = [TemperatureLogger(100, output = false), TrajectoryLogger(step = 100, output = false)]

    sys_2 = MDSys(
        n_atoms = n_atoms,
        atoms = atoms,
        boundary = boundary,
        interactions = interactions,
        loggers = loggers_2,
        simulator = simulator
    )
    simulate!(simulator, sys_2, info, 10000)

    @test info.running_step == 11000
end