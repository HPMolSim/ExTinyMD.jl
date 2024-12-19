@testset "test the rbe calculation" begin
    n_atoms = 300
    L_x = 100.0
    L_y = 100.0
    L_z = 100.0
    L = (L_x, L_y, L_z)
    s = 0.3
    p = 10
    α = 0.1
    boundary = Boundary((L_x, L_y, L_z), (1, 1, 1))
    atoms = create_atoms([(n_atoms ÷ 2, Atom(type = 1, mass = 1.0, charge = 1.0)), (n_atoms ÷ 2, Atom(type = 2, mass = 1.0, charge = - 1.0))])
    info = SimulationInfo(n_atoms, atoms, (0.0, L_x, 0.0, L_y, 0.0, L_z), boundary; min_r = 0.1, temp = 1.0)
    interactions = [
        (LennardJones(), CellList3D(info, 4.5, boundary, 100)),
        (RBEInteractions(α, L_x, L_y, L_z, p, s), CellList3D(info, 4.5, boundary, 100)),
    ]
    loggers = [TrajectoryLogger(;step = 10, trajectory_file="trajectory.txt",output = true), TemperatureLogger(100)]
    simulator = VerletProcess(dt = 0.001, thermostat = AndersenThermoStat(1.0, 0.05))
    sys = MDSys(
        n_atoms = n_atoms,
        atoms = atoms,
        boundary = boundary,
        interactions = interactions,
        loggers = loggers,
        simulator = simulator
    )
    
    simulate!(simulator, sys, info, 1000)
    @test info.running_step == 1000
end