using ExTinyMD, Plots

begin
    n_atoms = 1000
    n_atoms = Int64(round(n_atoms))
    L = 100.0
    boundary = CubicBoundary(L)
    atoms = Vector{Atom{Float64}}()

    for i in 1:n_atoms
        push!(atoms, Atom(type = 1, mass = 1.0, charge = 0.0))
    end

    info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, L), boundary; min_r = 0.1, temp = 1.0)

    interactions = [(LennardJones(), CellList3D(info, 4.5, boundary, 100))]
    loggers = [TempartureLogger(100, output = false), TrajectionLogger(info, 1000, output = false)]
    simulator = VerletProcess(dt = 0.001, thermostat = AndersenThermoStat(1.0, 0.05))

    sys = MDSys(
        n_atoms = n_atoms,
        atoms = atoms,
        boundary = boundary,
        interactions = interactions,
        loggers = loggers,
        simulator = simulator
    )

    simulate!(simulator, sys, info, 1000000)

    N = 10000
    bin_num = 100

    hist, volume, r, dr = hist_init(N, bin_num, 4.6)

    for i in 1:N
        simulate!(simulator, sys, info, 100)
        distance_hist!(hist, sys.interactions[1][2].neighbor_list, dr)
    end

    rdf = hist ./ (N .* volume)
    plot(r, 2 .* rdf, xlim = (0.0, 4.0), ylim = (0.0, 3.0))
    savefig("rdf_LJ_fluid.png")
end