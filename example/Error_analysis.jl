using ExTinyMD
using BenchmarkTools
using EwaldSummations
using LinearAlgebra

n_atoms = 300
L_x, L_y, L_z = 100.0, 100.0, 100.0
L = (L_x, L_y, L_z)
s = 0.3
p = 100
α = 0.1
boundary = Boundary((L_x, L_y, L_z), (1, 1, 1))

atoms = create_atoms([(n_atoms ÷ 2, Atom(type = 1, mass = 1.0, charge = 1.0)), (n_atoms ÷ 2, Atom(type = 2, mass = 1.0, charge = -1.0))])
info = SimulationInfo(n_atoms, atoms, (0.0, L_x, 0.0, L_y, 0.0, L_z), boundary; min_r = 0.1, temp = 1.0)

rbe_interaction = RBEInteractions(α, L_x, L_y, L_z, p, s)
neighborfinder_rbe = CellList3D(info, 4.5, boundary, 100)

ewald_interaction = Ewald3DInteraction(n_atoms, s, α, L)
neighborfinder_ewald = CellList3D(info, 4.5, boundary, 100)

function compute_rbe_forces!(rbe_interaction, neighborfinder, sys, info)
    update_acceleration!(rbe_interaction, neighborfinder, sys, info)
    return [info.particle_info[i].acceleration for i in 1:sys.n_atoms]
end

function compute_ewald_forces(ewald_interaction, neighborfinder, info, atoms)
    return force(ewald_interaction, neighborfinder, info, atoms)
end

function calculate_errors(rbe_forces, ewald_forces)
    n = length(rbe_forces)
    errors = [norm(rbe_forces[i] - ewald_forces[i]) for i in 1:n]
    mae = mean(errors)
    rmse = sqrt(mean(errors .^ 2))
    return mae, rmse
end

function save_errors(filename, step, error)
    open(filename, "a") do io
        println(io, "Step $step: $error")
    end
end

function simulate_with_comparison!(sys, info, num_steps)
    mae_list = []
    rmse_list = []

    for step in 1:num_steps
        rbe_forces = compute_rbe_forces!(rbe_interaction, neighborfinder_rbe, sys, info)
        ewald_forces = compute_ewald_forces(ewald_interaction, neighborfinder_ewald, info, atoms)

        mae, rmse = calculate_errors(rbe_forces, ewald_forces)
        push!(mae_list, mae)
        push!(rmse_list, rmse)

        save_errors("mae_errors.txt", step, mae)
        save_errors("rmse_errors.txt", step, rmse)

        simulate!(sys.simulator, sys, info, 1)
    end

    println("MAE over $num_steps steps: ", mean(mae_list))
    println("RMSE over $num_steps steps: ", mean(rmse_list))
end

loggers = [TrajectoryLogger(;step = 10, trajectory_file="trajectory.txt",output = true), TemperatureLogger(100)]

sys = MDSys(
    n_atoms = n_atoms,
    atoms = atoms,
    boundary = boundary,
    interactions = [(rbe_interaction, neighborfinder_rbe)],
    loggers = loggers,
    simulator = VerletProcess(0.001, AndersenThermoStat(1.0, 0.05))
)

simulator = VerletProcess(dt = 0.001, thermostat = AndersenThermoStat(1.0, 0.05))

num_steps = 10000

simulate_with_comparison!(sys, info, num_steps)
