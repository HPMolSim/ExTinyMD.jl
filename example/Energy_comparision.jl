using ExTinyMD
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

function compute_rbe_energy!(rbe_interaction, neighborfinder, sys, info)
    update_acceleration!(rbe_interaction, neighborfinder, sys, info)
    return rbe_interaction.total_energy[end]
end

function compute_ewald_energy(ewald_interaction, neighborfinder, info, atoms)
    return ExTinyMD.energy(ewald_interaction, neighborfinder, info, atoms)
end

function compare_and_save_energy!(sys, info, num_steps, filename)
    open(filename, "w") do io
        for step in 1:num_steps
            rbe_energy = compute_rbe_energy!(rbe_interaction, neighborfinder_rbe, sys, info)
            ewald_energy = compute_ewald_energy(ewald_interaction, neighborfinder_ewald, info, atoms)
            energy_diff = rbe_energy - ewald_energy
            
            println(io, "Step $step: Difference = $energy_diff")
            
            simulate!(sys.simulator, sys, info, 1)
        end
    end
end

loggers = [TrajectoryLogger(;step = 10, trajectory_file="trajectory.txt", output = true), TemperatureLogger(100)]

sys = MDSys(
    n_atoms = n_atoms,
    atoms = atoms,
    boundary = boundary,
    interactions = [(rbe_interaction, neighborfinder_rbe)],
    loggers = loggers,
    simulator = VerletProcess(0.001, AndersenThermoStat(1.0, 0.05))
)

simulator = VerletProcess(dt = 0.001, thermostat = AndersenThermoStat(1.0, 0.05))

num_steps = 1000
compare_and_save_energy!(sys, info, num_steps, "energy_comparison.txt")
