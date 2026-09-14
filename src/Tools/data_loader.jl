"""
    load_trajectory(filename::String) -> (trajectory_list, step_list)

Load a trajectory written by [`TrajectoryLogger`](@ref) from `filename`.
`trajectory_list[k]` is a vector of `[id, type, x, y, z, vx, vy, vz]` rows for
recorded step `step_list[k]`.
"""
function load_trajectory(filename::String)
    f = open(filename)

    trajectory_list = []
    step_list = []

    while !eof(f)
        step_string = split(readline(f), "=")
        step = parse(Int64, step_string[2])
        push!(step_list, step)

        readline(f)

        trajectory = []
        while (line = readline(f)) != ""
            data_string = split(line, ",")
            data_num = parse.(Float64, data_string)
            push!(trajectory, data_num)
        end

        push!(trajectory_list, trajectory)
    end

    close(f)

    return trajectory_list, step_list
end

"""
    data2info(data::Vector, boundary, masses, charges) -> (atoms, info)

Build `atoms` and a [`SimulationInfo`](@ref) from one frame of `data` (as
returned by [`load_trajectory`](@ref) or [`load_lammpstrj`](@ref)), looking up
each particle's mass and charge in `masses`/`charges` by its type.
"""
function data2info(data::Vector, boundary::Boundary{T}, masses::Vector{T}, charges::Vector{T}) where{T}
    L = boundary.length
	n_atoms = length(data)
	atoms = Vector{Atom{T}}()
	for i in 1:n_atoms
		type = data[i][2]
		q = charges[Int(type)]
        m = masses[Int(type)]
		push!(atoms, Atom(type = Int(type), mass = m, charge = q))
	end
	info = SimulationInfo(n_atoms, atoms, (0.0, L[1], 0.0, L[2], 0.0, L[3]), boundary; min_r = 0.0, temp = 1.0)
	for i in 1:n_atoms
		x, y, z = data[i][3], data[i][4], data[i][5]
		info.particle_info[i].position = Point(x, y, z)
	end
	return atoms, info
end

"""
    load_lammpstrj(filename::String) -> (trajectory, box, timestep, num_atoms)

Parse a LAMMPS dump file (`ITEM: TIMESTEP` / `ITEM: NUMBER OF ATOMS` /
`ITEM: BOX BOUNDS pp pp ff` / `ITEM: ATOMS id type x y z` blocks). Returns each
frame's `(id, type, x, y, z)` tuples in `trajectory`, the box bounds `box =
(xlo, xhi, ylo, yhi, zlo, zhi)`, the recorded `timestep`s, and `num_atoms`.
"""
function load_lammpstrj(filename::String)
    file = open(filename, "r")

    timestep = Vector{Int64}()
    num_atoms = 0
    box_bounds = [[0.0, 0.0] for i in 1:3]
    trajectory = Vector{Vector{Tuple{Int, Int, Float64, Float64, Float64}}}()

    for line in eachline(file)
        if line == "ITEM: TIMESTEP"
            push!(timestep, parse(Int, readline(file)))
        elseif line == "ITEM: NUMBER OF ATOMS"
            num_atoms = parse(Int, readline(file))
        elseif line == "ITEM: BOX BOUNDS pp pp ff"
            for i in 1:3
                bounds = split(readline(file))
                box_bounds[i] = [parse(Float64, bounds[1]), parse(Float64, bounds[2])]
            end
        elseif line == "ITEM: ATOMS id type x y z "
            data_step = Vector{Tuple{Int, Int, Float64, Float64, Float64}}()
            for _ in 1:num_atoms
                atom_line = readline(file)
                d = split(atom_line)
                atom_data = (parse(Int, d[1]), parse(Int, d[2]), parse(Float64, d[3]), parse(Float64, d[4]), parse(Float64, d[5]))
                push!(data_step, atom_data)
            end
            push!(trajectory, data_step)
        end
    end

    close(file)

    box = (box_bounds[1][1], box_bounds[1][2], box_bounds[2][1], box_bounds[2][2], box_bounds[3][1], box_bounds[3][2])

    return trajectory, box, timestep, num_atoms
end