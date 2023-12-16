struct TrajectoryLogger{TI} <: AbstractLogger 
    step::TI 
    trajectory_file::String
    output::Bool
end

function TrajectoryLogger(;step::TI, trajectory_file::String = "trajectory.txt", output::Bool = true) where{TI<:Integer}
    if output
        # if there is already a file used to store the result, refresh them
        f = open(trajectory_file, "a")
        close(f)
    end
    return TrajectoryLogger{TI}(step, trajectory_file, output)
end

function record!(logger::TrajectoryLogger{TI}, sys::MDSys{T}, info::SimulationInfo{T}) where {T<:Number, TI<:Integer}
    if iszero(info.running_step % logger.step) && logger.output
        IO_trajectory = open(logger.trajectory_file, "a")
        write(IO_trajectory, "step = $(info.running_step) \n")
        write(IO_trajectory, "id, atom_type, x, y, z, vx, vy, vz", "\n")
        for p_info in info.particle_info
            id = p_info.id
            atom_type = sys.atoms[id].type
            x, y, z = p_info.position
            vx, vy, vz = p_info.velocity
            data = "$id, $atom_type, $x, $y, $z, $vx, $vy, $vz"
            write(IO_trajectory, data, "\n")
        end
        write(IO_trajectory, "\n")
        close(IO_trajectory)
    end
    return nothing
end