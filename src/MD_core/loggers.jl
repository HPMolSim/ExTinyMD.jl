mutable struct TemperatureLogger{T, TI} <: AbstractLogger
    step::TI
    data::Vector{T}
    output::Bool
end

mutable struct PressureLogger{T, TI} <: AbstractLogger
    step::TI
    data::Vector{T}
    output::Bool
end

struct TrajectionLogger{TI} <: AbstractLogger 
    step::TI 
    trajection_file::String
    output::Bool
end


# TemperatureLogger(step::TI) where{TI<:Integer} = TemperatureLogger{Float64, TI}(step, Vector{Float64}())

function TemperatureLogger(step::TI; output::Bool = true) where{TI<:Integer}
    if output
        f = open("temperature.txt", "w")
        close(f)
    end
    return TemperatureLogger{Float64, TI}(step, Vector{Float64}(), output)
end
TemperatureLogger{T}(step::TI; output::Bool = true) where{TI<:Integer, T<:Number} = TemperatureLogger{T, TI}(step, Vector{T}(), output)

PressureLogger(step::TI; output::Bool = true) where{TI<:Integer} = PressureLogger{Float64}(step, Vector{Float64}(), output)
PressureLogger{T}(step::TI; output::Bool = true) where{TI<:Integer, T<:Number} = PressureLogger{T}(step, Vector{T}(), output)

function TrajectionLogger(;step::TI=100, trajection_file::String = "trajection.txt", output::Bool = true) where{TI<:Integer}
    if output
        # if there is already a file used to store the result, refresh them
        f = open(trajection_file, "w")
        close(f)
    end
    return TrajectionLogger{TI}(step, trajection_file, output)
end

@inbounds function temperature(sys::MDSys{T}, info::SimulationInfo{T}) where {T<:Number}
    Ek = zero(T)
    for p_info in info.particle_info
        id = p_info.id
        Ek += 0.5 * sys.atoms[id].mass * sum(abs2, p_info.velocity)
    end
    temperature = T(2.0 * Ek / (3.0 * sys.n_atoms))
    return temperature
end

# this function is used to record Temperature
function record!(logger::TemperatureLogger{T, TI}, sys::MDSys{T}, info::SimulationInfo{T}) where {T<:Number, TI<:Integer}
    if iszero(info.running_step % logger.step) && logger.output
        temp = temperature(sys, info)
        push!(logger.data, temp)

        IO_temperature = open("temperature.txt", "a")
        writedlm(IO_temperature, temp)
        close(IO_temperature)
    end
    return nothing
end

function record!(logger::TrajectionLogger{TI}, sys::MDSys{T}, info::SimulationInfo{T}) where {T<:Number, TI<:Integer}
    if iszero(info.running_step % logger.step) && logger.output
        IO_trajection = open(logger.trajection_file, "a")
        write(IO_trajection, "step = $(info.running_step) \n")
        write(IO_trajection, "id, atom_type, x, y, z, vx, vy, vz", "\n")
        for p_info in info.particle_info
            id = p_info.id
            atom_type = sys.atoms[id].type
            x, y, z = p_info.position
            vx, vy, vz = p_info.velocity
            data = "$id, $atom_type, $x, $y, $z, $vx, $vy, $vz"
            write(IO_trajection, data, "\n")
        end
        writedlm(IO_trajection, "\n")
        close(IO_trajection)
    end
    return nothing
end