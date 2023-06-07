export TempartureLogger, record!, TrajectionLogger

mutable struct TempartureLogger{T, TI} <: AbstractLogger
    step::TI
    data::Vector{T}
end

mutable struct PressureLogger{T, TI} <: AbstractLogger
    step::TI
    data::Vector{T}
end

struct TrajectionLogger{TI, T} <: AbstractLogger 
    step::TI 
    position_file::String
    velocity_file::String
    position::Vector{SVector{3, T}}
    velocity::Vector{SVector{3, T}}
end


# TempartureLogger(step::TI) where{TI<:Integer} = TempartureLogger{Float64, TI}(step, Vector{Float64}())

function TempartureLogger(step::TI) where{TI<:Integer}
    f = open("temparture.txt", "w")
    close(f)
    return TempartureLogger{Float64, TI}(step, Vector{Float64}())
end
TempartureLogger{T}(step::TI) where{TI<:Integer, T<:Number} = TempartureLogger{T, TI}(step, Vector{T}())

PressureLogger(step::TI) where{TI<:Integer} = PressureLogger{Float64}(step, Vector{Float64}())
PressureLogger{T}(step::TI) where{TI<:Integer, T<:Number} = PressureLogger{T}(step, Vector{T}())

function TrajectionLogger(info::SimulationInfo{T}, step::TI) where{T<:Number, TI<:Integer}
    position_file = "position.txt"
    velocity_file = "velocity.txt"

    # if there is already a file used to store the result, refresh them
    f = open(position_file, "w")
    close(f)
    f = open(velocity_file, "w")
    close(f)

    position = [SVector{3, T}(xi[1], xi[2], xi[3]) for xi in info.coords]
    velocity = [SVector{3, T}(xi[1], xi[2], xi[3]) for xi in info.velcoity]
    return TrajectionLogger{TI, T}(step, position_file, velocity_file, position, velocity)
end

# this function is used to record temparture
function record!(logger::TempartureLogger{T, TI}, sys::MDSys{T}, info::SimulationInfo{T}) where {T<:Number, TI<:Integer}
    if iszero(info.running_step % logger.step)
        Ek = zero(T)
        for i in 1:sys.n_atoms
            Ek += 0.5 * sys.atoms[i].mass * sum(abs2, info.velcoity[i])
        end
        temparture = T(2.0 * Ek / (3.0 * sys.n_atoms))
        push!(logger.data, temparture)

        IO_tempature = open("temparture.txt", "a")
        writedlm(IO_tempature, temparture)
        close(IO_tempature)
    end
    return nothing
end

function record!(logger::TrajectionLogger{TI}, sys::MDSys{T}, info::SimulationInfo{T}) where {T<:Number, TI<:Integer}
    if iszero(info.running_step % logger.step)
        PointToStaticArrays3D!(info.coords, logger.position)
        PointToStaticArrays3D!(info.velcoity, logger.velocity)

        IO_position = open(logger.position_file, "a")
        writedlm(IO_position, info.running_step)
        writedlm(IO_position, logger.position, ",")
        writedlm(IO_position, "\n")
        close(IO_position)

        IO_velocity = open(logger.velocity_file, "a")
        writedlm(IO_velocity, info.running_step)
        writedlm(IO_velocity, logger.velocity, ",")
        writedlm(IO_velocity, "\n")
        close(IO_velocity)
    end
    return nothing
end