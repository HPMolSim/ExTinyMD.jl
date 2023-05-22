export TempartureLogger

mutable struct TempartureLogger{T} <: AbstractLogger
    step::Integer
    data::Vector{T}
end

mutable struct PressureLogger{T} <: AbstractLogger
    step::Integer
    data::Vector{T}
end


TempartureLogger(step::TI) where{TI<:Integer} = TempartureLogger{Float64}(step, Vector{Float64}())
TempartureLogger{T}(step::TI) where{TI<:Integer, T<:Number} = TempartureLogger{T}(step, Vector{T}())

PressureLogger(step::TI) where{TI<:Integer} = PressureLogger{Float64}(step, Vector{Float64}())
PressureLogger{T}(step::TI) where{TI<:Integer, T<:Number} = PressureLogger{T}(step, Vector{T}())


# this function is used to record temparture
function record!(logger::TempartureLogger, sys::MDSys{T, T_INTERACTION, T_LOGGER, T_SIMULATOR}, info::SimulationInfo{T}) where {T, T_INTERACTION, T_LOGGER, T_SIMULATOR}
    if iszero(info.running_step % logger.step)
        Ek = sum([0.5 * sys.atoms[i].mass * sum(abs2, info.velcoity[i]) for i in 1:sys.n_atoms])
        temparture = 2.0 * Ek / (3.0 * sys.n_atoms)

        push!(logger.data, temparture)
    end
    return nothing
end