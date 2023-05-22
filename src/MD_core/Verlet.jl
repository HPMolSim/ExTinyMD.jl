export VerletProcess, info_update!

struct VerletProcess{T_DT, T_THERMO} <: AbstractSimulator
    dt::T_DT
    thermostat::T_THERMO
end

VerletProcess(; dt::T_DT, thermostat::T_THERMO = NoThermoStat()) where {T_DT<:Number, T_THERMO<:AbstractThermoStat} = VerletProcess{T_DT, T_THERMO}(dt, thermostat)


function info_update!(simulator::VerletProcess{T,T_THERMO}, sys::MDSys{T}, info::SimulationInfo{T}) where {T <: Number, T_THERMO<:AbstractThermoStat}

    dt = simulator.dt

    #cleanup the acceleration
    info.running_step += Int64(1)
    info.acceleration .*= zero(T)
    
    # step one v(t + dt/2) = v(t - dt/2) + a(t) dt
    for (interaction, neighborfinder) in sys.interactions
        update_acceleration!(interaction, neighborfinder, sys, info)
    end

    info.velcoity .+= info.acceleration .* dt

    # step two x(t + dt) = x(t) + v(t + dt/2) dt
    info.coords .+= info.velcoity .* dt

    thermostat_update!(simulator.thermostat, sys, info)

    return nothing
end