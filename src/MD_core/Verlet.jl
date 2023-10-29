struct VerletProcess{T_DT, T_THERMO} <: AbstractSimulator
    dt::T_DT
    thermostat::T_THERMO
end

VerletProcess(; dt::T_DT, thermostat::T_THERMO = NoThermoStat()) where {T_DT<:Number, T_THERMO<:AbstractThermoStat} = VerletProcess{T_DT, T_THERMO}(dt, thermostat)

@inbounds function erase_acceleration!(info::SimulationInfo{T}) where{T}
    for i in 1:length(info.particle_info)
        info.particle_info[i].acceleration = Point(zero(T), zero(T), zero(T))
    end
    return nothing
end

@inbounds function update_velocity!(info::SimulationInfo{T}, dt::T) where{T}
    for i in 1:length(info.particle_info)
        info.particle_info[i].velocity += info.particle_info[i].acceleration * dt
    end
    return nothing
end

@inbounds function update_position!(info::SimulationInfo{T}, dt::T) where{T}
    for i in 1:length(info.particle_info)
        info.particle_info[i].position += info.particle_info[i].velocity * dt
    end
    return nothing
end

@inbounds function update_dict!(info::SimulationInfo{T}) where{T}
    for i in 1:length(info.particle_info)
        info.id_dict[info.particle_info[i].id] = i
    end
    return nothing
end

@inbounds function info_update!(simulator::VerletProcess{T,T_THERMO}, sys::MDSys{T}, info::SimulationInfo{T}) where {T <: Number, T_THERMO<:AbstractThermoStat}

    dt = simulator.dt

    #cleanup the acceleration
    info.running_step += Int64(1)
    erase_acceleration!(info)
    
    # step one v(t + dt/2) = v(t - dt/2) + a(t) dt
    for (interaction, finder) in sys.interactions
        update_acceleration!(interaction, finder, sys, info)
    end

    update_velocity!(info, dt)
    # step two x(t + dt) = x(t) + v(t + dt/2) dt
    update_position!(info, dt)

    BoundaryCheck!(info, sys.boundary)
    update_dict!(info)

    thermostat_update!(simulator.thermostat, sys, info)

    return nothing
end

@inbounds function info_update!(simulator::VerletProcess{T,T_THERMO}, sys::MDSys{T}, info::SimulationInfo{T}, total_steps::TI) where {T <: Number, T_THERMO<:AbstractThermoStat, TI<:Integer}
    for step in 1:total_steps
        info_update!(simulator, sys, info)
    end
    return nothing
end