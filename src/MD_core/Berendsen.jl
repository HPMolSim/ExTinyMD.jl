struct BerendsenThermoStat{T} <: AbstractThermoStat
    temperature::T
    ν::T
end

BerendsenThermoStat(temperature::T, ν::T) where{T<:Number} = BerendsenThermoStat{T}(temperature, ν)

@inbounds function thermostat_update!(thermostat::BerendsenThermoStat{T}, sys::MDSys{T}, info::SimulationInfo{T}) where {T <: Number}

    λ = sqrt(1.0 + (sys.simulator.dt / thermostat.ν) * ((thermostat.temperature / temperature(sys, info)) - 1.0))

    for p_info in info.particle_info
        p_info.velocity *= λ
    end

    return nothing
end