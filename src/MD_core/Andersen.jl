struct AndersenThermoStat{T} <: AbstractThermoStat
    temperature::T
    ν::T
end

AndersenThermoStat(temperature::T, ν::T) where{T<:Number} = AndersenThermoStat{T}(temperature, ν)

@inbounds function thermostat_update!(thermostat::AndersenThermoStat{T}, sys::MDSys{T, T_INTERACTION, T_LOGGER, T_SIMULATOR}, info::SimulationInfo{T}) where {T <: Number, T_INTERACTION, T_LOGGER, T_SIMULATOR}
    for p_info in info.particle_info
        if rand(T) < thermostat.ν * sys.simulator.dt
            d = Normal{T}(zero(T), sqrt(thermostat.temperature / sys.atoms[p_info.id].mass))
            p_info.velocity = Point(rand(d), rand(d), rand(d))
        end
    end
    return nothing
end