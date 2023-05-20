export AndersenThermoStat, thermostat_update!

struct AndersenThermoStat{T} <: AbstractThermoStat
    tempature::T
    ν::T
end

AndersenThermoStat(tempature::T, ν::T) where{T<:Number} = AndersenThermoStat{T}(tempature, ν)

function thermostat_update!(thermostat::AndersenThermoStat{T}, sys::MDSys{T}, info::SimulationInfo{T}) where T <: Number
    sampling = rand(T, sys.n_atoms)
    for i in 1:sys.n_atoms
        if sampling[i] < thermostat.ν * sys.simulator.dt
            d = Normal{T}(zero(T), sqrt(thermostat.tempature / sys.atoms[i].mass))
            info.velcoity[i] = Point(rand(d), rand(d), rand(d))
        end
    end
    return nothing
end