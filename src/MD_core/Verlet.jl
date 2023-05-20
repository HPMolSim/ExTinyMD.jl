export VerletProcess, info_update!

struct VerletProcess{T} <: AbstractSimulator
    dt::T
    thermostat::AbstractThermoStat
end

VerletProcess(; dt::T, thermostat::AbstractThermoStat = NoThermoStat()) where T = VerletProcess{T}(dt, thermostat)


function info_update!(simulator::VerletProcess{T}, sys::MDSys{T}, info::SimulationInfo{T}) where {T <: Number}

    dt = simulator.dt
    
    # step one v(t + dt/2) = v(t - dt/2) + a(t) dt
    a  = sum.([acceleration(interaction, sys, info) for interaction in sys.interactions])
    info.velcoity .+= a .* dt

    # step two x(t + dt) = x(t) + v(t + dt/2) dt
    info.coords .+= info.velcoity .* dt

    thermostat_update!(simulator.thermostat, info.velcoity)

    return nothing
end