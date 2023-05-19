

struct VerletProcess{T} <: AbstractSimulator
    dt::T
    thermostat::AbstractThermoStat
end

Verlet(; dt::T, thermostat::AbstractThermoStat = NoThermoStat()) where T = VerletProcess{T}(dt, thermostat)
