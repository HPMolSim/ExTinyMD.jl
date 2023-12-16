# nose-hoover chain with M = 1. Follow the algorithm 30/31/32 from Frenkel & Smit's book:
# Understanding Molecular simulaton From Algorithms to Applications.

mutable struct NHVerletProcess{T <: Number} <: AbstractSimulator
    dt::T
    Q::T
    vξ::T
    ξ::T
    Ek::T
    temperature::T
end

NHVerletProcess(; dt::T, Q::T = 1.0, temperature::T) where {T <: Number} = NHVerletProcess{T}(dt, Q, zero(T), zero(T), zero(T), temperature)

@inbounds function Ek(sys::MDSys{T}, info::SimulationInfo{T}) where {T <: Number}
    Ek = zero(T)
    for p_info in info.particle_info
        id = p_info.id
        Ek += 0.5 * sys.atoms[id].mass * sum(abs2, p_info.velocity)
    end
    return Ek
end

function update_nhchain!(simulator::NHVerletProcess{T}, sys::MDSys{T}, info::SimulationInfo{T}) where {T <: Number}
    G1 = (2 * simulator.Ek - 3 * sys.n_atoms * simulator.temperature) / simulator.Q
    simulator.vξ += G1 * simulator.dt / 4
    simulator.ξ += simulator.vξ * simulator.dt / 2
    s = exp(-simulator.vξ * simulator.dt / 2)
    for p_info in info.particle_info
        p_info.velocity *= s
    end
    simulator.Ek *= s^2
    G1 = (2 * simulator.Ek - 3 * sys.n_atoms * simulator.temperature) / simulator.Q
    simulator.vξ += G1 * simulator.dt / 4
    return nothing
end

@inbounds function info_update!(simulator::NHVerletProcess{T}, sys::MDSys{T}, info::SimulationInfo{T}) where {T <: Number}

    dt = simulator.dt

    #cleanup the acceleration
    info.running_step += Int64(1)
    erase_acceleration!(info)

    # init the Ek in the first step
    if info.running_step == 1
        simulator.Ek = Ek(sys, info)
    end

    update_nhchain!(simulator, sys, info)

    update_position!(info, dt / 2)
    BoundaryCheck!(info, sys.boundary)

    for (interaction, finder) in sys.interactions
        update_acceleration!(interaction, finder, sys, info)
    end

    update_velocity!(info, dt)
    update_position!(info, dt / 2)

    simulator.Ek = Ek(sys, info)
    update_nhchain!(simulator, sys, info)

    BoundaryCheck!(info, sys.boundary)
    update_dict!(info)

    return nothing
end

@inbounds function info_update!(simulator::NHVerletProcess{T}, sys::MDSys{T}, info::SimulationInfo{T}, total_steps::TI) where {T <: Number, TI<:Integer}
    simulator.Ek = Ek(sys, info)
    for step in 1:total_steps
        info_update!(simulator, sys, info)
    end
    return nothing
end