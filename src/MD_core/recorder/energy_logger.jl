struct EnergyLogger{TI} <: AbstractLogger 
    step::TI 
    energy_file::String
    output::Bool
    interactions::Vector{Tuple{AbstractInteraction, AbstractNeighborFinder}}
end

Base.show(io::IO, logger::EnergyLogger) = print(io, "EnergyLogger")

function EnergyLogger(step::TI; interactions::Vector{Tuple{T_interaction, T_neighbor}}, energy_names::Vector{String}, energy_file::String = "energy.txt", output::Bool = true) where{TI<:Integer, T_interaction<:AbstractInteraction, T_neighbor<:AbstractNeighborFinder}
    @assert length(interactions) == length(energy_names)
    if output
        # if there is already a file used to store the result, refresh them
        f = open(energy_file, "w")
        names = "step, Ek, " * join(energy_names, ", ")
        write(f, names, "\n")
        close(f)
    end
    return EnergyLogger{TI}(step, energy_file, output, interactions)
end

@inbounds function kinetic_energy(sys::MDSys{T}, info::SimulationInfo{T}) where {T<:Number}
    Ek = zero(T)
    for p_info in info.particle_info
        id = p_info.id
        Ek += 0.5 * sys.atoms[id].mass * sum(abs2, p_info.velocity)
    end
    return Ek
end

function record!(logger::EnergyLogger{TI}, sys::MDSys{T}, info::SimulationInfo{T}) where {T<:Number, TI<:Integer}
    if iszero(info.running_step % logger.step) && logger.output
        IO_energy = open(logger.energy_file, "a")
        step = info.running_step
        Ek = kinetic_energy(sys, info)
        write(IO_energy, "$step, $Ek")
        for i in 1:length(logger.interactions)
            interaction, neighborfinder = logger.interactions[i]
            E = energy(interaction, neighborfinder, sys, info)
            write(IO_energy, ", $E")
        end
        write(IO_energy, "\n")
        close(IO_energy)
    end
    return nothing
end