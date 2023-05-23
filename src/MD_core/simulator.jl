export simulate!

function simulate!(simulator::T_SIMULATOR, sys::MDSys{T}, info::SimulationInfo{T}, total_steps::TI) where {T <: Number, T_SIMULATOR<:AbstractSimulator, TI<:Integer}

    for step in 1:total_steps
        info_update!(simulator, sys, info)
        for logger in sys.loggers
            record!(logger, sys, info)
        end
    end
    return nothing
end