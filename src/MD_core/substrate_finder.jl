export SubNeighborFinder, update_finder!

mutable struct SubNeighborFinder{T, TI} <: AbstractNeighborFinder
    cutoff::T
    update_steps::TI
    sub_down::T
    sub_up::T
    up_neighbor::Vector{TI}
    down_neighbor::Vector{TI}
end

function SubNeighborFinder(cutoff::T, info::SimulationInfo{T}, sub_down::T, sub_up::T; update_steps::TI = 100) where{T, TI}
    up_neighbor = Vector{T}()
    down_neighbor = Vector{T}()
    for p_info in info.particle_info
        z_i = p_info.position[3]
        if zero(T) < z_i - sub_down < cutoff
            push!(down_neighbor, p_info.id)
        end
        if zero(T) < sub_up - z_i < cutoff
            push!(up_neighbor, p_info.id)
        end
    end
    return SubNeighborFinder{T, TI}(cutoff, update_steps, sub_down, sub_up, up_neighbor, down_neighbor)
end

function update_finder!(neighborfinder::SubNeighborFinder{T, TI}, info::SimulationInfo{T}) where {T<:Number, TI<:Integer}
    if isone(info.running_step % neighborfinder.update_steps)
        neighborfinder.up_neighbor = Vector{T}()
        neighborfinder.down_neighbor = Vector{T}()
        for p_info in info.particle_info
            z_i = p_info.position[3]
            dz_down = z_i - neighborfinder.sub_down
            if zero(T) < dz_down < neighborfinder.cutoff
                push!(neighborfinder.down_neighbor, p_info.id)
            elseif dz_down < 0
                error("particle out of substrate")
            end
            dz_up = neighborfinder.sub_up - z_i
            if zero(T) < dz_up < neighborfinder.cutoff
                push!(neighborfinder.up_neighbor, p_info.id)
            elseif dz_up < 0
                error("particle out of substrate")
            end
        end
    end
    return nothing
end