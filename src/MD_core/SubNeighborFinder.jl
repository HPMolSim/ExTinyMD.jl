export SubNeighborFinder, update_finder!

mutable struct SubNeighborFinder{T, TI}
    cutoff::T
    update_steps::TI
    sub_down::T
    sub_up::T
    up_neighbor::Vector{TI}
    down_neighbor::Vector{TI}
end

function SubNeighborFinder(cutoff::T, coords::Vector{Point{3, T}}, sub_down::T, sub_up::T; update_steps::TI = 100) where{T, TI}
    n_atoms = length(coords)
    up_neighbor = Vector{T}()
    down_neighbor = Vector{T}()
    for i in 1:n_atoms
        z_i = coords[i][3]
        if zero(T) < z_i - sub_down < cutoff
            push!(down_neighbor, i)
        end
        if zero(T) < sub_up - z_i < cutoff
            push!(up_neighbor, i)
        end
    end
    return SubNeighborFinder{T, TI}(cutoff, update_steps, sub_down, sub_up, up_neighbor, down_neighbor)
end

function update_finder!(neighborfinder::SubNeighborFinder{T, TI}, info::SimulationInfo{T}) where {T<:Number, TI<:Integer}
    n_atoms = length(info.coords)
    if isone(info.running_step % neighborfinder.update_steps)
        neighborfinder.up_neighbor = Vector{T}()
        neighborfinder.down_neighbor = Vector{T}()
        for i in 1:n_atoms
            z_i = info.coords[i][3]
            if zero(T) < z_i - neighborfinder.sub_down < neighborfinder.cutoff
                push!(neighborfinder.down_neighbor, i)
            end
            if zero(T) < neighborfinder.sub_up - z_i < neighborfinder.cutoff
                push!(neighborfinder.up_neighbor, i)
            end
        end
    end
    return nothing
end