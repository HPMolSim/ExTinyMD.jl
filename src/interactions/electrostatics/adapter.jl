# Bridge between ExTinyMD's MD loop and the framework-free core API.
#
# Index convention: `info.particle_info` is indexed by storage slot, while
# `sys.atoms` is indexed by particle id. Positions are read in slot order and
# charges gathered to match, so core-layer index `i` consistently means "slot i"
# and forces come back in slot order.

"Gather charges in slot order into `buf`."
function gather_charges!(buf::Vector{T}, sys::MDSys{T}, info::SimulationInfo{T}) where {T}
    @inbounds for i in eachindex(info.particle_info)
        buf[i] = sys.atoms[info.particle_info[i].id].charge
    end
    return buf
end

"Gather positions in slot order into `buf`."
function gather_positions!(buf::Vector{SVector{3,T}},
                           info::SimulationInfo{T}) where {T}
    @inbounds for i in eachindex(info.particle_info)
        p = info.particle_info[i].position
        buf[i] = SVector{3,T}(p[1], p[2], p[3])
    end
    return buf
end

const ElectrostaticInteraction = Union{EwaldInteraction, ICM}

# Scratch for the gathered arrays, keyed on the interaction so repeated steps do
# not allocate. Stored on the interaction itself via these accessors.
_pos_scratch(inter::EwaldInteraction) = inter.pos_scratch
_charge_scratch(inter::EwaldInteraction) = inter.charge_scratch
_pos_scratch(inter::ICM) = inter.pos_scratch
_charge_scratch(inter::ICM) = inter.charge_scratch

# A NoNeighborFinder carries no usable list, so fall back to the interaction's
# own cell list in that case.
_finder_list(f::NoNeighborFinder) = nothing
_finder_list(f) = f.neighbor_list

function ExTinyMD.energy(inter::ElectrostaticInteraction, neighborfinder,
                         sys::MDSys{T}, info::SimulationInfo{T}) where {T}
    update_finder!(neighborfinder, info)
    poses   = gather_positions!(_pos_scratch(inter), info)
    charges = gather_charges!(_charge_scratch(inter), sys, info)
    return coulomb_energy(inter, poses, charges;
                          neighbor_list = _finder_list(neighborfinder))
end

function ExTinyMD.update_acceleration!(inter::ElectrostaticInteraction, neighborfinder,
                                       sys::MDSys{T}, info::SimulationInfo{T}) where {T}
    update_finder!(neighborfinder, info)
    poses   = gather_positions!(_pos_scratch(inter), info)
    charges = gather_charges!(_charge_scratch(inter), sys, info)

    F = coulomb_force!(inter.force_buffer, inter, poses, charges;
                       neighbor_list = _finder_list(neighborfinder))

    @inbounds for i in eachindex(info.particle_info)
        m = sys.atoms[info.particle_info[i].id].mass
        f = F[i]
        info.particle_info[i].acceleration += Point(f[1] / m, f[2] / m, f[3] / m)
    end
    return nothing
end
