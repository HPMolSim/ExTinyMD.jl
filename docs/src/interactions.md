```@meta
CurrentModule = ExTinyMD
```

# Interactions

An interaction is any object paired with a neighbour finder in `MDSys`'s
`interactions` list. Each step, `simulate!` calls
[`update_acceleration!`](@ref) on every `(interaction, finder)` pair to
accumulate forces, and a logger may separately call [`energy`](@ref) to record
that interaction's contribution.

## Built-in interactions

```@docs
LennardJones
SubLennardJones
ExternalField
NoInteraction
```

The electrostatic interactions (`Ewald3D`, `Ewald2D`, `ICMEwald2D`,
`ICMEwald3D`) are documented separately on the [Electrostatics](@ref) page —
they share the same `update_acceleration!`/`energy` contract described below,
but also work standalone, without any `MDSys`.

## The interaction contract

```@docs
update_acceleration!
energy
```

To add a new interaction, define a struct (optionally `<: AbstractInteraction`,
though this is not required — dispatch only needs the methods below) and
implement:

```julia
function ExTinyMD.update_acceleration!(
    interaction::YourInteraction,
    neighborfinder::YourFinder,
    sys::MDSys{T},
    info::SimulationInfo{T},
) where {T<:Number}
    update_finder!(neighborfinder, info)
    # accumulate into info.particle_info[i].acceleration for every particle i,
    # using neighborfinder.neighbor_list (or neighborfinder.up_neighbor /
    # down_neighbor for a SubNeighborFinder-style finder) to find pairs.
    return nothing
end
```

`update_acceleration!` **accumulates** — it must not overwrite
`info.particle_info[i].acceleration`, since `VerletProcess` zeroes it once per
step before calling every interaction in turn.

If the interaction should also be logged with [`EnergyLogger`](@ref), add a
matching method for `energy`:

```julia
function ExTinyMD.energy(
    interaction::YourInteraction,
    neighborfinder::YourFinder,
    sys::MDSys{T},
    info::SimulationInfo{T},
) where {T<:Number}
    # return the total energy, without mutating anything
end
```

[`LennardJones`](@ref) in `src/interactions/lennard_jones.jl` is the simplest
complete example of both methods; the electrostatics adapter in
`src/interactions/electrostatics/adapter.jl` is a more involved one, showing
how to bridge to a framework-free core API (see [Electrostatics](@ref)).
