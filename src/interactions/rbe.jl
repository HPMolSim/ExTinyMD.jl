struct RandomBatchEwald{T} <: AbstractInteraction
    ϵ::T
    ϵ_inf::T # boundary condition at infinity, default is Inf (conductor)
    α::T
    
    r_c::T
    k_c::T
    L::NTuple{3,T}
    k_set::Vector{Tuple{T, T, T, T}}
    weight::Weights{T, T, Vector{T}}
    S::T

    rbe::Bool
    p::Int
    weight_cut::T
end

Base.show(io::IO, interaction::RandomBatchEwald) = print(io, "RandomBatchEwald with ϵ = $(interaction.ϵ), ϵ_inf = $(interaction.ϵ_inf), α = $(interaction.α), r_c = $(interaction.r_c), k_c = $(interaction.k_c), L = $(interaction.L), rbe = $(interaction.rbe), p = $(interaction.p)")

function RandomBatchEwald(s::T, α::T, L::NTuple{3,T}; ϵ::T = one(T), ϵ_inf::T = Inf, rbe::Bool = true, p::Int = 30, weight_cut::T = 0.001) where{T}
    r_c = s / α
    k_c = 2 * α * s
    k_set, weight, S = Ewald_kset(k_c, L, α, weight_cut, rbe)
    return RandomBatchEwald(ϵ, ϵ_inf, α, r_c, k_c, L, k_set, weight, S, rbe, p, weight_cut)
end

function Ewald_kset(k_c::T, L::NTuple{3,T}, α::T, weight_cut::T, rbe::Bool) where{T}
    L_x, L_y, L_z = L

    mx_max = ceil(Int, k_c * L_x / 2π) + 1
    my_max = ceil(Int, k_c * L_y / 2π) + 1
    mz_max = ceil(Int, k_c * L_z / 2π) + 1

    k_set = Vector{NTuple{4, T}}()
    weight = Vector{T}()
    S = zero(T)

    for m_x in - mx_max : mx_max
        for m_y in - my_max : my_max
            for m_z in - mz_max : mz_max
                k_x = m_x * 2π / L_x
                k_y = m_y * 2π / L_y
                k_z = m_z * 2π / L_z
                k = sqrt(k_x^2 + k_y^2 + k_z^2)
                P = exp(-k^2 / (4 * α^2))
                S += P
                if 0 < k <= k_c && (!rbe || (rbe && P > weight_cut))
                    push!(k_set, (k_x, k_y, k_z, k))
                    push!(weight, P)
                end
            end
        end
    end

    return k_set, weights(weight), S
end

#################### energy ####################

function energy(interaction::RandomBatchEwald{T}, neighbor::CellList3D{T}, sys::MDSys{T}, info::SimulationInfo{T}) where{T}
    atoms = sys.atoms
    boundary = sys.boundary
    energy_long = Ewald3D_long_energy(interaction, info, atoms)
    energy_short = Ewald3D_short_energy(interaction, neighbor, info, atoms, boundary)

    @debug "Ewald3D El = $energy_long, Es = $energy_short"

    return energy_long + energy_short
end

function Ewald3D_long_energy(interaction::RandomBatchEwald{T}, info::SimulationInfo{T}, atoms::Vector{Atom{T}}) where{T}
    energy = zero(T)

    if interaction.ϵ_inf != Inf
        energy += Ewald3D_long_energy_k0(interaction, info, atoms)
    end

    if !(interaction.rbe)
        for K in interaction.k_set
            energy += Ewald3D_long_energy_k(K, interaction, info, atoms) * exp(- k[4]^2 / (4 * interaction.α^2))
        end
    else
        for _ in 1:p
            K = sample(interaction.k_set, interaction.weight)
            energy += interation.S / interaction.p * Ewald3D_long_energy_k(K, interaction, info, atoms)
        end
    end
    return energy / (4π * interaction.ϵ)
end

function Ewald3D_long_energy_k(K::Tuple{T, T, T, T}, interaction::RandomBatchEwald{T}, info::SimulationInfo{T}, atoms::Vector{Atom{T}}) where{T}
    k_x, k_y, k_z, k = K
    L_x, L_y, L_z = interaction.L
    α = interaction.α
    n_atoms = length(atoms)

    ρ_k = zero(ComplexF64)
    for j in 1:n_atoms
        x_j, y_j, z_j = info.particle_info[j].position
        qj = atoms[j].charge
        ρ_k += qj * exp(1.0im * (k_x * x_j + k_y * y_j + k_z * z_j))
    end

    sum_k = 2π / (L_x * L_y * L_z) * T(ρ_k' * ρ_k) / k^2

    return sum_k
end

function Ewald3D_long_energy_k0(interaction::RandomBatchEwald{T}, info::SimulationInfo{T}, atoms::Vector{Atom{T}}) where{T}

    L_x, L_y, L_z = interaction.L
    n_atoms = length(atoms)
    P = Point(zero(T), zero(T), zero(T))
    for i in 1:n_atoms
        P += atoms[i].charge * info.particle_info[i].position
    end

    energy_k0 = 2π * dist2(P, Point(zero(T), zero(T), zero(T))) / (2 * interaction.ϵ_inf + one(T)) / (L_x * L_y * L_z)

    return energy_k0
end

function Ewald3D_short_energy(interaction::RandomBatchEwald{T}, neighbor::CellList3D{T}, info::SimulationInfo{T}, atoms::Vector{Atom{T}}, boundary::Boundary) where{T}
    neighbor_list = neighbor.neighbor_list

    energy_short = zero(T)

    for (i, j, ρ) in neighbor_list
        position_i = info.particle_info[i].position
        position_j = info.particle_info[j].position
        coord_1, coord_2, r_sq = position_check3D(position_i, position_j, boundary, interaction.r_c)
        if iszero(r_sq)
            nothing
        else
            q_1 = atoms[i].charge
            q_2 = atoms[j].charge
            energy_short += Ewald3D_Es_pair(q_1, q_2, interaction.α, r_sq)
        end
    end

    for i in 1:interaction.n_atoms
        energy_short += Ewald3D_Es_self(q[i], interaction.α)
    end

    return energy_short / (4π * interaction.ϵ)
end

function Ewald3D_Es_pair(q_1::T, q_2::T, α::T, r_sq::T) where{T}
    return q_1 * q_2 * erfc(α * sqrt(r_sq)) / sqrt(r_sq)
end

function Ewald3D_Es_self(q::T, α::T) where{T}
    return - q^2 * α / sqrt(π)
end

#################### force ####################

function update_acceleration!(interaction::RandomBatchEwald{T}, neighbor::CellList3D{T}, sys::MDSys{T}, info::SimulationInfo{T}) where{T}
    atoms = sys.atoms
    boundary = sys.boundary
    Ewald3D_long_force!(interaction, info, atoms)
    Ewald3D_short_force!(interaction, neighbor, info, atoms, boundary)

    return nothing
end

function Ewald3D_short_force!(interaction::RandomBatchEwald{T}, neighbor::CellList3D{T}, info::SimulationInfo{T}, atoms::Vector{Atom{T}}, boundary::Boundary) where{T}
    neighbor_list = neighbor.neighbor_list

    for (i, j, ρ) in neighbor_list
        coord_1, coord_2, r_sq = position_check3D(info.particle_info[i].position, info.particle_info[j].position, boundary, interaction.r_c)
        if iszero(r_sq)
            nothing
        else
            q_1 = atoms[i].charge
            q_2 = atoms[j].charge
            F_ij = Ewald3D_Fs_pair(q_1, q_2, interaction.α, coord_1, coord_2)
            info.particle_info[i].acceleration += F_ij / (4π * interaction.ϵ) / atoms[i].mass
            info.particle_info[j].acceleration -= F_ij / (4π * interaction.ϵ) / atoms[j].mass
        end
    end

    nothing
end

function Ewald3D_Fs_pair(q_1::T, q_2::T, α::T, coord_1::Point{3, T}, coord_2::Point{3, T}) where{T}
    energy = r -> q_1 * q_2 * erfc(α * r) / r
    Δr = sqrt(dist2(coord_1, coord_2))
    if iszero(Δr)
        return Point(zero(T), zero(T), zero(T))
    else
        force = - ForwardDiff.derivative(energy, Δr)
        return T(force) * (coord_1 - coord_2) / Δr
    end
end

function Ewald3D_long_force!(interaction::RandomBatchEwald{T}, info::SimulationInfo{T}, atoms::Vector{Atom{T}}) where{T}
    if interaction.ϵ_inf != Inf
        Ewald3D_long_force_k0!(interaction, info, atoms)
    end

    if interaction.rbe
        for _ in 1:interaction.p
            K = sample(interaction.k_set, interaction.weight)
            Ewald3D_long_force_k!(K, interaction, info, atoms)
        end
    else
        for K in interaction.k_set
            Ewald3D_long_force_k!(K, interaction, info, atoms)
        end
    end

    nothing
end

function Ewald3D_long_force_k!(K::Tuple{T, T, T, T}, interaction::RandomBatchEwald{T}, info::SimulationInfo{T}, atoms::Vector{Atom{T}}) where{T}
    k_x, k_y, k_z, k = K
    L_x, L_y, L_z = interaction.L
    α = interaction.α
    n_atoms = length(atoms)

    ρ_k = zero(ComplexF64)
    for j in 1:n_atoms
        x_j, y_j, z_j = info.particle_info[j].position
        ρ_k += atoms[j].charge * exp(1.0im * (k_x * x_j + k_y * y_j + k_z * z_j))
    end

    for i in 1:n_atoms
        x_i, y_i, z_i = info.particle_info[i].position
        para = interaction.rbe ? interaction.S / interaction.p : exp(- k[4]^2 / (4 * α^2))
        info.particle_info[i].acceleration += - atoms[i].charge * T(imag(ρ_k * exp( - 1.0im * (k_x * x_i + k_y * y_i + k_z * z_i)))) / k^2 * para /interaction.ϵ * Point(k_x, k_y, k_z) / (L_x * L_y * L_z) / atoms[i].mass
    end

    return nothing
end

function Ewald3D_long_force_k0!(interaction::RandomBatchEwald{T}, info::SimulationInfo{T}, atoms::Vector{Atom{T}}) where{T}
    L_x, L_y, L_z = interaction.L

    P = Point(zero(T), zero(T), zero(T))
    n_atoms = lenght(info)
    for i in 1:n_atoms
        P += atoms[i].charge * info.particle_info[i].position
    end

    para = - one(T) / (2 * interaction.ϵ_inf + one(T)) / (2 * interaction.ϵ * L_x * L_y * L_z)

    for i in 1:interaction.n_atoms
        info.particle_info[i].acceleration += para * T(2) * P * atoms[i].charge / atoms[i].mass
    end

    return nothing
end