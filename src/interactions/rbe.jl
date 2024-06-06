struct RBEInteractions{T} <: AbstractInteraction
    α::T
    p::Int
    k_set::Vector{}
    prob::Vector{T}
    V::Float64
    cutoff::T
    k_c::T
end
LinearAlgebra.norm(p::Point) = sqrt(norm2(p))
norm2(p::Point) = sum(abs2, p.coo)
norm2(p::AbstractVector) = sum(abs2, p)
Base.zero(p::Point{N, T}) where {N, T} = zero(Point{N, T})
Base.zero(::Type{Point{N, T}}) where {N, T} = Point(ntuple(_->zero(T), N))
Base.zeros(::Type{Complex{Float64}}) = Complex{Float64}(0.0, 0.0)
Base.zeros(::Type{Complex{Float64}}, p::Int) = fill(Complex{Float64}(0.0, 0.0), p)

function RBEInteractions(α::T, L_x::Float64, L_y::Float64, L_z::Float64, p::Int, s::Float64) where T
    k_c = 2s/sqrt(α)
    cutoff = k_c / 2
    S = calculate_S(α, L_x, L_y, L_z)
    V = L_x * L_y * L_z
    k_set = generate_k_set(L_x, L_y, L_z, k_c)
    prob = compute_probabilities(k_set, α, S)
    return RBEInteractions{T}(α, p, k_set, prob, V, cutoff, k_c)
end

function generate_k_set(L_x::Float64, L_y::Float64, L_z::Float64, k_c::Float64)
    m_max_x = floor(Int, k_c * L_x / (2 * π))
    m_max_y = floor(Int, k_c * L_y / (2 * π))
    m_max_z = floor(Int, k_c * L_z / (2 * π))
    k_set = []
    for m1 in -m_max_x:m_max_x
        for m2 in -m_max_y:m_max_y
            for m3 in -m_max_z:m_max_z
                k_vector = [2 * π * m1 / L_x, 2 * π * m2 / L_y, 2 * π * m3 / L_z]
                if norm(k_vector) <= k_c && norm(k_vector) != 0
                    push!(k_set, k_vector)
                end
            end
        end
    end
    return k_set
end

function compute_probabilities(k_set, α::Float64, S::Float64)
    probabilities = []
    for k in k_set
        prob = exp(-norm(k)^2 / (4 * α)) / S
        push!(probabilities, prob)
    end
    return probabilities
end

function calculate_H(α::Float64, L_x::Float64, L_y::Float64, L_z::Float64)
    Hx = sqrt(α * L_x^2 / π) * (1 + 2 * exp(-α * L_x^2))
    Hy = sqrt(α * L_y^2 / π) * (1 + 2 * exp(-α * L_y^2))
    Hz = sqrt(α * L_z^2 / π) * (1 + 2 * exp(-α * L_z^2))
    return Hx * Hy * Hz
end
 
function calculate_S(α::Float64, L_x::Float64, L_y::Float64, L_z::Float64)
    H = calculate_H(α, L_x, L_y, L_z)
    return H^3 - 1
end

function calculate_F_long(i::Int, V, p::Int, rho_k::Vector{Complex{Float64}}, info::SimulationInfo{T}, samples::Vector, sys) where T<:Number
    Fi = zero(eltype(samples))
    for j in 1:p
        k2_ell = norm2(samples[j])
        exp_term = exp(-1im * dot(samples[j], info.particle_info[i].position))
        imag_part = imag(exp_term * rho_k[j])
        Fi += - (4 * π * samples[j] * sys.atoms[i].charge) / (V * k2_ell) * imag_part
    end
    return Fi
end

function update_rho_k(n_atoms::Int, p::Int, samples, sys, info,)
    rho_k = zeros(Complex{Float64}, p)
    for i in 1:p
        rho_k[i] = sum(sys.atoms[j].charge * exp(1im * dot(samples[i], info.particle_info[j].position)) for j in 1:n_atoms)
    end    
    return rho_k
end

function calculate_G(r::Float64, α::Float64)
    term1 = erfc(sqrt(α) * r) / r^2
    term2 = (2 * sqrt(α) * exp(-α * r^2)) / (sqrt(π) * r)
    return term1 + term2
end

function calculate_F_short(coord_1, coord_2, α::Float64, i::Int, j::Int, sys::MDSys{T}) where T<:Number
    q1 = sys.atoms[i].charge
    q2 = sys.atoms[j].charge
    r_ij = coord_2 - coord_1
    r = norm2(r_ij)
    G_r = calculate_G(r, α)
    F_ij = -q1 * q2 * G_r * r_ij / r
    return F_ij
end

function update_acceleration!(interaction::RBEInteractions{T}, neighborfinder::T_NIEGHBER, sys::MDSys{T}, info::SimulationInfo{T}) where {T<:Number ,T_NIEGHBER<:AbstractNeighborFinder}
    n_atoms = sys.n_atoms
    α = interaction.α
    p = interaction.p
    V = interaction.V
    boundary = sys.boundary
    samples = map(x->Point(T.(x)...), sample(interaction.k_set, weights(interaction.prob), p))
    rho_k = update_rho_k(n_atoms, p, samples, sys, info)
    update_finder!(neighborfinder, info)
    for i in 1:n_atoms
        force_long = calculate_F_long(i, V, p, rho_k, info, samples, sys)
        info.particle_info[i].acceleration += force_long / sys.atoms[i].mass
    end

    for (i, j, r) in neighborfinder.neighbor_list
        coord_1, coord_2, dist_sq = position_check3D(info.particle_info[i].position, info.particle_info[j].position, boundary, interaction.cutoff)
        if iszero(dist_sq)
            nothing
        else
            force_vector = calculate_F_short(coord_1, coord_2, α, i, j, sys)
            force_short = Point(force_vector...)
            info.particle_info[i].acceleration += force_short / sys.atoms[i].mass
            info.particle_info[j].acceleration -= force_short / sys.atoms[j].mass
        end
    end
end
