struct RBEInteractions{T} <: AbstractInteraction
    α::T
    p::Int
    k_set::Vector{}
    prob::Vector{T}
    V::Float64
    cutoff::Float64
end

Base.show(io::IO, interaction::RBEInteractions{}) = print(io, "RBEInteraction with α = $(interaction.α), p = $(interaction.p), k_set = $(interaction.k_set), prob = $(interaction.prob), V = $(interaction.V), cutoff = $(interaction.cutoff)")

LinearAlgebra.norm(p::Point) = sqrt(norm2(p))
norm2(p::Point) = sum(abs2, p.coo)
norm2(p::AbstractVector) = sum(abs2, p)
Base.zero(p::Point{N, T}) where {N, T} = zero(Point{N, T})
Base.zero(::Type{Point{N, T}}) where {N, T} = Point(ntuple(_->zero(T), N))
function initialize_complex_vector(n::Int)::Vector{Complex{Float64}}
    return fill(zero(Complex{Float64}), n)
end

function RBEInteractions(α::T, L::T, p::Int) where T
    α, p, k_set, prob, V, cutoff = sampling(α, L, p)
    return RBEInteractions{T}(α, p, k_set, prob, V, cutoff)
end

function acceptance_probability(m_star::Float64, α::Float64, L::Float64)
    if m_star == 0
        return erf(1 / (2 * sqrt(α * L^2 / π^2)))
    else
        sqrt_val = sqrt(α * L^2 / π^2)
        return 0.5 * (erf((abs(m_star) + 0.5) / sqrt_val) - erf((abs(m_star) - 0.5) / sqrt_val))
    end
end

function mh_sample(α::Float64, L::Float64)
    sample = 0.0
    while sample == 0.0
        x_star = rand(Normal(0, sqrt(α * L^2 / (2 * π^2))))
        m_star = Float64(round(Int, x_star))
        q = acceptance_probability(m_star, α, L)
        if rand() < q && m_star !=0
            sample =  m_star
            break
        end
    end
    return sample
end

function sampling(α, L, p)
    k_set = []
    V = L^3
    cutoff = 3.5
    for i in 0:2*L
        k_vector = generate_k_vector(α,L)
        push!(k_set, k_vector)
    end
    k_set = collect(k_set)
    z = sum(exp(-(norm(k)^2) / (4*α)) for k in k_set)
    prob = [exp(-(norm(k)^2) / (4*α)) / z for k in k_set]
    return α, p, k_set, prob, V, cutoff
end

function calculate_H(α, L)
  
    const_part = sqrt(α * L^2 / π)
    
    exp_term(m, α, L) = exp(-α * m^2 * L^2)
    
    sum_terms = 0.0
    for m in -1:1
        sum_terms += exp_term(m, α, L)
    end
    
    H = const_part * sum_terms
    
    return H
end
 
function calculate_S(α::Float64, L::Float64)
    H = calculate_H(α, L)
    return H^3 - 1
end

function calculate_probability(k::Vector{Float64}, α::Float64, S::Float64)
    k2 = sum(k .^ 2)
    return exp(-k2 / (4 * α)) / S
end

function generate_k_vector(α::Float64, L::Float64) 
    kx = mh_sample(α, L)
    ky = mh_sample(α, L)
    kz = mh_sample(α, L)
    return [kx, ky, kz]
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

function update_rho_k(n_atoms::Int, p, samples, sys, info,)
    rho_k = initialize_complex_vector(p)
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

function calculate_F_short(coord_1,coord_2, α, i, j, sys::MDSys{T}) where T<:Number
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
    samples = map(x->Point(T.(x)...), sample(interaction.k_set, weights(interaction.prob), p, replace = false))
    rho_k = update_rho_k(n_atoms, p, samples, sys, info)
    update_finder!(neighborfinder, info)
    for i in 1:n_atoms
        force_long = calculate_F_long(i, V, p, rho_k, info, samples, sys)
        info.particle_info[i].acceleration += Point(force_long...) / sys.atoms[i].mass
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