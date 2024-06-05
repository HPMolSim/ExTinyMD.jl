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

function sampling(α, L, n)
    k_set = []
    for i in 0:2*L
        k_vector = generate_k_vector(α,L)
        push!(k_set, k_vector)
    end
    k_set = collect(k_set)
    z = sum(exp(-(norm(k)^2) / (4*α)) for k in k_set)
    prob = [exp(-(norm(k)^2) / (4*α)) / z for k in k_set]
    samples = sample(k_set, weights(prob), n)

    return samples
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

function calculate_Fi(i::Int, p::Int, L::Float64, α::Float64, charges::Vector{Float64}, positions::Matrix{Tuple{Float64, Float64, Float64}}, rho_k::Vector{Complex{Float64}}, samples::Vector)
    V = L^3
    Fi = zeros(Float64, 3)
    qi = charges[i]
    ri = [positions[i]...]

    for j in 1:p
        k2_ell = norm(samples[j])^2
        exp_term = exp(-1im * dot(samples[j], ri))
        imag_part = imag(exp_term * rho_k[j])
        Fi += - (4 * π * samples[j] * qi) / (V * k2_ell) * imag_part
    end

    return Fi
end

function update_rho_k(p::Int, charges::Vector{Float64}, positions::Matrix{Tuple{Float64, Float64, Float64}}, samples::Vector{Any})
    n_atoms = length(charges)
    rho_k = zeros(Complex{Float64}, p)
    for i in 1:p
        rho_k[i] = sum(charges[j] * exp(1im * dot(samples[i], [positions[j]...])) for j in 1:n_atoms)
    end    
    return rho_k
end

function calculate_G(r::Float64, α::Float64)
    term1 = erfc(sqrt(α) * r) / r^2
    term2 = (2 * sqrt(α) * exp(-α * r^2)) / (sqrt(π) * r)
    return term1 + term2
end

function calculate_Fi_short(charges::Vector{Float64}, positions::Matrix{Point{3, Float64}}, neighborfinder::AbstractNeighborFinder, α::Float64, L::Float64)
    neighbor_list = neighborfinder.neighbor_list

    n_atoms = length(charges)
    force_short = [Point(zero(Float64), zero(Float64), zero(Float64)) for _=1:n_atoms]
    boundary = Boundary(L, (1, 1, 1))

    for (i, j, ρ) in neighbor_list
        coord_1, coord_2, r_sq = position_check3D(positions[i], positions[j], boundary, L)
        if iszero(r_sq)
            nothing
        else
            q_1 = charges[i]
            q_2 = charges[j]
            F_ij = calculate_Fs_pair(q_1, q_2, α, coord_1, coord_2)
            force_short[i] += F_ij
            force_short[j] -= F_ij
        end
    end

    return force_short
end

function calculate_Fs_pair(q1::Float64, q2::Float64, α::Float64, coord_1::Point{3, Float64}, coord_2::Point{3, Float64})::Point{3, Float64}
    rij = coord_2 - coord_1
    r = norm(rij)
    G_r = calculate_G(r, α)
    F_ij = -q1 * q2 * G_r * rij / r
    return F_ij
end


function update_acceleration!(interaction::RBEInteractions{T}, neighborfinder::T_NIEGHBER, sys::MDSys{T}, info::SimulationInfo{T}) where {T<:Number, T_NIEGHBER<:AbstractNeighborFinder}
    n_atoms = sys.n_atoms
    L = sys.boundary.length[1]
    α = interaction.α
    p = interaction.p
    samples = interaction.samples

    charges = [sys.atoms[i].charge for i in 1:n_atoms]
    positions = [Point{3, T}(info.particle_info[i].position) for i in 1:n_atoms]

    update_finder!(neighborfinder, info)
    rho_k = update_rho_k(p, charges, positions, samples)

    for i in 1:n_atoms
        Fi = calculate_Fi(i, p, L, α, charges, positions, rho_k, samples)
        Fi2 = calculate_Fi_short(charges, positions, neighborfinder, α, L)
        info.particle_info[i].acceleration += Point{3, T}(Tuple((Fi + Fi2) / sys.atoms[i].mass))
    end

    return nothing
end

