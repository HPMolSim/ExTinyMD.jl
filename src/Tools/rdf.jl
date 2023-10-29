function z_hist(z::Vector{T}, bin::Int, z_min::T, z_max::T, Lx::T, Ly::T, steps::Int) where{T}
	L = z_max - z_min
	dL = L / (bin - 1)

	z_list = [dL / 2 + dL * (i - 1) for i in 1:bin - 1]
	count = zeros(bin - 1)

	for zi in z
		l = zi - z_min
		id_i = Int(l ÷ dL) + 1
		count[id_i] += 1
	end

	return z_list, count / (steps * Lx * Ly * dL)
end

function hist_init(N::Integer, bin::Integer, r_c::T) where{T<:Number}
    dr = r_c / bin
    hist = zeros(T, bin)
    r = [dr * i for i in 1:bin]
    volume = 4π .* r.^2 * dr
    return hist, volume, r, dr
end

function distance_hist!(hist::Vector{T}, n_list::Vector{Tuple{Int64, Int64, T}}, d::T) where T
    for (i, j, r) in n_list
        bin = Int64(round(r / d)) + 1
        hist[bin] += 1
    end    
    return nothing
end