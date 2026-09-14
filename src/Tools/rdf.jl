"""
    z_hist(z, bin, z_min, z_max, Lx, Ly, steps) -> (z_list, density)

Histogram the z-coordinates in `z` into `bin - 1` bins spanning `[z_min,
z_max]`, normalised into a number density by the in-plane area `Lx * Ly`, the
bin width, and the number of sampled `steps`. Returns the bin centres
`z_list` alongside the density.
"""
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

"""
    hist_init(N, bin, r_c) -> (hist, volume, r, dr)

Allocate a zeroed radial histogram of `bin` bins spanning `(0, r_c]`, together
with the bin radii `r`, bin width `dr = r_c / bin`, and the spherical-shell
`volume` of each bin (for normalising into a radial distribution function).
`N` is accepted for interface symmetry with the sampling loop but does not
affect the allocation.
"""
function hist_init(N::Integer, bin::Integer, r_c::T) where{T<:Number}
    dr = r_c / bin
    hist = zeros(T, bin)
    r = [dr * i for i in 1:bin]
    volume = 4π .* r.^2 * dr
    return hist, volume, r, dr
end

"""
    distance_hist!(hist, n_list, d)

Bin the pair distances in neighbour list `n_list` (as `(i, j, r)` triples)
into `hist` using bin width `d`, accumulating in place across calls.
"""
function distance_hist!(hist::Vector{T}, n_list::Vector{Tuple{Int64, Int64, T}}, d::T) where T
    for (i, j, r) in n_list
        bin = Int64(round(r / d)) + 1
        hist[bin] += 1
    end    
    return nothing
end