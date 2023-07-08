export hist_init, distance_hist!

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