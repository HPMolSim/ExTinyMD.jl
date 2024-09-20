# notice that this function only support 3D systems
function random_position(n_atoms::TI, place::NTuple{6, T}, boundary::Boundary{T}; min_r=zero(T), max_attempts::TI=100) where {T, TI<:Integer}
    atoms_coords = Vector{Point{3, T}}()
    
    if iszero(min_r)
        for _ in 1:n_atoms
            push!(atoms_coords, random_coords(place::NTuple{6, T}))
        end
    else
        volume = (place[2] - place[1]) * (place[4] - place[3]) * (place[6] - place[5])
        max_atoms = volume / (4/3 * π * min_r^3)
        if n_atoms > max_atoms
            error("Too much atoms!")
        end
        min_r_sq = min_r^2
        failed_attempts = 0
        push!(atoms_coords, random_coords(place::NTuple{6, T}))

        while length(atoms_coords) < n_atoms
            new_coords = random_coords(place::NTuple{6, T})
            okay = true
            for old_coords in atoms_coords
                min_dist_sq = distance_boundary(new_coords, old_coords, boundary)
                if min_dist_sq < min_r_sq
                    failed_attempts += 1
                    okay = false
                    break
                end
            end
            if okay
                push!(atoms_coords, new_coords)
                failed_attempts = 0
            elseif failed_attempts >= max_attempts
                error("failed to place atoms")
            end
        end
    end

    return atoms_coords
end


function random_coords(place::NTuple{6, T}) where T
    return Point((place[2] - place[1]) * rand(T), (place[4] - place[3]) * rand(T), (place[6] - place[5]) * rand(T)) + Point(place[1], place[3], place[5])
end

@inbounds function distance_boundary(coord_1::Point{3, T}, coord_2::Point{3, T}, boundary::Boundary{T}) where T
    min_dist_sq = dist2(coord_1, coord_2)

    for mx = -boundary.period[1]:boundary.period[1]
        for my = -boundary.period[2]:boundary.period[2]
            for mz = -boundary.period[3]:boundary.period[3]
                new_dist_sq = dist2(coord_1 + Point(mx * boundary.length[1], my * boundary.length[2], mz * boundary.length[3]), coord_2)
                if new_dist_sq < min_dist_sq
                    min_dist_sq = new_dist_sq
                end
            end
        end
    end
    return min_dist_sq
end


function random_velocity(;temp::T, atoms::Vector{Atom{T}}, rng=Random.GLOBAL_RNG) where T
    atoms_velocity = Vector{Point{3, T}}()

    for i in 1:length(atoms)
        atom = atoms[i]
        σ = sqrt(temp / atom.mass)
        push!(atoms_velocity, Point(rand(rng, Normal(zero(T), σ)), rand(rng, Normal(zero(T), σ)), rand(rng, Normal(zero(T), σ))))
    end

    return atoms_velocity
end

function SimulationInfo(n_atoms::Int, atoms::Vector{Atom{T}}, place::NTuple{6, T}, boundary::Boundary{T}; min_r=zero(T), max_attempts::Int=100, rng=Random.GLOBAL_RNG, temp::T = 1.0) where {T}
    positions = random_position(n_atoms, place, boundary; min_r = min_r, max_attempts = max_attempts)
    velocities = random_velocity(;temp = temp, atoms = atoms, rng=rng)
    accelerations = [Point(zero(T), zero(T), zero(T)) for _=1:n_atoms]
    particle_info = [PatricleInfo(id, positions[id], velocities[id], accelerations[id]) for id in 1:n_atoms]

    dict_vec = Vector{Tuple{Int, Int}}()
    for i in 1:n_atoms
        id = particle_info[i].id
        push!(dict_vec, (id, i))
    end
    id_dict = Dict(dict_vec)

    return SimulationInfo{T}(zero(Int64), particle_info, id_dict)
end

function create_atoms(atoms_types::Vector{Tuple{Int, Atom{T}}}) where{T}
    atoms = Vector{Atom{T}}()
    for atom_type in atoms_types
        num, atom = atom_type
        for _=1:num
            push!(atoms, atom)
        end
    end
    return atoms
end