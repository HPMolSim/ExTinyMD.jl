function MSD(data::Vector, boundary::Boundary{T}; s::Int64 = 1) where{T}
    n_atoms = length(data[1])
    displacement = [Point(zero(T), zero(T), zero(T)) for j in 1:n_atoms]
    msd = [[zero(T)] for i in 1:3]

    L = boundary.length
    Lm = minimum(L)

    for i in s:length(data) - 1
        dx, dy, dz = zero(T), zero(T), zero(T)
        for j in 1:n_atoms
            coord_1 = Point(data[i][j][3], data[i][j][4], data[i][j][5])
            coord_2 = Point(data[i + 1][j][3], data[i + 1][j][4], data[i + 1][j][5])
            coord_2, coord_1, dist_sq = position_check3D(coord_2, coord_1, boundary, Lm / 2.1)
            displacement[j] += coord_2 - coord_1
        end

        dx = sum(displacement[i][1]^2 for i in 1:n_atoms)
        dy = sum(displacement[i][2]^2 for i in 1:n_atoms)
        dz = sum(displacement[i][3]^2 for i in 1:n_atoms)

        push!(msd[1], dx / n_atoms)
        push!(msd[2], dy / n_atoms)
        push!(msd[3], dz / n_atoms)
    end
    return msd
end