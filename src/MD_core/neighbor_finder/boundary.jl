function Boundary(L::NTuple{3, T}, period_set::NTuple{3, Char}) where T
    if period_set == ('p', 'p', 'p')
        return Boundary(L, (1, 1, 1))
    elseif period_set == ('p', 'p', 'f')
        return Boundary(L, (1, 1, 0))
    elseif period_set == ('p', 'f', 'f')
        return Boundary(L, (1, 0, 0))
    end
    error("Illegale Input!")
end

function Q2dBoundary(Lx::T, Ly::T, Lz::T) where T
    return Boundary((Lx, Ly, Lz), (1, 1, 0))
end

function CubicBoundary(L::T) where T
    return Boundary((L, L, L), (1, 1, 1))
end

function BoundaryCheck!(simulation_info::SimulationInfo, boundary::Boundary{T}) where T <: Number
    Lx, Ly, Lz = boundary.length
    px, py, pz = boundary.period
    for p_info in simulation_info.particle_info
        p_info.position -= Point(
            ((zero(T) < p_info.position[1] < Lx) || iszero(px)) ? zero(T) : Lx * div(p_info.position[1], Lx, RoundDown),
            ((zero(T) < p_info.position[2] < Ly) || iszero(py)) ? zero(T) : Ly * div(p_info.position[2], Ly, RoundDown),
            ((zero(T) < p_info.position[3] < Lz) || iszero(pz)) ? zero(T) : Lz * div(p_info.position[3], Lz, RoundDown)
        )
    end
    return nothing
end

@inbounds function position_check3D(coord_1::Point{3, T}, coord_2::Point{3, T}, boundary::Boundary{T}, cutoff::T) where T

    for mx = -boundary.period[1]:boundary.period[1]
        for my = -boundary.period[2]:boundary.period[2]
            for mz = -boundary.period[3]:boundary.period[3]
                dist_sq = dist2(coord_1 + Point(mx * boundary.length[1], my * boundary.length[2], mz * boundary.length[3]), coord_2)
                if dist_sq < abs2(cutoff)
                    return (coord_1 + Point(mx * boundary.length[1], my * boundary.length[2], mz * boundary.length[3]), coord_2, dist_sq)
                end
            end
        end
    end
    return (Point(zero(T), zero(T), zero(T)), Point(zero(T), zero(T), zero(T)), zero(T))
end

@inbounds function position_checkQ2D(coord_1::Point{3, T}, coord_2::Point{3, T}, boundary::Boundary{T}, cutoff::T) where T

    for mx = -boundary.period[1]:boundary.period[1]
        for my = -boundary.period[2]:boundary.period[2]
            dist_sq = abs2(coord_1[1] + mx * boundary.length[1] - coord_2[1]) + abs2(coord_1[2] + my * boundary.length[2] - coord_2[2])
            if dist_sq < abs2(cutoff)
                return (coord_1 + Point(mx * boundary.length[1], my * boundary.length[2], zero(T)), coord_2, dist_sq)
            end
        end
    end
    return (Point(zero(T), zero(T), zero(T)), Point(zero(T), zero(T), zero(T)), zero(T))
end