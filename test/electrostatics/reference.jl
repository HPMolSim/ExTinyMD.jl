using StaticArrays

"""
Direct lattice sum of the Coulomb energy, periodic in all three axes.
Deliberately naive and O(N² (2n_shell+1)³) — this is the oracle, so it is written
for obviousness, not speed. Converges slowly and only conditionally; use a neutral
configuration and a generous `n_shell`.
"""
function naive_energy_3D(poses, charges, L::NTuple{3,T}, n_shell::Int; ϵ::T = one(T)) where T
    E = zero(T)
    N = length(charges)
    for i in 1:N, j in 1:N
        qq = charges[i] * charges[j]
        for mx in -n_shell:n_shell, my in -n_shell:n_shell, mz in -n_shell:n_shell
            # skip only the i==j self term in the home cell; images of self do count
            (i == j && mx == 0 && my == 0 && mz == 0) && continue
            dx = poses[i][1] - poses[j][1] - mx * L[1]
            dy = poses[i][2] - poses[j][2] - my * L[2]
            dz = poses[i][3] - poses[j][3] - mz * L[3]
            E += qq / sqrt(dx^2 + dy^2 + dz^2)
        end
    end
    return E / (2 * 4π * ϵ)
end

"Direct lattice sum, periodic in x and y only (quasi-2D slab)."
function naive_energy_Q2D(poses, charges, L::NTuple{3,T}, n_shell::Int; ϵ::T = one(T)) where T
    E = zero(T)
    N = length(charges)
    for i in 1:N, j in 1:N
        qq = charges[i] * charges[j]
        for mx in -n_shell:n_shell, my in -n_shell:n_shell
            (i == j && mx == 0 && my == 0) && continue
            dx = poses[i][1] - poses[j][1] - mx * L[1]
            dy = poses[i][2] - poses[j][2] - my * L[2]
            dz = poses[i][3] - poses[j][3]
            E += qq / sqrt(dx^2 + dy^2 + dz^2)
        end
    end
    return E / (2 * 4π * ϵ)
end

"Central finite difference of `f(poses)` w.r.t. component `d` of particle `i`."
function fd_gradient(f, poses::Vector{SVector{3,T}}, i::Int, d::Int;
                     h::T = cbrt(eps(T))) where T
    shift = SVector{3,T}(ntuple(k -> k == d ? h : zero(T), 3))
    p = copy(poses)
    p[i] = poses[i] + shift
    fp = f(p)
    p[i] = poses[i] - shift
    fm = f(p)
    return (fp - fm) / (2h)
end

"""
Rock-salt (NaCl) lattice: `n_cells`³ conventional cells of edge `a`, two
interpenetrating FCC sublattices of opposite charge. Returns AoS positions,
charges, and the periodic box.
"""
function nacl_lattice(n_cells::Int, a::T) where T
    poses = SVector{3,T}[]
    charges = T[]
    h = a / 2
    for ix in 0:(2n_cells - 1), iy in 0:(2n_cells - 1), iz in 0:(2n_cells - 1)
        push!(poses, SVector{3,T}(ix * h, iy * h, iz * h))
        push!(charges, iseven(ix + iy + iz) ? one(T) : -one(T))
    end
    L = (T(n_cells) * a, T(n_cells) * a, T(n_cells) * a)
    return poses, charges, L
end
