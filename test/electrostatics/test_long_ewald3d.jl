@testset "Ewald3D total energy against direct lattice sum" begin
    # NaCl: neutral, symmetric, and with a known answer independent of our code
    a = 2.0
    poses, charges, L = nacl_lattice(2, a)   # 64 ions
    n = length(charges)

    α, s = 2.1, 4.0      # L = (4,4,4) so r_c must be < 2; s/α = 1.90
    short = EwaldShort(n, L; α = α, s = s)
    long  = Ewald3DLong(n, L; α = α, s = s)
    E_ewald = short_energy(short, poses, charges) + long_energy(long, poses, charges)

    # 8π, not 4π — see the oracle's Madelung testset for the derivation of the 2.
    M = -E_ewald * 8π * (a / 2) / n
    @test isapprox(M, 1.7475645946, rtol = 1e-4)
end

@testset "Ewald3D energy is independent of α" begin
    # The split point is arbitrary: short+long must be invariant. This is the
    # sharpest internal check on the relative normalisation of the two parts.
    Random.seed!(20260916)
    n = 20
    L = (12.0, 12.0, 12.0)
    poses = [SVector(rand() * L[1], rand() * L[2], rand() * L[3]) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]

    function total(α, s)
        short = EwaldShort(n, L; α = α, s = s)
        long  = Ewald3DLong(n, L; α = α, s = s)
        return short_energy(short, poses, charges) + long_energy(long, poses, charges)
    end

    # r_c = s/α = 5.71, 5.00, 4.44 — all < L/2 = 6. s = 4 caps accuracy near 1e-7,
    # so rtol is 1e-5 rather than 1e-6.
    E_ref = total(0.7, 4.0)
    @test isapprox(total(0.8, 4.0), E_ref, rtol = 1e-5)
    @test isapprox(total(0.9, 4.0), E_ref, rtol = 1e-5)
end

@testset "Ewald3D force matches -grad(energy)" begin
    Random.seed!(20260917)
    n = 10
    L = (12.0, 12.0, 12.0)
    poses = [SVector(rand() * L[1], rand() * L[2], rand() * L[3]) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]
    α, s = 0.8, 3.5      # r_c = 4.375 < 6
    short = EwaldShort(n, L; α = α, s = s)
    long  = Ewald3DLong(n, L; α = α, s = s)

    F = [zero(SVector{3,Float64}) for _ in 1:n]
    short_force!(F, short, poses, charges)
    long_force!(F, long, poses, charges)

    f = p -> short_energy(short, p, charges) + long_energy(long, p, charges)
    for i in 1:n, d in 1:3
        @test isapprox(F[i][d], -fd_gradient(f, poses, i, d; h = 1e-5),
                       rtol = 1e-4, atol = 1e-8)
    end
end

@testset "Ewald3D surface term responds to ϵ_inf" begin
    # A configuration with a net dipole: conducting (Inf) and vacuum (1.0) boundaries
    # must disagree, and the conducting case must drop the dipole term entirely.
    poses = [SVector(1.0, 4.0, 4.0), SVector(7.0, 4.0, 4.0)]
    charges = [1.0, -1.0]
    L = (8.0, 8.0, 8.0)
    l_cond = Ewald3DLong(2, L; α = 0.8, s = 3.0, ϵ_inf = Inf)   # r_c = 3.75 < 4
    l_vac  = Ewald3DLong(2, L; α = 0.8, s = 3.0, ϵ_inf = 1.0)
    @test !isapprox(long_energy(l_cond, poses, charges), long_energy(l_vac, poses, charges))

    # the difference is exactly the dipole term |P|²/(2Vϵ(2ϵ_inf+1))
    P = sum(charges[i] * poses[i] for i in 1:2)
    V = prod(L)
    @test isapprox(long_energy(l_vac, poses, charges) - long_energy(l_cond, poses, charges),
                   sum(abs2, P) / (2 * V * 1.0 * 3.0), rtol = 1e-10)
end

@testset "Ewald3DLong preserves Float32" begin
    L = (8.0f0, 8.0f0, 8.0f0)
    poses = [SVector(1.0f0, 2.0f0, 3.0f0), SVector(5.0f0, 6.0f0, 7.0f0)]
    charges = [1.0f0, -1.0f0]
    long = Ewald3DLong(2, L; α = 0.8f0, s = 3.0f0)   # r_c = 3.75 < 4
    @test long_energy(long, poses, charges) isa Float32
end

@testset "coulomb_energy/coulomb_force preserve Float32 through the composite" begin
    # FIX 2 regression. `long_energy` alone returning Float32 (the testset above)
    # passed even while the composite silently promoted to Float64 through the
    # `4π` literals in short.jl and icm.jl — that earlier test certified only the
    # one file that was right. Check the full plans that a user actually calls.
    L = (8.0f0, 8.0f0, 8.0f0)
    poses = [SVector(1.0f0, 2.0f0, 3.0f0), SVector(5.0f0, 6.0f0, 7.0f0)]
    charges = [1.0f0, -1.0f0]

    ewald3d = Ewald3D(2, L; α = 0.8f0, s = 3.0f0)   # r_c = 3.75 < min(L)/2 = 4
    @test coulomb_energy(ewald3d, poses, charges) isa Float32
    F3d = coulomb_force(ewald3d, poses, charges)
    @test F3d isa Vector{SVector{3,Float32}}

    # r_c = 3.75 < min(L[1], L[2])/2 = 4, as required for ICM/Ewald2D
    icm = ICMEwald2D(2, L; α = 0.8f0, s = 3.0f0, γ = (0.3f0, 0.3f0), N_image = 2)
    @test coulomb_energy(icm, poses, charges) isa Float32
    Ficm = coulomb_force(icm, poses, charges)
    @test Ficm isa Vector{SVector{3,Float32}}
end
