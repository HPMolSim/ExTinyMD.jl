@testset "Ewald2D energy is independent of α" begin
    Random.seed!(20260922)
    n = 16
    L = (6.0, 6.0, 20.0)
    # a slab: charges confined well inside the z extent
    poses = [SVector(rand() * L[1], rand() * L[2], 5.0 + 10.0 * rand()) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]

    total(α, s) = coulomb_energy(Ewald2D(n, L; α = α, s = s), poses, charges)
    # Only L[1], L[2] bound r_c under PeriodicQ2D: r_c = s/α = 2.86, 2.67, 2.50,
    # all < min(Lx,Ly)/2 = 3.
    E_ref = total(1.4, 4.0)
    @test isapprox(total(1.5, 4.0), E_ref, rtol = 1e-5)
    @test isapprox(total(1.6, 4.0), E_ref, rtol = 1e-5)
end

@testset "Ewald2D energy against quasi-2D direct sum" begin
    Random.seed!(20260923)
    n = 8
    L = (5.0, 5.0, 30.0)
    poses = [SVector(rand() * L[1], rand() * L[2], 10.0 + 10.0 * rand()) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]

    E_ewald  = coulomb_energy(Ewald2D(n, L; α = 1.7, s = 4.0), poses, charges)   # r_c = 2.35 < 2.5

    # The RAW 2D lattice sum is useless as a reference here: it converges as
    # 1/n_shell and is still 3.5% off at n_shell = 40. Use the extrapolated form,
    # which reaches ~3e-4 at (30, 60). Controller-measured, raw vs this Ewald value:
    #   n_shell = 10, 20, 30, 40, 60, 80  ->  13.5%, 6.9%, 4.6%, 3.5%, 2.3%, 1.7%
    E_direct = naive_energy_Q2D_extrap(poses, charges, L, 30, 60)
    @test isapprox(E_ewald, E_direct, rtol = 1e-3)

    # and confirm the extrapolation is doing the work, not a loose tolerance
    @test abs(E_ewald - E_direct) <
          abs(E_ewald - naive_energy_Q2D(poses, charges, L, 60))
end

@testset "Ewald2D force matches -grad(energy)" begin
    Random.seed!(20260924)
    n = 8
    L = (6.0, 6.0, 20.0)
    poses = [SVector(rand() * L[1], rand() * L[2], 5.0 + 10.0 * rand()) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]
    inter = Ewald2D(n, L; α = 1.3, s = 3.5)      # r_c = 2.69 < 3

    F = coulomb_force(inter, poses, charges)
    f = p -> coulomb_energy(inter, p, charges)
    for i in 1:n, d in 1:3
        @test isapprox(F[i][d], -fd_gradient(f, poses, i, d; h = 1e-5),
                       rtol = 1e-4, atol = 1e-8)
    end
end

@testset "Ewald2D net in-plane force vanishes" begin
    Random.seed!(20260925)
    n = 10
    L = (6.0, 6.0, 20.0)
    poses = [SVector(rand() * L[1], rand() * L[2], 5.0 + 10.0 * rand()) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]
    F = coulomb_force(Ewald2D(n, L; α = 1.3, s = 3.5), poses, charges)
    @test isapprox(sum(f -> f[1], F), 0.0, atol = 1e-9)
    @test isapprox(sum(f -> f[2], F), 0.0, atol = 1e-9)
    @test isapprox(sum(f -> f[3], F), 0.0, atol = 1e-9)
end

@testset "Ewald2DLong does not overflow for a tall box" begin
    # exp(k*z_ij) overflows for large k*z; the paired erfc underflows to zero at the
    # same time, so the product must be computed as zero rather than Inf*0 = NaN.
    L = (4.0, 4.0, 2000.0)
    poses = [SVector(1.0, 1.0, 10.0), SVector(2.0, 2.0, 1990.0)]
    charges = [1.0, -1.0]
    inter = Ewald2D(2, L; α = 1.6, s = 3.0)      # r_c = 1.875 < 2
    E = coulomb_energy(inter, poses, charges)
    @test isfinite(E)
    F = coulomb_force(inter, poses, charges)
    @test all(all(isfinite, f) for f in F)
end
