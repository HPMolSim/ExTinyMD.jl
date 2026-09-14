@testset "icm_reflect! layout and recurrence" begin
    T = Float64
    n, N_image = 2, 2
    L = (4.0, 4.0, 10.0)
    γ = (0.5, 0.25)
    poses = [SVector(1.0, 1.0, 2.0), SVector(2.0, 3.0, 7.0)]
    charges = [1.0, -2.0]

    m = n * (1 + 2N_image)
    rp = [zero(SVector{3,T}) for _ in 1:m]
    rq = zeros(T, m)
    @test ExTinyMD.icm_reflect!(rp, rq, γ, L, N_image, poses, charges) == m

    # real particles come first, unchanged
    @test rp[1] == poses[1] && rq[1] == charges[1]
    @test rp[2] == poses[2] && rq[2] == charges[2]

    # particle 1's images occupy slots 3..6, interleaved up, down, up, down
    z1 = poses[1][3]
    @test isapprox(rp[3][3], 2 * L[3] - z1)          # up 1
    @test isapprox(rq[3],    γ[1] * charges[1])
    @test isapprox(rp[4][3], -z1)                    # down 1
    @test isapprox(rq[4],    γ[2] * charges[1])
    @test isapprox(rp[5][3], 2 * L[3] - (-z1))       # up 2 = 2Lz - down1
    @test isapprox(rq[5],    γ[1] * γ[2] * charges[1])
    @test isapprox(rp[6][3], -(2 * L[3] - z1))       # down 2 = -up1
    @test isapprox(rq[6],    γ[2] * γ[1] * charges[1])

    # images keep the parent's x,y
    for s in 3:6
        @test rp[s][1] == poses[1][1] && rp[s][2] == poses[1][2]
    end

    # γ = 0 kills every image charge
    rq0 = zeros(T, m)
    ExTinyMD.icm_reflect!(rp, rq0, (0.0, 0.0), L, N_image, poses, charges)
    @test all(iszero, rq0[(n+1):end])
end

@testset "ICM with γ = 0 reduces to the un-imaged method" begin
    Random.seed!(20260926)
    n = 10
    L = (6.0, 6.0, 20.0)
    poses = [SVector(rand() * L[1], rand() * L[2], 5.0 + 10.0 * rand()) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]

    plain = Ewald2D(n, L; α = 1.3, s = 3.5)      # r_c = 2.69 < 3
    icm   = ICMEwald2D(n, L; α = 1.3, s = 3.5, γ = (0.0, 0.0), N_image = 3)
    @test isapprox(coulomb_energy(icm, poses, charges),
                   coulomb_energy(plain, poses, charges), rtol = 1e-9)
end

@testset "ICM image series converges in N_image" begin
    # For |γ| < 1 the image series is geometric, so the energy must converge as
    # N_image grows, and successive increments must shrink. This is a real
    # property of the method and it does not presuppose any particular
    # convention for how real-image pairs are weighted.
    Random.seed!(20260927)
    n = 6
    L = (5.0, 5.0, 25.0)
    γ = (0.4, 0.4)
    poses = [SVector(rand() * L[1], rand() * L[2], 8.0 + 9.0 * rand()) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]

    Es = [coulomb_energy(ICMEwald2D(n, L; α = 1.7, s = 4.0, γ = γ, N_image = m),
                         poses, charges) for m in 1:6]
    d = abs.(diff(Es))
    @test all(isfinite, Es)
    # increments shrink monotonically, and the tail is small relative to the total
    @test all(d[i + 1] < d[i] for i in 1:(length(d) - 1))
    @test d[end] < 1e-3 * abs(Es[end])
end

@testset "ICM+Ewald3D+ELC agrees with ICM+Ewald2D" begin
    # Different algorithms, same physics. This is the strongest check available on
    # the ELC slab correction and the z-padding bookkeeping.
    Random.seed!(20260928)
    n = 8
    L = (5.0, 5.0, 10.0)
    γ = (0.3, 0.3)
    poses = [SVector(rand() * L[1], rand() * L[2], 2.0 + 6.0 * rand()) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]

    # r_c = 2.35 < min(Lx,Ly)/2 = 2.5. ICMEwald3D builds its Ewald3DLong for the
    # z-padded box (5, 5, 50), whose smallest side is still 5, so the same bound holds.
    #
    # N_pad must be large enough that periodic replicas of the whole *image stack*
    # do not interact — not merely the real slab. Controller-measured disagreement
    # between the two routes at N_image = 3:
    #   N_pad = 1  ->  9.4e-3     (padding too small; this is NOT an algorithm error)
    #   N_pad = 2  ->  4.7e-8
    #   N_pad = 3  ->  4.8e-8
    # and at N_image = 5, N_pad = 2 still gives 8.5e-4, needing N_pad >= 3.
    # With adequate padding the two algorithms agree to near machine precision, so
    # the tolerance here is 1e-6 rather than the 1e-3 an earlier draft used.
    E_2d = coulomb_energy(ICMEwald2D(n, L; α = 1.7, s = 4.0, γ = γ, N_image = 3),
                          poses, charges)
    E_3d = coulomb_energy(ICMEwald3D(n, L; α = 1.7, s = 4.0, γ = γ, N_image = 3,
                                     N_pad = 2), poses, charges)
    @test isapprox(E_2d, E_3d, rtol = 1e-6)
end

@testset "ICM force matches -grad(energy)" begin
    # Spec §7 flagged the ICM force convention as a risk: the energy halves
    # real-image pairs while the force uses the full field, which looks inconsistent.
    # It is not. An image moves at twice the rate of its source, and that factor of 2
    # cancels the 1/2 exactly. The controller verified this against finite differences
    # in a standalone implementation before this task was dispatched: worst relative
    # error 1.0e-7 over all components. So this test is expected to PASS.
    #
    # If it nevertheless fails, STOP and report it — do not change the formula to make
    # it pass. A failure now means the port diverges from the verified convention,
    # which is a bug in the port, not a licence to adjust the physics.
    #
    # Note on what is NOT tested here: an earlier draft of this plan hand-rolled a
    # direct-sum ICM oracle that applied its own real-image weighting. That tests the
    # plan author's guess at the convention rather than the ported code, so it was
    # dropped. ICM correctness rests on four checks that do not beg the question:
    # the γ = 0 reduction to plain Ewald2D, convergence of the image series in
    # N_image, agreement between ICM+Ewald2D and ICM+Ewald3D+ELC (two different
    # algorithms), and this finite-difference force check.
    Random.seed!(20260929)
    n = 6
    L = (6.0, 6.0, 20.0)
    γ = (0.3, 0.3)
    inter = ICMEwald2D(n, L; α = 1.3, s = 3.5, γ = γ, N_image = 3)   # r_c = 2.69 < 3
    poses = [SVector(rand() * L[1], rand() * L[2], 6.0 + 8.0 * rand()) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]

    F = coulomb_force(inter, poses, charges)
    f = p -> coulomb_energy(inter, p, charges)
    for i in 1:n, d in 1:3
        @test isapprox(F[i][d], -fd_gradient(f, poses, i, d; h = 1e-5),
                       rtol = 1e-5, atol = 1e-8)
    end
end
