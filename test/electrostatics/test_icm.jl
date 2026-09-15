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
    # N_image grows. The GEOMETRY decides whether that is observable: with a thick
    # box and charges far from the walls the series is converged to machine
    # precision at N_image = 1, and the increments are then pure floating-point
    # noise with no ordering to assert. Controller-measured for a thick box
    # (L_z = 25, charges in z ∈ [8,17], γ = 0.4) the increments were
    # [2.2e-16, 0, 0, 3.3e-16, 2.2e-16] — an earlier draft of this test asserted
    # strict monotonic decrease on exactly that and failed on noise.
    #
    # So use a THIN slab with charges close to both walls and a strong dielectric
    # contrast, where the images genuinely contribute. Controller-measured
    # increments for the configuration below:
    #   [2.79e-5, 2.70e-7, 9.75e-10, 9.35e-12, 3.43e-14]
    Random.seed!(20260927)
    n = 6
    L = (5.0, 5.0, 4.0)
    γ = (0.9, 0.9)
    poses = [SVector(rand() * L[1], rand() * L[2], 0.8 + 2.4 * rand()) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]

    Es = [coulomb_energy(ICMEwald2D(n, L; α = 1.7, s = 4.0, γ = γ, N_image = m),
                         poses, charges) for m in 1:6]
    @test all(isfinite, Es)

    d = abs.(diff(Es))

    # the test must not be vacuous: the first image shell has to actually matter
    @test d[1] / abs(Es[end]) > 1e-5

    # increments shrink by at least 3x per shell, while they are above the noise
    # floor. Comparing noise against noise is what broke the earlier draft.
    noise = 1e-13 * abs(Es[end])
    for i in 1:(length(d) - 1)
        d[i] > noise || continue
        @test d[i + 1] < d[i] / 3
    end

    # and the series has converged by the last shell
    @test d[end] < 1e-9 * abs(Es[end])
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

@testset "ICM force matches -grad(energy) with real-image pairs in range" begin
    # The testset above uses a thick slab in which NO real particle sits within
    # r_c of a wall, so its reflected configuration produces zero real-image pairs
    # inside the cutoff (controller-verified: 3 real-real, 0 real-image). That
    # leaves the entire real-image branch of icm_short_energy / icm_short_force!
    # — the 1/2 weighting and the accumulate-on-the-real-index-only rule — unexercised.
    #
    # This configuration is deterministic rather than random, with two charges
    # placed inside r_c/2 of a wall so their own images are guaranteed inside the
    # cutoff on every run. Controller-measured: 2 real-real pairs, 4 real-image
    # pairs, worst finite-difference force error 6.2e-8.
    L = (6.0, 6.0, 5.0)
    γ = (0.8, 0.8)
    poses = [SVector(1.0, 1.0, 0.5), SVector(3.0, 1.5, 4.6), SVector(1.5, 3.5, 2.5),
             SVector(4.0, 4.0, 0.7), SVector(2.0, 4.5, 4.3), SVector(4.5, 2.0, 2.0)]
    charges = [1.0, -1.0, 1.0, -1.0, 1.0, -1.0]
    n = length(charges)

    inter = ICMEwald2D(n, L; α = 1.3, s = 3.5, γ = γ, N_image = 3)   # r_c = 2.692 < 3

    # Guard the intent: a charge at height z has its own image at -z, so the pair
    # separation is 2z. Assert that at least one real charge is close enough to a
    # wall for that pair to fall inside the cutoff, otherwise this test silently
    # degrades into a duplicate of the one above.
    r_c = inter.short.r_c
    @test 2 * minimum(p[3] for p in poses) < r_c

    F = coulomb_force(inter, poses, charges)
    f = p -> coulomb_energy(inter, p, charges)
    for i in 1:n, d in 1:3
        @test isapprox(F[i][d], -fd_gradient(f, poses, i, d; h = 1e-5),
                       rtol = 1e-5, atol = 1e-8)
    end
end

@testset "ICMEwald3D force matches -grad(energy) with real-image pairs in range" begin
    # FIX 3: both finite-difference testsets above use ICMEwald2D, so
    # `_elc_force!` and the `n_target < n_atoms` branch of `Ewald3DLong.long_force!`
    # (exercised whenever ICM's long-range solver is an `Ewald3DLong`, since real
    # particles are targets against all reflected charges as sources) are never
    # gradient-checked. Same thin-slab, real-image-pairs-in-range configuration
    # as the ICMEwald2D testset above, but through the ICMEwald3D + ELC route.
    #
    # Controller-measured: worst relative error 1.0e-5 on a force component of
    # magnitude 2.8e-5 (absolute error 3e-10), and 7.4e-8 elsewhere — hence the
    # generous atol, which matters only for near-zero components.
    L = (6.0, 6.0, 5.0)
    γ = (0.8, 0.8)
    poses = [SVector(1.0, 1.0, 0.5), SVector(3.0, 1.5, 4.6), SVector(1.5, 3.5, 2.5),
             SVector(4.0, 4.0, 0.7), SVector(2.0, 4.5, 4.3), SVector(4.5, 2.0, 2.0)]
    charges = [1.0, -1.0, 1.0, -1.0, 1.0, -1.0]
    n = length(charges)

    # r_c = 2.692 < min(Lx,Ly)/2 = 3. ICMEwald3D builds its Ewald3DLong for the
    # z-padded box (6, 6, 25) (N_pad = 2: (2*2+1)*5 = 25), whose smallest side is
    # still 6, so the same bound holds there too.
    inter = ICMEwald3D(n, L; α = 1.3, s = 3.5, γ = γ, N_image = 3, N_pad = 2)

    # Guard the intent: a charge at height z has its own image at -z, so the pair
    # separation is 2z; assert at least one real charge is close enough to a wall
    # for that pair to fall inside the cutoff.
    r_c = inter.short.r_c
    @test 2 * minimum(p[3] for p in poses) < r_c

    F = coulomb_force(inter, poses, charges)
    f = p -> coulomb_energy(inter, p, charges)
    for i in 1:n, d in 1:3
        @test isapprox(F[i][d], -fd_gradient(f, poses, i, d; h = 1e-5),
                       rtol = 1e-4, atol = 1e-8)
    end
end
