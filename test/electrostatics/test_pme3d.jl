using FINUFFT

@testset "PME3DLong energy equals Ewald3DLong exactly" begin
    # With the spherical D_k mask the two solvers sum over an identical k-set, so
    # this is an exact comparison rather than a convergence study. Controller
    # measured 6.4e-16 on this configuration.
    Random.seed!(4242)
    n = 40
    L = (12.0, 12.0, 12.0)
    α, s = 0.8, 4.0          # r_c = 5.0 < L/2 = 6
    poses = [SVector(rand() * L[1], rand() * L[2], rand() * L[3]) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]

    ewald = Ewald3DLong(n, L; α = α, s = s)
    pme   = PME3DLong(n, L; α = α, s = s)

    @test isapprox(long_energy(pme, poses, charges),
                   long_energy(ewald, poses, charges), rtol = 1e-12)
end

@testset "PME3DLong honours n_target" begin
    # ICMPME3D needs targets over real particles and sources over all reflected
    # charges. Controller measured 2e-16..8e-16 across these values.
    Random.seed!(99)
    n = 30
    L = (12.0, 12.0, 12.0)
    α, s = 0.8, 3.5
    poses = [SVector(rand() * L[1], rand() * L[2], rand() * L[3]) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]

    ewald = Ewald3DLong(n, L; α = α, s = s)
    pme   = PME3DLong(n, L; α = α, s = s)

    for nt in (n, 20, 10, 1)
        @test isapprox(long_energy(pme, poses, charges; n_target = nt),
                       long_energy(ewald, poses, charges; n_target = nt),
                       rtol = 1e-12)
    end
end

@testset "PME3DLong buffers survive interleaved calls" begin
    # Only masked-in grid points are written in some loops, and the plan reuses its
    # buffers across calls, so a stale buffer — or a skipped finufft_setpts! — would
    # contaminate a later result.
    #
    # The dirtying call must differ in point COUNT as well as in positions. An
    # earlier version of this test used the same n_target throughout, which meant a
    # reintroduced "skip setpts! when the count is unchanged" optimisation would
    # still pass: both B evaluations would land on the same original setpts state
    # and agree with each other while being wrong.
    Random.seed!(31337)
    n = 20
    L = (12.0, 12.0, 12.0)
    α, s = 0.8, 3.5
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]
    A = [SVector(rand() * L[1], rand() * L[2], rand() * L[3]) for _ in 1:n]
    B = [SVector(rand() * L[1], rand() * L[2], rand() * L[3]) for _ in 1:n]

    pme   = PME3DLong(n, L; α = α, s = s)
    ewald = Ewald3DLong(n, L; α = α, s = s)

    E_B_first = long_energy(pme, B, charges)
    long_energy(pme, A, charges; n_target = 7)   # different positions AND count
    E_B_after = long_energy(pme, B, charges)

    @test isapprox(E_B_first, E_B_after, rtol = 1e-14)
    # and it must still be the RIGHT answer, not merely self-consistent — two
    # equally-contaminated results would satisfy the assertion above
    @test isapprox(E_B_after, long_energy(ewald, B, charges), rtol = 1e-12)
    @test isfinite(E_B_after)
end

@testset "PME3D constructors error helpfully without FINUFFT" begin
    # Not testable in this file — FINUFFT is loaded. See Task 1 Step 6, which
    # checks the fallback message by other means.
    @test isdefined(ExTinyMD, :PME3D)
    @test isdefined(ExTinyMD, :ICMPME3D)
end
