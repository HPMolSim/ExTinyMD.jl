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
    # Verifies that a *reused* plan gives the correct answer after an interleaved
    # call at a different n_target, checked against an independent solver
    # (Ewald3DLong) rather than only self-consistency.
    #
    # The dirtying call must differ in point COUNT as well as in positions: if it
    # only changed positions, two equally-wrong results computed the same way could
    # still agree with each other while both being incorrect. Comparing against
    # Ewald3DLong at the end rules that out.
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

@testset "PME3DLong force equals Ewald3DLong force" begin
    # Controller measured 5.9e-15 relative / 2.6e-17 absolute on this shape.
    Random.seed!(7)
    n = 25
    L = (12.0, 12.0, 12.0)
    α, s = 0.8, 3.5
    poses = [SVector(rand() * L[1], rand() * L[2], rand() * L[3]) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]

    ewald = Ewald3DLong(n, L; α = α, s = s)
    pme   = PME3DLong(n, L; α = α, s = s)

    for nt in (n, 12)
        Fe = [zero(SVector{3,Float64}) for _ in 1:n]
        Fp = [zero(SVector{3,Float64}) for _ in 1:n]
        long_force!(Fe, ewald, poses, charges; n_target = nt)
        long_force!(Fp, pme,   poses, charges; n_target = nt)
        for i in 1:nt, d in 1:3
            @test isapprox(Fp[i][d], Fe[i][d], rtol = 1e-10, atol = 1e-16)
        end
    end
end

@testset "PME3DLong force accumulates, does not zero" begin
    # The composite coulomb_force! zeroes once and lets short and long accumulate.
    Random.seed!(8)
    n = 10
    L = (12.0, 12.0, 12.0)
    poses = [SVector(rand() * L[1], rand() * L[2], rand() * L[3]) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]
    pme = PME3DLong(n, L; α = 0.8, s = 3.5)

    F1 = [zero(SVector{3,Float64}) for _ in 1:n]
    long_force!(F1, pme, poses, charges)
    F2 = [zero(SVector{3,Float64}) for _ in 1:n]
    long_force!(F2, pme, poses, charges)
    long_force!(F2, pme, poses, charges)
    for i in 1:n
        @test isapprox(F2[i], 2 .* F1[i], rtol = 1e-12)
    end
end
