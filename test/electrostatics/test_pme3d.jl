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

@testset "PME3D and ICMPME3D are exported names" begin
    # This only checks that the stub declarations in long_pme3d.jl exist — it
    # passes whether or not FINUFFT (and therefore this extension) is loaded, and
    # would still pass with the entire fallback-error block deleted. It is not a
    # check that the "helpful error without FINUFFT" behaviour works; that is a
    # documented manual step — see Task 1 Step 6, which checks the fallback
    # message itself by other means (this file always has FINUFFT loaded).
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

@testset "PME3D composite matches Ewald3D" begin
    Random.seed!(1234)
    n = 30
    L = (12.0, 12.0, 12.0)
    α, s = 0.8, 3.5
    poses = [SVector(rand() * L[1], rand() * L[2], rand() * L[3]) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]

    e = Ewald3D(n, L; α = α, s = s)
    p = PME3D(n, L; α = α, s = s)
    @test p isa ExTinyMD.AbstractInteraction

    @test isapprox(coulomb_energy(p, poses, charges),
                   coulomb_energy(e, poses, charges), rtol = 1e-12)

    Fe = coulomb_force(e, poses, charges)
    Fp = coulomb_force(p, poses, charges)
    for i in 1:n, d in 1:3
        @test isapprox(Fp[i][d], Fe[i][d], rtol = 1e-10, atol = 1e-16)
    end
end

@testset "PME3D drives simulate! through the adapter" begin
    # The adapter and EwaldShort are reused unchanged; this confirms a PME3D
    # composite is a drop-in for Ewald3D in an actual MD run.
    Random.seed!(20260915)
    n, L = 30, 12.0
    boundary = Boundary((L, L, L), (1, 1, 1))
    atoms = Vector{Atom{Float64}}()
    for i in 1:(n ÷ 2); push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0)); end
    for i in (n ÷ 2 + 1):n; push!(atoms, Atom(type = 2, mass = 1.0, charge = -1.0)); end
    info = SimulationInfo(n, atoms, (0.0, L, 0.0, L, 0.0, L), boundary;
                          min_r = 1.0, temp = 1.0)
    info.running_step = 1

    inter = PME3D(n, (L, L, L); α = 0.8, s = 3.5)   # r_c = 4.375 < 6
    finder = CellList3D(info, inter.short.r_c, boundary, 1)
    sys = MDSys(n_atoms = n, atoms = atoms, boundary = boundary,
                interactions = [(inter, finder)],
                loggers = [TemperatureLogger(1000; output = false)],
                simulator = VerletProcess(dt = 1e-4))

    E0 = energy(inter, finder, sys, info)
    simulate!(sys.simulator, sys, info, 200)
    E1 = energy(inter, finder, sys, info)
    @test isfinite(E1)
    @test abs(E1 - E0) < 0.05 * max(abs(E0), 1.0)
end

@testset "PME3D force matches -grad(energy)" begin
    # Independent of the Ewald3D comparison: PME3D and Ewald3D share EwaldShort,
    # EwaldInteraction and the adapter, so agreement between them cannot detect a
    # defect in any of those. A finite-difference check constrains PME3D's own
    # energy against PME3D's own force and nothing else.
    Random.seed!(20260916)
    n = 12
    L = (12.0, 12.0, 12.0)
    poses = [SVector(rand() * L[1], rand() * L[2], rand() * L[3]) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]

    inter = PME3D(n, L; α = 0.8, s = 3.5)   # r_c = 4.375 < L/2 = 6

    F = coulomb_force(inter, poses, charges)
    f = p -> coulomb_energy(inter, p, charges)
    for i in 1:n, d in 1:3
        @test isapprox(F[i][d], -fd_gradient(f, poses, i, d; h = 1e-5),
                       rtol = 1e-4, atol = 1e-8)
    end
end

@testset "ICMPME3D matches ICMEwald3D" begin
    # Same physics, same k-set, different reciprocal-space engine.
    Random.seed!(20260928)
    n = 8
    L = (5.0, 5.0, 10.0)
    γ = (0.3, 0.3)
    poses = [SVector(rand() * L[1], rand() * L[2], 2.0 + 6.0 * rand()) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]

    # r_c = 2.35 < min(Lx,Ly)/2 = 2.5; N_image = 3 needs N_pad >= 2
    a = ICMEwald3D(n, L; α = 1.7, s = 4.0, γ = γ, N_image = 3, N_pad = 2)
    b = ICMPME3D(n, L;  α = 1.7, s = 4.0, γ = γ, N_image = 3, N_pad = 2)

    @test isapprox(coulomb_energy(b, poses, charges),
                   coulomb_energy(a, poses, charges), rtol = 1e-10)

    Fa = coulomb_force(a, poses, charges)
    Fb = coulomb_force(b, poses, charges)
    for i in 1:n, d in 1:3
        @test isapprox(Fb[i][d], Fa[i][d], rtol = 1e-8, atol = 1e-14)
    end
end

@testset "ICMPME3D force matches -grad(energy)" begin
    # Independent of the ICMEwald3D comparison above: a finite-difference check
    # constrains the PME force path on its own.
    Random.seed!(20260929)
    n = 6
    L = (6.0, 6.0, 5.0)
    γ = (0.8, 0.8)
    poses = [SVector(1.0, 1.0, 0.5), SVector(3.0, 1.5, 4.6), SVector(1.5, 3.5, 2.5),
             SVector(4.0, 4.0, 0.7), SVector(2.0, 4.5, 4.3), SVector(4.5, 2.0, 2.0)]
    charges = [1.0, -1.0, 1.0, -1.0, 1.0, -1.0]

    inter = ICMPME3D(n, L; α = 1.3, s = 3.5, γ = γ, N_image = 3, N_pad = 2)
    F = coulomb_force(inter, poses, charges)
    f = p -> coulomb_energy(inter, p, charges)
    for i in 1:n, d in 1:3
        @test isapprox(F[i][d], -fd_gradient(f, poses, i, d; h = 1e-5),
                       rtol = 1e-4, atol = 1e-8)
    end
end

@testset "PME3D and ICMPME3D forward ϵ" begin
    # Same check test_ewald.jl applies to Ewald3D. Without it, a dropped or
    # mis-forwarded ϵ would agree with a comparison partner that also defaulted
    # to ϵ = 1 and pass unnoticed.
    Random.seed!(20260917)
    n = 12
    L = (12.0, 12.0, 12.0)
    poses = [SVector(rand() * L[1], rand() * L[2], rand() * L[3]) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]

    E1 = coulomb_energy(PME3D(n, L; α = 0.8, s = 3.5, ϵ = 1.0), poses, charges)
    E2 = coulomb_energy(PME3D(n, L; α = 0.8, s = 3.5, ϵ = 2.0), poses, charges)
    @test isapprox(E2, E1 / 2, rtol = 1e-12)

    F1 = coulomb_force(PME3D(n, L; α = 0.8, s = 3.5, ϵ = 1.0), poses, charges)
    F2 = coulomb_force(PME3D(n, L; α = 0.8, s = 3.5, ϵ = 2.0), poses, charges)
    for i in 1:n
        @test isapprox(F2[i], F1[i] ./ 2, rtol = 1e-12)
    end

    # and the ICM route, on a confined slab
    Ls = (5.0, 5.0, 10.0)
    ps = [SVector(rand() * Ls[1], rand() * Ls[2], 2.0 + 6.0 * rand()) for _ in 1:n]
    kw = (; α = 1.7, s = 4.0, γ = (0.3, 0.3), N_image = 3, N_pad = 2)  # r_c = 2.353 < 2.5
    G1 = coulomb_energy(ICMPME3D(n, Ls; kw..., ϵ = 1.0), ps, charges)
    G2 = coulomb_energy(ICMPME3D(n, Ls; kw..., ϵ = 2.0), ps, charges)
    @test isapprox(G2, G1 / 2, rtol = 1e-12)
end

@testset "PME3DLong surface term responds to ϵ_inf" begin
    # FIX 4a coverage gap: every other testset in this file runs at the default
    # ϵ_inf = Inf, where 1/(2·Inf+1) makes the dipole surface term identically
    # zero in both long_energy and long_force! — deleting that code entirely
    # would not fail a single existing test. Mirrors
    # "Ewald3D surface term responds to ϵ_inf" in test_long_ewald3d.jl, but
    # checks PME3DLong (energy and force) against Ewald3DLong at finite ϵ_inf.
    poses = [SVector(1.0, 4.0, 4.0), SVector(7.0, 4.0, 4.0)]
    charges = [1.0, -1.0]
    L = (8.0, 8.0, 8.0)
    α, s = 0.8, 3.0   # r_c = 3.75 < 4

    for ϵ_inf in (1.0, 3.0)
        ewald = Ewald3DLong(2, L; α = α, s = s, ϵ_inf = ϵ_inf)
        pme   = PME3DLong(2, L;  α = α, s = s, ϵ_inf = ϵ_inf)
        @test isapprox(long_energy(pme, poses, charges),
                       long_energy(ewald, poses, charges), rtol = 1e-12)

        Fe = [zero(SVector{3,Float64}) for _ in 1:2]
        Fp = [zero(SVector{3,Float64}) for _ in 1:2]
        long_force!(Fe, ewald, poses, charges)
        long_force!(Fp, pme,   poses, charges)
        for i in 1:2, d in 1:3
            @test isapprox(Fp[i][d], Fe[i][d], rtol = 1e-10, atol = 1e-16)
        end
    end

    # and finite ϵ_inf must actually differ from the Inf (conducting) default —
    # otherwise the loop above could pass with the surface term dropped entirely
    # on both sides.
    cond = PME3DLong(2, L; α = α, s = s)   # ϵ_inf = Inf by default
    vac  = PME3DLong(2, L; α = α, s = s, ϵ_inf = 1.0)
    @test !isapprox(long_energy(cond, poses, charges), long_energy(vac, poses, charges))
end

@testset "PME3D preserves Float32 (and demonstrates FIX 3's tolerance fix)" begin
    # FIX 4b coverage gap: no Float32 test existed for the PME path. Mirrors
    # test_long_ewald3d.jl's "coulomb_energy/coulomb_force preserve Float32
    # through the composite" — `isfinite` is checked alongside the type, not
    # instead of it, since that earlier test's own comment records that a
    # type-only assertion once passed while the value was silently NaN.
    #
    # This also exercises FIX 3: before it, NUFFT_TOL = 1e-14 was below Float32
    # machine epsilon and FINUFFT warned on every call ("requested tolerance
    # epsilon too small to achieve", "increasing tol=1e-14 to eps_mach=1.19e-07").
    # `_nufft_tol(Float32)` now asks for 1f-6 instead, so this test should run
    # with no stderr warnings.
    L = (8.0f0, 8.0f0, 8.0f0)
    poses = [SVector(1.0f0, 2.0f0, 3.0f0), SVector(5.0f0, 6.0f0, 7.0f0)]
    charges = [1.0f0, -1.0f0]

    inter = PME3D(2, L; α = 0.8f0, s = 3.0f0)   # r_c = 3.75 < min(L)/2 = 4
    E = coulomb_energy(inter, poses, charges)
    @test E isa Float32
    @test isfinite(E)

    F = coulomb_force(inter, poses, charges)
    @test F isa Vector{SVector{3,Float32}}
    @test all(all(isfinite, f) for f in F)
end

@testset "PME3DLong long_energy and long_force! allocate nothing at steady state" begin
    # FIX 2 regression. `_structure_factor!` used to build a fresh
    # `Complex{T}[charges[j] for j in 1:m]` on every call — 16·n bytes, forever,
    # every timestep. Reviewer measured, at n = 500 with the default n_target
    # (the steady-state MD call): long_energy 8360 B, long_force! 8648 B, versus
    # 0 B for the equivalent Ewald3DLong calls. `PME3DLong` now carries a
    # preallocated `qs::Vector{Complex{T}}`, filled in place, and passes the
    # plan's own position/charge/output buffers straight through (rather than a
    # freshly built `view`) whenever the call is over the full particle set,
    # which is also what let this reach exactly 0 (see `_posview` in the
    # extension: a `SubArray` handed to `finufft_setpts!` gets boxed into the
    # plan's `AbstractVector{T}`-typed fields, at 288 B/call, unless the call is
    # skipped in favour of the concrete buffer itself).
    Random.seed!(555)
    n = 500
    L = (40.0, 40.0, 40.0)
    α, s = 0.4, 3.5      # r_c = 8.75 < 20
    poses = [SVector(rand() * L[1], rand() * L[2], rand() * L[3]) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]
    pme = PME3DLong(n, L; α = α, s = s)
    F = [zero(SVector{3,Float64}) for _ in 1:n]

    # Function barriers: an `@allocated` measurement taken directly at top level
    # would also see allocation from boxing non-const globals; wrapping the call
    # in a function with concrete argument types isolates the measurement to the
    # call itself. (No precedent for this exact idiom was found elsewhere in the
    # test suite at the time of writing; this is the standard `@allocated` +
    # function-barrier pattern.)
    _energy_call(long, poses, charges) = long_energy(long, poses, charges)
    _force_call!(F, long, poses, charges) = long_force!(F, long, poses, charges)

    _energy_call(pme, poses, charges)      # warm up: compile before measuring
    _force_call!(F, pme, poses, charges)

    @test @allocated(_energy_call(pme, poses, charges)) == 0
    @test @allocated(_force_call!(F, pme, poses, charges)) == 0
end

@testset "PME3DLong destroys its FINUFFT plans (does not leak them)" begin
    # FIX 1 regression (Major). `finufft_makeplan`'s C-side plan (FFTW plan,
    # sorted-point arrays, spreader workspace) gets no finalizer from FINUFFT.jl;
    # the only release path is `finufft_destroy!`, called here by a `finalizer`
    # attached to each plan in the constructor.
    #
    # An RSS-based version of this test was tried first (build ~25 PME3DLong
    # instances per batch at the ICMPME3D test's own grid — 56 sources,
    # 23×23×219 after masking — letting each go out of scope, `GC.gc(true)`,
    # compare batch-over-batch process RSS growth from /proc/self/statm). It
    # worked and clearly separated the two regimes when run standalone (with
    # the finalizer: -8, +14 MiB batch-over-batch across repeats; with it
    # removed: +37, +44 MiB) — but it flaked once (false failure) across
    # repeated runs of the full suite, where other tests' allocator activity
    # shares the same process. Per the brief, a leak test that fails
    # intermittently on CI is worse than none, so it was replaced with this
    # deterministic check instead: `Base.finalize` runs an object's registered
    # finalizers immediately, without needing it to be otherwise unreachable,
    # so it exercises exactly "is `finufft_destroy!` registered as a finalizer
    # on this plan" with no GC timing or OS memory-accounting involved. Verified
    # to fail (plan_ptr stays non-null) when the `finalizer(...)` calls in the
    # constructor are temporarily deleted.
    long = PME3DLong(10, (12.0, 12.0, 12.0); α = 0.8, s = 3.5)
    @test long.plan1.plan_ptr != C_NULL
    @test long.plan2.plan_ptr != C_NULL

    finalize(long.plan1)
    finalize(long.plan2)

    @test long.plan1.plan_ptr == C_NULL
    @test long.plan2.plan_ptr == C_NULL
end
