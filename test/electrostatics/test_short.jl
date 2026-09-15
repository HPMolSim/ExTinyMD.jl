@testset "EwaldShort energy against brute force" begin
    Random.seed!(20260914)
    n = 24
    L = (12.0, 12.0, 12.0)
    poses = [SVector(rand() * L[1], rand() * L[2], rand() * L[3]) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]

    α, s = 0.75, 4.0      # r_c = s/α = 5.33 < L/2 = 6
    short = EwaldShort(n, L; α = α, s = s)

    # brute force: same formula, every minimum-image pair inside r_c, plus self term
    function brute(poses, charges, L, α, r_c, ϵ)
        E = 0.0
        for i in 1:length(charges), j in (i+1):length(charges)
            d = ExTinyMD.min_image_disp(poses[i], poses[j], L, ExTinyMD.Periodic3D())
            r = sqrt(sum(abs2, d))
            r < r_c || continue
            E += charges[i] * charges[j] * erfc(α * r) / r
        end
        E -= α / sqrt(π) * sum(abs2, charges)
        return E / (4π * ϵ)
    end

    @test isapprox(short_energy(short, poses, charges),
                   brute(poses, charges, L, α, short.r_c, 1.0), rtol = 1e-12)
end

@testset "EwaldShort force matches -grad(energy)" begin
    Random.seed!(20260915)
    n = 12
    L = (12.0, 12.0, 12.0)
    poses = [SVector(rand() * L[1], rand() * L[2], rand() * L[3]) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]
    short = EwaldShort(n, L; α = 0.75, s = 3.5)   # r_c = 4.67 < 6

    F = [zero(SVector{3,Float64}) for _ in 1:n]
    short_force!(F, short, poses, charges)

    f = p -> short_energy(short, p, charges)
    for i in 1:n, d in 1:3
        @test isapprox(F[i][d], -fd_gradient(f, poses, i, d; h = 1e-5),
                       rtol = 1e-5, atol = 1e-9)
    end
end

@testset "EwaldShort Q2D convention does not wrap z" begin
    # two charges separated in z by more than L_z/2: Periodic3D sees the wrapped
    # image and Q2D sees the true separation, so the energies must differ
    L = (10.0, 10.0, 10.0)
    poses = [SVector(1.0, 1.0, 0.5), SVector(1.0, 1.0, 9.0)]
    charges = [1.0, -1.0]
    s3 = EwaldShort(2, L; α = 0.75, s = 3.0, convention = Periodic3D())
    sq = EwaldShort(2, L; α = 0.75, s = 3.0, convention = PeriodicQ2D())
    @test !isapprox(short_energy(s3, poses, charges), short_energy(sq, poses, charges))
end

@testset "EwaldShort accepts Point and NTuple positions" begin
    L = (10.0, 10.0, 10.0)
    sv = [SVector(1.0, 2.0, 3.0), SVector(2.0, 3.0, 4.0)]
    charges = [1.0, -1.0]
    short = EwaldShort(2, L; α = 0.75, s = 3.0)   # r_c = 4.0 < 5; pair at r = 1.73 counts
    E = short_energy(short, sv, charges)
    @test isapprox(short_energy(short, [Point(1.0, 2.0, 3.0), Point(2.0, 3.0, 4.0)],
                                charges), E)
    @test isapprox(short_energy(short, [(1.0, 2.0, 3.0), (2.0, 3.0, 4.0)], charges), E)
end
