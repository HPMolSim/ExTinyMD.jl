@testset "EwaldInteraction composes short and long" begin
    Random.seed!(20260918)
    n = 16
    L = (12.0, 12.0, 12.0)
    poses = [SVector(rand() * L[1], rand() * L[2], rand() * L[3]) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]

    α, s = 0.8, 4.0      # r_c = 5.0 < L/2 = 6
    inter = Ewald3D(n, L; α = α, s = s)
    @test inter isa ExTinyMD.AbstractInteraction

    short = EwaldShort(n, L; α = α, s = s)
    long  = Ewald3DLong(n, L; α = α, s = s)
    @test isapprox(coulomb_energy(inter, poses, charges),
                   short_energy(short, poses, charges) +
                   long_energy(long, poses, charges), rtol = 1e-12)
end

@testset "coulomb_force! zeroes its buffer" begin
    Random.seed!(20260919)
    n = 8
    L = (12.0, 12.0, 12.0)
    poses = [SVector(rand() * L[1], rand() * L[2], rand() * L[3]) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]
    inter = Ewald3D(n, L; α = 0.8, s = 4.0)

    F1 = coulomb_force(inter, poses, charges)
    # pre-filled buffer must not contaminate the result
    F2 = [SVector(99.0, 99.0, 99.0) for _ in 1:n]
    coulomb_force!(F2, inter, poses, charges)
    for i in 1:n
        @test isapprox(F1[i], F2[i], rtol = 1e-12)
    end
    # calling twice gives the same answer, i.e. no accumulation across calls
    coulomb_force!(F2, inter, poses, charges)
    for i in 1:n
        @test isapprox(F1[i], F2[i], rtol = 1e-12)
    end
end

@testset "Ewald3D total force matches -grad(energy)" begin
    Random.seed!(20260920)
    n = 10
    L = (12.0, 12.0, 12.0)
    poses = [SVector(rand() * L[1], rand() * L[2], rand() * L[3]) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]
    inter = Ewald3D(n, L; α = 0.8, s = 3.5)      # r_c = 4.375 < 6

    F = coulomb_force(inter, poses, charges)
    f = p -> coulomb_energy(inter, p, charges)
    for i in 1:n, d in 1:3
        @test isapprox(F[i][d], -fd_gradient(f, poses, i, d; h = 1e-5),
                       rtol = 1e-4, atol = 1e-8)
    end
end

@testset "Ewald3D net force vanishes" begin
    # Newton's third law: the total force on a periodic neutral system is zero
    Random.seed!(20260921)
    n = 14
    L = (12.0, 12.0, 12.0)
    poses = [SVector(rand() * L[1], rand() * L[2], rand() * L[3]) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]
    F = coulomb_force(Ewald3D(n, L; α = 0.8, s = 3.5), poses, charges)
    @test isapprox(sum(F), zero(SVector{3,Float64}), atol = 1e-10)
end
