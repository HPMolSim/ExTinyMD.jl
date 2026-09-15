@testset "ewald_cutoffs" begin
    r_c, k_c = ExTinyMD.ewald_cutoffs(4.0, 0.2)
    @test r_c ≈ 20.0
    @test k_c ≈ 1.6
end

@testset "k_set_3D" begin
    L = (10.0, 10.0, 10.0)
    k_c = 2.0
    ks = ExTinyMD.k_set_3D(k_c, L)
    @test !isempty(ks)
    # cutoff respected, k=0 excluded, stored magnitude consistent
    for (kx, ky, kz, k) in ks
        @test 0 < k <= k_c + 1e-12
        @test k ≈ sqrt(kx^2 + ky^2 + kz^2)
    end
    # ±k symmetry: the set is closed under negation. `canon0` folds -0.0 to 0.0
    # so the Set membership check (which uses isequal, and isequal(-0.0, 0.0)
    # is false) isn't tripped up by the sign of an exact-zero component that
    # arithmetic negation of 0.0 produces.
    canon0(x) = iszero(x) ? zero(x) : x
    s = Set((kx, ky, kz) for (kx, ky, kz, _) in ks)
    @test all(((canon0(-kx), canon0(-ky), canon0(-kz)) in s) for (kx, ky, kz) in s)
    # a cubic box gives a k-set invariant under axis permutation
    s2 = Set((ky, kz, kx) for (kx, ky, kz) in s)
    @test s == s2
end

@testset "k_set_2D" begin
    L = (10.0, 20.0, 5.0)
    k_c = 1.5
    ks = ExTinyMD.k_set_2D(k_c, L)
    @test !isempty(ks)
    for (kx, ky, k) in ks
        @test 0 < k <= k_c + 1e-12
        @test k ≈ sqrt(kx^2 + ky^2)
    end
    # L_z must not influence the 2D k-set
    @test ks == ExTinyMD.k_set_2D(k_c, (10.0, 20.0, 999.0))
end

@testset "check_neutrality" begin
    @test ExTinyMD.check_neutrality([1.0, -1.0, 2.0, -2.0]) ≈ 0.0
    @test_logs (:warn,) ExTinyMD.check_neutrality([1.0, 1.0])
    @test ExTinyMD.check_neutrality([1.0, 1.0]) ≈ 2.0
end

@testset "min_image_disp" begin
    L = (10.0, 10.0, 10.0)
    # nearest image wraps
    d = ExTinyMD.min_image_disp(SVector(9.5, 0.0, 0.0), SVector(0.5, 0.0, 0.0),
                                L, ExTinyMD.Periodic3D())
    @test isapprox(d[1], -1.0)
    # Q2D does not wrap z
    d2 = ExTinyMD.min_image_disp(SVector(0.0, 0.0, 9.5), SVector(0.0, 0.0, 0.5),
                                 L, ExTinyMD.PeriodicQ2D())
    @test isapprox(d2[3], 9.0)
    d3 = ExTinyMD.min_image_disp(SVector(0.0, 0.0, 9.5), SVector(0.0, 0.0, 0.5),
                                 L, ExTinyMD.Periodic3D())
    @test isapprox(d3[3], -1.0)
    # works on Point and NTuple too (generic AoS requirement)
    @test ExTinyMD.min_image_disp(Point(9.5, 0.0, 0.0), Point(0.5, 0.0, 0.0),
                                  L, ExTinyMD.Periodic3D())[1] ≈ -1.0
    @test ExTinyMD.min_image_disp((9.5, 0.0, 0.0), (0.5, 0.0, 0.0),
                                  L, ExTinyMD.Periodic3D())[1] ≈ -1.0
end
