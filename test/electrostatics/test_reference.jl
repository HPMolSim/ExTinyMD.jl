@testset "oracle: finite-difference helper" begin
    # gradient of a known scalar function of one particle's position
    poses = [SVector(0.3, 0.7, 1.1), SVector(2.0, 0.5, 0.25)]
    f = p -> 3 * p[1][1]^2 + 5 * p[1][2] * p[2][3] - p[1][3]^3
    @test isapprox(fd_gradient(f, poses, 1, 1), 6 * 0.3;          rtol = 1e-6)
    @test isapprox(fd_gradient(f, poses, 1, 2), 5 * 0.25;         rtol = 1e-6)
    @test isapprox(fd_gradient(f, poses, 1, 3), -3 * 1.1^2;       rtol = 1e-6)
    @test isapprox(fd_gradient(f, poses, 2, 3), 5 * 0.7;          rtol = 1e-6)
end

@testset "oracle: NaCl Madelung constant" begin
    # E_per_ion = -M q²/(4π ϵ0 r_nn), so M = -E * 8π * r_nn / N  with our 1/(4π)
    # convention: naive_energy_3D already applies the pair double-counting ½, so
    # recovering M from the total (already-halved) energy needs the 2 back.
    a = 2.0            # lattice constant; nearest-neighbour distance is a/2
    r_nn = a / 2
    poses, charges, L = nacl_lattice(1, a)
    @test length(charges) == 8
    @test sum(charges) == 0

    E = naive_energy_3D(poses, charges, L, 12)
    M = -E * 8π * r_nn / length(charges)
    # Cubic-shell truncation of the 8-ion cell (neutral, zero dipole, zero
    # quadrupole) converges fast, so the tolerance is tight.
    @test isapprox(M, 1.7475645946, atol = 1e-5)
end

@testset "oracle: Richardson extrapolation algebra" begin
    # A sequence that is exactly E_inf + c/n must be inverted exactly, from any
    # pair of shell counts. Accuracy on a real lattice sum is validated in Task 7,
    # where a converged Ewald2D reference exists; asserting it here would be
    # circular.
    E_inf, c = -0.25, 1.5
    fake(n) = E_inf + c / n
    ex(n1, n2) = (n2 * fake(n2) - n1 * fake(n1)) / (n2 - n1)
    @test isapprox(ex(10, 20), E_inf; rtol = 1e-12)
    @test isapprox(ex(30, 60), E_inf; rtol = 1e-12)
    @test isapprox(ex(40, 80), E_inf; rtol = 1e-12)
end

@testset "oracle: Q2D reduces to 3D for a tall box" begin
    # With one layer of charges and a box far taller than its width, the z-images
    # contribute negligibly, so the 3D and Q2D sums must agree.
    poses = [SVector(1.0, 1.0, 25.0), SVector(3.0, 2.0, 25.0),
             SVector(2.0, 3.5, 25.0), SVector(4.0, 4.5, 25.0)]
    charges = [1.0, -1.0, 1.0, -1.0]
    L = (6.0, 6.0, 400.0)
    @test isapprox(naive_energy_3D(poses, charges, L, 6),
                   naive_energy_Q2D(poses, charges, L, 6), rtol = 1e-3)
end
