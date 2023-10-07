@testset "test the cell list 3D" begin
    n_atoms = 100
    n_atoms = Int64(round(n_atoms))
    L = 50.0
    boundary = CubicBoundary(L)
    atoms = Vector{Atom{Float64}}()

    for i in 1:n_atoms
        push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0))
    end

    info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.5, L - 0.5), boundary; min_r = 2.0, temp = 1.0)
    info.running_step = 1

    r_c = 5.0
    cell_list_3D = CellList3D(info, r_c, boundary, 1)
    cell_list_dir_3D = CellListDir3D(info, r_c, boundary, 1)

    @test cell_list_3D.neighbor_list == cell_list_dir_3D.neighbor_list

    for _ in 1:10
        info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.5, L - 0.5), boundary; min_r = 2.0, temp = 1.0)
        update_finder!(cell_list_3D, info)
        update_finder!(cell_list_dir_3D, info)

        @test cell_list_3D.neighbor_list == cell_list_dir_3D.neighbor_list
    end
end

@testset "test the cell list Q2D" begin
    n_atoms = 100
    n_atoms = Int64(round(n_atoms))
    L = 50.0
    boundary = Q2dBoundary(L, L, L)
    atoms = Vector{Atom{Float64}}()

    for i in 1:n_atoms
        push!(atoms, Atom(type = 1, mass = 1.0, charge = 1.0))
    end

    info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.5, L - 0.5), boundary; min_r = 2.0, temp = 1.0)
    info.running_step = 1

    r_c = 5.0
    cell_list_Q2D = CellListQ2D(info, r_c, boundary, 1)
    cell_list_dir_Q2D = CellListDirQ2D(info, r_c, boundary, 1)

    @test cell_list_Q2D.neighbor_list == cell_list_dir_Q2D.neighbor_list

    for _ in 1:10
        info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.5, L - 0.5), boundary; min_r = 2.0, temp = 1.0)
        update_finder!(cell_list_Q2D, info)
        update_finder!(cell_list_dir_Q2D, info)

        @test cell_list_Q2D.neighbor_list == cell_list_dir_Q2D.neighbor_list
    end
end