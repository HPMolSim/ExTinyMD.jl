# ExTinyMD Electrostatics Stdlib — Phase 1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add Ewald3D, Ewald2D, ICM+Ewald2D and ICM+Ewald3D+ELC to ExTinyMD as a framework-free array API plus a thin MD adapter, with a self-contained test oracle and a Documenter site.

**Architecture:** Two layers. A numerical core taking a plan object plus AoS positions and charges — no ExTinyMD type is ever constructed to call it — and an adapter implementing `ExTinyMD.energy` / `ExTinyMD.update_acceleration!` on top. One shared real-space kernel (`EwaldShort`); interchangeable long-range solvers (`Ewald3DLong`, `Ewald2DLong`); ICM supplies its own split short-range kernel and decorates any long-range solver.

**Tech Stack:** Julia 1.10+, CellListMap 0.10, StaticArrays, SpecialFunctions (new dep), Documenter (docs env only).

**Spec:** `docs/superpowers/specs/2026-09-14-extinymd-electrostatics-design.md`

## Global Constraints

- `julia = "1.10"`. Do not use syntax newer than 1.10.
- `CellListMap = "0.10"`, `StaticArrays = "1"`. Add `SpecialFunctions = "2"` to `[deps]`.
- Do **not** add ForwardDiff, FINUFFT, EwaldSummations, or OhMyThreads to `[deps]`. FINUFFT is Phase 2 and enters as a `[weakdeps]` extension only.
- Do **not** add any dependency to the `test` target beyond `Test`. The oracle is written from scratch in the test suite (spec §8).
- Interchange type across the core API is `SVector{3,T}`. Core kernels must access positions only via `p[1]`, `p[2]`, `p[3]` so `Point{3,T}` and `NTuple{3,T}` also work.
- Never index a thread-local accumulator by `Threads.threadid()`. Partition work across tasks and give each task its own accumulator.
- Parameter convention, shared with every sibling package: `r_c = s/α`, `k_c = 2αs`.
- **`r_c` must be strictly less than half the smallest *periodic* box side.** CellListMap
  0.10 enforces this and raises
  `ArgumentError: UNIT CELL CHECK FAILED ... must be greater than 2*cutoff` otherwise —
  verified by probe, and the bound is strict, so `r_c == L/2` also fails. Since
  `r_c = s/α`, every choice of `α` and `s` is constrained by the box: `s/α < min(L)/2`.
  For `Periodic3D` all three sides count; for `PeriodicQ2D` only `L[1]` and `L[2]` do.
  Every parameter set in this plan has been checked against this. If you change one, re-check it.
- **CellListMap 0.10 in-place API uses keywords**: `InPlaceNeighborList(xpositions = ...)`
  and `update!(cl, xpositions = ...)`. The 0.9 spellings `InPlaceNeighborList(x = ...)` and
  positional `update!(cl, x)` raise `MethodError`. `neighborlist` likewise takes
  `xpositions =`. Follow `src/MD_core/neighbor_finder/cell_list.jl`, which is already on
  0.10; do **not** copy the idiom from `ParticleMeshEwald.jl`, which pins 0.9.
- Tests that build a `TemperatureLogger` must pass `output = false`. The default
  `output = true` opens and truncates `temperature.txt` in the working directory.
- Accumulate structure factors in `Complex{T}`, never a hardcoded `ComplexF64`.
- Prefactors are asymmetric and this is intentional: short-range carries `1/(4πϵ)`, Ewald2D's long-range carries `1/ϵ`. Ewald3D's long-range is written as `1/(2Vϵ)`, already folded. Do not "fix" this.
- Work on branch `electrostatics-stdlib`. Commit after every task.
- Run the full suite with `julia --project=. -e 'using Pkg; Pkg.test()'` from the package root.

---

## File Structure

| File | Responsibility |
|---|---|
| `src/interactions/electrostatics/common.jl` | k-set generation, neutrality check, minimum-image helpers, parameter conversion |
| `src/interactions/electrostatics/short.jl` | `EwaldShort` — shared real-space kernel, energy + force |
| `src/interactions/electrostatics/long_ewald3d.jl` | `Ewald3DLong` — direct k-sum + dipole term |
| `src/interactions/electrostatics/long_ewald2d.jl` | `Ewald2DLong` — 2D k-sum with z-kernel + k=0 term |
| `src/interactions/electrostatics/icm.jl` | `ICM_reflect!`, `ICMShort`, `ICM`, ELC slab term |
| `src/interactions/electrostatics/ewald.jl` | `EwaldInteraction` composite + the four public constructors |
| `src/interactions/electrostatics/adapter.jl` | `ExTinyMD.energy`, `ExTinyMD.update_acceleration!` |
| `test/electrostatics/reference.jl` | Self-contained oracle: naive lattice sum + finite-difference gradient |
| `test/electrostatics/*.jl` | One test file per source file above |
| `test/regression_finder.jl` | Regression test for the `neighbor_list` field bug |
| `docs/` | Documenter site |

Modified: `src/ExTinyMD.jl` (includes + exports), `src/types.jl` (field rename), `Project.toml` (SpecialFunctions), `test/runtests.jl`, `.github/workflows/CI.yml` (docs job).

---

### Task 1: Fix the `neighbor_list` field inconsistency

`AllNeighborFinder` and `NoNeighborFinder` name their field `neighborlist`; every interaction reads `.neighbor_list`. The shipped suite passes because it only exercises the `CellList*` finders, which use the correct name. This must be fixed first — Task 9's adapter tests use `NoNeighborFinder` for long-range-only interactions.

**Files:**
- Modify: `src/types.jl:86`, `src/types.jl:90`, `src/types.jl:97`
- Create: `test/regression_finder.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Consumes: nothing.
- Produces: `AllNeighborFinder{T}` and `NoNeighborFinder{T}` both expose field `neighbor_list::Vector{Tuple{Int64,Int64,T}}`.

- [ ] **Step 1: Write the failing test**

Create `test/regression_finder.jl`:

```julia
@testset "neighbor finders expose neighbor_list" begin
    # Regression: AllNeighborFinder/NoNeighborFinder named the field `neighborlist`
    # while every interaction reads `.neighbor_list`, so LJ + AllNeighborFinder threw.
    @test hasfield(ExTinyMD.AllNeighborFinder{Float64}, :neighbor_list)
    @test hasfield(ExTinyMD.NoNeighborFinder{Float64}, :neighbor_list)

    n_atoms = 20
    L = 20.0
    boundary = CubicBoundary(L)
    atoms = create_atoms([(n_atoms, Atom(type = 1, mass = 1.0))])
    info = SimulationInfo(n_atoms, atoms, (0.0, L, 0.0, L, 0.0, L), boundary;
                          min_r = 1.0, temp = 1.0)
    info.running_step = 1

    interaction = LennardJones(ϵ = 1.0, σ = 1.0, cutoff = 4.0)
    all_finder = AllNeighborFinder(n_atoms, Float64)

    # must not throw, and must produce a finite energy
    E = energy(interaction, all_finder, MDSys(
        n_atoms = n_atoms, atoms = atoms, boundary = boundary,
        interactions = [(interaction, all_finder)],
        loggers = [TemperatureLogger(100; output = false)],
        simulator = VerletProcess(dt = 0.001),
    ), info)
    @test isfinite(E)

    # NoNeighborFinder yields exactly zero pair energy for a cutoff-respecting interaction
    no_finder = NoNeighborFinder(Float64)
    @test isfinite(energy(interaction, no_finder, MDSys(
        n_atoms = n_atoms, atoms = atoms, boundary = boundary,
        interactions = [(interaction, no_finder)],
        loggers = [TemperatureLogger(100; output = false)],
        simulator = VerletProcess(dt = 0.001),
    ), info))
end
```

- [ ] **Step 2: Run it to verify it fails**

Add `include("regression_finder.jl")` to `test/runtests.jl` inside the top-level `@testset`, then run:

```
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: FAIL. The `hasfield` assertions return `false`, and the `energy` calls raise `FieldError: type AllNeighborFinder has no field neighbor_list`.

- [ ] **Step 3: Rename the field**

In `src/types.jl`, three edits:

```julia
struct AllNeighborFinder{T} <: AbstractNeighborFinder
    neighbor_list::Vector{Tuple{Int64, Int64, T}}
end
```

```julia
Base.show(io::IO, neighborfinder::AllNeighborFinder) = print(io, "AllNeighborFinder with $(length(neighborfinder.neighbor_list)) pairs")
```

```julia
struct NoNeighborFinder{T} <: AbstractNeighborFinder
    neighbor_list::Vector{Tuple{Int64, Int64, T}}
end
```

Then `grep -rn "neighborlist" src/` and confirm the only remaining hits are calls to CellListMap's `neighborlist` / `neighborlist!` functions in `cell_list.jl`, which are a different thing and must not be renamed.

- [ ] **Step 4: Run tests to verify they pass**

```
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS, and the pre-existing 25 tests still pass.

- [ ] **Step 5: Commit**

```bash
git add src/types.jl test/regression_finder.jl test/runtests.jl
git commit -m "fix: rename neighborlist field to neighbor_list on All/NoNeighborFinder

Every interaction reads .neighbor_list, so LennardJones + AllNeighborFinder
raised FieldError. The shipped suite missed it because it only exercised the
CellList finders, which already used the correct name."
```

---

### Task 2: Self-contained test oracle

Everything after this task is validated against this file, so it is built and verified first, against a literature constant rather than against our own code.

**Files:**
- Create: `test/electrostatics/reference.jl`
- Create: `test/electrostatics/test_reference.jl`
- Modify: `test/runtests.jl`

**Interfaces:**
- Produces:
  - `naive_energy_3D(poses, charges, L::NTuple{3,T}, n_shell::Int; ϵ=one(T))::T` — direct lattice sum over `±n_shell` images in all three axes, `1/(4πϵ)` prefactor, self-pair excluded within the home cell.
  - `naive_energy_Q2D(poses, charges, L::NTuple{3,T}, n_shell::Int; ϵ=one(T))::T` — same but images in x,y only.
  - `naive_energy_Q2D_extrap(poses, charges, L, n1::Int, n2::Int; ϵ=one(T))::T` — Richardson extrapolation of the quasi-2D sum, eliminating its `1/n_shell` tail.
  - `fd_gradient(f, poses, i::Int, d::Int; h)::T` — central difference of scalar `f(poses)` w.r.t. component `d` of particle `i`.
  - `nacl_lattice(n_cells::Int, a::T)` → `(poses::Vector{SVector{3,T}}, charges::Vector{T}, L::NTuple{3,T})`, a rock-salt configuration.

- [ ] **Step 1: Write the oracle and its self-test**

Create `test/electrostatics/reference.jl`:

```julia
using StaticArrays

"""
Direct lattice sum of the Coulomb energy, periodic in all three axes.
Deliberately naive and O(N² (2n_shell+1)³) — this is the oracle, so it is written
for obviousness, not speed. Converges slowly and only conditionally; use a neutral
configuration and a generous `n_shell`.
"""
function naive_energy_3D(poses, charges, L::NTuple{3,T}, n_shell::Int; ϵ::T = one(T)) where T
    E = zero(T)
    N = length(charges)
    for i in 1:N, j in 1:N
        qq = charges[i] * charges[j]
        for mx in -n_shell:n_shell, my in -n_shell:n_shell, mz in -n_shell:n_shell
            # skip only the i==j self term in the home cell; images of self do count
            (i == j && mx == 0 && my == 0 && mz == 0) && continue
            dx = poses[i][1] - poses[j][1] - mx * L[1]
            dy = poses[i][2] - poses[j][2] - my * L[2]
            dz = poses[i][3] - poses[j][3] - mz * L[3]
            E += qq / sqrt(dx^2 + dy^2 + dz^2)
        end
    end
    return E / (2 * 4π * ϵ)
end

"Direct lattice sum, periodic in x and y only (quasi-2D slab)."
function naive_energy_Q2D(poses, charges, L::NTuple{3,T}, n_shell::Int; ϵ::T = one(T)) where T
    E = zero(T)
    N = length(charges)
    for i in 1:N, j in 1:N
        qq = charges[i] * charges[j]
        for mx in -n_shell:n_shell, my in -n_shell:n_shell
            (i == j && mx == 0 && my == 0) && continue
            dx = poses[i][1] - poses[j][1] - mx * L[1]
            dy = poses[i][2] - poses[j][2] - my * L[2]
            dz = poses[i][3] - poses[j][3]
            E += qq / sqrt(dx^2 + dy^2 + dz^2)
        end
    end
    return E / (2 * 4π * ϵ)
end

"""
Richardson-extrapolated quasi-2D lattice sum.

The truncated 2D sum converges as `1/n_shell` — measured relative error 13.5%, 6.9%,
4.6%, 3.5% at `n_shell` = 10, 20, 30, 40 — which is far too slow to compare against an
Ewald result directly. Eliminating the `1/n` term with two shell counts reaches about
3e-4 at `(30, 60)`:

    E_inf ≈ (n2*E(n2) − n1*E(n1)) / (n2 − n1)

Use this, not the raw sum, whenever comparing against a converged method.
"""
function naive_energy_Q2D_extrap(poses, charges, L::NTuple{3,T}, n1::Int, n2::Int;
                                 ϵ::T = one(T)) where T
    E1 = naive_energy_Q2D(poses, charges, L, n1; ϵ = ϵ)
    E2 = naive_energy_Q2D(poses, charges, L, n2; ϵ = ϵ)
    return (n2 * E2 - n1 * E1) / (n2 - n1)
end

"Central finite difference of `f(poses)` w.r.t. component `d` of particle `i`."
function fd_gradient(f, poses::Vector{SVector{3,T}}, i::Int, d::Int;
                     h::T = cbrt(eps(T))) where T
    shift = SVector{3,T}(ntuple(k -> k == d ? h : zero(T), 3))
    p = copy(poses)
    p[i] = poses[i] + shift
    fp = f(p)
    p[i] = poses[i] - shift
    fm = f(p)
    return (fp - fm) / (2h)
end

"""
Rock-salt (NaCl) lattice: `n_cells`³ conventional cells of edge `a`, two
interpenetrating FCC sublattices of opposite charge. Returns AoS positions,
charges, and the periodic box.
"""
function nacl_lattice(n_cells::Int, a::T) where T
    poses = SVector{3,T}[]
    charges = T[]
    h = a / 2
    for ix in 0:(2n_cells - 1), iy in 0:(2n_cells - 1), iz in 0:(2n_cells - 1)
        push!(poses, SVector{3,T}(ix * h, iy * h, iz * h))
        push!(charges, iseven(ix + iy + iz) ? one(T) : -one(T))
    end
    L = (T(n_cells) * a, T(n_cells) * a, T(n_cells) * a)
    return poses, charges, L
end
```

Create `test/electrostatics/test_reference.jl`:

```julia
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
    # E_per_ion = -M q²/(4π ϵ0 r_nn), so M = -E * 4π * r_nn / N  with our 1/(4π) convention
    a = 2.0            # lattice constant; nearest-neighbour distance is a/2
    r_nn = a / 2
    poses, charges, L = nacl_lattice(1, a)
    @test length(charges) == 8
    @test sum(charges) == 0

    E = naive_energy_3D(poses, charges, L, 12)
    # E_total = -N*M / (2 * 4π * r_nn) in this unit convention. The 2 is the pair
    # double-counting factor naive_energy_3D already applies, so recovering M needs
    # it back — hence 8π, not 4π. Measured sequence (controller-verified against an
    # independent Ewald3D implementation):
    #   n_shell =  4  ->  M = 1.7475584843
    #   n_shell =  8  ->  M = 1.7475641146
    #   n_shell = 12  ->  M = 1.7475644920
    #   n_shell = 16  ->  M = 1.7475645609
    #   n_shell = 20  ->  M = 1.7475645804
    M = -E * 8π * r_nn / length(charges)
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
```

- [ ] **Step 2: Run to verify**

Add to `test/runtests.jl`. Note that `reference.jl` is included **once here**, ahead of
the testsets, not from inside `test_reference.jl` — Tasks 5, 7 and 8 all consume
`nacl_lattice`, `fd_gradient` and `naive_energy_Q2D`, so the oracle must be unambiguously
in scope for every electrostatics test file:

```julia
include("electrostatics/reference.jl")

@testset "electrostatics" begin
    include("electrostatics/test_reference.jl")
end
```

Run `julia --project=. -e 'using Pkg; Pkg.test()'`.

Expected: PASS.

If the Madelung test fails, the oracle is wrong and **must** be fixed before proceeding.
Check the sublattice charge assignment and the `i == j` home-cell exclusion first. Do not
loosen the tolerance on your own judgement: instead print the recovered constant at several
shell counts and put the sequence in your report, so the controller can adjudicate against
real numbers.

```julia
for n_shell in (8, 12, 16, 20)
    E = naive_energy_3D(poses, charges, L, n_shell)
    println("n_shell=", n_shell, "  M=", -E * 8π * (a/2) / length(charges))
end
```

The 8-ion cube is neutral with zero dipole and zero quadrupole, so cubic-shell truncation
converges quickly here — about 1e-7 by n_shell = 12. Report the sequence either way. If your
numbers differ materially from the ones quoted in the test above, say so rather than
adjusting anything: that would mean a real bug on one side or the other.

- [ ] **Step 3: Commit**

```bash
git add test/electrostatics/ test/runtests.jl
git commit -m "test: add self-contained electrostatics oracle

Naive direct lattice sums (3D and quasi-2D) plus a finite-difference gradient
helper, written from scratch so the suite takes no dependency on the packages
this stdlib supersedes. Validated against the NaCl Madelung constant."
```

---

### Task 3: `common.jl` — k-sets, neutrality, minimum image

**Files:**
- Create: `src/interactions/electrostatics/common.jl`
- Create: `test/electrostatics/test_common.jl`
- Modify: `src/ExTinyMD.jl`, `Project.toml`, `test/runtests.jl`

**Interfaces:**
- Consumes: nothing.
- Produces:
  - `k_set_3D(k_c::T, L::NTuple{3,T})::Vector{NTuple{4,T}}` — entries `(kx, ky, kz, k)`, `0 < k ≤ k_c`.
  - `k_set_2D(k_c::T, L::NTuple{3,T})::Vector{NTuple{3,T}}` — entries `(kx, ky, k)`.
  - `check_neutrality(charges; atol)::T` — returns `Σq`, warns once if `|Σq| > atol`.
  - `ewald_cutoffs(s::T, α::T)::NTuple{2,T}` — returns `(r_c, k_c) = (s/α, 2αs)`.
  - `Periodic3D`, `PeriodicQ2D` — singleton boundary-convention tags.
  - `min_image_disp(pi, pj, L::NTuple{3,T}, conv)::SVector{3,T}`.

- [ ] **Step 1: Write the failing tests**

Create `test/electrostatics/test_common.jl`:

```julia
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
    # ±k symmetry: the set is closed under negation.
    #
    # Fold -0.0 to 0.0 before any set membership. `Set` compares with `isequal`,
    # and `isequal(-0.0, 0.0)` is false even though `-0.0 == 0.0` is true, so
    # negating a k-vector with a zero component produces a key that is absent from
    # the set for reasons of floating-point sign, not of symmetry. Without the fold
    # this test fails deterministically against a perfectly symmetric k-set.
    fold(x::T) where {T} = x == 0 ? zero(T) : x
    key(kx, ky, kz) = (fold(kx), fold(ky), fold(kz))

    s = Set(key(kx, ky, kz) for (kx, ky, kz, _) in ks)
    # the fold must not merge distinct k-vectors, or the test above it goes vacuous
    @test length(s) == length(ks)
    @test all((key(-kx, -ky, -kz) in s) for (kx, ky, kz, _) in ks)
    # a cubic box gives a k-set invariant under axis permutation
    s2 = Set(key(ky, kz, kx) for (kx, ky, kz) in s)
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
```

- [ ] **Step 2: Run to verify it fails**

Add `include("electrostatics/test_common.jl")` to the `"electrostatics"` testset in `test/runtests.jl`. Run the suite.

Expected: FAIL with `UndefVarError: ewald_cutoffs not defined in ExTinyMD`.

- [ ] **Step 3: Add the SpecialFunctions dependency**

In `Project.toml`, add to `[deps]`:

```toml
SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"
```

and to `[compat]`:

```toml
SpecialFunctions = "2"
```

Then `julia --project=. -e 'using Pkg; Pkg.resolve()'`.

- [ ] **Step 4: Write the implementation**

Create `src/interactions/electrostatics/common.jl`:

```julia
# Boundary conventions. These select which axes wrap under the minimum-image
# convention, and nothing else.
abstract type AbstractBoundaryConvention end
struct Periodic3D  <: AbstractBoundaryConvention end
struct PeriodicQ2D <: AbstractBoundaryConvention end

"""
    ewald_cutoffs(s, α) -> (r_c, k_c)

Real- and reciprocal-space cutoffs from the splitting parameter `α` and the
dimensionless accuracy parameter `s`: `r_c = s/α`, `k_c = 2αs`. Larger `s` means
more accuracy and more work in both spaces. This convention is shared with
EwaldSummations, QuasiEwald, SoEwald2D and ParticleMeshEwald, so parameters
transfer between them unchanged.
"""
ewald_cutoffs(s::T, α::T) where {T} = (s / α, 2 * α * s)

"""
    k_set_3D(k_c, L) -> Vector{NTuple{4,T}}

Reciprocal lattice vectors `(k_x, k_y, k_z, |k|)` of the box `L` with
`0 < |k| ≤ k_c`. The set is closed under `k -> -k`; the long-range sums assume
that and therefore carry no factor-of-two correction.
"""
function k_set_3D(k_c::T, L::NTuple{3,T}) where {T}
    mx_max = ceil(Int, k_c * L[1] / 2π) + 1
    my_max = ceil(Int, k_c * L[2] / 2π) + 1
    mz_max = ceil(Int, k_c * L[3] / 2π) + 1

    k_set = Vector{NTuple{4,T}}()
    for m_x in -mx_max:mx_max, m_y in -my_max:my_max, m_z in -mz_max:mz_max
        k_x = m_x * 2π / L[1]
        k_y = m_y * 2π / L[2]
        k_z = m_z * 2π / L[3]
        k = sqrt(k_x^2 + k_y^2 + k_z^2)
        if 0 < k <= k_c
            push!(k_set, (k_x, k_y, k_z, k))
        end
    end
    return k_set
end

"""
    k_set_2D(k_c, L) -> Vector{NTuple{3,T}}

In-plane reciprocal lattice vectors `(k_x, k_y, |k|)` with `0 < |k| ≤ k_c`.
`L[3]` is ignored: the z axis is not periodic in the quasi-2D geometry.
"""
function k_set_2D(k_c::T, L::NTuple{3,T}) where {T}
    mx_max = ceil(Int, k_c * L[1] / 2π) + 1
    my_max = ceil(Int, k_c * L[2] / 2π) + 1

    k_set = Vector{NTuple{3,T}}()
    for m_x in -mx_max:mx_max, m_y in -my_max:my_max
        k_x = m_x * 2π / L[1]
        k_y = m_y * 2π / L[2]
        k = sqrt(k_x^2 + k_y^2)
        if 0 < k <= k_c
            push!(k_set, (k_x, k_y, k))
        end
    end
    return k_set
end

"""
    check_neutrality(charges; atol) -> Σq

Return the net charge, warning once when it is non-zero. Ewald summation of a
non-neutral system is conditionally convergent and picks up a box-volume-dependent
offset, which shows up as an energy that drifts with `L` rather than as an error.
"""
function check_neutrality(charges::AbstractVector{T}; atol::T = sqrt(eps(T))) where {T}
    net = sum(charges)
    if abs(net) > atol
        @warn "System is not charge neutral (Σq = $net); Ewald energies carry a " *
              "volume-dependent offset." maxlog = 1
    end
    return net
end

# Nearest-image displacement. `dx - L*round(dx/L)` is the true minimum image, unlike
# ExTinyMD's position_check3D, which returns the first image inside the cutoff and is
# equivalent only while r_c < L/2.
@inline _wrap(dx::T, L::T) where {T} = dx - L * round(dx / L)

@inline function min_image_disp(p_i, p_j, L::NTuple{3,T}, ::Periodic3D) where {T}
    return SVector{3,T}(_wrap(T(p_i[1]) - T(p_j[1]), L[1]),
                        _wrap(T(p_i[2]) - T(p_j[2]), L[2]),
                        _wrap(T(p_i[3]) - T(p_j[3]), L[3]))
end

@inline function min_image_disp(p_i, p_j, L::NTuple{3,T}, ::PeriodicQ2D) where {T}
    return SVector{3,T}(_wrap(T(p_i[1]) - T(p_j[1]), L[1]),
                        _wrap(T(p_i[2]) - T(p_j[2]), L[2]),
                        T(p_i[3]) - T(p_j[3]))
end
```

In `src/ExTinyMD.jl`, add `SpecialFunctions` to the `using` line and add the include after the existing interaction includes:

```julia
using LinearAlgebra, Random, Distributions, CellListMap, StaticArrays, DelimitedFiles, SpecialFunctions
```

```julia
# electrostatics standard library
include("interactions/electrostatics/common.jl")
```

Export the public names (not the internal helpers):

```julia
export Periodic3D, PeriodicQ2D, ewald_cutoffs, check_neutrality
```

- [ ] **Step 5: Run tests to verify they pass**

```
julia --project=. -e 'using Pkg; Pkg.test()'
```

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add Project.toml src/ExTinyMD.jl src/interactions/electrostatics/common.jl test/
git commit -m "feat: add electrostatics common utilities

k-set generation for 3D and 2D periodicity, charge-neutrality check, Ewald
cutoff convention, and true minimum-image displacement."
```

---

### Task 4: `EwaldShort` — the shared real-space kernel

**Files:**
- Create: `src/interactions/electrostatics/short.jl`
- Create: `test/electrostatics/test_short.jl`
- Modify: `src/ExTinyMD.jl`, `test/runtests.jl`

**Interfaces:**
- Consumes: `Periodic3D`, `PeriodicQ2D`, `min_image_disp`, `ewald_cutoffs` from Task 3.
- Produces:
  - `EwaldShort(n_atoms::Int, L::NTuple{3,T}; α, s, ϵ=one(T), convention=Periodic3D())`
  - `short_energy(short::EwaldShort{T}, poses, charges; neighbor_list=nothing)::T`
  - `short_force!(F::Vector{SVector{3,T}}, short, poses, charges; neighbor_list=nothing)` — **accumulates into `F`, does not zero it**
  - fields `α`, `r_c`, `k_c`, `ϵ`, `L`, `n_atoms`, `convention`

- [ ] **Step 1: Write the failing tests**

Create `test/electrostatics/test_short.jl`:

```julia
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
```

Add `using Random` and `using SpecialFunctions` at the top of `test/runtests.jl` (the test suite may use `erfc` directly; `SpecialFunctions` is a dependency of the package so it is available in the test env without adding it to the test target).

- [ ] **Step 2: Run to verify it fails**

Add the include; run the suite. Expected: FAIL with `UndefVarError: EwaldShort not defined`.

- [ ] **Step 3: Write the implementation**

Create `src/interactions/electrostatics/short.jl`:

```julia
"""
    EwaldShort(n_atoms, L; α, s, ϵ = 1.0, convention = Periodic3D())

Real-space part of an Ewald split, shared by every method in this library:

    E_s = 1/(4πϵ) [ Σ_{i<j, r<r_c} q_i q_j erfc(α r_ij)/r_ij − (α/√π) Σ_i q_i² ]

`convention` selects which axes wrap — [`Periodic3D`](@ref) for a triply periodic
box, [`PeriodicQ2D`](@ref) for a slab periodic in x and y only.

The struct owns a `CellListMap` neighbour list so that [`short_energy`](@ref) and
[`short_force!`](@ref) allocate nothing after construction. Pass
`neighbor_list = ...` to reuse a list maintained elsewhere, as the MD adapter does.
"""
mutable struct EwaldShort{T, C <: AbstractBoundaryConvention, TC}
    α::T
    r_c::T
    k_c::T
    ϵ::T
    L::NTuple{3,T}
    n_atoms::Int
    convention::C
    cell_list::TC
    pos_buffer::Vector{SVector{3,T}}
end

function EwaldShort(n_atoms::Int, L::NTuple{3,T}; α::T, s::T, ϵ::T = one(T),
                    convention::C = Periodic3D()) where {T, C <: AbstractBoundaryConvention}
    r_c, k_c = ewald_cutoffs(s, α)
    pos_buffer = [zero(SVector{3,T}) for _ in 1:n_atoms]
    cell_list = InPlaceNeighborList(xpositions = pos_buffer, cutoff = r_c,
                                    unitcell = _cell_unitcell(L, r_c, convention),
                                    parallel = true)
    return EwaldShort{T, C, typeof(cell_list)}(α, r_c, k_c, ϵ, L, n_atoms, convention,
                                               cell_list, pos_buffer)
end

Base.show(io::IO, s::EwaldShort) =
    print(io, "EwaldShort(α = $(s.α), r_c = $(s.r_c), ϵ = $(s.ϵ), $(s.convention))")

# For a non-periodic axis, inflate the unitcell so CellListMap finds no images along it.
_cell_unitcell(L::NTuple{3,T}, r_c::T, ::Periodic3D) where {T} =
    SVector{3,T}(L[1], L[2], L[3])
_cell_unitcell(L::NTuple{3,T}, r_c::T, ::PeriodicQ2D) where {T} =
    SVector{3,T}(L[1], L[2], max(L[3] + 2 * r_c, T(2) * r_c))

function _refresh_neighbors!(short::EwaldShort{T}, poses) where {T}
    @inbounds for i in 1:short.n_atoms
        p = poses[i]
        short.pos_buffer[i] = SVector{3,T}(T(p[1]), T(p[2]), T(p[3]))
    end
    update!(short.cell_list, xpositions = short.pos_buffer)
    return neighborlist!(short.cell_list)
end

"""
    short_energy(short, poses, charges; neighbor_list = nothing) -> T

Real-space Ewald energy. `poses` is AoS — any vector whose elements support
`p[1]`, `p[2]`, `p[3]`.
"""
function short_energy(short::EwaldShort{T}, poses, charges;
                      neighbor_list = nothing) where {T}
    nb = neighbor_list === nothing ? _refresh_neighbors!(short, poses) : neighbor_list
    α, r_c = short.α, short.r_c

    E = zero(T)
    @inbounds for (i, j, r) in nb
        (r < r_c && r > zero(T)) || continue
        E += charges[i] * charges[j] * erfc(α * r) / r
    end

    @inbounds for i in 1:short.n_atoms
        E -= charges[i]^2 * α / sqrt(T(π))
    end

    return E / (4π * short.ϵ)
end

# -dE/dr for E(r) = q_i q_j erfc(α r)/r
@inline function _short_pair_dEdr(q_i::T, q_j::T, α::T, r::T) where {T}
    return q_i * q_j * (erfc(α * r) / r^2 + 2α / sqrt(T(π)) * exp(-(α * r)^2) / r)
end

"""
    short_force!(F, short, poses, charges; neighbor_list = nothing)

Accumulate the real-space force into `F`. **Does not zero `F` first** — the
composite interaction zeroes once and lets short and long parts accumulate.
"""
function short_force!(F::Vector{SVector{3,T}}, short::EwaldShort{T}, poses, charges;
                      neighbor_list = nothing) where {T}
    nb = neighbor_list === nothing ? _refresh_neighbors!(short, poses) : neighbor_list
    α, r_c, ϵ = short.α, short.r_c, short.ϵ
    conv, L = short.convention, short.L
    pref = one(T) / (4π * ϵ)

    @inbounds for (i, j, r) in nb
        (r < r_c && r > zero(T)) || continue
        d = min_image_disp(poses[i], poses[j], L, conv)
        F_ij = _short_pair_dEdr(charges[i], charges[j], α, r) * d / r * pref
        F[i] += F_ij
        F[j] -= F_ij
    end
    return F
end
```

Note: the self-energy term is position-independent, so it contributes nothing to the force and is absent from `short_force!` by design.

Add to `src/ExTinyMD.jl`:

```julia
include("interactions/electrostatics/short.jl")
```
```julia
export EwaldShort, short_energy, short_force!
```

- [ ] **Step 4: Run tests to verify they pass**

Expected: PASS, including the finite-difference force check.

If the FD check fails only at a handful of `(i, d)` with small absolute force, the `atol` is doing its job and the failure is real — check the `F[j] -= F_ij` sign and that `d` runs from `i` to `j` (not `j` to `i`), since that sign error passes the energy test and fails only here.

- [ ] **Step 5: Commit**

```bash
git add src/ExTinyMD.jl src/interactions/electrostatics/short.jl test/
git commit -m "feat: add shared real-space Ewald kernel

One erfc pair kernel for every method, with an analytic pair derivative rather
than per-pair ForwardDiff. Validated against brute force and against the
finite-difference energy gradient."
```

---

### Task 5: `Ewald3DLong` — direct reciprocal-space sum

**Files:**
- Create: `src/interactions/electrostatics/long_ewald3d.jl`
- Create: `test/electrostatics/test_long_ewald3d.jl`
- Modify: `src/ExTinyMD.jl`, `test/runtests.jl`

**Interfaces:**
- Consumes: `k_set_3D`, `ewald_cutoffs` from Task 3.
- Produces:
  - `Ewald3DLong(n_atoms::Int, L::NTuple{3,T}; α, s, ϵ=one(T), ϵ_inf=T(Inf))`
  - `long_energy(long::Ewald3DLong{T}, poses, charges; n_target::Int = long.n_atoms)::T`
  - `long_force!(F::Vector{SVector{3,T}}, long::Ewald3DLong{T}, poses, charges; n_target::Int = long.n_atoms)` — accumulates

`n_target` exists for ICM (Task 8), which needs the energy summed over **real particles as
targets** against **all reflected charges as sources**. When `n_target == long.n_atoms`
(the default, and every non-ICM use) the target and source sets coincide and the
expressions reduce to the ordinary Ewald ones. Do not omit it: Task 8 depends on it, and
summing targets over the image charges too would add unphysical image-image
self-interaction. Controller-verified: real-only targets make ICM+Ewald2D and
ICM+Ewald3D+ELC agree to 4.7e-8, while all-reflected targets leave them 5.8% apart.
  - fields `α`, `k_c`, `r_c`, `ϵ`, `ϵ_inf`, `L`, `n_atoms`, `k_set`

- [ ] **Step 1: Write the failing tests**

Create `test/electrostatics/test_long_ewald3d.jl`:

```julia
@testset "Ewald3D total energy against direct lattice sum" begin
    # NaCl: neutral, symmetric, and with a known answer independent of our code
    a = 2.0
    poses, charges, L = nacl_lattice(2, a)   # 64 ions
    n = length(charges)

    α, s = 2.1, 4.0      # L = (4,4,4) so r_c must be < 2; s/α = 1.90
    short = EwaldShort(n, L; α = α, s = s)
    long  = Ewald3DLong(n, L; α = α, s = s)
    E_ewald = short_energy(short, poses, charges) + long_energy(long, poses, charges)

    # 8π, not 4π — see the oracle's Madelung testset for the derivation of the 2.
    M = -E_ewald * 8π * (a / 2) / n
    @test isapprox(M, 1.7475645946, rtol = 1e-4)
end

@testset "Ewald3D energy is independent of α" begin
    # The split point is arbitrary: short+long must be invariant. This is the
    # sharpest internal check on the relative normalisation of the two parts.
    Random.seed!(20260916)
    n = 20
    L = (12.0, 12.0, 12.0)
    poses = [SVector(rand() * L[1], rand() * L[2], rand() * L[3]) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]

    function total(α, s)
        short = EwaldShort(n, L; α = α, s = s)
        long  = Ewald3DLong(n, L; α = α, s = s)
        return short_energy(short, poses, charges) + long_energy(long, poses, charges)
    end

    # r_c = s/α = 5.71, 5.00, 4.44 — all < L/2 = 6. s = 4 caps accuracy near 1e-7,
    # so rtol is 1e-5 rather than 1e-6.
    E_ref = total(0.7, 4.0)
    @test isapprox(total(0.8, 4.0), E_ref, rtol = 1e-5)
    @test isapprox(total(0.9, 4.0), E_ref, rtol = 1e-5)
end

@testset "Ewald3D force matches -grad(energy)" begin
    Random.seed!(20260917)
    n = 10
    L = (12.0, 12.0, 12.0)
    poses = [SVector(rand() * L[1], rand() * L[2], rand() * L[3]) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]
    α, s = 0.8, 3.5      # r_c = 4.375 < 6
    short = EwaldShort(n, L; α = α, s = s)
    long  = Ewald3DLong(n, L; α = α, s = s)

    F = [zero(SVector{3,Float64}) for _ in 1:n]
    short_force!(F, short, poses, charges)
    long_force!(F, long, poses, charges)

    f = p -> short_energy(short, p, charges) + long_energy(long, p, charges)
    for i in 1:n, d in 1:3
        @test isapprox(F[i][d], -fd_gradient(f, poses, i, d; h = 1e-5),
                       rtol = 1e-4, atol = 1e-8)
    end
end

@testset "Ewald3D surface term responds to ϵ_inf" begin
    # A configuration with a net dipole: conducting (Inf) and vacuum (1.0) boundaries
    # must disagree, and the conducting case must drop the dipole term entirely.
    poses = [SVector(1.0, 4.0, 4.0), SVector(7.0, 4.0, 4.0)]
    charges = [1.0, -1.0]
    L = (8.0, 8.0, 8.0)
    l_cond = Ewald3DLong(2, L; α = 0.8, s = 3.0, ϵ_inf = Inf)   # r_c = 3.75 < 4
    l_vac  = Ewald3DLong(2, L; α = 0.8, s = 3.0, ϵ_inf = 1.0)
    @test !isapprox(long_energy(l_cond, poses, charges), long_energy(l_vac, poses, charges))

    # the difference is exactly the dipole term |P|²/(2Vϵ(2ϵ_inf+1))
    P = sum(charges[i] * poses[i] for i in 1:2)
    V = prod(L)
    @test isapprox(long_energy(l_vac, poses, charges) - long_energy(l_cond, poses, charges),
                   sum(abs2, P) / (2 * V * 1.0 * 3.0), rtol = 1e-10)
end

@testset "Ewald3DLong preserves Float32" begin
    L = (8.0f0, 8.0f0, 8.0f0)
    poses = [SVector(1.0f0, 2.0f0, 3.0f0), SVector(5.0f0, 6.0f0, 7.0f0)]
    charges = [1.0f0, -1.0f0]
    long = Ewald3DLong(2, L; α = 0.8f0, s = 3.0f0)   # r_c = 3.75 < 4
    @test long_energy(long, poses, charges) isa Float32
end
```

- [ ] **Step 2: Run to verify it fails**

Expected: FAIL with `UndefVarError: Ewald3DLong not defined`.

- [ ] **Step 3: Write the implementation**

Create `src/interactions/electrostatics/long_ewald3d.jl`:

```julia
"""
    Ewald3DLong(n_atoms, L; α, s, ϵ = 1.0, ϵ_inf = Inf)

Reciprocal-space part of the triply periodic Ewald sum, evaluated as a direct sum
over the k-set:

    E_l = 1/(2Vϵ) Σ_{0<|k|≤k_c} |ρ_k|² exp(−k²/4α²)/k²,   ρ_k = Σ_j q_j exp(i k·r_j)

plus the surface (dipole) term `|P|²/(2Vϵ(2ϵ_inf+1))` with `P = Σ_j q_j r_j`.
`ϵ_inf = Inf` is the conducting (tin-foil) boundary and drops the surface term;
finite `ϵ_inf` applies the correction for a medium of that permittivity at infinity.

Cost is `O(N·K)`. For large systems use the FINUFFT-backed `PME3D` instead, which
computes the same quantity with the same normalisation.
"""
struct Ewald3DLong{T}
    α::T
    r_c::T
    k_c::T
    ϵ::T
    ϵ_inf::T
    L::NTuple{3,T}
    n_atoms::Int
    k_set::Vector{NTuple{4,T}}
end

function Ewald3DLong(n_atoms::Int, L::NTuple{3,T}; α::T, s::T, ϵ::T = one(T),
                     ϵ_inf::T = T(Inf)) where {T}
    r_c, k_c = ewald_cutoffs(s, α)
    return Ewald3DLong{T}(α, r_c, k_c, ϵ, ϵ_inf, L, n_atoms, k_set_3D(k_c, L))
end

Base.show(io::IO, l::Ewald3DLong) =
    print(io, "Ewald3DLong(α = $(l.α), k_c = $(l.k_c), ϵ = $(l.ϵ), ϵ_inf = $(l.ϵ_inf), " *
              "$(length(l.k_set)) k-vectors)")

# Σ_j q_j exp(i k·r_j), accumulated in Complex{T} so Float32 and higher precisions survive.
@inline function _structure_factor(k_x::T, k_y::T, k_z::T, poses, charges,
                                   n::Int) where {T}
    ρ = zero(Complex{T})
    @inbounds for j in 1:n
        p = poses[j]
        ρ += charges[j] * cis(k_x * T(p[1]) + k_y * T(p[2]) + k_z * T(p[3]))
    end
    return ρ
end

@inline function _dipole(poses, charges, n::Int, ::Type{T}) where {T}
    P = zero(SVector{3,T})
    @inbounds for j in 1:n
        p = poses[j]
        P += charges[j] * SVector{3,T}(T(p[1]), T(p[2]), T(p[3]))
    end
    return P
end

function long_energy(long::Ewald3DLong{T}, poses, charges;
                     n_target::Int = long.n_atoms) where {T}
    α, n = long.α, long.n_atoms          # n = number of SOURCE charges
    V = long.L[1] * long.L[2] * long.L[3]

    E = zero(T)
    @inbounds for (k_x, k_y, k_z, k) in long.k_set
        ρ_src = _structure_factor(k_x, k_y, k_z, poses, charges, n)
        # real(conj(ρ_src) * ρ_tgt) collapses to abs2(ρ) when the sets coincide,
        # so there is one code path, not two.
        ρ_tgt = n_target == n ? ρ_src :
                _structure_factor(k_x, k_y, k_z, poses, charges, n_target)
        E += real(conj(ρ_src) * ρ_tgt) * exp(-k^2 / (4 * α^2)) / k^2
    end
    E /= (2 * V * long.ϵ)

    # Surface term. 1/(2*Inf+1) evaluates to zero, so the conducting case needs no branch.
    P_src = _dipole(poses, charges, n, T)
    P_tgt = n_target == n ? P_src : _dipole(poses, charges, n_target, T)
    E += dot(P_tgt, P_src) / (2 * V * long.ϵ * (2 * long.ϵ_inf + one(T)))

    return E
end

function long_force!(F::Vector{SVector{3,T}}, long::Ewald3DLong{T}, poses,
                     charges; n_target::Int = long.n_atoms) where {T}
    α, n = long.α, long.n_atoms          # n = number of SOURCE charges
    V = long.L[1] * long.L[2] * long.L[3]
    pref = one(T) / (V * long.ϵ)

    # Sources span all charges; the force is written only for the first n_target of
    # them. The coefficient is the full one — no 1/2 — matching the image-charge
    # convention, which is self-consistent because an image moves at twice the rate
    # of its source. Controller-verified against finite differences to 1e-7.
    @inbounds for (k_x, k_y, k_z, k) in long.k_set
        ρ_src = _structure_factor(k_x, k_y, k_z, poses, charges, n)
        D = exp(-k^2 / (4 * α^2)) / k^2
        kvec = SVector{3,T}(k_x, k_y, k_z)
        for i in 1:n_target
            p = poses[i]
            phase = cis(-(k_x * T(p[1]) + k_y * T(p[2]) + k_z * T(p[3])))
            F[i] -= (pref * charges[i] * D * imag(ρ_src * phase)) * kvec
        end
    end

    # Surface term force: F_i = -q_i P_src / (V ϵ (2ϵ_inf + 1))
    P_src = _dipole(poses, charges, n, T)
    surf = one(T) / (V * long.ϵ * (2 * long.ϵ_inf + one(T)))
    @inbounds for i in 1:n_target
        F[i] -= (surf * charges[i]) * P_src
    end

    return F
end
```

Add to `src/ExTinyMD.jl`:

```julia
include("interactions/electrostatics/long_ewald3d.jl")
```
```julia
export Ewald3DLong, long_energy, long_force!
```

- [ ] **Step 4: Run tests to verify they pass**

Expected: PASS. The α-independence test is the one to watch: if it fails, the relative normalisation of short and long is wrong, and the Madelung test may still pass at one particular α. Debug α-independence first.

- [ ] **Step 5: Commit**

```bash
git add src/ExTinyMD.jl src/interactions/electrostatics/long_ewald3d.jl test/
git commit -m "feat: add Ewald3D reciprocal-space solver

Direct k-sum with the dipole surface term. Validated by α-independence of the
total energy, the NaCl Madelung constant, and a finite-difference force check."
```

---

### Task 6: `EwaldInteraction` composite and the `Ewald3D` constructor

**Files:**
- Create: `src/interactions/electrostatics/ewald.jl`
- Create: `test/electrostatics/test_ewald.jl`
- Modify: `src/ExTinyMD.jl`, `test/runtests.jl`

**Interfaces:**
- Consumes: `EwaldShort`/`short_energy`/`short_force!` (Task 4), `Ewald3DLong`/`long_energy`/`long_force!` (Task 5).
- Produces:
  - `EwaldInteraction{T,S,L} <: AbstractInteraction` with fields `short::S`, `long::L`, `n_atoms::Int`, `force_buffer::Vector{SVector{3,T}}`
  - `coulomb_energy(inter, poses, charges; neighbor_list=nothing)::T`
  - `coulomb_force(inter, poses, charges; neighbor_list=nothing)::Vector{SVector{3,T}}`
  - `coulomb_force!(F, inter, poses, charges; neighbor_list=nothing)` — **zeroes `F` first**
  - `Ewald3D(n_atoms, L; α, s, ϵ=one(T), ϵ_inf=T(Inf))::EwaldInteraction`

Naming note: the core queries are `coulomb_energy` / `coulomb_force` rather than `energy` / `force`, because `ExTinyMD.energy` already exists with the four-argument MD signature and overloading it with a three-argument array form invites silent dispatch surprises. Task 9's adapter provides `ExTinyMD.energy`.

- [ ] **Step 1: Write the failing tests**

Create `test/electrostatics/test_ewald.jl`:

```julia
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
```

- [ ] **Step 2: Run to verify it fails**

Expected: FAIL with `UndefVarError: Ewald3D not defined`.

- [ ] **Step 3: Write the implementation**

Create `src/interactions/electrostatics/ewald.jl`:

```julia
"""
    EwaldInteraction(short, long, n_atoms)

An Ewald-split electrostatic interaction: a real-space kernel plus a
reciprocal-space solver. Composing them here rather than registering two separate
entries in `sys.interactions` keeps the pair together and hands the MD adapter a
single object.

Query it with [`coulomb_energy`](@ref) and [`coulomb_force`](@ref) for standalone
use, or add it to an `MDSys` and let the adapter drive it.
"""
struct EwaldInteraction{T, S, L} <: AbstractInteraction
    short::S
    long::L
    n_atoms::Int
    force_buffer::Vector{SVector{3,T}}
end

function EwaldInteraction(short::S, long::L, n_atoms::Int) where {S, L}
    T = typeof(short.α)
    return EwaldInteraction{T, S, L}(short, long, n_atoms,
                                     [zero(SVector{3,T}) for _ in 1:n_atoms])
end

Base.show(io::IO, i::EwaldInteraction) =
    print(io, "EwaldInteraction($(i.n_atoms) atoms)\n  short: $(i.short)\n  long:  $(i.long)")

"""
    coulomb_energy(interaction, poses, charges; neighbor_list = nothing) -> T

Total electrostatic energy. `poses` is AoS; no ExTinyMD type is required.
"""
function coulomb_energy(inter::EwaldInteraction{T}, poses, charges;
                        neighbor_list = nothing) where {T}
    return short_energy(inter.short, poses, charges; neighbor_list = neighbor_list) +
           long_energy(inter.long, poses, charges)
end

"""
    coulomb_force!(F, interaction, poses, charges; neighbor_list = nothing) -> F

Total electrostatic force, written into `F`. `F` is zeroed first, then the short-
and long-range parts accumulate into it.
"""
function coulomb_force!(F::Vector{SVector{3,T}}, inter::EwaldInteraction{T}, poses,
                        charges; neighbor_list = nothing) where {T}
    fill!(F, zero(SVector{3,T}))
    short_force!(F, inter.short, poses, charges; neighbor_list = neighbor_list)
    long_force!(F, inter.long, poses, charges)
    return F
end

"""
    coulomb_force(interaction, poses, charges; neighbor_list = nothing) -> Vector{SVector{3,T}}

Allocating form of [`coulomb_force!`](@ref). In an MD loop prefer the in-place
version, or let the adapter reuse `interaction.force_buffer`.
"""
function coulomb_force(inter::EwaldInteraction{T}, poses, charges;
                       neighbor_list = nothing) where {T}
    F = [zero(SVector{3,T}) for _ in 1:inter.n_atoms]
    return coulomb_force!(F, inter, poses, charges; neighbor_list = neighbor_list)
end

"""
    Ewald3D(n_atoms, L; α, s, ϵ = 1.0, ϵ_inf = Inf)

Standard Ewald summation for a triply periodic system. `α` splits real and
reciprocal space, `s` sets the accuracy (`r_c = s/α`, `k_c = 2αs`).

```jldoctest
julia> using StaticArrays

julia> inter = Ewald3D(2, (10.0, 10.0, 10.0); α = 1.0, s = 4.0);

julia> poses = [SVector(0.0, 0.0, 0.0), SVector(5.0, 0.0, 0.0)];

julia> round(coulomb_energy(inter, poses, [1.0, -1.0]); digits = 6)
-0.023873
```
"""
function Ewald3D(n_atoms::Int, L::NTuple{3,T}; α::T, s::T, ϵ::T = one(T),
                 ϵ_inf::T = T(Inf)) where {T}
    short = EwaldShort(n_atoms, L; α = α, s = s, ϵ = ϵ, convention = Periodic3D())
    long  = Ewald3DLong(n_atoms, L; α = α, s = s, ϵ = ϵ, ϵ_inf = ϵ_inf)
    return EwaldInteraction(short, long, n_atoms)
end
```

The `jldoctest` value above is a placeholder that **must be replaced with the real computed output**: run the snippet, paste the actual number, and confirm doctests pass (Task 10 wires doctests into CI). Do not invent it.

Add to `src/ExTinyMD.jl`:

```julia
include("interactions/electrostatics/ewald.jl")
```
```julia
export EwaldInteraction, coulomb_energy, coulomb_force, coulomb_force!, Ewald3D
```

- [ ] **Step 4: Run tests to verify they pass**

Expected: PASS.

- [ ] **Step 5: Fill in the real doctest output**

```
julia --project=. -e 'using ExTinyMD, StaticArrays; inter = Ewald3D(2, (10.0,10.0,10.0); α=1.0, s=4.0); poses=[SVector(0.0,0.0,0.0), SVector(5.0,0.0,0.0)]; println(round(coulomb_energy(inter, poses, [1.0,-1.0]); digits=6))'
```

Paste the printed value into the docstring.

- [ ] **Step 6: Commit**

```bash
git add src/ExTinyMD.jl src/interactions/electrostatics/ewald.jl test/
git commit -m "feat: add EwaldInteraction composite and Ewald3D constructor"
```

---

### Task 7: `Ewald2DLong` — quasi-2D reciprocal-space sum

**Files:**
- Create: `src/interactions/electrostatics/long_ewald2d.jl`
- Create: `test/electrostatics/test_long_ewald2d.jl`
- Modify: `src/ExTinyMD.jl`, `src/interactions/electrostatics/ewald.jl` (add `Ewald2D`), `test/runtests.jl`

**Interfaces:**
- Consumes: `k_set_2D` (Task 3), `EwaldShort` with `PeriodicQ2D` (Task 4), `EwaldInteraction` (Task 6).
- Produces:
  - `Ewald2DLong(n_atoms::Int, L::NTuple{3,T}; α, s, ϵ=one(T))`
  - `long_energy(long::Ewald2DLong{T}, poses, charges; n_target::Int = long.n_atoms)::T`
  - `long_force!(F, long::Ewald2DLong{T}, poses, charges; n_target::Int = long.n_atoms)` — accumulates

`n_target` has the same meaning as in Task 5: the `i` (target) loops run over
`1:n_target` while the `j` (source) loops run over all `1:long.n_atoms`. The default makes
them coincide. ICM (Task 8) is the only caller that passes a smaller value.
  - `Ewald2D(n_atoms, L; α, s, ϵ=one(T))::EwaldInteraction`

- [ ] **Step 1: Write the failing tests**

Create `test/electrostatics/test_long_ewald2d.jl`:

```julia
@testset "Ewald2D energy is independent of α" begin
    Random.seed!(20260922)
    n = 16
    L = (6.0, 6.0, 20.0)
    # a slab: charges confined well inside the z extent
    poses = [SVector(rand() * L[1], rand() * L[2], 5.0 + 10.0 * rand()) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]

    total(α, s) = coulomb_energy(Ewald2D(n, L; α = α, s = s), poses, charges)
    # Only L[1], L[2] bound r_c under PeriodicQ2D: r_c = s/α = 2.86, 2.67, 2.50,
    # all < min(Lx,Ly)/2 = 3.
    E_ref = total(1.4, 4.0)
    @test isapprox(total(1.5, 4.0), E_ref, rtol = 1e-5)
    @test isapprox(total(1.6, 4.0), E_ref, rtol = 1e-5)
end

@testset "Ewald2D energy against quasi-2D direct sum" begin
    Random.seed!(20260923)
    n = 8
    L = (5.0, 5.0, 30.0)
    poses = [SVector(rand() * L[1], rand() * L[2], 10.0 + 10.0 * rand()) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]

    E_ewald  = coulomb_energy(Ewald2D(n, L; α = 1.7, s = 4.0), poses, charges)   # r_c = 2.35 < 2.5

    # The RAW 2D lattice sum is useless as a reference here: it converges as
    # 1/n_shell and is still 3.5% off at n_shell = 40. Use the extrapolated form,
    # which reaches ~3e-4 at (30, 60). Controller-measured, raw vs this Ewald value:
    #   n_shell = 10, 20, 30, 40, 60, 80  ->  13.5%, 6.9%, 4.6%, 3.5%, 2.3%, 1.7%
    E_direct = naive_energy_Q2D_extrap(poses, charges, L, 30, 60)
    @test isapprox(E_ewald, E_direct, rtol = 1e-3)

    # and confirm the extrapolation is doing the work, not a loose tolerance
    @test abs(E_ewald - E_direct) <
          abs(E_ewald - naive_energy_Q2D(poses, charges, L, 60))
end

@testset "Ewald2D force matches -grad(energy)" begin
    Random.seed!(20260924)
    n = 8
    L = (6.0, 6.0, 20.0)
    poses = [SVector(rand() * L[1], rand() * L[2], 5.0 + 10.0 * rand()) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]
    inter = Ewald2D(n, L; α = 1.3, s = 3.5)      # r_c = 2.69 < 3

    F = coulomb_force(inter, poses, charges)
    f = p -> coulomb_energy(inter, p, charges)
    for i in 1:n, d in 1:3
        @test isapprox(F[i][d], -fd_gradient(f, poses, i, d; h = 1e-5),
                       rtol = 1e-4, atol = 1e-8)
    end
end

@testset "Ewald2D net in-plane force vanishes" begin
    Random.seed!(20260925)
    n = 10
    L = (6.0, 6.0, 20.0)
    poses = [SVector(rand() * L[1], rand() * L[2], 5.0 + 10.0 * rand()) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]
    F = coulomb_force(Ewald2D(n, L; α = 1.3, s = 3.5), poses, charges)
    @test isapprox(sum(f -> f[1], F), 0.0, atol = 1e-9)
    @test isapprox(sum(f -> f[2], F), 0.0, atol = 1e-9)
    @test isapprox(sum(f -> f[3], F), 0.0, atol = 1e-9)
end

@testset "Ewald2DLong does not overflow for a tall box" begin
    # exp(k*z_ij) overflows for large k*z; the paired erfc underflows to zero at the
    # same time, so the product must be computed as zero rather than Inf*0 = NaN.
    L = (4.0, 4.0, 2000.0)
    poses = [SVector(1.0, 1.0, 10.0), SVector(2.0, 2.0, 1990.0)]
    charges = [1.0, -1.0]
    inter = Ewald2D(2, L; α = 1.6, s = 3.0)      # r_c = 1.875 < 2
    E = coulomb_energy(inter, poses, charges)
    @test isfinite(E)
    F = coulomb_force(inter, poses, charges)
    @test all(all(isfinite, f) for f in F)
end
```

- [ ] **Step 2: Run to verify it fails**

Expected: FAIL with `UndefVarError: Ewald2D not defined`.

- [ ] **Step 3: Write the implementation**

Create `src/interactions/electrostatics/long_ewald2d.jl`:

```julia
"""
    Ewald2DLong(n_atoms, L; α, s, ϵ = 1.0)

Reciprocal-space part of the exact Ewald sum for a slab periodic in x and y and
free in z. For each in-plane wavevector `k`:

    E_k  = 1/ϵ Σ_i Σ_j q_i q_j cos(k·ρ_ij)
             [ e^{k z_ij} erfc(k/2α + α z_ij) + e^{−k z_ij} erfc(k/2α − α z_ij) ]
             / (8 L_x L_y k)

    E_k0 = −1/ϵ Σ_i Σ_j q_i q_j [ e^{−(α z_ij)²}/(α√π) + z_ij erf(α z_ij) ]
             / (4 L_x L_y)

!!! note "Prefactor"
    The long-range prefactor is `1/ϵ`, not `1/(4πϵ)` — the `4π` is already folded
    into the expressions above. The short-range part of the same method does carry
    `1/(4πϵ)`. This asymmetry is deliberate; α-independence of the total energy is
    the test that pins it.

Cost is `O(N²K)`: this is the *exact* 2D sum, intended as the accuracy reference
for quasi-2D systems rather than as a production method for large `N`.
"""
struct Ewald2DLong{T}
    α::T
    r_c::T
    k_c::T
    ϵ::T
    L::NTuple{3,T}
    n_atoms::Int
    k_set::Vector{NTuple{3,T}}
end

function Ewald2DLong(n_atoms::Int, L::NTuple{3,T}; α::T, s::T, ϵ::T = one(T)) where {T}
    r_c, k_c = ewald_cutoffs(s, α)
    return Ewald2DLong{T}(α, r_c, k_c, ϵ, L, n_atoms, k_set_2D(k_c, L))
end

Base.show(io::IO, l::Ewald2DLong) =
    print(io, "Ewald2DLong(α = $(l.α), k_c = $(l.k_c), ϵ = $(l.ϵ), " *
              "$(length(l.k_set)) k-vectors)")

# exp(±k z) erfc(k/2α ± α z), guarded. For large positive argument the exp overflows
# while the erfc underflows; the product tends to zero, so return zero rather than
# letting Inf * 0 produce NaN.
@inline function _exp_erfc(kz::T, arg::T) where {T}
    kz > T(600) && return zero(T)
    return exp(kz) * erfc(arg)
end

function long_energy(long::Ewald2DLong{T}, poses, charges;
                     n_target::Int = long.n_atoms) where {T}
    α, n = long.α, long.n_atoms          # n = number of SOURCE charges
    A = long.L[1] * long.L[2]

    E = zero(T)

    # k = 0 term. Targets i run to n_target, sources j over all n.
    @inbounds for i in 1:n_target, j in 1:n
        z = T(poses[i][3]) - T(poses[j][3])
        E -= charges[i] * charges[j] *
             (exp(-(α * z)^2) / (α * sqrt(T(π))) + z * erf(α * z)) / (4 * A)
    end

    # k != 0 terms
    @inbounds for (k_x, k_y, k) in long.k_set
        acc = zero(T)
        for i in 1:n_target
            p_i = poses[i]
            for j in 1:n
                p_j = poses[j]
                x = T(p_i[1]) - T(p_j[1])
                y = T(p_i[2]) - T(p_j[2])
                z = T(p_i[3]) - T(p_j[3])
                g = _exp_erfc(k * z, k / (2α) + α * z) +
                    _exp_erfc(-k * z, k / (2α) - α * z)
                acc += charges[i] * charges[j] * cos(k_x * x + k_y * y) * g
            end
        end
        E += acc / (8 * A * k)
    end

    return E / long.ϵ
end

function long_force!(F::Vector{SVector{3,T}}, long::Ewald2DLong{T}, poses,
                     charges; n_target::Int = long.n_atoms) where {T}
    α, n = long.α, long.n_atoms          # n = number of SOURCE charges
    A = long.L[1] * long.L[2]
    ϵ = long.ϵ

    # k = 0 term: F_i,z = q_i Σ_j q_j erf(α z_ij) / (2 L_x L_y ϵ)
    @inbounds for i in 1:n_target
        fz = zero(T)
        for j in 1:n
            z = T(poses[i][3]) - T(poses[j][3])
            fz += charges[j] * erf(α * z)
        end
        F[i] += SVector{3,T}(zero(T), zero(T), charges[i] * fz / (2 * A * ϵ))
    end

    # k != 0 terms
    @inbounds for (k_x, k_y, k) in long.k_set
        for i in 1:n_target
            p_i = poses[i]
            sx = zero(T); sy = zero(T); sz = zero(T)
            for j in 1:n
                p_j = poses[j]
                x = T(p_i[1]) - T(p_j[1])
                y = T(p_i[2]) - T(p_j[2])
                z = T(p_i[3]) - T(p_j[3])
                qq = charges[i] * charges[j]
                phase = k_x * x + k_y * y

                ee_p = _exp_erfc(k * z, k / (2α) + α * z)    # e^{kz}  erfc(k/2α + αz)
                ee_m = _exp_erfc(-k * z, k / (2α) - α * z)   # e^{-kz} erfc(k/2α - αz)

                sum_xy = -qq * sin(phase) * (ee_p + ee_m)
                sx += k_x * sum_xy / k
                sy += k_y * sum_xy / k

                # d/dz of (e^{kz} erfc(k/2α+αz) + e^{-kz} erfc(k/2α-αz))
                gauss_p = k * z > T(600) ? zero(T) :
                          exp(k * z) * exp(-(k / (2α) + α * z)^2)
                gauss_m = -k * z > T(600) ? zero(T) :
                          exp(-k * z) * exp(-(k / (2α) - α * z)^2)
                sz += qq * cos(phase) * (k * ee_p - k * ee_m -
                                         2α / sqrt(T(π)) * gauss_p +
                                         2α / sqrt(T(π)) * gauss_m) / k
            end
            F[i] -= SVector{3,T}(sx, sy, sz) / (4 * A * ϵ)
        end
    end

    return F
end
```

Then append the `Ewald2D` constructor to `src/interactions/electrostatics/ewald.jl`:

```julia
"""
    Ewald2D(n_atoms, L; α, s, ϵ = 1.0)

Exact Ewald summation for a slab periodic in x and y and free in z. `L[3]` is the
extent used for the real-space neighbour search, not a period.

`O(N²K)` — use it as the quasi-2D accuracy reference, and for production runs on
large systems reach for `QuasiEwald.jl` or `SoEwald2D.jl`.
"""
function Ewald2D(n_atoms::Int, L::NTuple{3,T}; α::T, s::T, ϵ::T = one(T)) where {T}
    short = EwaldShort(n_atoms, L; α = α, s = s, ϵ = ϵ, convention = PeriodicQ2D())
    long  = Ewald2DLong(n_atoms, L; α = α, s = s, ϵ = ϵ)
    return EwaldInteraction(short, long, n_atoms)
end
```

Add the include (before `ewald.jl`, since `ewald.jl` now references `Ewald2DLong`) and exports:

```julia
include("interactions/electrostatics/long_ewald2d.jl")
```
```julia
export Ewald2DLong, Ewald2D
```

- [ ] **Step 4: Run tests to verify they pass**

Expected: PASS. The α-independence test is again the key diagnostic — if only it fails, the `1/ϵ` vs `1/(4πϵ)` split is wrong somewhere.

- [ ] **Step 5: Commit**

```bash
git add src/ExTinyMD.jl src/interactions/electrostatics/ test/
git commit -m "feat: add exact Ewald2D reciprocal-space solver

Slab geometry periodic in x,y. Overflow-guarded exp*erfc products, forces
threaded over the outer loop rather than launching a thread region per particle
per wavevector as the reference implementation did."
```

---

### Task 8: ICM and ELC

**Files:**
- Create: `src/interactions/electrostatics/icm.jl`
- Create: `test/electrostatics/test_icm.jl`
- Modify: `src/ExTinyMD.jl`, `test/runtests.jl`

**Interfaces:**
- Consumes: `Ewald2DLong` (Task 7), `Ewald3DLong` (Task 5), `min_image_disp` (Task 3).
- Produces:
  - `icm_reflect!(ref_poses, ref_charges, γ::Tuple{T,T}, L, N_image, poses, charges)` → fills buffers, returns the number of entries written
  - `ICMShort(n_atoms, L; α, s, ϵ, N_image)` with `icm_short_energy` / `icm_short_force!`
  - `ICM{T,L} <: AbstractInteraction` with fields `long`, `short`, `γ`, `N_image`, `n_atoms`, `L`, `ref_poses`, `ref_charges`, `ref_force`, `elc::Bool`, `N_pad`
  - `coulomb_energy(icm::ICM, poses, charges; neighbor_list=nothing)` / `coulomb_force!` / `coulomb_force`
  - `ICMEwald2D(n_atoms, L; α, s, γ, N_image, ϵ=one(T))`
  - `ICMEwald3D(n_atoms, L; α, s, γ, N_image, N_pad, ϵ=one(T))`

Layout contract, which the force fold-back depends on: `ref_poses[1:n]` are the real particles in input order; entries `n + 2(i-1)*N_image + 1 : n + 2*i*N_image` are particle `i`'s images, interleaved up, down, up, down. Total length `n*(1 + 2*N_image)`.

- [ ] **Step 1: Write the failing tests**

Create `test/electrostatics/test_icm.jl`:

```julia
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
```

- [ ] **Step 2: Run to verify it fails**

Expected: FAIL with `UndefVarError: icm_reflect! not defined`.

- [ ] **Step 3: Implement `icm_reflect!`**

Create `src/interactions/electrostatics/icm.jl`, starting with the reflection:

```julia
"""
    icm_reflect!(ref_poses, ref_charges, γ, L, N_image, poses, charges) -> Int

Build the image-charge series for a slab confined between two dielectric walls,
writing into preallocated buffers and returning the number of entries written,
`n*(1 + 2*N_image)`.

Layout, which the force fold-back relies on:

- `ref_poses[1:n]` — the real particles, in input order.
- `ref_poses[n + 2(i-1)N_image + 1 : n + 2i·N_image]` — particle `i`'s images,
  interleaved up, down, up, down, …

`γ = (γ_up, γ_down)` are the dielectric contrast ratios at the two walls,
conventionally `(ϵ_mid − ϵ_out)/(ϵ_mid + ϵ_out)` for each; `γ = 0` is a matched
wall. The recurrence reflects each image off the opposite wall in turn:

    z_up[1]   = 2L_z − z          q_up[1]   = γ_up   q
    z_down[1] = −z                q_down[1] = γ_down q
    z_up[m]   = 2L_z − z_down[m−1]    q_up[m]   = γ_up   q_down[m−1]
    z_down[m] = −z_up[m−1]            q_down[m] = γ_down q_up[m−1]
"""
function icm_reflect!(ref_poses::Vector{SVector{3,T}}, ref_charges::Vector{T},
                      γ::Tuple{T,T}, L::NTuple{3,T}, N_image::Int,
                      poses, charges) where {T}
    n = length(charges)
    γ_up, γ_down = γ
    L_z = L[3]

    @inbounds for i in 1:n
        p = poses[i]
        ref_poses[i] = SVector{3,T}(T(p[1]), T(p[2]), T(p[3]))
        ref_charges[i] = charges[i]
    end

    idx = n
    @inbounds for i in 1:n
        p = poses[i]
        x, y, z = T(p[1]), T(p[2]), T(p[3])
        q = charges[i]

        z_up, z_down = 2 * L_z - z, -z
        q_up, q_down = γ_up * q, γ_down * q

        for m in 1:N_image
            if m > 1
                z_up, z_down = 2 * L_z - z_down, -z_up
                q_up, q_down = γ_up * q_down, γ_down * q_up
            end
            idx += 1
            ref_poses[idx] = SVector{3,T}(x, y, z_up)
            ref_charges[idx] = q_up
            idx += 1
            ref_poses[idx] = SVector{3,T}(x, y, z_down)
            ref_charges[idx] = q_down
        end
    end

    return idx
end
```

Careful: the recurrence updates `z_up`/`z_down` from the *previous* pair, so both must be
computed from the old values simultaneously. The tuple assignment above does that; a
sequential pair of assignments would not.

- [ ] **Step 4: Run the reflect tests**

Run only the reflect testset first. Expected: PASS for `icm_reflect!`, remaining ICM testsets still failing.

- [ ] **Step 5: Implement `ICMShort`, `ICM`, ELC, and the constructors**

Append to `src/interactions/electrostatics/icm.jl`:

```julia
"""
    ICMShort(n_atoms, L; α, s, ϵ = 1.0, N_image)

Real-space kernel for ICM. Unlike [`EwaldShort`](@ref) it must distinguish
real–real from real–image pairs: real–image energies carry a factor ½, real–image
forces accumulate on the real particle only, and the self-energy is summed over
real particles only.
"""
mutable struct ICMShort{T, TC}
    α::T
    r_c::T
    ϵ::T
    L::NTuple{3,T}
    n_atoms::Int
    N_image::Int
    cell_list::TC
    n_ref::Int
end

function ICMShort(n_atoms::Int, L::NTuple{3,T}; α::T, s::T, ϵ::T = one(T),
                  N_image::Int) where {T}
    r_c, _ = ewald_cutoffs(s, α)
    n_ref = n_atoms * (1 + 2 * N_image)
    buf = [zero(SVector{3,T}) for _ in 1:n_ref]
    # The reflected stack spans (2 N_image + 1) L_z in z; pad by 2 r_c so the
    # non-periodic z direction finds no spurious images.
    unitcell = SVector{3,T}(L[1], L[2], (2 * N_image + 1) * L[3] + 2 * r_c)
    cell_list = InPlaceNeighborList(xpositions = buf, cutoff = r_c, unitcell = unitcell,
                                    parallel = true)
    return ICMShort{T, typeof(cell_list)}(α, r_c, ϵ, L, n_atoms, N_image, cell_list, n_ref)
end

function icm_short_energy(short::ICMShort{T}, ref_poses::Vector{SVector{3,T}},
                          ref_charges::Vector{T}) where {T}
    n, α, r_c = short.n_atoms, short.α, short.r_c
    update!(short.cell_list, xpositions = ref_poses)
    nb = neighborlist!(short.cell_list)

    E = zero(T)
    @inbounds for (i, j, r) in nb
        (r < r_c && r > zero(T)) || continue
        both_real = (i <= n) && (j <= n)
        either_real = (i <= n) || (j <= n)
        either_real || continue
        w = both_real ? one(T) : T(0.5)
        E += w * ref_charges[i] * ref_charges[j] * erfc(α * r) / r
    end

    @inbounds for i in 1:n
        E -= ref_charges[i]^2 * α / sqrt(T(π))
    end

    return E / (4π * short.ϵ)
end

function icm_short_force!(F::Vector{SVector{3,T}}, short::ICMShort{T},
                          ref_poses::Vector{SVector{3,T}},
                          ref_charges::Vector{T}) where {T}
    n, α, r_c = short.n_atoms, short.α, short.r_c
    L, ϵ = short.L, short.ϵ
    pref = one(T) / (4π * ϵ)
    update!(short.cell_list, xpositions = ref_poses)
    nb = neighborlist!(short.cell_list)

    @inbounds for (i, j, r) in nb
        (r < r_c && r > zero(T)) || continue
        i_real = i <= n
        j_real = j <= n
        (i_real || j_real) || continue
        # Images are periodic in x,y and free in z, like the slab itself.
        d = min_image_disp(ref_poses[i], ref_poses[j], L, PeriodicQ2D())
        F_ij = _short_pair_dEdr(ref_charges[i], ref_charges[j], α, r) * d / r * pref
        # Real-image pairs push only on the real member (spec §7).
        i_real && (F[i] += F_ij)
        j_real && (F[j] -= F_ij)
    end
    return F
end

"""
    ICM(long, short, γ, N_image, n_atoms, L; elc = false, N_pad = 0)

Image-charge method for a slab between two dielectric walls. The long-range solver
`long` is evaluated on the reflected configuration and its forces are kept for the
real particles only; `short` is the split real-space kernel.

Set `elc = true` with `N_pad > 0` for the triply periodic + electrostatic layer
correction route, in which case `long` must be an [`Ewald3DLong`](@ref) built for
the z-padded box.

!!! warning "Force convention"
    Image positions depend on the real particles' z coordinates. The convention
    here — forces accumulated on real indices only — is ported verbatim from
    `EwaldSummations.jl` rather than re-derived, so published results remain
    reproducible. It is pinned by a finite-difference test.
"""
struct ICM{T, L, S} <: AbstractInteraction
    long::L
    short::S
    γ::Tuple{T,T}
    N_image::Int
    n_atoms::Int
    L::NTuple{3,T}
    elc::Bool
    N_pad::Int
    ref_poses::Vector{SVector{3,T}}
    ref_charges::Vector{T}
    ref_force::Vector{SVector{3,T}}
end

function ICM(long::L, short::S, γ::Tuple{T,T}, N_image::Int, n_atoms::Int,
             Lbox::NTuple{3,T}; elc::Bool = false, N_pad::Int = 0) where {T, L, S}
    n_ref = n_atoms * (1 + 2 * N_image)
    return ICM{T, L, S}(long, short, γ, N_image, n_atoms, Lbox, elc, N_pad,
                        [zero(SVector{3,T}) for _ in 1:n_ref], zeros(T, n_ref),
                        [zero(SVector{3,T}) for _ in 1:n_ref])
end

Base.show(io::IO, i::ICM) =
    print(io, "ICM($(i.n_atoms) atoms, γ = $(i.γ), N_image = $(i.N_image)" *
              (i.elc ? ", ELC with N_pad = $(i.N_pad)" : "") * ")\n  long: $(i.long)")

# ELC slab correction:
#   E = -1/(4πϵ) · π/(L_x L_y (2N_pad+1) L_z) Σ_i q_i Σ_j q_j (z_i - z_j)²
# summed over real i and all reflected j.
function _elc_energy(icm::ICM{T}, n_ref::Int) where {T}
    L = icm.L
    pref = -T(π) / (L[1] * L[2] * (2 * icm.N_pad + 1) * L[3]) / (4π * icm.long.ϵ)
    E = zero(T)
    @inbounds for i in 1:icm.n_atoms
        z_i = icm.ref_poses[i][3]
        t = zero(T)
        for j in 1:n_ref
            t += icm.ref_charges[j] * (z_i - icm.ref_poses[j][3])^2
        end
        E += icm.ref_charges[i] * t
    end
    return pref * E
end

function _elc_force!(F::Vector{SVector{3,T}}, icm::ICM{T}, n_ref::Int) where {T}
    L = icm.L
    pref = T(π) / (L[1] * L[2] * (2 * icm.N_pad + 1) * L[3]) / (4π * icm.long.ϵ)
    @inbounds for i in 1:icm.n_atoms
        z_i = icm.ref_poses[i][3]
        t = zero(T)
        for j in 1:n_ref
            t += 4 * icm.ref_charges[j] * (z_i - icm.ref_poses[j][3])
        end
        F[i] += SVector{3,T}(zero(T), zero(T), pref * icm.ref_charges[i] * t)
    end
    return F
end

function coulomb_energy(icm::ICM{T}, poses, charges; neighbor_list = nothing) where {T}
    n_ref = icm_reflect!(icm.ref_poses, icm.ref_charges, icm.γ, icm.L, icm.N_image,
                         poses, charges)
    # n_target = n_atoms is load-bearing: the long-range sum runs over REAL
    # particles as targets against ALL reflected charges as sources. Letting the
    # targets range over the images too would add unphysical image-image
    # self-interaction. Controller-verified: with real-only targets ICM+Ewald2D and
    # ICM+Ewald3D+ELC agree to 4.7e-8; with all-reflected targets they are 5.8% apart.
    E = icm_short_energy(icm.short, icm.ref_poses, icm.ref_charges) +
        long_energy(icm.long, icm.ref_poses, icm.ref_charges;
                    n_target = icm.n_atoms)
    icm.elc && (E += _elc_energy(icm, n_ref))
    return E
end

function coulomb_force!(F::Vector{SVector{3,T}}, icm::ICM{T}, poses, charges;
                        neighbor_list = nothing) where {T}
    n_ref = icm_reflect!(icm.ref_poses, icm.ref_charges, icm.γ, icm.L, icm.N_image,
                         poses, charges)
    fill!(icm.ref_force, zero(SVector{3,T}))
    icm_short_force!(icm.ref_force, icm.short, icm.ref_poses, icm.ref_charges)
    long_force!(icm.ref_force, icm.long, icm.ref_poses, icm.ref_charges;
                n_target = icm.n_atoms)
    icm.elc && _elc_force!(icm.ref_force, icm, n_ref)

    # Fold back: keep the real particles only.
    @inbounds for i in 1:icm.n_atoms
        F[i] = icm.ref_force[i]
    end
    return F
end

function coulomb_force(icm::ICM{T}, poses, charges; neighbor_list = nothing) where {T}
    F = [zero(SVector{3,T}) for _ in 1:icm.n_atoms]
    return coulomb_force!(F, icm, poses, charges; neighbor_list = neighbor_list)
end

"""
    ICMEwald2D(n_atoms, L; α, s, γ, N_image, ϵ = 1.0)

Image-charge method combined with the exact Ewald2D sum, for a slab confined by
two dielectric interfaces.
"""
function ICMEwald2D(n_atoms::Int, L::NTuple{3,T}; α::T, s::T, γ::Tuple{T,T},
                    N_image::Int, ϵ::T = one(T)) where {T}
    n_ref = n_atoms * (1 + 2 * N_image)
    long  = Ewald2DLong(n_ref, L; α = α, s = s, ϵ = ϵ)
    short = ICMShort(n_atoms, L; α = α, s = s, ϵ = ϵ, N_image = N_image)
    return ICM(long, short, γ, N_image, n_atoms, L)
end

"""
    ICMEwald3D(n_atoms, L; α, s, γ, N_image, N_pad, ϵ = 1.0)

Image-charge method combined with Ewald3D plus the electrostatic layer correction.
The slab is embedded in a box padded to `(2*N_pad + 1) * L[3]` in z and treated as
triply periodic; the ELC slab term removes the spurious interaction between
periodic replicas.
"""
function ICMEwald3D(n_atoms::Int, L::NTuple{3,T}; α::T, s::T, γ::Tuple{T,T},
                    N_image::Int, N_pad::Int, ϵ::T = one(T)) where {T}
    n_ref = n_atoms * (1 + 2 * N_image)
    L_pad = (L[1], L[2], (2 * N_pad + 1) * L[3])
    long  = Ewald3DLong(n_ref, L_pad; α = α, s = s, ϵ = ϵ, ϵ_inf = T(Inf))
    short = ICMShort(n_atoms, L; α = α, s = s, ϵ = ϵ, N_image = N_image)
    return ICM(long, short, γ, N_image, n_atoms, L; elc = true, N_pad = N_pad)
end
```

Note `long` is constructed with `n_ref`, not `n_atoms`: the long-range solver sees the whole
reflected configuration, and its `n_atoms` field is what its loops iterate over.

Add to `src/ExTinyMD.jl` (after `long_ewald2d.jl`, before `ewald.jl`):

```julia
include("interactions/electrostatics/icm.jl")
```
```julia
export ICM, ICMShort, ICMEwald2D, ICMEwald3D
```

- [ ] **Step 6: Run tests to verify they pass**

Expected: PASS. Two notes:

- The `γ = 0` reduction test is the first thing to check if others fail — it isolates the reflection and the real–real path from the image machinery.
- If only the ICM finite-difference test fails, that is the spec §7 risk materialising. **Stop and report it.** Record the observed relative discrepancy and which components disagree. Do not adjust the formula.

- [ ] **Step 7: Commit**

```bash
git add src/ExTinyMD.jl src/interactions/electrostatics/icm.jl test/
git commit -m "feat: add ICM and ELC for dielectrically confined slabs

Image-charge reflection with preallocated buffers, a split real-space kernel
shared by every ICM method, and the ELC slab correction. The force convention is
ported from EwaldSummations and pinned by a finite-difference test."
```

---

### Task 9: MD adapter

**Files:**
- Create: `src/interactions/electrostatics/adapter.jl`
- Create: `test/electrostatics/test_adapter.jl`
- Modify: `src/ExTinyMD.jl`, `test/runtests.jl`

**Interfaces:**
- Consumes: `EwaldInteraction`, `ICM`, `coulomb_energy`, `coulomb_force!`.
- Produces:
  - `ExTinyMD.energy(inter::Union{EwaldInteraction,ICM}, finder, sys::MDSys, info::SimulationInfo)::T`
  - `ExTinyMD.update_acceleration!(inter::Union{EwaldInteraction,ICM}, finder, sys, info)`
  - `gather_charges!(buf, sys, info)`, `gather_positions!(buf, info)`

- [ ] **Step 1: Write the failing tests**

Create `test/electrostatics/test_adapter.jl`:

```julia
function _charged_system(n, L; charge = 1.0)
    boundary = Boundary((L, L, L), (1, 1, 1))
    atoms = Vector{Atom{Float64}}()
    for i in 1:(n ÷ 2)
        push!(atoms, Atom(type = 1, mass = 1.0, charge = charge))
    end
    for i in (n ÷ 2 + 1):n
        push!(atoms, Atom(type = 2, mass = 1.0, charge = -charge))
    end
    info = SimulationInfo(n, atoms, (0.0, L, 0.0, L, 0.0, L), boundary;
                          min_r = 1.0, temp = 1.0)
    info.running_step = 1
    return boundary, atoms, info
end

@testset "adapter energy matches the core API" begin
    Random.seed!(20260930)
    n, L = 20, 10.0
    boundary, atoms, info = _charged_system(n, L)
    inter = Ewald3D(n, (L, L, L); α = 1.0, s = 4.0)   # r_c = 4.0 < L/2 = 5
    finder = CellList3D(info, inter.short.r_c, boundary, 1)

    sys = MDSys(n_atoms = n, atoms = atoms, boundary = boundary,
                interactions = [(inter, finder)],
                loggers = [TemperatureLogger(100; output = false)],
                simulator = VerletProcess(dt = 0.001))

    poses = [SVector(p.position[1], p.position[2], p.position[3])
             for p in info.particle_info]
    charges = [atoms[p.id].charge for p in info.particle_info]

    @test isapprox(energy(inter, finder, sys, info),
                   coulomb_energy(inter, poses, charges), rtol = 1e-10)
end

@testset "adapter writes acceleration = force/mass" begin
    Random.seed!(20260931)
    n, L = 16, 10.0
    boundary, atoms, info = _charged_system(n, L)
    # give the two species different masses so a missing division shows up
    atoms = [Atom(type = a.type, mass = a.type == 1 ? 1.0 : 4.0, charge = a.charge)
             for a in atoms]
    inter = Ewald3D(n, (L, L, L); α = 1.0, s = 4.0)   # r_c = 4.0 < L/2 = 5
    finder = CellList3D(info, inter.short.r_c, boundary, 1)
    sys = MDSys(n_atoms = n, atoms = atoms, boundary = boundary,
                interactions = [(inter, finder)],
                loggers = [TemperatureLogger(100; output = false)],
                simulator = VerletProcess(dt = 0.001))

    poses = [SVector(p.position[1], p.position[2], p.position[3])
             for p in info.particle_info]
    charges = [atoms[p.id].charge for p in info.particle_info]
    F = coulomb_force(inter, poses, charges)

    for p in info.particle_info
        p.acceleration = Point(0.0, 0.0, 0.0)
    end
    ExTinyMD.update_acceleration!(inter, finder, sys, info)

    for (slot, p) in enumerate(info.particle_info)
        m = atoms[p.id].mass
        for d in 1:3
            @test isapprox(p.acceleration[d], F[slot][d] / m, rtol = 1e-10)
        end
    end
end

@testset "adapter is correct when slot order differs from id order" begin
    # In stock ExTinyMD `particle_info[i].id == i`, so a gather that confuses slot
    # with id looks correct forever. Permute the mapping so the two differ and the
    # confusion becomes observable: charges must follow ids, positions must follow
    # slots. `substrate_lennard_jones.jl` already relies on this indirection via
    # `info.id_dict`, so it is a real invariant, not a hypothetical one.
    Random.seed!(20260934)
    n, L = 12, 10.0
    boundary, atoms, info = _charged_system(n, L)

    # give every id a distinct charge so a slot/id mix-up cannot cancel out
    atoms = [Atom(type = a.type, mass = 1.0 + 0.1 * i, charge = (isodd(i) ? 1.0 : -1.0) * (1 + 0.01 * i))
             for (i, a) in enumerate(atoms)]

    # reverse the slot order, keeping ids attached to their particles
    reverse!(info.particle_info)
    for i in eachindex(info.particle_info)
        info.id_dict[info.particle_info[i].id] = i
    end
    @test info.particle_info[1].id != 1        # the mapping really is permuted

    inter = Ewald3D(n, (L, L, L); α = 1.0, s = 4.0)
    finder = CellList3D(info, inter.short.r_c, boundary, 1)
    sys = MDSys(n_atoms = n, atoms = atoms, boundary = boundary,
                interactions = [(inter, finder)],
                loggers = [TemperatureLogger(100; output = false)],
                simulator = VerletProcess(dt = 0.001))

    poses   = [SVector(p.position[1], p.position[2], p.position[3])
               for p in info.particle_info]
    charges = [atoms[p.id].charge for p in info.particle_info]
    @test isapprox(energy(inter, finder, sys, info),
                   coulomb_energy(inter, poses, charges), rtol = 1e-10)

    # and the acceleration must land on the right particle
    F = coulomb_force(inter, poses, charges)
    for p in info.particle_info
        p.acceleration = Point(0.0, 0.0, 0.0)
    end
    ExTinyMD.update_acceleration!(inter, finder, sys, info)
    for (slot, p) in enumerate(info.particle_info)
        m = atoms[p.id].mass
        for d in 1:3
            @test isapprox(p.acceleration[d], F[slot][d] / m, rtol = 1e-10)
        end
    end
end

@testset "adapter runs inside simulate! with bounded energy drift" begin
    Random.seed!(20260932)
    n, L = 30, 12.0
    boundary, atoms, info = _charged_system(n, L)
    inter = Ewald3D(n, (L, L, L); α = 0.8, s = 3.5)   # r_c = 4.375 < L/2 = 6
    lj = LennardJones(ϵ = 1.0, σ = 1.0, cutoff = 3.0)
    finder = CellList3D(info, max(inter.short.r_c, 3.0), boundary, 1)

    sys = MDSys(n_atoms = n, atoms = atoms, boundary = boundary,
                interactions = [(lj, finder), (inter, finder)],
                loggers = [TemperatureLogger(1000; output = false)],
                simulator = VerletProcess(dt = 1e-4))

    E0 = energy(inter, finder, sys, info)
    simulate!(sys.simulator, sys, info, 200)
    E1 = energy(inter, finder, sys, info)

    @test isfinite(E1)
    # microcanonical Verlet at this dt should not let the electrostatic energy run
    # away; a sign error or a mass bug shows up as an unbounded value
    @test abs(E1 - E0) < 0.5 * max(abs(E0), 1.0)
end

@testset "adapter works for ICM with NoNeighborFinder" begin
    Random.seed!(20260933)
    n, L = 12, 8.0
    boundary, atoms, info = _charged_system(n, L)
    inter = ICMEwald2D(n, (L, L, L); α = 1.0, s = 3.0, γ = (0.3, 0.3), N_image = 3)   # r_c = 3.0 < 4
    finder = NoNeighborFinder(Float64)   # ICM keeps its own cell list
    sys = MDSys(n_atoms = n, atoms = atoms, boundary = boundary,
                interactions = [(inter, finder)],
                loggers = [TemperatureLogger(100; output = false)],
                simulator = VerletProcess(dt = 0.001))

    poses = [SVector(p.position[1], p.position[2], p.position[3])
             for p in info.particle_info]
    charges = [atoms[p.id].charge for p in info.particle_info]
    @test isapprox(energy(inter, finder, sys, info),
                   coulomb_energy(inter, poses, charges), rtol = 1e-10)
end
```

- [ ] **Step 2: Run to verify it fails**

Expected: FAIL — `MethodError: no method matching energy(::EwaldInteraction, ...)`.

- [ ] **Step 3: Write the implementation**

Create `src/interactions/electrostatics/adapter.jl`:

```julia
# Bridge between ExTinyMD's MD loop and the framework-free core API.
#
# Index convention: `info.particle_info` is indexed by storage slot, while
# `sys.atoms` is indexed by particle id. Positions are read in slot order and
# charges gathered to match, so core-layer index `i` consistently means "slot i"
# and forces come back in slot order.

"Gather charges in slot order into `buf`."
function gather_charges!(buf::Vector{T}, sys::MDSys{T}, info::SimulationInfo{T}) where {T}
    @inbounds for i in eachindex(info.particle_info)
        buf[i] = sys.atoms[info.particle_info[i].id].charge
    end
    return buf
end

"Gather positions in slot order into `buf`."
function gather_positions!(buf::Vector{SVector{3,T}},
                           info::SimulationInfo{T}) where {T}
    @inbounds for i in eachindex(info.particle_info)
        p = info.particle_info[i].position
        buf[i] = SVector{3,T}(p[1], p[2], p[3])
    end
    return buf
end

const ElectrostaticInteraction = Union{EwaldInteraction, ICM}

# Scratch for the gathered arrays, keyed on the interaction so repeated steps do
# not allocate. Stored on the interaction itself via these accessors.
_pos_scratch(inter::EwaldInteraction) = inter.pos_scratch
_charge_scratch(inter::EwaldInteraction) = inter.charge_scratch
_pos_scratch(inter::ICM) = inter.pos_scratch
_charge_scratch(inter::ICM) = inter.charge_scratch

function ExTinyMD.energy(inter::ElectrostaticInteraction, neighborfinder,
                         sys::MDSys{T}, info::SimulationInfo{T}) where {T}
    update_finder!(neighborfinder, info)
    poses   = gather_positions!(_pos_scratch(inter), info)
    charges = gather_charges!(_charge_scratch(inter), sys, info)
    return coulomb_energy(inter, poses, charges;
                          neighbor_list = _finder_list(neighborfinder))
end

function ExTinyMD.update_acceleration!(inter::ElectrostaticInteraction, neighborfinder,
                                       sys::MDSys{T}, info::SimulationInfo{T}) where {T}
    update_finder!(neighborfinder, info)
    poses   = gather_positions!(_pos_scratch(inter), info)
    charges = gather_charges!(_charge_scratch(inter), sys, info)

    F = coulomb_force!(inter.force_buffer, inter, poses, charges;
                       neighbor_list = _finder_list(neighborfinder))

    @inbounds for i in eachindex(info.particle_info)
        m = sys.atoms[info.particle_info[i].id].mass
        f = F[i]
        info.particle_info[i].acceleration += Point(f[1] / m, f[2] / m, f[3] / m)
    end
    return nothing
end

# A NoNeighborFinder carries a sentinel (0, 0, 0) entry rather than a usable list,
# so fall back to the interaction's own cell list in that case.
_finder_list(f::NoNeighborFinder) = nothing
_finder_list(f) = f.neighbor_list
```

This requires `pos_scratch` and `charge_scratch` fields. Add them to `EwaldInteraction` in `ewald.jl`:

```julia
struct EwaldInteraction{T, S, L} <: AbstractInteraction
    short::S
    long::L
    n_atoms::Int
    force_buffer::Vector{SVector{3,T}}
    pos_scratch::Vector{SVector{3,T}}
    charge_scratch::Vector{T}
end

function EwaldInteraction(short::S, long::L, n_atoms::Int) where {S, L}
    T = typeof(short.α)
    return EwaldInteraction{T, S, L}(short, long, n_atoms,
                                     [zero(SVector{3,T}) for _ in 1:n_atoms],
                                     [zero(SVector{3,T}) for _ in 1:n_atoms],
                                     zeros(T, n_atoms))
end
```

and to `ICM` in `icm.jl`, with `force_buffer`, `pos_scratch` and `charge_scratch` of length
`n_atoms` alongside the existing `ref_*` buffers of length `n_ref`:

```julia
struct ICM{T, L, S} <: AbstractInteraction
    long::L
    short::S
    γ::Tuple{T,T}
    N_image::Int
    n_atoms::Int
    L::NTuple{3,T}
    elc::Bool
    N_pad::Int
    ref_poses::Vector{SVector{3,T}}
    ref_charges::Vector{T}
    ref_force::Vector{SVector{3,T}}
    force_buffer::Vector{SVector{3,T}}
    pos_scratch::Vector{SVector{3,T}}
    charge_scratch::Vector{T}
end
```

Update both constructors to allocate the new fields, and re-run Tasks 6–8's tests to confirm
the struct change broke nothing.

Add to `src/ExTinyMD.jl` (last of the electrostatics includes):

```julia
include("interactions/electrostatics/adapter.jl")
```
```julia
export gather_charges!, gather_positions!
```

- [ ] **Step 4: Run tests to verify they pass**

Run the whole suite — the struct field additions touch Tasks 6–8.

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/ExTinyMD.jl src/interactions/electrostatics/ test/
git commit -m "feat: add MD adapter for electrostatic interactions

ExTinyMD.energy and update_acceleration! over the framework-free core, reusing
the simulation's own neighbour list and preallocated gather buffers."
```

---

### Task 10: Documentation

**Files:**
- Create: `docs/Project.toml`, `docs/make.jl`, `docs/src/index.md`, `docs/src/md_core.md`, `docs/src/interactions.md`, `docs/src/electrostatics.md`, `docs/src/api.md`
- Modify: `.github/workflows/CI.yml`, `README.md`, `.gitignore`

**Interfaces:**
- Consumes: every public name from Tasks 3–9.
- Produces: a Documenter site that builds with `--strict` and runs doctests.

- [ ] **Step 1: Create the docs environment**

`docs/Project.toml`:

```toml
[deps]
Documenter = "e30172f5-a6a5-5a46-863b-614d45cd2de4"
ExTinyMD = "fec76197-d59f-46dd-a0ed-76a83c21f7aa"
StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[compat]
Documenter = "1"
```

`docs/make.jl`:

```julia
using Documenter, ExTinyMD

DocMeta.setdocmeta!(ExTinyMD, :DocTestSetup, :(using ExTinyMD, StaticArrays);
                    recursive = true)

makedocs(
    sitename = "ExTinyMD.jl",
    modules  = [ExTinyMD],
    format   = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = [
        "Home"           => "index.md",
        "MD Core"        => "md_core.md",
        "Interactions"   => "interactions.md",
        "Electrostatics" => "electrostatics.md",
        "API Reference"  => "api.md",
    ],
    checkdocs = :exports,
    warnonly  = false,
)

deploydocs(repo = "github.com/HPMolSim/ExTinyMD.jl", devbranch = "main")
```

Add to `.gitignore`:

```
docs/build/
docs/Manifest.toml
```

- [ ] **Step 2: Write `docs/src/electrostatics.md`**

This is the substantive page. It must contain, as real prose and real runnable code:

1. **The two-layer API.** A worked standalone example and the same system inside `simulate!`.

````markdown
# Electrostatics

Every method here has two interfaces. The **core API** takes a plan object plus
plain arrays and requires no ExTinyMD type at all:

```julia
using ExTinyMD, StaticArrays

n, L = 100, (20.0, 20.0, 20.0)
poses   = [SVector(rand()*L[1], rand()*L[2], rand()*L[3]) for _ in 1:n]
charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]

inter = Ewald3D(n, L; α = 0.5, s = 4.0)   # r_c = s/α = 8.0, below L/2 = 10
E = coulomb_energy(inter, poses, charges)
F = coulomb_force(inter, poses, charges)
```

The page must state the `r_c < min(L)/2` rule explicitly at this point, because it is the
first thing a user hits: `α` and `s` are not independent of the box, and picking them
carelessly raises a `UNIT CELL CHECK FAILED` error from CellListMap rather than returning a
slightly wrong number.

The **MD interface** is the same object added to an `MDSys`, where
`update_acceleration!` is called for you each step.
````

2. **Method selection table.**

| Boundary condition | Dielectric walls | Method | Cost | Notes |
|---|---|---|---|---|
| Triply periodic | none | `Ewald3D` | `O(N·K)` | reference implementation |
| Triply periodic | none | `PME3D` | `O(N log N)` | Phase 2; needs `using FINUFFT` |
| Slab, periodic in x,y | none | `Ewald2D` | `O(N²K)` | exact; accuracy reference only |
| Slab, periodic in x,y | confined | `ICMEwald2D` | `O(N²K)` | exact for confined slabs |
| Slab, periodic in x,y | confined | `ICMEwald3D` | `O(N·K)` | Ewald3D + ELC; faster, approximate in `N_pad` |

3. **Parameter selection.** `r_c = s/α`, `k_c = 2αs`; `s` controls accuracy (error roughly `exp(-s²)`); `α` trades real- against reciprocal-space work and the total energy must be independent of it — a good self-check. State the α-independence test as the recommended way for a user to validate their own parameter choice.

4. **The prefactor asymmetry** of spec §5.4, stated plainly so a reader of the source is not misled.

5. **Error estimates.** Port the three functions from `EwaldSummations/src/error_estimate.jl` into the docs as guidance with the arXiv:2503.18126 citation. Do **not** add them to `src/` in this phase — they are not needed by any code path here.

- [ ] **Step 3: Write the remaining pages**

- `index.md`: what ExTinyMD is, installation, the LJ-fluid example from `example/LJ_fluid_rdf.jl` as a runnable walkthrough, and a pointer to the electrostatics page.
- `md_core.md`: `Boundary` / `Q2dBoundary` / `CubicBoundary`, `Atom`, `MDSys`, `SimulationInfo`, the four simulators, the three thermostats, the three loggers, and the neighbour finders — with `@docs` blocks.
- `interactions.md`: `LennardJones`, `SubLennardJones`, `ExternalField`, and how to write a new interaction (the `energy` / `update_acceleration!` contract).
- `api.md`: `@autodocs` over `ExTinyMD` to catch everything, so `checkdocs = :exports` has something to check against.

- [ ] **Step 4: Add docstrings to pre-existing public API**

`checkdocs = :exports` fails the build for any exported name without a docstring. ExTinyMD currently has essentially none, so this step is required for the build to pass. Add a docstring to every name in the `export` lists of `src/ExTinyMD.jl`. Keep them short — one sentence plus arguments for simple types.

- [ ] **Step 5: Build the docs locally**

```
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=".")); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

Expected: builds with no warnings, doctests pass. Fix any `checkdocs` failure by writing the missing docstring, never by relaxing `checkdocs`.

- [ ] **Step 6: Add the CI docs job**

Append to `.github/workflows/CI.yml`:

```yaml
  docs:
    name: Documentation
    runs-on: ubuntu-latest
    permissions:
      contents: write
      statuses: write
    steps:
      - uses: actions/checkout@v6
      - uses: julia-actions/setup-julia@v2
        with:
          version: '1'
      - uses: julia-actions/cache@v2
      - name: Install dependencies
        run: julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
      - name: Build and deploy
        env:
          GITHUB_TOKEN: ${{ secrets.GITHUB_TOKEN }}
          DOCUMENTER_KEY: ${{ secrets.DOCUMENTER_KEY }}
        run: julia --project=docs docs/make.jl
```

Check the existing `CI.yml` for the `actions/checkout` version already in use and match it rather than assuming v6.

- [ ] **Step 7: Update the README**

Add a documentation badge and a short electrostatics section pointing at the new page, listing the four available methods.

- [ ] **Step 8: Commit**

```bash
git add docs/ .github/workflows/CI.yml README.md .gitignore src/
git commit -m "docs: add Documenter site with electrostatics guide

Method selection table, parameter guidance, the two-layer API, and docstrings
for the previously undocumented public API."
```

---

## Self-Review

**Spec coverage.**

| Spec section | Task |
|---|---|
| §1 `neighbor_list` bug | 1 |
| §4.1 core layer | 4, 5, 6, 7, 8 |
| §4.2 adapter, id/slot convention | 9 |
| §5.1 common.jl | 3 |
| §5.2 EwaldShort, analytic force, task-partitioned threading | 4 |
| §5.3 Ewald3DLong, `Complex{T}` | 5 |
| §5.4 Ewald2DLong, overflow guard, threading restructure | 7 |
| §5.5 ICM, ELC, force convention | 8 |
| §5.6 composition, constructors | 6, 7, 8 |
| §5.7 PME3D | **Phase 2 — deliberately not in this plan** |
| §6 data flow | 9 |
| §7 ICM risk | 8 Step 6 |
| §8.1 unit | 2, 3, 4 |
| §8.2 cross-validation | 5, 7, 8 |
| §8.3 finite-difference forces | 4, 5, 6, 7, 8 |
| §8.4 integration | 9 |
| §8.5 regression | 1 |
| §9 docs | 10 |

Gaps accepted for this phase: §5.7 (PME3D) is Phase 2 by design. §8.2's
`test/validation/compare_ewaldsummations.jl` migration script is deferred to Phase 3, where
it gates the deletion of EwaldSummations' k-space code — writing it now would test code
nothing yet depends on.

**Type consistency.** `coulomb_energy` / `coulomb_force` / `coulomb_force!` are the core
queries throughout (Tasks 6, 8, 9); `short_energy` / `short_force!` and `long_energy` /
`long_force!` are the component queries (Tasks 4, 5, 7). `min_image_disp` is spelled
consistently from Task 3 onward. Accumulation contract: component `*_force!` functions
accumulate, composite `coulomb_force!` zeroes first — asserted by a test in Task 6. The
struct field additions in Task 9 modify structs defined in Tasks 6 and 8, which Task 9
Step 4 explicitly re-tests.

**Naming decision recorded:** the core queries are `coulomb_*`, not `energy` / `force`,
because `ExTinyMD.energy` already exists with a four-argument MD signature. This diverges
from the spec's §4.1 sketch, which wrote `energy(plan, poses, charges)`; the rename avoids
an overload that differs only in arity on an exported name.
