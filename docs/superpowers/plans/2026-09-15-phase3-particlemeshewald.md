# Phase 3a: ParticleMeshEwald — Decoupling Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `ParticleMeshEwald.jl` co-resolve with current ExTinyMD, give it an array-of-structs query API, and write the ExTinyMD extension its `[weakdeps]` has always implied but never had.

**Architecture:** The two-layer contract Phases 1–2 established: a framework-free core taking a plan object plus AoS positions and charges, and an ExTinyMD adapter in `ext/`. ExTinyMD moves from a declared-but-inert `[weakdeps]` entry to a working one.

**Tech Stack:** Julia 1.10+, CellListMap 0.10 (bumped from 0.9), FINUFFT, KernelAbstractions, SpecialFunctions.

**Spec:** `docs/superpowers/specs/2026-09-15-downstream-decoupling-design.md` (this is §5.1, the first of five packages)

**Repository:** `/mnt/home/xgao1/project/q2dmd/ParticleMeshEwald.jl` — **a different repository from the one this plan lives in.** ExTinyMD is at `../ExTinyMD.jl`.

## Baseline, measured

- ParticleMeshEwald's own suite **passes today**: 9 tests, on CellListMap 0.9.
- It **cannot co-resolve with ExTinyMD 0.3.** Verified:

```
Unsatisfiable requirements detected for package CellListMap [69e1c6dd]:
 ├─restricted to versions 0.9.14 - 0.9 by ParticleMeshEwald, leaving only 0.9.14 - 0.9.17
 └─restricted to versions 0.10 by ExTinyMD — no versions left
```

That conflict, not broken code, is the blocker. Everything else in this plan is the decoupling the user asked for plus four defects found along the way.

## Global Constraints

- `julia = "1.10"` minimum. No syntax newer than 1.10.
- **CellListMap 0.10 uses keyword arguments**: `InPlaceNeighborList(xpositions = ...)`, `update!(cl, xpositions = ...)`, `neighborlist(xpositions = ..., cutoff = ..., unitcell = ...)`. The 0.9 spellings (`x = ...`, positional `update!`) raise `MethodError`. ExTinyMD's `src/MD_core/neighbor_finder/cell_list.jl` is the reference for the 0.10 idiom.
- **`r_c` must be strictly less than half the smallest periodic box side.** CellListMap 0.10 enforces this with `ArgumentError: UNIT CELL CHECK FAILED ... must be greater than 2*cutoff`. Since `r_c = s/α`, every α/s choice is constrained by the box. Check every test configuration; this caught every parameter set in an earlier phase's first draft.
- **ExTinyMD goes in `[weakdeps]` and `[extensions]`, never `[deps]`.** It also goes in `[extras]` and the `test` target so `Pkg.test()` exercises the extension — that is the one allowed exception and does not make it a runtime dependency.
- **`energy`, `energy_short`, `energy_long` must NOT be exported** after this work. Five packages each exporting `energy` would make `using EwaldSummations, ParticleMeshEwald` ambiguous. Callers write `ParticleMeshEwald.energy(...)`. The `PME` constructor stays exported. This is a deliberate breaking change to a 0.1.0 package; note it in the README.
- **Never scale or otherwise mutate the caller's position arrays.** Scale into plan-owned buffers.
- **Never index a thread-local accumulator by `Threads.threadid()`** — unsound under task migration.
- **No forces.** ParticleMeshEwald has never had them, and ExTinyMD's `PME3D` (Phase 2) provides PME with forces already. Adding them here would duplicate that. The consequence is that this package's extension can supply `ExTinyMD.energy` but not `update_acceleration!`, so it cannot drive an MD run — document that and point users at `PME3D`.
- Work on a branch in the ParticleMeshEwald repository. Commit after every task.
- Suite: `julia --project=. -e 'using Pkg; Pkg.test()'` from that package's root.

## File Structure

| File | Responsibility |
|---|---|
| `src/types.jl` | `PME` plan struct; gains AoS scratch buffers, loses the SoA `pos` matrix if unused |
| `src/energy.jl` | AoS `energy`/`energy_short`/`energy_long`; no caller mutation; sound threading |
| `src/ParticleMeshEwald.jl` | module: exports, and removal of the KernelAbstractions examples include |
| `ext/ParticleMeshEwaldExTinyMDExt.jl` | **new** — `ExTinyMD.energy` adapter |
| `Project.toml` | CellListMap 0.10, `[extensions]`, ExTinyMD compat, stdlib compat fix |
| `test/energy.jl`, `test/runtests.jl` | AoS call sites, plus the new standalone and adapter tests |
| `README.md` | the AoS API, the un-exported names, the no-forces note |

---

### Task 1: Bump CellListMap to 0.10 and fix its API calls

The prerequisite for everything: without this the package cannot sit in the same environment as ExTinyMD.

**Files:** `Project.toml`, `src/types.jl`, `src/energy.jl`

- [ ] **Step 1: Record the baseline numbers**

Before changing anything, run the suite and **write down the energies its tests assert**, so you can prove the bump changed no physics:

```
julia --project=. -e 'using Pkg; Pkg.test()'
```

Put those values in your report. If any test only asserts a comparison rather than a number, capture the computed value by other means — you need something to compare against afterwards.

- [ ] **Step 2: Bump the compat entry**

In `Project.toml`, `CellListMap = "0.9.14"` becomes `CellListMap = "0.10"`.

While you are there, **fix two other compat problems**:

- `LinearAlgebra = "1.12.0"` pins a *stdlib* to a version that only exists on Julia 1.12+, contradicting `julia = "1.10, 1.11, 1.12"`. Stdlib compat entries should be loose (`"1"`) or absent. Same check for any other stdlib listed.
- `julia = "1.10, 1.11, 1.12"` excludes 1.13, which is current. Widen it to `"1.10"` (meaning 1.10 and up), matching ExTinyMD.

- [ ] **Step 3: Run it and watch it fail**

```
julia --project=. -e 'using Pkg; Pkg.resolve(); Pkg.test()'
```

Expected: `MethodError` on `InPlaceNeighborList` — 0.10 renamed the `x` keyword to `xpositions`, and `update!` no longer takes positions positionally.

- [ ] **Step 4: Fix the call sites**

`InPlaceNeighborList(x = pos, ...)` becomes `InPlaceNeighborList(xpositions = pos, ...)`; `update!(cl, pos)` becomes `update!(cl, xpositions = pos)`. Grep for every `InPlaceNeighborList`, `update!` and `neighborlist` in `src/`.

- [ ] **Step 5: Run and confirm the numbers are unchanged**

The suite must pass **and** produce the same energies you recorded in Step 1. Report both sets. A bump that changes a number is a finding, not a success — if one moves, stop and report it rather than updating the expected value.

- [ ] **Step 6: Confirm co-resolution now works**

```bash
julia --startup-file=no --project=$(mktemp -d) -e '
using Pkg
Pkg.develop(path="/mnt/home/xgao1/project/q2dmd/ExTinyMD.jl")
Pkg.develop(path="/mnt/home/xgao1/project/q2dmd/ParticleMeshEwald.jl")
Pkg.status("CellListMap")
println("co-resolved")'
```

Expected: succeeds, with one CellListMap 0.10.x. Paste the output — this is the whole point of the task.

- [ ] **Step 7: Commit**

```bash
git commit -m "compat: move to CellListMap 0.10 so the package can co-resolve with ExTinyMD

ParticleMeshEwald pinned CellListMap 0.9 while ExTinyMD 0.3 requires 0.10, so
the two could not sit in one environment. Also loosens a LinearAlgebra stdlib
compat that pinned 1.12.0 while claiming Julia 1.10 support, and widens the
julia compat to 1.10 and up."
```

---

### Task 2: Array-of-structs query API

**Files:** `src/energy.jl`, `src/types.jl`, `test/energy.jl`

Today the queries are structure-of-arrays and mutate their inputs:

```julia
energy(pme, x, y, z, q)      # three coordinate vectors plus complex charges
```

`energy_long` does `x .*= 2π/L[1]` and divides back afterwards, so it corrupts the caller's data if the transform throws in between.

**Interfaces to produce:**

- `ParticleMeshEwald.energy(pme, poses, charges)` — `poses` AoS, `charges::AbstractVector{<:Real}`
- `ParticleMeshEwald.energy_short(pme, poses, charges; neighbor_list = nothing)`
- `ParticleMeshEwald.energy_long(pme, poses, charges)`

`poses` elements are accessed only via `p[1]`, `p[2]`, `p[3]`, so `SVector{3,T}`, `NTuple{3,T}` and ExTinyMD's `Point{3,T}` all work with no conversion layer.

Note the charges become plain reals in the public API — the `Complex{T}` conversion FINUFFT needs is an internal detail and belongs in a plan-owned buffer, not in the caller's hands. Today callers must pass `ComplexF64.(q)`, which is an implementation detail leaking into the interface.

- [ ] **Step 1: Write the failing tests**

Add to `test/energy.jl` (adapt names to what is there):

```julia
@testset "AoS query API" begin
    Random.seed!(2026)
    n = 100
    L = (20.0, 20.0, 20.0)
    α, s = 0.5, 4.0            # r_c = s/α = 8.0 < min(L)/2 = 10
    poses = [SVector(rand() * L[1], rand() * L[2], rand() * L[3]) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]

    pme = PME(α, L, s, n)
    E = ParticleMeshEwald.energy(pme, poses, charges)
    @test isfinite(E)

    # the same positions as NTuple and as a 3-column read must agree
    @test ParticleMeshEwald.energy(pme, [Tuple(p) for p in poses], charges) ≈ E
end

@testset "energy does not mutate the caller's positions" begin
    Random.seed!(7)
    n = 50
    L = (20.0, 20.0, 20.0)
    poses = [SVector(rand() * L[1], rand() * L[2], rand() * L[3]) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]
    before = deepcopy(poses)

    pme = PME(0.5, L, 4.0, n)      # r_c = 8.0 < min(L)/2 = 10
    ParticleMeshEwald.energy(pme, poses, charges)
    @test poses == before
end
```

Note the inline `r_c` annotation on the `α, s` line. Every parameter set in an earlier phase's first draft violated the `r_c < min(L)/2` bound, and a documentation example I wrote myself shipped with a violation that raised the very error the page was warning about. Annotate every configuration you write the same way — it makes the check mechanical instead of remembered.

- [ ] **Step 2: Run and watch them fail** — `MethodError`, since `energy` takes five positional arguments today.

- [ ] **Step 3: Convert the queries**

Add plan-owned scaled-coordinate and charge buffers to `PME` in `src/types.jl` — `xs`, `ys`, `zs` of length `N` and `qs::Vector{Complex{T}}` — and have the queries scatter into them rather than touching `poses` or building new arrays. Keep the SoA methods only if something still needs them; if nothing does, delete them rather than leaving two paths.

- [ ] **Step 4: Run, confirm the recorded baseline energies still match, commit.**

The physics must not move. Report the before/after numbers.

---

### Task 3: Fix the unsound threading and the fragile include

**Files:** `src/energy.jl`, `src/ParticleMeshEwald.jl`

Two defects, both found by inspection during Phase 1's survey.

**3a. `energy_short_kernel!` indexes a thread-local accumulator by `Threads.threadid() - 1`.** When `threadid()` is 1 that is index 0, which is out of bounds in Julia; and `threadid()` is not a safe key for accumulation under a migrating task scheduler regardless. There is a separate single-threaded code path, which is presumably why the bug has not surfaced.

Replace with a task-partitioned reduction: split the neighbour list into chunks, give each chunk its own accumulator, reduce at the end. ExTinyMD's `src/interactions/electrostatics/short.jl` does the same job serially and is a good model for the arithmetic; the partitioning is yours.

**Then decide, and say which you chose:** whether the KernelAbstractions path earns its keep at all. It exists to allow a GPU backend, but the package has no GPU tests, `CPU()` is the only backend ever passed, and the kernel is the source of this bug. If removing it simplifies the package without losing a capability anyone uses, remove it and drop the KernelAbstractions dependency. If you keep it, fix the indexing properly. **This is a judgement call I am delegating — argue for what you chose.**

**3b. `src/ParticleMeshEwald.jl:7` reaches into a dependency's examples directory:**

```julia
include(joinpath(dirname(pathof(KernelAbstractions)), "../examples/utils.jl"))
```

That file is not part of KernelAbstractions' public API and can vanish in any patch release. Whatever it provides must be vendored into this package or removed along with the kernel. Find out what it is actually used for before deciding — it may well be unused.

- [ ] Write a test that exercises the multi-threaded short-range path (run under `-t 4` if needed) and asserts it agrees with the single-threaded one. Confirm it fails against the `threadid() - 1` version — if it cannot be made to fail, say so, because then the bug is unreachable and removing the kernel is the better answer.
- [ ] Fix, run, commit.

---

### Task 4: The ExTinyMD extension

**Files:** `ext/ParticleMeshEwaldExTinyMDExt.jl` (new), `Project.toml`, `test/`

`[weakdeps]` already lists ExTinyMD (and EwaldSummations, and Random), but there is no `[extensions]` entry and no `ext/` directory, so Julia has never loaded anything for them. They are inert declarations.

- [ ] **Step 1: Write the failing test**

```julia
@testset "ExTinyMD adapter" begin
    using ExTinyMD
    n, L = 60, 20.0
    boundary = Boundary((L, L, L), (1, 1, 1))
    atoms = Atom{Float64}[]
    for _ in 1:(n ÷ 2);      push!(atoms, Atom(type = 1, mass = 1.0, charge =  1.0)); end
    for _ in (n ÷ 2 + 1):n;  push!(atoms, Atom(type = 2, mass = 1.0, charge = -1.0)); end
    info = SimulationInfo(n, atoms, (0.0, L, 0.0, L, 0.0, L), boundary;
                          min_r = 1.0, temp = 1.0)
    info.running_step = 1

    pme = PME(0.5, (L, L, L), 4.0, n)        # r_c = 8.0 < 10
    finder = CellList3D(info, pme.r_c, boundary, 1)
    sys = MDSys(n_atoms = n, atoms = atoms, boundary = boundary,
                interactions = [(pme, finder)],
                loggers = [TemperatureLogger(100; output = false)],
                simulator = VerletProcess(dt = 1e-4))

    poses   = [SVector(p.position[1], p.position[2], p.position[3]) for p in info.particle_info]
    charges = [atoms[p.id].charge for p in info.particle_info]

    @test ExTinyMD.energy(pme, finder, sys, info) ≈
          ParticleMeshEwald.energy(pme, poses, charges)
end
```

`PME` will need to subtype `ExTinyMD.AbstractInteraction` for `MDSys` to accept it — but that type only exists when ExTinyMD is loaded, and `PME` is defined in `src/`. **Resolve this and explain your choice in the report.** The two options are: declare a local abstract type in `src/` and have the extension do nothing about the hierarchy (in which case `MDSys` will not accept it, and the test must be restructured to call `ExTinyMD.energy` directly without an `MDSys`), or restrict the adapter to the query functions and accept that `PME` cannot be placed in `sys.interactions`. Phase 1's ExTinyMD-internal types did not face this because they live inside ExTinyMD.

This is the most interesting design question in the task and the reason ParticleMeshEwald goes first: whatever you decide here sets the pattern for the four remaining packages. Take the time to get it right and write down why.

- [ ] **Step 2:** Add `[extensions] ParticleMeshEwaldExTinyMDExt = "ExTinyMD"`, ExTinyMD to `[extras]` and the `test` target, and `ExTinyMD = "0.3"` to `[compat]`. Remove the inert `EwaldSummations` and `Random` weakdeps unless you give them extensions too — a `[weakdeps]` entry with no `[extensions]` entry is dead weight that misleads the reader.

- [ ] **Step 3:** Write the adapter. It extracts charges and positions honouring ExTinyMD's id/slot indirection (`sys.atoms` is indexed by particle **id**, `info.particle_info` by storage **slot**), then calls the core. `../ExTinyMD.jl/src/interactions/electrostatics/adapter.jl` is the reference.

Provide `ExTinyMD.energy` only. **Do not** provide `update_acceleration!` — this package has no forces, and inventing them here would duplicate ExTinyMD's `PME3D`. Add a docstring saying so and pointing at `PME3D`.

- [ ] **Step 4:** Run, commit.

---

### Task 5: The standalone test, and documentation

**Files:** `test/`, `README.md`

- [ ] **Step 1: The test this whole phase exists for**

A test proving the core works with **ExTinyMD never loaded**. It must run in a session where ExTinyMD was not imported — a `@testset` inside the normal suite does not prove that, because `test/runtests.jl` loads ExTinyMD for the adapter test.

Run it as a subprocess with a clean environment:

```julia
@testset "core works without ExTinyMD" begin
    script = """
    using ParticleMeshEwald, StaticArrays
    @assert !haskey(Base.loaded_modules, Base.PkgId(
        Base.UUID("fec76197-d59f-46dd-a0ed-76a83c21f7aa"), "ExTinyMD"))
    n = 50; L = (20.0, 20.0, 20.0)
    poses = [SVector(rand()*L[1], rand()*L[2], rand()*L[3]) for _ in 1:n]
    charges = [isodd(i) ? 1.0 : -1.0 for i in 1:n]
    pme = PME(0.5, L, 4.0, n)
    E = ParticleMeshEwald.energy(pme, poses, charges)
    @assert isfinite(E)
    print("OK")
    """
    out = read(`\$(Base.julia_cmd()) --startup-file=no --project=\$(Base.active_project()) -e \$script`, String)
    @test out == "OK"
end
```

The `haskey(Base.loaded_modules, ...)` assertion is the load-bearing part: without it the test would pass even if ExTinyMD were loaded, and it would then be proving nothing. **Verify it fails if you add `using ExTinyMD` to the script** — otherwise you have written another test that passes for the wrong reason, which has happened ten times across these phases.

- [ ] **Step 2: README**

Document the AoS API, that `energy`/`energy_short`/`energy_long` are no longer exported and why, the `r_c < min(L)/2` rule, and that this package provides energy only — with a pointer to ExTinyMD's `PME3D` for anyone who needs forces or MD integration.

- [ ] **Step 3: Run the full suite, commit.**

---

## Self-Review

**Spec coverage.** §4.1 core layer → Task 2. §4.2 naming → Global Constraints and Task 5. §4.3 adapter → Task 4. §4.4 `Project.toml` → Tasks 1 and 4. §5.1's two named defects → Task 3. §6 testing: standalone smoke test → Task 5 Step 1; adapter test → Task 4; before/after numerical check → Tasks 1 and 2, where recording the baseline is an explicit step.

**What this plan does not do.** No forces (Global Constraints), so §6's `simulate!` adapter test is not applicable to this package — `ExTinyMD.energy` is testable but `update_acceleration!` does not exist to test. That is a real gap in the phase's success criterion 4 for *this* package only, and it is the direct consequence of the spec's §8 question about whether ParticleMeshEwald should continue to exist. Flag it rather than paper over it.

**Open questions delegated to the implementer, each requiring written justification:** whether to keep the KernelAbstractions kernel at all (Task 3a), what to do about the `examples/utils.jl` include (Task 3b), and how `PME` should relate to `ExTinyMD.AbstractInteraction` given it is defined in `src/` (Task 4 Step 1). The third sets the pattern for the four remaining packages.
