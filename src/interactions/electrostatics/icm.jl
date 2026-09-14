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
    # reflected configuration: length n_ref = n_atoms * (1 + 2*N_image)
    ref_poses::Vector{SVector{3,T}}
    ref_charges::Vector{T}
    ref_force::Vector{SVector{3,T}}
    # real particles only: length n_atoms. force_buffer is what the adapter hands
    # to coulomb_force!; pos_scratch/charge_scratch are its gather buffers.
    force_buffer::Vector{SVector{3,T}}
    pos_scratch::Vector{SVector{3,T}}
    charge_scratch::Vector{T}
end

function ICM(long::L, short::S, γ::Tuple{T,T}, N_image::Int, n_atoms::Int,
             Lbox::NTuple{3,T}; elc::Bool = false, N_pad::Int = 0) where {T, L, S}
    n_ref = n_atoms * (1 + 2 * N_image)
    return ICM{T, L, S}(long, short, γ, N_image, n_atoms, Lbox, elc, N_pad,
                        [zero(SVector{3,T}) for _ in 1:n_ref], zeros(T, n_ref),
                        [zero(SVector{3,T}) for _ in 1:n_ref],
                        [zero(SVector{3,T}) for _ in 1:n_atoms],
                        [zero(SVector{3,T}) for _ in 1:n_atoms],
                        zeros(T, n_atoms))
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
