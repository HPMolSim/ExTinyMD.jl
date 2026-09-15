"""
    PME3DLong(n_atoms, L; α, s, ϵ = 1.0, ϵ_inf = Inf)

Particle-mesh reciprocal-space solver for the triply periodic Ewald sum, computing
the same quantity as [`Ewald3DLong`](@ref) in `O(N log N)` instead of `O(N·K)`.

The structure factor is evaluated with a type-1 NUFFT onto a rectangular grid, and
the spherical cutoff `0 < |k| ≤ k_c` is recovered by zeroing the Green's function
`D_k` outside it. The two solvers therefore sum over an **identical** k-set and
agree to machine precision, rather than merely converging to the same limit.

Requires `using FINUFFT`; the implementation lives in a package extension so the
binary artifact stays off users who never ask for it.
"""
struct PME3DLong{T, P1, P2}
    α::T
    r_c::T
    k_c::T
    ϵ::T
    ϵ_inf::T
    L::NTuple{3,T}
    n_atoms::Int
    n_k::NTuple{3,Int}
    # D is the spherical-masked Green's function: exp(-k²/4α²)/k² inside the
    # cutoff, exactly zero outside it and at k = 0.
    D::Array{T,3}
    ρ_src::Array{Complex{T},3}
    ρ_tgt::Array{Complex{T},3}
    hx::Array{Complex{T},3}
    hy::Array{Complex{T},3}
    hz::Array{Complex{T},3}
    xs::Vector{T}
    ys::Vector{T}
    zs::Vector{T}
    ox::Vector{Complex{T}}
    oy::Vector{Complex{T}}
    oz::Vector{Complex{T}}
    qs::Vector{Complex{T}}
    plan1::P1
    plan2::P2
end

Base.show(io::IO, l::PME3DLong) =
    print(io, "PME3DLong(α = $(l.α), k_c = $(l.k_c), ϵ = $(l.ϵ), ϵ_inf = $(l.ϵ_inf), " *
              "grid = $(size(l.D)))")

"""
    PME3D(n_atoms, L; α, s, ϵ = 1.0, ϵ_inf = Inf)

Particle-mesh Ewald summation for a triply periodic system: [`EwaldShort`](@ref)
paired with [`PME3DLong`](@ref). Computes the same energy and forces as
[`Ewald3D`](@ref) at `O(N log N)` instead of `O(N·K)`.

Requires `using FINUFFT`.
"""
function PME3D end

"""
    ICMPME3D(n_atoms, L; α, s, γ, N_image, N_pad, ϵ = 1.0)

The image-charge method combined with particle-mesh Ewald and the electrostatic
layer correction — the `O(N log N)` counterpart of [`ICMEwald3D`](@ref).

Requires `using FINUFFT`.
"""
function ICMPME3D end

# Fallbacks so a user who forgot the extension gets a sentence rather than a
# MethodError listing zero candidates.
const _FINUFFT_HINT = "requires FINUFFT. Run `using FINUFFT` (and add it to your " *
                      "project) to load ExTinyMD's particle-mesh extension."

PME3D(args...; kwargs...)    = error("PME3D " * _FINUFFT_HINT)
ICMPME3D(args...; kwargs...) = error("ICMPME3D " * _FINUFFT_HINT)
PME3DLong(args...; kwargs...) = error("PME3DLong " * _FINUFFT_HINT)
