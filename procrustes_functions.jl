using Manopt, Manifolds, LinearAlgebra, RecursiveArrayTools
# Orthogonal Procrustes -> Frobenius norm -> closed form solution

# The following lines temporarily fix a Manopt bug 
import Manopt: mesh_adaptive_direct_search
function mesh_adaptive_direct_search(M::AbstractManifold, f, p=rand(M); kwargs...)
    mco = ManifoldCostObjective(f)
    return mesh_adaptive_direct_search(M, mco, p; kwargs...)
end

# The following lines temporarily fix a Manifolds bug 
import Manifolds: angles_4d_skew_sym_matrix
function angles_4d_skew_sym_matrix(A)
    @assert size(A) == (4, 4)
    @inbounds begin
        halfb = (A[1, 2]^2 + A[1, 3]^2 + A[2, 3]^2 + A[1, 4]^2 + A[2, 4]^2 + A[3, 4]^2) / 2
        c = (A[1, 2] * A[3, 4] - A[1, 3] * A[2, 4] + A[1, 4] * A[2, 3])^2
    end
    sqrtdisc = sqrt(max(halfb^2 - c,0))
    return sqrt(halfb + sqrtdisc), sqrt(halfb - sqrtdisc)
end

@doc raw"""
	Spectral_cost(M, p, A, B)

Compute the spectral norm ``\lVert p\cdot B- A\rVert_2``.
"""
spectral_procrustes_cost(M, p, A, B) = opnorm(p*B - A)

# ╔═╡ e7fae1c9-4800-4db4-bc66-24d0b89f24da
@doc raw"""
	Frobenius_cost(M, p)

Compute the Frobenious norm ``\lVert p\cdot B- A\rVert_{\mathrm{F}}``.
"""
orthogonal_procrustes_cost(M, p, A, B) = norm(p*B - A)

# ╔═╡ f5f795c2-a25f-48b3-86b8-ee43827eb719
@doc raw"""
	TwoOne_cost(M, p)

Compute the 2-1-norm ``\lVert p\cdot B- A\rVert_{2,1}``, i.e. the sum of the 2-norms of the columns.
"""
robust_procrustes_cost(M, p, A, B) = sum([norm(c, 2.0) for c in eachcol(p*B-A)])


"""
    orthogonal_procrustes(A,B[, p0]; kwargs...)

Compute the orthogonal procrustes, i.e. the minimizer of [`orthogonal_procrustes_cost`](@ref) which has a closed form solution.

The `kwargs...` and `p0` are just for compatibility to the other two calls and have no effect.
"""
function orthogonal_procrustes(A,B, p0=nothing; kwargs...)
  svd_AB = svd(A*B')
	return svd_AB.U*svd_AB.V'
end

"""
    spectral_procrustes(A,B[, p0]; kwargs...)

Compute the 2->1 norm (or robust) procrustes minimizer, see [`spectral_procrustes_cost`](@ref).
All kyword arguments are passed to the MADS solver.

You can provide an initial guess `p0` for the solver to start at, but default the identity is used.
"""
function spectral_procrustes(A,B, p0=nothing; kwargs...)
    n,m = size(A)
    M = Rotations(n)
    p0 = isnothing(p0) ? Matrix{Float64}(I,n,n) : p0
    p = mesh_adaptive_direct_search(M, (M,p) -> spectral_procrustes_cost(M, p, A, B), p0; kwargs...)
    return p
end

"""
    robust_procrustes(A,B[, p0]; kwargs...)

Compute the 2->1 norm (or robust) procrustes minimizer, see [`robust_procrustes_cost`](@ref).
All kyword arguments are passed to the MADS solver.

You can provide an initial guess `p0` for the solver to start at, but default the identity is used.
"""
function robust_procrustes(A, B, p0=nothing; kwargs...)
    n,m = size(A)
    M = Rotations(n)
    p0 = isnothing(p0) ? Matrix{Float64}(I,n,n) : p0
    p = mesh_adaptive_direct_search(M, (M,p) -> robust_procrustes_cost(M, p, A, B), p0; kwargs...)
    return p
 end
