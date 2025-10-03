# SpectralProcrustes
Codebase for "Procrustes Problems on Random Matrices" by Hajg Jasa, Ronny Bergmann, Christian Kümmerle, Avanti Athreya, and Zachary Lubberts.

**Abstract:** Meaningful comparison between sets of observations often necessitates alignment or registration between them, and the resulting optimization problems range in complexity from those admitting simple closed-form solutions to those requiring advanced and novel techniques. We compare different types of Procrustes problems, in which we align sets of points after various perturbations with respect to different choices of matrix norm, minimizing the difference after rotating one of these matrices by an orthogonal transformation. We highlight recent developments in nonsmooth Riemannian optimization and characterize which choices of norm work best for each perturbation. We show that in several applications, from low-dimensional alignments to more involved hypothesis testing for random networks, when Procrustes alignment with the spectral or robust norm is the appropriate choice, it is often feasible to replace the computationally more expensive spectral and robust minimizers with their closed-form Frobenius-norm counterpart. Our work reinforces the synergy between optimization, geometry, and statistics.

### How to use this code:

Besides for reproducibility, we want to provide this code for statisticians to easily compute the solutions to Procrustes problems under alternative choices of norm. To make use of either the spectral or robust Procrustes solution in your own R code, simply follow the procedure in any of the `powercurvesX.R` files. Here we compute the alignments based on Frobenius norm:

$$ p_0 = \mathrm{argmin}_{W\in \mathbb{R}^{d\times d}, W^\top W=I} \lVert \hat{X}-p_0 \hat{Y}\rVert_F. $$

Spectral, or 2-norm:

$$ p_S = \mathrm{argmin}_{W\in \mathbb{R}^{d\times d}, W^\top W=I} \lVert \hat{X}-p_0 \hat{Y}\rVert_2. $$

And "robust" norm, or sum of the row Euclidean norms $\lVert A\rVert_R:= \sum_{i=1}^n \lVert A_{i\cdot}\rVert_2$:

$$ p_R = \mathrm{argmin}_{W\in \mathbb{R}^{d\times d}, W^\top W=I} \lVert\hat{X}-p_0 \hat{Y}\rVert_R. $$

```
library("JuliaCall")
#julia_install_package_if_needed("Manopt")
#julia_install_package_if_needed("Manifolds")
#julia_install_package_if_needed("LinearAlgebra")
#julia_install_package_if_needed("RecursiveArrayTools")
julia_setup(JULIA_HOME=path_to_julia)
julia_source("procrustes_functions.jl")

procrustesFrob <- function(A,B) {
  #A, B are d x n
  #return W = argmin ||WB-A||_F
  tmp <- svd(A %*% t(B))
  return(tmp$u %*% t(tmp$v))
}

p0 <- procrustesFrob(t(Xhat),t(Yhat))
#Frobenius norm Procrustes solution
T1 <- norm(p0 %*% t(Yhat)-t(Xhat),type="F")
#Spectral norm with Frobenius norm Procrustes solution
T2 <- norm(p0 %*% t(Yhat)-t(Xhat),type="2")
#Robust norm with Frobenius norm Procrustes solution
T3 <- sum(apply(p0 %*% t(Yhat)-t(Xhat),2,function(x) sqrt(sum(x^2))))

#Starting the optimization from the Frobenius-norm solution generally gives good performance,
#but you can omit this if desired
pS <- julia_call("spectral_procrustes", t(Xhat), t(Yhat), p0)
#or pS <- julia_call("spectral_procrustes", t(Xhat), t(Yhat))
#Spectral norm with optimal minimizer
T3 <- norm(pS %*% t(Yhat)-t(Xhat),type="2")

#Robust norm with optimal minimizer
pR <- julia_call("robust_procrustes",t(Xhat),t(Yhat),p0)
T5 <- sum(apply(pR %*% t(Yhat)-t(Xhat),2,function(x) sqrt(sum(x^2))))
```

To compute the minimizer for another choice of norm, add the definition to `procrustes_functions.jl` (must be written in Julia). The following code provides the format of what must be added, but you will need to specify the appropriate definition of `myNew_norm`:

```
@doc raw"""
	myNew_cost(M, p, A, B)

Compute the MyNew norm ``\lVert p\cdot B- A\rVert_{MyNew}``.
"""
myNew_procrustes_cost(M, p, A, B) = myNew_norm(p*B - A)

"""
    myNew_procrustes(A,B[, p0]; kwargs...)

Compute the myNew procrustes minimizer, see [`myNew_procrustes_cost`](@ref).
All keyword arguments are passed to the MADS solver.

You can provide an initial guess `p0` for the solver to start at, but default the identity is used.
"""
function myNew_procrustes(A, B, p0=nothing; kwargs...)
    n,m = size(A)
    M = Rotations(n)
    p0 = isnothing(p0) ? Matrix{Float64}(I,n,n) : p0
    p = mesh_adaptive_direct_search(M, (M,p) -> myNew_procrustes_cost(M, p, A, B), p0; kwargs...)
    return p
end
```

You can then proceed as before, using either of the following lines in your R code:

```
pMy <- julia_call("myNew_procrustes",t(Xhat),t(Yhat),p0)
pMy <- julia_call("myNew_procrustes",t(Xhat),t(Yhat))
```

Happy aligning!
