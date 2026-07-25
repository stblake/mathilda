/* Mathilda — finite-difference stencils for the NDSolve method of lines.
 *
 * Weights are produced by Fornberg's algorithm (B. Fornberg, "Generation of
 * finite difference formulas on arbitrarily spaced grids", Math. Comp. 51
 * (1988) 699–706), which yields the coefficients for the m-th derivative at an
 * arbitrary evaluation point on an arbitrary node set — exact to the order the
 * node set allows.  On a uniform grid this gives central stencils of any even
 * order in the interior and one-sided stencils of the same order near the
 * boundaries. */
#ifndef NDSOLVE_STENCIL_H
#define NDSOLVE_STENCIL_H

/* Fornberg weights for the `m`-th derivative at point `z`, on the nodes
 * x[0..npts-1].  Writes w[0..npts-1] (the linear combination
 * sum_i w[i]*f(x[i]) approximates f^(m)(z)). */
void nd_fornberg_weights(int m, const double* x, int npts, double z, double* w);

/* Build a uniform-grid stencil for the `deriv`-th spatial derivative of
 * accuracy `order` at grid node `j` (0-based, of `nx` total nodes, spacing `h`).
 * The stencil is central where it fits and shifted (one-sided) near a boundary,
 * always keeping every node inside [0, nx-1] and preserving the requested order.
 * Writes the participating grid indices `idx[0..*n-1]` and their weights
 * `w[0..*n-1]` (already divided by h^deriv, so sum_k w[k]*u[idx[k]] ≈ ∂^deriv u).
 * `idx`/`w` must have room for at least deriv+order+1 entries; *n receives the
 * actual stencil width. */
void nd_stencil_build(int j, int nx, int deriv, int order, double h,
                      int* idx, double* w, int* n);

#endif /* NDSOLVE_STENCIL_H */
