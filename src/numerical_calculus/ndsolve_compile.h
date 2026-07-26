/* Mathilda — numeric compiler for the NDSolve reduced RHS.
 *
 * When the compiled linear operator does not apply (a nonlinear or otherwise
 * non-affine right-hand side), the reduced RHS f_i(t, Y) is still a fixed
 * expression tree in the reduced-state symbols NDSolve`w<k>, the time variable,
 * and numeric constants.  This module compiles each f_i into a compact
 * stack-machine program so the time integrator evaluates the RHS by running
 * bytecode over the numeric state vector — no symbol-table binding, no Expr
 * copying, no `numericalize` per call.  It is the nonlinear counterpart of the
 * linear operator fast path; the symbolic sampler remains the fallback for any
 * construct the compiler does not recognize (it simply returns NULL).
 *
 * The compiler also records, per component, which state indices it reads, so a
 * Jacobian can be formed by Curtis–Powell–Reid colored finite differences
 * (O(bandwidth) RHS evaluations instead of O(d)). */
#ifndef NDSOLVE_COMPILE_H
#define NDSOLVE_COMPILE_H

#include "ndsolve_common.h"

/* Opaque compiled RHS. */
typedef struct NdCompiled NdCompiled;

/* Compile every P->f[i] into bytecode.  Returns NULL if any component uses a
 * construct the VM does not support (caller keeps using the symbolic sampler).
 * Borrows P (reads P->f, P->ysym, P->tvar, P->d); the result is independent of
 * P afterwards and must be freed with nd_compiled_free. */
NdCompiled* nd_compile_rhs(const NdProblem* P);

/* out[0..d-1] = f(t, Y).  Returns false if any component evaluates to a
 * non-finite value (treated as a sampling failure, like the symbolic path). */
bool nd_compiled_eval(NdCompiled* C, double t, const double* Y, double* out);

/* Dense row-major Jacobian Jout[i*d + j] = d f_i/d y_j via colored central
 * finite differences over the compiled RHS.  Returns false on a non-finite
 * evaluation. */
bool nd_compiled_jacobian(NdCompiled* C, double t, const double* Y, double* Jout);

/* Number of CPR colors (== FD RHS-evaluation pairs per Jacobian); for reporting. */
int nd_compiled_ncolor(const NdCompiled* C);

void nd_compiled_free(NdCompiled* C);

#endif /* NDSOLVE_COMPILE_H */
