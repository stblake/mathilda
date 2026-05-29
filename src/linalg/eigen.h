#ifndef MATEIGEN_H
#define MATEIGEN_H

#include <stddef.h>
#include <stdint.h>
#include "expr.h"
#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* Eigenvalues / Eigenvectors of a square matrix.
 *
 * Symbolic path (default for non-numeric input):
 *   - Compute characteristic polynomial det(lambda*I - m) via
 *     Faddeev-Leverrier-Souriau for the ordinary case, or Laplace
 *     expansion of det(m - lambda*a) for the generalised case.
 *   - Route the resulting univariate polynomial through the public
 *     Solve builtin so its rationalise -> solve -> numericalize
 *     pipeline handles approximate (Real / MPFR) input automatically.
 *
 * Numeric Method dispatch (Phase 1: parsing + warning, real kernels
 * arrive in subsequent phases per lovely-roaming-diffie.md):
 *   - "Direct"   Hessenberg + implicit QR (general) / tridiag + symmetric
 *                QR (Hermitian).  Returns all eigenvalues.
 *   - "Arnoldi"  Krylov-subspace iteration; finds k extreme eigenvalues.
 *   - "Banded"   Hermitian only; reduces band to tridiagonal then QR.
 *   - "FEAST"    Hermitian only; eigenvalues in a user-specified interval.
 *   - Automatic  picks Direct unless k is small (-> Arnoldi) or matrix is
 *                Hermitian-banded (-> Banded).
 */
Expr* builtin_eigenvalues(Expr* res);
Expr* builtin_eigenvectors(Expr* res);

/* Registers Eigenvalues and Eigenvectors builtins and their attributes.
 * Must be called from core_init() (after symtab_init()). */
void  mateigen_init(void);

/* Method selector for the numeric eigenvalue dispatcher. */
typedef enum {
    MATEIGEN_AUTOMATIC = 0,
    MATEIGEN_DIRECT,
    MATEIGEN_ARNOLDI,
    MATEIGEN_BANDED,
    MATEIGEN_FEAST,
    MATEIGEN_METHOD_UNKNOWN  /* user-supplied string not recognised */
} MateigenMethod;

/* Parse a Method -> <value> rule into a MateigenMethod.  Accepts:
 *   - the symbol Automatic
 *   - the strings "Direct", "Arnoldi", "Banded", "FEAST"
 *   - any list whose head is one of the above strings (sub-options
 *     forms like {"Arnoldi", "MaxIterations" -> 100} -- only the head
 *     classifies the method here; sub-options are interpreted later).
 * Anything else maps to MATEIGEN_METHOD_UNKNOWN.  The argument must be
 * the right-hand side of the Rule (NOT the whole Rule expression). */
MateigenMethod mateigen_parse_method_value(Expr* method_value);

/* Flag telling a per-method kernel which outputs to populate.  Both
 * Eigenvalues and Eigenvectors route through the same kernel; the WANT
 * flag tells the kernel whether to skip the (potentially expensive)
 * eigenvector accumulation work.  This is the abstraction that LAPACK
 * exposes via its JOBV / JOBVL / JOBVR flags. */
typedef enum {
    MATEIGEN_WANT_VALUES  = 1u << 0,
    MATEIGEN_WANT_VECTORS = 1u << 1,
    MATEIGEN_WANT_BOTH    = MATEIGEN_WANT_VALUES | MATEIGEN_WANT_VECTORS
} MateigenWant;

#ifdef USE_MPFR
/* Compute all eigenvalues of a square n*n real matrix at MPFR precision.
 *
 * Black-box wrapper around the internal Hessenberg + Francis QR pipeline,
 * exposed for callers like N[Root[..]] that build a Frobenius companion
 * matrix and treat eigenvalue extraction as their all-roots backend.
 *
 * `A` is row-major, length n*n, with each cell mpfr_init2'd to `bits`.
 * The matrix is mutated. On return, eval_re[i] + I*eval_im[i] are the n
 * eigenvalues in arbitrary order — the caller is expected to sort them.
 *
 * Both eval_re and eval_im must be arrays of n already-mpfr_init2'd
 * cells at `bits` precision on entry. Returns 0 on success, non-zero on
 * QR non-convergence (exceedingly rare for well-conditioned matrices). */
int eigen_all_eigenvalues_real_mpfr(mpfr_t* A, size_t n, mpfr_prec_t bits,
                                    mpfr_t* eval_re, mpfr_t* eval_im);
#endif

#endif /* MATEIGEN_H */
