/* Mathilda — Compile[]: delegated ND dispatch tables.
 *
 * The heads here already have an NDArray entry point in the interpreter, and
 * each of those takes the WHOLE call and promises a result identical to the
 * equivalent List call.  So the compiled form is that same function, called
 * from the VM — which makes the compiled subset of these heads the interpreted
 * one by construction.  This file holds the three dispatch tables (ND_FNS /
 * ND_REDS / ND_FN2S — single-array structural, array→scalar reductions, and
 * two-array heads), the small adapter thunks that fit each interpreter entry
 * point to the table's `Expr* (*)(Expr*)` shape, and the lookup + static
 * result-type helpers the inference pass and emitter call.  The spec structs
 * and the lookup/result declarations live in compile_internal.h — the VM
 * dereferences a spec through an instruction immediate at run time. */
#include "compile_internal.h"      /* Nd*Spec structs, nd_*_lookup/result decls */
#include "../arithmetic.h"
#include "../ndarray.h"            /* EXPR_NDARRAY, builtin matrix producers */
#include "../ndreduce.h"           /* ndred_* reductions */
#include "../ndstruct.h"           /* ndstruct_* structural heads, ND_MOST/ND_REST */
#include "../linalg/ndlinalg.h"    /* ndla_* single/rank-2 linalg */
#include "../linalg/linalg.h"      /* builtin_matrixpower, nd_dot_machine */
#include "../list/transpose.h"     /* builtin_conjugate_transpose */
#include "../sort.h"               /* builtin_reverse_sort */
#include "../convolutions.h"       /* builtin_list_convolve / _correlate */
#include "../fourier.h"            /* fourier_compile / inverse_fourier_compile, DCT/DST */
#include <string.h>
#include <stddef.h>
#include <math.h>     /* INFINITY / NAN — the V_NORM p encoding */

/* ---- delegated structural heads -------------------------------------------
 *
 * These already have an NDArray entry point in the interpreter, and each of
 * those takes the WHOLE call and promises a result identical to the equivalent
 * List call (src/ndstruct.h:14, src/ndreduce.h:20).  So the compiled form is
 * that same function, called from the VM — which makes the compiled subset of
 * these heads the interpreted one by construction, exactly as A_PART does for
 * the general Part specs.  The win is not the operation (it was already a fast
 * buffer walk) but that a body CONTAINING one no longer bails wholesale.
 *
 * `rank_rule`: 0 = same rank as the operand, 1 = rank 1 (Flatten), 2 = rank 2
 * in and out (Transpose), 3 = rank 2 in -> rank 1 out (Diagonal: matrix ->
 * vector), 4 = rank 1 in and out (Normalize), 5 = rank 1 in -> rank 2 out
 * (the matrix producers DiagonalMatrix / HankelMatrix / ToeplitzMatrix /
 * VandermondeMatrix — COMPILE_MISSING.md §5).  Each entry preserves the element
 * type, restricted to the `elems` gate below (0 = any), unless `complex_result`

/* ndstruct_rotate / _head_tail / _take_extreme take a second C-level argument
 * that selects the direction, so each spelling needs its own one-line adapter
 * to fit the table's `Expr* (*)(Expr*)` shape. */
static Expr* nd_rotate_left(Expr* r)  { return ndstruct_rotate(r, true); }
static Expr* nd_rotate_right(Expr* r) { return ndstruct_rotate(r, false); }
static Expr* nd_most(Expr* r)         { return ndstruct_head_tail(r, ND_MOST); }
static Expr* nd_rest(Expr* r)         { return ndstruct_head_tail(r, ND_REST); }
static Expr* nd_take_largest(Expr* r) { return ndstruct_take_extreme(r, true); }
static Expr* nd_take_smallest(Expr* r){ return ndstruct_take_extreme(r, false); }

/* COMPILE_MISSING.md §2/§3 adapters.  Two shapes need one:
 *   - a DIRECT-args entry point (ndla_pseudoinverse_direct / _leastsquares_direct
 *     take the extracted matrix + Tolerance values, not the call);
 *   - a call-shaped entry point that PRINTS a diagnostic or runs an expensive
 *     delist-and-reeval on a structural mismatch its fast path cannot use.  A
 *     cheap dimension check declines (NULL) BEFORE the noisy path, so a compiled
 *     degrade stays silent — the interpreter fallback then reports it once.
 * Each fits the table's `Expr* (*)(Expr*)` shape, BORROWS the call (the VM frees
 * it) and returns a freshly-owned result or NULL/non-NDArray to decline. */
static bool nd_square_matrix(const Expr* a) {
    return a && a->type == EXPR_NDARRAY && a->data.ndarray.rank == 2
        && a->data.ndarray.dims[0] == a->data.ndarray.dims[1];
}
static bool nd_rank2(const Expr* a) {
    return a && a->type == EXPR_NDARRAY && a->data.ndarray.rank == 2;
}
static bool nd_vec3(const Expr* a) {
    return a && a->type == EXPR_NDARRAY && a->data.ndarray.rank == 1
        && a->data.ndarray.dims[0] == 3;
}
/* §2 single-array */
static Expr* nd_inverse(Expr* call) {          /* non-square declines silently */
    if (call->data.function.arg_count != 1
        || !nd_square_matrix(call->data.function.args[0])) return NULL;
    return ndla_inverse(call);
}
static Expr* nd_matrixpower(Expr* call) {      /* MatrixPower[m, e]: m square */
    if (call->data.function.arg_count != 2
        || !nd_square_matrix(call->data.function.args[0])) return NULL;
    return builtin_matrixpower(call);
}
static Expr* nd_pseudoinverse(Expr* call) {    /* direct-args, real float64 only */
    if (call->data.function.arg_count != 1) return NULL;
    const Expr* a = call->data.function.args[0];
    if (!nd_rank2(a)) return NULL;
    return ndla_pseudoinverse_direct(a, true, 0.0);
}
/* §3 two-array */
static Expr* nd_dot2(Expr* call) {             /* BLAS-first; scalar or matrix out */
    if (call->data.function.arg_count != 2) return NULL;
    return nd_dot_machine(call->data.function.args[0], call->data.function.args[1]);
}
static Expr* nd_linearsolve(Expr* call) {      /* non-conformable declines silently */
    if (call->data.function.arg_count != 2) return NULL;
    const Expr* a = call->data.function.args[0];
    const Expr* b = call->data.function.args[1];
    if (!nd_rank2(a) || !b || b->type != EXPR_NDARRAY
        || b->data.ndarray.dims[0] != a->data.ndarray.dims[0]) return NULL;
    return ndla_linearsolve(call);
}
static Expr* nd_cross(Expr* call) {            /* two length-3 vectors */
    if (call->data.function.arg_count != 2
        || !nd_vec3(call->data.function.args[0])
        || !nd_vec3(call->data.function.args[1])) return NULL;
    return ndla_cross(call);
}
static Expr* nd_leastsquares(Expr* call) {     /* direct-args, real float64 only */
    if (call->data.function.arg_count != 2) return NULL;
    const Expr* a = call->data.function.args[0];
    const Expr* b = call->data.function.args[1];
    if (!nd_rank2(a) || !b || b->type != EXPR_NDARRAY) return NULL;
    return ndla_leastsquares_direct(a, b, true, 0.0);
}

static const NdFnSpec ND_FNS[] = {
    { "Reverse",    ndstruct_reverse,   0, 0 },
    { "Sort",       ndstruct_sort,      0, 0 },   /* 1-arg only: a comparator is a
                                                   * function value the ND path
                                                   * cannot call back into */
    { "Ordering",   ndstruct_ordering,  0, 0, true },  /* argsort -> int64 positions;
                                                   * int_result forces the CT_INT
                                                   * element type (a permutation is
                                                   * integer for any input dtype) */
    { "Accumulate", ndred_accumulate,   0, 0 },
    { "Flatten",    ndstruct_flatten,   0, 1 },
    { "Transpose",  ndstruct_transpose, 0, 2 },
    /* Diagonal[m] and Diagonal[m, k]: a rank-2 buffer to its k-th diagonal (a
     * rank-1 vector).  Two rows for the two arities nd_fn_lookup keys on
     * (na == 1 + nextra); the optional k is the trailing integer argument. */
    { "Diagonal",   ndstruct_diagonal,  0, 3 },
    { "Diagonal",   ndstruct_diagonal,  1, 3 },
    { "Take",       ndstruct_take,      1, 0 },
    { "Drop",       ndstruct_drop,      1, 0 },
    /* Added 2026-08-02 by the sixth sweep (tools/compile_coverage.py).  Each
     * already had a buffer path in the interpreter and no lowering here, which
     * is the contradiction that audit exists to find: a head fast at the REPL
     * and slow inside the Compile[] meant to speed it up — and, because the
     * subset is a cliff, dragging every other head in the body down with it.
     *
     * All seven answer with an array of the SAME rank and element type, so
     * rank_rule 0 covers them.  Where a head is exact-only over integers
     * (Ratios of an integer vector is a list of Rationals) its ND entry point
     * already declines by handing back a nested List, which the VM reads as
     * "not an NDArray" and aborts to the interpreter — the same faithful
     * degrade the existing seven rely on. */
    { "Differences",   ndstruct_differences, 0, 0 },
    { "Ratios",        ndstruct_ratios,      0, 0 },
    { "Most",          nd_most,              0, 0 },
    { "Rest",          nd_rest,              0, 0 },
    { "Clip",          ndstruct_clip,        0, 0 },
    { "RotateLeft",    nd_rotate_left,       1, 0 },
    { "RotateRight",   nd_rotate_right,      1, 0 },
    { "MovingAverage", ndred_moving_average, 1, 0 },
    { "MovingMedian",  ndred_moving_median,  1, 0 },
    { "TakeLargest",   nd_take_largest,      1, 0 },
    { "TakeSmallest",  nd_take_smallest,     1, 0 },
    /* COMPILE_MISSING.md §2: single array -> array linalg heads, each delegating
     * to its interpreter entry point (an adapter where that entry takes direct
     * args or prints on a structural mismatch — see the adapter block above).
     * `elems` bars an operand dtype the fast path would CONVERT: Inverse /
     * Normalize / MatrixPower / PseudoInverse of an int operand are exact
     * Rationals/reals no int slot holds, so real (+complex where supported)
     * only; ReverseSort's buffer path declines complex (int+real); Conjugate-
     * Transpose preserves any dtype.  rank_rule 4 = require rank 1 -> rank 1. */
    { "Inverse",            nd_inverse,                   0, 2, false, NDF_REAL | NDF_CPLX },
    { "Normalize",          ndla_normalize,               0, 4, false, NDF_REAL | NDF_CPLX },
    { "MatrixPower",        nd_matrixpower,               1, 2, false, NDF_REAL | NDF_CPLX },
    { "ReverseSort",        builtin_reverse_sort,         0, 0, false, NDF_INT | NDF_REAL },
    { "ConjugateTranspose", builtin_conjugate_transpose,  0, 2, false, NDF_INT | NDF_REAL | NDF_CPLX },
    { "PseudoInverse",      nd_pseudoinverse,             0, 2, false, NDF_REAL },
    /* COMPILE_MISSING.md §4: the transforms.  Fourier / InverseFourier are
     * real-in complex-out — the ONLY heads that need complex_result (their
     * fourier_compile wrapper always builds NDT_COMPLEX64 so the buffer matches
     * the CT_COMPLEX promise; the interpreter's data-dependent real-collapse is
     * not statically typeable).  rank_rule 0 preserves rank, so 1-D and n-D FFT
     * both lower.  FourierDCT / FourierDST are real->real (deterministically for
     * a real operand), so they delegate straight to their builtins and are gated
     * NDF_REAL — an int operand would mis-promise CT_INT for a real buffer, a
     * complex operand would break the real result promise. */
    { "Fourier",            fourier_compile,          0, 0, false, NDF_INT | NDF_REAL | NDF_CPLX, true },
    { "InverseFourier",     inverse_fourier_compile,  0, 0, false, NDF_INT | NDF_REAL | NDF_CPLX, true },
    { "FourierDCT",         builtin_fourier_dct,      0, 0, false, NDF_REAL },
    { "FourierDST",         builtin_fourier_dst,      0, 0, false, NDF_REAL },
    /* COMPILE_MISSING.md §5: the matrix producers.  rank_rule 5 = rank 1 -> rank
     * 2; the delegate rebuilds H[ndarray] and the interpreter builtin returns a
     * packed rank-2 NDArray (DiagonalMatrix fills its buffer directly; Hankel /
     * Toeplitz / Vandermonde now do too — see their ndbuild_open branch).  This
     * needs no producing opcode: A_NDFN already stores whatever NDArray the
     * delegate returns, so the runtime shape is the delegate's; only the STATIC
     * rank+element rule lives here.  Gated NDF_INT|NDF_REAL: a complex operand
     * keeps the interpreter's exact delist path.  Single-vector form only (what
     * the coverage probe H[v] tests); the two-vector Hankel/Toeplitz[c, r] form
     * is a rank-1 x rank-1 -> rank-2 A_NDFN2 lowering, not done here. */
    { "DiagonalMatrix",     builtin_diagonalmatrix,    0, 5, false, NDF_INT | NDF_REAL },
    { "HankelMatrix",       builtin_hankelmatrix,      0, 5, false, NDF_INT | NDF_REAL },
    { "ToeplitzMatrix",     builtin_toeplitzmatrix,    0, 5, false, NDF_INT | NDF_REAL },
    { "VandermondeMatrix",  builtin_vandermondematrix, 0, 5, false, NDF_INT | NDF_REAL },
    /* Covariance[a] / Correlation[a]: the one-argument auto form, a rank-2 real
     * matrix to its p x p (co)variance matrix (rank_rule 2 = rank 2 in and out;
     * the n x p -> p x p dim change is fine, the register tracks only rank).  The
     * two-argument forms are the R2_DOT rows in ND_FN2S below.  NDF_REAL only: an
     * integer matrix's covariance is a matrix of Rationals no float slot holds, so
     * it declines to the interpreter's exact path (nd_covariance delists). */
    { "Covariance",         nd_covariance,             0, 2, false, NDF_REAL },
    { "Correlation",        nd_correlation,            0, 2, false, NDF_REAL },
};

/* ---- delegated array -> scalar reductions ---------------------------------
 *
 * The counterpart of ND_FNS for the heads that collapse an array to one
 * number.  Same contract and the same reason: each entry point takes the whole
 * call and promises the List call's answer, so the compiled subset is the
 * interpreted one and the ROUNDING is identical too — a compiled Mean that
 * summed in its own order would disagree with the interpreter in the last bits
 * for no benefit.
 *
 * Real element type only.  Mean, Variance, StandardDeviation, RootMeanSquare
 * and Median of an exact-integer vector are Rationals (Mean[{1, 2}] is 3/2,
 * Variance[Range[10]] is 55/6) and no machine slot holds one; the interpreter
 * is right to answer exactly, so the compiler declines rather than silently
 * rounding.  Max and Min are the exceptions — they SELECT an element, so an
 * integer vector gives an Integer — and carry int_ok to say so. */

static const NdRedSpec ND_REDS[] = {
    { "Mean",              ndred_mean,      false, 1, false },
    { "Median",            ndred_median,    false, 1, false },
    { "Variance",          ndred_variance,  false, 1, false },
    /* Standardized central moments (real vector -> Real scalar; an int vector's
     * value is a radical, so it declines to the interpreter, like Variance). */
    { "Skewness",          ndred_skewness,  false, 1, false },
    { "Kurtosis",          ndred_kurtosis,  false, 1, false },
    { "StandardDeviation", ndred_std,       false, 1, false },
    { "RootMeanSquare",    ndred_rms,       false, 1, false },
    { "Max",               ndred_max,       true,  1, false },
    { "Min",               ndred_min,       true,  1, false },
    /* Order statistics with a trailing rank (nextra 1): RankedMin[v, n] /
     * RankedMax[v, n] SELECT an element, so an int vector yields an Integer
     * (int_ok), like Max/Min.  The sign of n is a runtime register value, so it
     * does not affect the static result type. */
    { "RankedMin",         ndred_ranked_min, true, 1, false, 1 },
    { "RankedMax",         ndred_ranked_max, true, 1, false, 1 },
    /* Moment[v, r] (raw (1/n)Sum[x^r]) and CentralMoment[v, r] (moment about the
     * mean): the same array+trailing-int shape as Ranked* (nextra 1). A real
     * vector answers a Real; an int vector's exact answer is a Rational no machine
     * slot holds, so both decline (int_ok false) and the interpreter answers
     * exactly, exactly like Variance/Mean. */
    { "Moment",            ndred_moment,         false, 1, false, 1 },
    { "CentralMoment",     ndred_central_moment, false, 1, false, 1 },
    /* Linalg rank-2 -> scalar reductions (COMPILE_MISSING.md §1), each delegating
     * to its ndla_* interpreter entry point (src/linalg/ndlinalg.h).  Real operand
     * only: an exact-integer Tr / Det / Norm is an exact number no machine slot
     * holds, so int (and complex) operands decline to the interpreter, which is
     * right to answer exactly.  A real matrix already routes through these same
     * ndla_* functions in the interpreter (all four are on pack.c's AWARE list),
     * so the compiled answer is identical by construction.  MatrixRank answers an
     * Integer for a real matrix (int_result).  Norm accepts rank 1 or 2 (rank 0). */
    { "Tr",                ndla_tr,         false, 2, false },
    { "Det",               ndla_det,        false, 2, false },
    { "Norm",              ndla_norm,       false, 0, false },
    { "MatrixRank",        ndla_matrixrank, false, 2, true  },
};

/* Name accessors for the disassembler; see compile_internal.h. */
const char* nd_fn_head_name(const void* imm) {
    return imm ? ((const NdFnSpec*)imm)->head : NULL;
}
const char* nd_red_head_name(const void* imm) {
    return imm ? ((const NdRedSpec*)imm)->head : NULL;
}

const NdRedSpec* nd_red_lookup(const char* h, size_t na) {
    /* na == 1 + nextra: the plain reductions take just the array (Max[x, y] is
     * the scalar opcode, not this); RankedMin/Max take the array plus one int. */
    for (size_t i = 0; i < sizeof ND_REDS / sizeof ND_REDS[0]; i++)
        if (strcmp(h, ND_REDS[i].head) == 0 && na == 1u + (size_t)ND_REDS[i].nextra)
            return &ND_REDS[i];
    return NULL;
}

/* Scalar result type of a delegated reduction over an operand of type `ta`,
 * or CT_ERR.  The required rank is per-entry (`s->rank`): the statistics
 * reductions are rank 1 only — at rank 2 the interpreter's Mean is a VECTOR of
 * column means, not a number, and the two must not disagree — while the linalg
 * reductions (Tr / Det / MatrixRank) are rank 2 and Norm is either (rank 0). */
CompileType nd_red_result(const NdRedSpec* s, CompileType ta) {
    if (!CT_IS_ARRAY(ta)) return CT_ERR;
    int rk = CT_RANK(ta);
    if (s->rank == 0) { if (rk != 1 && rk != 2) return CT_ERR; }  /* Norm: 1 or 2 */
    else if (rk != s->rank) return CT_ERR;
    CompileType el = CT_ELEM(ta);
    if (s->int_result) return (el == CT_REAL) ? CT_INT : CT_ERR;  /* MatrixRank: real -> Integer */
    if (el == CT_REAL) return CT_REAL;
    if (el == CT_INT)  return s->int_ok ? CT_INT : CT_ERR;
    return CT_ERR;                       /* complex: no ND reduction promises it */
}

/* Norm[array, p]: is the two-arg form with this compile-time literal p, over an
 * operand of type `ta`, one that ndla_norm answers as a machine Real?  If so
 * return true and set *imm_r to the V_NORM encoding of p (opcode doc: +Inf for
 * Infinity, NaN for "Frobenius", else the positive value).  The validity mirrors
 * ndla_norm/ndla_matrix_norm_direct exactly (src/linalg/ndlinalg.c):
 *   - rank 1 (vector): any positive Integer/Real p, or Infinity.
 *   - rank 2 (matrix): p in {1, 2, Infinity, "Frobenius"} (an Integer, since the
 *     LAPACK path keys on Integer 1/2; a Real or other integer declines).
 * Real operand only — an exact-integer or complex Norm is not a machine real, so
 * it declines to the interpreter, exactly as the bare Norm[array] reduction does. */
bool nd_norm2_encode(const Expr* p, CompileType ta, double* imm_r) {
    if (!CT_IS_ARRAY(ta) || CT_ELEM(ta) != CT_REAL) return false;
    int rank = CT_RANK(ta);
    if (rank != 1 && rank != 2) return false;
    if (!p) return false;
    if (p->type == EXPR_SYMBOL && p->data.symbol.name
        && strcmp(p->data.symbol.name, "Infinity") == 0) { *imm_r = INFINITY; return true; }
    if (p->type == EXPR_STRING && p->data.string
        && strcmp(p->data.string, "Frobenius") == 0) {
        if (rank != 2) return false;                    /* ndla rank-1 rejects a string p */
        *imm_r = NAN; return true;
    }
    if (p->type == EXPR_INTEGER && p->data.integer > 0) {
        if (rank == 2 && p->data.integer != 1 && p->data.integer != 2) return false;
        *imm_r = (double)p->data.integer; return true;
    }
    if (p->type == EXPR_REAL && p->data.real > 0.0) {
        if (rank == 2) return false;                    /* ndla matrix norm: Integer p only */
        *imm_r = p->data.real; return true;
    }
    return false;
}

const NdFnSpec* nd_fn_lookup(const char* h, size_t na) {
    for (size_t i = 0; i < sizeof ND_FNS / sizeof ND_FNS[0]; i++)
        if (strcmp(h, ND_FNS[i].head) == 0 && na == 1u + (size_t)ND_FNS[i].nextra)
            return &ND_FNS[i];
    return NULL;
}

/* Result type of a delegated head over an operand of type `ta`, or CT_ERR.
 * `elems` (when non-zero) gates the operand's element type — a head is declined
 * for a dtype its fast path would convert, since A_NDFN cannot check the result
 * dtype against the promise.  `int_result` forces CT_INT (Ordering returns an
 * int64 permutation whatever the input dtype); complex_result forces CT_COMPLEX
 * (Fourier: real in, complex out); otherwise the operand's element type is
 * preserved.  rank_rule: 0 = same rank, 1 = rank 1, 2 = rank 2 in/out,
 * 3 = rank 2 -> rank 1 (Diagonal), 4 = rank 1 in/out (Normalize),
 * 5 = rank 1 -> rank 2 (matrix producers). */
CompileType nd_fn_result(const NdFnSpec* s, CompileType ta) {
    if (!CT_IS_ARRAY(ta)) return CT_ERR;
    if (s->elems && !((s->elems >> (unsigned)CT_ELEM(ta)) & 1u)) return CT_ERR;  /* operand gate */
    int rank = CT_RANK(ta);
    CompileType elem = s->complex_result ? CT_COMPLEX
                     : s->int_result     ? CT_INT
                     : CT_ELEM(ta);
    if (s->rank_rule == 1) return CT_ARRAY(elem, 1);
    if (s->rank_rule == 2) return rank == 2 ? CT_ARRAY(elem, 2) : CT_ERR;
    if (s->rank_rule == 3) return rank == 2 ? CT_ARRAY(elem, 1) : CT_ERR;  /* Diagonal */
    if (s->rank_rule == 4) return rank == 1 ? CT_ARRAY(elem, 1) : CT_ERR;  /* Normalize */
    if (s->rank_rule == 5) return rank == 1 ? CT_ARRAY(elem, 2) : CT_ERR;  /* matrix producers */
    return CT_ARRAY(elem, rank);
}

/* ---- delegated TWO-array heads (COMPILE_MISSING.md §3) --------------------
 *
 * The binary counterpart of ND_FNS: two machine arrays in, an array (opcode
 * A_NDFN2) or a scalar (V_NDFN2 — only Dot's vector.vector) out.  Same contract
 * as the single-array table — each `fn` takes the whole rebuilt call and either
 * returns the interpreter's own answer or declines — so the compiled subset is
 * the interpreted one by construction.  The scalar branch exists because a body
 * containing `a.b` (an inner product) would otherwise bail wholesale: the subset
 * is a cliff.
 *
 * `rank2`: the result rank as a function of the operand ranks (ra, rb):
 *   R2_DOT    ra + rb - 2  (matrix.matrix -> 2, matrix.vector -> 1,
 *             vector.vector -> 0 == a SCALAR, which the emitter routes to V_NDFN2)
 *   R2_SOLVE  result rank = rb, requires ra == 2  (LinearSolve, LeastSquares)
 *   R2_VEC    requires ra == rb == 1, result rank 1  (Cross, ListConvolve/Correlate)
 *   R2_JOIN   requires ra == rb, result rank = ra
 *   R2_MATRIX requires ra == rb == 1, result rank 2  (the two-vector matrix
 *             PRODUCERS HankelMatrix[c, r] / ToeplitzMatrix[c, r] — a rank-1
 *             column and a rank-1 row build a rank-2 matrix; COMPILE_MISSING.md §5)
 * `elem`: the result element type: PROMOTE (complex iff either operand is
 *   complex, else real) for the numeric contractions; SAME (require equal
 *   operand dtypes, preserve) for the structural Join and the matrix producers
 *   (int+int -> int64, real+real -> real; a mixed int/real pair declines to the
 *   interpreter, which coerces it exactly).  `elems` gates BOTH operands exactly
 *   as NdFnSpec.elems gates the one. */

static const NdFn2Spec ND_FN2S[] = {
    { "Dot",           nd_dot2,                NDF_REAL | NDF_CPLX,           R2_DOT,   NDF2_PROMOTE },
    { "LinearSolve",   nd_linearsolve,         NDF_REAL | NDF_CPLX,           R2_SOLVE, NDF2_PROMOTE },
    { "Cross",         nd_cross,               NDF_REAL | NDF_CPLX,           R2_VEC,   NDF2_PROMOTE },
    { "LeastSquares",  nd_leastsquares,        NDF_REAL,                      R2_SOLVE, NDF2_PROMOTE },
    { "ListConvolve",  builtin_list_convolve,  NDF_REAL | NDF_CPLX,           R2_VEC,   NDF2_PROMOTE },
    { "ListCorrelate", builtin_list_correlate, NDF_REAL | NDF_CPLX,           R2_VEC,   NDF2_PROMOTE },
    { "Join",          ndstruct_join,          NDF_INT | NDF_REAL | NDF_CPLX, R2_JOIN,  NDF2_SAME },
    /* COMPILE_MISSING.md §5, two-vector form: HankelMatrix[c, r] /
     * ToeplitzMatrix[c, r] build a rank-2 matrix from a rank-1 column and a
     * rank-1 row.  Delegated to the builtins (which delist the two NDArrays to
     * the two-vector form and write the result buffer directly — see hk_build /
     * tz_build).  Gated NDF_INT|NDF_REAL and NDF2_SAME to match the single-vector
     * rows: int+int -> int64, real+real -> real; a mixed int/real pair declines
     * to the interpreter (which coerces to real), and a complex pair likewise. */
    { "HankelMatrix",   builtin_hankelmatrix,   NDF_INT | NDF_REAL, R2_MATRIX, NDF2_SAME },
    { "ToeplitzMatrix", builtin_toeplitzmatrix, NDF_INT | NDF_REAL, R2_MATRIX, NDF2_SAME },
    /* Covariance[v, w] / Covariance[a, b] (and Correlation): the rank arithmetic
     * is exactly Dot's — two rank-1 vectors give a scalar (R2_DOT rr = 0 -> a
     * V_NDFN2 scalar), two rank-2 matrices give a rank-2 matrix (rr = 2 -> an
     * A_NDFN2 array), so R2_DOT covers both without a bespoke rank rule.  The
     * n x p -> p x q dim change is invisible to the register.  NDF_REAL only (a
     * complex or integer pair keeps the interpreter's exact/complex path); the fn
     * is the same nd_covariance / nd_correlation the REPL uses, so the compiled
     * answer is identical by construction. */
    { "Covariance",     nd_covariance,          NDF_REAL,           R2_DOT,    NDF2_PROMOTE },
    { "Correlation",    nd_correlation,         NDF_REAL,           R2_DOT,    NDF2_PROMOTE },
};

/* Name accessor for the disassembler; see compile_internal.h. */
const char* nd_fn2_head_name(const void* imm) {
    return imm ? ((const NdFn2Spec*)imm)->head : NULL;
}

const NdFn2Spec* nd_fn2_lookup(const char* h, size_t na) {
    if (na != 2) return NULL;
    for (size_t i = 0; i < sizeof ND_FN2S / sizeof ND_FN2S[0]; i++)
        if (strcmp(h, ND_FN2S[i].head) == 0) return &ND_FN2S[i];
    return NULL;
}

/* Result type of a delegated two-array head over operands `ta`, `tb`, or CT_ERR.
 * A rank-0 result is a SCALAR CompileType (Dot's vector.vector), which the emit
 * path recognises via !CT_IS_ARRAY to pick V_NDFN2 over A_NDFN2. */
CompileType nd_fn2_result(const NdFn2Spec* s, CompileType ta, CompileType tb) {
    if (!CT_IS_ARRAY(ta) || !CT_IS_ARRAY(tb)) return CT_ERR;
    unsigned ea = (unsigned)CT_ELEM(ta), eb = (unsigned)CT_ELEM(tb);
    if (!((s->elems >> ea) & 1u) || !((s->elems >> eb) & 1u)) return CT_ERR;  /* gate both */
    int ra = CT_RANK(ta), rb = CT_RANK(tb);

    CompileType elem;
    if (s->elem == NDF2_SAME) {
        if (ea != eb) return CT_ERR;                 /* Join needs equal operand dtypes */
        elem = CT_ELEM(ta);
    } else {                                          /* PROMOTE */
        elem = (ea == (unsigned)CT_COMPLEX || eb == (unsigned)CT_COMPLEX) ? CT_COMPLEX : CT_REAL;
    }

    int rr;
    switch (s->rank2) {
        case R2_DOT:   if (ra < 1 || rb < 1) return CT_ERR; rr = ra + rb - 2; break;
        case R2_SOLVE: if (ra != 2 || rb < 1 || rb > 2) return CT_ERR; rr = rb; break;
        case R2_VEC:   if (ra != 1 || rb != 1) return CT_ERR; rr = 1; break;
        case R2_JOIN:  if (ra != rb) return CT_ERR; rr = ra; break;
        case R2_MATRIX: if (ra != 1 || rb != 1) return CT_ERR; rr = 2; break;  /* Hankel/Toeplitz[c, r] */
        default: return CT_ERR;
    }
    if (rr < 0 || rr > CT_MAX_RANK) return CT_ERR;
    if (rr == 0) return elem;                         /* scalar -> V_NDFN2 */
    return CT_ARRAY(elem, rr);
}
