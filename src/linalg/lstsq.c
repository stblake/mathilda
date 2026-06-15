/* matlstsq.c
 *
 * LeastSquares[m, b] -- linear least-squares solver.
 *
 *   LeastSquares[m, b]
 *       Returns an x that minimises Norm[m . x - b].  When m has full
 *       column rank the minimiser is unique; when m is rank-deficient
 *       the result is the minimum-norm minimiser
 *       (PseudoInverse[m] . b).
 *
 *   LeastSquares[m, b, Method -> "<name>"]
 *   LeastSquares[m, b, Tolerance -> t]
 *       Optional Method and Tolerance arguments may appear in either
 *       order and both may be given.
 *
 * Method names (parsed case-sensitively; "Automatic" symbol also
 * accepted, matching the Mathematica grammar):
 *
 *   "Automatic"           -- alias for "Direct" (default).
 *   "Direct"              -- Moore-Penrose solve PseudoInverse[m] . b.
 *                            Works on every input family: exact
 *                            (Integer / Rational), symbolic, inexact
 *                            (Real / MPFR), and complex.  This is the
 *                            workhorse method; the LeastSquares ==
 *                            PseudoInverse . b identity is the
 *                            fundamental specification.
 *
 *   "IterativeRefinement" -- residual-correction loop on top of Direct.
 *                            Starting from x = PseudoInverse[m] . b we
 *                            repeatedly compute r = b - m . x,
 *                            dx = PseudoInverse[m] . r, x <- x + dx
 *                            until ||dx||^2 <= tol^2 or a 50-iteration
 *                            cap is hit.  For exact inputs the first
 *                            correction is exactly zero (the pseudoinverse
 *                            identity x = pinv b implies dx = pinv (I - m
 *                            pinv) b = 0), so the loop converges in a
 *                            single pass; for inexact inputs the loop
 *                            drives round-off down to Tolerance.
 *
 *   "LSQR"                -- Paige-Saunders LSQR.  Lanczos
 *                            bidiagonalisation of A combined with a
 *                            Givens-rotation update of the upper
 *                            triangular factor R, exactly the algorithm
 *                            from Paige & Saunders, ACM TOMS 1982.
 *                            Convergence test uses their |phi_bar *
 *                            alpha_{k+1}| estimate of ||A^T r||, scaled
 *                            against the initial gradient ||A^T b||.
 *                            Dispatches by input type: free symbols go
 *                            to Direct (the iteration's stopping test is
 *                            undecidable); exact (Integer / Rational)
 *                            and Complex inputs go to CGLS (Krylov is
 *                            mathematically equivalent and avoids the
 *                            square-root growth in exact arithmetic);
 *                            pure-real inexact inputs run the canonical
 *                            double-precision algorithm.
 *
 *   "Krylov"              -- Conjugate-Gradient-on-Least-Squares
 *                            (Hestenes-Stiefel CG applied to the normal
 *                            equations).  Iterates
 *                                q = A p,  alpha = |s|^2 / |q|^2,
 *                                x <- x + alpha p,  r <- r - alpha q,
 *                                s = A^T r,  beta = |s_new|^2 / |s|^2,
 *                                p <- s_new + beta p
 *                            with x_0 = 0, r_0 = b, s_0 = A^T b, p_0 = s_0.
 *                            Stops on |s|^2 <= tol^2 (or exact zero), a
 *                            null search direction |q|^2 = 0, or a
 *                            2 cols(A) + 10 iteration cap.  Restricted to
 *                            numeric inputs (Integer / Rational / Real /
 *                            Complex with numeric components); symbolic
 *                            inputs fall back to Direct to avoid running
 *                            the loop with an undecidable convergence
 *                            test.  Matrix RHS are solved column by
 *                            column and recombined via Transpose.
 *
 * Tolerance:
 *   Tolerance -> Automatic (default), or a non-negative number.
 *   Forwarded verbatim as the Tolerance option of the underlying
 *   PseudoInverse call so a future singular-value-truncation pass in
 *   PseudoInverse picks it up automatically.  The iterative methods
 *   will also use Tolerance as a convergence threshold once they have
 *   dedicated implementations.
 *
 * Memory ownership follows the standard builtin contract: this file
 * owns `res` on success and frees it; on failure (returning NULL) the
 * caller (the evaluator) retains ownership and the expression remains
 * unevaluated.  Every intermediate PseudoInverse / Dot / Plus / Times
 * call goes through eval_and_free so its argument tree is consumed and
 * its return value is owned by this function.
 */

#include "lstsq.h"
#include "linalg.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "print.h"
#include "sym_names.h"

#include <gmp.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef enum {
    LSTSQ_AUTOMATIC = 0,
    LSTSQ_DIRECT,
    LSTSQ_ITERREFINE,
    LSTSQ_LSQR,
    LSTSQ_KRYLOV,
    LSTSQ_INVALID
} LstsqMethod;

/* Default convergence tolerance when the user requests Tolerance ->
 * Automatic.  Picked an order of magnitude above sqrt(DBL_EPSILON) so
 * inexact iterations converge in a handful of steps without dropping
 * into the noise.  The squared-norm test compares against tol^2. */
#define LSTSQ_DEFAULT_TOL     1.0e-10

/* Hard cap on the IterativeRefinement loop.  Hit either when round-off
 * stalls at a level above Tolerance, or when the user passes a tiny
 * Tolerance unreachable in double precision.  Refinement is unconditionally
 * stable so any extra passes are wasted but not harmful. */
#define LSTSQ_REFINE_MAX_ITER 50

/* LSQR rank-deficiency / breakdown threshold.  When the next bidiagonal
 * alpha or beta drops below this fraction of its initial value the
 * Krylov subspace is exhausted; further iterations would divide by a
 * near-zero scalar and inflate the iterate.  1e-12 is a few orders of
 * magnitude above DBL_EPSILON, which empirically separates real rank
 * deficiency from honest round-off in well-conditioned systems. */
#define LSTSQ_LSQR_BREAKDOWN_EPS 1.0e-12

/* ------------------------------------------------------------------
 * Option parsing.
 *
 * LeastSquares accepts a flexible argument list -- the matrix m and
 * RHS b are required, after which any positional argument is parsed as
 * either a `Method -> ...` rule or a `Tolerance -> ...` rule.  Each
 * option may appear at most once; unknown options leave the call
 * unevaluated so the user sees a clean LeastSquares[...] back, matching
 * the PseudoInverse / LinearSolve grammar.
 * ------------------------------------------------------------------ */

static LstsqMethod parse_method_option(Expr* opt) {
    if (opt->type != EXPR_FUNCTION) return LSTSQ_INVALID;
    if (opt->data.function.head->type != EXPR_SYMBOL) return LSTSQ_INVALID;
    const char* hd = opt->data.function.head->data.symbol;
    if ((hd != SYM_Rule && hd != SYM_RuleDelayed)
        || opt->data.function.arg_count != 2) return LSTSQ_INVALID;
    Expr* lhs = opt->data.function.args[0];
    Expr* rhs = opt->data.function.args[1];
    if (lhs->type != EXPR_SYMBOL || lhs->data.symbol != SYM_Method)
        return LSTSQ_INVALID;

    /* Method -> Automatic (the symbol). */
    if (rhs->type == EXPR_SYMBOL && rhs->data.symbol == SYM_Automatic)
        return LSTSQ_AUTOMATIC;
    if (rhs->type != EXPR_STRING) return LSTSQ_INVALID;
    if (strcmp(rhs->data.string, "Automatic")           == 0) return LSTSQ_AUTOMATIC;
    if (strcmp(rhs->data.string, "Direct")              == 0) return LSTSQ_DIRECT;
    if (strcmp(rhs->data.string, "IterativeRefinement") == 0) return LSTSQ_ITERREFINE;
    if (strcmp(rhs->data.string, "LSQR")                == 0) return LSTSQ_LSQR;
    if (strcmp(rhs->data.string, "Krylov")              == 0) return LSTSQ_KRYLOV;
    return LSTSQ_INVALID;
}

/* Returns true iff `opt` is a Tolerance -> {Automatic | non-negative
 * number} rule.  Out-parameter `*is_automatic` is set to true for the
 * Automatic case (in which `*tol_out` is left as 0.0); a numeric value
 * sets `*is_automatic` to false and writes the value into `*tol_out`.
 *
 * Matches the grammar accepted by PseudoInverse so the same option
 * expression can be forwarded verbatim. */
static bool parse_tolerance_option(Expr* opt, bool* is_automatic_out,
                                   double* tol_out) {
    *is_automatic_out = true;
    *tol_out = 0.0;
    if (opt->type != EXPR_FUNCTION) return false;
    if (opt->data.function.head->type != EXPR_SYMBOL) return false;
    const char* hd = opt->data.function.head->data.symbol;
    if ((hd != SYM_Rule && hd != SYM_RuleDelayed)
        || opt->data.function.arg_count != 2) return false;
    Expr* lhs = opt->data.function.args[0];
    Expr* rhs = opt->data.function.args[1];
    if (lhs->type != EXPR_SYMBOL) return false;
    if (lhs->data.symbol != SYM_Tolerance) return false;

    if (rhs->type == EXPR_SYMBOL && rhs->data.symbol == SYM_Automatic) {
        *is_automatic_out = true;
        return true;
    }
    if (rhs->type == EXPR_INTEGER) {
        if (rhs->data.integer < 0) return false;
        *is_automatic_out = false;
        *tol_out = (double)rhs->data.integer;
        return true;
    }
    if (rhs->type == EXPR_REAL) {
        if (rhs->data.real < 0.0) return false;
        *is_automatic_out = false;
        *tol_out = rhs->data.real;
        return true;
    }
    if (rhs->type == EXPR_FUNCTION
        && rhs->data.function.head->type == EXPR_SYMBOL
        && rhs->data.function.head->data.symbol == SYM_Rational
        && rhs->data.function.arg_count == 2) {
        Expr* num = rhs->data.function.args[0];
        Expr* den = rhs->data.function.args[1];
        if (num->type != EXPR_INTEGER || den->type != EXPR_INTEGER) return false;
        if (num->data.integer < 0) return false;
        if (den->data.integer <= 0) return false;
        *is_automatic_out = false;
        *tol_out = (double)num->data.integer / (double)den->data.integer;
        return true;
    }
    return false;
}

/* ------------------------------------------------------------------
 * Shared numeric helpers used by every iterative method.
 *
 * - compute_norm_squared(v): returns Total[Flatten[v]^2] as an Expr.
 *   Works for both vector (rank-1) and matrix (rank-2) v; for a matrix
 *   it computes the Frobenius squared norm.  Implemented as
 *   Dot[Flatten[v], Flatten[v]] so it is correct for every numeric and
 *   symbolic type the evaluator knows about.
 *
 * - norm_squared_below(nsq, tol_sq): true iff the Expr nsq represents a
 *   non-negative value <= tol_sq.  Numeric atoms compare directly; an
 *   exact zero (Integer 0 or BigInt 0) is treated as converged
 *   unconditionally; anything else (Rational, symbolic) returns false so
 *   the caller keeps iterating until the cap.  This matches the
 *   semantics every iterative method needs: in exact arithmetic the
 *   residual eventually hits Integer 0 and we exit; in inexact
 *   arithmetic we exit when the residual passes below Tolerance.
 * ------------------------------------------------------------------ */
static Expr* compute_norm_squared(Expr* v) {
    Expr* flat = eval_and_free(expr_new_function(
        expr_new_symbol(SYM_Flatten),
        (Expr*[]){expr_copy(v)}, 1));
    if (!flat) return NULL;
    Expr* flat_copy = expr_copy(flat);
    Expr* nsq = eval_and_free(expr_new_function(
        expr_new_symbol("Dot"),
        (Expr*[]){flat, flat_copy}, 2));
    return nsq;
}

static bool norm_squared_below(Expr* nsq, double tol_sq) {
    if (!nsq) return false;
    if (nsq->type == EXPR_INTEGER) return nsq->data.integer == 0;
    if (nsq->type == EXPR_BIGINT)  return mpz_sgn(nsq->data.bigint) == 0;
    if (nsq->type == EXPR_REAL)    return nsq->data.real <= tol_sq;
    return false;
}

/* Numeric-leaf and numeric-tensor predicates.
 *
 * The iterative methods (CGLS / LSQR) have an undecidable convergence
 * test when the residual norm is a symbolic expression -- we cannot tell
 * whether 1/(a^2 + b) "is small".  We therefore restrict these methods
 * to inputs whose leaves are pure numbers (Integer, BigInt, Real,
 * Rational, or Complex of those).  Anything else routes to Direct,
 * which is happy with symbolic arithmetic via PseudoInverse. */
static bool is_numeric_atom(Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) return true;
    if (e->type == EXPR_BIGINT)  return true;
    if (e->type == EXPR_REAL)    return true;
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol;
        if (h == SYM_Rational || h == SYM_Complex) {
            for (size_t i = 0; i < e->data.function.arg_count; i++)
                if (!is_numeric_atom(e->data.function.args[i])) return false;
            return true;
        }
    }
    return false;
}

static bool is_numeric_tensor(Expr* e) {
    if (!e) return false;
    if (is_numeric_atom(e)) return true;
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_List) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (!is_numeric_tensor(e->data.function.args[i])) return false;
        return true;
    }
    return false;
}

/* Predicates used by the LSQR dispatcher.
 * - has_inexact_leaf: any Real value anywhere in the tree.
 * - has_complex_leaf: any Complex[...] subexpression anywhere in the tree.
 * Together they decide whether to run the fast double-precision LSQR
 * (pure real, with at least one Real entry forcing inexact arithmetic)
 * or to fall back to CGLS (exact rational or complex inputs). */
static bool has_inexact_leaf(Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
    if (e->type == EXPR_FUNCTION) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (has_inexact_leaf(e->data.function.args[i])) return true;
    }
    return false;
}
static bool has_complex_leaf(Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_FUNCTION) {
        const char* h = (e->data.function.head->type == EXPR_SYMBOL)
                       ? e->data.function.head->data.symbol : NULL;
        if (h == SYM_Complex) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (has_complex_leaf(e->data.function.args[i])) return true;
    }
    return false;
}

/* Convert a numeric atom (Integer / BigInt / Real / Rational of integers)
 * to a double.  Returns false for anything outside that grammar. */
static bool expr_as_double(Expr* e, double* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer; return true; }
    if (e->type == EXPR_REAL)    { *out = e->data.real;            return true; }
    if (e->type == EXPR_BIGINT)  { *out = mpz_get_d(e->data.bigint); return true; }
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Rational
        && e->data.function.arg_count == 2) {
        double n, d;
        if (!expr_as_double(e->data.function.args[0], &n)) return false;
        if (!expr_as_double(e->data.function.args[1], &d)) return false;
        if (d == 0.0) return false;
        *out = n / d;
        return true;
    }
    return false;
}

/* Row-major extraction of a rank-2 List[List[...], ...] matrix into a
 * heap-allocated double[].  Returns NULL on any non-numeric leaf or
 * non-rectangular shape.  Caller owns the returned buffer. */
static double* extract_matrix_doubles(Expr* m, int* rows_out, int* cols_out) {
    int64_t dims[64];
    if (get_tensor_dims(m, dims) != 2) return NULL;
    int rows = (int)dims[0], cols = (int)dims[1];
    if (rows <= 0 || cols <= 0) return NULL;
    if (m->type != EXPR_FUNCTION) return NULL;
    double* buf = (double*)malloc(sizeof(double) * (size_t)rows * (size_t)cols);
    if (!buf) return NULL;
    for (int i = 0; i < rows; i++) {
        Expr* row = m->data.function.args[i];
        if (row->type != EXPR_FUNCTION
            || (int)row->data.function.arg_count != cols) {
            free(buf); return NULL;
        }
        for (int j = 0; j < cols; j++) {
            if (!expr_as_double(row->data.function.args[j],
                                 &buf[i * cols + j])) {
                free(buf); return NULL;
            }
        }
    }
    *rows_out = rows; *cols_out = cols;
    return buf;
}
static double* extract_vector_doubles(Expr* v, int* len_out) {
    int64_t dims[64];
    if (get_tensor_dims(v, dims) != 1) return NULL;
    int n = (int)dims[0];
    if (n <= 0) return NULL;
    if (v->type != EXPR_FUNCTION) return NULL;
    double* buf = (double*)malloc(sizeof(double) * (size_t)n);
    if (!buf) return NULL;
    for (int i = 0; i < n; i++) {
        if (!expr_as_double(v->data.function.args[i], &buf[i])) {
            free(buf); return NULL;
        }
    }
    *len_out = n;
    return buf;
}
static Expr* doubles_to_vector_expr(const double* x, int n) {
    Expr** args = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    if (!args) return NULL;
    for (int i = 0; i < n; i++) args[i] = expr_new_real(x[i]);
    Expr* v = expr_new_function(expr_new_symbol(SYM_List), args, (size_t)n);
    free(args);
    return v;
}

/* Build {0, 0, ..., 0} of length n.  Used as the CGLS / LSQR starting
 * iterate; starting from the zero vector guarantees the minimum-norm
 * Moore-Penrose solution when m is rank-deficient. */
static Expr* make_zero_vector(int n) {
    Expr** args = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    if (!args) return NULL;
    for (int i = 0; i < n; i++) args[i] = expr_new_integer(0);
    Expr* v = expr_new_function(expr_new_symbol(SYM_List), args, (size_t)n);
    free(args);
    return v;
}

/* ------------------------------------------------------------------
 * Direct solver.
 *
 *   x = PseudoInverse[m, Tolerance -> tol_opt] . b
 *
 * If `tol_opt_or_null` is non-NULL we forward the user's exact Rule
 * expression so the value is preserved (including Tolerance ->
 * Automatic with the explicit symbol).  When NULL, PseudoInverse uses
 * its own default.
 *
 * Failure modes:
 *   - PseudoInverse returns unevaluated (e.g. on a degenerate shape
 *     we missed in our caller).  We detect that the result still has
 *     head PseudoInverse and propagate failure.
 *   - Dot evaluation produces a non-list (defensive); we return NULL.
 *
 * On success the returned Expr* is owned by the caller.
 * ------------------------------------------------------------------ */
static Expr* direct_solve(Expr* m, Expr* b, Expr* tol_opt_or_null) {
    Expr* pinv = NULL;
    if (tol_opt_or_null) {
        pinv = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_PseudoInverse),
            (Expr*[]){expr_copy(m), expr_copy(tol_opt_or_null)}, 2));
    } else {
        pinv = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_PseudoInverse),
            (Expr*[]){expr_copy(m)}, 1));
    }
    if (!pinv) return NULL;
    /* PseudoInverse failed (returned the unevaluated call). */
    if (pinv->type == EXPR_FUNCTION
        && pinv->data.function.head->type == EXPR_SYMBOL
        && pinv->data.function.head->data.symbol == SYM_PseudoInverse) {
        expr_free(pinv);
        return NULL;
    }

    Expr* result = eval_and_free(expr_new_function(
        expr_new_symbol("Dot"),
        (Expr*[]){pinv, expr_copy(b)}, 2));
    /* `pinv` was transferred into the Dot call above; eval_and_free
     * frees the call's input tree.  Ownership of `result` is ours. */
    return result;
}

/* ------------------------------------------------------------------
 * Iterative-refinement solver.
 *
 * Standard fixed-point refinement of an approximate LS solution:
 *
 *   x  <- PseudoInverse[m] . b
 *   loop:
 *       r  <- b - m . x
 *       dx <- PseudoInverse[m] . r
 *       x  <- x + dx
 *       break if ||dx||^2 <= tol^2  (numeric)
 *               or ||dx||^2 == 0    (exact)
 *               or iteration cap hit
 *               or the refinement direct_solve fails
 *
 * For an exact input the Moore-Penrose identity m pinv m = m gives
 *
 *     dx = pinv (b - m pinv b) = pinv b - pinv m pinv b = pinv b - pinv b = 0
 *
 * exactly on the first pass, so the loop exits after one iteration and
 * the result equals "Direct".  For inexact inputs each iteration roughly
 * squares the precision of the residual until round-off saturates, at
 * which point ||dx||^2 stops shrinking and we exit on either the
 * Tolerance test or the iteration cap.  If the inner direct_solve ever
 * fails we silently fall back to the best x computed so far -- never
 * worse than the previous iterate.
 *
 * The squared-norm test uses Total[Flatten[dx]^2] which is correct for
 * both vector b (Frobenius == 2-norm) and matrix b (Frobenius norm,
 * exactly the sum of column 2-norm^2).
 * ------------------------------------------------------------------ */
static Expr* iter_refine_solve(Expr* m, Expr* b, Expr* tol_opt_or_null,
                                bool tol_is_automatic, double tol_val) {
    double tol_sq;
    if (tol_is_automatic) {
        tol_sq = LSTSQ_DEFAULT_TOL * LSTSQ_DEFAULT_TOL;
    } else {
        tol_sq = tol_val * tol_val;
    }

    Expr* x = direct_solve(m, b, tol_opt_or_null);
    if (!x) return NULL;

    for (int k = 0; k < LSTSQ_REFINE_MAX_ITER; k++) {
        /* mx = m . x */
        Expr* mx = eval_and_free(expr_new_function(
            expr_new_symbol("Dot"),
            (Expr*[]){expr_copy(m), expr_copy(x)}, 2));
        /* neg_mx = -1 * mx */
        Expr* neg_mx = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Times),
            (Expr*[]){expr_new_integer(-1), mx}, 2));
        /* r = b - m . x */
        Expr* r = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Plus),
            (Expr*[]){expr_copy(b), neg_mx}, 2));

        /* dx = PseudoInverse[m] . r */
        Expr* dx = direct_solve(m, r, tol_opt_or_null);
        expr_free(r);
        if (!dx) break;  /* keep the current x */

        /* x <- x + dx (consumes x; we still need dx for the norm test) */
        Expr* dx_copy = expr_copy(dx);
        x = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Plus),
            (Expr*[]){x, dx_copy}, 2));

        /* Convergence: ||dx||^2 <= tol^2 (numeric) or exactly zero (exact). */
        Expr* nsq = compute_norm_squared(dx);
        expr_free(dx);
        bool converged = norm_squared_below(nsq, tol_sq);
        expr_free(nsq);
        if (converged) break;
    }
    return x;
}

/* ------------------------------------------------------------------
 * Krylov solver (CGLS / Hestenes-Stiefel on the normal equations).
 *
 * Single-RHS worker.  Starting from x_0 = 0 (so the iterate stays in
 * range(A^T) and converges to the minimum-norm LS solution when A is
 * rank-deficient), the loop maintains
 *
 *     r = b - A x       (residual in R^m)
 *     s = A^T r          (gradient in R^n)
 *     p                  (search direction in R^n)
 *     gamma = |s|^2      (squared gradient norm, monitored for exit)
 *
 * and at each step computes
 *
 *     q     = A p
 *     alpha = gamma / |q|^2
 *     x     <- x + alpha p
 *     r     <- r - alpha q
 *     s     = A^T r
 *     beta  = |s_new|^2 / gamma
 *     gamma <- |s_new|^2
 *     p     <- s_new + beta p
 *
 * Exits when gamma <= tol^2 (numeric) or gamma is exactly Integer 0
 * (exact), or when q^T q is exactly zero (null search direction --
 * happens at the LS solution of a rank-deficient system), or when the
 * iteration cap is hit.  The cap is 2 cols(A) + 10 -- conservatively
 * above the theoretical bound of cols(A) iterations for exact
 * arithmetic, with margin for round-off restarts in inexact arithmetic.
 *
 * Every intermediate Expr is tracked through the function so that the
 * cleanup at the bottom frees exactly the four live vectors r, s, p, and
 * gamma.  x is consumed by the caller's return path.
 * ------------------------------------------------------------------ */
static Expr* cgls_solve_single_rhs(Expr* m, Expr* mt, Expr* b,
                                    int max_iter,
                                    bool tol_is_automatic, double tol_val) {
    double tol_sq;
    if (tol_is_automatic) {
        tol_sq = LSTSQ_DEFAULT_TOL * LSTSQ_DEFAULT_TOL;
    } else {
        tol_sq = tol_val * tol_val;
    }

    int64_t m_dims[64];
    if (get_tensor_dims(m, m_dims) != 2) return NULL;
    int n = (int)m_dims[1];

    Expr* x = make_zero_vector(n);
    if (!x) return NULL;

    /* r = b (because x = 0 implies A x = 0). */
    Expr* r = expr_copy(b);

    /* s = A^T r */
    Expr* s = eval_and_free(expr_new_function(
        expr_new_symbol("Dot"),
        (Expr*[]){expr_copy(mt), expr_copy(r)}, 2));

    Expr* p = expr_copy(s);

    /* gamma = s . s */
    Expr* gamma = eval_and_free(expr_new_function(
        expr_new_symbol("Dot"),
        (Expr*[]){expr_copy(s), expr_copy(s)}, 2));

    for (int k = 0; k < max_iter; k++) {
        if (norm_squared_below(gamma, tol_sq)) break;

        /* q = A p */
        Expr* q = eval_and_free(expr_new_function(
            expr_new_symbol("Dot"),
            (Expr*[]){expr_copy(m), expr_copy(p)}, 2));

        /* qq = q . q */
        Expr* qq = eval_and_free(expr_new_function(
            expr_new_symbol("Dot"),
            (Expr*[]){expr_copy(q), expr_copy(q)}, 2));

        /* Null-direction guard: at the LS solution of a rank-deficient
         * system the next p lies in null(A) and q vanishes exactly. */
        if (qq->type == EXPR_INTEGER && qq->data.integer == 0) {
            expr_free(q); expr_free(qq); break;
        }
        if (qq->type == EXPR_BIGINT && mpz_sgn(qq->data.bigint) == 0) {
            expr_free(q); expr_free(qq); break;
        }
        if (qq->type == EXPR_REAL && qq->data.real == 0.0) {
            expr_free(q); expr_free(qq); break;
        }

        /* alpha = gamma / qq = gamma * Power[qq, -1] */
        Expr* inv_qq = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Power),
            (Expr*[]){qq, expr_new_integer(-1)}, 2));
        Expr* alpha = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Times),
            (Expr*[]){expr_copy(gamma), inv_qq}, 2));

        /* x <- x + alpha p */
        Expr* alpha_p = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Times),
            (Expr*[]){expr_copy(alpha), expr_copy(p)}, 2));
        x = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Plus),
            (Expr*[]){x, alpha_p}, 2));

        /* r <- r - alpha q  (alpha and q consumed here) */
        Expr* alpha_q = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Times),
            (Expr*[]){alpha, q}, 2));
        Expr* neg_alpha_q = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Times),
            (Expr*[]){expr_new_integer(-1), alpha_q}, 2));
        r = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Plus),
            (Expr*[]){r, neg_alpha_q}, 2));

        /* s_new = A^T r */
        Expr* s_new = eval_and_free(expr_new_function(
            expr_new_symbol("Dot"),
            (Expr*[]){expr_copy(mt), expr_copy(r)}, 2));

        /* gamma_new = s_new . s_new */
        Expr* gamma_new = eval_and_free(expr_new_function(
            expr_new_symbol("Dot"),
            (Expr*[]){expr_copy(s_new), expr_copy(s_new)}, 2));

        /* beta = gamma_new / gamma  (gamma consumed here) */
        Expr* inv_gamma = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Power),
            (Expr*[]){gamma, expr_new_integer(-1)}, 2));
        Expr* beta = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Times),
            (Expr*[]){expr_copy(gamma_new), inv_gamma}, 2));

        /* p <- s_new + beta p  (beta and old p consumed here) */
        Expr* beta_p = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Times),
            (Expr*[]){beta, p}, 2));
        p = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Plus),
            (Expr*[]){expr_copy(s_new), beta_p}, 2));

        /* Rotate s and gamma. */
        expr_free(s); s = s_new;
        gamma = gamma_new;
    }

    expr_free(r);
    expr_free(s);
    expr_free(p);
    expr_free(gamma);
    return x;
}

/* Multi-RHS wrapper: when b is a rows x k matrix we solve k single-RHS
 * problems column by column and recombine the answers via Transpose.
 * Single-vector b is forwarded straight to the worker.  m is converted
 * to Transpose[m] once and reused for every column. */
static Expr* cgls_solve(Expr* m, Expr* b,
                         bool tol_is_automatic, double tol_val) {
    int64_t m_dims[64];
    if (get_tensor_dims(m, m_dims) != 2) return NULL;
    int n = (int)m_dims[1];
    int max_iter = 2 * n + 10;

    /* mt = Transpose[m] shared across every RHS column. */
    Expr* mt = eval_and_free(expr_new_function(
        expr_new_symbol(SYM_Transpose),
        (Expr*[]){expr_copy(m)}, 1));
    if (!mt) return NULL;

    int64_t b_dims[64];
    int b_rank = get_tensor_dims(b, b_dims);
    if (b_rank == 1) {
        Expr* x = cgls_solve_single_rhs(m, mt, b, max_iter,
                                          tol_is_automatic, tol_val);
        expr_free(mt);
        return x;
    }
    if (b_rank != 2) { expr_free(mt); return NULL; }

    /* bT = Transpose[b]: rows of bT are columns of b. */
    Expr* bT = eval_and_free(expr_new_function(
        expr_new_symbol(SYM_Transpose),
        (Expr*[]){expr_copy(b)}, 1));
    if (!bT || bT->type != EXPR_FUNCTION
        || bT->data.function.head->type != EXPR_SYMBOL
        || bT->data.function.head->data.symbol != SYM_List) {
        if (bT) expr_free(bT);
        expr_free(mt);
        return NULL;
    }

    int k = (int)bT->data.function.arg_count;
    Expr** xT_rows = (Expr**)malloc(sizeof(Expr*) * (size_t)k);
    if (!xT_rows) { expr_free(bT); expr_free(mt); return NULL; }

    bool failed = false;
    for (int j = 0; j < k; j++) {
        Expr* bcol = bT->data.function.args[j];  /* borrowed */
        Expr* xcol = cgls_solve_single_rhs(m, mt, bcol, max_iter,
                                            tol_is_automatic, tol_val);
        if (!xcol) {
            failed = true;
            for (int i = 0; i < j; i++) expr_free(xT_rows[i]);
            break;
        }
        xT_rows[j] = xcol;
    }
    expr_free(bT);
    expr_free(mt);

    if (failed) { free(xT_rows); return NULL; }

    /* xT is k x n; transpose to get n x k. */
    Expr* xT = expr_new_function(expr_new_symbol(SYM_List), xT_rows, (size_t)k);
    free(xT_rows);  /* shallow free; row pointers transferred to xT */
    Expr* x = eval_and_free(expr_new_function(
        expr_new_symbol(SYM_Transpose),
        (Expr*[]){xT}, 1));
    return x;
}

/* ------------------------------------------------------------------
 * LSQR solver (Paige-Saunders 1982).
 *
 * Lanczos bidiagonalisation
 *
 *     beta_{k+1} u_{k+1} = A v_k     - alpha_k u_k
 *     alpha_{k+1} v_{k+1} = A^T u_{k+1} - beta_{k+1} v_k
 *
 * with initial step beta_1 u_1 = b, alpha_1 v_1 = A^T u_1, combined with
 * a Givens-rotation update of the resulting upper bidiagonal least-
 * squares problem.  Each iteration applies one Givens rotation
 *
 *     rho = sqrt(rho_bar^2 + beta_{k+1}^2),  c = rho_bar/rho, s = beta_{k+1}/rho
 *
 * and updates x and w as
 *
 *     x <- x + (phi/rho) w,  w <- v_{k+1} - (theta/rho) w
 *
 * The Paige-Saunders estimate ||A^T r_k|| ~= |phi_bar * alpha_{k+1}|
 * gives a cheap convergence test on the normal-equations residual
 * (Section 5 of the 1982 paper).  Initial gradient ||A^T b|| equals
 * beta_1 alpha_1, so the relative test is
 *
 *     |phi_bar alpha_{k+1}| <= tol * (beta_1 alpha_1).
 *
 * The cap of 2 cols(A) + 10 iterations matches Krylov.  All arithmetic
 * is plain C double, which is why the LSQR dispatcher only routes here
 * for pure-real inexact inputs -- exact rationals and complex entries
 * are pushed to CGLS in lsqr_solve.
 * ------------------------------------------------------------------ */
static bool lsqr_double(const double* A, int rows, int cols,
                         const double* b, double tol, int max_iter,
                         double* x) {
    double* u = (double*)malloc(sizeof(double) * (size_t)rows);
    double* v = (double*)malloc(sizeof(double) * (size_t)cols);
    double* w = (double*)malloc(sizeof(double) * (size_t)cols);
    double* tmp_r = (double*)malloc(sizeof(double) * (size_t)rows);
    double* tmp_c = (double*)malloc(sizeof(double) * (size_t)cols);
    if (!u || !v || !w || !tmp_r || !tmp_c) {
        free(u); free(v); free(w); free(tmp_r); free(tmp_c);
        return false;
    }

    for (int j = 0; j < cols; j++) x[j] = 0.0;

    /* Initial bidiagonal step: beta_1 u_1 = b, alpha_1 v_1 = A^T u_1. */
    double beta_1 = 0.0;
    for (int i = 0; i < rows; i++) beta_1 += b[i] * b[i];
    beta_1 = sqrt(beta_1);
    if (beta_1 == 0.0) goto cleanup;          /* b = 0  ->  x = 0 */
    for (int i = 0; i < rows; i++) u[i] = b[i] / beta_1;

    /* anorm is a cheap order-of-magnitude estimate of ||A||, used to
     * scale the breakdown thresholds: an alpha or beta that drops below
     * eps * anorm is rank deficiency, not honest progress.  Using the
     * matrix's max absolute entry is robust and avoids the trap of
     * scaling against alpha_init itself -- if A^T b is already zero
     * (the LS optimum at x = 0), alpha_init is machine-epsilon and any
     * threshold derived from it would never fire. */
    double anorm = 0.0;
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++) {
            double a = fabs(A[i * cols + j]);
            if (a > anorm) anorm = a;
        }
    double beta_break  = LSTSQ_LSQR_BREAKDOWN_EPS * anorm;
    double alpha_break = LSTSQ_LSQR_BREAKDOWN_EPS * anorm;

    for (int j = 0; j < cols; j++) {
        double s = 0.0;
        for (int i = 0; i < rows; i++) s += A[i * cols + j] * u[i];
        v[j] = s;
    }
    double alpha = 0.0;
    for (int j = 0; j < cols; j++) alpha += v[j] * v[j];
    alpha = sqrt(alpha);
    /* alpha <= eps * anorm  iff  A^T b is essentially zero -- the LS
     * optimum is x = 0 (any rank-deficient system where b lives in
     * null(A^T) hits this branch).  goto cleanup keeps x = 0. */
    if (alpha <= alpha_break) goto cleanup;
    for (int j = 0; j < cols; j++) v[j] /= alpha;
    for (int j = 0; j < cols; j++) w[j] = v[j];

    double phi_bar    = beta_1;
    double rho_bar    = alpha;
    double grad_init  = beta_1 * alpha;        /* ~= ||A^T b||         */
    double tol_grad   = tol * grad_init;

    for (int iter = 0; iter < max_iter; iter++) {
        /* p = A v - alpha u, beta_new = ||p|| */
        for (int i = 0; i < rows; i++) {
            double s = 0.0;
            for (int j = 0; j < cols; j++) s += A[i * cols + j] * v[j];
            tmp_r[i] = s - alpha * u[i];
        }
        double beta_new = 0.0;
        for (int i = 0; i < rows; i++) beta_new += tmp_r[i] * tmp_r[i];
        beta_new = sqrt(beta_new);

        /* End-of-bidiagonalisation: the next u would divide by ~0.  The
         * iteration k Givens collapses to (c = sign(rho_bar), s = 0):
         *     rho = |rho_bar|, phi/rho = phi_bar / rho_bar.
         * Applying this final x update incorporates the last Krylov
         * direction; without it, x is one step short of the LS answer
         * for full-rank near-converged systems. */
        if (beta_new <= beta_break) {
            if (rho_bar != 0.0) {
                double phi_over_rho = phi_bar / rho_bar;
                for (int j = 0; j < cols; j++) x[j] += phi_over_rho * w[j];
            }
            break;
        }
        for (int i = 0; i < rows; i++) u[i] = tmp_r[i] / beta_new;

        /* q = A^T u - beta_new v, alpha_new = ||q|| */
        for (int j = 0; j < cols; j++) {
            double s = 0.0;
            for (int i = 0; i < rows; i++) s += A[i * cols + j] * u[i];
            tmp_c[j] = s - beta_new * v[j];
        }
        double alpha_new = 0.0;
        for (int j = 0; j < cols; j++) alpha_new += tmp_c[j] * tmp_c[j];
        alpha_new = sqrt(alpha_new);
        bool alpha_breakdown = (alpha_new <= alpha_break);

        /* v <- q / alpha_new (only when alpha_new is non-trivial -- a
         * tiny alpha_new signals rank deficiency, dividing by it would
         * inflate v_{k+1} catastrophically). */
        if (!alpha_breakdown) {
            for (int j = 0; j < cols; j++) v[j] = tmp_c[j] / alpha_new;
        }

        /* Givens rotation on the upper bidiagonal factor.  With
         * alpha_new clamped to zero on breakdown, theta and the next
         * rho_bar vanish and the last x update is purely along the
         * just-computed Krylov direction. */
        double rho = sqrt(rho_bar * rho_bar + beta_new * beta_new);
        double cs  = rho_bar / rho;
        double sn  = beta_new / rho;
        double a_eff = alpha_breakdown ? 0.0 : alpha_new;
        double theta = sn * a_eff;
        double phi = cs * phi_bar;
        rho_bar    = -cs * a_eff;
        phi_bar    =  sn * phi_bar;

        /* Update x with the freshly computed phi/rho contribution. */
        double phi_over_rho   = phi / rho;
        double theta_over_rho = theta / rho;
        for (int j = 0; j < cols; j++) {
            x[j] += phi_over_rho * w[j];
            w[j] = alpha_breakdown
                 ? -theta_over_rho * w[j]
                 :  v[j] - theta_over_rho * w[j];
        }

        if (alpha_breakdown) break;
        alpha = alpha_new;

        /* Paige-Saunders normal-equations convergence test. */
        double est_grad = fabs(phi_bar) * alpha_new;
        if (est_grad <= tol_grad) break;
    }

cleanup:
    free(u); free(v); free(w); free(tmp_r); free(tmp_c);
    return true;
}

/* Wrapper: Expr in, Expr out.  Single-vector b runs LSQR once; matrix b
 * (rows x k) runs LSQR k times (one per column) and recombines via
 * Transpose, mirroring cgls_solve's shape handling. */
static Expr* lsqr_solve_double(Expr* m, Expr* b,
                                bool tol_is_automatic, double tol_val) {
    double tol = tol_is_automatic ? LSTSQ_DEFAULT_TOL : tol_val;

    int rows, cols;
    double* A = extract_matrix_doubles(m, &rows, &cols);
    if (!A) return NULL;
    int max_iter = 2 * cols + 10;

    int64_t b_dims[64];
    int b_rank = get_tensor_dims(b, b_dims);

    if (b_rank == 1) {
        int blen;
        double* bv = extract_vector_doubles(b, &blen);
        if (!bv) { free(A); return NULL; }
        double* xv = (double*)malloc(sizeof(double) * (size_t)cols);
        if (!xv) { free(A); free(bv); return NULL; }
        bool ok = lsqr_double(A, rows, cols, bv, tol, max_iter, xv);
        Expr* x = ok ? doubles_to_vector_expr(xv, cols) : NULL;
        free(A); free(bv); free(xv);
        return x;
    }
    if (b_rank != 2) { free(A); return NULL; }

    int b_cols = (int)b_dims[1];
    if (b_cols <= 0) { free(A); return NULL; }
    if (b->type != EXPR_FUNCTION) { free(A); return NULL; }

    /* Flat row-major buffer of b for column extraction. */
    double* B = (double*)malloc(sizeof(double) * (size_t)rows * (size_t)b_cols);
    if (!B) { free(A); return NULL; }
    for (int i = 0; i < rows; i++) {
        Expr* row = b->data.function.args[i];
        if (row->type != EXPR_FUNCTION
            || (int)row->data.function.arg_count != b_cols) {
            free(A); free(B); return NULL;
        }
        for (int j = 0; j < b_cols; j++) {
            if (!expr_as_double(row->data.function.args[j],
                                 &B[i * b_cols + j])) {
                free(A); free(B); return NULL;
            }
        }
    }

    double* col_b   = (double*)malloc(sizeof(double) * (size_t)rows);
    double* col_x   = (double*)malloc(sizeof(double) * (size_t)cols);
    Expr**  xT_rows = (Expr**)malloc(sizeof(Expr*) * (size_t)b_cols);
    if (!col_b || !col_x || !xT_rows) {
        free(A); free(B); free(col_b); free(col_x); free(xT_rows);
        return NULL;
    }
    bool failed = false;
    for (int j = 0; j < b_cols; j++) {
        for (int i = 0; i < rows; i++) col_b[i] = B[i * b_cols + j];
        if (!lsqr_double(A, rows, cols, col_b, tol, max_iter, col_x)) {
            failed = true;
            for (int k = 0; k < j; k++) expr_free(xT_rows[k]);
            break;
        }
        xT_rows[j] = doubles_to_vector_expr(col_x, cols);
    }
    free(A); free(B); free(col_b); free(col_x);
    if (failed) { free(xT_rows); return NULL; }

    Expr* xT = expr_new_function(expr_new_symbol(SYM_List),
                                  xT_rows, (size_t)b_cols);
    free(xT_rows);
    Expr* x = eval_and_free(expr_new_function(
        expr_new_symbol(SYM_Transpose),
        (Expr*[]){xT}, 1));
    return x;
}

/* Public LSQR dispatcher.  Routes by input grammar:
 *   - non-numeric (free symbols)  -> Direct
 *   - exact rational or complex   -> CGLS (equivalent without Sqrt)
 *   - pure-real with any Real     -> double-precision LSQR
 */
static Expr* lsqr_solve(Expr* m, Expr* b, Expr* tol_opt_or_null,
                         bool tol_is_automatic, double tol_val) {
    if (!is_numeric_tensor(m) || !is_numeric_tensor(b)) {
        return direct_solve(m, b, tol_opt_or_null);
    }
    bool complex_in = has_complex_leaf(m) || has_complex_leaf(b);
    bool inexact_in = has_inexact_leaf(m) || has_inexact_leaf(b);
    if (complex_in || !inexact_in) {
        return cgls_solve(m, b, tol_is_automatic, tol_val);
    }
    return lsqr_solve_double(m, b, tol_is_automatic, tol_val);
}

/* ------------------------------------------------------------------
 * Public dispatcher.
 *
 * The user-visible builtin.  After parsing options and validating the
 * shape of `m` and `b` we dispatch to one of the workers above.  LSQR
 * and Krylov currently dispatch to Direct (see the file header).
 * ------------------------------------------------------------------ */
Expr* builtin_leastsquares(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;

    Expr* m = res->data.function.args[0];
    Expr* b = res->data.function.args[1];

    LstsqMethod method   = LSTSQ_AUTOMATIC;
    bool   method_seen   = false;
    Expr*  tol_opt       = NULL;  /* forwarded verbatim to PseudoInverse */
    bool   tolerance_seen   = false;
    bool   tol_is_automatic = true;
    double tol_val          = 0.0;

    for (size_t i = 2; i < argc; i++) {
        Expr* opt = res->data.function.args[i];

        LstsqMethod m_try = parse_method_option(opt);
        if (m_try != LSTSQ_INVALID) {
            if (method_seen) return NULL;  /* duplicate Method -> ... */
            method = m_try;
            method_seen = true;
            continue;
        }

        bool   is_auto = true;
        double tv      = 0.0;
        if (parse_tolerance_option(opt, &is_auto, &tv)) {
            if (tolerance_seen) return NULL;  /* duplicate Tolerance -> ... */
            tol_opt          = opt;
            tol_is_automatic = is_auto;
            tol_val          = tv;
            tolerance_seen   = true;
            continue;
        }

        /* Unknown option -> stay unevaluated. */
        return NULL;
    }

    /* Shape validation -- shared by every method. */
    int64_t m_dims[64];
    int m_rank = get_tensor_dims(m, m_dims);
    if (m_rank != 2 || m_dims[0] == 0 || m_dims[1] == 0) {
        char* s = expr_to_string(m);
        fprintf(stderr,
                "LeastSquares::matrix: Argument %s at position 1 is "
                "not a non-empty rectangular matrix.\n", s);
        free(s);
        return NULL;
    }
    int rows = (int)m_dims[0];

    int64_t b_dims[64];
    int b_rank = get_tensor_dims(b, b_dims);
    if (b_rank != 1 && b_rank != 2) {
        char* s = expr_to_string(b);
        fprintf(stderr,
                "LeastSquares::lvec: %s is neither a vector nor a matrix.\n",
                s);
        free(s);
        return NULL;
    }
    if (b_dims[0] != rows) {
        char* ms = expr_to_string(m);
        char* bs = expr_to_string(b);
        fprintf(stderr,
                "LeastSquares::lvec1: Coefficient matrix and target "
                "vector %s . x == %s do not have the same dimensions.\n",
                ms, bs);
        free(ms);
        free(bs);
        return NULL;
    }

    /* For a rank-2 b, the inner-list length must be non-zero too --
     * a 3 x 0 b makes no sense as a multi-RHS and slips through the
     * leading-dim check above. */
    if (b_rank == 2 && b_dims[1] == 0) {
        char* s = expr_to_string(b);
        fprintf(stderr,
                "LeastSquares::lvec: %s is neither a vector nor a matrix.\n",
                s);
        free(s);
        return NULL;
    }

    Expr* result = NULL;
    switch (method) {
        case LSTSQ_AUTOMATIC:
        case LSTSQ_DIRECT:
            result = direct_solve(m, b, tol_opt);
            break;
        case LSTSQ_LSQR:
            result = lsqr_solve(m, b, tol_opt, tol_is_automatic, tol_val);
            break;
        case LSTSQ_KRYLOV:
            /* CGLS only converges for numeric inputs (its termination
             * test compares a residual norm against a tolerance).  Free
             * symbols make the test undecidable, so we fall back to
             * Direct -- the answer is the same. */
            if (is_numeric_tensor(m) && is_numeric_tensor(b)) {
                result = cgls_solve(m, b, tol_is_automatic, tol_val);
            } else {
                result = direct_solve(m, b, tol_opt);
            }
            break;
        case LSTSQ_ITERREFINE:
            result = iter_refine_solve(m, b, tol_opt,
                                        tol_is_automatic, tol_val);
            break;
        case LSTSQ_INVALID:
            return NULL;  /* unreachable: parse_method_option never sets this. */
    }
    return result;
}

/* ------------------------------------------------------------------
 * Registration.
 * ------------------------------------------------------------------ */
void matlstsq_init(void) {
    symtab_add_builtin("LeastSquares", builtin_leastsquares);
    symtab_get_def("LeastSquares")->attributes |= ATTR_PROTECTED;
}
