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
 *   "IterativeRefinement" -- numeric refinement of "Direct".  We
 *                            compute x = PseudoInverse[m] . b, the
 *                            residual r = b - m . x, the correction
 *                            dx = PseudoInverse[m] . r, and return
 *                            x + dx.  For exact inputs the correction
 *                            is exactly zero and the result equals
 *                            "Direct"; for inexact inputs the second
 *                            pass reduces round-off in the residual.
 *
 *   "LSQR"                -- Paige-Saunders LSQR iterative method.
 *                            Mathematica advertises this for sparse
 *                            machine-number matrices.  In this CAS the
 *                            method name is accepted by the parser so
 *                            user code can be written against the final
 *                            API; the current dispatch goes through
 *                            "Direct".  A dedicated implementation is a
 *                            future extension.
 *
 *   "Krylov"              -- iterative Krylov method (CGNR on the
 *                            normal equations).  As with "LSQR", the
 *                            current dispatch goes through "Direct".
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

#include "matlstsq.h"
#include "linalg.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "print.h"
#include "sym_names.h"

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
            expr_new_symbol("PseudoInverse"),
            (Expr*[]){expr_copy(m), expr_copy(tol_opt_or_null)}, 2));
    } else {
        pinv = eval_and_free(expr_new_function(
            expr_new_symbol("PseudoInverse"),
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
 *   x  <- PseudoInverse[m] . b
 *   r  <- b - m . x
 *   dx <- PseudoInverse[m] . r
 *   x  <- x + dx
 *
 * One refinement pass.  Reduces round-off in numerically tight numeric
 * systems; for exact inputs the residual is exactly zero so this
 * coincides with the Direct answer.  If the refinement step fails for
 * any reason, the initial Direct x is returned -- never worse than
 * Direct, often slightly better.
 * ------------------------------------------------------------------ */
static Expr* iter_refine_solve(Expr* m, Expr* b, Expr* tol_opt_or_null) {
    Expr* x = direct_solve(m, b, tol_opt_or_null);
    if (!x) return NULL;

    /* mx = m . x */
    Expr* mx = eval_and_free(expr_new_function(
        expr_new_symbol("Dot"),
        (Expr*[]){expr_copy(m), expr_copy(x)}, 2));

    /* neg_mx = -mx (a fresh copy, so mx is consumed). */
    Expr* neg_mx = eval_and_free(expr_new_function(
        expr_new_symbol("Times"),
        (Expr*[]){expr_new_integer(-1), mx}, 2));

    /* r = b + neg_mx = b - m.x */
    Expr* r = eval_and_free(expr_new_function(
        expr_new_symbol("Plus"),
        (Expr*[]){expr_copy(b), neg_mx}, 2));

    Expr* dx = direct_solve(m, r, tol_opt_or_null);
    expr_free(r);
    if (!dx) {
        return x;  /* refinement failed; fall back to the Direct answer */
    }

    /* refined = x + dx (x and dx consumed by the Plus). */
    Expr* refined = eval_and_free(expr_new_function(
        expr_new_symbol("Plus"),
        (Expr*[]){x, dx}, 2));
    return refined;
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
    bool   tolerance_seen = false;

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
            (void)tv; /* value used by PseudoInverse via the option Rule */
            if (tolerance_seen) return NULL;  /* duplicate Tolerance -> ... */
            tol_opt = opt;
            tolerance_seen = true;
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
        case LSTSQ_LSQR:    /* dispatched to Direct; see header */
        case LSTSQ_KRYLOV:  /* dispatched to Direct; see header */
            result = direct_solve(m, b, tol_opt);
            break;
        case LSTSQ_ITERREFINE:
            result = iter_refine_solve(m, b, tol_opt);
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
