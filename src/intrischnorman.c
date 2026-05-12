/* intrischnorman.c — parallel Risch / Risch-Norman heuristic
 * integrator.  Implements Bronstein's "Poor Man's Integrator"
 * (pmint, 2004).  Reference: parallel_risch/pmint.maple (99 lines)
 * and parallel_risch/bronstein 2004 - parallel integration.pdf.
 *
 * Phases of RISCH_NORMAN_PLAN.md:
 *   Phase 1: skeleton + dispatcher hook (always returns NULL).
 *   Phase 2: convert_to_tan + indet collection + substitution maps.
 *   Phase 3: vector field + splitFactor + deflation + monomials.
 *   Phase 4: candidate ansatz + linear-system extract + RowReduce.
 *   Phase 5: log-candidate sum + getSpecial + K=I retry.
 *   Phase 6: dispatcher polish + post-hoc verification.
 *
 * Memory contract: every public BuiltinFunc follows the picocas
 * convention — the caller (evaluator) owns `res`; the function
 * returns a freshly-allocated Expr* on success or NULL on failure.
 * Internal static helpers take borrowed input and return owned
 * output unless explicitly documented otherwise.
 *
 * Performance caps (enforced once Phase 4 lands):
 *   PMINT_MAX_MONOMIALS  — hard cap on the size of the candidate
 *                          monomial set; exceeding this returns NULL.
 *   PMINT_WALL_CLOCK_SEC — wall-clock budget per integrand.
 */

/* Expose POSIX symbols (sigaction, sigsetjmp/siglongjmp, sigjmp_buf,
 * setitimer) under -std=c99 on glibc; macOS exposes them implicitly,
 * which masks the cross-platform bug.  Must precede every header. */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "intrischnorman.h"

#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"
#include "attr.h"
#include "internal.h"
#include "poly.h"
#include "sym_names.h"
#include "sym_intern.h"

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <setjmp.h>
#include <signal.h>
#include <sys/time.h>

/* ------------------------------------------------------------------ */
/* Performance caps.                                                    */
/* ------------------------------------------------------------------ */

#define PMINT_MAX_MONOMIALS   5000
#define PMINT_WALL_CLOCK_SEC  30
#define PMINT_MAX_REWRITE_DEPTH 32      /* convert_to_tan giveup depth */
#define PMINT_MAX_INDETS      32        /* hard cap on indet set size  */

#ifndef PMINT_TRACE
#define PMINT_TRACE 0
#endif
#define PMTRACE(...) do { if (PMINT_TRACE) fprintf(stderr, "[pmint] " __VA_ARGS__); } while (0)

/* ------------------------------------------------------------------ */
/* Small utilities.                                                     */
/* ------------------------------------------------------------------ */

/* Build an integer Expr and return owned. */
static Expr* mk_int(int64_t v) {
    return expr_new_integer(v);
}

/* Build Plus[a, b] without evaluating (just structural). */
static Expr* mk_plus2(Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    return expr_new_function(expr_new_symbol(SYM_Plus), args, 2);
}

/* Build Times[a, b] without evaluating. */
static Expr* mk_times2(Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    return expr_new_function(expr_new_symbol(SYM_Times), args, 2);
}

/* Build Power[a, b] without evaluating. */
static Expr* mk_pow(Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    return expr_new_function(expr_new_symbol(SYM_Power), args, 2);
}

/* Build Times[a, b, c] for "a/b * c". */
static Expr* mk_div(Expr* a, Expr* b) {
    /* a / b = a * b^(-1) */
    return mk_times2(a, mk_pow(b, mk_int(-1)));
}

/* Build Function[head, arg] for unary heads (Sin, Cos, ...). */
static Expr* mk_unary(const char* head, Expr* arg) {
    Expr* args[1] = { arg };
    return expr_new_function(expr_new_symbol(head), args, 1);
}

/* Build Function[head, arg1, arg2] for binary heads (D, ...). */
static Expr* mk_binary(const char* head, Expr* a1, Expr* a2) {
    Expr* args[2] = { a1, a2 };
    return expr_new_function(expr_new_symbol(head), args, 2);
}

/* FreeQ[expr, sym]: returns true iff `sym` (as a symbol) does not
 * appear anywhere in `expr`. */
static bool expr_free_of_symbol(const Expr* expr, const char* sym_name) {
    if (!expr) return true;
    if (expr->type == EXPR_SYMBOL) return expr->data.symbol != sym_name;
    if (expr->type != EXPR_FUNCTION) return true;
    if (!expr_free_of_symbol(expr->data.function.head, sym_name)) return false;
    for (size_t i = 0; i < expr->data.function.arg_count; i++) {
        if (!expr_free_of_symbol(expr->data.function.args[i], sym_name)) return false;
    }
    return true;
}

/* Compute D[f, x] by building the Expr and evaluating.  Owned result. */
static Expr* call_d(Expr* f, Expr* x) {
    Expr* call = mk_binary("D", expr_copy(f), expr_copy(x));
    return eval_and_free(call);
}

/* True iff `e` evaluates structurally to the integer 0. */
static bool is_int_zero(const Expr* e) {
    return e && e->type == EXPR_INTEGER && e->data.integer == 0;
}

/* ------------------------------------------------------------------ */
/* Phase 2 — convert_to_tan: half-angle Tan / Tanh rewrite.            */
/*                                                                      */
/* pmint.maple line 5: `ff := eval(convert(f, tan))`.  Mirrors Maple's */
/* convert/tan: Sin / Cos / Sec / Csc rewritten to rational functions  */
/* of Tan[u/2]; Sinh / Cosh / Sech / Csch rewritten via Tanh[u/2].     */
/* Tan / Cot / Tanh / Coth are left alone (Cot ↦ 1/Tan, Coth ↦ 1/Tanh) */
/* — the Tan / Tanh atoms become their own field generators.           */
/*                                                                      */
/* Subexpressions free of x are left alone (FreeQ short-circuit).      */
/* Depth-capped at PMINT_MAX_REWRITE_DEPTH; returns NULL on giveup so  */
/* the top-level builtin can bubble back unevaluated.                  */
/* ------------------------------------------------------------------ */

static Expr* convert_to_tan_rec(Expr* f, const char* xname, int depth) {
    if (!f) return NULL;
    if (depth > PMINT_MAX_REWRITE_DEPTH) return NULL;

    /* Leaf or free of x: copy verbatim. */
    if (f->type != EXPR_FUNCTION || expr_free_of_symbol(f, xname)) {
        return expr_copy(f);
    }

    /* Recurse into the arguments first.  The head we leave alone since
     * pmint only rewrites known trig / hyperbolic heads. */
    size_t n = f->data.function.arg_count;
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    if (!new_args) return NULL;
    for (size_t i = 0; i < n; i++) {
        new_args[i] = convert_to_tan_rec(f->data.function.args[i], xname, depth + 1);
        if (!new_args[i]) {
            for (size_t j = 0; j < i; j++) expr_free(new_args[j]);
            free(new_args);
            return NULL;
        }
    }

    /* Identify the head as a trig / hyperbolic we know how to rewrite.
     * The argument u is new_args[0]. */
    Expr* head = f->data.function.head;
    Expr* result = NULL;

    if (n == 1 && head && head->type == EXPR_SYMBOL) {
        const char* hs = head->data.symbol;
        Expr* u = new_args[0];

        if (hs == SYM_Sin) {
            /* Sin[u] = 2 T / (1 + T^2),  T = Tan[u/2]. */
            Expr* T = mk_unary(SYM_Tan, mk_div(expr_copy(u), mk_int(2)));
            Expr* num = mk_times2(mk_int(2), expr_copy(T));
            Expr* den = mk_plus2(mk_int(1), mk_pow(T, mk_int(2)));
            result = mk_div(num, den);
        } else if (hs == SYM_Cos) {
            /* Cos[u] = (1 - T^2) / (1 + T^2),  T = Tan[u/2]. */
            Expr* T = mk_unary(SYM_Tan, mk_div(expr_copy(u), mk_int(2)));
            Expr* Tsq = mk_pow(expr_copy(T), mk_int(2));
            Expr* num = mk_plus2(mk_int(1), mk_times2(mk_int(-1), Tsq));
            Expr* den = mk_plus2(mk_int(1), mk_pow(T, mk_int(2)));
            result = mk_div(num, den);
        } else if (hs == SYM_Sec) {
            /* Sec[u] = (1 + T^2) / (1 - T^2),  T = Tan[u/2]. */
            Expr* T = mk_unary(SYM_Tan, mk_div(expr_copy(u), mk_int(2)));
            Expr* Tsq_n = mk_pow(expr_copy(T), mk_int(2));
            Expr* num = mk_plus2(mk_int(1), Tsq_n);
            Expr* den = mk_plus2(mk_int(1), mk_times2(mk_int(-1), mk_pow(T, mk_int(2))));
            result = mk_div(num, den);
        } else if (hs == SYM_Csc) {
            /* Csc[u] = (1 + T^2) / (2 T),  T = Tan[u/2]. */
            Expr* T = mk_unary(SYM_Tan, mk_div(expr_copy(u), mk_int(2)));
            Expr* Tsq = mk_pow(expr_copy(T), mk_int(2));
            Expr* num = mk_plus2(mk_int(1), Tsq);
            Expr* den = mk_times2(mk_int(2), T);
            result = mk_div(num, den);
        } else if (hs == SYM_Tan) {
            /* Tan[u] = 2 T / (1 - T^2),  T = Tan[u/2].
             * Full Weierstrass — keeps Tan[u/2] as the unique trig
             * generator and avoids the algebraic dependence between
             * Tan[u] and the Sin/Cos-derived Tan[u/2]. */
            Expr* T = mk_unary(SYM_Tan, mk_div(expr_copy(u), mk_int(2)));
            Expr* num = mk_times2(mk_int(2), expr_copy(T));
            Expr* den = mk_plus2(mk_int(1), mk_times2(mk_int(-1), mk_pow(T, mk_int(2))));
            result = mk_div(num, den);
        } else if (hs == SYM_Cot) {
            /* Cot[u] = (1 - T^2) / (2 T),  T = Tan[u/2]. */
            Expr* T = mk_unary(SYM_Tan, mk_div(expr_copy(u), mk_int(2)));
            Expr* Tsq = mk_pow(expr_copy(T), mk_int(2));
            Expr* num = mk_plus2(mk_int(1), mk_times2(mk_int(-1), Tsq));
            Expr* den = mk_times2(mk_int(2), T);
            result = mk_div(num, den);
        } else if (hs == SYM_Sinh) {
            /* Sinh[u] = 2 T / (1 - T^2),  T = Tanh[u/2]. */
            Expr* T = mk_unary(SYM_Tanh, mk_div(expr_copy(u), mk_int(2)));
            Expr* num = mk_times2(mk_int(2), expr_copy(T));
            Expr* den = mk_plus2(mk_int(1), mk_times2(mk_int(-1), mk_pow(T, mk_int(2))));
            result = mk_div(num, den);
        } else if (hs == SYM_Cosh) {
            /* Cosh[u] = (1 + T^2) / (1 - T^2),  T = Tanh[u/2]. */
            Expr* T = mk_unary(SYM_Tanh, mk_div(expr_copy(u), mk_int(2)));
            Expr* Tsq = mk_pow(expr_copy(T), mk_int(2));
            Expr* num = mk_plus2(mk_int(1), Tsq);
            Expr* den = mk_plus2(mk_int(1), mk_times2(mk_int(-1), mk_pow(T, mk_int(2))));
            result = mk_div(num, den);
        } else if (hs == SYM_Sech) {
            /* Sech[u] = (1 - T^2) / (1 + T^2),  T = Tanh[u/2]. */
            Expr* T = mk_unary(SYM_Tanh, mk_div(expr_copy(u), mk_int(2)));
            Expr* Tsq = mk_pow(expr_copy(T), mk_int(2));
            Expr* num = mk_plus2(mk_int(1), mk_times2(mk_int(-1), Tsq));
            Expr* den = mk_plus2(mk_int(1), mk_pow(T, mk_int(2)));
            result = mk_div(num, den);
        } else if (hs == SYM_Csch) {
            /* Csch[u] = (1 - T^2) / (2 T),  T = Tanh[u/2]. */
            Expr* T = mk_unary(SYM_Tanh, mk_div(expr_copy(u), mk_int(2)));
            Expr* Tsq = mk_pow(expr_copy(T), mk_int(2));
            Expr* num = mk_plus2(mk_int(1), mk_times2(mk_int(-1), Tsq));
            Expr* den = mk_times2(mk_int(2), T);
            result = mk_div(num, den);
        } else if (hs == SYM_Tanh) {
            /* Tanh[u] = 2 T / (1 + T^2),  T = Tanh[u/2]. */
            Expr* T = mk_unary(SYM_Tanh, mk_div(expr_copy(u), mk_int(2)));
            Expr* num = mk_times2(mk_int(2), expr_copy(T));
            Expr* den = mk_plus2(mk_int(1), mk_pow(T, mk_int(2)));
            result = mk_div(num, den);
        } else if (hs == SYM_Coth) {
            /* Coth[u] = (1 + T^2) / (2 T),  T = Tanh[u/2]. */
            Expr* T = mk_unary(SYM_Tanh, mk_div(expr_copy(u), mk_int(2)));
            Expr* Tsq = mk_pow(expr_copy(T), mk_int(2));
            Expr* num = mk_plus2(mk_int(1), Tsq);
            Expr* den = mk_times2(mk_int(2), T);
            result = mk_div(num, den);
        }
    }

    if (result) {
        /* We have replacements; the original new_args contents were
         * already consumed (we expr_copy'd from them inside the
         * formulas).  Free the new_args wrappers. */
        for (size_t i = 0; i < n; i++) expr_free(new_args[i]);
        free(new_args);
        return result;
    }

    /* Not a head we rewrite: rebuild the function with new args. */
    Expr* newf = expr_new_function(expr_copy(head), new_args, n);
    free(new_args);
    return newf;
}

/* Pythagorean rewrite: walks the tree and replaces:
 *   Sec[u]^2 -> 1 + Tan[u]^2
 *   Csc[u]^2 -> 1 + Cot[u]^2
 *   Sech[u]^2 -> 1 - Tanh[u]^2
 *   Csch[u]^2 -> Coth[u]^2 - 1
 *
 * Used inside the closure walk in collect_indets_closed: derivatives
 * of Tan[u]/Tanh[u] produce Sec[u]^2/Sech[u]^2 which we must rewrite
 * back into the Tan / Tanh field generators so the closure doesn't
 * accidentally collect Sec / Csc / Sech / Csch as separate (and
 * algebraically dependent) atoms.
 *
 * Does NOT recursively simplify after rewrite — returns a structural
 * substitution.  Caller may evaluate if appropriate. */
static Expr* pythagorean_rewrite(Expr* f) {
    if (!f) return NULL;
    if (f->type != EXPR_FUNCTION) return expr_copy(f);

    Expr* h = f->data.function.head;
    size_t n = f->data.function.arg_count;

    /* Detect Power[Sec[u], 2] etc. */
    if (h && h->type == EXPR_SYMBOL && h->data.symbol == SYM_Power
        && n == 2
        && f->data.function.args[0]
        && f->data.function.args[0]->type == EXPR_FUNCTION
        && f->data.function.args[0]->data.function.arg_count == 1
        && f->data.function.args[1]
        && f->data.function.args[1]->type == EXPR_INTEGER
        && f->data.function.args[1]->data.integer == 2) {
        Expr* base = f->data.function.args[0];
        Expr* bh = base->data.function.head;
        Expr* u = base->data.function.args[0];
        if (bh && bh->type == EXPR_SYMBOL) {
            const char* bhs = bh->data.symbol;
            Expr* urew = pythagorean_rewrite(u);
            if (bhs == SYM_Sec) {
                /* 1 + Tan[u]^2 */
                return mk_plus2(mk_int(1),
                                mk_pow(mk_unary(SYM_Tan, urew), mk_int(2)));
            } else if (bhs == SYM_Csc) {
                /* 1 + Cot[u]^2 */
                return mk_plus2(mk_int(1),
                                mk_pow(mk_unary(SYM_Cot, urew), mk_int(2)));
            } else if (bhs == SYM_Sech) {
                /* 1 - Tanh[u]^2 */
                return mk_plus2(mk_int(1),
                                mk_times2(mk_int(-1),
                                          mk_pow(mk_unary(SYM_Tanh, urew), mk_int(2))));
            } else if (bhs == SYM_Csch) {
                /* Coth[u]^2 - 1 */
                return mk_plus2(mk_int(-1),
                                mk_pow(mk_unary(SYM_Coth, urew), mk_int(2)));
            }
            expr_free(urew);
        }
    }

    /* Generic recursion. */
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    if (!new_args) return NULL;
    for (size_t i = 0; i < n; i++) {
        new_args[i] = pythagorean_rewrite(f->data.function.args[i]);
        if (!new_args[i]) {
            for (size_t j = 0; j < i; j++) expr_free(new_args[j]);
            free(new_args);
            return NULL;
        }
    }
    Expr* result = expr_new_function(expr_copy(h), new_args, n);
    free(new_args);
    return result;
}

/* Decot: post-eval walk that rewrites Cot[u] back to Tan[u]^(-1)
 * (and Coth[u] to Tanh[u]^(-1)).  Picocas's evaluator collapses
 * Power[Tan[u], -1] to Cot[u], which would break pmint's field-
 * generator invariant (the collected indet set should contain
 * Tan[u] as a single atom; Cot[u] is algebraically dependent and
 * must not become a second generator).
 *
 * After this rewrite the result is structurally rational in
 * Tan[u/2] / Tanh[u/2] alone.  We do NOT re-evaluate — the inverse
 * Power form must survive into pmint's substitution stage. */
static Expr* decot_rec(Expr* f) {
    if (!f) return NULL;

    if (f->type == EXPR_FUNCTION) {
        Expr* h = f->data.function.head;
        if (h && h->type == EXPR_SYMBOL
            && f->data.function.arg_count == 1) {
            const char* hs = h->data.symbol;
            if (hs == SYM_Cot) {
                Expr* inner = decot_rec(f->data.function.args[0]);
                return mk_pow(mk_unary(SYM_Tan, inner), mk_int(-1));
            }
            if (hs == SYM_Coth) {
                Expr* inner = decot_rec(f->data.function.args[0]);
                return mk_pow(mk_unary(SYM_Tanh, inner), mk_int(-1));
            }
            if (hs == SYM_Sec) {
                Expr* inner = decot_rec(f->data.function.args[0]);
                return mk_pow(mk_unary(SYM_Cos, inner), mk_int(-1));
            }
            if (hs == SYM_Csc) {
                Expr* inner = decot_rec(f->data.function.args[0]);
                return mk_pow(mk_unary(SYM_Sin, inner), mk_int(-1));
            }
            if (hs == SYM_Sech) {
                Expr* inner = decot_rec(f->data.function.args[0]);
                return mk_pow(mk_unary(SYM_Cosh, inner), mk_int(-1));
            }
            if (hs == SYM_Csch) {
                Expr* inner = decot_rec(f->data.function.args[0]);
                return mk_pow(mk_unary(SYM_Sinh, inner), mk_int(-1));
            }
        }

        size_t n = f->data.function.arg_count;
        Expr** new_args = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
        if (!new_args) return NULL;
        for (size_t i = 0; i < n; i++) {
            new_args[i] = decot_rec(f->data.function.args[i]);
            if (!new_args[i]) {
                for (size_t j = 0; j < i; j++) expr_free(new_args[j]);
                free(new_args);
                return NULL;
            }
        }
        Expr* result = expr_new_function(expr_copy(h), new_args, n);
        free(new_args);
        return result;
    }

    return expr_copy(f);
}

/* convert_sincos_to_tan_rec — partial Weierstrass: rewrite only
 * Sin[u] / Cos[u] / Sec[u] / Csc[u] and their hyperbolic siblings.
 * Tan / Cot / Tanh / Coth are left alone.  Used by the post-hoc
 * verifier so that pmint's Tan-form output and the original
 * Sin/Cos-form integrand can be unified by Cancel[Together[...]]
 * without nesting Tan[x/4] into Tan[x/2]'s. */
static Expr* convert_sincos_to_tan_rec(Expr* f, const char* xname, int depth) {
    if (!f) return NULL;
    if (depth > PMINT_MAX_REWRITE_DEPTH) return NULL;
    if (f->type != EXPR_FUNCTION || expr_free_of_symbol(f, xname)) {
        return expr_copy(f);
    }

    size_t n = f->data.function.arg_count;
    Expr** new_args = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    if (!new_args) return NULL;
    for (size_t i = 0; i < n; i++) {
        new_args[i] = convert_sincos_to_tan_rec(f->data.function.args[i], xname, depth + 1);
        if (!new_args[i]) {
            for (size_t j = 0; j < i; j++) expr_free(new_args[j]);
            free(new_args);
            return NULL;
        }
    }

    Expr* head = f->data.function.head;
    Expr* result = NULL;
    if (n == 1 && head && head->type == EXPR_SYMBOL) {
        const char* hs = head->data.symbol;
        Expr* u = new_args[0];
        if (hs == SYM_Sin) {
            Expr* T = mk_unary(SYM_Tan, mk_div(expr_copy(u), mk_int(2)));
            result = mk_div(mk_times2(mk_int(2), expr_copy(T)),
                             mk_plus2(mk_int(1), mk_pow(T, mk_int(2))));
        } else if (hs == SYM_Cos) {
            Expr* T = mk_unary(SYM_Tan, mk_div(expr_copy(u), mk_int(2)));
            Expr* Tsq = mk_pow(expr_copy(T), mk_int(2));
            Expr* num = mk_plus2(mk_int(1), mk_times2(mk_int(-1), Tsq));
            Expr* den = mk_plus2(mk_int(1), mk_pow(T, mk_int(2)));
            result = mk_div(num, den);
        } else if (hs == SYM_Sec) {
            Expr* T = mk_unary(SYM_Tan, mk_div(expr_copy(u), mk_int(2)));
            Expr* Tsq_n = mk_pow(expr_copy(T), mk_int(2));
            result = mk_div(mk_plus2(mk_int(1), Tsq_n),
                             mk_plus2(mk_int(1),
                                      mk_times2(mk_int(-1), mk_pow(T, mk_int(2)))));
        } else if (hs == SYM_Csc) {
            Expr* T = mk_unary(SYM_Tan, mk_div(expr_copy(u), mk_int(2)));
            result = mk_div(mk_plus2(mk_int(1), mk_pow(expr_copy(T), mk_int(2))),
                             mk_times2(mk_int(2), T));
        } else if (hs == SYM_Sinh) {
            Expr* T = mk_unary(SYM_Tanh, mk_div(expr_copy(u), mk_int(2)));
            result = mk_div(mk_times2(mk_int(2), expr_copy(T)),
                             mk_plus2(mk_int(1),
                                      mk_times2(mk_int(-1), mk_pow(T, mk_int(2)))));
        } else if (hs == SYM_Cosh) {
            Expr* T = mk_unary(SYM_Tanh, mk_div(expr_copy(u), mk_int(2)));
            Expr* Tsq = mk_pow(expr_copy(T), mk_int(2));
            result = mk_div(mk_plus2(mk_int(1), Tsq),
                             mk_plus2(mk_int(1),
                                      mk_times2(mk_int(-1), mk_pow(T, mk_int(2)))));
        } else if (hs == SYM_Sech) {
            Expr* T = mk_unary(SYM_Tanh, mk_div(expr_copy(u), mk_int(2)));
            Expr* Tsq = mk_pow(expr_copy(T), mk_int(2));
            result = mk_div(mk_plus2(mk_int(1), mk_times2(mk_int(-1), Tsq)),
                             mk_plus2(mk_int(1), mk_pow(T, mk_int(2))));
        } else if (hs == SYM_Csch) {
            Expr* T = mk_unary(SYM_Tanh, mk_div(expr_copy(u), mk_int(2)));
            Expr* Tsq = mk_pow(expr_copy(T), mk_int(2));
            result = mk_div(mk_plus2(mk_int(1), mk_times2(mk_int(-1), Tsq)),
                             mk_times2(mk_int(2), T));
        }
    }

    if (result) {
        for (size_t i = 0; i < n; i++) expr_free(new_args[i]);
        free(new_args);
        return result;
    }
    Expr* newf = expr_new_function(expr_copy(head), new_args, n);
    free(new_args);
    return newf;
}

static Expr* convert_sincos_to_tan(Expr* f, Expr* x) {
    if (!f || !x || x->type != EXPR_SYMBOL) return NULL;
    Expr* r = convert_sincos_to_tan_rec(f, x->data.symbol, 0);
    if (!r) return NULL;
    Expr* normalised = eval_and_free(r);
    if (!normalised) return NULL;
    Expr* decot = decot_rec(normalised);
    expr_free(normalised);
    return decot;
}

static Expr* builtin_pm_sincos_to_tan(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)
        return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;
    return convert_sincos_to_tan(f, x);
}

/* Integrate`Helpers`PMPythagoreanRewrite[f] — rewrites Sec[u]^2,
 * Csc[u]^2, Sech[u]^2, Csch[u]^2 via the Pythagorean identities so the
 * result stays in the {Tan, Cot, Tanh, Coth} field. */
static Expr* builtin_pm_pythagorean(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION || res->data.function.arg_count != 1)
        return NULL;
    Expr* f = res->data.function.args[0];
    return pythagorean_rewrite(f);
}

/* Public-ish convert_to_tan: top-level entry.  The pipeline is:
 *   (1) recursive half-angle rewrite (convert_to_tan_rec)
 *   (2) full evaluate to normalise the algebra (this may collapse
 *       1/Tan to Cot, which we undo in step 3)
 *   (3) decot pass to put any residual Cot / Coth / Sec / ... back
 *       into Power[Tan, -1] / Power[Tanh, -1] form so the indet set
 *       is exactly { Tan[u/2], Tanh[u/2], ... }.
 *
 * Step 3 produces an UNEVALUATED tree — re-evaluating would undo the
 * decot.  Downstream pmint stages must work with the un-normalised
 * Power form. */
static Expr* convert_to_tan(Expr* f, Expr* x) {
    if (!f || !x || x->type != EXPR_SYMBOL) return NULL;
    Expr* rewritten = convert_to_tan_rec(f, x->data.symbol, 0);
    if (!rewritten) return NULL;
    Expr* normalised = eval_and_free(rewritten);
    if (!normalised) return NULL;
    Expr* decot = decot_rec(normalised);
    expr_free(normalised);
    return decot;
}

/* ------------------------------------------------------------------ */
/* Phase 2 — collect_indets_closed: gather transcendental atoms in    */
/* `ff` whose derivative wrt x is non-zero, closed under one round of */
/* differentiation.  pmint.maple lines 6-8.                            */
/* ------------------------------------------------------------------ */

/* A growing array of borrowed Expr* pointers (atoms still own their
 * memory in the caller's tree; we copy on output). */
typedef struct {
    Expr** items;
    size_t n;
    size_t cap;
} ExprList;

static void elist_init(ExprList* el) { el->items = NULL; el->n = 0; el->cap = 0; }

__attribute__((unused))
static void elist_free_items(ExprList* el) {
    for (size_t i = 0; i < el->n; i++) expr_free(el->items[i]);
    free(el->items);
    el->items = NULL; el->n = 0; el->cap = 0;
}

/* Append (taking ownership of e) iff no existing item is structurally
 * equal.  Returns true on append, false if already present (then frees
 * e); returns false also on capacity overflow (frees e). */
static bool elist_add_unique(ExprList* el, Expr* e) {
    if (!e) return false;
    for (size_t i = 0; i < el->n; i++) {
        if (expr_eq(el->items[i], e)) {
            expr_free(e);
            return false;
        }
    }
    if (el->n >= PMINT_MAX_INDETS) {
        expr_free(e);
        return false;
    }
    if (el->n >= el->cap) {
        size_t ncap = el->cap ? el->cap * 2 : 8;
        Expr** ni = (Expr**)realloc(el->items, sizeof(Expr*) * ncap);
        if (!ni) { expr_free(e); return false; }
        el->items = ni;
        el->cap = ncap;
    }
    el->items[el->n++] = e;
    return true;
}

/* Heads we collect as transcendental atoms.  Power[a, b] is collected
 * only when b is non-integer or when a or b mention x — pure x^n with
 * integer n is polynomial-in-x territory.
 *
 * Sin, Cos, Sec, Csc and their hyperbolic siblings are nominally
 * eliminated by convert_to_tan (rewritten to half-angle Tan / Tanh
 * rationals), but picocas's evaluator may rationalise them back to
 * Cot / Sec / Csc / ... forms when 1/Tan collapses.  We still
 * collect them so PMCollectIndets is robust to either form.  pmint's
 * downstream linear-system stage will fail gracefully (return NULL)
 * if the resulting set is algebraically over-determined. */
static bool head_is_transcendental_atom(const Expr* f, const char* xname) {
    if (!f || f->type != EXPR_FUNCTION) return false;
    Expr* h = f->data.function.head;
    if (!h || h->type != EXPR_SYMBOL) return false;
    const char* hs = h->data.symbol;
    /* Trig / hyperbolic field generators.  Tan / Tanh are pmint's
     * canonical generators (per convert_to_tan); the others get
     * collected only if they appear post-eval as evaluator artifacts. */
    if (hs == SYM_Tan || hs == SYM_Tanh) return true;
    if (hs == SYM_Cot || hs == SYM_Coth) return true;
    if (hs == SYM_Sin || hs == SYM_Cos || hs == SYM_Sec || hs == SYM_Csc) return true;
    if (hs == SYM_Sinh || hs == SYM_Cosh || hs == SYM_Sech || hs == SYM_Csch) return true;
    if (hs == SYM_Log) return true;
    if (hs == SYM_ArcSin || hs == SYM_ArcCos || hs == SYM_ArcTan
        || hs == SYM_ArcSinh || hs == SYM_ArcCosh || hs == SYM_ArcTanh) return true;
    /* LambertW if present in the build. */
    if (strcmp(hs, "LambertW") == 0) return true;
    /* Power[E, u] (i.e. Exp[u]) where u mentions x. */
    if (hs == SYM_Power
        && f->data.function.arg_count == 2
        && f->data.function.args[0]
        && f->data.function.args[0]->type == EXPR_SYMBOL
        && f->data.function.args[0]->data.symbol == SYM_E
        && !expr_free_of_symbol(f->data.function.args[1], xname)) {
        return true;
    }
    /* Power[base, exponent] where exponent is non-integer or symbolic
     * and mentions x — captures x^(1/2), x^a, etc. */
    if (hs == SYM_Power
        && f->data.function.arg_count == 2) {
        const Expr* base = f->data.function.args[0];
        const Expr* exp  = f->data.function.args[1];
        bool exp_is_int = exp && (exp->type == EXPR_INTEGER || exp->type == EXPR_BIGINT);
        if (!exp_is_int
            && (!expr_free_of_symbol(base, xname) || !expr_free_of_symbol(exp, xname))) {
            return true;
        }
    }
    return false;
}

/* Canonicalise the head of `expr` so algebraically-dependent
 * generators collapse to a single canonical atom.  In particular:
 *   Cot[u]  -> Tan[u]
 *   Coth[u] -> Tanh[u]
 * pmint's substitution map then carries an extra reciprocal rule so
 * the original Cot / Coth are still substituted correctly. */
static Expr* canonicalize_atom_head(Expr* expr) {
    if (!expr || expr->type != EXPR_FUNCTION) return expr_copy(expr);
    Expr* h = expr->data.function.head;
    if (!h || h->type != EXPR_SYMBOL) return expr_copy(expr);
    if (expr->data.function.arg_count != 1) return expr_copy(expr);
    const char* hs = h->data.symbol;
    if (hs == SYM_Cot) {
        return mk_unary(SYM_Tan, expr_copy(expr->data.function.args[0]));
    }
    if (hs == SYM_Coth) {
        return mk_unary(SYM_Tanh, expr_copy(expr->data.function.args[0]));
    }
    return expr_copy(expr);
}

/* Walk `expr` adding every transcendental-atom subexpression whose
 * derivative wrt x is non-zero. */
static void collect_atoms_rec(Expr* expr, const char* xname,
                              Expr* x_sym, ExprList* out) {
    if (!expr) return;

    if (head_is_transcendental_atom(expr, xname)
        && !expr_free_of_symbol(expr, xname)) {
        /* D[atom, x] != 0 guard. */
        Expr* d = call_d(expr, x_sym);
        bool nonzero = !is_int_zero(d);
        expr_free(d);
        if (nonzero) {
            /* Canonicalise the head so Cot[u] / Coth[u] collapse to
             * Tan[u] / Tanh[u] — they're algebraically dependent. */
            elist_add_unique(out, canonicalize_atom_head(expr));
            /* Continue descent — atoms can contain sub-atoms. */
        }
    }

    if (expr->type != EXPR_FUNCTION) return;
    collect_atoms_rec(expr->data.function.head, xname, x_sym, out);
    for (size_t i = 0; i < expr->data.function.arg_count; i++) {
        collect_atoms_rec(expr->data.function.args[i], xname, x_sym, out);
    }
}

/* Collect transcendental atoms in `ff` with D[a,x] ≠ 0, then close
 * under one round of differentiation (any new atoms appearing in the
 * derivatives are added too).  The lone `x` symbol is always the first
 * element of the output set.
 *
 * Returns 0 on success; writes owned Expr* array to *out_si, owned by
 * caller (free each item then the array).  Returns 1 on any failure
 * (cap exceeded, allocation error). */
static int collect_indets_closed(Expr* ff, Expr* x,
                                  Expr*** out_si, size_t* out_n) {
    *out_si = NULL;
    *out_n = 0;
    if (!ff || !x || x->type != EXPR_SYMBOL) return 1;

    const char* xname = x->data.symbol;
    ExprList el;
    elist_init(&el);

    /* Always include x first.  pmint relies on having x in the
     * substitution map so the resulting field is K(x, atoms). */
    if (!expr_free_of_symbol(ff, xname)) {
        elist_add_unique(&el, expr_copy(x));
    }

    collect_atoms_rec(ff, xname, x, &el);

    /* Closure: differentiate each atom in `el`, apply the Pythagorean
     * rewrite (so Sec[u]^2 / Csc[u]^2 / Sech[u]^2 / Csch[u]^2 turn
     * back into Tan / Cot / Tanh / Coth — keeping the generator set
     * to {x, Tan/Tanh/Log/Exp/...}), then collect atoms in the rewrite.
     * One pass is enough because pmint only needs the generators of
     * the differential field — any subsequent derivative is already a
     * polynomial in what we've collected. */
    size_t initial = el.n;
    for (size_t i = 1; i < initial; i++) {   /* skip x (index 0) */
        Expr* d = call_d(el.items[i], x);
        if (d) {
            Expr* rewritten = pythagorean_rewrite(d);
            expr_free(d);
            if (rewritten) {
                collect_atoms_rec(rewritten, xname, x, &el);
                expr_free(rewritten);
            }
        }
    }

    /* Move ownership out. */
    *out_si = el.items;
    *out_n = el.n;
    return 0;
}

/* ------------------------------------------------------------------ */
/* Phase 2 — build_substitution_maps: allocate pmint$v_i fresh        */
/* symbols and build the forward / reverse rule lists.                 */
/* ------------------------------------------------------------------ */

typedef struct {
    Expr* lhs;  /* original term     (owned) */
    Expr* rhs;  /* fresh-var symbol  (owned) */
} PMSubEntry;

typedef struct {
    PMSubEntry* items;
    size_t n;
} PMSubMap;

static void pm_sub_map_init(PMSubMap* m) { m->items = NULL; m->n = 0; }

static void pm_sub_map_free(PMSubMap* m) {
    if (!m) return;
    for (size_t i = 0; i < m->n; i++) {
        expr_free(m->items[i].lhs);
        expr_free(m->items[i].rhs);
    }
    free(m->items);
    m->items = NULL;
    m->n = 0;
}

/* Materialise a PMSubMap as List[Rule[lhs, rhs], ...] for use with
 * ReplaceAll.  When a rule maps Tan[u] -> v_k, also emits a
 * companion Cot[u] -> 1/v_k rule so the algebraically-dependent
 * generator is normalised through the same fresh variable.  Returns
 * owned Expr*. */
static Expr* pm_sub_map_to_rule_list(const PMSubMap* m) {
    /* Allocate worst-case 2× size: one extra rule per Tan / Tanh atom. */
    size_t cap = m->n * 2 + 1;
    Expr** items = (Expr**)malloc(sizeof(Expr*) * cap);
    if (!items) return NULL;
    size_t out_n = 0;
    for (size_t i = 0; i < m->n; i++) {
        items[out_n++] = mk_binary(SYM_Rule,
                                     expr_copy(m->items[i].lhs),
                                     expr_copy(m->items[i].rhs));
        /* Add reciprocal rule if lhs is Tan[u] or Tanh[u]. */
        Expr* lhs = m->items[i].lhs;
        if (lhs && lhs->type == EXPR_FUNCTION
            && lhs->data.function.arg_count == 1
            && lhs->data.function.head
            && lhs->data.function.head->type == EXPR_SYMBOL) {
            const char* hs = lhs->data.function.head->data.symbol;
            const char* recip_head = NULL;
            if (hs == SYM_Tan)  recip_head = SYM_Cot;
            else if (hs == SYM_Tanh) recip_head = SYM_Coth;
            if (recip_head) {
                Expr* recip_lhs = mk_unary(recip_head,
                                            expr_copy(lhs->data.function.args[0]));
                Expr* recip_rhs = mk_pow(expr_copy(m->items[i].rhs),
                                          mk_int(-1));
                items[out_n++] = mk_binary(SYM_Rule, recip_lhs, recip_rhs);
            }
        }
    }
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), items, out_n);
    free(items);
    return list;
}

/* Build lin (terms -> fresh vars), lout (fresh vars -> terms), and
 * vars (fresh symbols as Expr* per indet).  All outputs owned.
 *
 * The fresh variables are named pmint$v1, pmint$v2, ...  The leading
 * `pmint$` prefix gates them out of the user namespace and lets the
 * Phase 4 free-variable sweep zero out any unbound `_A`/`_B` /
 * `pmint$` symbols that survived the solver.
 *
 * The `gen_offset` parameter lets a caller request a name-suffix base
 * different from 1 if other pmint state has already allocated names.
 *
 * Returns 0 on success, 1 on failure. */
static int build_substitution_maps(Expr** si, size_t n,
                                    PMSubMap* out_lin, PMSubMap* out_lout,
                                    Expr*** out_vars, size_t gen_offset) {
    if (!si || !out_lin || !out_lout || !out_vars) return 1;
    pm_sub_map_init(out_lin);
    pm_sub_map_init(out_lout);
    *out_vars = NULL;

    if (n == 0) return 0;
    out_lin->items  = (PMSubEntry*)calloc(n, sizeof(PMSubEntry));
    out_lout->items = (PMSubEntry*)calloc(n, sizeof(PMSubEntry));
    *out_vars = (Expr**)calloc(n, sizeof(Expr*));
    if (!out_lin->items || !out_lout->items || !*out_vars) {
        free(out_lin->items);  out_lin->items  = NULL;
        free(out_lout->items); out_lout->items = NULL;
        free(*out_vars);       *out_vars       = NULL;
        return 1;
    }
    out_lin->n  = n;
    out_lout->n = n;

    for (size_t i = 0; i < n; i++) {
        char buf[32];
        snprintf(buf, sizeof(buf), "pmint$v%zu", gen_offset + i);
        const char* interned = intern_symbol(buf);

        Expr* fresh = expr_new_symbol(interned);

        out_lin->items[i].lhs  = expr_copy(si[i]);
        out_lin->items[i].rhs  = expr_copy(fresh);
        out_lout->items[i].lhs = expr_copy(fresh);
        out_lout->items[i].rhs = expr_copy(si[i]);
        (*out_vars)[i] = fresh;
    }
    return 0;
}

/* ------------------------------------------------------------------ */
/* Testable surfaces — Integrate`Helpers`PM*                            */
/* ------------------------------------------------------------------ */

/* Integrate`Helpers`PMConvertToTan[f, x] — exposes convert_to_tan
 * directly for REPL / unit testing. */
static Expr* builtin_pm_convert_to_tan(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)
        return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;
    return convert_to_tan(f, x);
}

/* Integrate`Helpers`PMCollectIndets[f, x] — returns the list of
 * collected transcendental atoms after the closure pass.  The first
 * element is always x. */
static Expr* builtin_pm_collect_indets(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)
        return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;

    Expr** si = NULL;
    size_t n = 0;
    if (collect_indets_closed(f, x, &si, &n) != 0) {
        for (size_t i = 0; i < n; i++) expr_free(si[i]);
        free(si);
        return NULL;
    }
    Expr* list = expr_new_function(expr_new_symbol(SYM_List), si, n);
    free(si);
    return list;
}

/* ------------------------------------------------------------------ */
/* Phase 3 — total-derivation vector field, splitFactor, deflation,    */
/* monomial enumeration.                                                */
/* ------------------------------------------------------------------ */

/* Apply ReplaceAll[expr, rules] via evaluate.  rules is owned by us
 * and consumed; expr is borrowed.  Returns owned result. */
static Expr* eval_replace_all(Expr* expr, Expr* rules) {
    return eval_and_free(mk_binary("ReplaceAll", expr_copy(expr), rules));
}

/* Evaluate Together[f]. */
static Expr* eval_together(Expr* f) {
    return eval_and_free(mk_unary("Together", expr_copy(f)));
}

/* Evaluate Expand[f]. */
static Expr* eval_expand(Expr* f) {
    return eval_and_free(mk_unary("Expand", expr_copy(f)));
}

/* Evaluate Numerator[f]. */
static Expr* eval_numer(Expr* f) {
    return eval_and_free(mk_unary("Numerator", expr_copy(f)));
}

/* Evaluate Denominator[f]. */
static Expr* eval_denom(Expr* f) {
    return eval_and_free(mk_unary("Denominator", expr_copy(f)));
}

/* Evaluate Cancel[f]. */
static Expr* eval_cancel(Expr* f) {
    return eval_and_free(mk_unary("Cancel", expr_copy(f)));
}

/* Evaluate PolynomialGCD[a, b]. */
static Expr* eval_poly_gcd(Expr* a, Expr* b) {
    return eval_and_free(mk_binary("PolynomialGCD", expr_copy(a), expr_copy(b)));
}

/* Evaluate Coefficient[p, x, n]. */
static Expr* eval_coefficient(Expr* p, Expr* x, int64_t deg) {
    Expr* args[3] = { expr_copy(p), expr_copy(x), mk_int(deg) };
    Expr* call = expr_new_function(expr_new_symbol("Coefficient"), args, 3);
    return eval_and_free(call);
}

/* Evaluate the degree of p in x.  picocas doesn't ship `Exponent`; we
 * compute it as Length[CoefficientList[p, x]] - 1.  Returns 0 when p
 * is free of x or constant. */
static int64_t eval_degree(Expr* p, Expr* x) {
    Expr* call = mk_binary("CoefficientList", expr_copy(p), expr_copy(x));
    Expr* clist = eval_and_free(call);
    int64_t d = 0;
    if (clist && clist->type == EXPR_FUNCTION
        && clist->data.function.head
        && clist->data.function.head->type == EXPR_SYMBOL
        && clist->data.function.head->data.symbol == SYM_List
        && clist->data.function.arg_count > 0) {
        d = (int64_t)clist->data.function.arg_count - 1;
    }
    expr_free(clist);
    return d;
}

/* Evaluate Variables[p]. Returns the List Expr* (owned). */
__attribute__((unused))
static Expr* eval_variables(Expr* p) {
    Expr* call = mk_unary("Variables", expr_copy(p));
    return eval_and_free(call);
}

/* Test whether `e` evaluates to zero (numeric, symbolic, or via
 * Cancel/Together).  Cheap path first. */
__attribute__((unused))
static bool is_zero_after_cancel(Expr* e) {
    if (!e) return true;
    if (e->type == EXPR_INTEGER && e->data.integer == 0) return true;
    Expr* c = eval_cancel(e);
    bool z = c && c->type == EXPR_INTEGER && c->data.integer == 0;
    expr_free(c);
    return z;
}

/* Build the total-derivation vector field for an integrand whose
 * indets are `si` (length n, including x at index 0).  Outputs:
 *   *out_l[k]  = q * Together[D[si[k], x] /. lin_rules]  (the scaled
 *                derivative of fresh-var vars[k], in fresh vars)
 *   *out_q     = lcm of denominators of the unscaled derivatives.
 *
 * Pmint.maple lines 11-15.  All outputs are caller-owned.
 *
 * `vars` is the list of fresh symbols (one per atom).  `lin` is the
 * forward substitution map. */
static int build_vector_field(Expr** si, size_t n,
                               PMSubMap* lin,
                               Expr* x,
                               Expr*** out_l,
                               Expr** out_q) {
    *out_l = NULL; *out_q = NULL;
    if (n == 0) return 1;

    Expr* lin_rules = pm_sub_map_to_rule_list(lin);
    Expr** numers = (Expr**)calloc(n, sizeof(Expr*));
    Expr** denoms = (Expr**)calloc(n, sizeof(Expr*));
    if (!numers || !denoms || !lin_rules) {
        free(numers); free(denoms); if (lin_rules) expr_free(lin_rules);
        return 1;
    }

    for (size_t k = 0; k < n; k++) {
        /* D[si[k], x] in original (un-substituted) form, then
         * Pythagorean rewrite (Sec^2 -> 1+Tan^2 etc.), then substitute
         * atoms -> fresh vars. */
        Expr* d = call_d(si[k], x);
        Expr* rew = pythagorean_rewrite(d);
        expr_free(d);
        Expr* sub = eval_replace_all(rew, expr_copy(lin_rules));
        expr_free(rew);
        Expr* tog = eval_together(sub);
        expr_free(sub);
        numers[k] = eval_numer(tog);
        denoms[k] = eval_denom(tog);
        expr_free(tog);
    }
    expr_free(lin_rules);

    /* q = lcm of all denoms.  lcm(a, b) = a*b / gcd(a, b), via Cancel. */
    Expr* q = expr_copy(denoms[0]);
    for (size_t k = 1; k < n; k++) {
        Expr* g = eval_poly_gcd(q, denoms[k]);
        Expr* prod = eval_expand(mk_times2(expr_copy(q), expr_copy(denoms[k])));
        Expr* div = mk_div(prod, g);
        expr_free(q);
        q = eval_cancel(div);
        expr_free(div);
    }

    /* l[k] = (q * numers[k]) / denoms[k], Expanded after Cancel. */
    Expr** lvec = (Expr**)calloc(n, sizeof(Expr*));
    if (!lvec) {
        for (size_t k = 0; k < n; k++) { expr_free(numers[k]); expr_free(denoms[k]); }
        free(numers); free(denoms);
        expr_free(q);
        return 1;
    }
    for (size_t k = 0; k < n; k++) {
        Expr* prod = mk_times2(expr_copy(q), expr_copy(numers[k]));
        Expr* ratio = mk_div(prod, expr_copy(denoms[k]));
        Expr* canc = eval_cancel(ratio);
        expr_free(ratio);
        lvec[k] = eval_expand(canc);
        expr_free(canc);
    }
    for (size_t k = 0; k < n; k++) { expr_free(numers[k]); expr_free(denoms[k]); }
    free(numers); free(denoms);

    *out_l = lvec;
    *out_q = q;
    return 0;
}

/* apply_d(f) — total derivation: Σ l[k] · D[f, vars[k]].
 *
 * Returns an unevaluated Plus[...] tree (no Together / Expand) so the
 * outer pipeline can decide when to normalise.  The caller is
 * responsible for any subsequent Together / Numerator. */
static Expr* apply_d(Expr* f, Expr** vars, Expr** l, size_t n) {
    Expr* sum = mk_int(0);
    for (size_t k = 0; k < n; k++) {
        Expr* dfk = call_d(f, vars[k]);
        if (!dfk) { expr_free(sum); return NULL; }
        Expr* term = mk_times2(expr_copy(l[k]), dfk);
        Expr* acc = mk_plus2(sum, term);
        sum = acc;
    }
    /* Evaluate just once to canonicalise.  Avoids Together (which on a
     * 15-unknown candidate with rational sub-terms is the bottleneck). */
    return eval_and_free(sum);
}

/* `vars_in_poly(p, vars, n)` — returns the subset of vars[i] that
 * actually appear in p (i.e., not FreeQ).  Output is an array of
 * borrowed pointers into the input vars; length written to *out_n.
 * `out_idx[i]` is the index in `vars` of the i-th surviving var. */
static void vars_in_poly(Expr* p, Expr** vars, size_t n,
                          size_t* out_idx, size_t* out_n) {
    *out_n = 0;
    for (size_t i = 0; i < n; i++) {
        if (vars[i] && vars[i]->type == EXPR_SYMBOL
            && !expr_free_of_symbol(p, vars[i]->data.symbol)) {
            out_idx[(*out_n)++] = i;
        }
    }
}

/* split_factor(p) — pmint.maple lines 80-90.
 *
 * Returns [s, h] such that p = s · h, where s is the "normal" part
 * (gcd-of-derivation factors removed) and h is the "special" part.
 *
 * Recursion: pick a var x in p with d(x)≠0; split content c (recurse)
 * and primitive part q.  Compute s = PolynomialGCD(q, d(q)) /
 * PolynomialGCD(q, ∂q/∂x).  If deg(s, x) == 0 return [spl_c[0], q ·
 * spl_c[1]]; else recurse on q/s.
 *
 * On entry `vars` are the fresh vars currently considered; `l` are
 * their scaled derivatives. */
static int split_factor(Expr* p,
                         Expr** vars, Expr** l, size_t n,
                         Expr** out_s, Expr** out_h) {
    *out_s = NULL; *out_h = NULL;

    /* Find a var x in p. */
    size_t idx_buf[PMINT_MAX_INDETS];
    size_t idxn = 0;
    vars_in_poly(p, vars, n, idx_buf, &idxn);
    if (idxn == 0) {
        *out_s = mk_int(1);
        *out_h = expr_copy(p);
        return 0;
    }
    Expr* xv = vars[idx_buf[0]];

    /* Content (gcd of coefficient list in xv) and primitive part. */
    Expr* p_exp = eval_expand(p);
    int64_t deg_p = eval_degree(p_exp, xv);

    /* Build content via PolynomialGCD fold of CoefficientList. */
    Expr* content = NULL;
    for (int64_t k = 0; k <= deg_p; k++) {
        Expr* ck = eval_coefficient(p_exp, xv, k);
        if (is_int_zero(ck)) { expr_free(ck); continue; }
        if (!content) {
            content = ck;
        } else {
            Expr* g = eval_poly_gcd(content, ck);
            expr_free(content); expr_free(ck);
            content = g;
        }
    }
    if (!content) content = mk_int(1);

    /* Primitive part q = Cancel[p / content]. */
    Expr* q = eval_cancel(mk_div(expr_copy(p_exp), expr_copy(content)));
    expr_free(p_exp);

    /* Recurse on content. */
    Expr* spl_c_s = NULL;
    Expr* spl_c_h = NULL;
    if (split_factor(content, vars, l, n, &spl_c_s, &spl_c_h) != 0) {
        expr_free(content); expr_free(q);
        return 1;
    }
    expr_free(content);

    /* d(q) and ∂q/∂xv. */
    Expr* dq    = apply_d(q, vars, l, n);
    Expr* pdq   = call_d(q, xv);
    Expr* dqrew = pythagorean_rewrite(pdq);
    expr_free(pdq);

    /* s = PolynomialGCD(q, dq) / PolynomialGCD(q, dqrew), Cancel'd.
     * Special case: deg(s, xv) == 0 → take [spl_c_s, q * spl_c_h]. */
    Expr* g1 = eval_poly_gcd(q, dq);
    Expr* g2 = eval_poly_gcd(q, dqrew);
    expr_free(dq); expr_free(dqrew);
    Expr* s_raw = mk_div(expr_copy(g1), expr_copy(g2));
    expr_free(g1); expr_free(g2);
    Expr* s = eval_cancel(s_raw);
    expr_free(s_raw);

    if (eval_degree(s, xv) == 0) {
        /* s contributes nothing to the split at this var. */
        expr_free(s);
        Expr* h = eval_expand(mk_times2(expr_copy(q), expr_copy(spl_c_h)));
        expr_free(spl_c_h); expr_free(q);
        *out_s = spl_c_s;
        *out_h = h;
        return 0;
    }

    /* Recurse on q / s. */
    Expr* q_over_s = eval_cancel(mk_div(expr_copy(q), expr_copy(s)));
    expr_free(q);
    Expr* splh_s = NULL;
    Expr* splh_h = NULL;
    if (split_factor(q_over_s, vars, l, n, &splh_s, &splh_h) != 0) {
        expr_free(q_over_s); expr_free(s);
        expr_free(spl_c_s); expr_free(spl_c_h);
        return 1;
    }
    expr_free(q_over_s);

    /* [spl_c[0] * splh[0] * s, spl_c[1] * splh[1]]. */
    Expr* s_final = eval_expand(
        mk_times2(mk_times2(expr_copy(spl_c_s), expr_copy(splh_s)),
                  expr_copy(s)));
    Expr* h_final = eval_expand(mk_times2(expr_copy(spl_c_h), expr_copy(splh_h)));
    expr_free(spl_c_s); expr_free(spl_c_h);
    expr_free(splh_s);  expr_free(splh_h);
    expr_free(s);

    *out_s = s_final;
    *out_h = h_final;
    return 0;
}

/* deflation(p) — pmint.maple lines 92-98.
 *
 * Returns deflation(content) * PolynomialGCD(q, ∂q/∂x). */
static Expr* deflation(Expr* p,
                        Expr** vars, Expr** l, size_t n) {
    (void)l;
    size_t idx_buf[PMINT_MAX_INDETS];
    size_t idxn = 0;
    vars_in_poly(p, vars, n, idx_buf, &idxn);
    if (idxn == 0) return expr_copy(p);

    Expr* xv = vars[idx_buf[0]];
    Expr* p_exp = eval_expand(p);
    int64_t deg_p = eval_degree(p_exp, xv);

    Expr* content = NULL;
    for (int64_t k = 0; k <= deg_p; k++) {
        Expr* ck = eval_coefficient(p_exp, xv, k);
        if (is_int_zero(ck)) { expr_free(ck); continue; }
        if (!content) content = ck;
        else {
            Expr* g = eval_poly_gcd(content, ck);
            expr_free(content); expr_free(ck);
            content = g;
        }
    }
    if (!content) content = mk_int(1);

    Expr* q = eval_cancel(mk_div(expr_copy(p_exp), expr_copy(content)));
    expr_free(p_exp);

    Expr* defl_c = deflation(content, vars, l, n);
    expr_free(content);

    Expr* dq_xv = call_d(q, xv);
    Expr* dq_rew = pythagorean_rewrite(dq_xv);
    expr_free(dq_xv);
    Expr* g = eval_poly_gcd(q, dq_rew);
    expr_free(q); expr_free(dq_rew);

    Expr* result = eval_expand(mk_times2(defl_c, g));
    return result;
}

/* enumerate_monomials(vars, nv, total_degree) — recursive build of
 * the set of monomials in vars[] of total degree ≤ total_degree.
 *
 * pmint.maple lines 69-78.
 *
 * Cap at PMINT_MAX_MONOMIALS; on overflow set *out_n = 0 and return 1
 * (caller bubbles back with ::overlarge). */
static int enumerate_monomials(Expr** vars, size_t nv, int total_degree,
                                Expr*** out_monoms, size_t* out_n) {
    *out_monoms = NULL;
    *out_n = 0;
    if (total_degree < 0) return 0;

    if (nv == 0) {
        /* enumerate empty set => {1}. */
        Expr** items = (Expr**)malloc(sizeof(Expr*));
        if (!items) return 1;
        items[0] = mk_int(1);
        *out_monoms = items;
        *out_n = 1;
        return 0;
    }

    /* Recursion: enumerate({vars[:-1]}, d) ∪
     *            {vars[-1]^i * w : 1 ≤ i ≤ d, w ∈ enumerate({vars[:-1]}, d-i)} */
    Expr* x = vars[nv - 1];

    Expr** lower = NULL;
    size_t lower_n = 0;
    if (enumerate_monomials(vars, nv - 1, total_degree, &lower, &lower_n) != 0) {
        return 1;
    }

    size_t cap = 16;
    Expr** acc = (Expr**)malloc(sizeof(Expr*) * cap);
    if (!acc) {
        for (size_t i = 0; i < lower_n; i++) expr_free(lower[i]);
        free(lower);
        return 1;
    }
    size_t accn = 0;

    /* d == 0 case from lower. */
    for (size_t i = 0; i < lower_n; i++) {
        if (accn >= PMINT_MAX_MONOMIALS) {
            for (size_t j = 0; j < accn; j++) expr_free(acc[j]);
            for (size_t j = 0; j < lower_n; j++) expr_free(lower[j]);
            free(acc); free(lower);
            return 1;
        }
        if (accn == cap) {
            cap *= 2;
            Expr** na = (Expr**)realloc(acc, sizeof(Expr*) * cap);
            if (!na) {
                for (size_t j = 0; j < accn; j++) expr_free(acc[j]);
                for (size_t j = 0; j < lower_n; j++) expr_free(lower[j]);
                free(acc); free(lower);
                return 1;
            }
            acc = na;
        }
        acc[accn++] = lower[i];
        lower[i] = NULL;  /* moved */
    }
    free(lower);

    /* i = 1..total_degree. */
    for (int i = 1; i <= total_degree; i++) {
        Expr** lower_i = NULL;
        size_t lower_i_n = 0;
        if (enumerate_monomials(vars, nv - 1, total_degree - i,
                                 &lower_i, &lower_i_n) != 0) {
            for (size_t j = 0; j < accn; j++) expr_free(acc[j]);
            free(acc);
            return 1;
        }
        for (size_t k = 0; k < lower_i_n; k++) {
            if (accn >= PMINT_MAX_MONOMIALS) {
                for (size_t j = 0; j < accn; j++) expr_free(acc[j]);
                for (size_t j = 0; j < lower_i_n; j++) expr_free(lower_i[j]);
                free(acc); free(lower_i);
                return 1;
            }
            if (accn == cap) {
                cap *= 2;
                Expr** na = (Expr**)realloc(acc, sizeof(Expr*) * cap);
                if (!na) {
                    for (size_t j = 0; j < accn; j++) expr_free(acc[j]);
                    for (size_t j = 0; j < lower_i_n; j++) expr_free(lower_i[j]);
                    free(acc); free(lower_i);
                    return 1;
                }
                acc = na;
            }
            Expr* xi = mk_pow(expr_copy(x), mk_int(i));
            Expr* prod = mk_times2(xi, lower_i[k]);
            acc[accn++] = eval_expand(prod);
            expr_free(prod);
        }
        free(lower_i);
    }

    *out_monoms = acc;
    *out_n = accn;
    return 0;
}

/* ------------------------------------------------------------------ */
/* Phase 4 — candidate ansatz, linear-system extraction, RowReduce.    */
/* ------------------------------------------------------------------ */

/* Build the candidate antiderivative's polynomial numerator
 *   num = Σ A_i · monomials[i]
 * (without the /cden division — the caller carries cden separately
 * and assembles cand = num/cden at the back-substitution step).
 *
 * Allocates fresh symbol names pmint$A1, pmint$A2, ...  Outputs:
 *   *out_num       = the unknown-bearing polynomial in fresh vars   (owned)
 *   *out_A_names   = interned-name pointers (interned strings live forever)
 *   *out_unknowns  = owned Expr* copies, one per A_i, for use as
 *                    substitution LHS in the linear-system solver. */
static int build_candidate(Expr** monomials, size_t nm,
                            Expr** out_num,
                            const char*** out_A_names, size_t* out_nA,
                            Expr*** out_unknowns) {
    *out_num = NULL;
    *out_A_names = NULL;
    *out_nA = 0;
    *out_unknowns = NULL;
    if (nm == 0) return 1;

    const char** names = (const char**)calloc(nm, sizeof(char*));
    Expr** unknowns = (Expr**)calloc(nm, sizeof(Expr*));
    if (!names || !unknowns) { free(names); free(unknowns); return 1; }

    Expr* sum = mk_int(0);
    for (size_t i = 0; i < nm; i++) {
        char buf[32];
        snprintf(buf, sizeof(buf), "pmint$A%zu", i + 1);
        names[i] = intern_symbol(buf);
        unknowns[i] = expr_new_symbol(names[i]);
        Expr* term = mk_times2(expr_copy(unknowns[i]), expr_copy(monomials[i]));
        sum = mk_plus2(sum, term);
    }
    Expr* sum_eval = eval_expand(sum);
    expr_free(sum);

    *out_num = sum_eval;
    *out_A_names = names;
    *out_nA = nm;
    *out_unknowns = unknowns;
    return 0;
}

/* Recursively walk a nested-list CoefficientList output, calling
 * `cb(exp_vec, coeff)` for each non-zero leaf.  `depth` tracks the
 * current variable index. */
typedef void (*coeff_visitor)(const int64_t* exp_vec, size_t nv,
                               Expr* coeff, void* userdata);

static void walk_coefficient_table(Expr* node, size_t depth, size_t nv,
                                    int64_t* exp_vec,
                                    coeff_visitor cb, void* user) {
    if (!node) return;
    if (depth == nv) {
        /* Leaf — the coefficient of v_1^e_1 ... v_nv^e_nv. */
        cb(exp_vec, nv, node, user);
        return;
    }
    /* Must be a List of length d+1. */
    if (node->type != EXPR_FUNCTION
        || !node->data.function.head
        || node->data.function.head->type != EXPR_SYMBOL
        || node->data.function.head->data.symbol != SYM_List) {
        /* Scalar at non-leaf depth — pad zeros for the remaining
         * dimensions. */
        for (size_t i = depth; i < nv; i++) exp_vec[i] = 0;
        cb(exp_vec, nv, node, user);
        return;
    }
    for (size_t i = 0; i < node->data.function.arg_count; i++) {
        exp_vec[depth] = (int64_t)i;
        walk_coefficient_table(node->data.function.args[i],
                                depth + 1, nv, exp_vec, cb, user);
    }
}

/* Row builder context for the linear-system extractor. */
typedef struct {
    Expr*** rows;           /* rows of mpq-as-Expr augmenter */
    size_t* row_caps;
    size_t nrows;
    size_t cap_rows;
    size_t ncols;           /* nunknowns + 1 (last col = -RHS) */
    Expr** unknowns;
    size_t nunknowns;
    bool ok;
} LinearBuilder;

static void linear_builder_add_row(LinearBuilder* lb, Expr* coeff_expr) {
    if (!lb->ok) return;
    if (is_zero_poly(coeff_expr)) return;   /* zero row, skip */

    /* Extract the row.  coeff_expr is supposed to be LINEAR in
     * lb->unknowns.  We compute:
     *   row[j] = Coefficient[coeff_expr, unknowns[j], 1]  for j < nunknowns
     *   row[nunknowns] = - (coeff_expr /. {unknown_k -> 0})  (negated RHS) */
    Expr* row = (Expr*)0;
    Expr** entries = (Expr**)calloc(lb->ncols, sizeof(Expr*));
    if (!entries) { lb->ok = false; return; }

    for (size_t j = 0; j < lb->nunknowns; j++) {
        Expr* c = eval_coefficient(coeff_expr, lb->unknowns[j], 1);
        entries[j] = c;
    }

    /* Constant term = expr with all unknowns set to 0. */
    Expr** zeros = (Expr**)malloc(sizeof(Expr*) * lb->nunknowns);
    for (size_t j = 0; j < lb->nunknowns; j++) {
        zeros[j] = mk_binary(SYM_Rule, expr_copy(lb->unknowns[j]), mk_int(0));
    }
    Expr* zero_rules = expr_new_function(expr_new_symbol(SYM_List),
                                          zeros, lb->nunknowns);
    free(zeros);
    Expr* const_part = eval_replace_all(coeff_expr, zero_rules);
    Expr* neg_const = eval_and_free(
        mk_times2(mk_int(-1), const_part));
    entries[lb->nunknowns] = neg_const;
    (void)row;

    /* Push the row. */
    if (lb->nrows >= lb->cap_rows) {
        size_t ncap = lb->cap_rows ? lb->cap_rows * 2 : 16;
        Expr*** nrows = (Expr***)realloc(lb->rows, sizeof(Expr**) * ncap);
        if (!nrows) {
            for (size_t j = 0; j < lb->ncols; j++) expr_free(entries[j]);
            free(entries);
            lb->ok = false;
            return;
        }
        lb->rows = nrows;
        lb->cap_rows = ncap;
    }
    lb->rows[lb->nrows++] = entries;
}

static void coeff_cb(const int64_t* exp_vec, size_t nv,
                      Expr* coeff, void* user) {
    (void)exp_vec; (void)nv;
    linear_builder_add_row((LinearBuilder*)user, coeff);
}

/* solve_linear_undet — extract the per-monomial coefficient rows of
 * equation_numer (a polynomial in vars, linear in unknowns), assemble
 * an augmented matrix [M | -b], call RowReduce, decode the solution.
 *
 * Outputs:
 *   *out_solution_rules  = List[Rule[unknown_j, value_j], ...]   owned
 *   *out_status          = 0 unique, 1 free vars (zeroed), -1 inconsistent
 *
 * Returns 0 on success (status set), 1 on internal error. */
static int solve_linear_undet(Expr* equation_numer,
                               Expr** vars, size_t nvars,
                               Expr** unknowns, size_t nunknowns,
                               Expr** out_solution_rules,
                               int* out_status) {
    *out_solution_rules = NULL;
    *out_status = -1;
    if (nunknowns == 0) return 1;

    /* Get CoefficientList[equation_numer, {vars...}] as a nested list. */
    Expr* vars_list = NULL;
    {
        Expr** items = (Expr**)malloc(sizeof(Expr*) * nvars);
        for (size_t i = 0; i < nvars; i++) items[i] = expr_copy(vars[i]);
        vars_list = expr_new_function(expr_new_symbol(SYM_List), items, nvars);
        free(items);
    }
    Expr* eq_exp = eval_expand(equation_numer);
    Expr* clist = eval_and_free(mk_binary("CoefficientList",
                                           eq_exp, vars_list));

    /* Walk the table, building linear rows. */
    LinearBuilder lb = {0};
    lb.rows = NULL;
    lb.row_caps = NULL;
    lb.nrows = 0;
    lb.cap_rows = 0;
    lb.ncols = nunknowns + 1;
    lb.unknowns = unknowns;
    lb.nunknowns = nunknowns;
    lb.ok = true;
    int64_t* exp_vec = (int64_t*)calloc(nvars, sizeof(int64_t));
    walk_coefficient_table(clist, 0, nvars, exp_vec, coeff_cb, &lb);
    free(exp_vec);
    expr_free(clist);

    if (!lb.ok) {
        for (size_t i = 0; i < lb.nrows; i++) {
            for (size_t j = 0; j < lb.ncols; j++) expr_free(lb.rows[i][j]);
            free(lb.rows[i]);
        }
        free(lb.rows);
        return 1;
    }

    if (lb.nrows == 0) {
        /* No constraints — every unknown is free.  Set them all to 0. */
        Expr** items = (Expr**)malloc(sizeof(Expr*) * nunknowns);
        for (size_t j = 0; j < nunknowns; j++) {
            items[j] = mk_binary(SYM_Rule, expr_copy(unknowns[j]), mk_int(0));
        }
        *out_solution_rules = expr_new_function(expr_new_symbol(SYM_List),
                                                  items, nunknowns);
        free(items);
        *out_status = 1;
        free(lb.rows);
        return 0;
    }

    /* Build augmented matrix as List[List[...]]. */
    Expr** mat_rows = (Expr**)malloc(sizeof(Expr*) * lb.nrows);
    for (size_t i = 0; i < lb.nrows; i++) {
        Expr** entries = (Expr**)malloc(sizeof(Expr*) * lb.ncols);
        for (size_t j = 0; j < lb.ncols; j++) entries[j] = lb.rows[i][j];
        free(lb.rows[i]);
        mat_rows[i] = expr_new_function(expr_new_symbol(SYM_List),
                                          entries, lb.ncols);
        free(entries);
    }
    free(lb.rows);
    Expr* matrix = expr_new_function(expr_new_symbol(SYM_List),
                                       mat_rows, lb.nrows);
    free(mat_rows);

    /* Call RowReduce.  Result is also a nested List. */
    Expr* rr = eval_and_free(mk_unary("RowReduce", matrix));
    if (!rr || rr->type != EXPR_FUNCTION
        || rr->data.function.head->type != EXPR_SYMBOL
        || rr->data.function.head->data.symbol != SYM_List) {
        if (rr) expr_free(rr);
        return 1;
    }

    /* Decode RREF.  Each row: walk for first non-zero entry (the pivot
     * column).  If pivot column < nunknowns: that unknown == rhs
     * (assuming all other entries in row × free vars = 0; free vars
     * pinned to 0 per pmint convention).  If pivot column == nunknowns
     * (and entry != 0): inconsistent.  Build solution: every unknown
     * with a pivot gets its rhs; unknowns without pivots are free → 0. */
    Expr** sol_vals = (Expr**)calloc(nunknowns, sizeof(Expr*));
    bool* pivoted   = (bool*)calloc(nunknowns, sizeof(bool));
    bool inconsistent = false;
    for (size_t r = 0; r < rr->data.function.arg_count; r++) {
        Expr* row_expr = rr->data.function.args[r];
        if (row_expr->type != EXPR_FUNCTION
            || row_expr->data.function.arg_count != lb.ncols) {
            inconsistent = true;
            break;
        }
        int64_t pivot_col = -1;
        for (size_t j = 0; j < lb.ncols; j++) {
            if (!is_zero_poly(row_expr->data.function.args[j])) {
                pivot_col = (int64_t)j;
                break;
            }
        }
        if (pivot_col < 0) continue;       /* zero row */
        if ((size_t)pivot_col == nunknowns) {
            /* All-zero coefficient with non-zero RHS → infeasible. */
            inconsistent = true;
            break;
        }
        /* Pivot in unknown column.  RHS sits in column nunknowns; the
         * pivot has been normalised to 1 by RowReduce.  Free-var
         * convention: pmint zeros every unknown without a pivot, so
         * the entries to the right (with j > pivot_col, j < nunknowns)
         * become 0 in the final solution — net effect is just RHS. */
        Expr* rhs = expr_copy(row_expr->data.function.args[nunknowns]);
        sol_vals[pivot_col] = rhs;
        pivoted[pivot_col] = true;
    }

    if (inconsistent) {
        for (size_t i = 0; i < nunknowns; i++) {
            if (sol_vals[i]) expr_free(sol_vals[i]);
        }
        free(sol_vals); free(pivoted);
        expr_free(rr);
        *out_status = -1;
        return 0;
    }

    bool any_free = false;
    Expr** items = (Expr**)malloc(sizeof(Expr*) * nunknowns);
    for (size_t j = 0; j < nunknowns; j++) {
        if (!pivoted[j]) {
            sol_vals[j] = mk_int(0);
            any_free = true;
        }
        items[j] = mk_binary(SYM_Rule, expr_copy(unknowns[j]), sol_vals[j]);
    }
    *out_solution_rules = expr_new_function(expr_new_symbol(SYM_List),
                                              items, nunknowns);
    free(items); free(sol_vals); free(pivoted);
    expr_free(rr);
    *out_status = any_free ? 1 : 0;
    return 0;
}

/* ------------------------------------------------------------------ */
/* Phase 5 — log candidates, getSpecial Darboux polys, K=I retry.      */
/* ------------------------------------------------------------------ */

/* Factor `p` into irreducible factors.  When `over_Qi` is non-zero,
 * use picocas's `Factor[p, Extension -> I]` (Trager algebraic
 * factoring); otherwise just `Factor[p]` over Z.
 *
 * Writes a newly-allocated array of owned Expr* factors into
 * *out_factors and the count into *out_n.  Unit factors (numeric
 * constants) and factors free of any vars[k] are excluded.
 *
 * Returns 0 on success, 1 on failure. */
static int my_factors(Expr* p, int over_Qi,
                       Expr** vars, size_t nv,
                       Expr*** out_factors, size_t* out_n) {
    *out_factors = NULL;
    *out_n = 0;

    Expr* factored;
    if (over_Qi) {
        Expr* ext_rule = mk_binary(SYM_Rule,
                                     expr_new_symbol("Extension"),
                                     expr_new_symbol(SYM_I));
        Expr* args[2] = { expr_copy(p), ext_rule };
        factored = eval_and_free(
            expr_new_function(expr_new_symbol("Factor"), args, 2));
    } else {
        Expr* arg = expr_copy(p);
        factored = eval_and_free(
            expr_new_function(expr_new_symbol("Factor"), &arg, 1));
    }
    if (!factored) return 1;

    /* Walk Factor's result, collecting irreducible factors.  The
     * structure is one of:
     *   - a single irreducible polynomial
     *   - Times[c, f1, f2, ...] (with c numeric)
     *   - Power[f, k] (treat as f)
     *   - Times[c, Power[f1, k1], f2, Power[f3, k3], ...] */
    Expr** acc = NULL;
    size_t cap = 0, count = 0;

    Expr* nodes_to_process[64];
    size_t np = 0;
    nodes_to_process[np++] = factored;

    while (np > 0 && count < 64) {
        Expr* node = nodes_to_process[--np];
        if (!node) continue;
        /* Skip purely numeric units. */
        if (node->type == EXPR_INTEGER || node->type == EXPR_REAL
            || node->type == EXPR_BIGINT) continue;
        if (node->type == EXPR_FUNCTION
            && node->data.function.head
            && node->data.function.head->type == EXPR_SYMBOL) {
            const char* hs = node->data.function.head->data.symbol;
            if (hs == SYM_Times) {
                for (size_t i = 0; i < node->data.function.arg_count && np < 60; i++) {
                    nodes_to_process[np++] = node->data.function.args[i];
                }
                continue;
            }
            if (hs == SYM_Power && node->data.function.arg_count == 2) {
                /* Factor's Power[base, k] denotes a multiplicity — we
                 * only need the base once. */
                nodes_to_process[np++] = node->data.function.args[0];
                continue;
            }
            if (hs == SYM_Rational) continue;   /* numeric */
        }
        /* Check that this factor actually mentions one of vars[]. */
        bool mentions_var = false;
        for (size_t k = 0; k < nv; k++) {
            if (vars[k] && vars[k]->type == EXPR_SYMBOL
                && !expr_free_of_symbol(node, vars[k]->data.symbol)) {
                mentions_var = true;
                break;
            }
        }
        if (!mentions_var) continue;

        /* Dedupe by structural equality. */
        bool dup = false;
        for (size_t i = 0; i < count; i++) {
            if (expr_eq(acc[i], node)) { dup = true; break; }
        }
        if (dup) continue;
        if (count >= cap) {
            cap = cap ? cap * 2 : 8;
            Expr** na = (Expr**)realloc(acc, sizeof(Expr*) * cap);
            if (!na) {
                for (size_t i = 0; i < count; i++) expr_free(acc[i]);
                free(acc);
                expr_free(factored);
                return 1;
            }
            acc = na;
        }
        acc[count++] = expr_copy(node);
    }
    expr_free(factored);

    *out_factors = acc;
    *out_n = count;
    return 0;
}

static void get_special_free(Expr** polys, bool* flags, size_t n) {
    if (polys) {
        for (size_t i = 0; i < n; i++) expr_free(polys[i]);
        free(polys);
    }
    free(flags);
}

/* getSpecial — pmint.maple lines 22-28.  For each atom in si[], if it
 * has a known Darboux polynomial under the total derivation, emit it
 * (in fresh-var form, so the result is a polynomial in vars[]).
 *
 * Output structure: parallel arrays of (poly, is_integral_flag).
 * `is_integral_flag` true means the polynomial should also be
 * multiplied into s in cden = s * spl[0] * deflation(spl[1]).
 * Returns 0 on success. */
static int get_special_all(Expr** si, Expr** vars, size_t n,
                            Expr*** out_darboux, bool** out_flags,
                            size_t* out_count) {
    *out_darboux = NULL;
    *out_flags = NULL;
    *out_count = 0;

    Expr** acc_p = NULL;
    bool* acc_f = NULL;
    size_t cap = 0, count = 0;

    for (size_t k = 0; k < n; k++) {
        if (!si[k] || si[k]->type != EXPR_FUNCTION) continue;
        Expr* head = si[k]->data.function.head;
        if (!head || head->type != EXPR_SYMBOL) continue;
        const char* hs = head->data.symbol;

        Expr* darboux1 = NULL;
        Expr* darboux2 = NULL;
        bool integral_flag = false;

        if (hs == SYM_Tan) {
            /* Darboux poly for Tan[u]: 1 + Tan[u]^2 = 1 + v_k^2. */
            darboux1 = mk_plus2(mk_int(1),
                                 mk_pow(expr_copy(vars[k]), mk_int(2)));
            integral_flag = false;
        } else if (hs == SYM_Tanh) {
            /* Two Darboux polys: 1 + Tanh[u] and 1 - Tanh[u]. */
            darboux1 = mk_plus2(mk_int(1), expr_copy(vars[k]));
            darboux2 = mk_plus2(mk_int(1),
                                 mk_times2(mk_int(-1), expr_copy(vars[k])));
            integral_flag = false;
        } else if (strcmp(hs, "LambertW") == 0) {
            /* Darboux poly for LambertW[u]: LambertW[u] = v_k, integral. */
            darboux1 = expr_copy(vars[k]);
            integral_flag = true;
        }

        Expr* tmp[2] = { darboux1, darboux2 };
        for (int j = 0; j < 2; j++) {
            if (!tmp[j]) continue;
            if (count >= cap) {
                cap = cap ? cap * 2 : 4;
                Expr** np = (Expr**)realloc(acc_p, sizeof(Expr*) * cap);
                if (!np) {
                    for (size_t i = 0; i < count; i++) expr_free(acc_p[i]);
                    free(acc_p); free(acc_f);
                    expr_free(darboux1);
                    if (darboux2) expr_free(darboux2);
                    return 1;
                }
                acc_p = np;
                bool* nf = (bool*)realloc(acc_f, sizeof(bool) * cap);
                if (!nf) {
                    for (size_t i = 0; i < count; i++) expr_free(acc_p[i]);
                    free(acc_p); free(acc_f);
                    expr_free(darboux1);
                    if (darboux2) expr_free(darboux2);
                    return 1;
                }
                acc_f = nf;
            }
            acc_p[count] = tmp[j];
            acc_f[count] = integral_flag;
            count++;
        }
    }

    *out_darboux = acc_p;
    *out_flags = acc_f;
    *out_count = count;
    return 0;
}

/* try_integral — Phase 4 inner loop (no log candidates yet).
 *
 *   ff_numer / ff_denom = ff_fresh   (already substituted to fresh
 *                                     vars, written as a numer/denom
 *                                     pair so all arithmetic stays
 *                                     polynomial)
 *   cand_num            = Σ A_i · monomials[i]   (polynomial in
 *                                                 fresh vars, linear
 *                                                 in unknowns)
 *   cden                = candidate denominator (polynomial)
 *
 * The pmint equation f − d(cand_num/cden)/q = 0 becomes, after
 * clearing all denominators (q, cden, ff_denom) and using the
 * quotient rule apply_d(u/v) = (apply_d(u)*v − u*apply_d(v))/v^2:
 *
 *   E :=  ff_numer * cden^2 * q
 *       − ff_denom * (apply_d(cand_num) * cden − cand_num * apply_d(cden))
 *
 * E is a pure polynomial in fresh vars, linear in the unknowns A_i.
 * Setting every monomial coefficient (in fresh vars) to zero gives
 * the linear system the solver tackles.
 *
 * On success writes the integrated result (back-substituted to atom
 * form) into *out_result and returns 0.  On infeasible returns 0 with
 * *out_result = NULL.  On error returns 1. */
static int try_integral_full(Expr* ff_numer, Expr* ff_denom,
                              Expr* cand_num, Expr* cden,
                              Expr** A_unknowns, size_t nA,
                              Expr** candlogs, size_t nlog,
                              Expr** vars, Expr** l, size_t nv,
                              Expr* q,
                              PMSubMap* lout,
                              Expr** out_result) {
    *out_result = NULL;

    /* Allocate B_j unknowns. */
    Expr** B_unknowns = (Expr**)calloc(nlog ? nlog : 1, sizeof(Expr*));
    for (size_t j = 0; j < nlog; j++) {
        char buf[32];
        snprintf(buf, sizeof(buf), "pmint$B%zu", j + 1);
        B_unknowns[j] = expr_new_symbol(intern_symbol(buf));
    }

    /* Combined unknowns array. */
    size_t nunk = nA + nlog;
    Expr** unknowns = (Expr**)calloc(nunk ? nunk : 1, sizeof(Expr*));
    for (size_t i = 0; i < nA; i++) unknowns[i] = expr_copy(A_unknowns[i]);
    for (size_t j = 0; j < nlog; j++) unknowns[nA + j] = expr_copy(B_unknowns[j]);

    PMTRACE("try_integral_full: nA=%zu nlog=%zu nv=%zu\n", nA, nlog, nv);

    /* apply_d on polynomials only.  Log candidates contribute
     * apply_d(g_j)/g_j after the divide-by-q, but we keep things as
     * polynomials by multiplying through.
     *
     * Equation (cleared of all denominators):
     *
     *   E :=  ff_numer * Π_j(g_j) * cden^2 * q
     *       - ff_denom * Π_j(g_j) * (apply_d(num) cden − num apply_d(cden))
     *       - ff_denom * cden^2 * Σ_j B_j * apply_d(g_j) * Π_{k≠j}(g_k)
     *
     * setting E ≡ 0 gives the linear system in {A_i, B_j}. */
    Expr* d_num  = apply_d(cand_num, vars, l, nv);
    if (!d_num) { goto err; }
    Expr* d_cden = apply_d(cden, vars, l, nv);
    if (!d_cden) { expr_free(d_num); goto err; }

    /* d_logs[j] = apply_d(g_j). */
    Expr** d_logs = (Expr**)calloc(nlog ? nlog : 1, sizeof(Expr*));
    for (size_t j = 0; j < nlog; j++) {
        d_logs[j] = apply_d(candlogs[j], vars, l, nv);
        if (!d_logs[j]) {
            for (size_t k = 0; k < j; k++) expr_free(d_logs[k]);
            free(d_logs);
            expr_free(d_num); expr_free(d_cden);
            goto err;
        }
    }

    /* prod_g = Π g_j. */
    Expr* prod_g = mk_int(1);
    for (size_t j = 0; j < nlog; j++) {
        prod_g = eval_expand(mk_times2(prod_g, expr_copy(candlogs[j])));
    }

    Expr* cden_sq = eval_expand(mk_times2(expr_copy(cden), expr_copy(cden)));

    /* term1 = ff_numer * prod_g * cden^2 * q   (no unknowns) */
    Expr* term1 = eval_expand(
        mk_times2(mk_times2(mk_times2(expr_copy(ff_numer),
                                        expr_copy(prod_g)),
                              expr_copy(cden_sq)),
                  expr_copy(q)));

    /* term2 = ff_denom * prod_g * (d_num cden − cand_num d_cden) */
    Expr* dpart_a = eval_expand(mk_times2(expr_copy(d_num), expr_copy(cden)));
    Expr* dpart_b = eval_expand(mk_times2(expr_copy(cand_num), expr_copy(d_cden)));
    Expr* dpart   = eval_expand(
        mk_plus2(dpart_a, mk_times2(mk_int(-1), dpart_b)));
    Expr* term2 = eval_expand(
        mk_times2(mk_times2(expr_copy(ff_denom), expr_copy(prod_g)), dpart));

    /* term3 = ff_denom * cden^2 * Σ_j B_j apply_d(g_j) Π_{k≠j}(g_k) */
    Expr* term3_sum = mk_int(0);
    for (size_t j = 0; j < nlog; j++) {
        /* g_others = Π_{k≠j}(g_k) — compute via prod_g / g_j (clean
         * since prod_g is just the product). */
        Expr* g_others = eval_cancel(mk_div(expr_copy(prod_g),
                                             expr_copy(candlogs[j])));
        Expr* contrib = mk_times2(
            mk_times2(expr_copy(B_unknowns[j]), expr_copy(d_logs[j])),
            g_others);
        term3_sum = mk_plus2(term3_sum, contrib);
    }
    Expr* term3 = eval_expand(
        mk_times2(mk_times2(expr_copy(ff_denom), expr_copy(cden_sq)), term3_sum));

    /* E = term1 − term2 − term3. */
    Expr* equation_numer = eval_expand(
        mk_plus2(term1,
                  mk_plus2(mk_times2(mk_int(-1), term2),
                            mk_times2(mk_int(-1), term3))));
    PMTRACE("equation built nrows≈?\n");

    expr_free(d_num); expr_free(d_cden);
    expr_free(prod_g); expr_free(cden_sq);
    for (size_t j = 0; j < nlog; j++) expr_free(d_logs[j]);
    free(d_logs);

    Expr* sol = NULL;
    int status = -1;
    int err = solve_linear_undet(equation_numer, vars, nv, unknowns, nunk,
                                  &sol, &status);
    expr_free(equation_numer);
    PMTRACE("solve returned err=%d status=%d\n", err, status);
    if (err != 0) goto err;
    if (status < 0) {
        if (sol) expr_free(sol);
        goto ret_none;
    }

    /* Build cand = cand_num/cden + Σ B_j Log[g_j], substitute solution,
     * then substitute fresh vars back to atoms. */
    Expr* log_sum = mk_int(0);
    for (size_t j = 0; j < nlog; j++) {
        Expr* logterm = mk_times2(expr_copy(B_unknowns[j]),
                                    mk_unary(SYM_Log, expr_copy(candlogs[j])));
        log_sum = mk_plus2(log_sum, logterm);
    }
    Expr* cand_full = mk_plus2(
        mk_div(expr_copy(cand_num), expr_copy(cden)),
        log_sum);
    Expr* cand_solved = eval_replace_all(cand_full, sol);
    expr_free(cand_full);

    Expr* lout_rules = pm_sub_map_to_rule_list(lout);
    Expr* result = eval_replace_all(cand_solved, lout_rules);
    expr_free(cand_solved);

    *out_result = result;
    goto cleanup;

err:
    for (size_t i = 0; i < nunk; i++) expr_free(unknowns[i]);
    free(unknowns);
    for (size_t j = 0; j < nlog; j++) expr_free(B_unknowns[j]);
    free(B_unknowns);
    return 1;

ret_none:
cleanup:
    for (size_t i = 0; i < nunk; i++) expr_free(unknowns[i]);
    free(unknowns);
    for (size_t j = 0; j < nlog; j++) expr_free(B_unknowns[j]);
    free(B_unknowns);
    return 0;
}

/* ------------------------------------------------------------------ */
/* Public pmint pipeline entry — replaces the Phase 1 stub.            */
/* ------------------------------------------------------------------ */

/* Total degree of p in vars (sum-of-individual-degrees upper bound;
 * conservative).  Returns 0 when p is free of all vars. */
static int64_t total_degree(Expr* p, Expr** vars, size_t n) {
    int64_t total = 0;
    for (size_t i = 0; i < n; i++) {
        total += eval_degree(p, vars[i]);
    }
    return total;
}

static Expr* rischnorman_integrate(Expr* f, Expr* x) {
    PMTRACE("entry\n");
    /* Phase 4 pipeline (no log candidates yet): */
    /*   1. convert_to_tan(f, x)                                     */
    /*   2. collect_indets_closed                                    */
    /*   3. build_substitution_maps / build_vector_field             */
    /*   4. splitFactor(q, d) → splq; splitFactor(denom(ff), d) → spl */
    /*   5. cden = splq[0] * spl[0] * deflation(spl[1])              */
    /*   6. dg = 1 + total_deg(splq[0])                              */
    /*           + max(total_deg(numer ff), total_deg(denom ff))     */
    /*   7. enumerate_monomials(vars, dg)                            */
    /*   8. build_candidate                                          */
    /*   9. try_integral                                             */

    Expr* ff_atoms = convert_to_tan(f, x);
    if (!ff_atoms) return NULL;
    PMTRACE("after convert_to_tan, n_atoms=?\n");

    Expr** si = NULL; size_t n = 0;
    if (collect_indets_closed(ff_atoms, x, &si, &n) != 0 || n == 0) {
        expr_free(ff_atoms);
        for (size_t i = 0; i < n; i++) expr_free(si[i]);
        free(si);
        return NULL;
    }

    PMSubMap lin, lout;
    Expr** vars = NULL;
    if (build_substitution_maps(si, n, &lin, &lout, &vars, 1) != 0) {
        for (size_t i = 0; i < n; i++) expr_free(si[i]);
        free(si);
        expr_free(ff_atoms);
        return NULL;
    }

    Expr** lvec = NULL; Expr* q = NULL;
    if (build_vector_field(si, n, &lin, x, &lvec, &q) != 0) {
        for (size_t i = 0; i < n; i++) expr_free(si[i]);
        free(si);
        pm_sub_map_free(&lin); pm_sub_map_free(&lout);
        for (size_t i = 0; i < n; i++) expr_free(vars[i]);
        free(vars);
        expr_free(ff_atoms);
        return NULL;
    }

    /* Substitute ff_atoms to fresh-var form. */
    Expr* lin_rules = pm_sub_map_to_rule_list(&lin);
    Expr* ff_fresh = eval_expand(eval_replace_all(ff_atoms, lin_rules));
    expr_free(ff_atoms);

    PMTRACE("about to splitFactor(q)\n");
    /* splitFactor(q, dx). */
    Expr* splq_s = NULL, *splq_h = NULL;
    if (split_factor(q, vars, lvec, n, &splq_s, &splq_h) != 0) {
        expr_free(ff_fresh); expr_free(q);
        for (size_t i = 0; i < n; i++) { expr_free(si[i]); expr_free(vars[i]); expr_free(lvec[i]); }
        free(si); free(vars); free(lvec);
        pm_sub_map_free(&lin); pm_sub_map_free(&lout);
        return NULL;
    }

    PMTRACE("after splitFactor(q)\n");
    /* ff = Together[f_fresh], df = Denominator. */
    Expr* ff_together = eval_together(ff_fresh);
    Expr* ff_numer = eval_numer(ff_together);
    Expr* ff_denom = eval_denom(ff_together);
    expr_free(ff_together);

    /* splitFactor(df, dx). */
    Expr* spl_s = NULL, *spl_h = NULL;
    if (split_factor(ff_denom, vars, lvec, n, &spl_s, &spl_h) != 0) {
        expr_free(ff_fresh); expr_free(ff_numer); expr_free(ff_denom);
        expr_free(splq_s); expr_free(splq_h);
        expr_free(q);
        for (size_t i = 0; i < n; i++) { expr_free(si[i]); expr_free(vars[i]); expr_free(lvec[i]); }
        free(si); free(vars); free(lvec);
        pm_sub_map_free(&lin); pm_sub_map_free(&lout);
        return NULL;
    }

    /* Compute Darboux specials.  Their integral-flagged polys go into
     * the `s` part of cden; non-integral ones become log candidates. */
    Expr** darboux_polys = NULL;
    bool* darboux_flags = NULL;
    size_t n_darboux = 0;
    get_special_all(si, vars, n, &darboux_polys, &darboux_flags, &n_darboux);

    /* s_part = splq_s * (product of integral-flagged Darboux). */
    Expr* s_part = expr_copy(splq_s);
    for (size_t i = 0; i < n_darboux; i++) {
        if (!darboux_flags[i]) continue;
        s_part = eval_expand(mk_times2(s_part, expr_copy(darboux_polys[i])));
    }

    /* cden = s_part * spl_s * deflation(spl_h). */
    Expr* defl_spl_h = deflation(spl_h, vars, lvec, n);
    Expr* cden = eval_expand(
        mk_times2(mk_times2(expr_copy(s_part), expr_copy(spl_s)),
                  defl_spl_h));
    expr_free(s_part);
    /* Avoid zero cden. */
    if (is_zero_poly(cden)) {
        expr_free(cden);
        cden = mk_int(1);
    }

    /* dg = 1 + deg(splq_s) + max(deg(numer), deg(denom)). */
    int64_t deg_splq_s = total_degree(splq_s, vars, n);
    int64_t deg_num    = total_degree(ff_numer, vars, n);
    int64_t deg_den    = total_degree(ff_denom, vars, n);
    int64_t dg = 1 + deg_splq_s + (deg_num > deg_den ? deg_num : deg_den);
    if (dg < 0) dg = 0;
    expr_free(ff_numer); expr_free(ff_denom);

    /* Keep splq_s, spl_s, spl_h around for the candlog factor lists in
     * the K=0/K=I retry loop below — they get factored over Q (then
     * Q[i]). */
    Expr* splq_s_keep = splq_s;   splq_s = NULL;
    Expr* spl_s_keep  = spl_s;    spl_s = NULL;
    Expr* spl_h_keep  = spl_h;    spl_h = NULL;
    expr_free(splq_h);

    /* Enumerate monomials. */
    Expr** monoms = NULL; size_t nmon = 0;
    if (enumerate_monomials(vars, n, (int)dg, &monoms, &nmon) != 0
        || nmon == 0) {
        expr_free(ff_fresh); expr_free(cden); expr_free(q);
        for (size_t i = 0; i < n; i++) { expr_free(si[i]); expr_free(vars[i]); expr_free(lvec[i]); }
        free(si); free(vars); free(lvec);
        pm_sub_map_free(&lin); pm_sub_map_free(&lout);
        for (size_t i = 0; i < nmon; i++) expr_free(monoms[i]);
        free(monoms);
        return NULL;
    }

    PMTRACE("dg=%lld nmon=%zu\n", (long long)dg, nmon);
    /* Build candidate. */
    Expr* cand_num = NULL;
    const char** A_names = NULL; size_t nA = 0;
    Expr** unknowns_A = NULL;
    if (build_candidate(monoms, nmon, &cand_num, &A_names, &nA, &unknowns_A) != 0) {
        for (size_t i = 0; i < nmon; i++) expr_free(monoms[i]);
        free(monoms);
        expr_free(ff_fresh); expr_free(cden); expr_free(q);
        expr_free(splq_s_keep); expr_free(spl_s_keep); expr_free(spl_h_keep);
        get_special_free(darboux_polys, darboux_flags, n_darboux);
        for (size_t i = 0; i < n; i++) { expr_free(si[i]); expr_free(vars[i]); expr_free(lvec[i]); }
        free(si); free(vars); free(lvec);
        pm_sub_map_free(&lin); pm_sub_map_free(&lout);
        return NULL;
    }
    for (size_t i = 0; i < nmon; i++) expr_free(monoms[i]);
    free(monoms);

    /* Decompose ff_fresh = ff_n / ff_d so the equation can stay
     * polynomial (no Together inside try_integral). */
    Expr* ff_t = eval_together(ff_fresh);
    Expr* ff_n = eval_numer(ff_t);
    Expr* ff_d = eval_denom(ff_t);
    expr_free(ff_t);
    expr_free(ff_fresh);

    PMTRACE("about to try_integral_full\n");
    Expr* result = NULL;

    for (int K_iteration = 0; K_iteration < 2 && !result; K_iteration++) {
        int over_Qi = K_iteration;

        /* Build candlog list: factors of splq_s, spl_s, spl_h
         * (over Q or Q[i] depending on K_iteration), plus non-integral
         * Darboux polys. */
        Expr** candlogs = NULL;
        size_t ncandlogs = 0;

        Expr** facs_splq_s = NULL; size_t nf_splq_s = 0;
        Expr** facs_spl_s  = NULL; size_t nf_spl_s = 0;
        Expr** facs_spl_h  = NULL; size_t nf_spl_h = 0;
        my_factors(splq_s_keep, over_Qi, vars, n, &facs_splq_s, &nf_splq_s);
        my_factors(spl_s_keep, over_Qi, vars, n, &facs_spl_s, &nf_spl_s);
        my_factors(spl_h_keep, over_Qi, vars, n, &facs_spl_h, &nf_spl_h);

        size_t total_candidates = nf_splq_s + nf_spl_s + nf_spl_h + n_darboux;
        candlogs = (Expr**)calloc(total_candidates ? total_candidates : 1,
                                    sizeof(Expr*));
        for (size_t i = 0; i < nf_splq_s; i++) {
            bool dup = false;
            for (size_t k = 0; k < ncandlogs; k++)
                if (expr_eq(candlogs[k], facs_splq_s[i])) { dup = true; break; }
            if (!dup) candlogs[ncandlogs++] = expr_copy(facs_splq_s[i]);
        }
        for (size_t i = 0; i < nf_spl_s; i++) {
            bool dup = false;
            for (size_t k = 0; k < ncandlogs; k++)
                if (expr_eq(candlogs[k], facs_spl_s[i])) { dup = true; break; }
            if (!dup) candlogs[ncandlogs++] = expr_copy(facs_spl_s[i]);
        }
        for (size_t i = 0; i < nf_spl_h; i++) {
            bool dup = false;
            for (size_t k = 0; k < ncandlogs; k++)
                if (expr_eq(candlogs[k], facs_spl_h[i])) { dup = true; break; }
            if (!dup) candlogs[ncandlogs++] = expr_copy(facs_spl_h[i]);
        }
        /* Non-integral Darboux polys join the candlog list. */
        for (size_t i = 0; i < n_darboux; i++) {
            if (darboux_flags[i]) continue;
            bool dup = false;
            for (size_t k = 0; k < ncandlogs; k++)
                if (expr_eq(candlogs[k], darboux_polys[i])) { dup = true; break; }
            if (!dup) candlogs[ncandlogs++] = expr_copy(darboux_polys[i]);
        }

        for (size_t i = 0; i < nf_splq_s; i++) expr_free(facs_splq_s[i]);
        free(facs_splq_s);
        for (size_t i = 0; i < nf_spl_s; i++) expr_free(facs_spl_s[i]);
        free(facs_spl_s);
        for (size_t i = 0; i < nf_spl_h; i++) expr_free(facs_spl_h[i]);
        free(facs_spl_h);

        PMTRACE("K_iteration=%d ncandlogs=%zu\n", K_iteration, ncandlogs);

        int err = try_integral_full(ff_n, ff_d, cand_num, cden,
                                     unknowns_A, nA,
                                     candlogs, ncandlogs,
                                     vars, lvec, n, q, &lout, &result);

        for (size_t i = 0; i < ncandlogs; i++) expr_free(candlogs[i]);
        free(candlogs);
        if (err != 0) { result = NULL; break; }
    }

    expr_free(cand_num);
    expr_free(ff_n); expr_free(ff_d);
    expr_free(cden);
    expr_free(splq_s_keep); expr_free(spl_s_keep); expr_free(spl_h_keep);
    get_special_free(darboux_polys, darboux_flags, n_darboux);
    for (size_t i = 0; i < nA; i++) expr_free(unknowns_A[i]);
    free(unknowns_A); free(A_names);
    expr_free(q);
    for (size_t i = 0; i < n; i++) { expr_free(si[i]); expr_free(vars[i]); expr_free(lvec[i]); }
    free(si); free(vars); free(lvec);
    pm_sub_map_free(&lin); pm_sub_map_free(&lout);

    if (!result) return NULL;

    /* ----------------------------------------------------------- */
    /* Output cleanup: rewrite Tan[u/2] / Cot[u/2] / Tanh[u/2] /   */
    /* Coth[u/2] half-angle terms back into Sin[u] / Cos[u] /      */
    /* Sinh[u] / Cosh[u] using the closed-form Weierstrass         */
    /* inverses, then per-term Cancel to clean up.                  */
    /*                                                              */
    /* Identities used (no Simplify — Cancel-only path):           */
    /*   Tan[u/2]^2  -> (1 - Cos[u]) / (1 + Cos[u])                 */
    /*   Tan[u/2]    -> Sin[u] / (1 + Cos[u])                       */
    /*   Cot[u/2]^2  -> (1 + Cos[u]) / (1 - Cos[u])                 */
    /*   Cot[u/2]    -> (1 + Cos[u]) / Sin[u]                       */
    /*   Tanh[u/2]^2 -> (Cosh[u] - 1) / (Cosh[u] + 1)               */
    /*   Tanh[u/2]   -> Sinh[u] / (1 + Cosh[u])                     */
    /*   Coth[u/2]^2 -> (Cosh[u] + 1) / (Cosh[u] - 1)               */
    /*   Coth[u/2]   -> (Cosh[u] + 1) / Sinh[u]                     */
    /*                                                              */
    /* Then Together to combine fractions, and per-term Cancel if  */
    /* the top-level head is Plus.  Avoids the expensive Simplify  */
    /* pass that would be needed for fuller trig reduction. */
    {
        Expr* rules_list = parse_expression(
            "{Tan[u_/2]^2 -> (1 - Cos[u])/(1 + Cos[u]),"
            " Tan[u_/2] -> Sin[u]/(1 + Cos[u]),"
            " Cot[u_/2]^2 -> (1 + Cos[u])/(1 - Cos[u]),"
            " Cot[u_/2] -> (1 + Cos[u])/Sin[u],"
            " Tanh[u_/2]^2 -> (Cosh[u] - 1)/(Cosh[u] + 1),"
            " Tanh[u_/2] -> Sinh[u]/(1 + Cosh[u]),"
            " Coth[u_/2]^2 -> (Cosh[u] + 1)/(Cosh[u] - 1),"
            " Coth[u_/2] -> (Cosh[u] + 1)/Sinh[u]}");
        /* Secondary cleanup rules: Log[a/b] -> Log[a] - Log[b] so
         * fractional log arguments produced by the Tan-half-angle
         * substitutions can fold into each other.  Valid up to a
         * branch-cut constant; pmint's antiderivatives are defined
         * modulo an arbitrary constant. */
        Expr* log_rules = parse_expression(
            "{Log[a_/b_] :> Log[a] - Log[b]}");
        if (rules_list) {
            Expr* rewritten = eval_replace_all(result, rules_list);
            expr_free(result);
            /* Together inside Log arguments: walk the tree, apply
             * Together to each Log's argument so Log[1 + (1-c)/(1+c)]
             * collapses to Log[2/(1+c)] before the Log[a/b] expansion
             * fires. */
            Expr* log_together_rule = parse_expression(
                "{Log[a_] :> Log[Together[a]]}");
            if (log_together_rule) {
                Expr* tightened = eval_replace_all(rewritten, log_together_rule);
                expr_free(rewritten);
                rewritten = tightened;
            }
            if (log_rules) {
                Expr* with_logs = eval_replace_all(rewritten, log_rules);
                expr_free(rewritten);
                rewritten = with_logs;
            }
            Expr* combined = eval_together(rewritten);
            expr_free(rewritten);

            /* Fresh-symbol PolynomialGCD reduction.  picocas's Cancel
             * treats Cos[x] / Sin[x] / E^x as opaque transcendental
             * atoms, so common factors like (1 + Cos[x]) between
             * numerator and denominator don't simplify.  We substitute
             * each common transcendental to a fresh polynomial-style
             * symbol, factor out the GCD, then substitute back.
             *
             * This handles e.g. Sin[x] Exp[x]'s output, where Cancel
             * alone leaves `4 E^x (1+Cos[x])(Sin[x]-Cos[x]) / (8(1+Cos[x]))`
             * uncancelled. */
            const char* x_name = (x && x->type == EXPR_SYMBOL)
                ? x->data.symbol : NULL;
            if (x_name && combined) {
                char buf[1024];
                snprintf(buf, sizeof(buf),
                    "{Cos[%s] -> pmint$cc, Sin[%s] -> pmint$ss,"
                    " Cosh[%s] -> pmint$ch, Sinh[%s] -> pmint$sh,"
                    " Log[%s] -> pmint$lg, E^%s -> pmint$ee}",
                    x_name, x_name, x_name, x_name, x_name, x_name);
                Expr* sub_fwd = parse_expression(buf);

                snprintf(buf, sizeof(buf),
                    "{pmint$cc -> Cos[%s], pmint$ss -> Sin[%s],"
                    " pmint$ch -> Cosh[%s], pmint$sh -> Sinh[%s],"
                    " pmint$lg -> Log[%s], pmint$ee -> E^%s}",
                    x_name, x_name, x_name, x_name, x_name, x_name);
                Expr* sub_back = parse_expression(buf);

                if (sub_fwd && sub_back) {
                    Expr* with_fresh = eval_replace_all(combined, sub_fwd);
                    Expr* fresh_cancelled = eval_cancel(with_fresh);
                    expr_free(with_fresh);
                    Expr* restored = eval_replace_all(fresh_cancelled, sub_back);
                    expr_free(fresh_cancelled);
                    expr_free(combined);
                    combined = restored;
                } else {
                    if (sub_fwd) expr_free(sub_fwd);
                    if (sub_back) expr_free(sub_back);
                }
            }

            /* If the top-level head is Plus, Cancel each summand
             * separately — picocas's Cancel is more effective on
             * simpler subexpressions than on a single Together'd
             * fraction with multiple variables. */
            if (combined && combined->type == EXPR_FUNCTION
                && combined->data.function.head
                && combined->data.function.head->type == EXPR_SYMBOL
                && combined->data.function.head->data.symbol == SYM_Plus) {
                size_t k = combined->data.function.arg_count;
                Expr** new_args = (Expr**)malloc(sizeof(Expr*) * (k ? k : 1));
                for (size_t i = 0; i < k; i++) {
                    new_args[i] = eval_cancel(combined->data.function.args[i]);
                }
                Expr* head_copy = expr_copy(combined->data.function.head);
                expr_free(combined);
                result = expr_new_function(head_copy, new_args, k);
                free(new_args);
                /* Evaluate once more so picocas re-canonicalises the Plus. */
                result = eval_and_free(result);
            } else {
                result = eval_cancel(combined);
                expr_free(combined);
            }
        }
    }

    return result;
}

/* ------------------------------------------------------------------ */
/* Testable surfaces for Phase 3.                                      */
/*                                                                      */
/* All take user-level expressions (in original atoms) and return     */
/* outputs in user-level form.  Internally they go through the pmint  */
/* preprocessing (convert_to_tan, collect_indets, build_subst_maps,   */
/* build_vector_field) so the test inputs read naturally.             */
/* ------------------------------------------------------------------ */

/* Build the full pmint-pipeline preprocessing state for an integrand.
 * Output: si (atoms), lin / lout (subst maps), vars (fresh symbols),
 * l (scaled derivatives), q (lcm denom).  All caller-owned.
 *
 * Returns 0 on success. */
static int build_pipeline_state(Expr* f, Expr* x,
                                 Expr*** out_si, size_t* out_n,
                                 PMSubMap* out_lin, PMSubMap* out_lout,
                                 Expr*** out_vars,
                                 Expr*** out_l, Expr** out_q) {
    Expr* ff = convert_to_tan(f, x);
    if (!ff) return 1;
    if (collect_indets_closed(ff, x, out_si, out_n) != 0) {
        expr_free(ff);
        return 1;
    }
    expr_free(ff);
    if (build_substitution_maps(*out_si, *out_n,
                                 out_lin, out_lout, out_vars, 1) != 0) {
        for (size_t i = 0; i < *out_n; i++) expr_free((*out_si)[i]);
        free(*out_si); *out_si = NULL;
        return 1;
    }
    if (build_vector_field(*out_si, *out_n, out_lin, x, out_l, out_q) != 0) {
        for (size_t i = 0; i < *out_n; i++) expr_free((*out_si)[i]);
        free(*out_si); *out_si = NULL;
        pm_sub_map_free(out_lin);
        pm_sub_map_free(out_lout);
        for (size_t i = 0; i < *out_n; i++) expr_free((*out_vars)[i]);
        free(*out_vars); *out_vars = NULL;
        return 1;
    }
    return 0;
}

static void free_pipeline_state(Expr** si, size_t n,
                                 PMSubMap* lin, PMSubMap* lout,
                                 Expr** vars,
                                 Expr** l, Expr* q) {
    for (size_t i = 0; i < n; i++) expr_free(si[i]);
    free(si);
    pm_sub_map_free(lin);
    pm_sub_map_free(lout);
    if (vars) {
        for (size_t i = 0; i < n; i++) expr_free(vars[i]);
        free(vars);
    }
    if (l) {
        for (size_t i = 0; i < n; i++) expr_free(l[i]);
        free(l);
    }
    if (q) expr_free(q);
}

/* Integrate`Helpers`PMVectorField[f, x] — returns
 * {indets, lin_rules, l_vec, q} where indets is the list of fresh-var
 * substitutions, lin is the substitution rule list, l_vec is the
 * scaled derivative vector, and q is the common denominator. */
static Expr* builtin_pm_vector_field(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)
        return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;

    Expr** si = NULL; size_t n = 0;
    PMSubMap lin, lout;
    Expr** vars = NULL;
    Expr** l = NULL; Expr* q = NULL;

    if (build_pipeline_state(f, x, &si, &n, &lin, &lout,
                              &vars, &l, &q) != 0) {
        return NULL;
    }

    /* Build {vars, l_list, q}. */
    Expr** var_items = (Expr**)malloc(sizeof(Expr*) * n);
    Expr** l_items   = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        var_items[i] = expr_copy(vars[i]);
        l_items[i]   = expr_copy(l[i]);
    }
    Expr* vars_list = expr_new_function(expr_new_symbol(SYM_List), var_items, n);
    Expr* l_list    = expr_new_function(expr_new_symbol(SYM_List), l_items, n);
    free(var_items); free(l_items);
    Expr* q_copy = expr_copy(q);

    Expr* result_items[3] = { vars_list, l_list, q_copy };
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), result_items, 3);

    free_pipeline_state(si, n, &lin, &lout, vars, l, q);
    return result;
}

/* Integrate`Helpers`PMApplyD[f_in_atoms, x] — applies the total
 * derivation to an integrand-level expression.  The expression is
 * first substituted to fresh vars; apply_d is called; the result is
 * substituted back via lout.  Returns owned Expr*. */
static Expr* builtin_pm_apply_d(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)
        return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;

    Expr** si = NULL; size_t n = 0;
    PMSubMap lin, lout;
    Expr** vars = NULL;
    Expr** l = NULL; Expr* q = NULL;

    if (build_pipeline_state(f, x, &si, &n, &lin, &lout,
                              &vars, &l, &q) != 0) {
        return NULL;
    }

    /* Substitute f to fresh-var representation. */
    Expr* ff = convert_to_tan(f, x);
    if (!ff) {
        free_pipeline_state(si, n, &lin, &lout, vars, l, q);
        return NULL;
    }
    Expr* lin_rules = pm_sub_map_to_rule_list(&lin);
    Expr* f_subbed = eval_replace_all(ff, lin_rules);
    expr_free(ff);

    Expr* df = apply_d(f_subbed, vars, l, n);
    expr_free(f_subbed);

    /* Substitute back. */
    Expr* lout_rules = pm_sub_map_to_rule_list(&lout);
    Expr* result = eval_replace_all(df, lout_rules);
    expr_free(df);

    free_pipeline_state(si, n, &lin, &lout, vars, l, q);
    return result;
}

/* Integrate`Helpers`PMSplitFactor[p_in_atoms, x] — returns {s, h}
 * such that p = s · h, after substituting p to fresh-var form. */
static Expr* builtin_pm_split_factor(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)
        return NULL;
    Expr* p = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;

    Expr** si = NULL; size_t n = 0;
    PMSubMap lin, lout;
    Expr** vars = NULL;
    Expr** l = NULL; Expr* q = NULL;
    if (build_pipeline_state(p, x, &si, &n, &lin, &lout,
                              &vars, &l, &q) != 0) {
        return NULL;
    }

    Expr* p_conv = convert_to_tan(p, x);
    if (!p_conv) {
        free_pipeline_state(si, n, &lin, &lout, vars, l, q);
        return NULL;
    }
    Expr* lin_rules = pm_sub_map_to_rule_list(&lin);
    Expr* p_sub = eval_expand(eval_replace_all(p_conv, lin_rules));
    expr_free(p_conv);

    Expr* s = NULL;
    Expr* h = NULL;
    int ok = split_factor(p_sub, vars, l, n, &s, &h);
    expr_free(p_sub);
    if (ok != 0 || !s || !h) {
        if (s) expr_free(s);
        if (h) expr_free(h);
        free_pipeline_state(si, n, &lin, &lout, vars, l, q);
        return NULL;
    }

    Expr* lout_rules = pm_sub_map_to_rule_list(&lout);
    Expr* s_back = eval_replace_all(s, expr_copy(lout_rules));
    Expr* h_back = eval_replace_all(h, lout_rules);
    expr_free(s); expr_free(h);

    Expr* result_items[2] = { s_back, h_back };
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), result_items, 2);
    free_pipeline_state(si, n, &lin, &lout, vars, l, q);
    return result;
}

/* Integrate`Helpers`PMDeflation[p_in_atoms, x] — returns deflation(p)
 * in atom form. */
static Expr* builtin_pm_deflation(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)
        return NULL;
    Expr* p = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;

    Expr** si = NULL; size_t n = 0;
    PMSubMap lin, lout;
    Expr** vars = NULL;
    Expr** l = NULL; Expr* q = NULL;
    if (build_pipeline_state(p, x, &si, &n, &lin, &lout,
                              &vars, &l, &q) != 0) {
        return NULL;
    }

    Expr* p_conv = convert_to_tan(p, x);
    Expr* lin_rules = pm_sub_map_to_rule_list(&lin);
    Expr* p_sub = eval_expand(eval_replace_all(p_conv, lin_rules));
    expr_free(p_conv);

    Expr* defl = deflation(p_sub, vars, l, n);
    expr_free(p_sub);

    Expr* lout_rules = pm_sub_map_to_rule_list(&lout);
    Expr* result = eval_replace_all(defl, lout_rules);
    expr_free(defl);
    free_pipeline_state(si, n, &lin, &lout, vars, l, q);
    return result;
}

/* Integrate`Helpers`PMEnumerateMonoms[{v1, v2, ...}, d] — list of
 * monomials of total degree ≤ d in the given variables. */
static Expr* builtin_pm_enumerate_monoms(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)
        return NULL;
    Expr* vlist = res->data.function.args[0];
    Expr* dexpr = res->data.function.args[1];
    if (vlist->type != EXPR_FUNCTION
        || vlist->data.function.head->type != EXPR_SYMBOL
        || vlist->data.function.head->data.symbol != SYM_List)
        return NULL;
    if (dexpr->type != EXPR_INTEGER) return NULL;
    int total_degree = (int)dexpr->data.integer;
    if (total_degree < 0) return NULL;

    size_t nv = vlist->data.function.arg_count;
    Expr** vars = (Expr**)malloc(sizeof(Expr*) * (nv ? nv : 1));
    for (size_t i = 0; i < nv; i++) vars[i] = vlist->data.function.args[i];

    Expr** monoms = NULL;
    size_t n = 0;
    int err = enumerate_monomials(vars, nv, total_degree, &monoms, &n);
    free(vars);
    if (err) return NULL;

    Expr* list = expr_new_function(expr_new_symbol(SYM_List), monoms, n);
    free(monoms);
    return list;
}

/* Integrate`Helpers`PMSubstMap[f, x] — returns the forward
 * substitution list `lin` (term -> fresh var) as a List of Rules.
 * Convenient for verifying the substitution round-trip in tests. */
static Expr* builtin_pm_subst_map(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)
        return NULL;
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;

    Expr** si = NULL;
    size_t n = 0;
    if (collect_indets_closed(f, x, &si, &n) != 0) {
        for (size_t i = 0; i < n; i++) expr_free(si[i]);
        free(si);
        return NULL;
    }
    PMSubMap lin, lout;
    Expr** vars = NULL;
    if (build_substitution_maps(si, n, &lin, &lout, &vars, 1) != 0) {
        for (size_t i = 0; i < n; i++) expr_free(si[i]);
        free(si);
        return NULL;
    }
    Expr* result = pm_sub_map_to_rule_list(&lin);
    pm_sub_map_free(&lin);
    pm_sub_map_free(&lout);
    for (size_t i = 0; i < n; i++) expr_free(vars[i]);
    free(vars);
    for (size_t i = 0; i < n; i++) expr_free(si[i]);
    free(si);
    return result;
}

/* ------------------------------------------------------------------ */
/* Public entry: Integrate`RischNorman[f, x].                           */
/*                                                                      */
/* Phase 4: wires the full pmint pipeline (no log candidates yet).     */
/* ------------------------------------------------------------------ */

/* Wall-clock timeout machinery.  pmint's candidate space can grow
 * fast; we cap each invocation at PMINT_WALL_CLOCK_SEC seconds and
 * bail with NULL if it overshoots.  Uses SIGALRM + sigsetjmp so the
 * deep recursion in apply_d / split_factor / etc. can be unwound
 * safely. */
static sigjmp_buf pmint_timeout_env;
static volatile sig_atomic_t pmint_timed_out = 0;

static void pmint_timeout_handler(int sig) {
    (void)sig;
    pmint_timed_out = 1;
    siglongjmp(pmint_timeout_env, 1);
}

Expr* builtin_rischnorman(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) {
        return NULL;
    }
    Expr* f = res->data.function.args[0];
    Expr* x = res->data.function.args[1];
    if (x->type != EXPR_SYMBOL) return NULL;

    /* Wall-clock budget.  Install a SIGALRM handler, arm a one-shot
     * setitimer, and run rischnorman_integrate inside sigsetjmp.  On
     * timeout we jump back and return NULL — the dispatcher then
     * bubbles back unevaluated. */
    struct sigaction sa, old_sa;
    struct itimerval old_it;
    memset(&sa, 0, sizeof(sa));
    sa.sa_handler = pmint_timeout_handler;
    sigemptyset(&sa.sa_mask);
    sigaction(SIGALRM, &sa, &old_sa);

    struct itimerval new_it;
    memset(&new_it, 0, sizeof(new_it));
    new_it.it_value.tv_sec = PMINT_WALL_CLOCK_SEC;
    setitimer(ITIMER_REAL, &new_it, &old_it);

    pmint_timed_out = 0;
    Expr* result = NULL;
    if (sigsetjmp(pmint_timeout_env, 1) == 0) {
        result = rischnorman_integrate(f, x);
    }

    /* Disarm timer + restore handler. */
    struct itimerval disarm;
    memset(&disarm, 0, sizeof(disarm));
    setitimer(ITIMER_REAL, &disarm, NULL);
    setitimer(ITIMER_REAL, &old_it, NULL);
    sigaction(SIGALRM, &old_sa, NULL);

    if (pmint_timed_out) {
        /* Result may be partially constructed inside Maple-callees —
         * we don't have a safe way to clean up, but the next gc cycle
         * will catch the unreachable allocations.  Bubble back NULL. */
        return NULL;
    }
    return result;
}

/* ------------------------------------------------------------------ */
/* Init.                                                                */
/* ------------------------------------------------------------------ */

static void install(const char* name, Expr* (*fn)(Expr*), const char* docstring) {
    symtab_add_builtin(name, fn);
    symtab_get_def(name)->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    if (docstring) symtab_set_docstring(name, docstring);
}

void intrischnorman_init(void) {
    install("Integrate`RischNorman",
            builtin_rischnorman,
            "Integrate`RischNorman[f, x] applies the parallel-Risch /\n"
            "Risch-Norman heuristic (Bronstein's Poor Man's Integrator,\n"
            "2004) to find an elementary antiderivative of f in x.\n"
            "Returns the antiderivative on success or the call unevaluated\n"
            "on failure.  Called by the Integrate dispatcher as a\n"
            "fall-through after Integrate`IntegrateRational declines a\n"
            "non-rational integrand.");

    install("Integrate`Helpers`PMConvertToTan",
            builtin_pm_convert_to_tan,
            "Integrate`Helpers`PMConvertToTan[f, x] rewrites Sin / Cos /\n"
            "Sec / Csc and their hyperbolic siblings as rational functions\n"
            "of Tan[u/2] / Tanh[u/2].  Mirrors Maple's convert(_, tan).\n"
            "pmint.maple line 5.");

    install("Integrate`Helpers`PMPythagoreanRewrite",
            builtin_pm_pythagorean,
            "Integrate`Helpers`PMPythagoreanRewrite[f] rewrites Sec[u]^2,\n"
            "Csc[u]^2, Sech[u]^2, Csch[u]^2 via the Pythagorean\n"
            "identities so the result stays in the {Tan, Cot, Tanh,\n"
            "Coth} field generators.");

    install("Integrate`Helpers`PMSincosToTan",
            builtin_pm_sincos_to_tan,
            "Integrate`Helpers`PMSincosToTan[f, x] is the partial\n"
            "Weierstrass rewrite — only Sin / Cos / Sec / Csc and their\n"
            "hyperbolic siblings are converted to Tan[u/2] / Tanh[u/2]\n"
            "rationals; Tan / Cot / Tanh / Coth are left intact.  Used\n"
            "by the post-hoc differentiate-and-cancel verifier so that\n"
            "pmint's Tan-form output can be reconciled with the\n"
            "original Sin/Cos integrand without nesting Tan[x/4] etc.");

    install("Integrate`Helpers`PMCollectIndets",
            builtin_pm_collect_indets,
            "Integrate`Helpers`PMCollectIndets[f, x] returns the list of\n"
            "transcendental atoms in f (Tan, Tanh, Log, Exp, ArcTan, ...)\n"
            "whose derivative wrt x is non-zero, closed under one round\n"
            "of differentiation.  The first element is x itself.\n"
            "pmint.maple lines 6-8.");

    install("Integrate`Helpers`PMSubstMap",
            builtin_pm_subst_map,
            "Integrate`Helpers`PMSubstMap[f, x] returns the forward\n"
            "substitution rules `term -> pmint$v_i` as a List of Rule[].\n"
            "Used by pmint to lift the differential field to a polynomial\n"
            "ring over fresh variables.  pmint.maple lines 9-10.");

    /* Phase 3 — vector field / split_factor / deflation / monomials. */
    install("Integrate`Helpers`PMVectorField",
            builtin_pm_vector_field,
            "Integrate`Helpers`PMVectorField[f, x] returns\n"
            "{vars, l_vec, q} for the integrand f in x: vars is the list\n"
            "of fresh symbols, l_vec is the scaled derivative vector\n"
            "(l[k] = q · D[atom_k, x] /. lin), and q is the common\n"
            "denominator.  pmint.maple lines 11-15.");

    install("Integrate`Helpers`PMApplyD",
            builtin_pm_apply_d,
            "Integrate`Helpers`PMApplyD[f, x] applies the total\n"
            "derivation operator d (defined by the integrand's\n"
            "differential field) to f and returns the result back in\n"
            "original-atom form.  d(f) = Σ l[k] · ∂f/∂vars[k], where\n"
            "{vars, l} = Integrate`Helpers`PMVectorField[f, x].");

    install("Integrate`Helpers`PMSplitFactor",
            builtin_pm_split_factor,
            "Integrate`Helpers`PMSplitFactor[p, x] returns {s, h} such\n"
            "that p = s · h, where s captures the factors whose total\n"
            "derivation has non-zero gcd with the factor itself (the\n"
            "'normal' part) and h is the 'special' part.  pmint.maple\n"
            "lines 80-90.");

    install("Integrate`Helpers`PMDeflation",
            builtin_pm_deflation,
            "Integrate`Helpers`PMDeflation[p, x] returns the deflation\n"
            "of p under the total-derivation operator: the recursive\n"
            "PolynomialGCD of p with its ordinary derivative through\n"
            "content / primitive recursion.  pmint.maple lines 92-98.");

    install("Integrate`Helpers`PMEnumerateMonoms",
            builtin_pm_enumerate_monoms,
            "Integrate`Helpers`PMEnumerateMonoms[{v1, v2, ...}, d]\n"
            "enumerates all monomials of total degree ≤ d in the given\n"
            "variables.  Capped at PMINT_MAX_MONOMIALS (= 5000).\n"
            "pmint.maple lines 69-78.");
}
