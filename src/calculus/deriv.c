/*
 * deriv.c -- Native C implementation of Mathematica-style differentiation.
 *
 * This module replaces the fragile rule-based bootstrap in
 * src/internal/deriv.m with a direct, dispatch-driven implementation.
 *
 * Overview
 * --------
 * The two key entry points are the builtins D (partial derivative) and
 * Dt (total derivative). Both ultimately funnel through a single
 * recursive core, ``compute_deriv``, parameterised by an optional
 * differentiation variable. When the variable is non-NULL we compute a
 * partial derivative treating everything else as constant (using a fast
 * FreeQ-style walk to short-circuit constant sub-trees). When the
 * variable is NULL we compute a total derivative -- unknown symbols
 * then participate as ``Dt[sym]`` terms.
 *
 * Why this is faster than the rule-based implementation
 * -----------------------------------------------------
 * The old deriv.m relied on ~60 DownValues. Each call to D[f, x] would:
 *   * scan the DownValues list for D linearly,
 *   * attempt pattern matching against every rule head (Plus, Times,
 *     Power, every elementary function, ...),
 *   * run ``/;`` side-conditions such as FreeQ,
 *   * perform attempt-evaluate/backtrack cycles in the matcher,
 *   * recursively re-evaluate the result through the full rule engine.
 *
 * In contrast, this module performs a single head-symbol strcmp dispatch
 * per call, constructs the derivative expression directly, and lets the
 * outer evaluator simplify arithmetic. Crucially, the constant-detection
 * step uses a tailored structural traversal (expr_free_of) that avoids
 * calling out to the generic FreeQ builtin.
 *
 * Returned expressions
 * --------------------
 * Every builder below produces plain un-reduced expression trees (e.g.
 * Plus[0, x] or Times[1, x]). The outer Mathilda evaluator runs a full
 * fixed-point reduction on the value we return, so Plus[0, ...],
 * Times[1, ...], and all subsequent chain-rule simplifications fold
 * automatically. This keeps the code readable and avoids duplicating
 * the arithmetic simplifier.
 *
 * Memory ownership
 * ----------------
 * Every helper that returns an ``Expr*`` returns a freshly allocated
 * tree owned by the caller. Input expressions are never mutated;
 * sub-expressions that need to be reused are always deep-copied.
 */

#include "deriv.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "sym_intern.h"
#include "print.h"
#include "series.h"

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parse.h"   /* parse_expression: the Jacobi derivative templates */

/* ---------------------------------------------------------------------- */
/* Tiny expression builders                                                */
/* ---------------------------------------------------------------------- */

static Expr* mk_int(int64_t v) { return expr_new_integer(v); }
static Expr* mk_sym(const char* s) { return expr_new_symbol(s); }

/* Build a Function expression with the given head (already owned) and
 * a stack-provided argument array. Arguments are moved (ownership
 * transferred) into the new node. */
static Expr* mk_fn_take(Expr* head, Expr** args, size_t n) {
    /* expr_new_function memcpys the pointer slice, so a local stack
     * buffer is enough. */
    return expr_new_function(head, args, n);
}

static Expr* mk_fn1(const char* name, Expr* a) {
    Expr* args[1] = { a };
    return mk_fn_take(mk_sym(name), args, 1);
}

static Expr* mk_fn2(const char* name, Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    return mk_fn_take(mk_sym(name), args, 2);
}

/* Take ownership of `items` pointers and wrap them as name[items...].
 * The `items` array itself is freed here. */
static Expr* mk_fnN_adopt(const char* name, Expr** items, size_t n) {
    Expr* r = expr_new_function(mk_sym(name), items, n);
    free(items);
    return r;
}

/* Like mk_fnN_adopt but the head is an arbitrary pre-built expression
 * (e.g. Derivative[k][f] when producing Derivative[k][f][g]). The
 * `items` pointer must be a heap-allocated array owned by the caller;
 * this helper frees it. */
static Expr* mk_fn_head_adopt(Expr* head, Expr** items, size_t n) {
    Expr* r = expr_new_function(head, items, n);
    free(items);
    return r;
}

/* Wrap a single expression `arg` under an arbitrary pre-built head,
 * producing `head[arg]`. Used to build things like Derivative[1][f]
 * and Derivative[1][f][g]. */
static Expr* mk_fn_head1(Expr* head, Expr* arg) {
    Expr* items[1] = { arg };
    return expr_new_function(head, items, 1);
}

/* Negation helper: build Times[-1, a]. */
static Expr* mk_neg(Expr* a) {
    return mk_fn2("Times", mk_int(-1), a);
}

/* ---------------------------------------------------------------------- */
/* Predicates                                                              */
/* ---------------------------------------------------------------------- */

static bool is_sym(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, name) == 0;
}

static bool is_lit_zero(const Expr* e) {
    return e && e->type == EXPR_INTEGER && e->data.integer == 0;
}

/* True if ``f`` contains no subexpression structurally equal to ``x``.
 *
 * This is a hot path: it is called once per recursive deriv step to
 * short-circuit constant sub-trees. We avoid expr_hash and the generic
 * FreeQ builtin to keep the inner loop allocation-free. */
static bool expr_free_of(const Expr* f, const Expr* x) {
    if (expr_eq(f, x)) return false;
    if (f->type == EXPR_FUNCTION) {
        if (!expr_free_of(f->data.function.head, x)) return false;
        for (size_t i = 0; i < f->data.function.arg_count; i++) {
            if (!expr_free_of(f->data.function.args[i], x)) return false;
        }
    }
    return true;
}

/* True if ``sym`` (an EXPR_SYMBOL) is listed in ``nonconsts`` (a
 * List[...] of symbols, or NULL). Used to recognise variables that the
 * caller has declared implicitly dependent on the differentiation
 * variable via the NonConstants option. */
static bool nonconsts_contains(const Expr* nonconsts, const Expr* sym) {
    if (!nonconsts || !sym) return false;
    if (sym->type != EXPR_SYMBOL) return false;
    if (nonconsts->type != EXPR_FUNCTION) return false;
    if (!is_sym(nonconsts->data.function.head, "List")) return false;
    size_t n = nonconsts->data.function.arg_count;
    for (size_t i = 0; i < n; i++) {
        const Expr* v = nonconsts->data.function.args[i];
        if (v && v->type == EXPR_SYMBOL && v->data.symbol.name == sym->data.symbol.name) {
            return true;
        }
    }
    return false;
}

/* True if ``f`` contains any symbol listed in ``nonconsts``. With
 * NonConstants set we must NOT short-circuit constant subtrees that
 * happen to mention an implicitly-dependent symbol -- otherwise
 * D[2 y, x, NonConstants -> {y}] would wrongly fold to 0. */
static bool expr_contains_nonconst(const Expr* f, const Expr* nonconsts) {
    if (!nonconsts || nonconsts->type != EXPR_FUNCTION) return false;
    if (!is_sym(nonconsts->data.function.head, "List")) return false;
    size_t n = nonconsts->data.function.arg_count;
    for (size_t i = 0; i < n; i++) {
        if (!expr_free_of(f, nonconsts->data.function.args[i])) return true;
    }
    return false;
}

/* Build the canonical unevaluated form D[sym, x, NonConstants -> {...}].
 * ``nonconsts`` is the already-canonicalised List of implicit-dependence
 * symbols; we deep-copy it into the result so the caller retains
 * ownership of the original. */
static Expr* build_unevaluated_d_nonconsts(const Expr* sym, const Expr* x,
                                           const Expr* nonconsts) {
    Expr* opt = mk_fn2("Rule",
                       mk_sym("NonConstants"),
                       expr_copy((Expr*)nonconsts));
    Expr** args = malloc(sizeof(Expr*) * 3);
    args[0] = expr_copy((Expr*)sym);
    args[1] = expr_copy((Expr*)x);
    args[2] = opt;
    Expr* r = expr_new_function(mk_sym("D"), args, 3);
    free(args);
    return r;
}

/* Recognise the symbols that Dt treats as absolute constants. These
 * mirror the rule set in the original deriv.m plus a few additional
 * Mathematica-compatible constants. */
static bool symbol_has_constant_attr(const char* name) {
    SymbolDef* d = symtab_lookup(name);
    return d && (d->attributes & ATTR_CONSTANT);
}

static bool is_dt_constant_symbol(const Expr* e) {
    if (!e || e->type != EXPR_SYMBOL) return false;
    const char* s = e->data.symbol.name;
    return s == SYM_Pi || s == SYM_E || s == SYM_I ||
           s == SYM_Infinity || s == SYM_ComplexInfinity ||
           s == SYM_EulerGamma || s == SYM_Catalan ||
           s == SYM_GoldenRatio || s == SYM_Degree ||
           s == SYM_GoldenAngle || s == SYM_Glaisher ||
           s == SYM_Khinchin ||
           s == SYM_True || s == SYM_False ||
           /* A user symbol given the Constant attribute
            * (SetAttributes[R, Constant]) is a constant for Dt: Dt[R] = 0,
            * Dt[R x] = R Dt[x]. */
           symbol_has_constant_attr(s);
}

/* ---------------------------------------------------------------------- */
/* Derivative of an elementary unary function                              */
/* ---------------------------------------------------------------------- */

/* Given the head name of a known elementary single-argument function F
 * and its argument g, produce F'(g). Returns NULL if F is not
 * recognised.
 *
 * The argument `g` is NOT consumed -- we copy it as needed.
 */
static Expr* elementary_fprime(const char* name, Expr* g) {
    /* --- trigonometric --- */
    if (name == SYM_Sin) return mk_fn1("Cos", expr_copy(g));
    if (name == SYM_Cos) return mk_neg(mk_fn1("Sin", expr_copy(g)));
    if (name == SYM_Tan) return mk_fn2("Power", mk_fn1("Sec", expr_copy(g)), mk_int(2));
    if (name == SYM_Cot) return mk_neg(mk_fn2("Power", mk_fn1("Csc", expr_copy(g)), mk_int(2)));
    if (name == SYM_Sec) return mk_fn2("Times", mk_fn1("Sec", expr_copy(g)), mk_fn1("Tan", expr_copy(g)));
    if (name == SYM_Csc) return mk_neg(mk_fn2("Times", mk_fn1("Csc", expr_copy(g)), mk_fn1("Cot", expr_copy(g))));

    /* --- inverse trigonometric --- */
    /* D[ArcSin[g], g] = 1 / Sqrt[1 - g^2]; we encode 1/Sqrt[u] as Power[u, -1/2]. */
    if (name == SYM_ArcSin) {
        Expr* one_minus_g2 = mk_fn2("Plus", mk_int(1),
                                    mk_neg(mk_fn2("Power", expr_copy(g), mk_int(2))));
        return mk_fn2("Power", mk_fn1("Sqrt", one_minus_g2), mk_int(-1));
    }
    if (name == SYM_ArcCos) {
        Expr* one_minus_g2 = mk_fn2("Plus", mk_int(1),
                                    mk_neg(mk_fn2("Power", expr_copy(g), mk_int(2))));
        return mk_neg(mk_fn2("Power", mk_fn1("Sqrt", one_minus_g2), mk_int(-1)));
    }
    if (name == SYM_ArcTan) {
        Expr* one_plus_g2 = mk_fn2("Plus", mk_int(1), mk_fn2("Power", expr_copy(g), mk_int(2)));
        return mk_fn2("Power", one_plus_g2, mk_int(-1));
    }
    if (name == SYM_ArcCot) {
        Expr* one_plus_g2 = mk_fn2("Plus", mk_int(1), mk_fn2("Power", expr_copy(g), mk_int(2)));
        return mk_neg(mk_fn2("Power", one_plus_g2, mk_int(-1)));
    }
    if (name == SYM_ArcSec) {
        /* 1 / (g^2 * Sqrt[1 - 1/g^2]) */
        Expr* g2 = mk_fn2("Power", expr_copy(g), mk_int(2));
        Expr* inv_g2 = mk_fn2("Power", expr_copy(g), mk_int(-2));
        Expr* radicand = mk_fn2("Plus", mk_int(1), mk_neg(inv_g2));
        Expr* denom = mk_fn2("Times", g2, mk_fn1("Sqrt", radicand));
        return mk_fn2("Power", denom, mk_int(-1));
    }
    if (name == SYM_ArcCsc) {
        Expr* g2 = mk_fn2("Power", expr_copy(g), mk_int(2));
        Expr* inv_g2 = mk_fn2("Power", expr_copy(g), mk_int(-2));
        Expr* radicand = mk_fn2("Plus", mk_int(1), mk_neg(inv_g2));
        Expr* denom = mk_fn2("Times", g2, mk_fn1("Sqrt", radicand));
        return mk_neg(mk_fn2("Power", denom, mk_int(-1)));
    }

    /* --- hyperbolic --- */
    if (name == SYM_Sinh) return mk_fn1("Cosh", expr_copy(g));
    if (name == SYM_Cosh) return mk_fn1("Sinh", expr_copy(g));
    if (name == SYM_Tanh) return mk_fn2("Power", mk_fn1("Sech", expr_copy(g)), mk_int(2));
    if (name == SYM_Coth) return mk_neg(mk_fn2("Power", mk_fn1("Csch", expr_copy(g)), mk_int(2)));
    if (name == SYM_Sech) return mk_neg(mk_fn2("Times", mk_fn1("Sech", expr_copy(g)), mk_fn1("Tanh", expr_copy(g))));
    if (name == SYM_Csch) return mk_neg(mk_fn2("Times", mk_fn1("Csch", expr_copy(g)), mk_fn1("Coth", expr_copy(g))));

    /* --- inverse hyperbolic --- */
    if (name == SYM_ArcSinh) {
        Expr* radicand = mk_fn2("Plus", mk_int(1), mk_fn2("Power", expr_copy(g), mk_int(2)));
        return mk_fn2("Power", mk_fn1("Sqrt", radicand), mk_int(-1));
    }
    if (name == SYM_ArcCosh) {
        Expr* radicand = mk_fn2("Plus", mk_fn2("Power", expr_copy(g), mk_int(2)), mk_int(-1));
        return mk_fn2("Power", mk_fn1("Sqrt", radicand), mk_int(-1));
    }
    if (name == SYM_ArcTanh) {
        Expr* denom = mk_fn2("Plus", mk_int(1), mk_neg(mk_fn2("Power", expr_copy(g), mk_int(2))));
        return mk_fn2("Power", denom, mk_int(-1));
    }
    if (name == SYM_ArcCoth) {
        Expr* denom = mk_fn2("Plus", mk_int(1), mk_neg(mk_fn2("Power", expr_copy(g), mk_int(2))));
        return mk_fn2("Power", denom, mk_int(-1));
    }
    if (name == SYM_ArcSech) {
        Expr* g2 = mk_fn2("Power", expr_copy(g), mk_int(2));
        Expr* inv_g2 = mk_fn2("Power", expr_copy(g), mk_int(-2));
        Expr* radicand = mk_fn2("Plus", inv_g2, mk_int(-1));
        Expr* denom = mk_fn2("Times", g2, mk_fn1("Sqrt", radicand));
        return mk_neg(mk_fn2("Power", denom, mk_int(-1)));
    }
    if (name == SYM_ArcCsch) {
        Expr* g2 = mk_fn2("Power", expr_copy(g), mk_int(2));
        Expr* inv_g2 = mk_fn2("Power", expr_copy(g), mk_int(-2));
        Expr* radicand = mk_fn2("Plus", mk_int(1), inv_g2);
        Expr* denom = mk_fn2("Times", g2, mk_fn1("Sqrt", radicand));
        return mk_neg(mk_fn2("Power", denom, mk_int(-1)));
    }

    /* --- Airy Ai: d/dg AiryAi[g] = AiryAiPrime[g]. --- */
    if (name == SYM_AiryAi) {
        return mk_fn1("AiryAiPrime", expr_copy(g));
    }

    /* --- Airy Ai': d/dg AiryAiPrime[g] = g AiryAi[g]  (from Ai'' = z Ai). --- */
    if (name == SYM_AiryAiPrime) {
        return mk_fn2("Times", expr_copy(g), mk_fn1("AiryAi", expr_copy(g)));
    }

    /* --- Airy Bi: d/dg AiryBi[g] = AiryBiPrime[g]. --- */
    if (name == SYM_AiryBi) {
        return mk_fn1("AiryBiPrime", expr_copy(g));
    }

    /* --- Airy Bi': d/dg AiryBiPrime[g] = g AiryBi[g]  (from Bi'' = z Bi). --- */
    if (name == SYM_AiryBiPrime) {
        return mk_fn2("Times", expr_copy(g), mk_fn1("AiryBi", expr_copy(g)));
    }

    /* --- Lambert W: d/dg ProductLog[g] = ProductLog[g] / (g (1 + ProductLog[g])). --- */
    if (name == SYM_ProductLog) {
        Expr* w = mk_fn1("ProductLog", expr_copy(g));
        Expr* one_plus_w = mk_fn2("Plus", mk_int(1), w);
        Expr* denom = mk_fn2("Times", expr_copy(g), one_plus_w);
        return mk_fn2("Times", mk_fn1("ProductLog", expr_copy(g)),
                      mk_fn2("Power", denom, mk_int(-1)));
    }

    /* --- sine integral: d/dg SinIntegral[g] = Sinc[g] = Sin[g]/g. --- */
    if (name == SYM_SinIntegral) {
        return mk_fn1("Sinc", expr_copy(g));
    }

    /* --- cosine integral: d/dg CosIntegral[g] = Cos[g]/g. --- */
    if (name == SYM_CosIntegral) {
        return mk_fn2("Times", mk_fn1("Cos", expr_copy(g)),
                      mk_fn2("Power", expr_copy(g), mk_int(-1)));
    }

    /* --- hyperbolic sine integral: d/dg SinhIntegral[g] = Sinh[g]/g.
     *     (No Sinhc builtin, so emit Sinh[g]/g directly.) --- */
    if (name == SYM_SinhIntegral) {
        return mk_fn2("Times", mk_fn1("Sinh", expr_copy(g)),
                      mk_fn2("Power", expr_copy(g), mk_int(-1)));
    }

    /* --- hyperbolic cosine integral: d/dg CoshIntegral[g] = Cosh[g]/g. --- */
    if (name == SYM_CoshIntegral) {
        return mk_fn2("Times", mk_fn1("Cosh", expr_copy(g)),
                      mk_fn2("Power", expr_copy(g), mk_int(-1)));
    }

    /* --- Fresnel integrals: d/dg FresnelC[g] = Cos[Pi g^2/2],
     *     d/dg FresnelS[g] = Sin[Pi g^2/2]. (Pi g^2/2 as Times[1/2, Pi, g^2].) */
    if (name == SYM_FresnelC || name == SYM_FresnelS) {
        Expr* half = mk_fn2("Power", mk_int(2), mk_int(-1));   /* 1/2 */
        Expr* g2   = mk_fn2("Power", expr_copy(g), mk_int(2));
        Expr* arg  = mk_fn2("Times", half, mk_fn2("Times", mk_sym("Pi"), g2));
        return mk_fn1(name == SYM_FresnelC ? "Cos" : "Sin", arg);
    }

    /* --- cardinal sine: d/dg Sinc[g] = Cos[g]/g - Sin[g]/g^2. --- */
    if (name == SYM_Sinc) {
        Expr* t1 = mk_fn2("Times", mk_fn1("Cos", expr_copy(g)),
                          mk_fn2("Power", expr_copy(g), mk_int(-1)));
        Expr* t2 = mk_neg(mk_fn2("Times", mk_fn1("Sin", expr_copy(g)),
                          mk_fn2("Power", expr_copy(g), mk_int(-2))));
        return mk_fn2("Plus", t1, t2);
    }

    /* --- complete elliptic integrals of the first and second kind, in the
     * PARAMETER convention m = k^2 (g is m here):
     *   d/dm EllipticK[m] = (EllipticE[m] - (1 - m) EllipticK[m]) / (2 m (1 - m))
     *   d/dm EllipticE[m] = (EllipticE[m] - EllipticK[m]) / (2 m)
     * Only the one-argument spellings reach this table, so EllipticE[phi, m]
     * (and EllipticF, EllipticPi) are handled in the multi-argument block. --- */
    if (name == SYM_EllipticK) {
        Expr* one_minus_m = mk_fn2("Plus", mk_int(1), mk_neg(expr_copy(g)));
        Expr* num = mk_fn2("Plus", mk_fn1("EllipticE", expr_copy(g)),
                        mk_neg(mk_fn2("Times", expr_copy(one_minus_m),
                                      mk_fn1("EllipticK", expr_copy(g)))));
        Expr* den = mk_fn2("Times", mk_fn2("Times", mk_int(2), expr_copy(g)),
                           one_minus_m);
        return mk_fn2("Times", num, mk_fn2("Power", den, mk_int(-1)));
    }
    if (name == SYM_EllipticE) {
        Expr* num = mk_fn2("Plus", mk_fn1("EllipticE", expr_copy(g)),
                        mk_neg(mk_fn1("EllipticK", expr_copy(g))));
        Expr* den = mk_fn2("Times", mk_int(2), expr_copy(g));
        return mk_fn2("Times", num, mk_fn2("Power", den, mk_int(-1)));
    }

    /* --- error function: d/dg Erf[g] = (2/Sqrt[Pi]) E^(-g^2). --- */
    if (name == SYM_Erf) {
        Expr* coeff = mk_fn2("Times", mk_int(2),
                             mk_fn2("Power", mk_fn1("Sqrt", mk_sym("Pi")), mk_int(-1)));
        Expr* gaussian = mk_fn1("Exp",
                             mk_neg(mk_fn2("Power", expr_copy(g), mk_int(2))));
        return mk_fn2("Times", coeff, gaussian);
    }

    /* --- imaginary error function: d/dg Erfi[g] = (2/Sqrt[Pi]) E^(g^2). --- */
    if (name == SYM_Erfi) {
        Expr* coeff = mk_fn2("Times", mk_int(2),
                             mk_fn2("Power", mk_fn1("Sqrt", mk_sym("Pi")), mk_int(-1)));
        Expr* gaussian = mk_fn1("Exp",
                             mk_fn2("Power", expr_copy(g), mk_int(2)));
        return mk_fn2("Times", coeff, gaussian);
    }

    /* --- Heaviside step: d/dg HeavisideTheta[g] = DiracDelta[g].  Needed so a
     * Green's-function / impulse solution differentiates cleanly (e.g. the IC
     * fit of a causal impulse response), matching Mathematica. --- */
    if (name == SYM_HeavisideTheta) return mk_fn1("DiracDelta", expr_copy(g));

    /* --- complementary error function: d/dg Erfc[g] = -(2/Sqrt[Pi]) E^(-g^2). --- */
    if (name == SYM_Erfc) {
        Expr* coeff = mk_fn2("Times", mk_int(-2),
                             mk_fn2("Power", mk_fn1("Sqrt", mk_sym("Pi")), mk_int(-1)));
        Expr* gaussian = mk_fn1("Exp",
                             mk_neg(mk_fn2("Power", expr_copy(g), mk_int(2))));
        return mk_fn2("Times", coeff, gaussian);
    }

    /* --- inverse error function: d/dg InverseErf[g] = (Sqrt[Pi]/2) E^(InverseErf[g]^2). --- */
    if (name == SYM_InverseErf) {
        Expr* coeff = mk_fn2("Times", mk_fn1("Sqrt", mk_sym("Pi")),
                             mk_fn2("Power", mk_int(2), mk_int(-1)));
        Expr* sq = mk_fn2("Power", mk_fn1("InverseErf", expr_copy(g)), mk_int(2));
        return mk_fn2("Times", coeff, mk_fn1("Exp", sq));
    }

    /* --- inverse complementary error function:
     *     d/dg InverseErfc[g] = -(Sqrt[Pi]/2) E^(InverseErfc[g]^2). --- */
    if (name == SYM_InverseErfc) {
        Expr* coeff = mk_fn2("Times", mk_int(-1),
                             mk_fn2("Times", mk_fn1("Sqrt", mk_sym("Pi")),
                                    mk_fn2("Power", mk_int(2), mk_int(-1))));
        Expr* sq = mk_fn2("Power", mk_fn1("InverseErfc", expr_copy(g)), mk_int(2));
        return mk_fn2("Times", coeff, mk_fn1("Exp", sq));
    }

    /* --- exponential integral: d/dg ExpIntegralEi[g] = E^g / g. --- */
    if (name == SYM_ExpIntegralEi) {
        return mk_fn2("Times", mk_fn1("Exp", expr_copy(g)),
                      mk_fn2("Power", expr_copy(g), mk_int(-1)));
    }

    /* --- logarithmic integral: d/dg LogIntegral[g] = 1 / Log[g]. --- */
    if (name == SYM_LogIntegral) {
        return mk_fn2("Power", mk_fn1("Log", expr_copy(g)), mk_int(-1));
    }

    /* --- exp/log and sqrt --- */
    if (name == SYM_Exp) return mk_fn1("Exp", expr_copy(g));
    if (name == SYM_Log) return mk_fn2("Power", expr_copy(g), mk_int(-1));
    if (name == SYM_Sqrt) {
        /* d/dg Sqrt[g] = 1/(2 Sqrt[g]) */
        return mk_fn2("Power", mk_fn2("Times", mk_int(2), mk_fn1("Sqrt", expr_copy(g))), mk_int(-1));
    }

    /* --- Fibonacci number: derivative w.r.t. its single argument. ---
     *   d/dn Fibonacci[n] = (1/Sqrt[5]) GoldenRatio^-n
     *       (GoldenRatio^(2 n) Log[GoldenRatio] + Cos[n Pi] Log[GoldenRatio]
     *        + Pi Sin[n Pi]). */
    if (name == SYM_Fibonacci) {
        Expr* logphi = mk_fn1("Log", mk_sym("GoldenRatio"));
        Expr* t1 = mk_fn2("Times",
                          mk_fn2("Power", mk_sym("GoldenRatio"),
                                 mk_fn2("Times", mk_int(2), expr_copy(g))),
                          expr_copy(logphi));
        Expr* t2 = mk_fn2("Times",
                          mk_fn1("Cos", mk_fn2("Times", expr_copy(g), mk_sym("Pi"))),
                          logphi);                          /* consumes logphi */
        Expr* t3 = mk_fn2("Times", mk_sym("Pi"),
                          mk_fn1("Sin", mk_fn2("Times", expr_copy(g), mk_sym("Pi"))));
        Expr** terms = malloc(sizeof(Expr*) * 3);
        terms[0] = t1; terms[1] = t2; terms[2] = t3;
        Expr* inner = mk_fnN_adopt("Plus", terms, 3);
        Expr* pref = mk_fn2("Times",
                            mk_fn2("Power", mk_fn1("Sqrt", mk_int(5)), mk_int(-1)),
                            mk_fn2("Power", mk_sym("GoldenRatio"),
                                   mk_fn2("Times", mk_int(-1), expr_copy(g))));
        return mk_fn2("Times", pref, inner);
    }

    /* --- LucasL number: derivative w.r.t. its single argument. ---
     *   d/dn LucasL[n] = GoldenRatio^-n (GoldenRatio^(2 n) Log[GoldenRatio]
     *       - Cos[n Pi] Log[GoldenRatio] - Pi Sin[n Pi]). */
    if (name == SYM_LucasL) {
        Expr* logphi = mk_fn1("Log", mk_sym("GoldenRatio"));
        Expr* t1 = mk_fn2("Times",
                          mk_fn2("Power", mk_sym("GoldenRatio"),
                                 mk_fn2("Times", mk_int(2), expr_copy(g))),
                          expr_copy(logphi));
        Expr* t2 = mk_neg(mk_fn2("Times",
                          mk_fn1("Cos", mk_fn2("Times", expr_copy(g), mk_sym("Pi"))),
                          logphi));                         /* consumes logphi */
        Expr* t3 = mk_neg(mk_fn2("Times", mk_sym("Pi"),
                          mk_fn1("Sin", mk_fn2("Times", expr_copy(g), mk_sym("Pi")))));
        Expr** terms = malloc(sizeof(Expr*) * 3);
        terms[0] = t1; terms[1] = t2; terms[2] = t3;
        Expr* inner = mk_fnN_adopt("Plus", terms, 3);
        Expr* pref = mk_fn2("Power", mk_sym("GoldenRatio"),
                            mk_fn2("Times", mk_int(-1), expr_copy(g)));
        return mk_fn2("Times", pref, inner);
    }

    return NULL; /* not a recognised elementary unary */
}

/* ---------------------------------------------------------------------- */
/* Core recursive derivative                                               */
/* ---------------------------------------------------------------------- */

/* Total-derivative mode for ``Dt[f, x]``. When set, ``compute_deriv``
 * (always called with a non-NULL ``x``) treats every symbol other than
 * ``x`` and the distinguished constants as an *implicit* function of
 * ``x``: such a symbol ``s`` contributes the unevaluated total
 * derivative ``Dt[s, x]`` rather than ``0`` (the ``D`` convention).
 * Numbers and Dt-constants still vanish. The flag is set only for the
 * duration of a ``Dt[f, x]`` evaluation and saved/restored so nested
 * re-entrant differentiation behaves correctly. */
static bool g_dt_total_mode = false;

/* Forward declarations. ``nonconsts`` may be NULL (no NonConstants
 * option in effect) or a fully canonicalised List[sym1, ...] of
 * symbols treated as implicit functions of x. */
static Expr* compute_deriv(Expr* f, Expr* x, Expr* nonconsts);

/* Shortcut: return D[g, x] as a fresh tree. When x is NULL, return
 * Dt[g] instead. For symbols and numeric atoms the answer is folded
 * immediately; for compound expressions we recurse through
 * compute_deriv so constants short-circuit there too.
 *
 * NEVER returns NULL. compute_deriv returns NULL to mean "I cannot
 * differentiate this", which is a valid answer only at the TOP level,
 * where builtin_d / builtin_dt keep the whole call unevaluated and
 * higher_order_partial wraps the residue as an inert D[...]. Every
 * caller of THIS function, by contrast, stores the result straight into
 * an argument slot of a larger tree (the sum rule, the product rule,
 * the chain rule, ...) -- and a NULL in an argument array is
 * dereferenced unconditionally by evaluate_step, i.e. a segfault rather
 * than an unevaluated result. There are 82 such call sites and none of
 * them checks, so the guarantee belongs here, once, instead of being
 * restated 82 times. Emitting the inert derivative also gives the
 * better answer: D[Inactive[Integrate][f[u], Sqrt[x]] + 1, x] becomes
 * D[Inactive[Integrate][f[u], Sqrt[x]], x] -- the inert derivative of
 * the inert part, with the differentiable part still differentiated --
 * which is what Mathematica returns. */
static Expr* deriv_of(Expr* g, Expr* x, Expr* nonconsts) {
    Expr* d = compute_deriv(g, x, nonconsts);
    if (d) return d;
    return x ? mk_fn2("D", expr_copy(g), expr_copy(x))
             : mk_fn1("Dt", expr_copy(g));
}

/* The chain rule applied to an ``f[g1, g2, ..., gn]`` expression whose
 * head is not one of the explicitly-handled elementary heads. Produces
 *
 *   Sum_{k : D[gk, x] != 0} Derivative[0..1_k..0][f][g1..gn] * D[gk, x]
 *
 * If all partials vanish the result is 0. If exactly one contributes,
 * the outer Plus is elided. */
static Expr* chain_rule_unknown(Expr* f, Expr* x, Expr* nonconsts) {
    Expr* head = f->data.function.head;
    size_t n = f->data.function.arg_count;
    Expr** args = f->data.function.args;

    /* Collect (nonzero-derivative, partial-index) terms. */
    Expr** terms = malloc(sizeof(Expr*) * n);
    size_t nterms = 0;
    for (size_t k = 0; k < n; k++) {
        Expr* dk = deriv_of(args[k], x, nonconsts);
        if (is_lit_zero(dk)) { expr_free(dk); continue; }

        /* Build the indicator Derivative[0,..,1,..,0]. */
        Expr** idx = malloc(sizeof(Expr*) * n);
        for (size_t j = 0; j < n; j++) idx[j] = mk_int(j == k ? 1 : 0);
        Expr* op = mk_fnN_adopt("Derivative", idx, n);          /* Derivative[...]   */
        Expr* op_f = mk_fn_head1(op, expr_copy(head));          /* Derivative[...][f] */

        /* Apply to the original argument list. */
        Expr** gcopy = malloc(sizeof(Expr*) * n);
        for (size_t j = 0; j < n; j++) gcopy[j] = expr_copy(args[j]);
        Expr* applied = mk_fn_head_adopt(op_f, gcopy, n);

        terms[nterms++] = mk_fn2("Times", applied, dk);
    }

    if (nterms == 0) { free(terms); return mk_int(0); }
    if (nterms == 1) { Expr* r = terms[0]; free(terms); return r; }
    return mk_fnN_adopt("Plus", terms, nterms);
}

/* Differentiate an expression of the form
 *
 *     Derivative[n1, ..., nm][f] [ g1, ..., gn ]
 *
 * Each nonzero partial derivative advances the corresponding index in
 * the Derivative head. Returns NULL if the expression doesn't match
 * this shape (so the caller can try other dispatch rules). */
static Expr* deriv_of_derivative_form(Expr* f, Expr* x, Expr* nonconsts) {
    Expr* head = f->data.function.head;
    if (head->type != EXPR_FUNCTION) return NULL;
    if (head->data.function.arg_count != 1) return NULL;   /* need exactly one f */

    Expr* outer = head->data.function.head;                /* Derivative[n1..nm] */
    if (outer->type != EXPR_FUNCTION) return NULL;
    if (!is_sym(outer->data.function.head, "Derivative")) return NULL;

    size_t m = outer->data.function.arg_count;
    size_t n = f->data.function.arg_count;
    if (m != n) return NULL;                               /* arity mismatch */

    Expr* f0 = head->data.function.args[0];                /* the underlying f */
    Expr** args = f->data.function.args;

    Expr** terms = malloc(sizeof(Expr*) * n);
    size_t nterms = 0;
    for (size_t k = 0; k < n; k++) {
        Expr* dk = deriv_of(args[k], x, nonconsts);
        if (is_lit_zero(dk)) { expr_free(dk); continue; }

        /* Advance the k-th order index. If it's a plain integer, bump it;
         * otherwise wrap it in a Plus so the evaluator can simplify. */
        Expr** new_idx = malloc(sizeof(Expr*) * m);
        for (size_t j = 0; j < m; j++) {
            Expr* orig = outer->data.function.args[j];
            if (j == k) {
                if (orig->type == EXPR_INTEGER) {
                    new_idx[j] = mk_int(orig->data.integer + 1);
                } else {
                    new_idx[j] = mk_fn2("Plus", expr_copy(orig), mk_int(1));
                }
            } else {
                new_idx[j] = expr_copy(orig);
            }
        }
        Expr* new_op = mk_fnN_adopt("Derivative", new_idx, m);
        Expr* new_op_f = mk_fn_head1(new_op, expr_copy(f0));

        Expr** gcopy = malloc(sizeof(Expr*) * n);
        for (size_t j = 0; j < n; j++) gcopy[j] = expr_copy(args[j]);
        Expr* applied = mk_fn_head_adopt(new_op_f, gcopy, n);

        terms[nterms++] = mk_fn2("Times", applied, dk);
    }

    if (nterms == 0) { free(terms); return mk_int(0); }
    if (nterms == 1) { Expr* r = terms[0]; free(terms); return r; }
    return mk_fnN_adopt("Plus", terms, nterms);
}

/* Core dispatch. If `x` is non-NULL we are computing D[f, x]; otherwise
 * we are computing Dt[f]. The two paths share 95% of their logic --
 * only constant-handling and the "unknown symbol" base-case differ.
 * ``nonconsts`` (may be NULL) names symbols that should be treated as
 * implicit functions of x; these short-circuit neither expr_free_of
 * nor the symbol-base-case. */
static Expr* compute_deriv(Expr* f, Expr* x, Expr* nonconsts) {
    /* --- Partial-derivative base cases. --- */
    if (x) {
        if (expr_eq(f, x)) return mk_int(1);
        /* D[Indeterminate, x] = Indeterminate: an indeterminate quantity has
         * no well-defined derivative. Matches Mathematica and keeps higher
         * derivatives of UnitStep (whose first derivative carries an
         * Indeterminate Piecewise value) stable rather than collapsing the
         * value to 0. */
        if (f->type == EXPR_SYMBOL && f->data.symbol.name == SYM_Indeterminate)
            return expr_copy(f);
        /* A symbol declared as NonConstant produces the canonical
         * unevaluated implicit-derivative form. */
        if (f->type == EXPR_SYMBOL && nonconsts_contains(nonconsts, f)) {
            return build_unevaluated_d_nonconsts(f, x, nonconsts);
        }
        /* Dt[f, x] total-derivative mode: every symbol other than x and
         * the distinguished constants is an implicit function of x, so it
         * does NOT vanish -- its contribution is the unevaluated total
         * derivative Dt[sym, x]. Numbers and Dt-constants still vanish.
         * Compound expressions fall through to the dispatch below; we must
         * NOT apply the free-of-x short-circuit, or dependent sub-symbols
         * would be dropped. */
        if (g_dt_total_mode) {
            if (f->type == EXPR_INTEGER || f->type == EXPR_REAL ||
                f->type == EXPR_BIGINT)
                return mk_int(0);
#ifdef USE_MPFR
            if (f->type == EXPR_MPFR) return mk_int(0);
#endif
            if (f->type == EXPR_SYMBOL) {
                if (is_dt_constant_symbol(f)) return mk_int(0);
                /* f is some other symbol (and not x, handled above). */
                return mk_fn2("Dt", expr_copy(f), expr_copy(x));
            }
            /* Compound: skip the free-of-x short-circuit, fall through. */
        } else {
        /* Don't short-circuit `Equal[...]` or `Inequality[...]` even when
         * they have no x: the dispatch below distributes D over each
         * value slot and lets the outer evaluator simplify the residue,
         * matching Mathematica's semantics for D[a == b, x] and
         * D[a < b < c, x]. */
        bool f_is_equal_or_ineq = (f->type == EXPR_FUNCTION &&
                                   f->data.function.head->type == EXPR_SYMBOL &&
                                   (f->data.function.head->data.symbol.name == SYM_Equal
                                    || f->data.function.head->data.symbol.name == SYM_Inequality));
        if (!f_is_equal_or_ineq &&
            expr_free_of(f, x) && !expr_contains_nonconst(f, nonconsts)) {
            return mk_int(0);
        }
        }
    } else {
        /* Dt mode: numeric atoms and distinguished constants vanish. */
        if (f->type == EXPR_INTEGER || f->type == EXPR_REAL || f->type == EXPR_BIGINT) {
            return mk_int(0);
        }
        if (is_dt_constant_symbol(f)) return mk_int(0);
        if (f->type == EXPR_SYMBOL) {
            /* Unknown symbols are treated as implicit functions of an
             * implicit parameter -- their total derivative is Dt[x]. */
            return mk_fn1("Dt", expr_copy(f));
        }
    }

    /* At this point f must be compound. */
    if (f->type != EXPR_FUNCTION) {
        /* Strings and other atoms we cannot differentiate -- bail. */
        return NULL;
    }

    Expr* head = f->data.function.head;
    size_t n = f->data.function.arg_count;
    Expr** args = f->data.function.args;

    /* ------------------------------------------------------------------ */
    /* Head is a plain symbol -- dispatch on the head name.               */
    /* ------------------------------------------------------------------ */
    if (head->type == EXPR_SYMBOL) {
        const char* h = head->data.symbol.name;

        /* --- Plus: sum of derivatives. --- */
        if (h == SYM_Plus) {
            Expr** ts = malloc(sizeof(Expr*) * n);
            for (size_t i = 0; i < n; i++) ts[i] = deriv_of(args[i], x, nonconsts);
            return mk_fnN_adopt("Plus", ts, n);
        }

        /* --- Equal (a == b [== c ...]): distribute the derivative
         *     through each side. Mathematica returns Equal[D[a,x],
         *     D[b,x], ...], collapsing arms whose derivative is 0
         *     via the outer evaluator. The same convention applies
         *     for partial, total, and NonConstants-tagged derivatives. */
        if (h == SYM_Equal) {
            Expr** ts = malloc(sizeof(Expr*) * n);
            for (size_t i = 0; i < n; i++) ts[i] = deriv_of(args[i], x, nonconsts);
            return mk_fnN_adopt("Equal", ts, n);
        }

        /* --- Inequality[v0, op0, v1, op1, ...]: differentiate each value
         *     slot, keep the operator symbols verbatim. The outer
         *     evaluator may then collapse the residual (e.g. if all
         *     derivatives are 0, it becomes Inequality[0, op, 0, op, 0]
         *     and reduces in builtin_inequality). */
        if (h == SYM_Inequality && (n & 1u) == 1) {
            Expr** ts = malloc(sizeof(Expr*) * n);
            for (size_t i = 0; i < n; i++) {
                if ((i & 1u) == 1) ts[i] = expr_copy(args[i]);    /* op symbol */
                else               ts[i] = deriv_of(args[i], x, nonconsts);
            }
            return mk_fnN_adopt("Inequality", ts, n);
        }

        /* --- Times: general product rule. For n factors this is
         *      sum_i (d/dx arg_i) * prod_{j!=i} arg_j.
         * We build each term directly; the evaluator handles
         * simplification of Times[0, ...] and trivial identities. --- */
        if (h == SYM_Times) {
            Expr** ts = malloc(sizeof(Expr*) * n);
            for (size_t i = 0; i < n; i++) {
                Expr** factors = malloc(sizeof(Expr*) * n);
                for (size_t j = 0; j < n; j++) {
                    factors[j] = (i == j) ? deriv_of(args[j], x, nonconsts)
                                          : expr_copy(args[j]);
                }
                ts[i] = mk_fnN_adopt("Times", factors, n);
            }
            return mk_fnN_adopt("Plus", ts, n);
        }

        /* --- Dot, NonCommutativeMultiply, Cross: multilinear, so the
         * product rule keeps the factor order:
         *     d/dx a.b = a'.b + a.b'
         * (the generic chain rule gave Derivative[1, 0][Dot][a, b] a', as
         * the reference never does; FeynCalc's FunctionalD/FeynRule
         * differentiate DOT products of fields).  A factor whose
         * derivative is 0 contributes no term (D[a . b[x], x] = a . b'[x]). */
        if ((h == SYM_Dot || h == SYM_NonCommutativeMultiply || h == SYM_Cross) && n >= 1) {
            Expr** ts = malloc(sizeof(Expr*) * n);
            size_t nt = 0;
            for (size_t i = 0; i < n; i++) {
                Expr* di = deriv_of(args[i], x, nonconsts);
                if (is_lit_zero(di)) { expr_free(di); continue; }
                Expr** fs = malloc(sizeof(Expr*) * n);
                for (size_t j = 0; j < n; j++)
                    fs[j] = (i == j) ? di : expr_copy(args[j]);
                ts[nt++] = mk_fnN_adopt(h, fs, n);
            }
            if (nt == 0) { free(ts); return mk_int(0); }
            if (nt == 1) { Expr* t = ts[0]; free(ts); return t; }
            return mk_fnN_adopt("Plus", ts, nt);
        }

        /* --- Power: d/dx f^g = f^(g-1) * (g*f' + f*Log[f]*g'). ---
         *
         * When g turns out to be independent of x we collapse to the
         * usual g * f^(g-1) * f'. This avoids introducing spurious
         * Log[f] terms that the outer evaluator would then have to
         * fight to remove (and which are problematic for non-positive
         * bases). We detect "g is constant" by computing d g / d x once
         * and checking for the literal 0. */
        if (h == SYM_Power && n == 2) {
            Expr* a = args[0];
            Expr* b = args[1];
            Expr* da = deriv_of(a, x, nonconsts);
            Expr* db = deriv_of(b, x, nonconsts);

            if (is_lit_zero(db)) {
                expr_free(db);
                Expr* bm1 = mk_fn2("Plus", expr_copy(b), mk_int(-1));
                Expr* f_pow = mk_fn2("Power", expr_copy(a), bm1);
                Expr** factors = malloc(sizeof(Expr*) * 3);
                factors[0] = expr_copy(b);
                factors[1] = f_pow;
                factors[2] = da;
                return mk_fnN_adopt("Times", factors, 3);
            }

            Expr* bm1 = mk_fn2("Plus", expr_copy(b), mk_int(-1));
            Expr* f_pow = mk_fn2("Power", expr_copy(a), bm1);
            Expr* t1 = mk_fn2("Times", expr_copy(b), da);
            Expr* t2_inner = mk_fn2("Times", expr_copy(a), mk_fn1("Log", expr_copy(a)));
            Expr* t2 = mk_fn2("Times", t2_inner, db);
            Expr* inner = mk_fn2("Plus", t1, t2);
            return mk_fn2("Times", f_pow, inner);
        }

        /* --- List: thread element-wise. --- */
        if (h == SYM_List) {
            Expr** ts = malloc(sizeof(Expr*) * n);
            for (size_t i = 0; i < n; i++) ts[i] = deriv_of(args[i], x, nonconsts);
            return mk_fnN_adopt("List", ts, n);
        }

        /* --- RootSum[Function[t, p], Function[t, body]]: thread the
         *     derivative through the body Function.  The bound variable
         *     `t` ranges over roots of p[t] == 0, which are independent
         *     of the outer differentiation variable x; only the body
         *     can depend on x.  See src/root.c for the construct's
         *     definition. */
        if (h == SYM_RootSum && n == 2) {
            Expr* fn1 = args[0];
            Expr* fn2 = args[1];
            if (fn2->type == EXPR_FUNCTION
                && fn2->data.function.head->type == EXPR_SYMBOL
                && fn2->data.function.head->data.symbol.name == SYM_Function) {
                /* Two accepted body forms (the bound variable is
                 * always independent of x, so D threads through):
                 *   2-arg  Function[t, body]      — named bound var
                 *   1-arg  Function[body]         — Slot[_]-bound
                 *                                   ((... &) syntax)
                 * Both are produced by root_make_rootsum at different
                 * historical points; we accept either. */
                size_t fa = fn2->data.function.arg_count;
                if (fa == 2) {
                    Expr* bvar = fn2->data.function.args[0];
                    Expr* body = fn2->data.function.args[1];
                    Expr* dbody = deriv_of(body, x, nonconsts);
                    Expr* dbody_eval = eval_and_free(dbody);
                    Expr* new_fn2 = mk_fn2("Function",
                        expr_copy(bvar), dbody_eval);
                    return mk_fn2("RootSum", expr_copy(fn1), new_fn2);
                }
                if (fa == 1) {
                    Expr* body = fn2->data.function.args[0];
                    Expr* dbody = deriv_of(body, x, nonconsts);
                    Expr* dbody_eval = eval_and_free(dbody);
                    Expr* new_fn2 = mk_fn1("Function", dbody_eval);
                    return mk_fn2("RootSum", expr_copy(fn1), new_fn2);
                }
            }
            /* Unrecognised body form: leave RootSum unevaluated. */
            return NULL;
        }

        /* --- Differentiation through an (unevaluated) Sum / Product:
         *   D[Sum[f, s1, ..., sm], x]     = Sum[D[f, x], s1, ..., sm]
         *   D[Product[f, s1, ..., sm], x] = Product[f, s..] Sum[D[f, x]/f, s..]
         * valid when every iterator spec is {k, ...} with a symbol index
         * k != x (and not a NonConstant) and the rest of the spec (limits,
         * step, or explicit value list) free of x.  Without this rule the
         * generic chain rule treated Sum as an opaque two-argument function
         * and emitted Derivative[1, 0][Sum][f, {k, 1, n}] D[f, x] -- with the
         * bound index k leaked into the outer factor, a meaningless result. */
        if (x && !g_dt_total_mode && (h == SYM_Sum || h == SYM_Product) &&
            n >= 2) {
            bool ok = true;
            for (size_t s = 1; s < n && ok; s++) {
                const Expr* sp = args[s];
                const Expr* idx = sp;          /* indefinite form: a bare index */
                if (sp->type == EXPR_FUNCTION) {
                    if (sp->data.function.head->type != EXPR_SYMBOL ||
                        sp->data.function.head->data.symbol.name != SYM_List ||
                        sp->data.function.arg_count < 1 ||
                        sp->data.function.arg_count > 4) {
                        ok = false;
                        break;
                    }
                    idx = sp->data.function.args[0];
                    for (size_t j = 1; j < sp->data.function.arg_count; j++) {
                        const Expr* lim = sp->data.function.args[j];
                        if (!expr_free_of(lim, x) ||
                            expr_contains_nonconst(lim, nonconsts)) {
                            ok = false;
                            break;
                        }
                    }
                }
                if (idx->type != EXPR_SYMBOL || expr_eq(idx, x) ||
                    nonconsts_contains(nonconsts, idx))
                    ok = false;
            }
            /* A limit (or the index itself) depending on x: no term-wise
             * rule applies, and the generic chain rule would leak the bound
             * index -- leave D[Sum[...], x] unevaluated instead. */
            if (!ok) return NULL;
            /* Sum and Product are HoldAll: the term-wise derivative must be
             * evaluated here, or it would stay as a raw, unsimplified tree
             * (0 x^k + k x^(k - 1) 1 ...) inside the held sum. */
            Expr* dterm = eval_and_free(deriv_of(args[0], x, nonconsts));
            Expr** items = malloc(n * sizeof(Expr*));
            if (h == SYM_Sum) {
                items[0] = dterm;
            } else {
                items[0] = eval_and_free(mk_fn2("Times", dterm,
                                  mk_fn2("Power", expr_copy(args[0]), mk_int(-1))));
            }
            for (size_t s = 1; s < n; s++) items[s] = expr_copy(args[s]);
            Expr* dsum = mk_fnN_adopt(SYM_Sum, items, n);
            if (h == SYM_Sum) return dsum;
            return mk_fn2("Times", expr_copy(f), dsum);
        }

        /* --- Fundamental theorem of calculus: D[Integrate[f, x], x] = f. ---
         * Only the indefinite form with integration variable identical to the
         * differentiation variable is the FTC.  A definite/iterated integral
         * (Integrate[f, {x, a, b}], 2nd arg a List) or a different integration
         * variable is left to fall through (differentiation under the integral
         * sign is not done here). */
        if (h == SYM_Integrate && n == 2 && expr_eq(args[1], x))
            return expr_copy(args[0]);

        /* --- Leibniz rule for a variable-limit definite integral:
         *   D[Integrate[e, {u, a, b}], x]
         *     = (e/.u->b) D[b,x] - (e/.u->a) D[a,x] + Integrate[D[e,x], {u,a,b}].
         * The bound variable u must be a symbol distinct from x.  This replaces
         * the generic fallback, which would emit an inert Derivative[_,_][
         * Integrate] with the bound variable leaked into the result. */
        if (h == SYM_Integrate && n == 2 &&
            args[1]->type == EXPR_FUNCTION &&
            args[1]->data.function.head->type == EXPR_SYMBOL &&
            args[1]->data.function.head->data.symbol.name == SYM_List &&
            args[1]->data.function.arg_count == 3 &&
            args[1]->data.function.args[0]->type == EXPR_SYMBOL &&
            !expr_eq(args[1]->data.function.args[0], x)) {
            Expr* e  = args[0];
            Expr* u  = args[1]->data.function.args[0];
            Expr* lo = args[1]->data.function.args[1];
            Expr* hi = args[1]->data.function.args[2];
            Expr* e_hi = eval_and_free(mk_fn2("ReplaceAll", expr_copy(e),
                             mk_fn2("Rule", expr_copy(u), expr_copy(hi))));
            Expr* e_lo = eval_and_free(mk_fn2("ReplaceAll", expr_copy(e),
                             mk_fn2("Rule", expr_copy(u), expr_copy(lo))));
            Expr* dhi = eval_and_free(deriv_of(hi, x, nonconsts));
            Expr* dlo = eval_and_free(deriv_of(lo, x, nonconsts));
            Expr* de  = eval_and_free(deriv_of(e, x, nonconsts));
            Expr** listitems = malloc(3 * sizeof(Expr*));
            listitems[0] = expr_copy(u);
            listitems[1] = expr_copy(lo);
            listitems[2] = expr_copy(hi);
            Expr* newspec = mk_fnN_adopt(SYM_List, listitems, 3);
            Expr** terms = malloc(3 * sizeof(Expr*));
            terms[0] = mk_fn2("Times", e_hi, dhi);
            terms[1] = mk_neg(mk_fn2("Times", e_lo, dlo));
            terms[2] = mk_fn2("Integrate", de, newspec);
            return mk_fnN_adopt(SYM_Plus, terms, 3);
        }

        /* --- Log[b, f]: reduce to Log[f]/Log[b]. --- */
        if (h == SYM_Log && n == 2) {
            Expr* lb = mk_fn1("Log", expr_copy(args[0]));
            Expr* lf = mk_fn1("Log", expr_copy(args[1]));
            Expr* quot = mk_fn2("Times", lf, mk_fn2("Power", lb, mk_int(-1)));
            Expr* r = compute_deriv(quot, x, nonconsts);
            expr_free(quot);
            return r;
        }

        /* --- ArcTan[u, v] (two-argument, = Arg[u + I v] = atan2(v, u)):
         *   D = (u D[v] - v D[u]) / (u^2 + v^2). --- */
        if (h == SYM_ArcTan && n == 2) {
            Expr* u = args[0];
            Expr* v = args[1];
            Expr* du = deriv_of(u, x, nonconsts);
            Expr* dv = deriv_of(v, x, nonconsts);
            Expr* num = mk_fn2("Plus",
                mk_fn2("Times", expr_copy(u), dv),
                mk_neg(mk_fn2("Times", expr_copy(v), du)));
            Expr* den = mk_fn2("Plus",
                mk_fn2("Power", expr_copy(u), mk_int(2)),
                mk_fn2("Power", expr_copy(v), mk_int(2)));
            return mk_fn2("Times", num, mk_fn2("Power", den, mk_int(-1)));
        }

        /* --- MeijerG[{{a..},{a'..}}, {{b..},{b'..}}, z(, r)]: derivative w.r.t.
         * z (parameters free of x), times D[z, x]:
         *   3 args: MeijerG[{{-1, a-1..},{a'-1..}}, {{b-1..},{0, b'-1..}}, z]
         *           (Mathematica's form: z d/dz G = G^{m,n+1}_{p+1,q+1} with the
         *           extra 0 / 1 pair, then z^-1 G(a; b) = G(a-1; b-1));
         *   4 args: (1/(r z)) MeijerG[{{0, a..},{a'..}}, {{b..},{1, b'..}}, z, r],
         *           since z d/dz z^-s = -s z^-s and -s = (1/r)(-r s) is the ratio
         *           Gamma(1-0-rs)/Gamma(1-1-rs). */
        if (h == SYM_MeijerG && (n == 3 || n == 4)
            && expr_free_of(args[0], x) && expr_free_of(args[1], x)
            && (n == 3 || expr_free_of(args[3], x))) {
            const Expr* g[4];
            bool shape = true;
            for (int s2 = 0; s2 < 2 && shape; s2++) {
                const Expr* L = args[s2];
                if (L->type != EXPR_FUNCTION || L->data.function.head->type != EXPR_SYMBOL
                    || L->data.function.head->data.symbol.name != SYM_List
                    || L->data.function.arg_count != 2) { shape = false; break; }
                for (int t = 0; t < 2; t++) {
                    const Expr* G = L->data.function.args[t];
                    if (G->type != EXPR_FUNCTION || G->data.function.head->type != EXPR_SYMBOL
                        || G->data.function.head->data.symbol.name != SYM_List) { shape = false; break; }
                    g[2 * s2 + t] = G;
                }
            }
            if (shape) {
                Expr* z  = args[2];
                Expr* dz = deriv_of(z, x, nonconsts);
                if (is_lit_zero(dz)) { expr_free(dz); return mk_int(0); }
                long shift = (n == 3) ? -1 : 0;
                /* group k gets an extra leading (k == 0) or trailing-first (k == 3)
                 * element: -1 / 0 (3 args) or 0 / 1 (4 args). */
                Expr* grp[4];
                for (int k = 0; k < 4; k++) {
                    size_t len = g[k]->data.function.arg_count;
                    bool extra = (k == 0 || k == 3);
                    Expr** it = malloc(sizeof(Expr*) * (len + 1));
                    size_t w = 0;
                    if (extra) it[w++] = mk_int(k == 0 ? (n == 3 ? -1 : 0) : (n == 3 ? 0 : 1));
                    for (size_t i = 0; i < len; i++) {
                        Expr* e = expr_copy(g[k]->data.function.args[i]);
                        it[w++] = shift ? mk_fn2("Plus", e, mk_int(shift)) : e;
                    }
                    grp[k] = expr_new_function(mk_sym("List"), it, w);
                    free(it);
                }
                Expr* al = mk_fn2("List", grp[0], grp[1]);
                Expr* bl = mk_fn2("List", grp[2], grp[3]);
                Expr* G;
                if (n == 3) {
                    G = expr_new_function(mk_sym("MeijerG"),
                            (Expr*[]){ al, bl, expr_copy(z) }, 3);
                } else {
                    G = expr_new_function(mk_sym("MeijerG"),
                            (Expr*[]){ al, bl, expr_copy(z), expr_copy(args[3]) }, 4);
                    G = mk_fn2("Times", G,
                            mk_fn2("Power", mk_fn2("Times", expr_copy(args[3]), expr_copy(z)),
                                   mk_int(-1)));
                }
                return mk_fn2("Times", G, dz);
            }
        }

        /* --- HypergeometricPFQ[{a..}, {b..}, z]: derivative w.r.t. z is
         *   (prod a_i / prod b_j) pFq[{a_i+1}, {b_j+1}, z], times D[z, x].
         * Only the z-derivative is handled here; if the parameter lists
         * depend on x, defer to the generic chain rule (return NULL). */
        if (h == SYM_HypergeometricPFQ && n == 3
            && args[0]->type == EXPR_FUNCTION && args[1]->type == EXPR_FUNCTION
            && expr_free_of(args[0], x) && expr_free_of(args[1], x)) {
            Expr* al = args[0];
            Expr* bl = args[1];
            Expr* z  = args[2];
            Expr* dz = deriv_of(z, x, nonconsts);
            if (is_lit_zero(dz)) { expr_free(dz); return mk_int(0); }
            size_t p = al->data.function.arg_count;
            size_t q = bl->data.function.arg_count;

            /* prefactor = prod a_i / prod b_j */
            Expr* pref;
            {
                size_t w = 0;
                Expr** pf = malloc(sizeof(Expr*) * (p + 1));
                for (size_t i = 0; i < p; i++) pf[w++] = expr_copy(al->data.function.args[i]);
                if (q > 0) {
                    Expr** db = malloc(sizeof(Expr*) * q);
                    for (size_t j = 0; j < q; j++) db[j] = expr_copy(bl->data.function.args[j]);
                    Expr* bprod = mk_fnN_adopt("Times", db, q);
                    pf[w++] = mk_fn2("Power", bprod, mk_int(-1));
                }
                if (w == 0) { free(pf); pref = mk_int(1); }
                else pref = mk_fnN_adopt("Times", pf, w);
            }

            /* shifted parameter lists {a_i + 1}, {b_j + 1} */
            Expr** na = p ? malloc(sizeof(Expr*) * p) : NULL;
            for (size_t i = 0; i < p; i++)
                na[i] = mk_fn2("Plus", expr_copy(al->data.function.args[i]), mk_int(1));
            Expr* nal = expr_new_function(mk_sym("List"), na, p);
            free(na);
            Expr** nb = q ? malloc(sizeof(Expr*) * q) : NULL;
            for (size_t j = 0; j < q; j++)
                nb[j] = mk_fn2("Plus", expr_copy(bl->data.function.args[j]), mk_int(1));
            Expr* nbl = expr_new_function(mk_sym("List"), nb, q);
            free(nb);

            Expr* newp = expr_new_function(mk_sym("HypergeometricPFQ"),
                             (Expr*[]){ nal, nbl, expr_copy(z) }, 3);
            Expr** ff = malloc(sizeof(Expr*) * 3);
            ff[0] = pref; ff[1] = newp; ff[2] = dz;
            return mk_fnN_adopt("Times", ff, 3);
        }

        /* --- Fibonacci[A, B]: chain rule through both arguments.
         *   dF/dB = (2 A F[A-1,B] + (A-1) B F[A,B]) / (4 + B^2)
         *   dF/dA = (Pi + 2 ArcSinh[B/2] Cot[A Pi]) Csc[A Pi] F[-A,B]
         *           + (Pi Cot[A Pi]
         *              + ArcSinh[B/2] (Cot[A Pi]^2 + Csc[A Pi]^2)) F[A,B]
         * The total derivative is dF/dA D[A,x] + dF/dB D[B,x]; zero-derivative
         * arms are dropped (so D[F[n,x],x] keeps only the dF/dB term, and
         * D[F[n,x],n] keeps only dF/dA). */
        if (h == SYM_Fibonacci && n == 2) {
            Expr* A = args[0];
            Expr* B = args[1];
            Expr* dA = deriv_of(A, x, nonconsts);
            Expr* dB = deriv_of(B, x, nonconsts);
            Expr* terms[2];
            size_t nt = 0;

            if (!is_lit_zero(dB)) {
                Expr* fAm1 = mk_fn2("Fibonacci",
                                    mk_fn2("Plus", expr_copy(A), mk_int(-1)),
                                    expr_copy(B));
                Expr* fA = mk_fn2("Fibonacci", expr_copy(A), expr_copy(B));
                Expr* t1 = mk_fn2("Times",
                                  mk_fn2("Times", mk_int(2), expr_copy(A)), fAm1);
                Expr* t2 = mk_fn2("Times",
                                  mk_fn2("Times",
                                         mk_fn2("Plus", expr_copy(A), mk_int(-1)),
                                         expr_copy(B)),
                                  fA);
                Expr* numer = mk_fn2("Plus", t1, t2);
                Expr* denom = mk_fn2("Plus", mk_int(4),
                                     mk_fn2("Power", expr_copy(B), mk_int(2)));
                Expr* dFdB = mk_fn2("Times", numer,
                                    mk_fn2("Power", denom, mk_int(-1)));
                terms[nt++] = mk_fn2("Times", dFdB, dB);
            } else {
                expr_free(dB);
            }

            if (!is_lit_zero(dA)) {
                Expr* asinh = mk_fn1("ArcSinh",
                                     mk_fn2("Times", expr_copy(B),
                                            mk_fn2("Power", mk_int(2), mk_int(-1))));
                Expr* cot = mk_fn1("Cot",
                                   mk_fn2("Times", expr_copy(A), mk_sym("Pi")));
                Expr* csc = mk_fn1("Csc",
                                   mk_fn2("Times", expr_copy(A), mk_sym("Pi")));
                Expr* fnegA = mk_fn2("Fibonacci",
                                     mk_fn2("Times", mk_int(-1), expr_copy(A)),
                                     expr_copy(B));
                Expr* fA = mk_fn2("Fibonacci", expr_copy(A), expr_copy(B));

                /* coef1 = (Pi + 2 ArcSinh[B/2] Cot[A Pi]) Csc[A Pi] */
                Expr* c1inner = mk_fn2("Plus", mk_sym("Pi"),
                                       mk_fn2("Times",
                                              mk_fn2("Times", mk_int(2), expr_copy(asinh)),
                                              expr_copy(cot)));
                Expr* coef1 = mk_fn2("Times", c1inner, expr_copy(csc));

                /* coef2 = Pi Cot[A Pi] + ArcSinh[B/2] (Cot[A Pi]^2 + Csc[A Pi]^2) */
                Expr* cot2 = mk_fn2("Power", expr_copy(cot), mk_int(2));
                Expr* csc2 = mk_fn2("Power", csc, mk_int(2));
                Expr* c2b = mk_fn2("Times", asinh, mk_fn2("Plus", cot2, csc2));
                Expr* c2a = mk_fn2("Times", mk_sym("Pi"), cot);
                Expr* coef2 = mk_fn2("Plus", c2a, c2b);

                Expr* dFdA = mk_fn2("Plus",
                                    mk_fn2("Times", coef1, fnegA),
                                    mk_fn2("Times", coef2, fA));
                terms[nt++] = mk_fn2("Times", dFdA, dA);
            } else {
                expr_free(dA);
            }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return mk_fn2("Plus", terms[0], terms[1]);
        }

        /* --- LucasL[A, B]: chain rule through both arguments.
         *   dL/dB = (A (2 L[A-1,B] + B L[A,B])) / (4 + B^2)
         *   dL/dA = L[A+1,B] (2 ArcSinh[B/2] + Pi Tan[A Pi]) / Sqrt[4+B^2]
         *           + L[A,B] (-2 B ArcSinh[B/2]
         *                     - Pi (B + Sqrt[4+B^2]) Tan[A Pi]) / (2 Sqrt[4+B^2])
         * The total derivative is dL/dA D[A,x] + dL/dB D[B,x]; zero-derivative
         * arms are dropped (so D[L[n,x],x] keeps only the dL/dB term, and
         * D[L[n,x],n] keeps only dL/dA). */
        if (h == SYM_LucasL && n == 2) {
            Expr* A = args[0];
            Expr* B = args[1];
            Expr* dA = deriv_of(A, x, nonconsts);
            Expr* dB = deriv_of(B, x, nonconsts);
            Expr* terms[2];
            size_t nt = 0;

            if (!is_lit_zero(dB)) {
                Expr* lAm1 = mk_fn2("LucasL",
                                    mk_fn2("Plus", expr_copy(A), mk_int(-1)),
                                    expr_copy(B));
                Expr* lA = mk_fn2("LucasL", expr_copy(A), expr_copy(B));
                Expr* inner = mk_fn2("Plus",
                                     mk_fn2("Times", mk_int(2), lAm1),
                                     mk_fn2("Times", expr_copy(B), lA));
                Expr* numer = mk_fn2("Times", expr_copy(A), inner);
                Expr* denom = mk_fn2("Plus", mk_int(4),
                                     mk_fn2("Power", expr_copy(B), mk_int(2)));
                Expr* dLdB = mk_fn2("Times", numer,
                                    mk_fn2("Power", denom, mk_int(-1)));
                terms[nt++] = mk_fn2("Times", dLdB, dB);
            } else {
                expr_free(dB);
            }

            if (!is_lit_zero(dA)) {
                Expr* S = mk_fn1("Sqrt",
                                 mk_fn2("Plus", mk_int(4),
                                        mk_fn2("Power", expr_copy(B), mk_int(2))));
                Expr* asinh = mk_fn1("ArcSinh",
                                     mk_fn2("Times", expr_copy(B),
                                            mk_fn2("Power", mk_int(2), mk_int(-1))));
                Expr* tan = mk_fn1("Tan",
                                   mk_fn2("Times", expr_copy(A), mk_sym("Pi")));
                Expr* lAp1 = mk_fn2("LucasL",
                                    mk_fn2("Plus", expr_copy(A), mk_int(1)),
                                    expr_copy(B));
                Expr* lA = mk_fn2("LucasL", expr_copy(A), expr_copy(B));

                /* term1 = L[A+1,B] (2 ArcSinh[B/2] + Pi Tan[A Pi]) / Sqrt[4+B^2] */
                Expr* c1 = mk_fn2("Plus",
                                  mk_fn2("Times", mk_int(2), expr_copy(asinh)),
                                  mk_fn2("Times", mk_sym("Pi"), expr_copy(tan)));
                Expr* term1 = mk_fn2("Times",
                                     mk_fn2("Times", lAp1, c1),
                                     mk_fn2("Power", expr_copy(S), mk_int(-1)));

                /* term2 = L[A,B] (-2 B ArcSinh[B/2] - Pi (B + S) Tan[A Pi])
                 *         / (2 Sqrt[4+B^2]) */
                Expr* c2a = mk_neg(mk_fn2("Times",
                                   mk_fn2("Times", mk_int(2), expr_copy(B)),
                                   asinh));                  /* consumes asinh */
                Expr* c2b = mk_neg(mk_fn2("Times",
                                   mk_fn2("Times", mk_sym("Pi"),
                                          mk_fn2("Plus", expr_copy(B), S)),  /* consumes S */
                                   tan));                    /* consumes tan */
                Expr* c2 = mk_fn2("Plus", c2a, c2b);
                Expr* denom2 = mk_fn2("Times", mk_int(2),
                                      mk_fn1("Sqrt",
                                             mk_fn2("Plus", mk_int(4),
                                                    mk_fn2("Power", expr_copy(B),
                                                           mk_int(2)))));
                Expr* term2 = mk_fn2("Times",
                                     mk_fn2("Times", lA, c2),
                                     mk_fn2("Power", denom2, mk_int(-1)));

                Expr* dLdA = mk_fn2("Plus", term1, term2);
                terms[nt++] = mk_fn2("Times", dLdA, dA);
            } else {
                expr_free(dA);
            }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return mk_fn2("Plus", terms[0], terms[1]);
        }

        /* --- Zeta[S, A] (Hurwitz zeta): chain rule on both args.
         *   d/dA Zeta[S, A] = -S Zeta[1+S, A]   (elementary)
         *   d/dS Zeta[S, A] = Derivative[1,0][Zeta][S,A]  (generic; the
         *               derivative with respect to the order s has no
         *               elementary closed form).
         * Zero-derivative arms are dropped, so D[Zeta[s,a],a] keeps only the
         * dA term, matching Mathematica's -s Zeta[1+s, a]. */
        if (h == SYM_Zeta && n == 2) {
            Expr* S = args[0];
            Expr* A = args[1];
            Expr* dS = deriv_of(S, x, nonconsts);
            Expr* dA = deriv_of(A, x, nonconsts);
            Expr* terms[2];
            size_t nt = 0;

            if (!is_lit_zero(dA)) {
                /* -S Zeta[1+S, A] * dA */
                Expr* s1 = mk_fn2("Plus", mk_int(1), expr_copy(S));
                Expr* z  = mk_fn2("Zeta", s1, expr_copy(A));
                Expr* dZdA = mk_neg(mk_fn2("Times", expr_copy(S), z));
                terms[nt++] = mk_fn2("Times", dZdA, dA);
            } else {
                expr_free(dA);
            }

            if (!is_lit_zero(dS)) {
                /* Derivative[1, 0][Zeta][S, A] * dS */
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(1), mk_int(0) }, 2);
                Expr* op_g = mk_fn_head1(op, mk_sym("Zeta"));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(S), expr_copy(A) }, 2);
                terms[nt++] = mk_fn2("Times", applied, dS);
            } else {
                expr_free(dS);
            }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return mk_fn2("Plus", terms[0], terms[1]);
        }

        /* --- HurwitzZeta[S, A]: chain rule on both args.
         *   d/dA HurwitzZeta[S, A] = -S HurwitzZeta[1+S, A]   (elementary)
         *   d/dS HurwitzZeta[S, A] = Derivative[1,0][HurwitzZeta][S,A] (generic;
         *               the derivative in the order s has no elementary form).
         * Zero-derivative arms are dropped, so D[HurwitzZeta[s,a],a] keeps only
         * the dA term, matching Mathematica's -s HurwitzZeta[1+s, a]. */
        if (h == SYM_HurwitzZeta && n == 2) {
            Expr* S = args[0];
            Expr* A = args[1];
            Expr* dS = deriv_of(S, x, nonconsts);
            Expr* dA = deriv_of(A, x, nonconsts);
            Expr* terms[2];
            size_t nt = 0;

            if (!is_lit_zero(dA)) {
                /* -S HurwitzZeta[1+S, A] * dA */
                Expr* s1 = mk_fn2("Plus", mk_int(1), expr_copy(S));
                Expr* z  = mk_fn2("HurwitzZeta", s1, expr_copy(A));
                Expr* dZdA = mk_neg(mk_fn2("Times", expr_copy(S), z));
                terms[nt++] = mk_fn2("Times", dZdA, dA);
            } else {
                expr_free(dA);
            }

            if (!is_lit_zero(dS)) {
                /* Derivative[1, 0][HurwitzZeta][S, A] * dS */
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(1), mk_int(0) }, 2);
                Expr* op_g = mk_fn_head1(op, mk_sym("HurwitzZeta"));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(S), expr_copy(A) }, 2);
                terms[nt++] = mk_fn2("Times", applied, dS);
            } else {
                expr_free(dS);
            }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return mk_fn2("Plus", terms[0], terms[1]);
        }

        /* --- HarmonicNumber[N] and HarmonicNumber[N, R]: derivative w.r.t. the
         *   continuous index N.
         *   d/dN H_N     = Zeta[2] - HarmonicNumber[N, 2]
         *   d/dN H_N^(R) = R (Zeta[R+1] - HarmonicNumber[N, R+1])
         *   d/dR H_N^(R) = Derivative[0,1][HarmonicNumber][N, R] (no elementary
         *                  closed form in the order R). */
        if (h == SYM_HarmonicNumber && n == 1) {
            Expr* N = args[0];
            Expr* dN = deriv_of(N, x, nonconsts);
            if (is_lit_zero(dN)) { expr_free(dN); return mk_int(0); }
            Expr* dHdN = mk_fn2("Plus", mk_fn1("Zeta", mk_int(2)),
                                mk_neg(mk_fn2("HarmonicNumber",
                                              expr_copy(N), mk_int(2))));
            return mk_fn2("Times", dHdN, dN);
        }
        if (h == SYM_HarmonicNumber && n == 2) {
            Expr* N = args[0];
            Expr* R = args[1];
            Expr* dN = deriv_of(N, x, nonconsts);
            Expr* dR = deriv_of(R, x, nonconsts);
            Expr* terms[2];
            size_t nt = 0;

            if (!is_lit_zero(dN)) {
                /* R (Zeta[R+1] - HarmonicNumber[N, R+1]) * dN */
                Expr* zr = mk_fn1("Zeta", mk_fn2("Plus", expr_copy(R), mk_int(1)));
                Expr* hn = mk_fn2("HarmonicNumber", expr_copy(N),
                                  mk_fn2("Plus", expr_copy(R), mk_int(1)));
                Expr* inner = mk_fn2("Plus", zr, mk_neg(hn));
                Expr* dHdN = mk_fn2("Times", expr_copy(R), inner);
                terms[nt++] = mk_fn2("Times", dHdN, dN);
            } else {
                expr_free(dN);
            }

            if (!is_lit_zero(dR)) {
                /* Derivative[0, 1][HarmonicNumber][N, R] * dR */
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(0), mk_int(1) }, 2);
                Expr* op_g = mk_fn_head1(op, mk_sym("HarmonicNumber"));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(N), expr_copy(R) }, 2);
                terms[nt++] = mk_fn2("Times", applied, dR);
            } else {
                expr_free(dR);
            }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return mk_fn2("Plus", terms[0], terms[1]);
        }

        /* --- Zeta[S] (one-argument Riemann zeta): d/dS Zeta[S] has no
         * elementary closed form, so emit Derivative[1][Zeta][S] * D[S, x]. */
        if (h == SYM_Zeta && n == 1) {
            Expr* S = args[0];
            Expr* dS = deriv_of(S, x, nonconsts);
            if (is_lit_zero(dS)) { expr_free(dS); return mk_int(0); }
            Expr* op = expr_new_function(mk_sym("Derivative"),
                          (Expr*[]){ mk_int(1) }, 1);
            Expr* op_g = mk_fn_head1(op, mk_sym("Zeta"));
            Expr* applied = expr_new_function(op_g, (Expr*[]){ expr_copy(S) }, 1);
            return mk_fn2("Times", applied, dS);
        }

        /* --- PolyLog[N, Z] (polylogarithm): chain rule on both args.
         *   d/dZ PolyLog[N, Z] = PolyLog[N-1, Z]/Z   (elementary)
         *   d/dN PolyLog[N, Z] = Derivative[1,0][PolyLog][N, Z]  (no elementary
         *               closed form in the order N).
         * Zero-derivative arms are dropped, so D[PolyLog[n,z],z] keeps only the
         * dZ term, matching Mathematica's PolyLog[n-1, z]/z. */
        if (h == SYM_PolyLog && n == 2) {
            Expr* N = args[0];
            Expr* Z = args[1];
            Expr* dN = deriv_of(N, x, nonconsts);
            Expr* dZ = deriv_of(Z, x, nonconsts);
            Expr* terms[2];
            size_t nt = 0;

            if (!is_lit_zero(dZ)) {
                /* PolyLog[N-1, Z] Z^-1 * dZ */
                Expr* nm1 = mk_fn2("Plus", expr_copy(N), mk_int(-1));
                Expr* pl  = mk_fn2("PolyLog", nm1, expr_copy(Z));
                Expr* dPdZ = mk_fn2("Times", pl,
                                    mk_fn2("Power", expr_copy(Z), mk_int(-1)));
                terms[nt++] = mk_fn2("Times", dPdZ, dZ);
            } else {
                expr_free(dZ);
            }

            if (!is_lit_zero(dN)) {
                /* Derivative[1, 0][PolyLog][N, Z] * dN */
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(1), mk_int(0) }, 2);
                Expr* op_g = mk_fn_head1(op, mk_sym("PolyLog"));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(N), expr_copy(Z) }, 2);
                terms[nt++] = mk_fn2("Times", applied, dN);
            } else {
                expr_free(dN);
            }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return mk_fn2("Plus", terms[0], terms[1]);
        }

        /* --- PolyLog[N, P, Z] (Nielsen generalized polylogarithm): chain rule.
         *   d/dZ S_{N,P}(Z) = S_{N-1,P}(Z)/Z  (from the integral representation;
         *               the recurrence that defines S_{N,P} in N)
         *   d/dN, d/dP = Derivative[1,0,0] / Derivative[0,1,0][PolyLog][N,P,Z]. */
        if (h == SYM_PolyLog && n == 3) {
            Expr* N = args[0];
            Expr* P = args[1];
            Expr* Z = args[2];
            Expr* dN = deriv_of(N, x, nonconsts);
            Expr* dP = deriv_of(P, x, nonconsts);
            Expr* dZ = deriv_of(Z, x, nonconsts);
            Expr* terms[3];
            size_t nt = 0;

            if (!is_lit_zero(dZ)) {
                Expr* pl = expr_new_function(mk_sym("PolyLog"), (Expr*[]){
                    mk_fn2("Plus", mk_int(-1), expr_copy(N)), expr_copy(P),
                    expr_copy(Z) }, 3);
                Expr* dPdZ = mk_fn2("Times", pl,
                                    mk_fn2("Power", expr_copy(Z), mk_int(-1)));
                terms[nt++] = mk_fn2("Times", dPdZ, dZ);
            } else {
                expr_free(dZ);
            }

            Expr* dvs[2] = { dN, dP };
            for (int k = 0; k < 2; k++) {
                if (is_lit_zero(dvs[k])) { expr_free(dvs[k]); continue; }
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(k == 0), mk_int(k == 1), mk_int(0) }, 3);
                Expr* op_g = mk_fn_head1(op, mk_sym("PolyLog"));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(N), expr_copy(P), expr_copy(Z) }, 3);
                terms[nt++] = mk_fn2("Times", applied, dvs[k]);
            }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return expr_new_function(mk_sym("Plus"), terms, nt);
        }

        /* --- LerchPhi[Z, S, A] (Lerch transcendent): chain rule on all three args.
         *   d/dZ Phi = (LerchPhi[Z,-1+S,A] - A LerchPhi[Z,S,A]) / Z   (elementary)
         *   d/dA Phi = -S LerchPhi[Z,1+S,A]                           (elementary)
         *   d/dS Phi = Derivative[0,1,0][LerchPhi][Z,S,A]  (no elementary form in
         *               the order S).
         * Zero-derivative arms are dropped. */
        if (h == SYM_LerchPhi && n == 3) {
            Expr* Z = args[0];
            Expr* S = args[1];
            Expr* A = args[2];
            Expr* dZ = deriv_of(Z, x, nonconsts);
            Expr* dS = deriv_of(S, x, nonconsts);
            Expr* dA = deriv_of(A, x, nonconsts);
            Expr* terms[3];
            size_t nt = 0;

            if (!is_lit_zero(dZ)) {
                /* (LerchPhi[Z,-1+S,A] - A LerchPhi[Z,S,A]) Z^-1 * dZ */
                Expr* sm1 = mk_fn2("Plus", mk_int(-1), expr_copy(S));
                Expr* lp1 = expr_new_function(mk_sym("LerchPhi"),
                              (Expr*[]){ expr_copy(Z), sm1, expr_copy(A) }, 3);
                Expr* lp0 = expr_new_function(mk_sym("LerchPhi"),
                              (Expr*[]){ expr_copy(Z), expr_copy(S), expr_copy(A) }, 3);
                Expr* num = mk_fn2("Plus", lp1,
                                   mk_neg(mk_fn2("Times", expr_copy(A), lp0)));
                Expr* dPdZ = mk_fn2("Times", num,
                                    mk_fn2("Power", expr_copy(Z), mk_int(-1)));
                terms[nt++] = mk_fn2("Times", dPdZ, dZ);
            } else {
                expr_free(dZ);
            }

            if (!is_lit_zero(dA)) {
                /* -S LerchPhi[Z, 1+S, A] * dA */
                Expr* sp1 = mk_fn2("Plus", mk_int(1), expr_copy(S));
                Expr* lp = expr_new_function(mk_sym("LerchPhi"),
                              (Expr*[]){ expr_copy(Z), sp1, expr_copy(A) }, 3);
                Expr* dPdA = mk_neg(mk_fn2("Times", expr_copy(S), lp));
                terms[nt++] = mk_fn2("Times", dPdA, dA);
            } else {
                expr_free(dA);
            }

            if (!is_lit_zero(dS)) {
                /* Derivative[0, 1, 0][LerchPhi][Z, S, A] * dS */
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(0), mk_int(1), mk_int(0) }, 3);
                Expr* op_g = mk_fn_head1(op, mk_sym("LerchPhi"));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(Z), expr_copy(S), expr_copy(A) }, 3);
                terms[nt++] = mk_fn2("Times", applied, dS);
            } else {
                expr_free(dS);
            }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            if (nt == 2) return mk_fn2("Plus", terms[0], terms[1]);
            return expr_new_function(mk_sym("Plus"), terms, 3);
        }

        /* --- BesselJ[N, Z]: chain rule on both args.
         *   d/dZ BesselJ[N, Z] = (BesselJ[N-1, Z] - BesselJ[N+1, Z]) / 2  (DLMF 10.6.1)
         *   d/dN BesselJ[N, Z] = Derivative[1,0][BesselJ][N,Z]  (the derivative
         *               with respect to the order has no elementary form).
         * Zero-derivative arms are dropped, so D[BesselJ[n,x],x] keeps only the
         * dZ term, matching Mathematica's 1/2 (BesselJ[-1+n,x]-BesselJ[1+n,x]). */
        /* --- ProductLog[K, Z] (Lambert W, k-th branch): the branch index K is
         *   discrete, so only the Z-derivative contributes.
         *   d/dZ ProductLog[K, Z] = ProductLog[K, Z] / (Z (1 + ProductLog[K, Z])). */
        if (h == SYM_ProductLog && n == 2) {
            Expr* K = args[0];
            Expr* Z = args[1];
            Expr* dZ = deriv_of(Z, x, nonconsts);
            if (is_lit_zero(dZ)) { expr_free(dZ); return mk_int(0); }
            Expr* w = mk_fn2("ProductLog", expr_copy(K), expr_copy(Z));
            Expr* one_plus_w = mk_fn2("Plus", mk_int(1), w);
            Expr* denom = mk_fn2("Times", expr_copy(Z), one_plus_w);
            Expr* dWdZ = mk_fn2("Times",
                              mk_fn2("ProductLog", expr_copy(K), expr_copy(Z)),
                              mk_fn2("Power", denom, mk_int(-1)));
            return mk_fn2("Times", dWdZ, dZ);
        }

        /* --- The Legendre elliptic integrals. Second (F, E) or third (Pi)
         * argument is the PARAMETER m = k^2, not the modulus.
         *
         *   d/dphi EllipticF[phi, m]     = 1 / Sqrt[1 - m Sin[phi]^2]
         *   d/dphi EllipticE[phi, m]     =     Sqrt[1 - m Sin[phi]^2]
         *   d/dphi EllipticPi[n, phi, m] = 1 / ((1 - n Sin[phi]^2) Sqrt[1 - m Sin[phi]^2])
         *   d/dm   EllipticK[m]          = (EllipticE[m] - (1 - m) EllipticK[m]) / (2 m (1 - m))
         *   d/dm   EllipticE[m]          = (EllipticE[m] - EllipticK[m]) / (2 m)
         *   d/dm   EllipticE[phi, m]     = (EllipticE[phi, m] - EllipticF[phi, m]) / (2 m)
         *
         *   d/dm   EllipticF[phi, m]     = -E(phi,m)/(2(m-1)m) - F(phi,m)/(2m)
         *                                   + Sin[2 phi]/(4(m-1) Sqrt[1 - m Sin[phi]^2])
         *   d/dn   EllipticPi[n, m]      = (n E + (m-n) K + (n^2-m) Pi) / (2(m-n)(n-1)n)
         *   d/dm   EllipticPi[n, m]      = (E/(m-1) + Pi) / (2(n-m))
         *
         * The phi-derivatives are just the integrands, which is why a numeric
         * verify of an antiderivative built from these kernels closes.
         *
         * The three m-/n- forms above were inert Derivative[...] placeholders
         * until each was checked against a central difference at 30 digits
         * (agreement ~1e-30), on the principle that an inert derivative is
         * honest where a guessed one corrupts every caller silently.
         *
         * The n- and m-derivatives of the INCOMPLETE EllipticPi stay inert, and
         * not for want of trying: a least-squares fit over the natural candidate
         * basis (E, F, Pi over (n-m), (m-1), (n-1), and the Sin[2 phi]/Sqrt[...]
         * boundary term) does NOT recover them -- residual 1.5 relative, with no
         * simple rational coefficients -- so the closed form involves terms that
         * basis does not span. Guessing it is exactly the failure mode the
         * original note warned about. (BesselJ treats its order the same way.)
         *
         * Consequence worth knowing: Interval[] threading for the incomplete
         * EllipticPi in its m slot therefore stays symbolic, since the general
         * certifier differentiates to get a sign. --- */
        if (h == SYM_EllipticPi && n == 3) {
            Expr* N   = args[0];
            Expr* PHI = args[1];
            Expr* M   = args[2];
            Expr* dPHI = deriv_of(PHI, x, nonconsts);
            Expr* dN   = deriv_of(N, x, nonconsts);
            Expr* dM   = deriv_of(M, x, nonconsts);
            Expr* terms[3];
            size_t nt = 0;
            if (!is_lit_zero(dPHI)) {
                Expr* s2 = mk_fn2("Power", mk_fn1("Sin", expr_copy(PHI)), mk_int(2));
                Expr* rad = mk_fn2("Plus", mk_int(1),
                                mk_neg(mk_fn2("Times", expr_copy(M), expr_copy(s2))));
                Expr* circ = mk_fn2("Plus", mk_int(1),
                                mk_neg(mk_fn2("Times", expr_copy(N), s2)));
                Expr* den = mk_fn2("Times", circ, mk_fn1("Sqrt", rad));
                terms[nt++] = mk_fn2("Times", mk_fn2("Power", den, mk_int(-1)), dPHI);
            } else expr_free(dPHI);
            if (!is_lit_zero(dN)) {
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(1), mk_int(0), mk_int(0) }, 3);
                Expr* og = mk_fn_head1(op, mk_sym("EllipticPi"));
                terms[nt++] = mk_fn2("Times",
                    expr_new_function(og, (Expr*[]){ expr_copy(N), expr_copy(PHI), expr_copy(M) }, 3), dN);
            } else expr_free(dN);
            if (!is_lit_zero(dM)) {
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(0), mk_int(0), mk_int(1) }, 3);
                Expr* og = mk_fn_head1(op, mk_sym("EllipticPi"));
                terms[nt++] = mk_fn2("Times",
                    expr_new_function(og, (Expr*[]){ expr_copy(N), expr_copy(PHI), expr_copy(M) }, 3), dM);
            } else expr_free(dM);
            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            Expr* acc = mk_fn2("Plus", terms[0], terms[1]);
            for (size_t i = 2; i < nt; i++) acc = mk_fn2("Plus", acc, terms[i]);
            return acc;
        }

        if ((h == SYM_EllipticF || h == SYM_EllipticE) && n == 2) {
            bool is_F = (h == SYM_EllipticF);
            Expr* PHI = args[0];
            Expr* M   = args[1];
            Expr* dPHI = deriv_of(PHI, x, nonconsts);
            Expr* dM   = deriv_of(M, x, nonconsts);
            Expr* terms[2];
            size_t nt = 0;
            if (!is_lit_zero(dPHI)) {
                Expr* s2 = mk_fn2("Power", mk_fn1("Sin", expr_copy(PHI)), mk_int(2));
                Expr* rad = mk_fn2("Plus", mk_int(1),
                                mk_neg(mk_fn2("Times", expr_copy(M), s2)));
                Expr* root = mk_fn1("Sqrt", rad);
                /* F: 1/Sqrt[...]   E: Sqrt[...] */
                Expr* dfdphi = is_F ? mk_fn2("Power", root, mk_int(-1)) : root;
                terms[nt++] = mk_fn2("Times", dfdphi, dPHI);
            } else expr_free(dPHI);
            if (!is_lit_zero(dM)) {
                if (is_F) {
                    /* d/dm F(phi|m) =
                     *     -E(phi,m)/(2(m-1)m) - F(phi,m)/(2m)
                     *   + Sin[2 phi] / (4 (m-1) Sqrt[1 - m Sin[phi]^2])
                     * Verified against a central difference at 30 digits to
                     * ~1e-30 at (phi, m) = (0.7, 0.3), (1.1, 0.6), (0.4, 0.85).
                     * The boundary term is what distinguishes this from the
                     * complete case; it vanishes at phi = Pi/2, where the
                     * identity reduces to d/dm K's. */
                    Expr* s2  = mk_fn2("Power", mk_fn1("Sin", expr_copy(PHI)), mk_int(2));
                    Expr* rad = mk_fn2("Plus", mk_int(1),
                                    mk_neg(mk_fn2("Times", expr_copy(M), s2)));
                    Expr* mm1 = mk_fn2("Plus", expr_copy(M), mk_int(-1));
                    Expr* t1 = mk_neg(mk_fn2("Times",
                        mk_fn2("EllipticE", expr_copy(PHI), expr_copy(M)),
                        mk_fn2("Power", mk_fn2("Times", mk_int(2),
                            mk_fn2("Times", expr_copy(mm1), expr_copy(M))), mk_int(-1))));
                    Expr* t2 = mk_neg(mk_fn2("Times",
                        mk_fn2("EllipticF", expr_copy(PHI), expr_copy(M)),
                        mk_fn2("Power", mk_fn2("Times", mk_int(2), expr_copy(M)), mk_int(-1))));
                    Expr* t3 = mk_fn2("Times",
                        mk_fn1("Sin", mk_fn2("Times", mk_int(2), expr_copy(PHI))),
                        mk_fn2("Power", mk_fn2("Times", mk_int(4),
                            mk_fn2("Times", mm1, mk_fn1("Sqrt", rad))), mk_int(-1)));
                    Expr* sum = mk_fn2("Plus", t1, mk_fn2("Plus", t2, t3));
                    terms[nt++] = mk_fn2("Times", sum, dM);
                } else {
                    /* (E[phi,m] - F[phi,m]) / (2 m) */
                    Expr* diff = mk_fn2("Plus",
                        mk_fn2("EllipticE", expr_copy(PHI), expr_copy(M)),
                        mk_neg(mk_fn2("EllipticF", expr_copy(PHI), expr_copy(M))));
                    Expr* den = mk_fn2("Times", mk_int(2), expr_copy(M));
                    terms[nt++] = mk_fn2("Times",
                        mk_fn2("Times", diff, mk_fn2("Power", den, mk_int(-1))), dM);
                }
            } else expr_free(dM);
            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return mk_fn2("Plus", terms[0], terms[1]);
        }

        if (h == SYM_EllipticPi && n == 2) {   /* complete: both closed form */
            /* d/dn Pi(n|m) = (n E(m) + (m-n) K(m) + (n^2-m) Pi(n|m))
             *                / (2 (m-n) (n-1) n)
             * d/dm Pi(n|m) = (E(m)/(m-1) + Pi(n|m)) / (2 (n-m))
             * Both verified against a central difference at 30 digits to ~1e-31
             * at (n, m) = (0.3, 0.5), (-2, 0.4), (0.8, 0.2). */
            Expr* N = args[0];
            Expr* M = args[1];
            Expr* dN = deriv_of(N, x, nonconsts);
            Expr* dM = deriv_of(M, x, nonconsts);
            Expr* terms[2];
            size_t nt = 0;
            if (!is_lit_zero(dN)) {
                Expr* mn  = mk_fn2("Plus", expr_copy(M), mk_neg(expr_copy(N)));
                Expr* num = mk_fn2("Plus",
                    mk_fn2("Times", expr_copy(N), mk_fn1("EllipticE", expr_copy(M))),
                    mk_fn2("Plus",
                        mk_fn2("Times", expr_copy(mn), mk_fn1("EllipticK", expr_copy(M))),
                        mk_fn2("Times",
                            mk_fn2("Plus", mk_fn2("Power", expr_copy(N), mk_int(2)),
                                           mk_neg(expr_copy(M))),
                            mk_fn2("EllipticPi", expr_copy(N), expr_copy(M)))));
                Expr* den = mk_fn2("Times", mk_int(2),
                    mk_fn2("Times", mn,
                        mk_fn2("Times", mk_fn2("Plus", expr_copy(N), mk_int(-1)),
                                        expr_copy(N))));
                terms[nt++] = mk_fn2("Times",
                    mk_fn2("Times", num, mk_fn2("Power", den, mk_int(-1))), dN);
            } else expr_free(dN);
            if (!is_lit_zero(dM)) {
                Expr* num = mk_fn2("Plus",
                    mk_fn2("Times", mk_fn1("EllipticE", expr_copy(M)),
                        mk_fn2("Power", mk_fn2("Plus", expr_copy(M), mk_int(-1)), mk_int(-1))),
                    mk_fn2("EllipticPi", expr_copy(N), expr_copy(M)));
                Expr* den = mk_fn2("Times", mk_int(2),
                    mk_fn2("Plus", expr_copy(N), mk_neg(expr_copy(M))));
                terms[nt++] = mk_fn2("Times",
                    mk_fn2("Times", num, mk_fn2("Power", den, mk_int(-1))), dM);
            } else expr_free(dM);
            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return mk_fn2("Plus", terms[0], terms[1]);
        }

        if (h == SYM_BesselJ && n == 2) {
            Expr* N = args[0];
            Expr* Z = args[1];
            Expr* dN = deriv_of(N, x, nonconsts);
            Expr* dZ = deriv_of(Z, x, nonconsts);
            Expr* terms[2];
            size_t nt = 0;

            if (!is_lit_zero(dZ)) {
                /* (BesselJ[N-1, Z] - BesselJ[N+1, Z]) / 2 * dZ */
                Expr* bm1 = mk_fn2("BesselJ",
                              mk_fn2("Plus", expr_copy(N), mk_int(-1)), expr_copy(Z));
                Expr* bp1 = mk_fn2("BesselJ",
                              mk_fn2("Plus", expr_copy(N), mk_int(1)), expr_copy(Z));
                Expr* diff = mk_fn2("Plus", bm1, mk_neg(bp1));
                Expr* dBdZ = mk_fn2("Times",
                              mk_fn2("Power", mk_int(2), mk_int(-1)), diff);
                terms[nt++] = mk_fn2("Times", dBdZ, dZ);
            } else {
                expr_free(dZ);
            }

            if (!is_lit_zero(dN)) {
                /* Derivative[1,0][BesselJ][N,Z] * dN */
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(1), mk_int(0) }, 2);
                Expr* op_g = mk_fn_head1(op, mk_sym("BesselJ"));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(N), expr_copy(Z) }, 2);
                terms[nt++] = mk_fn2("Times", applied, dN);
            } else {
                expr_free(dN);
            }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return mk_fn2("Plus", terms[0], terms[1]);
        }

        /* --- AngerJ / WeberE / KelvinBer / KelvinBei / KelvinKer / KelvinKei
         * [N, Z]: chain rule on both args. In Z:
         *   d/dZ J_N = (J_(N-1) - J_(N+1))/2,  d/dZ E_N = (E_(N-1) - E_(N+1))/2
         *     (DLMF 11.10.32-11.10.33 combined)
         *   with R = ber|ker, I = bei|kei (DLMF 10.63.1 halved and recombined):
         *   d/dZ R_N = (R_(N+1) - R_(N-1) + I_(N+1) - I_(N-1)) / (2 Sqrt[2])
         *   d/dZ I_N = (I_(N+1) - I_(N-1) - R_(N+1) + R_(N-1)) / (2 Sqrt[2])
         * In N: Derivative[1,0][head][N, Z] (no closed form in the order). */
        if ((h == SYM_AngerJ || h == SYM_WeberE || h == SYM_KelvinBer ||
             h == SYM_KelvinBei || h == SYM_KelvinKer || h == SYM_KelvinKei) && n == 2) {
            const char* hn = h == SYM_AngerJ ? "AngerJ" : h == SYM_WeberE ? "WeberE"
                           : h == SYM_KelvinBer ? "KelvinBer" : h == SYM_KelvinBei ? "KelvinBei"
                           : h == SYM_KelvinKer ? "KelvinKer" : "KelvinKei";
            Expr* N = args[0];
            Expr* Z = args[1];
            Expr* dN = deriv_of(N, x, nonconsts);
            Expr* dZ = deriv_of(Z, x, nonconsts);
            Expr* terms[2];
            size_t nt = 0;

            if (!is_lit_zero(dZ)) {
                Expr* dSdZ;
                if (h == SYM_AngerJ || h == SYM_WeberE) {
                    Expr* lo = mk_fn2(hn, mk_fn2("Plus", mk_int(-1), expr_copy(N)), expr_copy(Z));
                    Expr* hi = mk_fn2(hn, mk_fn2("Plus", mk_int(1), expr_copy(N)), expr_copy(Z));
                    dSdZ = mk_fn2("Times", mk_fn2("Rational", mk_int(1), mk_int(2)),
                                  mk_fn2("Plus", lo, mk_neg(hi)));
                } else {
                    bool ker = h == SYM_KelvinKer || h == SYM_KelvinKei;
                    bool imag = h == SYM_KelvinBei || h == SYM_KelvinKei;
                    const char* R = ker ? "KelvinKer" : "KelvinBer";
                    const char* Ih = ker ? "KelvinKei" : "KelvinBei";
                    Expr* Rhi = mk_fn2(R, mk_fn2("Plus", mk_int(1), expr_copy(N)), expr_copy(Z));
                    Expr* Rlo = mk_fn2(R, mk_fn2("Plus", mk_int(-1), expr_copy(N)), expr_copy(Z));
                    Expr* Ihi = mk_fn2(Ih, mk_fn2("Plus", mk_int(1), expr_copy(N)), expr_copy(Z));
                    Expr* Ilo = mk_fn2(Ih, mk_fn2("Plus", mk_int(-1), expr_copy(N)), expr_copy(Z));
                    Expr* sum = imag
                        ? expr_new_function(mk_sym("Plus"), (Expr*[]){
                              Ihi, mk_neg(Ilo), mk_neg(Rhi), Rlo }, 4)
                        : expr_new_function(mk_sym("Plus"), (Expr*[]){
                              Rhi, mk_neg(Rlo), Ihi, mk_neg(Ilo) }, 4);
                    /* 1/(2 Sqrt[2]) = 2^(-3/2) */
                    dSdZ = mk_fn2("Times", mk_fn2("Power", mk_int(2),
                                  mk_fn2("Rational", mk_int(-3), mk_int(2))), sum);
                }
                terms[nt++] = mk_fn2("Times", dSdZ, dZ);
            } else {
                expr_free(dZ);
            }

            if (!is_lit_zero(dN)) {
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(1), mk_int(0) }, 2);
                Expr* op_g = mk_fn_head1(op, mk_sym(hn));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(N), expr_copy(Z) }, 2);
                terms[nt++] = mk_fn2("Times", applied, dN);
            } else {
                expr_free(dN);
            }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return mk_fn2("Plus", terms[0], terms[1]);
        }

        /* --- StruveH[N, Z] / StruveL[N, Z]: chain rule on both args.
         *   d/dZ H_N(Z) = ((Z/2)^N / (Sqrt[Pi] Gamma[N + 3/2])
         *                  + H_(N-1)(Z) - H_(N+1)(Z)) / 2
         *   d/dZ L_N(Z) = ((Z/2)^N / (Sqrt[Pi] Gamma[N + 3/2])
         *                  + L_(N-1)(Z) + L_(N+1)(Z)) / 2
         * (half the sum of the two first-order recurrences, DLMF 11.4.23-26).
         *   d/dN = Derivative[1,0][head][N, Z]  (no closed form in the order).
         * Zero-derivative arms are dropped, so D[StruveH[nu,x],x] keeps only the
         * dZ term, the form the reference system returns. */
        if ((h == SYM_StruveH || h == SYM_StruveL) && n == 2) {
            const char* hn = h == SYM_StruveH ? "StruveH" : "StruveL";
            Expr* N = args[0];
            Expr* Z = args[1];
            Expr* dN = deriv_of(N, x, nonconsts);
            Expr* dZ = deriv_of(Z, x, nonconsts);
            Expr* terms[2];
            size_t nt = 0;

            if (!is_lit_zero(dZ)) {
                /* 2^-N Z^N Pi^(-1/2) / Gamma[N + 3/2] */
                Expr* src = expr_new_function(mk_sym("Times"), (Expr*[]){
                    mk_fn2("Power", mk_int(2), mk_neg(expr_copy(N))),
                    mk_fn2("Power", expr_copy(Z), expr_copy(N)),
                    mk_fn2("Power", mk_sym("Pi"),
                           mk_fn2("Rational", mk_int(-1), mk_int(2))),
                    mk_fn2("Power",
                           mk_fn1("Gamma", mk_fn2("Plus",
                               mk_fn2("Rational", mk_int(3), mk_int(2)), expr_copy(N))),
                           mk_int(-1)) }, 4);
                Expr* sm1 = mk_fn2(hn, mk_fn2("Plus", mk_int(-1), expr_copy(N)),
                                   expr_copy(Z));
                Expr* sp1 = mk_fn2(hn, mk_fn2("Plus", mk_int(1), expr_copy(N)),
                                   expr_copy(Z));
                if (h == SYM_StruveH) sp1 = mk_neg(sp1);
                Expr* sum = expr_new_function(mk_sym("Plus"),
                                              (Expr*[]){ src, sm1, sp1 }, 3);
                Expr* dSdZ = mk_fn2("Times",
                              mk_fn2("Rational", mk_int(1), mk_int(2)), sum);
                terms[nt++] = mk_fn2("Times", dSdZ, dZ);
            } else {
                expr_free(dZ);
            }

            if (!is_lit_zero(dN)) {
                /* Derivative[1,0][head][N,Z] * dN */
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(1), mk_int(0) }, 2);
                Expr* op_g = mk_fn_head1(op, mk_sym(hn));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(N), expr_copy(Z) }, 2);
                terms[nt++] = mk_fn2("Times", applied, dN);
            } else {
                expr_free(dN);
            }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return mk_fn2("Plus", terms[0], terms[1]);
        }

        /* --- HypergeometricU[A, B, Z]: chain rule on the argument Z (DLMF 13.3.22,
         *   d/dz U(a,b,z) = -a U(a+1, b+1, z)).
         *   d/dA, d/dB have no elementary form -> inert
         *     Derivative[1,0,0][HypergeometricU][A,B,Z] / Derivative[0,1,0][...].
         * Zero-derivative arms are dropped, so D[HypergeometricU[a,b,z],z] keeps
         * only the dZ term. */
        if (h == SYM_HypergeometricU && n == 3) {
            Expr* A = args[0];
            Expr* B = args[1];
            Expr* Z = args[2];
            Expr* dA = deriv_of(A, x, nonconsts);
            Expr* dB = deriv_of(B, x, nonconsts);
            Expr* dZ = deriv_of(Z, x, nonconsts);
            Expr* terms[3];
            size_t nt = 0;

            if (!is_lit_zero(dZ)) {
                /* -A U(A+1, B+1, Z) */
                Expr* u1 = expr_new_function(mk_sym("HypergeometricU"),
                              (Expr*[]){ mk_fn2("Plus", mk_int(1), expr_copy(A)),
                                         mk_fn2("Plus", mk_int(1), expr_copy(B)),
                                         expr_copy(Z) }, 3);
                Expr* dUdZ = mk_fn2("Times", mk_neg(expr_copy(A)), u1);
                terms[nt++] = mk_fn2("Times", dUdZ, dZ);
            } else {
                expr_free(dZ);
            }

            if (!is_lit_zero(dA)) {
                /* Derivative[1,0,0][HypergeometricU][A,B,Z] * dA */
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(1), mk_int(0), mk_int(0) }, 3);
                Expr* op_g = mk_fn_head1(op, mk_sym("HypergeometricU"));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(A), expr_copy(B), expr_copy(Z) }, 3);
                terms[nt++] = mk_fn2("Times", applied, dA);
            } else {
                expr_free(dA);
            }

            if (!is_lit_zero(dB)) {
                /* Derivative[0,1,0][HypergeometricU][A,B,Z] * dB */
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(0), mk_int(1), mk_int(0) }, 3);
                Expr* op_g = mk_fn_head1(op, mk_sym("HypergeometricU"));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(A), expr_copy(B), expr_copy(Z) }, 3);
                terms[nt++] = mk_fn2("Times", applied, dB);
            } else {
                expr_free(dB);
            }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            if (nt == 2) return mk_fn2("Plus", terms[0], terms[1]);
            return expr_new_function(mk_sym("Plus"), terms, 3);
        }

        /* --- MathieuC[A, Q, Z] / MathieuCPrime[A, Q, Z]: the Mathieu equation
         *   ce'' = (2 Q Cos[2 Z] - A) ce, so
         *   d/dZ MathieuC      = MathieuCPrime,
         *   d/dZ MathieuCPrime = (-A + 2 Q Cos[2 Z]) MathieuC.
         *   d/dA, d/dQ have no closed form -> inert
         *     Derivative[1,0,0][head][A,Q,Z] / Derivative[0,1,0][head][A,Q,Z].
         * Zero-derivative arms are dropped, so D[MathieuC[a,q,x],x] keeps only
         * the dZ term (and Series about any point follows by Taylor). */
        /*   MathieuS / MathieuSPrime: the same equation for se,
         *   d/dZ MathieuS = MathieuSPrime,
         *   d/dZ MathieuSPrime = (-A + 2 Q Cos[2 Z]) MathieuS. */
        if ((h == SYM_MathieuC || h == SYM_MathieuCPrime ||
             h == SYM_MathieuS || h == SYM_MathieuSPrime) && n == 3) {
            const int odd = (h == SYM_MathieuS || h == SYM_MathieuSPrime);
            const int is_prime = (h == SYM_MathieuCPrime || h == SYM_MathieuSPrime);
            const char* hn = h == SYM_MathieuC ? "MathieuC" :
                             h == SYM_MathieuCPrime ? "MathieuCPrime" :
                             h == SYM_MathieuS ? "MathieuS" : "MathieuSPrime";
            const char* base = odd ? "MathieuS" : "MathieuC";
            const char* prime = odd ? "MathieuSPrime" : "MathieuCPrime";
            Expr* A = args[0];
            Expr* Q = args[1];
            Expr* Z = args[2];
            Expr* dA = deriv_of(A, x, nonconsts);
            Expr* dQ = deriv_of(Q, x, nonconsts);
            Expr* dZ = deriv_of(Z, x, nonconsts);
            Expr* terms[3];
            size_t nt = 0;

            if (!is_lit_zero(dZ)) {
                Expr* dfdz;
                if (!is_prime) {
                    dfdz = expr_new_function(mk_sym(prime),
                              (Expr*[]){ expr_copy(A), expr_copy(Q), expr_copy(Z) }, 3);
                } else {
                    /* (-A + 2 Q Cos[2 Z]) MathieuC[A, Q, Z] */
                    Expr* coef = mk_fn2("Plus", mk_neg(expr_copy(A)),
                                   expr_new_function(mk_sym("Times"), (Expr*[]){
                                       mk_int(2), expr_copy(Q),
                                       mk_fn1("Cos", mk_fn2("Times", mk_int(2), expr_copy(Z)))
                                   }, 3));
                    Expr* ce = expr_new_function(mk_sym(base),
                              (Expr*[]){ expr_copy(A), expr_copy(Q), expr_copy(Z) }, 3);
                    dfdz = mk_fn2("Times", coef, ce);
                }
                terms[nt++] = mk_fn2("Times", dfdz, dZ);
            } else {
                expr_free(dZ);
            }
            for (int which = 0; which < 2; which++) {
                Expr* dP = which == 0 ? dA : dQ;
                if (is_lit_zero(dP)) { expr_free(dP); continue; }
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(which == 0 ? 1 : 0), mk_int(which == 1 ? 1 : 0),
                                         mk_int(0) }, 3);
                Expr* op_g = mk_fn_head1(op, mk_sym(hn));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(A), expr_copy(Q), expr_copy(Z) }, 3);
                terms[nt++] = mk_fn2("Times", applied, dP);
            }
            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            if (nt == 2) return mk_fn2("Plus", terms[0], terms[1]);
            return expr_new_function(mk_sym("Plus"), terms, 3);
        }

        /* --- WhittakerM[K, M, Z]: chain rule on the argument Z (DLMF 13.4.27,
         *   z M' = (z/2 - k) M + (k + m + 1/2) M_{k+1,m}):
         *   d/dZ M_{k,m}(Z) = (1/2 - k/Z) M_{k,m}(Z) + ((1/2+k+m)/Z) M_{k+1,m}(Z).
         *   d/dK, d/dM have no elementary form -> inert
         *     Derivative[1,0,0][WhittakerM][K,M,Z] / Derivative[0,1,0][...].
         * Zero-derivative arms are dropped, so D[WhittakerM[k,m,z],z] keeps only
         * the dZ term. */
        if (h == SYM_WhittakerM && n == 3) {
            Expr* K = args[0];
            Expr* M = args[1];
            Expr* Z = args[2];
            Expr* dK = deriv_of(K, x, nonconsts);
            Expr* dM = deriv_of(M, x, nonconsts);
            Expr* dZ = deriv_of(Z, x, nonconsts);
            Expr* terms[3];
            size_t nt = 0;

            if (!is_lit_zero(dZ)) {
                Expr* invZ = mk_fn2("Power", expr_copy(Z), mk_int(-1));   /* 1/Z */
                /* (1/2 - k/Z) M_{k,m}(Z) */
                Expr* c0 = mk_fn2("Plus",
                              mk_fn2("Rational", mk_int(1), mk_int(2)),
                              mk_neg(mk_fn2("Times", expr_copy(K), expr_copy(invZ))));
                Expr* wkm = expr_new_function(mk_sym("WhittakerM"),
                              (Expr*[]){ expr_copy(K), expr_copy(M), expr_copy(Z) }, 3);
                Expr* t0 = mk_fn2("Times", c0, wkm);
                /* ((1/2 + k + m)/Z) M_{k+1,m}(Z) */
                Expr* c1 = mk_fn2("Times",
                              mk_fn2("Plus", mk_fn2("Rational", mk_int(1), mk_int(2)),
                                     mk_fn2("Plus", expr_copy(K), expr_copy(M))),
                              invZ);                           /* original invZ adopted */
                Expr* wk1 = expr_new_function(mk_sym("WhittakerM"),
                              (Expr*[]){ mk_fn2("Plus", mk_int(1), expr_copy(K)),
                                         expr_copy(M), expr_copy(Z) }, 3);
                Expr* t1 = mk_fn2("Times", c1, wk1);
                Expr* dMdZ = mk_fn2("Plus", t0, t1);
                terms[nt++] = mk_fn2("Times", dMdZ, dZ);
            } else {
                expr_free(dZ);
            }

            if (!is_lit_zero(dK)) {
                /* Derivative[1,0,0][WhittakerM][K,M,Z] * dK */
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(1), mk_int(0), mk_int(0) }, 3);
                Expr* op_g = mk_fn_head1(op, mk_sym("WhittakerM"));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(K), expr_copy(M), expr_copy(Z) }, 3);
                terms[nt++] = mk_fn2("Times", applied, dK);
            } else {
                expr_free(dK);
            }

            if (!is_lit_zero(dM)) {
                /* Derivative[0,1,0][WhittakerM][K,M,Z] * dM */
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(0), mk_int(1), mk_int(0) }, 3);
                Expr* op_g = mk_fn_head1(op, mk_sym("WhittakerM"));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(K), expr_copy(M), expr_copy(Z) }, 3);
                terms[nt++] = mk_fn2("Times", applied, dM);
            } else {
                expr_free(dM);
            }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            if (nt == 2) return mk_fn2("Plus", terms[0], terms[1]);
            return expr_new_function(mk_sym("Plus"), terms, 3);
        }

        /* --- WhittakerW[K, M, Z]: chain rule on the argument Z (DLMF 13.14.31,
         *   d/dz W_{k,m}(z) = (1/2 - k/z) W_{k,m}(z) - (1/z) W_{k+1,m}(z)).
         *   d/dK, d/dM have no elementary form -> inert
         *     Derivative[1,0,0][WhittakerW][K,M,Z] / Derivative[0,1,0][...]. */
        if (h == SYM_WhittakerW && n == 3) {
            Expr* K = args[0];
            Expr* M = args[1];
            Expr* Z = args[2];
            Expr* dK = deriv_of(K, x, nonconsts);
            Expr* dM = deriv_of(M, x, nonconsts);
            Expr* dZ = deriv_of(Z, x, nonconsts);
            Expr* terms[3];
            size_t nt = 0;

            if (!is_lit_zero(dZ)) {
                Expr* invZ = mk_fn2("Power", expr_copy(Z), mk_int(-1));   /* 1/Z */
                /* (1/2 - k/Z) W_{k,m}(Z) */
                Expr* c0 = mk_fn2("Plus",
                              mk_fn2("Rational", mk_int(1), mk_int(2)),
                              mk_neg(mk_fn2("Times", expr_copy(K), expr_copy(invZ))));
                Expr* wkm = expr_new_function(mk_sym("WhittakerW"),
                              (Expr*[]){ expr_copy(K), expr_copy(M), expr_copy(Z) }, 3);
                Expr* t0 = mk_fn2("Times", c0, wkm);
                /* -(1/Z) W_{k+1,m}(Z) */
                Expr* wk1 = expr_new_function(mk_sym("WhittakerW"),
                              (Expr*[]){ mk_fn2("Plus", mk_int(1), expr_copy(K)),
                                         expr_copy(M), expr_copy(Z) }, 3);
                Expr* t1 = mk_fn2("Times", mk_neg(invZ), wk1);   /* original invZ adopted */
                Expr* dWdZ = mk_fn2("Plus", t0, t1);
                terms[nt++] = mk_fn2("Times", dWdZ, dZ);
            } else {
                expr_free(dZ);
            }

            if (!is_lit_zero(dK)) {
                /* Derivative[1,0,0][WhittakerW][K,M,Z] * dK */
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(1), mk_int(0), mk_int(0) }, 3);
                Expr* op_g = mk_fn_head1(op, mk_sym("WhittakerW"));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(K), expr_copy(M), expr_copy(Z) }, 3);
                terms[nt++] = mk_fn2("Times", applied, dK);
            } else {
                expr_free(dK);
            }

            if (!is_lit_zero(dM)) {
                /* Derivative[0,1,0][WhittakerW][K,M,Z] * dM */
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(0), mk_int(1), mk_int(0) }, 3);
                Expr* op_g = mk_fn_head1(op, mk_sym("WhittakerW"));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(K), expr_copy(M), expr_copy(Z) }, 3);
                terms[nt++] = mk_fn2("Times", applied, dM);
            } else {
                expr_free(dM);
            }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            if (nt == 2) return mk_fn2("Plus", terms[0], terms[1]);
            return expr_new_function(mk_sym("Plus"), terms, 3);
        }

        /* --- LegendreQ[N, Z]: chain rule on the argument only.
         *   d/dZ Q_N(Z) = ((-1-N) Z Q_N(Z) + (1+N) Q_{N+1}(Z)) / (Z^2 - 1).
         *   d/dN Q_N(Z) = Derivative[1,0][LegendreQ][N,Z]  (no elementary form).
         * Zero-derivative arms are dropped, so D[LegendreQ[n,x],x] keeps only the
         * dZ term (LegendreP has the identical recurrence). */
        if (h == SYM_LegendreQ && n == 2) {
            Expr* N = args[0];
            Expr* Z = args[1];
            Expr* dN = deriv_of(N, x, nonconsts);
            Expr* dZ = deriv_of(Z, x, nonconsts);
            Expr* terms[2];
            size_t nt = 0;

            if (!is_lit_zero(dZ)) {
                /* (-1-N) Z Q_N */
                Expr* c0 = mk_fn2("Plus", mk_int(-1), mk_neg(expr_copy(N)));
                Expr* qN = mk_fn2("LegendreQ", expr_copy(N), expr_copy(Z));
                Expr* t0 = mk_fn2("Times", c0, mk_fn2("Times", expr_copy(Z), qN));
                /* (1+N) Q_{N+1} */
                Expr* c1 = mk_fn2("Plus", mk_int(1), expr_copy(N));
                Expr* qN1 = mk_fn2("LegendreQ",
                                mk_fn2("Plus", mk_int(1), expr_copy(N)), expr_copy(Z));
                Expr* t1 = mk_fn2("Times", c1, qN1);
                Expr* num = mk_fn2("Plus", t0, t1);
                /* Z^2 - 1 */
                Expr* den = mk_fn2("Plus", mk_int(-1),
                                mk_fn2("Power", expr_copy(Z), mk_int(2)));
                Expr* dQdZ = mk_fn2("Times", num, mk_fn2("Power", den, mk_int(-1)));
                terms[nt++] = mk_fn2("Times", dQdZ, dZ);
            } else {
                expr_free(dZ);
            }

            if (!is_lit_zero(dN)) {
                /* Derivative[1,0][LegendreQ][N,Z] * dN */
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(1), mk_int(0) }, 2);
                Expr* op_g = mk_fn_head1(op, mk_sym("LegendreQ"));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(N), expr_copy(Z) }, 2);
                terms[nt++] = mk_fn2("Times", applied, dN);
            } else {
                expr_free(dN);
            }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return mk_fn2("Plus", terms[0], terms[1]);
        }

        /* --- LegendreP[N, Z]: identical recurrence to LegendreQ.
         *   d/dZ P_N(Z) = ((-1-N) Z P_N(Z) + (1+N) P_{N+1}(Z)) / (Z^2 - 1).
         *   d/dN P_N(Z) = Derivative[1,0][LegendreP][N,Z]  (no elementary form).
         * Needed so the Legendre-equation solutions DSolve emits back-substitution
         * verify (LegendreQ already had this; LegendreP did not). */
        if (h == SYM_LegendreP && n == 2) {
            Expr* N = args[0];
            Expr* Z = args[1];
            Expr* dN = deriv_of(N, x, nonconsts);
            Expr* dZ = deriv_of(Z, x, nonconsts);
            Expr* terms[2];
            size_t nt = 0;

            if (!is_lit_zero(dZ)) {
                /* (-1-N) Z P_N */
                Expr* c0 = mk_fn2("Plus", mk_int(-1), mk_neg(expr_copy(N)));
                Expr* pN = mk_fn2("LegendreP", expr_copy(N), expr_copy(Z));
                Expr* t0 = mk_fn2("Times", c0, mk_fn2("Times", expr_copy(Z), pN));
                /* (1+N) P_{N+1} */
                Expr* c1 = mk_fn2("Plus", mk_int(1), expr_copy(N));
                Expr* pN1 = mk_fn2("LegendreP",
                                mk_fn2("Plus", mk_int(1), expr_copy(N)), expr_copy(Z));
                Expr* t1 = mk_fn2("Times", c1, pN1);
                Expr* num = mk_fn2("Plus", t0, t1);
                /* Z^2 - 1 */
                Expr* den = mk_fn2("Plus", mk_int(-1),
                                mk_fn2("Power", expr_copy(Z), mk_int(2)));
                Expr* dPdZ = mk_fn2("Times", num, mk_fn2("Power", den, mk_int(-1)));
                terms[nt++] = mk_fn2("Times", dPdZ, dZ);
            } else {
                expr_free(dZ);
            }

            if (!is_lit_zero(dN)) {
                /* Derivative[1,0][LegendreP][N,Z] * dN */
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(1), mk_int(0) }, 2);
                Expr* op_g = mk_fn_head1(op, mk_sym("LegendreP"));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(N), expr_copy(Z) }, 2);
                terms[nt++] = mk_fn2("Times", applied, dN);
            } else {
                expr_free(dN);
            }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return mk_fn2("Plus", terms[0], terms[1]);
        }

        /* --- ChebyshevU[N, Z]:
         *   d/dZ U_N(Z) = ((-1-N) U_(N-1)(Z) + N Z U_N(Z)) / (Z^2 - 1).
         *   d/dN U_N(Z) = Derivative[1,0][ChebyshevU][N,Z]  (no elementary form).
         * Enables D[] and, via it, the generic Interval threading in eval.c. */
        if (h == SYM_ChebyshevU && n == 2) {
            Expr* N = args[0];
            Expr* Z = args[1];
            Expr* dN = deriv_of(N, x, nonconsts);
            Expr* dZ = deriv_of(Z, x, nonconsts);
            Expr* terms[2];
            size_t nt = 0;

            if (!is_lit_zero(dZ)) {
                /* (-1-N) U_{N-1}(Z) */
                Expr* c0 = mk_fn2("Plus", mk_int(-1), mk_neg(expr_copy(N)));
                Expr* uNm1 = mk_fn2("ChebyshevU",
                                mk_fn2("Plus", expr_copy(N), mk_int(-1)), expr_copy(Z));
                Expr* t0 = mk_fn2("Times", c0, uNm1);
                /* N Z U_N(Z) */
                Expr* uN = mk_fn2("ChebyshevU", expr_copy(N), expr_copy(Z));
                Expr* t1 = mk_fn2("Times", expr_copy(N),
                                mk_fn2("Times", expr_copy(Z), uN));
                Expr* num = mk_fn2("Plus", t0, t1);
                /* Z^2 - 1 */
                Expr* den = mk_fn2("Plus", mk_int(-1),
                                mk_fn2("Power", expr_copy(Z), mk_int(2)));
                Expr* dUdZ = mk_fn2("Times", num, mk_fn2("Power", den, mk_int(-1)));
                terms[nt++] = mk_fn2("Times", dUdZ, dZ);
            } else {
                expr_free(dZ);
            }

            if (!is_lit_zero(dN)) {
                /* Derivative[1,0][ChebyshevU][N,Z] * dN */
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(1), mk_int(0) }, 2);
                Expr* op_g = mk_fn_head1(op, mk_sym("ChebyshevU"));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(N), expr_copy(Z) }, 2);
                terms[nt++] = mk_fn2("Times", applied, dN);
            } else {
                expr_free(dN);
            }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return mk_fn2("Plus", terms[0], terms[1]);
        }

        /* --- ChebyshevT[N, Z]:
         *   d/dZ T_N(Z) = N U_(N-1)(Z)  (first kind differentiates to the
         *                                second kind).
         *   d/dN T_N(Z) = Derivative[1,0][ChebyshevT][N,Z]  (no elementary form).
         * Enables D[] and, via it, the generic Interval threading in eval.c. */
        if (h == SYM_ChebyshevT && n == 2) {
            Expr* N = args[0];
            Expr* Z = args[1];
            Expr* dN = deriv_of(N, x, nonconsts);
            Expr* dZ = deriv_of(Z, x, nonconsts);
            Expr* terms[2];
            size_t nt = 0;

            if (!is_lit_zero(dZ)) {
                /* N U_{N-1}(Z) */
                Expr* uNm1 = mk_fn2("ChebyshevU",
                                mk_fn2("Plus", expr_copy(N), mk_int(-1)), expr_copy(Z));
                Expr* dTdZ = mk_fn2("Times", expr_copy(N), uNm1);
                terms[nt++] = mk_fn2("Times", dTdZ, dZ);
            } else {
                expr_free(dZ);
            }

            if (!is_lit_zero(dN)) {
                /* Derivative[1,0][ChebyshevT][N,Z] * dN */
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(1), mk_int(0) }, 2);
                Expr* op_g = mk_fn_head1(op, mk_sym("ChebyshevT"));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(N), expr_copy(Z) }, 2);
                terms[nt++] = mk_fn2("Times", applied, dN);
            } else {
                expr_free(dN);
            }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return mk_fn2("Plus", terms[0], terms[1]);
        }

        /* --- HermiteH[N, Z]:
         *   d/dZ H_N(Z) = 2 N HermiteH[N-1, Z].
         *   d/dN H_N(Z) = Derivative[1,0][HermiteH][N,Z]  (no elementary form).
         * Enables D[], the generic Interval threading in eval.c, and -- via the
         * Taylor-via-D path in series.c -- the finite-point Series expansion. */
        /* --- ParabolicCylinderD[N, Z] (DLMF 12.8.2):
         *   d/dZ D_N(Z) = Z/2 D_N(Z) - D_{N+1}(Z);
         *   d/dN is the inert Derivative[1,0][ParabolicCylinderD][N, Z]. */
        if (h == SYM_ParabolicCylinderD && n == 2) {
            Expr* N = args[0];
            Expr* Z = args[1];
            Expr* dN = deriv_of(N, x, nonconsts);
            Expr* dZ = deriv_of(Z, x, nonconsts);
            Expr* terms[2];
            size_t nt = 0;
            if (!is_lit_zero(dZ)) {
                Expr* d0 = mk_fn2("ParabolicCylinderD", expr_copy(N), expr_copy(Z));
                Expr* d1 = mk_fn2("ParabolicCylinderD",
                                  mk_fn2("Plus", expr_copy(N), mk_int(1)), expr_copy(Z));
                Expr* half = expr_new_function(mk_sym("Times"),
                                (Expr*[]){ mk_fn2("Rational", mk_int(1), mk_int(2)),
                                           expr_copy(Z), d0 }, 3);
                Expr* dd = mk_fn2("Plus", half, mk_fn2("Times", mk_int(-1), d1));
                terms[nt++] = mk_fn2("Times", dd, dZ);
            } else {
                expr_free(dZ);
            }
            if (!is_lit_zero(dN)) {
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(1), mk_int(0) }, 2);
                Expr* op_g = mk_fn_head1(op, mk_sym("ParabolicCylinderD"));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(N), expr_copy(Z) }, 2);
                terms[nt++] = mk_fn2("Times", applied, dN);
            } else {
                expr_free(dN);
            }
            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return mk_fn2("Plus", terms[0], terms[1]);
        }

        /* --- BetaRegularized[Z, A, B] and [Z0, Z1, A, B] (DLMF 8.17.1): the
         * z-derivative z^(a-1) (1-z)^(b-1)/Beta[a, b]; the a-derivative
         *   -z^a 3F2(a, a, 1-b; a+1, a+1; z)/(a^2 Beta[a, b])
         *     + I_z(a, b) (Log[z] - PolyGamma[a] + PolyGamma[a + b]),
         * and the b-derivative from I_z(a, b) = 1 - I_(1-z)(b, a). The
         * four-argument form differentiates as I_z1 - I_z0. */
        if ((n == 3 || n == 4) && h && strcmp(h, "BetaRegularized") == 0) {
            static const char* const TB[3] = {
                "((1 - #1)^(#3 - 1) #1^(#2 - 1)/Beta[#2, #3])&",
                "(-#1^#2 HypergeometricPFQ[{#2, #2, 1 - #3}, {1 + #2, 1 + #2}, #1]/(#2^2 Beta[#2, #3]) + "
                    "BetaRegularized[#1, #2, #3] (Log[#1] - PolyGamma[0, #2] + PolyGamma[0, #2 + #3]))&",
                "((1 - #1)^#3 HypergeometricPFQ[{#3, #3, 1 - #2}, {1 + #3, 1 + #3}, 1 - #1]/(#3^2 Beta[#2, #3]) - "
                    "(1 - BetaRegularized[#1, #2, #3]) (Log[1 - #1] - PolyGamma[0, #3] + PolyGamma[0, #2 + #3]))&" };
            Expr* A = args[n - 2];
            Expr* B = args[n - 1];
            Expr* terms[6];
            size_t nt = 0;
            for (int end = (n == 4 ? 0 : 1); end < 2; end++) {
                /* end 1: the upper (or only) endpoint, + ; end 0: Z0, - */
                Expr* Z = (n == 4) ? args[end] : args[0];
                double sg = end == 1 ? 1.0 : -1.0;
                Expr* d[3] = { deriv_of(Z, x, nonconsts), deriv_of(A, x, nonconsts),
                               deriv_of(B, x, nonconsts) };
                for (int j = 0; j < 3; j++) {
                    if (is_lit_zero(d[j])) { expr_free(d[j]); continue; }
                    Expr* f = parse_expression(TB[j]);
                    Expr* g = expr_new_function(f, (Expr*[]){ expr_copy(Z), expr_copy(A),
                                                              expr_copy(B) }, 3);
                    Expr* t = mk_fn2("Times", g, d[j]);
                    if (sg < 0) t = mk_fn2("Times", mk_int(-1), t);
                    terms[nt++] = t;
                }
            }
            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return expr_new_function(mk_sym("Plus"), terms, nt);
        }
        /* --- WeierstrassP / PPrime / Zeta / Sigma [U, {G2, G3}] (DLMF 23.3.10,
         * 23.2.8-9; the invariant derivatives from the homogeneity relations
         * of wp in (u, g2, g3), DLMF 23.10), templates in #1 = u, #2 = g2,
         * #3 = g3. */
        if (n == 2 && h && strncmp(h, "Weierstrass", 11) == 0 &&
            args[1]->type == EXPR_FUNCTION && args[1]->data.function.head->type == EXPR_SYMBOL &&
            strcmp(args[1]->data.function.head->data.symbol.name, "List") == 0 &&
            args[1]->data.function.arg_count == 2) {
            int k = strcmp(h, "WeierstrassP") == 0 ? 0 : strcmp(h, "WeierstrassPPrime") == 0 ? 1
                  : strcmp(h, "WeierstrassZeta") == 0 ? 2 : strcmp(h, "WeierstrassSigma") == 0 ? 3 : -1;
            if (k >= 0) {
                #define WPP "WeierstrassP[#1, {#2, #3}]"
                #define WPD "WeierstrassPPrime[#1, {#2, #3}]"
                #define WPZ "WeierstrassZeta[#1, {#2, #3}]"
                #define WPS "WeierstrassSigma[#1, {#2, #3}]"
                #define WPDEL "(#2^3 - 27 #3^2)"
                static const char* const T[4][3] = {
                    { "(" WPD ")&",
                      "((6 #2 #3 + 2 #2^2 " WPP " - 36 #3 " WPP "^2 + " WPD " (#2^2 #1 - 18 #3 " WPZ "))/(4 " WPDEL "))&",
                      "((-2 #2^2 - 18 #3 " WPP " + 12 #2 " WPP "^2 + " WPD " (-9 #3 #1 + 6 #2 " WPZ "))/(2 " WPDEL "))&" },
                    { "(6 " WPP "^2 - #2/2)&",
                      "((-(#2^3 #1) + 12 #2^2 #1 " WPP "^2 + 6 (#2^2 - 18 #3 " WPP ") " WPD " + 18 #3 (#2 - 12 " WPP "^2) " WPZ ")/(8 " WPDEL "))&",
                      "((9 #2 #3 #1 - 108 #3 #1 " WPP "^2 + (-54 #3 + 36 #2 " WPP ") " WPD " + (-6 #2^2 + 72 #2 " WPP "^2) " WPZ ")/(4 " WPDEL "))&" },
                    { "(-" WPP ")&",
                      "((-(#2 #1 (3 #3 + 2 #2 " WPP ")) + 18 #3 " WPD " + 2 (#2^2 + 18 #3 " WPP ") " WPZ ")/(8 " WPDEL "))&",
                      "((#1 (#2^2 + 18 #3 " WPP ") - 6 #2 " WPD " - 6 (3 #3 + 2 #2 " WPP ") " WPZ ")/(4 " WPDEL "))&" },
                    { "(" WPS " " WPZ ")&",
                      "((" WPS " (-3 #2 #3 #1^2 + 4 #2^2 (-1 + #1 " WPZ ") + 36 #3 (" WPP " - " WPZ "^2)))/(16 " WPDEL "))&",
                      "((" WPS " (#2^2 #1^2 - 12 #2 " WPP " + 12 #2 " WPZ "^2 - 36 #3 (-1 + #1 " WPZ ")))/(8 " WPDEL "))&" } };
                #undef WPP
                #undef WPD
                #undef WPZ
                #undef WPS
                #undef WPDEL
                Expr* U = args[0];
                Expr* G2 = args[1]->data.function.args[0];
                Expr* G3 = args[1]->data.function.args[1];
                Expr* d[3] = { deriv_of(U, x, nonconsts), deriv_of(G2, x, nonconsts),
                               deriv_of(G3, x, nonconsts) };
                Expr* terms[3];
                size_t nt = 0;
                for (int j = 0; j < 3; j++) {
                    if (is_lit_zero(d[j])) { expr_free(d[j]); continue; }
                    Expr* f = parse_expression(T[k][j]);
                    Expr* g = expr_new_function(f, (Expr*[]){ expr_copy(U), expr_copy(G2),
                                                              expr_copy(G3) }, 3);
                    terms[nt++] = mk_fn2("Times", g, d[j]);
                }
                if (nt == 0) return mk_int(0);
                if (nt == 1) return terms[0];
                return expr_new_function(mk_sym("Plus"), terms, nt);
            }
        }
        /* --- EllipticTheta[A, U, Q] / EllipticThetaPrime[A, U, Q]:
         *   d/dU EllipticTheta = EllipticThetaPrime; every other partial is the
         *   inert Derivative form (the reference's). A is a discrete index. */
        if (n == 3 && h && (strcmp(h, "EllipticTheta") == 0 || strcmp(h, "EllipticThetaPrime") == 0)) {
            bool pr = h[13] == 'P';
            Expr* A = args[0];
            Expr* U = args[1];
            Expr* Q = args[2];
            Expr* dU = deriv_of(U, x, nonconsts);
            Expr* dQ = deriv_of(Q, x, nonconsts);
            Expr* terms[2];
            size_t nt = 0;
            if (!is_lit_zero(dU)) {
                Expr* g;
                if (!pr) {
                    g = expr_new_function(mk_sym("EllipticThetaPrime"),
                        (Expr*[]){ expr_copy(A), expr_copy(U), expr_copy(Q) }, 3);
                } else {
                    Expr* op = expr_new_function(mk_sym("Derivative"),
                                  (Expr*[]){ mk_int(0), mk_int(1), mk_int(0) }, 3);
                    g = expr_new_function(mk_fn_head1(op, mk_sym(h)),
                        (Expr*[]){ expr_copy(A), expr_copy(U), expr_copy(Q) }, 3);
                }
                terms[nt++] = mk_fn2("Times", g, dU);
            } else {
                expr_free(dU);
            }
            if (!is_lit_zero(dQ)) {
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(0), mk_int(0), mk_int(1) }, 3);
                Expr* g = expr_new_function(mk_fn_head1(op, mk_sym(h)),
                    (Expr*[]){ expr_copy(A), expr_copy(U), expr_copy(Q) }, 3);
                terms[nt++] = mk_fn2("Times", g, dQ);
            } else {
                expr_free(dQ);
            }
            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return mk_fn2("Plus", terms[0], terms[1]);
        }
        /* --- HankelH1/H2 and the spherical Bessel family f[N, Z] (DLMF 10.6.1,
         * 10.51.2): H' = (H_{N-1} - H_{N+1})/2, f' = -f/(2Z) + (f_{N-1} - f_{N+1})/2;
         * the order derivative stays inert. */
        if (n == 2 && h && (strcmp(h, "HankelH1") == 0 || strcmp(h, "HankelH2") == 0 ||
                            strcmp(h, "SphericalBesselJ") == 0 || strcmp(h, "SphericalBesselY") == 0 ||
                            strcmp(h, "SphericalHankelH1") == 0 || strcmp(h, "SphericalHankelH2") == 0)) {
            bool sph = h[0] == 'S';
            {
                char tpl[256];
                if (sph)
                    snprintf(tpl, sizeof tpl, "(-%s[#1, #2]/(2 #2) + (%s[#1 - 1, #2] - %s[#1 + 1, #2])/2)&", h, h, h);
                else
                    snprintf(tpl, sizeof tpl, "((%s[#1 - 1, #2] - %s[#1 + 1, #2])/2)&", h, h);
                Expr* N = args[0];
                Expr* Z = args[1];
                Expr* dN = deriv_of(N, x, nonconsts);
                Expr* dZ = deriv_of(Z, x, nonconsts);
                Expr* terms[2];
                size_t nt = 0;
                if (!is_lit_zero(dZ)) {
                    Expr* f = parse_expression(tpl);
                    Expr* g = expr_new_function(f, (Expr*[]){ expr_copy(N), expr_copy(Z) }, 2);
                    terms[nt++] = mk_fn2("Times", g, dZ);
                } else {
                    expr_free(dZ);
                }
                if (!is_lit_zero(dN)) {
                    Expr* op = expr_new_function(mk_sym("Derivative"),
                                  (Expr*[]){ mk_int(1), mk_int(0) }, 2);
                    Expr* op_g = mk_fn_head1(op, mk_sym(h));
                    Expr* applied = expr_new_function(op_g,
                                  (Expr*[]){ expr_copy(N), expr_copy(Z) }, 2);
                    terms[nt++] = mk_fn2("Times", applied, dN);
                } else {
                    expr_free(dN);
                }
                if (nt == 0) return mk_int(0);
                if (nt == 1) return terms[0];
                return mk_fn2("Plus", terms[0], terms[1]);
            }
        }
        /* --- The Jacobi elliptic family, f[U, M]: d/dU from DLMF 22.13 (and the
         * EllipticF integrands for the inverses), as pure-function templates in
         * (U, M); d/dM stays the inert Derivative[0, 1][f][U, M]. */
        if (n == 2 && h && (strncmp(h, "Jacobi", 6) == 0 || strncmp(h, "InverseJacobi", 13) == 0)) {
            static const char* const JD[][2] = {
                { "JacobiSN", "(JacobiCN[#1, #2] JacobiDN[#1, #2])&" },
                { "JacobiCN", "(-JacobiSN[#1, #2] JacobiDN[#1, #2])&" },
                { "JacobiDN", "(-#2 JacobiSN[#1, #2] JacobiCN[#1, #2])&" },
                { "JacobiCD", "((-1 + #2) JacobiND[#1, #2] JacobiSD[#1, #2])&" },
                { "JacobiCS", "(-JacobiDS[#1, #2] JacobiNS[#1, #2])&" },
                { "JacobiDC", "((1 - #2) JacobiNC[#1, #2] JacobiSC[#1, #2])&" },
                { "JacobiDS", "(-JacobiCS[#1, #2] JacobiNS[#1, #2])&" },
                { "JacobiNC", "(JacobiDC[#1, #2] JacobiSC[#1, #2])&" },
                { "JacobiND", "(#2 JacobiCD[#1, #2] JacobiSD[#1, #2])&" },
                { "JacobiNS", "(-JacobiCS[#1, #2] JacobiDS[#1, #2])&" },
                { "JacobiSC", "(JacobiDC[#1, #2] JacobiNC[#1, #2])&" },
                { "JacobiSD", "(JacobiCD[#1, #2] JacobiND[#1, #2])&" },
                { "JacobiAmplitude", "JacobiDN[#1, #2]&" },
                { "InverseJacobiSN", "(JacobiCD[InverseJacobiSN[#1, #2], #2]/(1 - #1^2))&" },
                { "InverseJacobiCN", "(-JacobiDS[InverseJacobiCN[#1, #2], #2]/(1 - #2 + #2 #1^2))&" },
                { "InverseJacobiDN", "(-JacobiCS[InverseJacobiDN[#1, #2], #2]/(-1 + #2 + #1^2))&" },
            };
            for (size_t j = 0; j < sizeof(JD) / sizeof(JD[0]); j++) {
                if (strcmp(h, JD[j][0]) != 0) continue;
                Expr* U = args[0];
                Expr* M = args[1];
                Expr* dU = deriv_of(U, x, nonconsts);
                Expr* dM = deriv_of(M, x, nonconsts);
                Expr* terms[2];
                size_t nt = 0;
                if (!is_lit_zero(dU)) {
                    Expr* f = parse_expression(JD[j][1]);
                    Expr* g = expr_new_function(f, (Expr*[]){ expr_copy(U), expr_copy(M) }, 2);
                    terms[nt++] = mk_fn2("Times", g, dU);
                } else {
                    expr_free(dU);
                }
                if (!is_lit_zero(dM)) {
                    Expr* op = expr_new_function(mk_sym("Derivative"),
                                  (Expr*[]){ mk_int(0), mk_int(1) }, 2);
                    Expr* op_g = mk_fn_head1(op, mk_sym(JD[j][0]));
                    Expr* applied = expr_new_function(op_g,
                                  (Expr*[]){ expr_copy(U), expr_copy(M) }, 2);
                    terms[nt++] = mk_fn2("Times", applied, dM);
                } else {
                    expr_free(dM);
                }
                if (nt == 0) return mk_int(0);
                if (nt == 1) return terms[0];
                return mk_fn2("Plus", terms[0], terms[1]);
            }
        }
        /* --- GegenbauerC[N, L, Z] (DLMF 18.9.19):
         *   d/dZ C_N^L(Z) = 2 L C_{N-1}^{L+1}(Z);
         *   d/dN, d/dL are inert Derivative[1,0,0] / Derivative[0,1,0] forms.
         * (The two-argument form rewrites to ChebyshevT on evaluation.) */
        if (h == SYM_GegenbauerC && n == 3) {
            Expr* N = args[0];
            Expr* L = args[1];
            Expr* Z = args[2];
            Expr* dN = deriv_of(N, x, nonconsts);
            Expr* dL = deriv_of(L, x, nonconsts);
            Expr* dZ = deriv_of(Z, x, nonconsts);
            Expr* terms[3];
            size_t nt = 0;
            if (!is_lit_zero(dZ)) {
                Expr* c = expr_new_function(mk_sym("GegenbauerC"),
                    (Expr*[]){ mk_fn2("Plus", expr_copy(N), mk_int(-1)),
                               mk_fn2("Plus", expr_copy(L), mk_int(1)),
                               expr_copy(Z) }, 3);
                Expr* d = expr_new_function(mk_sym("Times"),
                    (Expr*[]){ mk_int(2), expr_copy(L), c }, 3);
                terms[nt++] = mk_fn2("Times", d, dZ);
            } else {
                expr_free(dZ);
            }
            Expr* dd[2] = { dN, dL };
            for (int s = 0; s < 2; s++) {
                if (is_lit_zero(dd[s])) { expr_free(dd[s]); continue; }
                Expr* op = expr_new_function(mk_sym("Derivative"),
                    (Expr*[]){ mk_int(s == 0), mk_int(s == 1), mk_int(0) }, 3);
                Expr* op_g = mk_fn_head1(op, mk_sym("GegenbauerC"));
                Expr* applied = expr_new_function(op_g,
                    (Expr*[]){ expr_copy(N), expr_copy(L), expr_copy(Z) }, 3);
                terms[nt++] = mk_fn2("Times", applied, dd[s]);
            }
            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return expr_new_function(mk_sym("Plus"), terms, nt);
        }

        if (h == SYM_HermiteH && n == 2) {
            Expr* N = args[0];
            Expr* Z = args[1];
            Expr* dN = deriv_of(N, x, nonconsts);
            Expr* dZ = deriv_of(Z, x, nonconsts);
            Expr* terms[2];
            size_t nt = 0;

            if (!is_lit_zero(dZ)) {
                /* 2 N HermiteH[N-1, Z] */
                Expr* hNm1 = mk_fn2("HermiteH",
                                mk_fn2("Plus", expr_copy(N), mk_int(-1)), expr_copy(Z));
                Expr* dHdZ = expr_new_function(mk_sym("Times"),
                                (Expr*[]){ mk_int(2), expr_copy(N), hNm1 }, 3);
                terms[nt++] = mk_fn2("Times", dHdZ, dZ);
            } else {
                expr_free(dZ);
            }

            if (!is_lit_zero(dN)) {
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(1), mk_int(0) }, 2);
                Expr* op_g = mk_fn_head1(op, mk_sym("HermiteH"));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(N), expr_copy(Z) }, 2);
                terms[nt++] = mk_fn2("Times", applied, dN);
            } else {
                expr_free(dN);
            }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return mk_fn2("Plus", terms[0], terms[1]);
        }

        /* --- LaguerreL[N, Z]:
         *   d/dZ L_N(Z) = -LaguerreL[N-1, 1, Z]  (the generalized Laguerre).
         *   d/dN L_N(Z) = Derivative[1,0][LaguerreL][N,Z]  (no elementary form).
         * Enables D[] and, via it, the generic Interval threading in eval.c. */
        if (h == SYM_LaguerreL && n == 2) {
            Expr* N = args[0];
            Expr* Z = args[1];
            Expr* dN = deriv_of(N, x, nonconsts);
            Expr* dZ = deriv_of(Z, x, nonconsts);
            Expr* terms[2];
            size_t nt = 0;

            if (!is_lit_zero(dZ)) {
                /* -LaguerreL[N-1, 1, Z] */
                Expr* lg = expr_new_function(mk_sym("LaguerreL"),
                              (Expr*[]){ mk_fn2("Plus", expr_copy(N), mk_int(-1)),
                                         mk_int(1), expr_copy(Z) }, 3);
                terms[nt++] = mk_fn2("Times", mk_neg(lg), dZ);
            } else {
                expr_free(dZ);
            }

            if (!is_lit_zero(dN)) {
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(1), mk_int(0) }, 2);
                Expr* op_g = mk_fn_head1(op, mk_sym("LaguerreL"));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(N), expr_copy(Z) }, 2);
                terms[nt++] = mk_fn2("Times", applied, dN);
            } else {
                expr_free(dN);
            }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return mk_fn2("Plus", terms[0], terms[1]);
        }

        /* --- LaguerreL[N, A, Z] (generalized):
         *   d/dZ L_N^A(Z) = -LaguerreL[N-1, A+1, Z].
         *   d/dN, d/dA have no elementary form (kept inert via Derivative). */
        if (h == SYM_LaguerreL && n == 3) {
            Expr* N = args[0];
            Expr* A = args[1];
            Expr* Z = args[2];
            Expr* dN = deriv_of(N, x, nonconsts);
            Expr* dA = deriv_of(A, x, nonconsts);
            Expr* dZ = deriv_of(Z, x, nonconsts);
            Expr* terms[3];
            size_t nt = 0;

            if (!is_lit_zero(dZ)) {
                /* -LaguerreL[N-1, A+1, Z] */
                Expr* lg = expr_new_function(mk_sym("LaguerreL"),
                              (Expr*[]){ mk_fn2("Plus", expr_copy(N), mk_int(-1)),
                                         mk_fn2("Plus", expr_copy(A), mk_int(1)),
                                         expr_copy(Z) }, 3);
                terms[nt++] = mk_fn2("Times", mk_neg(lg), dZ);
            } else {
                expr_free(dZ);
            }

            if (!is_lit_zero(dN)) {
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(1), mk_int(0), mk_int(0) }, 3);
                Expr* op_g = mk_fn_head1(op, mk_sym("LaguerreL"));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(N), expr_copy(A), expr_copy(Z) }, 3);
                terms[nt++] = mk_fn2("Times", applied, dN);
            } else {
                expr_free(dN);
            }

            if (!is_lit_zero(dA)) {
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(0), mk_int(1), mk_int(0) }, 3);
                Expr* op_g = mk_fn_head1(op, mk_sym("LaguerreL"));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(N), expr_copy(A), expr_copy(Z) }, 3);
                terms[nt++] = mk_fn2("Times", applied, dA);
            } else {
                expr_free(dA);
            }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return expr_new_function(mk_sym("Plus"), terms, nt);
        }

        /* --- JacobiP[N, A, B, Z]:
         *   d/dZ P_N^(A,B)(Z) = (N+A+B+1)/2 JacobiP[N-1, A+1, B+1, Z]  (DLMF 18.9.15).
         *   d/dN, d/dA, d/dB have no elementary form (kept inert via Derivative). */
        if (h == SYM_JacobiP && n == 4) {
            Expr* terms[4];
            size_t nt = 0;
            for (int slot = 0; slot < 4; slot++) {
                Expr* ds = deriv_of(args[slot], x, nonconsts);
                if (is_lit_zero(ds)) { expr_free(ds); continue; }
                Expr* piece;
                if (slot == 3) {
                    Expr* jp = expr_new_function(mk_sym("JacobiP"), (Expr*[]){
                        mk_fn2("Plus", expr_copy(args[0]), mk_int(-1)),
                        mk_fn2("Plus", expr_copy(args[1]), mk_int(1)),
                        mk_fn2("Plus", expr_copy(args[2]), mk_int(1)),
                        expr_copy(args[3]) }, 4);
                    Expr* c = expr_new_function(mk_sym("Plus"), (Expr*[]){
                        expr_copy(args[0]), expr_copy(args[1]), expr_copy(args[2]),
                        mk_int(1) }, 4);
                    piece = expr_new_function(mk_sym("Times"), (Expr*[]){
                        mk_fn2("Rational", mk_int(1), mk_int(2)), c, jp }, 3);
                } else {
                    Expr* op = expr_new_function(mk_sym("Derivative"), (Expr*[]){
                        mk_int(slot == 0), mk_int(slot == 1), mk_int(slot == 2), mk_int(0) }, 4);
                    Expr* op_g = mk_fn_head1(op, mk_sym("JacobiP"));
                    piece = expr_new_function(op_g, (Expr*[]){
                        expr_copy(args[0]), expr_copy(args[1]), expr_copy(args[2]),
                        expr_copy(args[3]) }, 4);
                }
                terms[nt++] = mk_fn2("Times", piece, ds);
            }
            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return expr_new_function(mk_sym("Plus"), terms, nt);
        }

        /* --- ExpIntegralE[N, Z]: d/dZ E_N(Z) = -E_(N-1)(Z) (DLMF 8.19.13);
         *     d/dN has no elementary form (kept inert via Derivative[1,0]). */
        if (h == SYM_ExpIntegralE && n == 2) {
            Expr* N = args[0];
            Expr* Z = args[1];
            Expr* dN = deriv_of(N, x, nonconsts);
            Expr* dZ = deriv_of(Z, x, nonconsts);
            Expr* terms[2];
            size_t nt = 0;

            if (!is_lit_zero(dZ)) {
                /* -ExpIntegralE[N-1, Z] * dZ */
                Expr* em1 = mk_fn2("ExpIntegralE",
                              mk_fn2("Plus", expr_copy(N), mk_int(-1)), expr_copy(Z));
                terms[nt++] = mk_fn2("Times", mk_fn2("Times", mk_int(-1), em1), dZ);
            } else {
                expr_free(dZ);
            }

            if (!is_lit_zero(dN)) {
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(1), mk_int(0) }, 2);
                Expr* op_g = mk_fn_head1(op, mk_sym("ExpIntegralE"));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(N), expr_copy(Z) }, 2);
                terms[nt++] = mk_fn2("Times", applied, dN);
            } else {
                expr_free(dN);
            }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return mk_fn2("Plus", terms[0], terms[1]);
        }

        /* --- BesselK[N, Z]: chain rule on both args.
         *   d/dZ BesselK[N, Z] = -(BesselK[N-1, Z] + BesselK[N+1, Z]) / 2  (DLMF 10.29.5)
         *   d/dN BesselK[N, Z] = Derivative[1,0][BesselK][N,Z]  (no elementary form).
         * Note the sign differs from BesselJ. */
        if (h == SYM_BesselK && n == 2) {
            Expr* N = args[0];
            Expr* Z = args[1];
            Expr* dN = deriv_of(N, x, nonconsts);
            Expr* dZ = deriv_of(Z, x, nonconsts);
            Expr* terms[2];
            size_t nt = 0;

            if (!is_lit_zero(dZ)) {
                /* -(BesselK[N-1, Z] + BesselK[N+1, Z]) / 2 * dZ */
                Expr* bm1 = mk_fn2("BesselK",
                              mk_fn2("Plus", expr_copy(N), mk_int(-1)), expr_copy(Z));
                Expr* bp1 = mk_fn2("BesselK",
                              mk_fn2("Plus", expr_copy(N), mk_int(1)), expr_copy(Z));
                Expr* sum = mk_fn2("Plus", bm1, bp1);
                Expr* dBdZ = mk_fn2("Times",
                              mk_fn2("Power", mk_int(-2), mk_int(-1)), sum);
                terms[nt++] = mk_fn2("Times", dBdZ, dZ);
            } else {
                expr_free(dZ);
            }

            if (!is_lit_zero(dN)) {
                /* Derivative[1,0][BesselK][N,Z] * dN */
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(1), mk_int(0) }, 2);
                Expr* op_g = mk_fn_head1(op, mk_sym("BesselK"));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(N), expr_copy(Z) }, 2);
                terms[nt++] = mk_fn2("Times", applied, dN);
            } else {
                expr_free(dN);
            }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return mk_fn2("Plus", terms[0], terms[1]);
        }

        /* --- BesselI[N, Z]: chain rule on both args.
         *   d/dZ BesselI[N, Z] = (BesselI[N-1, Z] + BesselI[N+1, Z]) / 2  (DLMF 10.29.5)
         *   d/dN BesselI[N, Z] = Derivative[1,0][BesselI][N,Z]  (no elementary form).
         * Like BesselK the two-term sum carries a '+', but (like BesselJ) the
         * overall coefficient is +1/2 rather than -1/2. */
        if (h == SYM_BesselI && n == 2) {
            Expr* N = args[0];
            Expr* Z = args[1];
            Expr* dN = deriv_of(N, x, nonconsts);
            Expr* dZ = deriv_of(Z, x, nonconsts);
            Expr* terms[2];
            size_t nt = 0;

            if (!is_lit_zero(dZ)) {
                /* (BesselI[N-1, Z] + BesselI[N+1, Z]) / 2 * dZ */
                Expr* bm1 = mk_fn2("BesselI",
                              mk_fn2("Plus", expr_copy(N), mk_int(-1)), expr_copy(Z));
                Expr* bp1 = mk_fn2("BesselI",
                              mk_fn2("Plus", expr_copy(N), mk_int(1)), expr_copy(Z));
                Expr* sum = mk_fn2("Plus", bm1, bp1);
                Expr* dBdZ = mk_fn2("Times",
                              mk_fn2("Power", mk_int(2), mk_int(-1)), sum);
                terms[nt++] = mk_fn2("Times", dBdZ, dZ);
            } else {
                expr_free(dZ);
            }

            if (!is_lit_zero(dN)) {
                /* Derivative[1,0][BesselI][N,Z] * dN */
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(1), mk_int(0) }, 2);
                Expr* op_g = mk_fn_head1(op, mk_sym("BesselI"));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(N), expr_copy(Z) }, 2);
                terms[nt++] = mk_fn2("Times", applied, dN);
            } else {
                expr_free(dN);
            }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return mk_fn2("Plus", terms[0], terms[1]);
        }

        /* --- BesselY[N, Z]: chain rule on both args.
         *   d/dZ BesselY[N, Z] = (BesselY[N-1, Z] - BesselY[N+1, Z]) / 2  (DLMF 10.6.1)
         *   d/dN BesselY[N, Z] = Derivative[1,0][BesselY][N,Z]  (no elementary form).
         * Identical in shape to BesselJ (Y satisfies the same recurrence). */
        if (h == SYM_BesselY && n == 2) {
            Expr* N = args[0];
            Expr* Z = args[1];
            Expr* dN = deriv_of(N, x, nonconsts);
            Expr* dZ = deriv_of(Z, x, nonconsts);
            Expr* terms[2];
            size_t nt = 0;

            if (!is_lit_zero(dZ)) {
                /* (BesselY[N-1, Z] - BesselY[N+1, Z]) / 2 * dZ */
                Expr* bm1 = mk_fn2("BesselY",
                              mk_fn2("Plus", expr_copy(N), mk_int(-1)), expr_copy(Z));
                Expr* bp1 = mk_fn2("BesselY",
                              mk_fn2("Plus", expr_copy(N), mk_int(1)), expr_copy(Z));
                Expr* diff = mk_fn2("Plus", bm1, mk_neg(bp1));
                Expr* dBdZ = mk_fn2("Times",
                              mk_fn2("Power", mk_int(2), mk_int(-1)), diff);
                terms[nt++] = mk_fn2("Times", dBdZ, dZ);
            } else {
                expr_free(dZ);
            }

            if (!is_lit_zero(dN)) {
                /* Derivative[1,0][BesselY][N,Z] * dN */
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(1), mk_int(0) }, 2);
                Expr* op_g = mk_fn_head1(op, mk_sym("BesselY"));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(N), expr_copy(Z) }, 2);
                terms[nt++] = mk_fn2("Times", applied, dN);
            } else {
                expr_free(dN);
            }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return mk_fn2("Plus", terms[0], terms[1]);
        }

        /* --- Gamma[A, Z] (upper incomplete gamma): chain rule on both args.
         *   dGamma/dZ = -Z^(A-1) E^-Z    (elementary)
         *   dGamma/dA = Derivative[1,0][Gamma][A,Z]  (no closed form without
         *               PolyGamma; left as the generic partial).
         * Zero-derivative arms are dropped, so D[Gamma[a,z],z] keeps only the
         * dZ term, matching Mathematica's -E^-z z^(a-1). */
        if (h == SYM_Gamma && n == 2) {
            Expr* A = args[0];
            Expr* Z = args[1];
            Expr* dA = deriv_of(A, x, nonconsts);
            Expr* dZ = deriv_of(Z, x, nonconsts);
            Expr* terms[2];
            size_t nt = 0;

            if (!is_lit_zero(dZ)) {
                /* -Z^(A-1) E^-Z * dZ */
                Expr* zpow = mk_fn2("Power", expr_copy(Z),
                                    mk_fn2("Plus", expr_copy(A), mk_int(-1)));
                Expr* enz  = mk_fn1("Exp", mk_neg(expr_copy(Z)));
                Expr* dGdZ = mk_neg(mk_fn2("Times", zpow, enz));
                terms[nt++] = mk_fn2("Times", dGdZ, dZ);
            } else {
                expr_free(dZ);
            }

            if (!is_lit_zero(dA)) {
                /* Derivative[1,0][Gamma][A,Z] * dA */
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(1), mk_int(0) }, 2);
                Expr* op_g = mk_fn_head1(op, mk_sym("Gamma"));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(A), expr_copy(Z) }, 2);
                terms[nt++] = mk_fn2("Times", applied, dA);
            } else {
                expr_free(dA);
            }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return mk_fn2("Plus", terms[0], terms[1]);
        }

        /* --- PolyGamma[N, Z] (the polygamma family): chain rule on both args.
         *   d/dZ PolyGamma[N, Z] = PolyGamma[N+1, Z]   (raises the order)
         *   d/dN PolyGamma[N, Z] = Derivative[1,0][PolyGamma][N,Z]  (generic; the
         *               derivative with respect to the order has no elementary form).
         * Zero-derivative arms are dropped, so D[PolyGamma[n,x],x] keeps only the
         * dZ term, matching Mathematica's PolyGamma[1+n, x]. */
        if (h == SYM_PolyGamma && n == 2) {
            Expr* N = args[0];
            Expr* Z = args[1];
            Expr* dN = deriv_of(N, x, nonconsts);
            Expr* dZ = deriv_of(Z, x, nonconsts);
            Expr* terms[2];
            size_t nt = 0;

            if (!is_lit_zero(dZ)) {
                /* PolyGamma[N+1, Z] * dZ */
                Expr* np1 = mk_fn2("Plus", expr_copy(N), mk_int(1));
                Expr* pg  = mk_fn2("PolyGamma", np1, expr_copy(Z));
                terms[nt++] = mk_fn2("Times", pg, dZ);
            } else {
                expr_free(dZ);
            }

            if (!is_lit_zero(dN)) {
                /* Derivative[1, 0][PolyGamma][N, Z] * dN */
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(1), mk_int(0) }, 2);
                Expr* op_g = mk_fn_head1(op, mk_sym("PolyGamma"));
                Expr* applied = expr_new_function(op_g,
                              (Expr*[]){ expr_copy(N), expr_copy(Z) }, 2);
                terms[nt++] = mk_fn2("Times", applied, dN);
            } else {
                expr_free(dN);
            }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return mk_fn2("Plus", terms[0], terms[1]);
        }

        /* --- Beta[A, B] (Euler beta): chain rule on both arguments.
         *   d/dA Beta[A,B] = Beta[A,B] (PolyGamma[0,A] - PolyGamma[0,A+B])
         *   d/dB Beta[A,B] = Beta[A,B] (PolyGamma[0,B] - PolyGamma[0,A+B])
         * Higher derivatives compose automatically by re-differentiating the
         * Beta and PolyGamma factors. Zero-derivative arms are dropped. */
        if (h == SYM_Beta && n == 2) {
            Expr* A = args[0];
            Expr* B = args[1];
            Expr* dA = deriv_of(A, x, nonconsts);
            Expr* dB = deriv_of(B, x, nonconsts);
            Expr* terms[2];
            size_t nt = 0;

            if (!is_lit_zero(dA)) {
                Expr* psiA = mk_fn2("PolyGamma", mk_int(0), expr_copy(A));
                Expr* psiS = mk_fn2("PolyGamma", mk_int(0),
                                    mk_fn2("Plus", expr_copy(A), expr_copy(B)));
                Expr* fac  = mk_fn2("Times",
                                    mk_fn2("Beta", expr_copy(A), expr_copy(B)),
                                    mk_fn2("Subtract", psiA, psiS));
                terms[nt++] = mk_fn2("Times", fac, dA);
            } else { expr_free(dA); }

            if (!is_lit_zero(dB)) {
                Expr* psiB = mk_fn2("PolyGamma", mk_int(0), expr_copy(B));
                Expr* psiS = mk_fn2("PolyGamma", mk_int(0),
                                    mk_fn2("Plus", expr_copy(A), expr_copy(B)));
                Expr* fac  = mk_fn2("Times",
                                    mk_fn2("Beta", expr_copy(A), expr_copy(B)),
                                    mk_fn2("Subtract", psiB, psiS));
                terms[nt++] = mk_fn2("Times", fac, dB);
            } else { expr_free(dB); }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return mk_fn2("Plus", terms[0], terms[1]);
        }

        /* --- Beta[Z, A, B] (incomplete beta): chain rule on all three args.
         *   d/dZ Beta[Z,A,B] = Z^(A-1) (1-Z)^(B-1)   (the integrand; elementary)
         *   d/dA, d/dB        = generic Derivative[0,1,0]/[0,0,1][Beta][Z,A,B]
         *               (no elementary form). Zero-derivative arms are dropped,
         * so D[Beta[z,a,b],z] keeps only Z^(a-1)(1-z)^(b-1). */
        if (h == SYM_Beta && n == 3) {
            Expr* Z = args[0];
            Expr* A = args[1];
            Expr* B = args[2];
            Expr* dZ = deriv_of(Z, x, nonconsts);
            Expr* dA = deriv_of(A, x, nonconsts);
            Expr* dB = deriv_of(B, x, nonconsts);
            Expr* terms[3];
            size_t nt = 0;

            if (!is_lit_zero(dZ)) {
                Expr* za  = mk_fn2("Power", expr_copy(Z),
                                   mk_fn2("Plus", expr_copy(A), mk_int(-1)));
                Expr* omz = mk_fn2("Power",
                                   mk_fn2("Subtract", mk_int(1), expr_copy(Z)),
                                   mk_fn2("Plus", expr_copy(B), mk_int(-1)));
                terms[nt++] = mk_fn2("Times", mk_fn2("Times", za, omz), dZ);
            } else { expr_free(dZ); }

            if (!is_lit_zero(dA)) {
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(0), mk_int(1), mk_int(0) }, 3);
                Expr* opb = mk_fn_head1(op, mk_sym("Beta"));
                Expr* applied = expr_new_function(opb,
                              (Expr*[]){ expr_copy(Z), expr_copy(A), expr_copy(B) }, 3);
                terms[nt++] = mk_fn2("Times", applied, dA);
            } else { expr_free(dA); }

            if (!is_lit_zero(dB)) {
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(0), mk_int(0), mk_int(1) }, 3);
                Expr* opb = mk_fn_head1(op, mk_sym("Beta"));
                Expr* applied = expr_new_function(opb,
                              (Expr*[]){ expr_copy(Z), expr_copy(A), expr_copy(B) }, 3);
                terms[nt++] = mk_fn2("Times", applied, dB);
            } else { expr_free(dB); }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return expr_new_function(mk_sym("Plus"), terms, nt);
        }

        /* --- Beta[Z0, Z1, A, B] (generalized incomplete beta).
         *   d/dZ1 = Z1^(A-1) (1-Z1)^(B-1),  d/dZ0 = -Z0^(A-1) (1-Z0)^(B-1)
         *   d/dA, d/dB = generic Derivative[0,0,1,0]/[0,0,0,1][Beta][...]. */
        if (h == SYM_Beta && n == 4) {
            Expr* Z0 = args[0];
            Expr* Z1 = args[1];
            Expr* A  = args[2];
            Expr* B  = args[3];
            Expr* dZ0 = deriv_of(Z0, x, nonconsts);
            Expr* dZ1 = deriv_of(Z1, x, nonconsts);
            Expr* dA  = deriv_of(A,  x, nonconsts);
            Expr* dB  = deriv_of(B,  x, nonconsts);
            Expr* terms[4];
            size_t nt = 0;

            if (!is_lit_zero(dZ1)) {
                Expr* za  = mk_fn2("Power", expr_copy(Z1),
                                   mk_fn2("Plus", expr_copy(A), mk_int(-1)));
                Expr* omz = mk_fn2("Power",
                                   mk_fn2("Subtract", mk_int(1), expr_copy(Z1)),
                                   mk_fn2("Plus", expr_copy(B), mk_int(-1)));
                terms[nt++] = mk_fn2("Times", mk_fn2("Times", za, omz), dZ1);
            } else { expr_free(dZ1); }

            if (!is_lit_zero(dZ0)) {
                Expr* za  = mk_fn2("Power", expr_copy(Z0),
                                   mk_fn2("Plus", expr_copy(A), mk_int(-1)));
                Expr* omz = mk_fn2("Power",
                                   mk_fn2("Subtract", mk_int(1), expr_copy(Z0)),
                                   mk_fn2("Plus", expr_copy(B), mk_int(-1)));
                terms[nt++] = mk_fn2("Times", mk_neg(mk_fn2("Times", za, omz)), dZ0);
            } else { expr_free(dZ0); }

            if (!is_lit_zero(dA)) {
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(0), mk_int(0), mk_int(1), mk_int(0) }, 4);
                Expr* opb = mk_fn_head1(op, mk_sym("Beta"));
                Expr* applied = expr_new_function(opb,
                              (Expr*[]){ expr_copy(Z0), expr_copy(Z1),
                                         expr_copy(A), expr_copy(B) }, 4);
                terms[nt++] = mk_fn2("Times", applied, dA);
            } else { expr_free(dA); }

            if (!is_lit_zero(dB)) {
                Expr* op = expr_new_function(mk_sym("Derivative"),
                              (Expr*[]){ mk_int(0), mk_int(0), mk_int(0), mk_int(1) }, 4);
                Expr* opb = mk_fn_head1(op, mk_sym("Beta"));
                Expr* applied = expr_new_function(opb,
                              (Expr*[]){ expr_copy(Z0), expr_copy(Z1),
                                         expr_copy(A), expr_copy(B) }, 4);
                terms[nt++] = mk_fn2("Times", applied, dB);
            } else { expr_free(dB); }

            if (nt == 0) return mk_int(0);
            if (nt == 1) return terms[0];
            return expr_new_function(mk_sym("Plus"), terms, nt);
        }

        /* --- Gamma[A] (one-argument): d/dA Gamma[A] = Gamma[A] PolyGamma[0, A].
         * Repeated differentiation then composes via the product rule, e.g.
         * D[Gamma[z], {z, 2}] = Gamma[z] PolyGamma[0,z]^2 + Gamma[z] PolyGamma[1,z]. */
        if (h == SYM_Gamma && n == 1) {
            Expr* A = args[0];
            Expr* dA = deriv_of(A, x, nonconsts);
            Expr* g  = mk_fn1("Gamma", expr_copy(A));
            Expr* pg = mk_fn2("PolyGamma", mk_int(0), expr_copy(A));
            Expr* dGdA = mk_fn2("Times", g, pg);
            return mk_fn2("Times", dGdA, dA);
        }

        /* --- LogGamma[A]: d/dA LogGamma[A] = PolyGamma[0, A] (the digamma).
         * Higher derivatives then raise the PolyGamma order automatically. */
        if (h == SYM_LogGamma && n == 1) {
            Expr* A = args[0];
            Expr* dA = deriv_of(A, x, nonconsts);
            Expr* pg = mk_fn2("PolyGamma", mk_int(0), expr_copy(A));
            return mk_fn2("Times", pg, dA);
        }

        /* --- Piecewise[{{v1, c1}, ...}, default]: differentiate each value
         *     expression and the default; the conditions ride through
         *     unchanged. D[Piecewise[{{vi, ci}}, d], x]
         *       = Piecewise[{{D[vi, x], ci}}, D[d, x]].
         *     Without this rule a second derivative of UnitStep (whose first
         *     derivative is a Piecewise) falls through to the generic chain
         *     rule and produces garbage. The default is 0 when omitted. --- */
        if (h == SYM_Piecewise && n >= 1 &&
            args[0]->type == EXPR_FUNCTION &&
            args[0]->data.function.head->type == EXPR_SYMBOL &&
            args[0]->data.function.head->data.symbol.name == SYM_List) {
            Expr* pairs = args[0];
            size_t np = pairs->data.function.arg_count;
            Expr** new_pairs = malloc(sizeof(Expr*) * (np ? np : 1));
            bool ok = true;
            size_t built = 0;
            for (size_t i = 0; i < np; i++) {
                Expr* pr = pairs->data.function.args[i];
                if (!(pr->type == EXPR_FUNCTION &&
                      pr->data.function.head->type == EXPR_SYMBOL &&
                      pr->data.function.head->data.symbol.name == SYM_List &&
                      pr->data.function.arg_count == 2)) { ok = false; break; }
                /* Piecewise has HoldAll, so its branch values are not
                 * re-evaluated by the outer evaluator; reduce each derivative
                 * here (e.g. 2 x^(2-1) -> 2 x). Conditions ride through
                 * verbatim. */
                Expr* dv   = eval_and_free(deriv_of(pr->data.function.args[0], x, nonconsts));
                Expr* cond = expr_copy(pr->data.function.args[1]);
                new_pairs[built++] = mk_fn2("List", dv, cond);
            }
            if (ok) {
                Expr* pairs_list = mk_fnN_adopt("List", new_pairs, np);
                Expr* ddefault = (n >= 2) ? eval_and_free(deriv_of(args[1], x, nonconsts))
                                          : mk_int(0);
                return mk_fn2("Piecewise", pairs_list, ddefault);
            }
            for (size_t i = 0; i < built; i++) expr_free(new_pairs[i]);
            free(new_pairs);
        }

        /* --- UnitStep[a1, ..., an]: acts as the product of the per-argument
         *     UnitStep[ai], so the derivative follows the product rule. Each
         *     factor contributes UnitStep'[ai] = Piecewise[{{Indeterminate,
         *     ai == 0}}, 0], multiplied by D[ai, x] and the UnitStep of the
         *     remaining arguments. For n == 1 the trailing UnitStep factor is
         *     absent and the result collapses to the single Piecewise. --- */
        if (h == SYM_UnitStep && n >= 1) {
            Expr** terms = malloc(sizeof(Expr*) * n);
            for (size_t i = 0; i < n; i++) {
                Expr* cond  = mk_fn2("Equal", expr_copy(args[i]), mk_int(0));
                Expr* pair  = mk_fn2("List", mk_sym("Indeterminate"), cond);
                Expr* pairs = mk_fn1("List", pair);
                Expr* pw    = mk_fn2("Piecewise", pairs, mk_int(0));
                Expr* dai   = deriv_of(args[i], x, nonconsts);

                size_t fc = 0;
                Expr** factors = malloc(sizeof(Expr*) * 3);
                factors[fc++] = pw;
                factors[fc++] = dai;
                if (n > 1) {
                    Expr** rest = malloc(sizeof(Expr*) * (n - 1));
                    size_t rc = 0;
                    for (size_t j = 0; j < n; j++)
                        if (j != i) rest[rc++] = expr_copy(args[j]);
                    factors[fc++] = mk_fnN_adopt("UnitStep", rest, n - 1);
                }
                terms[i] = mk_fnN_adopt("Times", factors, fc);
            }
            return mk_fnN_adopt("Plus", terms, n);
        }

        /* --- Floor / Ceiling / Round / IntegerPart / FractionalPart: the
         *     integer-rounding functions are piecewise-constant, so their
         *     derivative is a Piecewise that is 0 (or 1, for FractionalPart)
         *     off the jump set and Indeterminate on it, multiplied by D[g, x]
         *     (the chain rule, via deriv_of). Matches Mathematica; modeled on
         *     the UnitStep handler above (Times[pw, 1] collapses to pw for the
         *     bare g == x case). --- */
        if (h == SYM_Floor && n == 1) {
            Expr* cond = mk_fn2("Greater", expr_copy(args[0]),
                                mk_fn1("Floor", expr_copy(args[0])));
            Expr* pw = mk_fn2("Piecewise",
                              mk_fn1("List", mk_fn2("List", mk_int(0), cond)),
                              mk_sym("Indeterminate"));
            return mk_fn2("Times", pw, deriv_of(args[0], x, nonconsts));
        }
        if (h == SYM_Ceiling && n == 1) {
            Expr* cond = mk_fn2("Less", expr_copy(args[0]),
                                mk_fn1("Ceiling", expr_copy(args[0])));
            Expr* pw = mk_fn2("Piecewise",
                              mk_fn1("List", mk_fn2("List", mk_int(0), cond)),
                              mk_sym("Indeterminate"));
            return mk_fn2("Times", pw, deriv_of(args[0], x, nonconsts));
        }
        if (h == SYM_Round && n == 1) {
            /* zero off the half-integer jump set: both Re[g]-1/2 and Im[g]-1/2
             * non-integers. */
            Expr* reCond = mk_fn2("NotElement",
                mk_fn2("Plus", mk_fn2("Rational", mk_int(-1), mk_int(2)),
                       mk_fn1("Re", expr_copy(args[0]))),
                mk_sym("Integers"));
            Expr* imCond = mk_fn2("NotElement",
                mk_fn2("Plus", mk_fn2("Rational", mk_int(-1), mk_int(2)),
                       mk_fn1("Im", expr_copy(args[0]))),
                mk_sym("Integers"));
            Expr* pw = mk_fn2("Piecewise",
                mk_fn1("List", mk_fn2("List", mk_int(0),
                                      mk_fn2("And", reCond, imCond))),
                mk_sym("Indeterminate"));
            return mk_fn2("Times", pw, deriv_of(args[0], x, nonconsts));
        }
        if ((h == SYM_IntegerPart || h == SYM_FractionalPart) && n == 1) {
            /* IntegerPart' = 0, FractionalPart' = 1, off the integer jump set:
             * each of Re[g], Im[g] either zero or a non-integer. */
            Expr* reCond = mk_fn2("Or",
                mk_fn2("Equal", mk_fn1("Re", expr_copy(args[0])), mk_int(0)),
                mk_fn2("NotElement", mk_fn1("Re", expr_copy(args[0])),
                       mk_sym("Integers")));
            Expr* imCond = mk_fn2("Or",
                mk_fn2("Equal", mk_fn1("Im", expr_copy(args[0])), mk_int(0)),
                mk_fn2("NotElement", mk_fn1("Im", expr_copy(args[0])),
                       mk_sym("Integers")));
            Expr* val = mk_int(h == SYM_FractionalPart ? 1 : 0);
            Expr* pw = mk_fn2("Piecewise",
                mk_fn1("List", mk_fn2("List", val,
                                      mk_fn2("And", reCond, imCond))),
                mk_sym("Indeterminate"));
            return mk_fn2("Times", pw, deriv_of(args[0], x, nonconsts));
        }

        /* --- SeriesData: term-by-term differentiation. --- */
        if (h == SYM_SeriesData && n == 6) {
            Expr* r = series_differentiate(f, x);
            if (r) return r;
        }

        /* A held expression is not a function of its contents: D[Hold[x], x]
         * stays inert, as in the reference (and Sin[Hold[x]] gives
         * Cos[Hold[x]] D[Hold[x], x]). */
        if (h == SYM_Hold || h == SYM_HoldForm || h == SYM_HoldComplete ||
            h == SYM_HoldPattern)
            return NULL;

        /* --- Known elementary unary function: F'(g) * D[g, x]. --- */
        if (n == 1) {
            Expr* fp = elementary_fprime(h, args[0]);
            if (fp) {
                Expr* dg = deriv_of(args[0], x, nonconsts);
                return mk_fn2("Times", fp, dg);
            }
            /* Unknown single-argument function -- standard chain rule:
             *     Derivative[1][f][g] * D[g, x]. */
            Expr* op = mk_fn1("Derivative", mk_int(1));         /* Derivative[1]      */
            Expr* op_f = mk_fn_head1(op, expr_copy(head));      /* Derivative[1][f]   */
            Expr* applied = mk_fn_head1(op_f, expr_copy(args[0])); /* Derivative[1][f][g] */
            Expr* dg = deriv_of(args[0], x, nonconsts);
            return mk_fn2("Times", applied, dg);
        }

        /* --- Unknown multi-argument function: full chain rule. --- */
        return chain_rule_unknown(f, x, nonconsts);
    }

    /* ------------------------------------------------------------------ */
    /* Head is itself a function. The only pattern we can reduce is       */
    /* Derivative[...][f][args...].                                       */
    /* ------------------------------------------------------------------ */
    if (head->type == EXPR_FUNCTION) {
        /* Fundamental theorem of calculus for an INACTIVE integral:
         *   D[Inactive[Integrate][f, u], u] = f   (indefinite, u == diff var).
         * This lets an inert integral (built to dodge Integrate's cost) verify by
         * the implicit-function rule without ever running the integration cascade.
         * A different differentiation variable is already handled by the free-of-x
         * short-circuit above (an x-free inactive integral -> 0). */
        Expr* ih = head->data.function.head;
        bool inact_int = (head->data.function.arg_count == 1
            && ih->type == EXPR_SYMBOL && ih->data.symbol.name == SYM_Inactive
            && head->data.function.args[0]->type == EXPR_SYMBOL
            && head->data.function.args[0]->data.symbol.name == SYM_Integrate);
        if (x && inact_int && n == 2 && expr_eq(args[1], x))
            return expr_copy(args[0]);
        /* An inactive integral over a non-list "variable" v (Sqrt[x], or a
         * symbol other than x): the reference's Leibniz form
         *   D[Inactive[Integrate][e, v], x] = e D[v, x] + Inactive[Integrate][D[e, x], v],
         * so D[Inactive[Integrate][f[u], Sqrt[x]], x] is f[u]/(2 Sqrt[x]). */
        bool v_is_list = n == 2 && args[1]->type == EXPR_FUNCTION &&
                         is_sym(args[1]->data.function.head, "List");
        if (x && inact_int && n == 2 && !v_is_list) {
            Expr* dv = deriv_of(args[1], x, nonconsts);
            Expr* de = deriv_of(args[0], x, nonconsts);
            if (!dv || !de) { if (dv) expr_free(dv); if (de) expr_free(de); return NULL; }
            dv = eval_and_free(dv);
            de = eval_and_free(de);
            Expr** terms = malloc(2 * sizeof(Expr*));
            size_t nt = 0;
            if (!(dv->type == EXPR_INTEGER && dv->data.integer == 0))
                terms[nt++] = mk_fn2("Times", expr_copy(args[0]), dv);
            else expr_free(dv);
            if (!(de->type == EXPR_INTEGER && de->data.integer == 0)) {
                Expr* inact = expr_new_function(expr_new_symbol(SYM_Inactive),
                                  (Expr*[]){ expr_new_symbol(SYM_Integrate) }, 1);
                terms[nt++] = expr_new_function(inact, (Expr*[]){ de, expr_copy(args[1]) }, 2);
            } else expr_free(de);
            return nt == 0 ? (free(terms), expr_new_integer(0)) : mk_fnN_adopt(SYM_Plus, terms, nt);
        }
        /* The DEFINITE companion — the same Leibniz rule the ACTIVE head gets above,
         * for Inactive[Integrate][e, {u, a, b}].  It is what lets a first integral
         * whose constant has been fitted at a point be written in the honest definite
         * form Integrate[1/p, {t, y0, y[x]}] == x - x0 and still verify through the
         * implicit-function rule; without it the fitted relation had to be built by
         * substituting y[x] -> y0 across the whole expression, which also rewrote the
         * integration-variable slot and produced a corrupt answer.  The inner
         * Integrate[D[e,x], ...] term stays INERT, so no integration is attempted. */
        if (x && inact_int && n == 2
            && args[1]->type == EXPR_FUNCTION
            && args[1]->data.function.head->type == EXPR_SYMBOL
            && args[1]->data.function.head->data.symbol.name == SYM_List
            && args[1]->data.function.arg_count == 3
            && args[1]->data.function.args[0]->type == EXPR_SYMBOL
            && !expr_eq(args[1]->data.function.args[0], x)) {
            Expr* e  = args[0];
            Expr* u  = args[1]->data.function.args[0];
            Expr* lo = args[1]->data.function.args[1];
            Expr* hi = args[1]->data.function.args[2];
            Expr* e_hi = eval_and_free(mk_fn2("ReplaceAll", expr_copy(e),
                             mk_fn2("Rule", expr_copy(u), expr_copy(hi))));
            Expr* e_lo = eval_and_free(mk_fn2("ReplaceAll", expr_copy(e),
                             mk_fn2("Rule", expr_copy(u), expr_copy(lo))));
            Expr* dhi = eval_and_free(deriv_of(hi, x, nonconsts));
            Expr* dlo = eval_and_free(deriv_of(lo, x, nonconsts));
            Expr* de  = eval_and_free(deriv_of(e, x, nonconsts));
            Expr** listitems = malloc(3 * sizeof(Expr*));
            listitems[0] = expr_copy(u);
            listitems[1] = expr_copy(lo);
            listitems[2] = expr_copy(hi);
            Expr* newspec = mk_fnN_adopt(SYM_List, listitems, 3);
            /* An x-free integrand contributes nothing; dropping the term keeps the
             * common case exactly f[b] b' instead of trailing an inert
             * Inactive[Integrate][0, ...] that no simplifier will remove. */
            bool de_zero = (de->type == EXPR_INTEGER && de->data.integer == 0);
            Expr** terms = malloc(3 * sizeof(Expr*));
            size_t nt = 0;
            terms[nt++] = mk_fn2("Times", e_hi, dhi);
            terms[nt++] = mk_neg(mk_fn2("Times", e_lo, dlo));
            if (de_zero) { expr_free(de); expr_free(newspec); }
            else {
                Expr* inact = expr_new_function(expr_new_symbol(SYM_Inactive),
                                  (Expr*[]){ expr_new_symbol(SYM_Integrate) }, 1);
                terms[nt++] = expr_new_function(inact, (Expr*[]){ de, newspec }, 2);
            }
            return mk_fnN_adopt(SYM_Plus, terms, nt);
        }

        /* Any other inactive integral (a list spec the rules above do not
         * cover, e.g. {u, Sqrt[x]}) stays an inert D[...], as in the reference;
         * it is not an unknown function of its arguments. */
        if (inact_int) return NULL;

        Expr* r = deriv_of_derivative_form(f, x, nonconsts);
        if (r) return r;

        /* Applied InterpolatingFunction object: differentiate by the generic
         * chain rule, producing Derivative[0..1_k..0][InterpolatingFunction[...]]
         * [g1..gn] * D[gk, x]. The evaluator then reduces the Derivative head
         * to a derivative-annotated InterpolatingFunction (see interp.c), so
         * D[ifun[x], x] and ifun'[x] agree. */
        if (is_sym(head->data.function.head, "InterpolatingFunction")) {
            return chain_rule_unknown(f, x, nonconsts);
        }

        /* Any other compound head h[a][g1..gn] is an unknown function h[a] of
         * its arguments, as in the reference: D[g[a][x], x] is
         * Derivative[1][g[a]][x], by the same generic chain rule. A head that
         * itself depends on x adds the derivative through the head, applied to
         * the arguments: D[g[x][x], x] is Derivative[1][g[x]][x] +
         * Derivative[1][g][x][x]. (FeynCalc's FCLoopGLIDifferentiate
         * differentiates GLI[topo, {n..}][m] this way.) */
        Expr* dhead = deriv_of(head, x, nonconsts);
        if (!dhead) return NULL;
        dhead = eval_and_free(dhead);
        Expr* through_args = chain_rule_unknown(f, x, nonconsts);
        if (!through_args) { expr_free(dhead); return NULL; }
        if (dhead->type == EXPR_INTEGER && dhead->data.integer == 0) {
            expr_free(dhead);
            return through_args;
        }
        Expr** a = malloc((n ? n : 1) * sizeof(Expr*));
        for (size_t i = 0; i < n; i++) a[i] = expr_copy(args[i]);
        Expr* through_head = expr_new_function(dhead, a, n);
        free(a);
        return mk_fn2("Plus", through_args, through_head);
    }

    /* Give up -- caller keeps the expression unevaluated. */
    return NULL;
}

/* ---------------------------------------------------------------------- */
/* Higher-order: D[f, {x, n}]                                              */
/* ---------------------------------------------------------------------- */

/* Returns a new expression representing the n-th partial derivative of
 * f with respect to x. For n == 0 the input is deep-copied. */
static Expr* higher_order_partial(Expr* f, Expr* x, int64_t order,
                                  Expr* nonconsts) {
    if (order <= 0) return expr_copy(f);
    Expr* current = expr_copy(f);
    for (int64_t i = 0; i < order; i++) {
        Expr* nxt = compute_deriv(current, x, nonconsts);
        if (!nxt) {
            /* Cannot differentiate further -- leave as D[current, x]
             * and wrap any remaining orders as nested D calls. In
             * practice this only occurs for truly opaque heads. */
            Expr* wrap = mk_fn2("D", current, expr_copy(x));
            for (int64_t k = i + 1; k < order; k++) {
                wrap = mk_fn2("D", wrap, expr_copy(x));
            }
            return wrap;
        }
        /* The outer evaluator will also reduce the final expression,
         * but we need a simplified form between iterations so that
         * expr_free_of short-circuits work on the next pass.  Note that
         * evaluate() returns a new tree without freeing its argument;
         * we must explicitly free `nxt`. */
        Expr* reduced = evaluate(nxt);
        expr_free(nxt);
        expr_free(current);
        current = reduced;
    }
    return current;
}

/* ---------------------------------------------------------------------- */
/* Builtin: D                                                              */
/* ---------------------------------------------------------------------- */

/* Parse a second-argument spec which may be either
 *     x            -- a bare variable, order 1
 *     {x, n}       -- an integer-order spec
 *     {x, k}       -- a symbolic-order spec (k an Expr, not an Integer)
 *     {{x1,...,xN}}      -- an array (vector) spec, order 1 (gradient)
 *     {{x1,...,xN}, n}   -- an array spec, integer order n
 * On success, *var is set to a NON-owned pointer into the spec.
 *   - For an integer order, *order is set to the integer and *order_expr
 *     is NULL.
 *   - For a symbolic order, *order is set to the sentinel -1 and
 *     *order_expr is a NON-owned pointer to the spec's symbolic order
 *     argument.
 *   - For an array spec, *is_array is set to true and *var points at the
 *     inner array List expression (so its args are the variables).
 * Returns false if the spec cannot be interpreted. */
static bool parse_var_spec(Expr* spec, Expr** var, int64_t* order,
                           Expr** order_expr, bool* is_array) {
    if (order_expr) *order_expr = NULL;
    if (is_array) *is_array = false;
    if (spec->type == EXPR_FUNCTION &&
        is_sym(spec->data.function.head, "List")) {
        size_t ac = spec->data.function.arg_count;
        Expr* first = (ac >= 1) ? spec->data.function.args[0] : NULL;
        bool first_is_list = (first && first->type == EXPR_FUNCTION &&
                              is_sym(first->data.function.head, "List"));

        /* Array specs: {{x1,...,xN}} or {{x1,...,xN}, n}. */
        if (first_is_list) {
            if (ac == 1) {
                if (is_array) *is_array = true;
                *var = first;
                *order = 1;
                return true;
            }
            if (ac == 2) {
                Expr* k = spec->data.function.args[1];
                if (k->type == EXPR_INTEGER) {
                    if (k->data.integer < 0) return false;
                    if (is_array) *is_array = true;
                    *var = first;
                    *order = k->data.integer;
                    return true;
                }
                /* Symbolic / non-integer order with array spec: not handled. */
                return false;
            }
            return false;
        }

        /* Scalar specs: {var, n} or {var, k}. Any other shape
         * (e.g., {x}, {x, y, z}) is a malformed multiple-derivative
         * specifier; reject so the caller can emit D::dvar and leave
         * the call unevaluated (matching Mathematica). */
        if (ac != 2) {
            return false;
        }
        Expr* k = spec->data.function.args[1];
        if (k->type == EXPR_INTEGER) {
            *var = spec->data.function.args[0];
            *order = k->data.integer;
            return *order >= 0;
        }
        /* Reject non-integer numeric orders (Real, Rational, Complex);
         * MMA's D::dvar requires n to be symbolic or non-negative integer. */
        if (k->type == EXPR_REAL || k->type == EXPR_BIGINT ||
            (k->type == EXPR_FUNCTION &&
             k->data.function.head->type == EXPR_SYMBOL &&
             (k->data.function.head->data.symbol.name == SYM_Rational ||
              k->data.function.head->data.symbol.name == SYM_Complex))) {
            return false;
        }
        /* Symbolic order. */
        *var = spec->data.function.args[0];
        *order = -1;
        if (order_expr) *order_expr = k;
        return true;
    }
    *var = spec;
    *order = 1;
    return true;
}

/* ---------------------------------------------------------------------- */
/* Array (Outer-style) derivatives                                         */
/* ---------------------------------------------------------------------- */

/* Forward declaration. */
static Expr* outer_d_step(Expr* f, Expr* array, Expr* nonconsts);

/* For a scalar (non-List) ``leaf_of_f`` and an ``array`` that may itself
 * be a nested List, build the tensor whose leaves are D[leaf_of_f, v]
 * for each leaf v of ``array``. Preserves ``array``'s nested structure.
 * When the inner D cannot be reduced, we leave a literal D[leaf, v]
 * wrapper so the outer evaluator can re-attempt later. */
static Expr* d_at_array_leaves(Expr* leaf_of_f, Expr* array, Expr* nonconsts) {
    if (array->type == EXPR_FUNCTION &&
        is_sym(array->data.function.head, "List")) {
        size_t n = array->data.function.arg_count;
        Expr** items = malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) {
            items[i] = d_at_array_leaves(leaf_of_f,
                                         array->data.function.args[i],
                                         nonconsts);
        }
        return mk_fnN_adopt("List", items, n);
    }
    Expr* d = compute_deriv(leaf_of_f, array, nonconsts);
    if (!d) {
        /* Fallback: emit an unevaluated D[leaf, var] so the outer
         * evaluator can take another pass. Matches MMA's behaviour for
         * heads with no known derivative rule. */
        return mk_fn2("D", expr_copy(leaf_of_f), expr_copy(array));
    }
    return d;
}

/* Compute one array-derivative step: First[Outer[D, {f}, array]].
 *
 * This is equivalent to Outer[D, f, array] (the {f} wrapper plus First
 * cancel cleanly), so we recurse into f at each List level, and at every
 * non-List leaf we lay down a copy of ``array``'s skeleton populated
 * with D[leaf, v_i] entries. The result therefore has shape
 *     shape(f) ++ shape(array)
 * matching Mathematica's definition. */
static Expr* outer_d_step(Expr* f, Expr* array, Expr* nonconsts) {
    if (f->type == EXPR_FUNCTION &&
        is_sym(f->data.function.head, "List")) {
        size_t n = f->data.function.arg_count;
        Expr** items = malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) {
            items[i] = outer_d_step(f->data.function.args[i], array, nonconsts);
        }
        return mk_fnN_adopt("List", items, n);
    }
    return d_at_array_leaves(f, array, nonconsts);
}

/* Repeated array derivative: D[f, {array, n}] == apply outer_d_step n
 * times. Between iterations we evaluate so that the next pass sees a
 * simplified expression (mirrors higher_order_partial). */
static Expr* array_higher_order(Expr* f, Expr* array, int64_t order,
                                Expr* nonconsts) {
    if (order <= 0) return expr_copy(f);
    Expr* current = expr_copy(f);
    for (int64_t i = 0; i < order; i++) {
        Expr* nxt = outer_d_step(current, array, nonconsts);
        Expr* reduced = evaluate(nxt);
        expr_free(nxt);
        expr_free(current);
        current = reduced;
    }
    return current;
}

/* Emit MMA-style D::dvar message for malformed multiple-derivative
 * specifiers. Only called when parse_var_spec rejects a List spec. */
static void emit_dvar_message(Expr* spec) {
    if (spec->type != EXPR_FUNCTION ||
        !is_sym(spec->data.function.head, "List")) {
        return;
    }
    char* s = expr_to_string(spec);
    printf("D::dvar: Multiple derivative specifier %s does not have "
           "the form {variable, n}, where n is symbolic or a "
           "non-negative integer.\n", s ? s : "?");
    free(s);
}

/* Compute D[f, {var, k}] for a symbolic order k. The algorithm pattern-
 * matches the small set of forms with closed-form symbolic-order
 * derivatives, recurses through additive structure, pulls var-free
 * factors out of products, and otherwise returns NULL (caller falls
 * back to leaving D[f, {var, k}] unevaluated).
 *
 *   D[c, {x, k}]              = 0                          (c free of x)
 *   D[c f, {x, k}]            = c D[f, {x, k}]              (c free of x)
 *   D[a + b + ..., {x, k}]    = D[a, {x, k}] + ...
 *   D[Power[x, n], {x, k}]    = FactorialPower[n, k] x^(n-k)  (n free of x)
 *   D[Power[b, x], {x, k}]    = b^x Log[b]^k                  (b free of x)
 *
 * Sound for non-negative integer k (MMA's Piecewise conventions).
 * Returns a fresh Expr* on success, NULL when no symbolic-k closed form
 * applies. */
static Expr* compute_deriv_symbolic_order(Expr* f, Expr* var, Expr* k);

static Expr* compute_deriv_symbolic_order(Expr* f, Expr* var, Expr* k) {
    if (!f || !var || !k) return NULL;

    /* Var-free constant: D[c, {x, k}] = 0 for k >= 1. MMA also gives 0;
     * we match. (Strictly Piecewise[{{c, k==0}}, 0], but our convention
     * matches the unsubscripted MMA result.) */
    if (expr_free_of(f, var)) return mk_int(0);

    /* The variable itself: D[x, {x, k}] = If[k == 1, 1, 0] but MMA
     * returns Piecewise; for the symbolic-k path we leave unevaluated
     * (no general single-Expr closed form). */
    if (expr_eq(f, var)) return NULL;

    /* Plus: distribute. */
    if (f->type == EXPR_FUNCTION && is_sym(f->data.function.head, "Plus")) {
        size_t n = f->data.function.arg_count;
        Expr** parts = malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) {
            Expr* d = compute_deriv_symbolic_order(f->data.function.args[i], var, k);
            if (!d) {
                for (size_t j = 0; j < i; j++) expr_free(parts[j]);
                free(parts);
                return NULL;
            }
            parts[i] = d;
        }
        Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), parts, n);
        free(parts);
        return eval_and_free(sum);
    }

    /* Times: pull out var-free factors. If all but one factor is var-
     * free, the var-bearing factor's symbolic-k derivative scales by
     * the product of the var-free factors. */
    if (f->type == EXPR_FUNCTION && is_sym(f->data.function.head, "Times")) {
        size_t n = f->data.function.arg_count;
        Expr** free_factors = malloc(sizeof(Expr*) * n);
        size_t free_count = 0;
        Expr* var_factor = NULL;
        size_t var_factor_count = 0;
        for (size_t i = 0; i < n; i++) {
            if (expr_free_of(f->data.function.args[i], var)) {
                free_factors[free_count++] = expr_copy(f->data.function.args[i]);
            } else {
                var_factor = f->data.function.args[i];
                var_factor_count++;
            }
        }
        if (var_factor_count == 0 || var_factor_count > 1) {
            /* All var-free already handled above; >1 var-bearing factors
             * have no general symbolic-k form (Leibniz produces a
             * binomial sum that the user can request via Sum). */
            for (size_t i = 0; i < free_count; i++) expr_free(free_factors[i]);
            free(free_factors);
            return NULL;
        }
        Expr* d_var = compute_deriv_symbolic_order(var_factor, var, k);
        if (!d_var) {
            for (size_t i = 0; i < free_count; i++) expr_free(free_factors[i]);
            free(free_factors);
            return NULL;
        }
        if (free_count == 0) { free(free_factors); return d_var; }
        Expr** product_args = malloc(sizeof(Expr*) * (free_count + 1));
        for (size_t i = 0; i < free_count; i++) product_args[i] = free_factors[i];
        product_args[free_count] = d_var;
        Expr* prod = expr_new_function(expr_new_symbol(SYM_Times),
                                       product_args, free_count + 1);
        free(product_args);
        free(free_factors);
        return eval_and_free(prod);
    }

    /* Power forms with symbolic-k closed-form derivatives. */
    if (f->type == EXPR_FUNCTION && is_sym(f->data.function.head, "Power") &&
        f->data.function.arg_count == 2) {
        Expr* base = f->data.function.args[0];
        Expr* exp  = f->data.function.args[1];

        /* Power[var, n] with n free of var:
         *   D[x^n, {x, k}] = FactorialPower[n, k] * Power[x, n - k]. */
        if (expr_eq(base, var) && expr_free_of(exp, var)) {
            Expr* fp = mk_fn2("FactorialPower", expr_copy(exp), expr_copy(k));
            Expr* new_exp = expr_new_function(expr_new_symbol(SYM_Plus),
                (Expr*[]){ expr_copy(exp), mk_neg(expr_copy(k)) }, 2);
            Expr* new_pow = mk_fn2("Power", expr_copy(var), new_exp);
            Expr* product = expr_new_function(expr_new_symbol(SYM_Times),
                (Expr*[]){ fp, new_pow }, 2);
            return eval_and_free(product);
        }

        /* Power[b, var] with b free of var:
         *   D[b^x, {x, k}] = b^x * Log[b]^k. */
        if (expr_eq(exp, var) && expr_free_of(base, var)) {
            Expr* same  = mk_fn2("Power", expr_copy(base), expr_copy(var));
            Expr* logb  = mk_fn1("Log", expr_copy(base));
            Expr* logbk = mk_fn2("Power", logb, expr_copy(k));
            Expr* product = expr_new_function(expr_new_symbol(SYM_Times),
                (Expr*[]){ same, logbk }, 2);
            return eval_and_free(product);
        }
    }

    /* Sin/Cos/Sinh/Cosh of a linear-in-var argument.
     *
     * For u(var) = a*var + b with a, b free of var, the chain rule gives
     *   D^k/dvar^k F[u] = a^k * F^(k)[u]
     * with closed-form k-th derivatives matching Mathematica:
     *   Sin^(k)[t]  = Sin[k Pi/2 + t]
     *   Cos^(k)[t]  = Cos[k Pi/2 + t]
     *   Cosh^(k)[t] = (-I)^k Cos[k Pi/2 - I t]
     *   Sinh^(k)[t] = I (-I)^k Sin[k Pi/2 - I t]
     *
     * Linearity is detected by computing du/dvar; if that derivative is
     * free of var then u is linear in var and the residual is exactly a.
     * Non-linear arguments fall through (caller leaves D unevaluated). */
    if (f->type == EXPR_FUNCTION && f->data.function.arg_count == 1 &&
        f->data.function.head->type == EXPR_SYMBOL) {
        const char* h = f->data.function.head->data.symbol.name;
        bool is_trig  = (h == SYM_Sin) || (h == SYM_Cos);
        bool is_hyper = (h == SYM_Sinh) || (h == SYM_Cosh);
        if (is_trig || is_hyper) {
            Expr* u = f->data.function.args[0];
            Expr* du_raw = compute_deriv(u, var, NULL);
            if (du_raw) {
                Expr* a = evaluate(du_raw);
                expr_free(du_raw);
                if (expr_free_of(a, var)) {
                    /* half_pi_k = k * Pi / 2 */
                    Expr* half_pi_k = expr_new_function(expr_new_symbol(SYM_Times),
                        (Expr*[]){ expr_copy(k),
                                   expr_new_symbol(SYM_Pi),
                                   mk_fn2("Power", mk_int(2), mk_int(-1)) }, 3);
                    Expr* ak = mk_fn2("Power", a, expr_copy(k));
                    Expr* result;
                    if (is_trig) {
                        /* a^k * F[k*Pi/2 + u] */
                        Expr* shifted = mk_fn2("Plus", half_pi_k, expr_copy(u));
                        Expr* trig    = mk_fn1(h, shifted);
                        result = mk_fn2("Times", ak, trig);
                    } else {
                        /* a^k * (-I)^k * Cos[k*Pi/2 - I u]      (Cosh)
                         * a^k * I * (-I)^k * Sin[k*Pi/2 - I u]  (Sinh) */
                        Expr* minus_I_u = mk_fn2("Times",
                            mk_fn2("Times", mk_int(-1), expr_new_symbol(SYM_I)),
                            expr_copy(u));
                        Expr* arg     = mk_fn2("Plus", half_pi_k, minus_I_u);
                        const char* trig_name = (h == SYM_Cosh) ? "Cos" : "Sin";
                        Expr* trig    = mk_fn1(trig_name, arg);
                        Expr* neg_I   = mk_fn2("Times", mk_int(-1),
                                               expr_new_symbol(SYM_I));
                        Expr* neg_I_k = mk_fn2("Power", neg_I, expr_copy(k));
                        if (h == SYM_Cosh) {
                            result = expr_new_function(expr_new_symbol(SYM_Times),
                                (Expr*[]){ ak, neg_I_k, trig }, 3);
                        } else {
                            result = expr_new_function(expr_new_symbol(SYM_Times),
                                (Expr*[]){ ak, expr_new_symbol(SYM_I),
                                           neg_I_k, trig }, 4);
                        }
                    }
                    return eval_and_free(result);
                }
                expr_free(a);
            }
        }
    }

    /* HurwitzZeta[s, a] differentiated to symbolic order k in its second
     * argument a (s free of a = var):
     *   D[HurwitzZeta[s, a], {a, k}] = (-1)^k Pochhammer[s, k] HurwitzZeta[s+k, a].
     */
    if (f->type == EXPR_FUNCTION && is_sym(f->data.function.head, "HurwitzZeta") &&
        f->data.function.arg_count == 2) {
        Expr* s = f->data.function.args[0];
        Expr* a = f->data.function.args[1];
        if (expr_eq(a, var) && expr_free_of(s, var)) {
            Expr* sign = mk_fn2("Power", mk_int(-1), expr_copy(k));
            Expr* poch = mk_fn2("Pochhammer", expr_copy(s), expr_copy(k));
            Expr* spk  = mk_fn2("Plus", expr_copy(s), expr_copy(k));
            Expr* hz   = mk_fn2("HurwitzZeta", spk, expr_copy(a));
            Expr* prod = expr_new_function(expr_new_symbol(SYM_Times),
                (Expr*[]){ sign, poch, hz }, 3);
            return eval_and_free(prod);
        }
    }

    /* No symbolic-k closed form recognised for this shape. */
    return NULL;
}

/* Build an unevaluated D[f, {var, k}] tree. Used as the fallback when
 * compute_deriv_symbolic_order can't reduce the expression. */
static Expr* build_unevaluated_d(Expr* f, Expr* var, Expr* k) {
    Expr* spec = expr_new_function(expr_new_symbol(SYM_List),
        (Expr*[]){ expr_copy(var), expr_copy(k) }, 2);
    return expr_new_function(expr_new_symbol(SYM_D),
        (Expr*[]){ expr_copy(f), spec }, 2);
}

/* Detect Rule[NonConstants, val] / RuleDelayed[NonConstants, val].
 * Returns true when ``a`` is such an option, and writes a fresh
 * canonical List[sym1, ...] to ``*out`` (caller owns; we wrap a bare
 * single symbol into a singleton List to match Mathematica's output
 * canonicalisation). */
static bool parse_nonconsts_option(Expr* a, Expr** out) {
    if (!a || a->type != EXPR_FUNCTION) return false;
    if (a->data.function.arg_count != 2) return false;
    Expr* head = a->data.function.head;
    if (head->type != EXPR_SYMBOL) return false;
    if (head->data.symbol.name != SYM_Rule && head->data.symbol.name != SYM_RuleDelayed) {
        return false;
    }
    Expr* lhs = a->data.function.args[0];
    Expr* rhs = a->data.function.args[1];
    if (lhs->type != EXPR_SYMBOL || lhs->data.symbol.name != SYM_NonConstants) {
        return false;
    }
    if (rhs->type == EXPR_FUNCTION && is_sym(rhs->data.function.head, "List")) {
        *out = expr_copy(rhs);
    } else {
        *out = mk_fn1("List", expr_copy(rhs));
    }
    return true;
}

Expr* builtin_d(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2) return NULL;                 /* D[f] stays unevaluated */

    /* D is always a partial derivative. Pin total-derivative mode off for
     * the duration so a Dt-in-progress flag can never leak into a D node
     * evaluated re-entrantly (e.g. via eval_and_free inside the core). */
    bool saved_dt_mode = g_dt_total_mode;
    g_dt_total_mode = false;

    Expr* f = res->data.function.args[0];

    /* Split trailing arguments into option Rules and var-specs. Options
     * may appear interleaved (Mathematica conventionally places them at
     * the end, but we accept either order to match how users write
     * D[expr, x, NonConstants -> y]). */
    Expr** specs = malloc(sizeof(Expr*) * (argc - 1));
    size_t nspecs = 0;
    Expr* nonconsts = NULL;
    for (size_t i = 1; i < argc; i++) {
        Expr* a = res->data.function.args[i];
        Expr* parsed = NULL;
        if (parse_nonconsts_option(a, &parsed)) {
            if (nonconsts) expr_free(nonconsts);
            nonconsts = parsed;             /* later wins, mirroring MMA */
            continue;
        }
        specs[nspecs++] = a;
    }
    if (nspecs == 0) {
        if (nonconsts) expr_free(nonconsts);
        free(specs);
        g_dt_total_mode = saved_dt_mode;
        return NULL;                        /* D[f, opts...] stays unevaluated */
    }

    /* Fixed-point form D[sym, var, NonConstants -> {... sym ...}].
     * The internal compute_deriv emits exactly this shape when it hits
     * an implicit-dependence symbol; we return the canonical
     * List-wrapped form so that NonConstants -> y normalises to
     * NonConstants -> {y}. If the input is already canonical,
     * returning a structurally-equal result lets the outer fixed-point
     * detector stop without further work. */
    if (nonconsts && nspecs == 1 && f->type == EXPR_SYMBOL &&
        nonconsts_contains(nonconsts, f) &&
        specs[0]->type == EXPR_SYMBOL && !expr_eq(f, specs[0])) {
        Expr* r = build_unevaluated_d_nonconsts(f, specs[0], nonconsts);
        expr_free(nonconsts);
        free(specs);
        g_dt_total_mode = saved_dt_mode;
        return r;
    }

    /* Sequentially apply each spec (handling mixed partials D[f, x, y, ...]). */
    Expr* current = expr_copy(f);
    for (size_t i = 0; i < nspecs; i++) {
        Expr* var = NULL;
        int64_t order = 1;
        Expr* order_expr = NULL;
        bool is_array = false;
        if (!parse_var_spec(specs[i], &var, &order, &order_expr, &is_array)) {
            emit_dvar_message(specs[i]);       /* MMA-style D::dvar */
            if (nonconsts) expr_free(nonconsts);
            free(specs);
            expr_free(current);
            g_dt_total_mode = saved_dt_mode;
            return NULL;                       /* malformed spec */
        }
        Expr* stepped = NULL;
        if (is_array) {
            /* {{x1,...,xN}} or {{x1,...,xN}, n}: array (Outer) derivative. */
            stepped = array_higher_order(current, var, order, nonconsts);
        } else if (order >= 0) {
            stepped = higher_order_partial(current, var, order, nonconsts);
        } else {
            /* Symbolic-order spec. NonConstants is not threaded through
             * the closed-form symbolic-k path; it falls through to the
             * unevaluated form for any mixed case. */
            Expr* sym = compute_deriv_symbolic_order(current, var, order_expr);
            if (sym) stepped = sym;
            else stepped = build_unevaluated_d(current, var, order_expr);
        }
        expr_free(current);
        current = stepped;
    }

    if (nonconsts) expr_free(nonconsts);
    free(specs);
    g_dt_total_mode = saved_dt_mode;
    return current;
}

/* ---------------------------------------------------------------------- */
/* Builtin: Dt                                                             */
/* ---------------------------------------------------------------------- */

Expr* builtin_dt(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc == 0) return NULL;

    Expr* f = res->data.function.args[0];

    /* Dt[f]: total derivative. */
    if (argc == 1) {
        return compute_deriv(f, NULL, NULL);  /* may be NULL to stay unevaluated */
    }

    /* Terminal case: the total derivative of a bare symbol with respect
     * to a single bare variable is its own normal form. compute_deriv (in
     * g_dt_total_mode) re-emits Dt[s, x] for such a symbol, and
     * higher_order_partial then evaluates that result -- which would
     * re-enter this builtin and recurse without bound. Resolve the closed
     * cases here and leave the genuinely irreducible form unevaluated. */
    if (argc == 2 && f->type == EXPR_SYMBOL) {
        Expr* v = res->data.function.args[1];
        if (is_dt_constant_symbol(f)) return mk_int(0);   /* Dt[Pi, x] = 0 */
        if (v->type == EXPR_SYMBOL) {
            if (f->data.symbol.name == v->data.symbol.name) return mk_int(1);  /* Dt[x, x] */
            return NULL;                                  /* Dt[s, x] stays */
        }
    }

    /* Dt[f, x, ...] is the *total* derivative with respect to each
     * variable: unlike D, every symbol other than the differentiation
     * variable and the distinguished constants is treated as an implicit
     * function of that variable. We reuse the shared partial-derivative
     * machinery with g_dt_total_mode enabled so free symbols contribute
     * Dt[sym, x] terms instead of vanishing. The flag is saved/restored to
     * stay correct under nested (re-entrant) differentiation. */
    bool saved_dt_mode = g_dt_total_mode;
    g_dt_total_mode = true;
    Expr* current = expr_copy(f);
    for (size_t i = 1; i < argc; i++) {
        Expr* var = NULL;
        int64_t order = 1;
        Expr* order_expr = NULL;
        bool is_array = false;
        if (!parse_var_spec(res->data.function.args[i], &var, &order, &order_expr,
                            &is_array)) {
            emit_dvar_message(res->data.function.args[i]);
            expr_free(current);
            g_dt_total_mode = saved_dt_mode;
            return NULL;
        }
        Expr* stepped = NULL;
        if (is_array) {
            stepped = array_higher_order(current, var, order, NULL);
        } else if (order >= 0) {
            stepped = higher_order_partial(current, var, order, NULL);
        } else {
            Expr* sym = compute_deriv_symbolic_order(current, var, order_expr);
            stepped = sym ? sym : build_unevaluated_d(current, var, order_expr);
        }
        expr_free(current);
        current = stepped;
    }

    g_dt_total_mode = saved_dt_mode;
    return current;
}

/* ---------------------------------------------------------------------- */
/* Reduce Derivative[n1,...,nm][Function[...]] to a plain Function[body']  */
/* ---------------------------------------------------------------------- */

/* Build a fresh Slot[i] expression. */
static Expr* make_slot(int64_t i) {
    Expr* idx[1] = { mk_int(i) };
    return expr_new_function(mk_sym("Slot"), idx, 1);
}

/* See deriv.h for the contract. */
Expr* derivative_of_pure_function(Expr* deriv_head, Expr* pure_fn) {
    if (!deriv_head || !pure_fn) return NULL;
    if (deriv_head->type != EXPR_FUNCTION) return NULL;
    if (!is_sym(deriv_head->data.function.head, "Derivative")) return NULL;
    size_t m = deriv_head->data.function.arg_count;
    if (m == 0) return NULL;

    /* All derivative orders must be nonnegative integers. */
    int64_t* orders = malloc(sizeof(int64_t) * m);
    for (size_t i = 0; i < m; i++) {
        Expr* oi = deriv_head->data.function.args[i];
        if (oi->type != EXPR_INTEGER || oi->data.integer < 0) {
            free(orders);
            return NULL;
        }
        orders[i] = oi->data.integer;
    }

    if (pure_fn->type != EXPR_FUNCTION ||
        !is_sym(pure_fn->data.function.head, "Function")) {
        free(orders);
        return NULL;
    }

    size_t fargc = pure_fn->data.function.arg_count;
    if (fargc == 0) { free(orders); return NULL; }

    /* Determine body and the variables to differentiate with respect to.
     * For each supported Function signature we identify m variables
     * (matching the arity of Derivative) and remember whether the rebuilt
     * Function preserves the original parameter/attribute shape. */
    Expr* body = NULL;
    Expr** vars = malloc(sizeof(Expr*) * m);
    bool slot_form = false;

    if (fargc == 1) {
        /* Function[body]: slot form. Differentiate wrt Slot[1..m]. */
        body = pure_fn->data.function.args[0];
        for (size_t i = 0; i < m; i++) vars[i] = make_slot((int64_t)(i + 1));
        slot_form = true;
    } else {
        /* Function[params, body, ...]. */
        Expr* params = pure_fn->data.function.args[0];
        body = pure_fn->data.function.args[1];

        if (params->type == EXPR_SYMBOL &&
            params->data.symbol.name == SYM_Null) {
            /* Function[Null, body, attrs]: slot form with attributes. */
            for (size_t i = 0; i < m; i++) vars[i] = make_slot((int64_t)(i + 1));
        } else if (params->type == EXPR_SYMBOL) {
            /* Function[x, body]: single named parameter. Only Derivative[n]
             * (m == 1) is meaningful here. */
            if (m != 1) { free(vars); free(orders); return NULL; }
            vars[0] = expr_copy(params);
        } else if (params->type == EXPR_FUNCTION &&
                   is_sym(params->data.function.head, "List")) {
            /* Function[{x1,...,xk}, body]: arity of Derivative must match. */
            size_t k = params->data.function.arg_count;
            if (k != m) { free(vars); free(orders); return NULL; }
            for (size_t i = 0; i < m; i++) {
                vars[i] = expr_copy(params->data.function.args[i]);
            }
        } else {
            free(vars);
            free(orders);
            return NULL;
        }
    }

    /* Apply each partial derivative sequentially, simplifying between
     * iterations so subsequent passes can see a reduced expression. */
    Expr* current = expr_copy(body);
    for (size_t i = 0; i < m; i++) {
        for (int64_t k = 0; k < orders[i]; k++) {
            Expr* d = compute_deriv(current, vars[i], NULL);
            expr_free(current);
            if (!d) {
                for (size_t j = 0; j < m; j++) expr_free(vars[j]);
                free(vars);
                free(orders);
                return NULL;
            }
            current = evaluate(d);
            expr_free(d);
        }
    }

    for (size_t j = 0; j < m; j++) expr_free(vars[j]);
    free(vars);
    free(orders);

    /* Rebuild a Function with the original parameter/attribute shape but
     * the differentiated body. */
    if (slot_form) {
        Expr* fn_args[1] = { current };
        return expr_new_function(mk_sym("Function"), fn_args, 1);
    }

    size_t new_argc = fargc;
    Expr** new_args = malloc(sizeof(Expr*) * new_argc);
    new_args[0] = expr_copy(pure_fn->data.function.args[0]);
    new_args[1] = current;
    for (size_t i = 2; i < new_argc; i++) {
        new_args[i] = expr_copy(pure_fn->data.function.args[i]);
    }
    Expr* out = expr_new_function(mk_sym("Function"), new_args, new_argc);
    free(new_args);
    return out;
}

/* ---------------------------------------------------------------------- */
/* Reduce Derivative[n1,...,nm][f] for a symbol f with DownValues          */
/* ---------------------------------------------------------------------- */

/* Counter used to mint unique temporary symbol names when synthesising
 * the Function[{t1,...,tm}, f[t1,...,tm]] wrapper. The temporaries never
 * have OwnValues installed and are not registered with the symbol table,
 * so they cannot collide with anything in the body produced by the
 * DownValue rewrite (which would have to literally contain the same
 * unique name to clash, and the names are generated fresh per call). */
static int64_t deriv_temp_counter = 0;

/* See deriv.h for the contract. */
Expr* derivative_of_symbol(Expr* deriv_head, Expr* fsym) {
    if (!deriv_head || !fsym) return NULL;
    if (deriv_head->type != EXPR_FUNCTION) return NULL;
    if (!is_sym(deriv_head->data.function.head, "Derivative")) return NULL;
    if (fsym->type != EXPR_SYMBOL) return NULL;

    size_t m = deriv_head->data.function.arg_count;
    if (m == 0) return NULL;

    /* All derivative orders must be nonnegative integers; otherwise the
     * downstream derivative_of_pure_function would reject them and we
     * would have done substantial work for nothing. */
    for (size_t i = 0; i < m; i++) {
        Expr* oi = deriv_head->data.function.args[i];
        if (oi->type != EXPR_INTEGER || oi->data.integer < 0) return NULL;
    }

    /* Only proceed when the symbol carries at least one DownValue. For a
     * bare symbol with no rules, attempting to substitute would just give
     * us back f[t1,...,tm] and we would correctly fall through; the
     * early exit here is purely a fast path that avoids the call/eval
     * roundtrip in the common no-definition case. */
    if (!symtab_get_down_values(fsym->data.symbol.name)) return NULL;

    /* Mint fresh temporary variable symbols. We use a static counter so
     * names are globally unique per process invocation; this avoids the
     * need to register temporaries (and clean them up) in the symbol
     * table just to guarantee freshness. */
    int64_t my_id = ++deriv_temp_counter;

    Expr** vars = malloc(sizeof(Expr*) * m);
    for (size_t i = 0; i < m; i++) {
        char namebuf[64];
        snprintf(namebuf, sizeof(namebuf),
                 "Derivative$%lld$%zu", (long long)my_id, i + 1);
        vars[i] = mk_sym(namebuf);
    }

    /* Build f[t1,...,tm] and evaluate it. The DownValue lookup happens
     * during evaluation, so we get back the substituted body when a
     * single-symbol-argument rule like f[x_] := body matches. */
    Expr** call_args = malloc(sizeof(Expr*) * m);
    for (size_t i = 0; i < m; i++) call_args[i] = expr_copy(vars[i]);
    Expr* call = expr_new_function(expr_copy(fsym), call_args, m);
    free(call_args);

    Expr* body = evaluate(call);
    expr_free(call);

    /* If the call did not rewrite (no matching DownValue), abort. We
     * detect this by checking that the result is still structurally
     * f[t1,...,tm]. */
    bool unchanged = (body->type == EXPR_FUNCTION
                      && body->data.function.head->type == EXPR_SYMBOL
                      && body->data.function.head->data.symbol.name == fsym->data.symbol.name
                      && body->data.function.arg_count == m);
    if (unchanged) {
        for (size_t i = 0; i < m; i++) {
            if (!expr_eq(body->data.function.args[i], vars[i])) {
                unchanged = false;
                break;
            }
        }
    }
    if (unchanged) {
        expr_free(body);
        for (size_t i = 0; i < m; i++) expr_free(vars[i]);
        free(vars);
        return NULL;
    }

    /* Wrap the substituted body in Function[{t1,...,tm}, body] so we can
     * reuse the existing pure-function differentiation pipeline. */
    Expr** list_args = malloc(sizeof(Expr*) * m);
    for (size_t i = 0; i < m; i++) list_args[i] = expr_copy(vars[i]);
    Expr* var_list = mk_fnN_adopt("List", list_args, m);

    Expr* fn_args[2] = { var_list, body };
    Expr* pure_fn = expr_new_function(mk_sym("Function"), fn_args, 2);

    Expr* result = derivative_of_pure_function(deriv_head, pure_fn);
    expr_free(pure_fn);

    for (size_t i = 0; i < m; i++) expr_free(vars[i]);
    free(vars);

    return result;
}

/* ---------------------------------------------------------------------- */
/* Builtin: Derivative                                                     */
/* ---------------------------------------------------------------------- */

/*
 * The Derivative symbol primarily serves as a tag; the actual
 * differentiation happens inside the D dispatch. However, there is one
 * useful simplification we can perform at the head level:
 *
 *     Derivative[0, 0, ..., 0]  -->  Identity-like behaviour via the
 *     calling evaluator (the expression Derivative[0,..,0][f][args]
 *     is recognised by compute_deriv and does not need a builtin
 *     rewrite here).
 *
 * We therefore simply return NULL, leaving ``Derivative[n]`` in its
 * canonical unevaluated form. The value of this builtin is that it
 * lets us register attributes on the Derivative symbol.
 */
Expr* builtin_derivative(Expr* res) {
    (void)res;
    return NULL;
}

/* Registry of numeric evaluators for inert derivatives (see deriv.h).  A
 * handful of entries at most, so a linear scan of interned pointers. */
#define DERIV_NUMERIC_MAX 32
static struct { const char* fname; DerivNumericFn fn; size_t nargs; } deriv_numeric_tab[DERIV_NUMERIC_MAX];
static size_t deriv_numeric_n = 0;

void deriv_register_numeric(const char* fname, DerivNumericFn fn, size_t nargs) {
    for (size_t i = 0; i < deriv_numeric_n; i++)
        if (deriv_numeric_tab[i].fname == fname && deriv_numeric_tab[i].nargs == nargs) {
            deriv_numeric_tab[i].fn = fn;
            return;
        }
    if (deriv_numeric_n < DERIV_NUMERIC_MAX) {
        deriv_numeric_tab[deriv_numeric_n].fname = fname;
        deriv_numeric_tab[deriv_numeric_n].fn = fn;
        deriv_numeric_tab[deriv_numeric_n].nargs = nargs;
        deriv_numeric_n++;
    }
}

DerivNumericFn deriv_numeric_for(const Expr* dop) {
    if (deriv_numeric_n == 0 || !dop || dop->type != EXPR_FUNCTION ||
        dop->data.function.arg_count != 1) return NULL;
    const Expr* d = dop->data.function.head;            /* Derivative[n..] */
    const Expr* f = dop->data.function.args[0];
    if (!d || d->type != EXPR_FUNCTION || !d->data.function.head ||
        d->data.function.head->type != EXPR_SYMBOL ||
        d->data.function.head->data.symbol.name != SYM_Derivative ||
        f->type != EXPR_SYMBOL) return NULL;
    for (size_t i = 0; i < deriv_numeric_n; i++)
        if (deriv_numeric_tab[i].fname == f->data.symbol.name
            && deriv_numeric_tab[i].nargs == d->data.function.arg_count)
            return deriv_numeric_tab[i].fn;
    return NULL;
}

/* True when e mentions an inert derivative of the symbol fname
 * (Derivative[..][fname]) or an unevaluated D[...]: the closed form below is
 * then not a closed form. */
static bool deriv_mentions_inert(const Expr* e, const char* fname) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* hd = e->data.function.head;
    if (hd->type == EXPR_SYMBOL && hd->data.symbol.name == SYM_D) return true;
    if (hd->type == EXPR_FUNCTION && hd->data.function.arg_count == 1
        && hd->data.function.head->type == EXPR_FUNCTION
        && hd->data.function.head->data.function.head->type == EXPR_SYMBOL
        && hd->data.function.head->data.function.head->data.symbol.name == SYM_Derivative
        && hd->data.function.args[0]->type == EXPR_SYMBOL
        && hd->data.function.args[0]->data.symbol.name == fname)
        return true;
    if (deriv_mentions_inert(hd, fname)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (deriv_mentions_inert(e->data.function.args[i], fname)) return true;
    return false;
}

/* Closed forms of Derivative[n1..nm][f] for a built-in numeric function f,
 * memoised per operator: the derivative of f[t1..tm] in placeholder slots, or
 * NULL when D has no closed form for it (D[BesselJ[n, z], n] stays inert). */
#define DERIV_CF_MEMO 64
static struct { Expr* dop; Expr* body; } deriv_cf_memo[DERIV_CF_MEMO];
static int deriv_cf_memo_n = 0, deriv_cf_memo_next = 0;
#define DERIV_CF_MAXARGS 8

/* Derivative[n1..nm][f][a1..am] for a built-in NumericFunction f (Sin, Gamma,
 * BesselJ, ...): the value D gives, Derivative[1][Sin][x] -> Cos[x] and
 * Derivative[0, 1][BesselJ][n, z] -> (BesselJ[-1 + n, z] - BesselJ[1 + n, z])/2,
 * as in the reference.  Without it every applied derivative of a built-in
 * stayed inert -- N[Derivative[0, 1][BesselJ][3/8, 5/8]] and
 * Derivative[1][Gamma][0.5] were left unevaluated.  Declines (NULL) for a
 * user symbol (no builtin: the inert form is the answer DSolve and friends
 * build on), a non-numeric builtin (evaluating f[t] could have side effects),
 * and a derivative with no closed form. */
static Expr* deriv_closed_form(const Expr* call) {
    const Expr* dop = call->data.function.head;
    if (dop->type != EXPR_FUNCTION || dop->data.function.arg_count != 1) return NULL;
    const Expr* d = dop->data.function.head;
    const Expr* f = dop->data.function.args[0];
    if (d->type != EXPR_FUNCTION || d->data.function.head->type != EXPR_SYMBOL
        || d->data.function.head->data.symbol.name != SYM_Derivative
        || f->type != EXPR_SYMBOL) return NULL;
    size_t m = d->data.function.arg_count;
    if (m == 0 || m > DERIV_CF_MAXARGS || call->data.function.arg_count != m) return NULL;
    int64_t total = 0;
    for (size_t i = 0; i < m; i++) {
        const Expr* oi = d->data.function.args[i];
        if (oi->type != EXPR_INTEGER || oi->data.integer < 0) return NULL;
        total += oi->data.integer;
    }
    if (total <= 0 || total > 16) return NULL;
    SymbolDef* def = f->data.symbol.def ? f->data.symbol.def
                                        : symtab_lookup(f->data.symbol.name);
    if (!def || !def->builtin_func || !(def->attributes & ATTR_NUMERICFUNCTION))
        return NULL;

    static const char* slot[DERIV_CF_MAXARGS];
    if (!slot[0])
        for (int i = 0; i < DERIV_CF_MAXARGS; i++) {
            char nb[32];
            snprintf(nb, sizeof nb, "Derivative$Slot%d", i + 1);
            slot[i] = intern_symbol(nb);
        }

    Expr* body = NULL;
    int hit = -1;
    for (int i = 0; i < deriv_cf_memo_n; i++)
        if (expr_eq(deriv_cf_memo[i].dop, dop)) { hit = i; break; }
    if (hit >= 0) {
        body = deriv_cf_memo[hit].body;
    } else {
        /* Reserve the entry as a decline BEFORE differentiating: D of f[t..]
         * may itself produce Derivative[..][f][t..] (no closed form), whose
         * evaluation re-enters here and must find "no closed form" rather
         * than recurse without bound. */
        int at;
        if (deriv_cf_memo_n < DERIV_CF_MEMO) at = deriv_cf_memo_n++;
        else {
            at = deriv_cf_memo_next;
            deriv_cf_memo_next = (deriv_cf_memo_next + 1) % DERIV_CF_MEMO;
            expr_free(deriv_cf_memo[at].dop);
            expr_free(deriv_cf_memo[at].body);
        }
        Expr* key = expr_copy((Expr*)dop);
        deriv_cf_memo[at].dop = key;
        deriv_cf_memo[at].body = NULL;
        Expr** fa = malloc(m * sizeof(Expr*));
        for (size_t i = 0; i < m; i++) fa[i] = expr_new_symbol(slot[i]);
        Expr* g = expr_new_function(expr_copy((Expr*)f), fa, m);
        free(fa);
        size_t nspec = 0;
        Expr** da = malloc((m + 1) * sizeof(Expr*));
        da[nspec++] = g;
        for (size_t i = 0; i < m; i++) {
            int64_t k = d->data.function.args[i]->data.integer;
            if (k == 0) continue;
            da[nspec++] = mk_fn2("List", expr_new_symbol(slot[i]), mk_int(k));
        }
        Expr* r = eval_and_free(expr_new_function(mk_sym("D"), da, nspec));
        free(da);
        if (r && deriv_mentions_inert(r, f->data.symbol.name)) { expr_free(r); r = NULL; }
        body = r;
        /* The slot may have been recycled by a nested insertion meanwhile. */
        if (deriv_cf_memo[at].dop == key && !deriv_cf_memo[at].body)
            deriv_cf_memo[at].body = body ? expr_copy(body) : NULL;
        if (body) {
            Expr** rules0 = malloc(m * sizeof(Expr*));
            for (size_t i = 0; i < m; i++)
                rules0[i] = mk_fn2("Rule", expr_new_symbol(slot[i]),
                                   expr_copy(call->data.function.args[i]));
            Expr* rl0 = mk_fnN_adopt("List", rules0, m);
            return eval_and_free(mk_fn2("ReplaceAll", body, rl0));
        }
        return NULL;
    }
    if (!body) return NULL;
    Expr** rules = malloc(m * sizeof(Expr*));
    for (size_t i = 0; i < m; i++)
        rules[i] = mk_fn2("Rule", expr_new_symbol(slot[i]),
                          expr_copy(call->data.function.args[i]));
    Expr* rl = mk_fnN_adopt("List", rules, m);
    return eval_and_free(mk_fn2("ReplaceAll", expr_copy(body), rl));
}

Expr* deriv_apply_numeric(const Expr* call) {
    if (!call || call->type != EXPR_FUNCTION) return NULL;
    Expr* cf = deriv_closed_form(call);
    if (cf) return cf;
    if (deriv_numeric_n == 0) return NULL;
    DerivNumericFn fn = deriv_numeric_for(call->data.function.head);
    return fn ? fn(call) : NULL;
}

/* ---------------------------------------------------------------------- */
/* Registration                                                            */
/* ---------------------------------------------------------------------- */

void deriv_init(void) {
    symtab_add_builtin("D", builtin_d);
    symtab_add_builtin("Dt", builtin_dt);
    symtab_add_builtin("Derivative", builtin_derivative);

    /* Match the original deriv.m attribute set:
     *     SetAttributes[D, {Protected}]
     *     SetAttributes[Dt, {Protected}]
     *     SetAttributes[Derivative, {Protected}] */
    symtab_get_def("D")->attributes          |= ATTR_PROTECTED;
    symtab_get_def("Dt")->attributes         |= ATTR_PROTECTED;
    symtab_get_def("Derivative")->attributes |= ATTR_NHOLDALL;   /* not Protected: packages define Derivative[...][f][x] := ... SubValues */
}
