/*
 * complex_expand.c  --  ComplexExpand
 *
 * See complex_expand.h for the user-facing contract.
 *
 * Architecture
 * ------------
 * The engine is one recursive routine, cx_decompose(e, ctx, &re, &im), that
 * writes real-valued expressions re, im with  e == re + I*im  under the
 * decomposition context ctx (which symbols are complex, and the output
 * TargetFunctions basis).  Every other operation is a thin wrapper: the
 * builtin front-end evaluates its argument, threads over lists / relations,
 * runs cx_decompose, and assembles Expand[re + I*im] (or pulls out one
 * component for a Re/Im/Abs/Arg/... wrapper, which cx_decompose already
 * handles as ordinary nodes).
 *
 * TargetFunctions -> {Re, Im} and {Abs, Arg} both flow through the (re, im)
 * engine; they differ only in how a *complex atom* is decomposed (the single
 * substitution point cx_atom_reim()).  TargetFunctions -> Conjugate is a
 * separate, simpler path: it conjugates the whole expression (I -> -I,
 * z -> Conjugate[z]) and averages, which reproduces the z^2/2 + Conjugate[z]^2/2
 * family directly.
 *
 * Memory
 * ------
 * All cx_* helpers BORROW their Expr* arguments and return freshly-owned,
 * evaluated Expr*.  The builtin never frees `res` (the evaluator owns it).
 */

#include "complex_expand.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "symtab.h"
#include "eval.h"
#include "attr.h"
#include "arithmetic.h"
#include "common.h"
#include "sym_names.h"
#include "match.h"

/* ------------------------------------------------------------------ *
 *  Expression-building shorthand.  mk_* take ownership of pointer      *
 *  arguments (no deep copy); callers pass fresh or expr_copy()'d nodes.*
 * ------------------------------------------------------------------ */

static Expr* mk_int(int64_t v)          { return expr_new_integer(v); }
static Expr* mk_sym(const char* s)      { return expr_new_symbol(s); }
static Expr* mk_fn1(const char* h, Expr* a)          { return expr_new_function(mk_sym(h), (Expr*[]){ a }, 1); }
static Expr* mk_fn2(const char* h, Expr* a, Expr* b) { return expr_new_function(mk_sym(h), (Expr*[]){ a, b }, 2); }
static Expr* mk_pow(Expr* b, Expr* e)   { return mk_fn2("Power", b, e); }
static Expr* mk_half(void)              { return mk_fn2("Rational", mk_int(1), mk_int(2)); }
static Expr* mk_I(void)                 { return make_complex(mk_int(0), mk_int(1)); }

/* Evaluate and free (consume) e. */
static Expr* ev(Expr* e) { return eval_and_free(e); }

/* Build head[args...] copying nothing (adopts args[]); helper for N args. */
static Expr* mk_fnN(const char* h, Expr** args, size_t n) {
    return expr_new_function(mk_sym(h), args, n);
}

/* ------------------------------------------------------------------ *
 *  Small numeric predicates.                                          *
 * ------------------------------------------------------------------ */

/* Structural zero: Integer 0 or Real 0.0. */
static bool cx_is_zero(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) return e->data.integer == 0;
    if (e->type == EXPR_REAL)    return e->data.real == 0.0;
    return false;
}

/* True iff `e` is a positive numeric literal (Integer > 0, Real > 0,
 * Rational[p,q] with p,q > 0).  Used to keep real^real real. */
static bool cx_is_pos_numeric(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) return e->data.integer > 0;
    if (e->type == EXPR_REAL)    return e->data.real > 0.0;
    int64_t n, d;
    if (is_rational(e, &n, &d)) return (n > 0 && d > 0) || (n < 0 && d < 0);
    if (e->type == EXPR_SYMBOL) {
        /* Pi, E, EulerGamma, GoldenRatio, Degree, Catalan are positive reals. */
        const char* s = e->data.symbol.name;
        if (s == SYM_Pi || s == SYM_E) return true;
    }
    return false;
}

/* True iff `e` is a plain (small or big) integer literal, writing it to *n
 * when it fits in int64. */
static bool cx_is_int_literal(const Expr* e, long* n) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *n = (long)e->data.integer; return true; }
    return false;
}

/* ------------------------------------------------------------------ *
 *  Decomposition context.                                             *
 * ------------------------------------------------------------------ */

typedef enum { TF_REIM, TF_ABSARG, TF_CONJUGATE } TargetMode;

typedef struct {
    Expr** cvars;    /* borrowed complex-variable patterns */
    size_t ncvars;
    TargetMode target;
} CxCtx;

/* True iff atom `e` matches one of the complex-variable patterns in ctx. */
static bool cx_is_complex_atom(const Expr* e, const CxCtx* ctx) {
    for (size_t i = 0; i < ctx->ncvars; i++) {
        MatchEnv* env = env_new();
        bool ok = match((Expr*)e, ctx->cvars[i], env);
        env_free(env);
        if (ok) return true;
    }
    return false;
}

/* ------------------------------------------------------------------ *
 *  Complex arithmetic on (re, im) pairs.  Borrow inputs, own outputs. *
 * ------------------------------------------------------------------ */

static Expr* cx_add2(Expr* a, Expr* b) {
    return ev(mk_fn2("Plus", expr_copy(a), expr_copy(b)));
}
static Expr* cx_mul2(Expr* a, Expr* b) {
    return ev(mk_fn2("Times", expr_copy(a), expr_copy(b)));
}
static Expr* cx_sub2(Expr* a, Expr* b) {
    return ev(mk_fn2("Plus", expr_copy(a), mk_fn2("Times", mk_int(-1), expr_copy(b))));
}

/* (ar + I ai) * (br + I bi) -> (rr, ri). */
static void cx_cmul(Expr* ar, Expr* ai, Expr* br, Expr* bi, Expr** rr, Expr** ri) {
    Expr* ac = cx_mul2(ar, br);
    Expr* bd = cx_mul2(ai, bi);
    Expr* ad = cx_mul2(ar, bi);
    Expr* bc = cx_mul2(ai, br);
    *rr = cx_sub2(ac, bd);
    *ri = cx_add2(ad, bc);
    expr_free(ac); expr_free(bd); expr_free(ad); expr_free(bc);
}

/* 1 / (ar + I ai) -> (rr, ri) = (ar/D, -ai/D), D = ar^2 + ai^2. */
static void cx_crecip(Expr* ar, Expr* ai, Expr** rr, Expr** ri) {
    Expr* a2 = ev(mk_pow(expr_copy(ar), mk_int(2)));
    Expr* b2 = ev(mk_pow(expr_copy(ai), mk_int(2)));
    Expr* den = cx_add2(a2, b2);
    expr_free(a2); expr_free(b2);
    *rr = ev(mk_fn2("Times", expr_copy(ar), mk_pow(expr_copy(den), mk_int(-1))));
    *ri = ev(mk_fn2("Times", mk_fn2("Times", mk_int(-1), expr_copy(ai)),
                    mk_pow(expr_copy(den), mk_int(-1))));
    expr_free(den);
}

/* (ar + I ai) / (br + I bi) -> (rr, ri). */
static void cx_cdiv(Expr* ar, Expr* ai, Expr* br, Expr* bi, Expr** rr, Expr** ri) {
    Expr *irr, *iri;
    cx_crecip(br, bi, &irr, &iri);
    cx_cmul(ar, ai, irr, iri, rr, ri);
    expr_free(irr); expr_free(iri);
}

/* Reconstruct the (possibly complex) value re + I*im. */
static Expr* cx_recon(Expr* re, Expr* im) {
    if (cx_is_zero(im)) return expr_copy(re);
    return ev(mk_fn2("Plus", expr_copy(re),
                     mk_fn2("Times", mk_I(), expr_copy(im))));
}

/* Sum of squares re^2 + im^2, expanded. */
static Expr* cx_sumsq(Expr* re, Expr* im) {
    Expr* r2 = ev(mk_pow(expr_copy(re), mk_int(2)));
    Expr* i2 = ev(mk_pow(expr_copy(im), mk_int(2)));
    Expr* s  = cx_add2(r2, i2);
    expr_free(r2); expr_free(i2);
    Expr* xs = ev(mk_fn1("Expand", s));
    return xs;
}

/* ------------------------------------------------------------------ *
 *  The engine.                                                        *
 * ------------------------------------------------------------------ */

static void cx_decompose(const Expr* e, const CxCtx* ctx, Expr** re, Expr** im);

/* Decompose an atom (symbol / number) that is *not* a Complex literal. */
static void cx_atom(const Expr* e, const CxCtx* ctx, Expr** re, Expr** im) {
    if (e->type == EXPR_SYMBOL && cx_is_complex_atom(e, ctx)) {
        if (ctx->target == TF_ABSARG) {
            Expr* ab = mk_fn1("Abs", expr_copy((Expr*)e));
            Expr* ar = mk_fn1("Arg", expr_copy((Expr*)e));
            *re = ev(mk_fn2("Times", expr_copy(ab), mk_fn1("Cos", expr_copy(ar))));
            *im = ev(mk_fn2("Times", ab, mk_fn1("Sin", ar)));
        } else { /* TF_REIM */
            *re = mk_fn1("Re", expr_copy((Expr*)e));
            *im = mk_fn1("Im", expr_copy((Expr*)e));
        }
        return;
    }
    /* Real (default): numbers, real constants, unlisted symbols. */
    *re = expr_copy((Expr*)e);
    *im = mk_int(0);
}

/* Apply an elementary-function formula given the decomposed argument (u, v). */
static void cx_apply_elem(const char* head, Expr* u, Expr* v, Expr** re, Expr** im) {
    /* Circular */
    if (strcmp(head, "Sin") == 0) {
        *re = ev(mk_fn2("Times", mk_fn1("Sin", expr_copy(u)), mk_fn1("Cosh", expr_copy(v))));
        *im = ev(mk_fn2("Times", mk_fn1("Cos", expr_copy(u)), mk_fn1("Sinh", expr_copy(v))));
        return;
    }
    if (strcmp(head, "Cos") == 0) {
        *re = ev(mk_fn2("Times", mk_fn1("Cos", expr_copy(u)), mk_fn1("Cosh", expr_copy(v))));
        *im = ev(mk_fn2("Times", mk_fn2("Times", mk_int(-1), mk_fn1("Sin", expr_copy(u))),
                         mk_fn1("Sinh", expr_copy(v))));
        return;
    }
    /* Hyperbolic */
    if (strcmp(head, "Sinh") == 0) {
        *re = ev(mk_fn2("Times", mk_fn1("Sinh", expr_copy(u)), mk_fn1("Cos", expr_copy(v))));
        *im = ev(mk_fn2("Times", mk_fn1("Cosh", expr_copy(u)), mk_fn1("Sin", expr_copy(v))));
        return;
    }
    if (strcmp(head, "Cosh") == 0) {
        *re = ev(mk_fn2("Times", mk_fn1("Cosh", expr_copy(u)), mk_fn1("Cos", expr_copy(v))));
        *im = ev(mk_fn2("Times", mk_fn1("Sinh", expr_copy(u)), mk_fn1("Sin", expr_copy(v))));
        return;
    }
    /* Tan / Cot / Tanh / Coth via double-angle rational forms. */
    if (strcmp(head, "Tan") == 0 || strcmp(head, "Cot") == 0) {
        Expr* two_u = ev(mk_fn2("Times", mk_int(2), expr_copy(u)));
        Expr* two_v = ev(mk_fn2("Times", mk_int(2), expr_copy(v)));
        Expr* den = ev(mk_fn2("Plus", mk_fn1("Cos", expr_copy(two_u)), mk_fn1("Cosh", expr_copy(two_v))));
        Expr* s2u = mk_fn1("Sin", expr_copy(two_u));
        Expr* sh2v = mk_fn1("Sinh", expr_copy(two_v));
        if (strcmp(head, "Tan") == 0) {
            *re = ev(mk_fn2("Times", s2u, mk_pow(expr_copy(den), mk_int(-1))));
            *im = ev(mk_fn2("Times", sh2v, mk_pow(expr_copy(den), mk_int(-1))));
        } else { /* Cot: (Sin2u - I Sinh2v)/(Cosh2v - Cos2u) */
            Expr* den2 = ev(mk_fn2("Plus", mk_fn1("Cosh", expr_copy(two_v)),
                                   mk_fn2("Times", mk_int(-1), mk_fn1("Cos", expr_copy(two_u)))));
            *re = ev(mk_fn2("Times", s2u, mk_pow(expr_copy(den2), mk_int(-1))));
            *im = ev(mk_fn2("Times", mk_fn2("Times", mk_int(-1), sh2v),
                            mk_pow(expr_copy(den2), mk_int(-1))));
            expr_free(den2);
        }
        expr_free(two_u); expr_free(two_v); expr_free(den);
        return;
    }
    if (strcmp(head, "Tanh") == 0 || strcmp(head, "Coth") == 0) {
        Expr* two_u = ev(mk_fn2("Times", mk_int(2), expr_copy(u)));
        Expr* two_v = ev(mk_fn2("Times", mk_int(2), expr_copy(v)));
        Expr* den = ev(mk_fn2("Plus", mk_fn1("Cosh", expr_copy(two_u)), mk_fn1("Cos", expr_copy(two_v))));
        Expr* sh2u = mk_fn1("Sinh", expr_copy(two_u));
        Expr* s2v = mk_fn1("Sin", expr_copy(two_v));
        if (strcmp(head, "Tanh") == 0) {
            *re = ev(mk_fn2("Times", sh2u, mk_pow(expr_copy(den), mk_int(-1))));
            *im = ev(mk_fn2("Times", s2v, mk_pow(expr_copy(den), mk_int(-1))));
        } else { /* Coth: (Sinh2u - I Sin2v)/(Cosh2u - Cos2v) */
            Expr* den2 = ev(mk_fn2("Plus", mk_fn1("Cosh", expr_copy(two_u)),
                                   mk_fn2("Times", mk_int(-1), mk_fn1("Cos", expr_copy(two_v)))));
            *re = ev(mk_fn2("Times", sh2u, mk_pow(expr_copy(den2), mk_int(-1))));
            *im = ev(mk_fn2("Times", mk_fn2("Times", mk_int(-1), s2v),
                            mk_pow(expr_copy(den2), mk_int(-1))));
            expr_free(den2);
        }
        expr_free(two_u); expr_free(two_v); expr_free(den);
        return;
    }
    /* Sec / Csc / Sech / Csch via reciprocal of Cos / Sin / Cosh / Sinh. */
    if (strcmp(head, "Sec") == 0 || strcmp(head, "Csc") == 0 ||
        strcmp(head, "Sech") == 0 || strcmp(head, "Csch") == 0) {
        const char* base = strcmp(head, "Sec") == 0 ? "Cos"
                         : strcmp(head, "Csc") == 0 ? "Sin"
                         : strcmp(head, "Sech") == 0 ? "Cosh" : "Sinh";
        Expr *br, *bi;
        cx_apply_elem(base, u, v, &br, &bi);
        cx_crecip(br, bi, re, im);
        expr_free(br); expr_free(bi);
        return;
    }
    /* Unreachable. */
    *re = mk_int(0); *im = mk_int(0);
}

/* Rewrite an inverse trig / hyperbolic head into logarithmic form, given the
 * (already-copied) argument expression `w`.  Returns a fresh Expr (owns w). */
static Expr* cx_inverse_to_log(const char* head, Expr* w) {
    #define W       (expr_copy(w))
    #define NEG(e)  (mk_fn2("Times", mk_int(-1), (e)))
    #define SQRT(e) (mk_pow((e), mk_half()))
    #define RECIP(e)(mk_pow((e), mk_int(-1)))
    Expr* out = NULL;
    if (strcmp(head, "ArcSin") == 0) {
        /* -I Log[I w + Sqrt[1 - w^2]] */
        Expr* inner = mk_fn2("Plus", mk_fn2("Times", mk_I(), W),
                             SQRT(mk_fn2("Plus", mk_int(1), NEG(mk_pow(W, mk_int(2))))));
        out = mk_fn2("Times", NEG(mk_I()), mk_fn1("Log", inner));
    } else if (strcmp(head, "ArcCos") == 0) {
        /* Pi/2 + I Log[I w + Sqrt[1 - w^2]] */
        Expr* inner = mk_fn2("Plus", mk_fn2("Times", mk_I(), W),
                             SQRT(mk_fn2("Plus", mk_int(1), NEG(mk_pow(W, mk_int(2))))));
        out = mk_fn2("Plus", mk_fn2("Times", mk_half(), mk_sym("Pi")),
                     mk_fn2("Times", mk_I(), mk_fn1("Log", inner)));
    } else if (strcmp(head, "ArcTan") == 0) {
        /* (I/2)(Log[1 - I w] - Log[1 + I w]) */
        Expr* la = mk_fn1("Log", mk_fn2("Plus", mk_int(1), NEG(mk_fn2("Times", mk_I(), W))));
        Expr* lb = mk_fn1("Log", mk_fn2("Plus", mk_int(1), mk_fn2("Times", mk_I(), W)));
        out = mk_fn2("Times", mk_fn2("Times", mk_I(), mk_half()),
                     mk_fn2("Plus", la, NEG(lb)));
    } else if (strcmp(head, "ArcCot") == 0) {
        /* (I/2)(Log[1 - I/w] - Log[1 + I/w]) */
        Expr* la = mk_fn1("Log", mk_fn2("Plus", mk_int(1), NEG(mk_fn2("Times", mk_I(), RECIP(W)))));
        Expr* lb = mk_fn1("Log", mk_fn2("Plus", mk_int(1), mk_fn2("Times", mk_I(), RECIP(W))));
        out = mk_fn2("Times", mk_fn2("Times", mk_I(), mk_half()),
                     mk_fn2("Plus", la, NEG(lb)));
    } else if (strcmp(head, "ArcCsc") == 0) {
        /* ArcSin[1/w] = -I Log[I/w + Sqrt[1 - 1/w^2]] */
        Expr* inner = mk_fn2("Plus", mk_fn2("Times", mk_I(), RECIP(W)),
                             SQRT(mk_fn2("Plus", mk_int(1), NEG(mk_pow(RECIP(W), mk_int(2))))));
        out = mk_fn2("Times", NEG(mk_I()), mk_fn1("Log", inner));
    } else if (strcmp(head, "ArcSec") == 0) {
        /* ArcCos[1/w] = Pi/2 + I Log[I/w + Sqrt[1 - 1/w^2]] */
        Expr* inner = mk_fn2("Plus", mk_fn2("Times", mk_I(), RECIP(W)),
                             SQRT(mk_fn2("Plus", mk_int(1), NEG(mk_pow(RECIP(W), mk_int(2))))));
        out = mk_fn2("Plus", mk_fn2("Times", mk_half(), mk_sym("Pi")),
                     mk_fn2("Times", mk_I(), mk_fn1("Log", inner)));
    } else if (strcmp(head, "ArcSinh") == 0) {
        /* Log[w + Sqrt[1 + w^2]] */
        out = mk_fn1("Log", mk_fn2("Plus", W, SQRT(mk_fn2("Plus", mk_int(1), mk_pow(W, mk_int(2))))));
    } else if (strcmp(head, "ArcCosh") == 0) {
        /* Log[w + Sqrt[w - 1] Sqrt[w + 1]] */
        Expr* rad = mk_fn2("Times", SQRT(mk_fn2("Plus", W, mk_int(-1))),
                           SQRT(mk_fn2("Plus", W, mk_int(1))));
        out = mk_fn1("Log", mk_fn2("Plus", W, rad));
    } else if (strcmp(head, "ArcCsch") == 0) {
        /* ArcSinh[1/w] = Log[1/w + Sqrt[1 + 1/w^2]] */
        out = mk_fn1("Log", mk_fn2("Plus", RECIP(W),
                     SQRT(mk_fn2("Plus", mk_int(1), mk_pow(RECIP(W), mk_int(2))))));
    } else if (strcmp(head, "ArcTanh") == 0) {
        /* (1/2)(Log[1 + w] - Log[1 - w]) */
        Expr* la = mk_fn1("Log", mk_fn2("Plus", mk_int(1), W));
        Expr* lb = mk_fn1("Log", mk_fn2("Plus", mk_int(1), NEG(W)));
        out = mk_fn2("Times", mk_half(), mk_fn2("Plus", la, NEG(lb)));
    } else if (strcmp(head, "ArcCoth") == 0) {
        /* (1/2)(Log[1 + 1/w] - Log[1 - 1/w]) */
        Expr* la = mk_fn1("Log", mk_fn2("Plus", mk_int(1), RECIP(W)));
        Expr* lb = mk_fn1("Log", mk_fn2("Plus", mk_int(1), NEG(RECIP(W))));
        out = mk_fn2("Times", mk_half(), mk_fn2("Plus", la, NEG(lb)));
    } else if (strcmp(head, "ArcSech") == 0) {
        /* ArcCosh[1/w] = Log[1/w + Sqrt[1/w - 1] Sqrt[1/w + 1]] */
        Expr* rad = mk_fn2("Times", SQRT(mk_fn2("Plus", RECIP(W), mk_int(-1))),
                           SQRT(mk_fn2("Plus", RECIP(W), mk_int(1))));
        out = mk_fn1("Log", mk_fn2("Plus", RECIP(W), rad));
    }
    #undef W
    #undef NEG
    #undef SQRT
    #undef RECIP
    expr_free(w);
    return out;
}

static bool cx_is_inverse_head(const char* s) {
    static const char* names[] = {
        "ArcSin","ArcCos","ArcTan","ArcCot","ArcSec","ArcCsc",
        "ArcSinh","ArcCosh","ArcTanh","ArcCoth","ArcSech","ArcCsch"
    };
    for (size_t i = 0; i < sizeof(names)/sizeof(names[0]); i++)
        if (strcmp(s, names[i]) == 0) return true;
    return false;
}

static bool cx_is_elem_head(const char* s) {
    static const char* names[] = {
        "Sin","Cos","Tan","Cot","Sec","Csc",
        "Sinh","Cosh","Tanh","Coth","Sech","Csch"
    };
    for (size_t i = 0; i < sizeof(names)/sizeof(names[0]); i++)
        if (strcmp(s, names[i]) == 0) return true;
    return false;
}

/* Power[b, e] decomposition. */
static void cx_power(const Expr* base, const Expr* exp, const CxCtx* ctx,
                     Expr** re, Expr** im) {
    Expr *bu, *bv, *eu, *ev_;
    cx_decompose(base, ctx, &bu, &bv);
    cx_decompose(exp, ctx, &eu, &ev_);

    bool base_real = cx_is_zero(bv);
    bool exp_real  = cx_is_zero(ev_);

    /* E^z  (base symbol E) -> Exp decomposition. */
    if (base->type == EXPR_SYMBOL && base->data.symbol.name == SYM_E) {
        Expr* mod = ev(mk_fn1("Exp", expr_copy(eu)));           /* Exp[Re] */
        *re = ev(mk_fn2("Times", expr_copy(mod), mk_fn1("Cos", expr_copy(ev_))));
        *im = ev(mk_fn2("Times", mod, mk_fn1("Sin", expr_copy(ev_))));
        goto done;
    }

    /* Real base, real integer exponent -> real. */
    long n;
    bool exp_int = exp_real && cx_is_int_literal(eu, &n);
    if (base_real && exp_int) {
        *re = ev(mk_pow(expr_copy(bu), mk_int(n)));
        *im = mk_int(0);
        goto done;
    }

    /* Complex base, integer exponent -> Cartesian by repeated complex
     * multiplication (square-and-multiply).  We deliberately do NOT rebuild
     * and re-decompose the value: bu/bv may already contain decomposed atoms
     * (Abs[z], Cos[Arg z], ...) whose re-decomposition is not idempotent. */
    if (exp_int) {
        long m = n < 0 ? -n : n;
        Expr* rr = mk_int(1);           /* accumulator = 1 + 0 I */
        Expr* ri = mk_int(0);
        Expr* pr = expr_copy(bu);       /* running power of the base */
        Expr* pi = expr_copy(bv);
        while (m > 0) {
            if (m & 1) {
                Expr *nr, *ni;
                cx_cmul(rr, ri, pr, pi, &nr, &ni);
                expr_free(rr); expr_free(ri); rr = nr; ri = ni;
            }
            m >>= 1;
            if (m > 0) {
                Expr *sr, *si;
                cx_cmul(pr, pi, pr, pi, &sr, &si);
                expr_free(pr); expr_free(pi); pr = sr; pi = si;
            }
        }
        expr_free(pr); expr_free(pi);
        if (n < 0) {
            Expr *ir, *ii;
            cx_crecip(rr, ri, &ir, &ii);
            expr_free(rr); expr_free(ri); rr = ir; ri = ii;
        }
        *re = rr; *im = ri;
        goto done;
    }

    /* Real positive base, real exponent -> real. */
    if (base_real && exp_real && cx_is_pos_numeric(bu)) {
        *re = ev(mk_pow(expr_copy(bu), expr_copy(eu)));
        *im = mk_int(0);
        goto done;
    }

    /* General polar master formula.
     *   r = Abs[base] = Sqrt[bu^2 + bv^2],   theta = Arg[base]
     *   modulus = r^eu * Exp[-ev*theta]
     *   phase   = eu*theta + ev*Log[r]
     *   re = modulus Cos[phase], im = modulus Sin[phase]
     * theta uses the *original* base so a bare complex atom stays Arg[z]. */
    {
        Expr* ss = cx_sumsq(bu, bv);
        Expr* r  = ev(mk_pow(ss, mk_half()));                    /* Sqrt[bu^2+bv^2] */
        Expr* logr = mk_fn1("Log", expr_copy(r));
        Expr* theta = ev(mk_fn1("Arg", expr_copy((Expr*)base)));

        Expr* modulus = ev(mk_fn2("Times",
                             mk_pow(expr_copy(r), expr_copy(eu)),
                             mk_fn1("Exp", mk_fn2("Times", mk_int(-1),
                                     mk_fn2("Times", expr_copy(ev_), expr_copy(theta))))));
        Expr* phase = ev(mk_fn2("Plus",
                          mk_fn2("Times", expr_copy(eu), expr_copy(theta)),
                          mk_fn2("Times", expr_copy(ev_), logr)));
        *re = ev(mk_fn2("Times", expr_copy(modulus), mk_fn1("Cos", expr_copy(phase))));
        *im = ev(mk_fn2("Times", modulus, mk_fn1("Sin", phase)));
        expr_free(r); expr_free(theta);
    }

done:
    expr_free(bu); expr_free(bv); expr_free(eu); expr_free(ev_);
}

/* Core recursion. */
static void cx_decompose(const Expr* e, const CxCtx* ctx, Expr** re, Expr** im) {
    if (!e) { *re = mk_int(0); *im = mk_int(0); return; }

    /* Complex[a, b] literal (covers I = Complex[0,1]). */
    Expr *cre, *cim;
    if (is_complex((Expr*)e, &cre, &cim)) {
        *re = expr_copy(cre);
        *im = expr_copy(cim);
        return;
    }

    if (e->type != EXPR_FUNCTION) { cx_atom(e, ctx, re, im); return; }

    const Expr* head = e->data.function.head;
    if (head->type != EXPR_SYMBOL) { cx_atom(e, ctx, re, im); return; }
    const char* h = head->data.symbol.name;
    Expr** args = e->data.function.args;
    size_t nargs = e->data.function.arg_count;

    /* Plus: componentwise sum. */
    if (h == SYM_Plus) {
        Expr* sr = mk_int(0);
        Expr* si = mk_int(0);
        for (size_t k = 0; k < nargs; k++) {
            Expr *tr, *ti;
            cx_decompose(args[k], ctx, &tr, &ti);
            Expr* nr = cx_add2(sr, tr);
            Expr* ni = cx_add2(si, ti);
            expr_free(sr); expr_free(si); expr_free(tr); expr_free(ti);
            sr = nr; si = ni;
        }
        *re = sr; *im = si;
        return;
    }

    /* Times: complex-multiply fold. */
    if (h == SYM_Times) {
        Expr* ar = mk_int(1);
        Expr* ai = mk_int(0);
        for (size_t k = 0; k < nargs; k++) {
            Expr *fr, *fi;
            cx_decompose(args[k], ctx, &fr, &fi);
            Expr *nr, *ni;
            cx_cmul(ar, ai, fr, fi, &nr, &ni);
            expr_free(ar); expr_free(ai); expr_free(fr); expr_free(fi);
            ar = nr; ai = ni;
        }
        *re = ar; *im = ai;
        return;
    }

    /* Power. */
    if (h == SYM_Power && nargs == 2) {
        cx_power(args[0], args[1], ctx, re, im);
        return;
    }

    /* Exp[z] = Power[E, z] normally, but handle a literal Exp head too. */
    if (strcmp(h, "Exp") == 0 && nargs == 1) {
        Expr *u, *v;
        cx_decompose(args[0], ctx, &u, &v);
        Expr* mod = ev(mk_fn1("Exp", expr_copy(u)));
        *re = ev(mk_fn2("Times", expr_copy(mod), mk_fn1("Cos", expr_copy(v))));
        *im = ev(mk_fn2("Times", mod, mk_fn1("Sin", expr_copy(v))));
        expr_free(u); expr_free(v);
        return;
    }

    /* Log[w] = (1/2) Log[Expand[u^2+v^2]] + I Arg[w]. */
    if (h == SYM_Log && nargs == 1) {
        Expr *u, *v;
        cx_decompose(args[0], ctx, &u, &v);
        if (cx_is_zero(v) && cx_is_pos_numeric(u)) {
            *re = mk_fn1("Log", expr_copy(u));
            *im = mk_int(0);
        } else {
            Expr* ss = cx_sumsq(u, v);
            *re = ev(mk_fn2("Times", mk_half(), mk_fn1("Log", ss)));
            *im = ev(mk_fn1("Arg", expr_copy(args[0])));
        }
        expr_free(u); expr_free(v);
        return;
    }

    /* Circular / hyperbolic. */
    if (nargs == 1 && cx_is_elem_head(h)) {
        Expr *u, *v;
        cx_decompose(args[0], ctx, &u, &v);
        cx_apply_elem(h, u, v, re, im);
        expr_free(u); expr_free(v);
        return;
    }

    /* Inverse circular / hyperbolic: rewrite to Log form and recurse. */
    if (nargs == 1 && cx_is_inverse_head(h)) {
        Expr* logform = cx_inverse_to_log(h, expr_copy(args[0]));
        if (logform) {
            cx_decompose(logform, ctx, re, im);
            expr_free(logform);
            return;
        }
    }

    /* Re[w] / Im[w]: real and imaginary components of w. */
    if (h == SYM_Re && nargs == 1) {
        Expr *u, *v;
        cx_decompose(args[0], ctx, &u, &v);
        *re = u; *im = mk_int(0);
        expr_free(v);
        return;
    }
    if (h == SYM_Im && nargs == 1) {
        Expr *u, *v;
        cx_decompose(args[0], ctx, &u, &v);
        *re = v; *im = mk_int(0);
        expr_free(u);
        return;
    }

    /* Abs[w] = Sqrt[u^2 + v^2]  (real-valued). */
    if (h == SYM_Abs && nargs == 1) {
        Expr *u, *v;
        cx_decompose(args[0], ctx, &u, &v);
        if (cx_is_zero(v)) {
            *re = ev(mk_fn1("Abs", expr_copy(u)));
        } else {
            Expr* ss = cx_sumsq(u, v);
            *re = ev(mk_pow(ss, mk_half()));
        }
        *im = mk_int(0);
        expr_free(u); expr_free(v);
        return;
    }

    /* Arg[w]: kept symbolic on the original expression (real-valued). */
    if (h == SYM_Arg && nargs == 1) {
        *re = ev(mk_fn1("Arg", expr_copy(args[0])));
        *im = mk_int(0);
        return;
    }

    /* Conjugate[w] = u - I v. */
    if (h == SYM_Conjugate && nargs == 1) {
        Expr *u, *v;
        cx_decompose(args[0], ctx, &u, &v);
        *re = u;
        *im = ev(mk_fn2("Times", mk_int(-1), v));
        return;
    }

    /* Sign[w]: real -> Sign[u]; complex -> w / Abs[w]. */
    if (h == SYM_Sign && nargs == 1) {
        Expr *u, *v;
        cx_decompose(args[0], ctx, &u, &v);
        if (cx_is_zero(v)) {
            *re = ev(mk_fn1("Sign", expr_copy(u)));
            *im = mk_int(0);
        } else {
            Expr* ss = cx_sumsq(u, v);
            Expr* absv = ev(mk_pow(ss, mk_half()));
            Expr* zero = mk_int(0);
            cx_cdiv(u, v, absv, zero, re, im);  /* cx_cdiv copies its args */
            expr_free(absv); expr_free(zero);
        }
        expr_free(u); expr_free(v);
        return;
    }

    /* Unknown head. */
    {
        /* Decompose every argument; if all are real, assume f is real-valued. */
        bool all_real = true;
        Expr** recon = malloc(sizeof(Expr*) * (nargs ? nargs : 1));
        for (size_t k = 0; k < nargs; k++) {
            Expr *u, *v;
            cx_decompose(args[k], ctx, &u, &v);
            if (!cx_is_zero(v)) all_real = false;
            recon[k] = cx_recon(u, v);
            expr_free(u); expr_free(v);
        }
        Expr* rebuilt = ev(mk_fnN(h, recon, nargs)); /* adopts recon[] */
        free(recon);
        if (all_real) {
            *re = rebuilt;
            *im = mk_int(0);
        } else {
            *re = mk_fn1("Re", expr_copy(rebuilt));
            *im = mk_fn1("Im", rebuilt);
        }
        return;
    }
}

/* ------------------------------------------------------------------ *
 *  TargetFunctions -> Conjugate path.                                 *
 * ------------------------------------------------------------------ */

/* Conjugate expression: I -> -I (Complex[a,b] -> Complex[a,-b]) and each
 * listed complex variable z -> Conjugate[z]; real symbols unchanged. */
static Expr* cx_conjugate_of(const Expr* e, const CxCtx* ctx) {
    Expr *cre, *cim;
    if (is_complex((Expr*)e, &cre, &cim)) {
        return make_complex(expr_copy(cre), ev(mk_fn2("Times", mk_int(-1), expr_copy(cim))));
    }
    if (e->type == EXPR_SYMBOL) {
        if (cx_is_complex_atom(e, ctx)) return mk_fn1("Conjugate", expr_copy((Expr*)e));
        return expr_copy((Expr*)e);
    }
    if (e->type == EXPR_FUNCTION) {
        size_t n = e->data.function.arg_count;
        Expr* head = cx_conjugate_of(e->data.function.head, ctx);
        Expr** a = malloc(sizeof(Expr*) * (n ? n : 1));
        for (size_t k = 0; k < n; k++)
            a[k] = cx_conjugate_of(e->data.function.args[k], ctx);
        Expr* out = expr_new_function(head, a, n);
        free(a);
        return out;
    }
    return expr_copy((Expr*)e);
}

/* ------------------------------------------------------------------ *
 *  Front-end helpers.                                                 *
 * ------------------------------------------------------------------ */

/* True iff `e` is f[...] with head symbol named exactly `name`. */
static bool head_name_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           strcmp(e->data.function.head->data.symbol.name, name) == 0;
}

/* Heads over which ComplexExpand threads (besides List, handled separately). */
static bool cx_threads_over(const char* s) {
    static const char* names[] = {
        "Equal","Unequal","Less","LessEqual","Greater","GreaterEqual",
        "SameQ","UnsameQ","And","Or","Not","Xor","Implies","Rule","RuleDelayed"
    };
    for (size_t i = 0; i < sizeof(names)/sizeof(names[0]); i++)
        if (strcmp(s, names[i]) == 0) return true;
    return false;
}

/* Assemble the final result for a single (non-threaded) expression. */
static Expr* cx_expand_one(const Expr* arg, const CxCtx* ctx);

/* Whole-expression decomposition dispatch (handles ReIm and Conjugate mode). */
static Expr* cx_expand_one(const Expr* arg, const CxCtx* ctx) {
    /* ReIm[w] -> {Re-expansion, Im-expansion}. */
    if (head_name_is(arg, "ReIm") && arg->data.function.arg_count == 1) {
        const Expr* inner = arg->data.function.args[0];
        Expr* rr = mk_fn1("Re", expr_copy((Expr*)inner));
        Expr* ii = mk_fn1("Im", expr_copy((Expr*)inner));
        Expr* er = cx_expand_one(rr, ctx);
        Expr* ei = cx_expand_one(ii, ctx);
        expr_free(rr); expr_free(ii);
        return ev(mk_fn2("List", er, ei));
    }

    if (ctx->target == TF_CONJUGATE) {
        /* Peel a single real-valued / conjugation wrapper so we split the
         * inner expression, not the (already-real) wrapper. */
        const char* wrap = NULL;
        const Expr* inner = arg;
        if (arg->type == EXPR_FUNCTION && arg->data.function.arg_count == 1) {
            const char* wh = head_is(arg, SYM_Re) ? "Re"
                           : head_is(arg, SYM_Im) ? "Im"
                           : head_is(arg, SYM_Abs) ? "Abs"
                           : head_is(arg, SYM_Arg) ? "Arg"
                           : head_is(arg, SYM_Conjugate) ? "Conjugate" : NULL;
            if (wh) { wrap = wh; inner = arg->data.function.args[0]; }
        }
        Expr* wc = cx_conjugate_of(inner, ctx);
        /* re = (inner + wc)/2, im = (inner - wc)/(2 I) = -I (inner - wc)/2. */
        Expr* re = ev(mk_fn2("Times", mk_half(), mk_fn2("Plus", expr_copy((Expr*)inner), expr_copy(wc))));
        Expr* im = ev(mk_fn2("Times", mk_fn2("Times", mk_int(-1), mk_half()),
                        mk_fn2("Times", mk_I(), mk_fn2("Plus", expr_copy((Expr*)inner),
                               mk_fn2("Times", mk_int(-1), expr_copy(wc))))));
        Expr* out;
        if (!wrap) {
            out = ev(mk_fn2("Plus", re, mk_fn2("Times", mk_I(), im)));
        } else if (strcmp(wrap, "Re") == 0)  { out = re; re = NULL; }
        else if (strcmp(wrap, "Im") == 0)    { out = im; im = NULL; }
        else if (strcmp(wrap, "Arg") == 0)   { out = ev(mk_fn1("Arg", expr_copy((Expr*)inner))); }
        else if (strcmp(wrap, "Conjugate") == 0) { out = expr_copy(wc); }
        else /* Abs */ {
            Expr* ss = cx_sumsq(re, im);
            out = ev(mk_pow(ss, mk_half()));
        }
        if (re) expr_free(re);
        if (im) expr_free(im);
        expr_free(wc);
        return ev(mk_fn1("Simplify", out));
    }

    Expr *re, *im;
    cx_decompose(arg, ctx, &re, &im);
    Expr* combined = mk_fn2("Plus", re, mk_fn2("Times", mk_I(), im));
    return ev(mk_fn1("Expand", combined));
}

/* Recursively thread over List and relation/logic heads, else expand. */
static Expr* cx_expand_thread(const Expr* arg, const CxCtx* ctx) {
    if (arg->type == EXPR_FUNCTION && arg->data.function.head->type == EXPR_SYMBOL) {
        const char* h = arg->data.function.head->data.symbol.name;
        bool is_list = (h == SYM_List);
        if (is_list || cx_threads_over(h)) {
            size_t n = arg->data.function.arg_count;
            Expr** out = malloc(sizeof(Expr*) * (n ? n : 1));
            for (size_t k = 0; k < n; k++)
                out[k] = cx_expand_thread(arg->data.function.args[k], ctx);
            Expr* rebuilt = expr_new_function(expr_copy(arg->data.function.head), out, n);
            free(out);
            Expr* evd = ev(rebuilt);
            /* For relations, try to collapse to True/False. */
            if (!is_list && (strcmp(h, "Equal") == 0 || strcmp(h, "Unequal") == 0)) {
                Expr* simp = ev(mk_fn1("Simplify", expr_copy(evd)));
                if (simp->type == EXPR_SYMBOL &&
                    (simp->data.symbol.name == SYM_True || simp->data.symbol.name == SYM_False)) {
                    expr_free(evd);
                    return simp;
                }
                expr_free(simp);
            }
            return evd;
        }
    }
    return cx_expand_one(arg, ctx);
}

/* ------------------------------------------------------------------ *
 *  Builtin front-end.                                                 *
 * ------------------------------------------------------------------ */

Expr* builtin_complex_expand(Expr* res) {
    size_t argc = res->data.function.arg_count;
    Expr** argv = res->data.function.args;

    /* Split trailing TargetFunctions -> ... option Rules from positionals. */
    TargetMode target = TF_REIM;
    size_t npos = argc;
    for (size_t k = argc; k > 0; k--) {
        Expr* a = argv[k - 1];
        if (a->type == EXPR_FUNCTION && a->data.function.arg_count == 2 &&
            a->data.function.head->type == EXPR_SYMBOL &&
            (a->data.function.head->data.symbol.name == SYM_Rule ||
             a->data.function.head->data.symbol.name == SYM_RuleDelayed)) {
            Expr* lhs = a->data.function.args[0];
            if (lhs->type == EXPR_SYMBOL && strcmp(lhs->data.symbol.name, "TargetFunctions") == 0) {
                Expr* rhs = a->data.function.args[1];
                /* Detect Conjugate / {Abs,Arg}; default Re,Im otherwise. */
                bool has_conj = false, has_abs = false, has_arg = false;
                Expr** items; size_t ni;
                if (head_is(rhs, SYM_List)) { items = rhs->data.function.args; ni = rhs->data.function.arg_count; }
                else { items = &rhs; ni = 1; }
                for (size_t j = 0; j < ni; j++) {
                    if (items[j]->type == EXPR_SYMBOL) {
                        const char* nm = items[j]->data.symbol.name;
                        if (strcmp(nm, "Conjugate") == 0) has_conj = true;
                        else if (strcmp(nm, "Abs") == 0) has_abs = true;
                        else if (strcmp(nm, "Arg") == 0) has_arg = true;
                    }
                }
                if (has_conj) target = TF_CONJUGATE;
                else if (has_abs || has_arg) target = TF_ABSARG;
                else target = TF_REIM;
                npos--;   /* consume this option */
                continue;
            }
        }
        break;  /* first non-option positional from the right: stop. */
    }

    if (npos < 1 || npos > 2) {
        fprintf(stderr, "General::argct: ComplexExpand called with %zu arguments.\n", argc);
        return NULL;
    }

    /* Build complex-variable list from the optional second positional. */
    Expr** cvars = NULL;
    size_t ncvars = 0;
    if (npos == 2) {
        Expr* spec = argv[1];
        if (head_is(spec, SYM_List)) {
            ncvars = spec->data.function.arg_count;
            cvars = ncvars ? malloc(sizeof(Expr*) * ncvars) : NULL;
            for (size_t k = 0; k < ncvars; k++) cvars[k] = spec->data.function.args[k];
        } else {
            ncvars = 1;
            cvars = malloc(sizeof(Expr*));
            cvars[0] = spec;
        }
    }

    CxCtx ctx = { cvars, ncvars, target };
    Expr* out = cx_expand_thread(argv[0], &ctx);
    if (cvars) free(cvars);
    return out;
}

/* ------------------------------------------------------------------ *
 *  Registration.                                                      *
 * ------------------------------------------------------------------ */

void complex_expand_init(void) {
    symtab_add_builtin("ComplexExpand", builtin_complex_expand);
    symtab_get_def("ComplexExpand")->attributes |= ATTR_PROTECTED;
}
