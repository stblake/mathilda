/*
 * expand_power.c -- PowerExpand.
 *
 * PowerExpand distributes powers over products and collapses nested powers
 * and logarithms:
 *
 *     (a b)^c   -> a^c b^c
 *     (a^b)^c   -> a^(b c)
 *     Log[a b]  -> Log[a] + Log[b]
 *     Log[a^b]  -> b Log[a]
 *     Arg[a b]  -> Arg[a] + Arg[b]
 *
 * Sqrt is stored as Power[x, 1/2] and Log[1/z] as Log[Power[z, -1]], so the
 * rules above cover Sqrt[x y] -> Sqrt[x] Sqrt[y], Sqrt[z^2] -> z and
 * Log[1/z] -> -Log[z] without any special-casing.
 *
 * The rewrites are applied top-down to a fixed point (mirroring
 * Mathematica's ReplaceRepeated semantics): a rule fires at the outermost
 * matching node and the transformed result is reprocessed.  This is why, e.g.
 * Log[(a b)^c] becomes c (Log[a] + Log[b]) rather than expanding the inner
 * power first.
 *
 * Three modes, selected by the Assumptions option:
 *
 *   - Automatic (default): the textbook transforms above.  Correct when the
 *     bases are positive reals (or the exponents integers); branch issues are
 *     ignored.
 *   - Assumptions -> True: emit the universally-correct formulas, attaching a
 *     branch-correction term built from Floor / Arg / Im / E / I / Pi.
 *   - Assumptions -> assum: emit the True-mode formula and then refine the
 *     correction terms under the assumptions (Arg / Im of known-sign reals,
 *     Floor over assumption-bounded intervals).  Faithful on the documented
 *     examples; where the assumptions fall outside this reasoning it degrades
 *     gracefully to the symbolic True-mode form rather than a wrong value.
 *
 * PowerExpand threads over List, equations, inequalities and logic functions,
 * and supports the variable-restricted form PowerExpand[expr, {x1, ...}].
 */

#include "expand_power.h"
#include "eval.h"
#include "symtab.h"
#include "sym_names.h"
#include "simp/simp.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ------------------------------------------------------------------ */
/* Context                                                            */
/* ------------------------------------------------------------------ */

typedef enum { PE_AUTOMATIC, PE_TRUE, PE_ASSUME } PEMode;

typedef struct {
    PEMode      mode;
    Expr**      vars;   /* variable restriction set (NULL/0 = unrestricted) */
    size_t      nvars;
    AssumeCtx*  ctx;    /* non-NULL only in PE_ASSUME */
} PECtx;

/* ------------------------------------------------------------------ */
/* Small constructors (each takes ownership of its Expr* arguments)   */
/* ------------------------------------------------------------------ */

static Expr* fn1(const char* head, Expr* a) {
    Expr* args[1] = { a };
    return expr_new_function(expr_new_symbol(head), args, 1);
}
static Expr* fn2(const char* head, Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    return expr_new_function(expr_new_symbol(head), args, 2);
}
static Expr* fnN(const char* head, Expr** args, size_t n) {
    return expr_new_function(expr_new_symbol(head), args, n);
}

static Expr* pe_int(int64_t n)  { return expr_new_integer(n); }
static Expr* pe_pi(void)        { return expr_new_symbol(SYM_Pi); }
static Expr* pe_E(void)         { return expr_new_symbol(SYM_E); }
static Expr* pe_I(void)         { return fn2("Complex", pe_int(0), pe_int(1)); }
static Expr* pe_half(void)      { return fn2("Rational", pe_int(1), pe_int(2)); }
static Expr* pe_recip_2pi(void) { return fn2("Power", fn2("Times", pe_int(2), pe_pi()), pe_int(-1)); }

/* ------------------------------------------------------------------ */
/* Head / structure helpers                                           */
/* ------------------------------------------------------------------ */

static const char* head_name(const Expr* e) {
    if (e->type != EXPR_FUNCTION) return "";
    if (e->data.function.head->type != EXPR_SYMBOL) return "";
    return e->data.function.head->data.symbol;
}

static bool is_head(const Expr* e, const char* h) {
    return e->type == EXPR_FUNCTION && strcmp(head_name(e), h) == 0;
}

static bool is_thread_head(const char* h) {
    return strcmp(h, "List") == 0 || strcmp(h, "Equal") == 0 ||
           strcmp(h, "Unequal") == 0 || strcmp(h, "Less") == 0 ||
           strcmp(h, "LessEqual") == 0 || strcmp(h, "Greater") == 0 ||
           strcmp(h, "GreaterEqual") == 0 || strcmp(h, "Inequality") == 0 ||
           strcmp(h, "And") == 0 || strcmp(h, "Or") == 0 || strcmp(h, "Not") == 0;
}

/* True if e mentions any of the restriction variables (or if unrestricted). */
static bool expr_mentions(const Expr* e, const Expr* v) {
    if (expr_eq((Expr*)e, (Expr*)v)) return true;
    if (e->type == EXPR_FUNCTION) {
        if (expr_mentions(e->data.function.head, v)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (expr_mentions(e->data.function.args[i], v)) return true;
    }
    return false;
}
static bool var_gate(const Expr* src, const PECtx* ctx) {
    if (ctx->nvars == 0) return true;            /* unrestricted */
    for (size_t i = 0; i < ctx->nvars; i++)
        if (expr_mentions(src, ctx->vars[i])) return true;
    return false;
}

/* ------------------------------------------------------------------ */
/* Correction-term builders                                           */
/*                                                                    */
/* All return NULL in Automatic mode (no correction).  In True/Assume */
/* mode they build the full symbolic branch-correction term.          */
/* ------------------------------------------------------------------ */

/* K_prod = Floor[1/2 - (Arg[z1] + ... + Arg[zn]) / (2 Pi)] */
static Expr* k_prod(Expr** zs, size_t n) {
    Expr** argterms = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) argterms[i] = fn1("Arg", expr_copy(zs[i]));
    Expr* argsum = fnN("Plus", argterms, n);
    free(argterms);
    Expr* frac = fn2("Times", pe_recip_2pi(), argsum);            /* sum/(2Pi) */
    Expr* inside = fn2("Plus", pe_half(), fn2("Times", pe_int(-1), frac));
    return fn1("Floor", inside);
}

/* K_pow = Floor[1/2 - Im[x Log[a]] / (2 Pi)] for base a, inner exponent x. */
static Expr* k_pow(Expr* a, Expr* x) {
    Expr* xloga = fn2("Times", expr_copy(x), fn1("Log", expr_copy(a)));
    Expr* im = fn1("Im", xloga);
    Expr* frac = fn2("Times", pe_recip_2pi(), im);
    Expr* inside = fn2("Plus", pe_half(), fn2("Times", pe_int(-1), frac));
    return fn1("Floor", inside);
}

/* (z1...zn)^c correction: E^(2 I c Pi K_prod) */
static Expr* corr_powprod(const PECtx* ctx, Expr** zs, size_t n, Expr* c) {
    if (ctx->mode == PE_AUTOMATIC) return NULL;
    Expr* exps[5] = { pe_int(2), pe_I(), expr_copy(c), pe_pi(), k_prod(zs, n) };
    return fn2("Power", pe_E(), fnN("Times", exps, 5));
}
/* Log[z1...zn] correction: 2 I Pi K_prod */
static Expr* corr_logprod(const PECtx* ctx, Expr** zs, size_t n) {
    if (ctx->mode == PE_AUTOMATIC) return NULL;
    Expr* t[4] = { pe_int(2), pe_I(), pe_pi(), k_prod(zs, n) };
    return fnN("Times", t, 4);
}
/* Arg[z1...zn] correction: 2 Pi K_prod */
static Expr* corr_argprod(const PECtx* ctx, Expr** zs, size_t n) {
    if (ctx->mode == PE_AUTOMATIC) return NULL;
    Expr* t[3] = { pe_int(2), pe_pi(), k_prod(zs, n) };
    return fnN("Times", t, 3);
}
/* (a^x)^y correction: E^(2 I Pi y K_pow) */
static Expr* corr_powpow(const PECtx* ctx, Expr* a, Expr* x, Expr* y) {
    if (ctx->mode == PE_AUTOMATIC) return NULL;
    Expr* exps[5] = { pe_int(2), pe_I(), pe_pi(), expr_copy(y), k_pow(a, x) };
    return fn2("Power", pe_E(), fnN("Times", exps, 5));
}
/* Log[a^x] correction: 2 I Pi K_pow */
static Expr* corr_logpow(const PECtx* ctx, Expr* a, Expr* x) {
    if (ctx->mode == PE_AUTOMATIC) return NULL;
    Expr* t[4] = { pe_int(2), pe_I(), pe_pi(), k_pow(a, x) };
    return fnN("Times", t, 4);
}

/* ------------------------------------------------------------------ */
/* Inverse-trig                                                        */
/* ------------------------------------------------------------------ */

/* Returns the direct-function head that head h inverts, or NULL. */
static const char* inverse_of(const char* h) {
    static const char* pairs[][2] = {
        {"ArcSin","Sin"}, {"ArcCos","Cos"}, {"ArcTan","Tan"}, {"ArcCot","Cot"},
        {"ArcSec","Sec"}, {"ArcCsc","Csc"},
        {"ArcSinh","Sinh"}, {"ArcCosh","Cosh"}, {"ArcTanh","Tanh"},
        {"ArcCoth","Coth"}, {"ArcSech","Sech"}, {"ArcCsch","Csch"},
    };
    for (size_t i = 0; i < sizeof(pairs)/sizeof(pairs[0]); i++)
        if (strcmp(h, pairs[i][0]) == 0) return pairs[i][1];
    return NULL;
}

/* The universally-correct ArcTan[Tan[x]] formula:
 *   x - Pi Floor[1/2 + Re[x]/Pi]
 *     + 1/2 (1 + E^(I Pi (Floor[1/2 - Re[x]/Pi] + Floor[-1/2 + Re[x]/Pi])))
 *       Pi UnitStep[Im[x]]                                                  */
static Expr* arctan_tan_universal(Expr* x) {
    Expr* rex_over_pi = fn2("Times", fn1("Re", expr_copy(x)), fn2("Power", pe_pi(), pe_int(-1)));

    Expr* t1 = expr_copy(x);

    /* -Pi Floor[1/2 + Re[x]/Pi] */
    Expr* f2 = fn1("Floor", fn2("Plus", pe_half(), expr_copy(rex_over_pi)));
    Expr* t2 = fnN("Times", (Expr*[]){ pe_int(-1), pe_pi(), f2 }, 3);

    /* Floor[1/2 - Re[x]/Pi] + Floor[-1/2 + Re[x]/Pi] */
    Expr* fa = fn1("Floor", fn2("Plus", pe_half(),
                    fn2("Times", pe_int(-1), expr_copy(rex_over_pi))));
    Expr* fb = fn1("Floor", fn2("Plus", fn2("Rational", pe_int(-1), pe_int(2)),
                    expr_copy(rex_over_pi)));
    Expr* expo = fnN("Times", (Expr*[]){ pe_I(), pe_pi(), fn2("Plus", fa, fb) }, 3);
    Expr* ph = fn2("Power", pe_E(), expo);
    Expr* bracket = fn2("Plus", pe_int(1), ph);
    Expr* t3 = fnN("Times", (Expr*[]){ pe_half(), bracket, pe_pi(),
                    fn1("UnitStep", fn1("Im", expr_copy(x))) }, 4);

    expr_free(rex_over_pi);
    return fnN("Plus", (Expr*[]){ t1, t2, t3 }, 3);
}

/* ------------------------------------------------------------------ */
/* Rule application (one rewrite at the current node)                  */
/* ------------------------------------------------------------------ */

static bool pe_to_double(const Expr* e, double* out);   /* defined below */

/* True when c is a numeric, non-integer exponent (i.e. a genuine root such as
 * 1/2 or 1/3) -- the case where a negative base would otherwise turn imaginary. */
static bool pe_exp_nonint(const Expr* c) {
    double d;
    return pe_to_double(c, &d) && d != floor(d);
}

/* Negate a factor, distributing -1 through a leading Plus so that, e.g.,
 * -(-1 + u) becomes 1 - u rather than the unevaluated Times[-1, -1 + u]. */
static Expr* pe_negate(const Expr* z) {
    if (is_head(z, "Plus")) {
        size_t k = z->data.function.arg_count;
        Expr** t = malloc(sizeof(Expr*) * k);
        for (size_t i = 0; i < k; i++)
            t[i] = fn2("Times", pe_int(-1), expr_copy(z->data.function.args[i]));
        Expr* sum = fnN("Plus", t, k);
        free(t);
        return eval_and_free(sum);
    }
    return eval_and_free(fn2("Times", pe_int(-1), expr_copy((Expr*)z)));
}

/* Build Times[powers...] (* corr if non-NULL) for (z1..zn)^c.
 *
 * When c is a non-integer numeric exponent (a root) and the product carries a
 * negative numeric coefficient together with a Plus factor, fold the sign of
 * the coefficient into that Plus factor: (-k ... (p))^c -> k^c ... (-p)^c.
 * This keeps the root real where possible, turning e.g.
 *   Sqrt[-4 Dt[u]^2 (-1 + u)]  ->  2 Dt[u] Sqrt[1 - u]
 * instead of 2 I Dt[u] Sqrt[-1 + u].  The rewrite is value-preserving under
 * PowerExpand's branch-cut-ignoring semantics, since (-1)^c (p)^c = (-p)^c. */
static Expr* split_power_product(Expr** zs, size_t n, Expr* c, Expr* corr) {
    size_t neg_i = 0, plus_i = 0;
    bool has_neg = false, has_plus = false;
    if (pe_exp_nonint(c)) {
        for (size_t i = 0; i < n; i++) {
            double d;
            if (!has_neg && pe_to_double(zs[i], &d) && d < 0) { neg_i = i; has_neg = true; }
            if (!has_plus && is_head(zs[i], "Plus"))          { plus_i = i; has_plus = true; }
        }
    }
    bool fold = has_neg && has_plus;

    size_t m = n + (corr ? 1 : 0);
    Expr** factors = malloc(sizeof(Expr*) * m);
    for (size_t i = 0; i < n; i++) {
        Expr* z = (fold && (i == neg_i || i == plus_i))
                      ? pe_negate(zs[i])            /* -k -> k, p -> -p */
                      : expr_copy(zs[i]);
        factors[i] = fn2("Power", z, expr_copy(c));
    }
    if (corr) factors[n] = corr;
    Expr* res = fnN("Times", factors, m);
    free(factors);
    return res;
}

/* Build Plus[heads[zi]...] (+ corr if non-NULL) for Log/Arg of a product. */
static Expr* split_additive_product(const char* h, Expr** zs, size_t n, Expr* corr) {
    size_t m = n + (corr ? 1 : 0);
    Expr** terms = malloc(sizeof(Expr*) * m);
    for (size_t i = 0; i < n; i++) terms[i] = fn1(h, expr_copy(zs[i]));
    if (corr) terms[n] = corr;
    Expr* res = fnN("Plus", terms, m);
    free(terms);
    return res;
}

/*
 * Attempt a single rewrite at e. On success returns true and sets:
 *   *main      the rewritten expression that should be reprocessed (it may
 *              contain further-expandable structure, e.g. Log of a product);
 *   *corr      an inert branch-correction term (NULL in Automatic mode) that
 *              must NOT be reprocessed -- it carries Floor/Arg/Im sub-terms,
 *              including 1/(2 Pi) = Power[Times[2,Pi],-1], that would otherwise
 *              re-trigger the product-power rule forever;
 *   *additive  true when corr combines with Plus (Log/Arg), false for Times.
 */
static bool pe_apply_rule(Expr* e, const PECtx* ctx,
                          Expr** main, Expr** corr, bool* additive) {
    *main = NULL; *corr = NULL; *additive = false;
    if (e->type != EXPR_FUNCTION) return false;
    const char* h = head_name(e);
    Expr** a = e->data.function.args;
    size_t ac = e->data.function.arg_count;

    /* ---- Power ---- */
    if (strcmp(h, "Power") == 0 && ac == 2) {
        Expr* base = a[0];
        Expr* exp  = a[1];
        if (is_head(base, "Times")) {                 /* (z1..zn)^c */
            if (!var_gate(base, ctx)) return false;
            Expr** zs = base->data.function.args;
            size_t n  = base->data.function.arg_count;
            *main = split_power_product(zs, n, exp, NULL);
            *corr = corr_powprod(ctx, zs, n, exp);
            return true;
        }
        if (is_head(base, "Power") && base->data.function.arg_count == 2) {
            if (!var_gate(base, ctx)) return false;    /* (a^b)^c */
            Expr* ib = base->data.function.args[0];
            Expr* ie = base->data.function.args[1];
            *main = fn2("Power", expr_copy(ib),
                        fn2("Times", expr_copy(ie), expr_copy(exp)));
            *corr = corr_powpow(ctx, ib, ie, exp);
            return true;
        }
        return false;
    }

    /* ---- Log ---- */
    if (strcmp(h, "Log") == 0 && ac == 1) {
        Expr* arg = a[0];
        if (is_head(arg, "Power") && arg->data.function.arg_count == 2) {
            if (!var_gate(arg, ctx)) return false;     /* Log[a^b] */
            Expr* base = arg->data.function.args[0];
            Expr* exp  = arg->data.function.args[1];
            *main = fn2("Times", expr_copy(exp), fn1("Log", expr_copy(base)));
            *corr = corr_logpow(ctx, base, exp);
            *additive = true;
            return true;
        }
        if (is_head(arg, "Times")) {                   /* Log[z1..zn] */
            if (!var_gate(arg, ctx)) return false;
            Expr** zs = arg->data.function.args;
            size_t n  = arg->data.function.arg_count;
            *main = split_additive_product("Log", zs, n, NULL);
            *corr = corr_logprod(ctx, zs, n);
            *additive = true;
            return true;
        }
        return false;
    }

    /* ---- Arg ---- */
    if (strcmp(h, "Arg") == 0 && ac == 1) {
        Expr* arg = a[0];
        if (is_head(arg, "Times")) {
            if (!var_gate(arg, ctx)) return false;
            Expr** zs = arg->data.function.args;
            size_t n  = arg->data.function.arg_count;
            *main = split_additive_product("Arg", zs, n, NULL);
            *corr = corr_argprod(ctx, zs, n);
            *additive = true;
            return true;
        }
        return false;
    }

    /* ---- Inverse trig: f^-1[f[x]] ---- */
    {
        const char* inner = inverse_of(h);
        if (inner && ac == 1 && is_head(a[0], inner) &&
            a[0]->data.function.arg_count == 1) {
            Expr* x = a[0]->data.function.args[0];
            if (!var_gate(x, ctx)) return false;
            if (ctx->mode != PE_AUTOMATIC && strcmp(h, "ArcTan") == 0)
                *main = arctan_tan_universal(x);
            else
                *main = expr_copy(x);                  /* principal branch */
            return true;
        }
    }

    return false;
}

/* ------------------------------------------------------------------ */
/* Top-down recursive driver                                          */
/* ------------------------------------------------------------------ */

static Expr* pe_rec(Expr* e, const PECtx* ctx) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy(e);

    const char* h = head_name(e);

    /* Thread over lists / equations / inequalities / logic. */
    if (is_thread_head(h)) {
        bool is_ineq = (strcmp(h, "Inequality") == 0);
        size_t n = e->data.function.arg_count;
        Expr** args = malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) {
            if (is_ineq && (i & 1u) == 1)              /* operator slot */
                args[i] = expr_copy(e->data.function.args[i]);
            else
                args[i] = pe_rec(e->data.function.args[i], ctx);
        }
        Expr* res = fnN(h, args, n);
        free(args);
        return res;
    }

    /* Top-down: try a rule at this node. Reprocess the rewritten "main" part
     * (it may expand further) but leave the correction term inert. */
    Expr* main = NULL; Expr* corr = NULL; bool additive = false;
    if (pe_apply_rule(e, ctx, &main, &corr, &additive)) {
        Expr* main2 = pe_rec(main, ctx);
        expr_free(main);
        if (!corr) return main2;
        return additive ? fn2("Plus", main2, corr) : fn2("Times", main2, corr);
    }

    /* No rule here: recurse into head and arguments. */
    size_t n = e->data.function.arg_count;
    Expr** args = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++)
        args[i] = pe_rec(e->data.function.args[i], ctx);
    Expr* res = expr_new_function(expr_copy(e->data.function.head), args, n);
    free(args);
    return res;
}

/* ================================================================== */
/* Assumption-based refinement (PE_ASSUME)                            */
/* ================================================================== */

/* Read a literal numeric expression (Integer / Real / BigInt / Rational). */
static bool pe_to_double(const Expr* e, double* out) {
    switch (e->type) {
        case EXPR_INTEGER: *out = (double)e->data.integer; return true;
        case EXPR_REAL:    *out = e->data.real;            return true;
        case EXPR_BIGINT:  *out = mpz_get_d(e->data.bigint); return true;
#ifdef USE_MPFR
        case EXPR_MPFR:    *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true;
#endif
        case EXPR_FUNCTION:
            if (is_head(e, "Rational") && e->data.function.arg_count == 2) {
                double p, q;
                if (pe_to_double(e->data.function.args[0], &p) &&
                    pe_to_double(e->data.function.args[1], &q) && q != 0.0) {
                    *out = p / q; return true;
                }
            }
            return false;
        default: return false;
    }
}

/* Extract numeric bounds on variable v from the assumption facts. Recognises
 * binary comparisons and Inequality chains with a numeric side. */
static void pe_var_bounds(const Expr* v, const PECtx* ctx,
                          double* lo, bool* haslo, double* hi, bool* hashi) {
    *haslo = *hashi = false;
    if (!ctx->ctx) return;

    for (size_t i = 0; i < ctx->ctx->count; i++) {
        Expr* f = ctx->ctx->facts[i];
        if (f->type != EXPR_FUNCTION) continue;
        const char* fh = head_name(f);

        /* Binary comparison: op[lhs, rhs]. */
        if (f->data.function.arg_count == 2 &&
            (strcmp(fh,"Less")==0 || strcmp(fh,"LessEqual")==0 ||
             strcmp(fh,"Greater")==0 || strcmp(fh,"GreaterEqual")==0)) {
            Expr* L = f->data.function.args[0];
            Expr* R = f->data.function.args[1];
            double num;
            if (expr_eq(L,(Expr*)v) && pe_to_double(R,&num)) {     /* v op num */
                if (strcmp(fh,"Less")==0 || strcmp(fh,"LessEqual")==0) {
                    if (!*hashi || num < *hi) { *hi=num; *hashi=true; }
                } else {
                    if (!*haslo || num > *lo) { *lo=num; *haslo=true; }
                }
            } else if (expr_eq(R,(Expr*)v) && pe_to_double(L,&num)) { /* num op v */
                if (strcmp(fh,"Less")==0 || strcmp(fh,"LessEqual")==0) {
                    if (!*haslo || num > *lo) { *lo=num; *haslo=true; }
                } else {
                    if (!*hashi || num < *hi) { *hi=num; *hashi=true; }
                }
            }
            continue;
        }

        /* Inequality chain: Inequality[e0, op0, e1, op1, e2, ...]. */
        if (strcmp(fh,"Inequality")==0 && (f->data.function.arg_count & 1u)) {
            size_t na = f->data.function.arg_count;
            for (size_t k = 0; k + 2 < na; k += 2) {
                Expr* L  = f->data.function.args[k];
                Expr* op = f->data.function.args[k+1];
                Expr* R  = f->data.function.args[k+2];
                const char* o = (op->type==EXPR_SYMBOL)? op->data.symbol : "";
                double num;
                bool lessop = (strcmp(o,"Less")==0 || strcmp(o,"LessEqual")==0);
                if (expr_eq(L,(Expr*)v) && pe_to_double(R,&num)) {       /* v < num */
                    if (lessop) { if (!*hashi || num<*hi){*hi=num;*hashi=true;} }
                    else        { if (!*haslo || num>*lo){*lo=num;*haslo=true;} }
                } else if (expr_eq(R,(Expr*)v) && pe_to_double(L,&num)) { /* num < v */
                    if (lessop) { if (!*haslo || num>*lo){*lo=num;*haslo=true;} }
                    else        { if (!*hashi || num<*hi){*hi=num;*hashi=true;} }
                }
            }
        }
    }
}

/* Sign of a real expression under the assumptions: +1, -1, or 0 (unknown). */
static int pe_real_sign(const Expr* e, const PECtx* ctx) {
    double num;
    if (pe_to_double(e, &num)) return num > 0 ? 1 : (num < 0 ? -1 : 0);

    if (e->type == EXPR_SYMBOL) {
        if (e->data.symbol == SYM_Pi || e->data.symbol == SYM_E) return 1;
        double lo, hi; bool hl, hh;
        pe_var_bounds(e, ctx, &lo, &hl, &hi, &hh);
        if (hl && lo >= 0) return 1;
        if (hh && hi <= 0) return -1;
        return 0;
    }
    if (is_head(e, "Times")) {                 /* product of factor signs */
        int s = 1;
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            int fi = pe_real_sign(e->data.function.args[i], ctx);
            if (fi == 0) return 0;
            s *= fi;
        }
        return s;
    }
    return 0;
}

/* Is e known to be real under the assumptions (good enough to factor out)? */
static bool pe_is_real(const Expr* e, const PECtx* ctx) {
    double num;
    if (pe_to_double(e, &num)) return true;
    if (e->type == EXPR_SYMBOL) {
        if (e->data.symbol == SYM_Pi || e->data.symbol == SYM_E) return true;
        if (ctx->ctx && assume_known_real(ctx->ctx, (Expr*)e)) return true;
        double lo, hi; bool hl, hh;
        pe_var_bounds(e, ctx, &lo, &hl, &hi, &hh);
        return hl || hh;                       /* a numeric bound implies real */
    }
    if (is_head(e, "Arg")) return true;
    return false;
}

/* Arg[e] under assumptions: 0 (positive), Pi (negative), else Arg[e]. */
static Expr* pe_arg_under(const Expr* e, const PECtx* ctx) {
    int s = pe_real_sign(e, ctx);
    if (s > 0) return pe_int(0);
    if (s < 0) return pe_pi();
    return fn1("Arg", expr_copy((Expr*)e));
}

/* Imaginary part of e under assumptions. Returns a refined expr (caller owns).
 * Handles the forms produced by the correction terms: real constants, Pi/E,
 * Log[v] (-> Arg[v]), products with real factors, and sums. */
static Expr* pe_im_under(const Expr* e, const PECtx* ctx) {
    double num;
    if (pe_to_double(e, &num)) return pe_int(0);
    if (e->type == EXPR_SYMBOL) {
        if (e->data.symbol == SYM_Pi || e->data.symbol == SYM_E) return pe_int(0);
        if (pe_is_real(e, ctx)) return pe_int(0);
        return fn1("Im", expr_copy((Expr*)e));
    }
    if (is_head(e, "Log") && e->data.function.arg_count == 1)
        return pe_arg_under(e->data.function.args[0], ctx);   /* Im[Log[v]] = Arg[v] */

    if (is_head(e, "Plus")) {
        size_t n = e->data.function.arg_count;
        Expr** t = malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) t[i] = pe_im_under(e->data.function.args[i], ctx);
        return fnN("Plus", t, n);
    }
    if (is_head(e, "Times")) {
        /* Pull out the real factors: Im[(real) w] = (real) Im[w]. */
        size_t n = e->data.function.arg_count;
        Expr** realf = malloc(sizeof(Expr*) * (n + 1));
        Expr** restf = malloc(sizeof(Expr*) * n);
        size_t nr = 0, nrest = 0;
        for (size_t i = 0; i < n; i++) {
            Expr* fi = e->data.function.args[i];
            if (pe_is_real(fi, ctx)) realf[nr++] = expr_copy(fi);
            else                     restf[nrest++] = expr_copy(fi);
        }
        Expr* result;
        if (nrest == 0) {                       /* fully real -> Im = 0 */
            for (size_t i = 0; i < nr; i++) expr_free(realf[i]);
            result = pe_int(0);
        } else if (nr == 0) {                   /* nothing real to pull out */
            for (size_t i = 0; i < nrest; i++) expr_free(restf[i]);
            result = fn1("Im", expr_copy((Expr*)e));
        } else {
            Expr* rest = (nrest == 1) ? restf[0] : fnN("Times", restf, nrest);
            Expr* imrest = pe_im_under(rest, ctx);   /* (real) Im[rest] */
            expr_free(rest);                          /* frees the rest factors */
            realf[nr++] = imrest;
            result = fnN("Times", realf, nr);
        }
        free(realf);
        free(restf);
        return result;
    }
    return fn1("Im", expr_copy((Expr*)e));
}

/* Interval [lo,hi] of e over assumption-bounded atoms. Arg[..] residuals are
 * treated as bounded in [-Pi, Pi]. Returns false if e is not reducible. */
static bool pe_interval(const Expr* e, const PECtx* ctx, double* lo, double* hi) {
    double num;
    if (pe_to_double(e, &num)) { *lo = *hi = num; return true; }

    if (e->type == EXPR_SYMBOL) {
        if (e->data.symbol == SYM_Pi) { *lo = *hi = M_PI; return true; }
        double l, h; bool hl, hh;
        pe_var_bounds(e, ctx, &l, &hl, &h, &hh);
        if (hl && hh) { *lo = l; *hi = h; return true; }
        return false;
    }
    if (is_head(e, "Arg")) { *lo = -M_PI; *hi = M_PI; return true; }

    if (is_head(e, "Plus")) {
        double s_lo = 0, s_hi = 0;
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            double l, h;
            if (!pe_interval(e->data.function.args[i], ctx, &l, &h)) return false;
            s_lo += l; s_hi += h;
        }
        *lo = s_lo; *hi = s_hi; return true;
    }
    if (is_head(e, "Times")) {
        double p_lo = 1, p_hi = 1;
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            double l, h;
            if (!pe_interval(e->data.function.args[i], ctx, &l, &h)) return false;
            double c1 = p_lo*l, c2 = p_lo*h, c3 = p_hi*l, c4 = p_hi*h;
            double mn = c1, mx = c1;
            if (c2<mn) mn=c2; if (c3<mn) mn=c3; if (c4<mn) mn=c4;
            if (c2>mx) mx=c2; if (c3>mx) mx=c3; if (c4>mx) mx=c4;
            p_lo = mn; p_hi = mx;
        }
        *lo = p_lo; *hi = p_hi; return true;
    }
    if (is_head(e, "Power") && e->data.function.arg_count == 2) {
        /* Only the simple reciprocal of a constant (e.g. 1/(2 Pi)) is needed. */
        double bl, bh, el, eh;
        if (pe_interval(e->data.function.args[0], ctx, &bl, &bh) &&
            pe_interval(e->data.function.args[1], ctx, &el, &eh) &&
            el == eh && el == -1.0 && bl == bh && bl != 0.0) {
            *lo = *hi = 1.0 / bl; return true;
        }
        return false;
    }
    return false;
}

/* If Floor[arg] is constant over the assumption-bounded interval, return that
 * integer; otherwise NULL. */
static Expr* pe_floor_const(const Expr* arg, const PECtx* ctx) {
    double lo, hi;
    if (!pe_interval(arg, ctx, &lo, &hi)) return NULL;
    if (lo > hi) { double t = lo; lo = hi; hi = t; }
    double mid = 0.5 * (lo + hi);
    double k = floor(mid);
    const double eps = 1e-9;
    if (k <= lo + eps && hi <= k + 1.0 + eps)
        return pe_int((int64_t)k);
    return NULL;
}

/* Collapse E^(n I Pi) to (-1)^n for integer n. After Floor terms refine to
 * integers the multiplicative corrections become E^(integer I Pi); reducing
 * them yields results like Sqrt[z^2] -> -z under z < 0. */
static Expr* pe_reduce_ipi(Expr* e) {
    if (e->type != EXPR_FUNCTION) return expr_copy(e);

    size_t n = e->data.function.arg_count;
    Expr** args = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) args[i] = pe_reduce_ipi(e->data.function.args[i]);
    Expr* e2 = expr_new_function(expr_copy(e->data.function.head), args, n);
    free(args);

    if (is_head(e2, "Power") && e2->data.function.arg_count == 2 &&
        e2->data.function.args[0]->type == EXPR_SYMBOL &&
        strcmp(e2->data.function.args[0]->data.symbol, "E") == 0) {
        Expr* X = e2->data.function.args[1];
        Expr* ratio = eval_and_free(
            fn2("Times", expr_copy(X),
                fn2("Power", fn2("Times", pe_I(), pe_pi()), pe_int(-1))));
        if (ratio->type == EXPR_INTEGER) {
            Expr* r = eval_and_free(fn2("Power", pe_int(-1), expr_copy(ratio)));
            expr_free(ratio);
            expr_free(e2);
            return r;
        }
        expr_free(ratio);
    }
    return e2;
}

/* Walk the expression refining Arg / Im / Floor under the assumptions. */
static Expr* pe_refine(Expr* e, const PECtx* ctx) {
    if (e->type != EXPR_FUNCTION) return expr_copy(e);
    const char* h = head_name(e);

    if (strcmp(h, "Arg") == 0 && e->data.function.arg_count == 1) {
        Expr* inner = pe_refine(e->data.function.args[0], ctx);
        Expr* r = pe_arg_under(inner, ctx);
        expr_free(inner);
        return r;
    }
    if (strcmp(h, "Im") == 0 && e->data.function.arg_count == 1) {
        Expr* inner = pe_refine(e->data.function.args[0], ctx);
        Expr* r = pe_im_under(inner, ctx);
        expr_free(inner);
        return r;
    }
    if (strcmp(h, "Floor") == 0 && e->data.function.arg_count == 1) {
        Expr* inner = pe_refine(e->data.function.args[0], ctx);
        Expr* k = pe_floor_const(inner, ctx);
        if (k) { expr_free(inner); return k; }
        return fn1("Floor", inner);
    }

    size_t n = e->data.function.arg_count;
    Expr** args = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) args[i] = pe_refine(e->data.function.args[i], ctx);
    Expr* res = expr_new_function(expr_copy(e->data.function.head), args, n);
    free(args);
    return res;
}

/* ------------------------------------------------------------------ */
/* Argument parsing + builtin                                          */
/* ------------------------------------------------------------------ */

extern bool is_rule_with_lhs(const Expr* e, const char* lhs_symbol);

Expr* builtin_powerexpand(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    if (argc < 1 || argc > 3) {
        fprintf(stderr,
                "PowerExpand::argt: PowerExpand called with %zu arguments; "
                "1 or 2 arguments are expected.\n", argc);
        return NULL;
    }

    Expr** av = res->data.function.args;
    Expr* expr = av[0];
    Expr* varlist = NULL;          /* {x1, ...} restriction */
    Expr* assum = NULL;            /* Assumptions -> value */

    for (size_t i = 1; i < argc; i++) {
        Expr* a = av[i];
        if (is_rule_with_lhs(a, "Assumptions")) {
            assum = a->data.function.args[1];
        } else if (is_head(a, "List")) {
            varlist = a;
        } else if (varlist == NULL) {
            varlist = a;           /* a bare symbol restricts to that variable */
        }
    }

    PECtx ctx;
    ctx.mode  = PE_AUTOMATIC;
    ctx.vars  = NULL;
    ctx.nvars = 0;
    ctx.ctx   = NULL;

    if (varlist) {
        if (is_head(varlist, "List")) {
            ctx.vars  = varlist->data.function.args;
            ctx.nvars = varlist->data.function.arg_count;
        } else {
            ctx.vars  = &varlist;  /* single variable */
            ctx.nvars = 1;
        }
    }

    /* Assumptions option. Automatic keeps the default behaviour; True emits
     * the universal formulas; anything else triggers refinement. */
    bool is_automatic = (!assum) ||
        (assum->type == EXPR_SYMBOL && strcmp(assum->data.symbol, "Automatic") == 0);
    bool is_true = (assum && assum->type == EXPR_SYMBOL &&
                    strcmp(assum->data.symbol, "True") == 0);
    if (!is_automatic) ctx.mode = is_true ? PE_TRUE : PE_ASSUME;

    if (ctx.mode == PE_ASSUME) ctx.ctx = assume_ctx_from_expr(assum);

    Expr* out = pe_rec(expr, &ctx);
    if (ctx.mode == PE_ASSUME) {
        Expr* refined = pe_refine(out, &ctx);
        expr_free(out);
        Expr* reduced = pe_reduce_ipi(refined);   /* E^(n I Pi) -> (-1)^n */
        expr_free(refined);
        out = reduced;
    }
    if (ctx.ctx) assume_ctx_free(ctx.ctx);
    return out;
}

void expand_power_init(void) {
    symtab_add_builtin("PowerExpand", builtin_powerexpand);
    symtab_get_def("PowerExpand")->attributes |= ATTR_PROTECTED;
}
