/*
 * reduce_realfn.c
 *
 * Elementary-real-function support for `Reduce` over the Reals (Phase 9).  See
 * reduce_realfn.h.  Three services:
 *
 *   - reduce_stmt_has_realfn:      does the statement need this path at all?
 *   - reduce_real_domain_collect:  the head->domain table (constraints/boundaries)
 *   - reduce_realfn_preprocess:    Mod->Floor, integer-part isolation, Abs split
 *
 * Everything here is a pure Expr->Expr / Expr->data transformation; the actual
 * sign-diagram solving is in reduce_realdiag.c.
 */
#include "reduce_realfn.h"

#include "eval.h"
#include "sym_names.h"

#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ *
 *  Small helpers                                                      *
 * ------------------------------------------------------------------ */

static bool is_head(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == name;
}

static bool contains_symbol(const Expr* e, const char* name) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == name;
    if (e->type == EXPR_FUNCTION) {
        if (contains_symbol(e->data.function.head, name)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (contains_symbol(e->data.function.args[i], name)) return true;
    }
    return false;
}

static bool contains_x(const Expr* e, const Expr* x) {
    return x && x->type == EXPR_SYMBOL && contains_symbol(e, x->data.symbol.name);
}

/* True iff `e` contains no Abs node anywhere. */
static bool is_absfree(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return true;
    if (is_head(e, SYM_Abs)) return false;
    if (!is_absfree(e->data.function.head)) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!is_absfree(e->data.function.args[i])) return false;
    return true;
}

/* Structural replace: a fresh copy of `e` with every subtree equal to `target`
 * replaced by a copy of `repl`. */
static Expr* subst_expr(const Expr* e, const Expr* target, const Expr* repl) {
    if (expr_eq((Expr*)e, (Expr*)target)) return expr_copy((Expr*)repl);
    if (!e || e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    Expr* head = subst_expr(e->data.function.head, target, repl);
    size_t n = e->data.function.arg_count;
    Expr** args = (n > 0) ? malloc(n * sizeof(Expr*)) : NULL;
    for (size_t i = 0; i < n; i++)
        args[i] = subst_expr(e->data.function.args[i], target, repl);
    Expr* out = expr_new_function(head, args, n);
    free(args);
    return out;
}

/* A two-argument relation head[a,b], UNEVALUATED (adopts a, b). */
static Expr* rel2(const char* head, Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(head), (Expr*[]){ a, b }, 2);
}

/* Evaluate a scalar helper expression to normal form (adopts e). */
static Expr* enorm(Expr* e) { return eval_and_free(e); }

/* u + k  (k an integer), evaluated. */
static Expr* plus_int(const Expr* u, long k) {
    return enorm(expr_new_function(expr_new_symbol(SYM_Plus),
        (Expr*[]){ expr_copy((Expr*)u), expr_new_integer(k) }, 2));
}
/* k - u  (k an integer), evaluated. */
static Expr* int_minus(long k, const Expr* u) {
    return enorm(expr_new_function(expr_new_symbol(SYM_Plus),
        (Expr*[]){ expr_new_integer(k),
                   expr_new_function(expr_new_symbol(SYM_Times),
                       (Expr*[]){ expr_new_integer(-1), expr_copy((Expr*)u) }, 2) }, 2));
}

/* Is `e` a positive numeric constant (via N)?  Handles Pi, 2 Pi, Sqrt[2], ... */
static bool positive_const(const Expr* e) {
    Expr* v = eval_and_free(expr_new_function(expr_new_symbol(SYM_N),
        (Expr*[]){ expr_copy((Expr*)e) }, 1));
    double d = 0.0; bool ok = true;
    if (v->type == EXPR_REAL)         d = v->data.real;
    else if (v->type == EXPR_INTEGER) d = (double)v->data.integer;
    else                              ok = false;
    expr_free(v);
    return ok && d > 0.0;
}

/* ------------------------------------------------------------------ *
 *  Head -> real-domain table                                          *
 * ------------------------------------------------------------------ */

static void rdcon_push(RDomCon** cons, int* n, int* cap, Expr* poly, bool strict) {
    if (*n == *cap) { *cap = *cap ? *cap * 2 : 8; *cons = realloc(*cons, (size_t)*cap * sizeof(RDomCon)); }
    (*cons)[*n].poly = poly;
    (*cons)[*n].strict = strict;
    (*n)++;
}

/* Even-order rational-power radical Power[u, p/q] (q even) needs u>=0 (u>0 when
 * the exponent is negative).  Returns 1 if it added a constraint, 0 otherwise. */
static int power_domain(const Expr* e, const Expr* x, RDomCon** cons, int* n, int* cap) {
    if (!is_head(e, SYM_Power) || e->data.function.arg_count != 2) return 0;
    const Expr* u   = e->data.function.args[0];
    const Expr* exp = e->data.function.args[1];
    if (!is_head(exp, SYM_Rational) || exp->data.function.arg_count != 2) return 0;
    const Expr* pnum = exp->data.function.args[0];
    const Expr* pden = exp->data.function.args[1];
    if (pden->type != EXPR_INTEGER || (pden->data.integer % 2) != 0) return 0;   /* q even */
    if (!contains_x(u, x)) return 0;
    bool neg_exp = (pnum->type == EXPR_INTEGER && pnum->data.integer < 0);
    rdcon_push(cons, n, cap, expr_copy((Expr*)u), neg_exp);   /* u>=0, or u>0 */
    return 1;
}

void reduce_real_domain_collect(const Expr* e, const Expr* x,
                                RDomCon** cons, int* n, int* cap) {
    if (!e || e->type != EXPR_FUNCTION) return;

    /* Classify this node. */
    const char* h = e->data.function.head->type == EXPR_SYMBOL
                  ? e->data.function.head->data.symbol.name : NULL;
    size_t argc = e->data.function.arg_count;
    const Expr* a0 = argc >= 1 ? e->data.function.args[0] : NULL;

    if (h == SYM_Power) {
        power_domain(e, x, cons, n, cap);
    } else if (h == SYM_Sqrt && argc == 1 && contains_x(a0, x)) {
        rdcon_push(cons, n, cap, expr_copy((Expr*)a0), false);        /* u>=0 */
    } else if (h == SYM_Log && contains_x(argc == 2 ? e->data.function.args[1] : a0, x)) {
        const Expr* u = (argc == 2) ? e->data.function.args[1] : a0;  /* Log[u] or Log[b,u] */
        rdcon_push(cons, n, cap, expr_copy((Expr*)u), true);          /* u>0 */
    } else if ((h == SYM_ArcSin || h == SYM_ArcCos) && argc == 1 && contains_x(a0, x)) {
        rdcon_push(cons, n, cap, plus_int(a0, 1), false);             /* u+1>=0 */
        rdcon_push(cons, n, cap, int_minus(1, a0), false);            /* 1-u>=0 */
    } else if (h == SYM_ArcTanh && argc == 1 && contains_x(a0, x)) {
        rdcon_push(cons, n, cap, plus_int(a0, 1), true);              /* u+1>0 */
        rdcon_push(cons, n, cap, int_minus(1, a0), true);             /* 1-u>0 */
    } else if (h == SYM_ArcCosh && argc == 1 && contains_x(a0, x)) {
        rdcon_push(cons, n, cap, plus_int(a0, -1), false);           /* u-1>=0 */
    } else if (h == SYM_ArcSech && argc == 1 && contains_x(a0, x)) {
        rdcon_push(cons, n, cap, expr_copy((Expr*)a0), true);         /* u>0 */
        rdcon_push(cons, n, cap, int_minus(1, a0), false);           /* 1-u>=0 */
    }

    /* Recurse into every argument (nested radicals / logs). */
    for (size_t i = 0; i < argc; i++)
        reduce_real_domain_collect(e->data.function.args[i], x, cons, n, cap);
}

/* ------------------------------------------------------------------ *
 *  Detector                                                           *
 * ------------------------------------------------------------------ */

static bool node_is_realfn(const Expr* e, const Expr* x) {
    if (!e || e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    size_t argc = e->data.function.arg_count;
    const Expr* a0 = argc >= 1 ? e->data.function.args[0] : NULL;

    if (h == SYM_Abs || h == SYM_Floor || h == SYM_Ceiling || h == SYM_Round
        || h == SYM_IntegerPart || h == SYM_Mod || h == SYM_Log
        || h == SYM_ArcSin || h == SYM_ArcCos || h == SYM_ArcTanh
        || h == SYM_ArcCosh || h == SYM_ArcSech || h == SYM_Sqrt
        || h == SYM_Max || h == SYM_Min
        /* piecewise selectors (case-split away in preprocessing) */
        || h == SYM_Piecewise || h == SYM_Sign || h == SYM_UnitStep
        || h == SYM_Ramp || h == SYM_Clip || h == SYM_HeavisideTheta
        || h == SYM_Boole || h == SYM_UnitBox || h == SYM_FractionalPart) {
        /* any of these with x under it is enough */
        for (size_t i = 0; i < argc; i++)
            if (contains_x(e->data.function.args[i], x)) return true;
    }
    /* Non-integer rational power of an x-expression (radical). */
    if (h == SYM_Power && argc == 2) {
        const Expr* exp = e->data.function.args[1];
        if (is_head(exp, SYM_Rational) && contains_x(a0, x)) return true;
    }
    return false;
}

bool reduce_stmt_has_realfn(const Expr* e, const Expr* x) {
    if (!e) return false;
    if (node_is_realfn(e, x)) return true;
    if (e->type == EXPR_FUNCTION) {
        if (reduce_stmt_has_realfn(e->data.function.head, x)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (reduce_stmt_has_realfn(e->data.function.args[i], x)) return true;
    }
    return false;
}

/* True iff `e` mentions any of the `nv` reduce variables. */
static bool contains_any_var(const Expr* e, Expr** vars, int nv) {
    for (int i = 0; i < nv; i++) if (contains_x(e, vars[i])) return true;
    return false;
}

/* A selector head (Abs/Max/Min + the fully-polynomial-reducible piecewise
 * builtins) with a reduce variable under it -- the multivariate dispatch's cue
 * to run reduce_piecewise_preprocess before the FM/CAD engines.  (IntegerPart /
 * FractionalPart / Mod are excluded: they reduce to Floor/Ceiling, which are
 * univariate-only, so a multivariate occurrence is left to decline.) */
bool reduce_stmt_has_piecewise(const Expr* e, Expr** vars, int nv) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        if ((h == SYM_Abs || h == SYM_Max || h == SYM_Min || h == SYM_Piecewise
             || h == SYM_Sign || h == SYM_UnitStep || h == SYM_Ramp || h == SYM_Clip
             || h == SYM_HeavisideTheta || h == SYM_Boole || h == SYM_UnitBox)
            && contains_any_var(e, vars, nv))
            return true;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (reduce_stmt_has_piecewise(e->data.function.args[i], vars, nv)) return true;
    return false;
}

static bool is_sqrt_radical(const Expr* e);   /* defined with the rationalizer below */

/* A square-root radical Power[u, 1/2] with a reduce variable in its radicand --
 * the multivariate dispatch's cue to rationalize radicals into polynomial
 * constraints before FM/CAD. */
bool reduce_stmt_has_radical(const Expr* e, Expr** vars, int nv) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (is_sqrt_radical(e) && contains_any_var(e->data.function.args[0], vars, nv)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (reduce_stmt_has_radical(e->data.function.args[i], vars, nv)) return true;
    return false;
}

/* ------------------------------------------------------------------ *
 *  Preprocessing 1: Mod -> Floor                                      *
 * ------------------------------------------------------------------ */

/* Replace every Mod[u,m] (m a positive constant) by u - m*Floor[u/m]. */
static Expr* subst_mod(const Expr* e, bool* changed) {
    if (!e || e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t n = e->data.function.arg_count;
    Expr* head = subst_mod(e->data.function.head, changed);
    Expr** args = (n > 0) ? malloc(n * sizeof(Expr*)) : NULL;
    for (size_t i = 0; i < n; i++) args[i] = subst_mod(e->data.function.args[i], changed);
    Expr* rebuilt = expr_new_function(head, args, n);
    free(args);

    if (is_head(rebuilt, SYM_Mod) && rebuilt->data.function.arg_count == 2) {
        Expr* u = rebuilt->data.function.args[0];
        Expr* m = rebuilt->data.function.args[1];
        if (positive_const(m)) {
            *changed = true;
            /* u/m */
            Expr* um = expr_new_function(expr_new_symbol(SYM_Times),
                (Expr*[]){ expr_copy(u),
                           expr_new_function(expr_new_symbol(SYM_Power),
                               (Expr*[]){ expr_copy(m), expr_new_integer(-1) }, 2) }, 2);
            Expr* fl = expr_new_function(expr_new_symbol(SYM_Floor), (Expr*[]){ um }, 1);
            /* u - m*Floor[u/m] */
            Expr* out = expr_new_function(expr_new_symbol(SYM_Plus),
                (Expr*[]){ expr_copy(u),
                           expr_new_function(expr_new_symbol(SYM_Times),
                               (Expr*[]){ expr_new_integer(-1), expr_copy(m), fl }, 3) }, 2);
            expr_free(rebuilt);
            return out;
        }
    }
    return rebuilt;
}

/* ------------------------------------------------------------------ *
 *  Preprocessing 2: integer-part relation isolation                   *
 * ------------------------------------------------------------------ */

/* Find the single distinct Floor/Ceiling/Round/IntegerPart node in `e`; returns
 * it (borrowed) and sets *iphead to its head name, or NULL if there is not
 * exactly one distinct such node. */
static const Expr* find_single_ip(const Expr* e, const char** iphead, const Expr** found) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;
    const char* h = e->data.function.head->type == EXPR_SYMBOL
                  ? e->data.function.head->data.symbol.name : NULL;
    if ((h == SYM_Floor || h == SYM_Ceiling || h == SYM_Round || h == SYM_IntegerPart)
        && e->data.function.arg_count == 1) {
        if (*found && !expr_eq((Expr*)*found, (Expr*)e)) { *iphead = NULL; return NULL; }
        *found = e; *iphead = h;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        find_single_ip(e->data.function.args[i], iphead, found);
    return *found;
}

/* IntegerQ[c] via evaluate. */
static bool is_integer_expr(const Expr* c) {
    return c->type == EXPR_INTEGER || c->type == EXPR_BIGINT;
}

/* Emit the defining relation for  IP[v]  REL  c  (c a numeric constant), where
 * IP is Floor/Ceiling/Round.  Returns a freshly-owned statement Expr, or NULL to
 * decline. */
static Expr* ip_defining(const char* iphead, const Expr* v, const char* rel, const Expr* c) {
    /* floor(c), ceil(c) as integer Exprs */
    Expr* fc = enorm(expr_new_function(expr_new_symbol(SYM_Floor), (Expr*[]){ expr_copy((Expr*)c) }, 1));
    Expr* cc = enorm(expr_new_function(expr_new_symbol(SYM_Ceiling), (Expr*[]){ expr_copy((Expr*)c) }, 1));
    bool c_int = is_integer_expr(c);
    Expr* out = NULL;

    if (iphead == SYM_Floor) {
        if (rel == SYM_Equal) {
            if (!c_int) { out = expr_new_symbol(SYM_False); }
            else out = expr_new_function(expr_new_symbol(SYM_Inequality),
                (Expr*[]){ expr_copy((Expr*)c), expr_new_symbol(SYM_LessEqual),
                           expr_copy((Expr*)v), expr_new_symbol(SYM_Less),
                           plus_int(c, 1) }, 5);
        } else if (rel == SYM_Unequal) {
            if (!c_int) { out = expr_new_symbol(SYM_True); }
            else out = expr_new_function(expr_new_symbol(SYM_Or), (Expr*[]){
                rel2(SYM_Less, expr_copy((Expr*)v), expr_copy((Expr*)c)),
                rel2(SYM_GreaterEqual, expr_copy((Expr*)v), plus_int(c, 1)) }, 2);
        } else if (rel == SYM_LessEqual) {      /* Floor[v]<=c ⟺ v < floor(c)+1 */
            out = rel2(SYM_Less, expr_copy((Expr*)v), plus_int(fc, 1));
        } else if (rel == SYM_Less) {           /* Floor[v]<c ⟺ v < ceil(c) */
            out = rel2(SYM_Less, expr_copy((Expr*)v), expr_copy(cc));
        } else if (rel == SYM_GreaterEqual) {   /* Floor[v]>=c ⟺ v >= ceil(c) */
            out = rel2(SYM_GreaterEqual, expr_copy((Expr*)v), expr_copy(cc));
        } else if (rel == SYM_Greater) {        /* Floor[v]>c ⟺ v >= floor(c)+1 */
            out = rel2(SYM_GreaterEqual, expr_copy((Expr*)v), plus_int(fc, 1));
        }
    } else if (iphead == SYM_Ceiling) {
        /* Ceiling[v]=m is the least integer >= v  (m-1 < v <= m). */
        if (rel == SYM_Equal) {
            if (!c_int) { out = expr_new_symbol(SYM_False); }
            else out = expr_new_function(expr_new_symbol(SYM_Inequality),
                (Expr*[]){ plus_int(c, -1), expr_new_symbol(SYM_Less),
                           expr_copy((Expr*)v), expr_new_symbol(SYM_LessEqual),
                           expr_copy((Expr*)c) }, 5);
        } else if (rel == SYM_Unequal) {
            if (!c_int) { out = expr_new_symbol(SYM_True); }
            else out = expr_new_function(expr_new_symbol(SYM_Or), (Expr*[]){
                rel2(SYM_LessEqual, expr_copy((Expr*)v), plus_int(c, -1)),
                rel2(SYM_Greater,   expr_copy((Expr*)v), expr_copy((Expr*)c)) }, 2);
        } else if (rel == SYM_LessEqual) {      /* Ceiling[v]<=c ⟺ v <= floor(c)   */
            out = rel2(SYM_LessEqual, expr_copy((Expr*)v), expr_copy(fc));
        } else if (rel == SYM_Less) {           /* Ceiling[v]<c  ⟺ v <= ceil(c)-1  */
            out = rel2(SYM_LessEqual, expr_copy((Expr*)v), plus_int(cc, -1));
        } else if (rel == SYM_GreaterEqual) {   /* Ceiling[v]>=c ⟺ v > ceil(c)-1   */
            out = rel2(SYM_Greater, expr_copy((Expr*)v), plus_int(cc, -1));
        } else if (rel == SYM_Greater) {        /* Ceiling[v]>c  ⟺ v > floor(c)    */
            out = rel2(SYM_Greater, expr_copy((Expr*)v), expr_copy(fc));
        }
    } else if (iphead == SYM_Round) {
        if (rel == SYM_Equal && c_int) {        /* n-1/2 <= v < n+1/2 (half-up) */
            Expr* lo = enorm(expr_new_function(expr_new_symbol(SYM_Plus),
                (Expr*[]){ expr_copy((Expr*)c),
                           expr_new_function(expr_new_symbol(SYM_Rational),
                               (Expr*[]){ expr_new_integer(-1), expr_new_integer(2) }, 2) }, 2));
            Expr* hi = enorm(expr_new_function(expr_new_symbol(SYM_Plus),
                (Expr*[]){ expr_copy((Expr*)c),
                           expr_new_function(expr_new_symbol(SYM_Rational),
                               (Expr*[]){ expr_new_integer(1), expr_new_integer(2) }, 2) }, 2));
            out = expr_new_function(expr_new_symbol(SYM_Inequality),
                (Expr*[]){ lo, expr_new_symbol(SYM_LessEqual),
                           expr_copy((Expr*)v), expr_new_symbol(SYM_Less), hi }, 5);
        }
    }
    expr_free(fc); expr_free(cc);
    return out;
}

/* Try to rewrite a relational leaf REL[L,R] that is linear in a single
 * integer-part node.  Returns the rewritten statement, or NULL to leave as-is. */
static Expr* try_ip_isolate(const Expr* rel, bool* changed) {
    if (!rel || rel->type != EXPR_FUNCTION || rel->data.function.arg_count != 2) return NULL;
    const char* h = rel->data.function.head->type == EXPR_SYMBOL
                  ? rel->data.function.head->data.symbol.name : NULL;
    if (h != SYM_Equal && h != SYM_Unequal && h != SYM_Less && h != SYM_LessEqual
        && h != SYM_Greater && h != SYM_GreaterEqual) return NULL;

    Expr* L = rel->data.function.args[0];
    Expr* R = rel->data.function.args[1];
    Expr* diff = enorm(expr_new_function(expr_new_symbol(SYM_Plus),
        (Expr*[]){ expr_copy(L),
                   expr_new_function(expr_new_symbol(SYM_Times),
                       (Expr*[]){ expr_new_integer(-1), expr_copy(R) }, 2) }, 2));

    const char* iphead = NULL; const Expr* ip = NULL;
    find_single_ip(diff, &iphead, &ip);
    if (!ip || !iphead) { expr_free(diff); return NULL; }
    const Expr* v = ip->data.function.args[0];

    /* Extract B = diff|_{ip=0}, A = diff|_{ip=1} - B; verify linearity at ip=2. */
    Expr* zero = expr_new_integer(0), *one = expr_new_integer(1), *two = expr_new_integer(2);
    Expr* B  = enorm(subst_expr(diff, ip, zero));
    Expr* d1 = enorm(subst_expr(diff, ip, one));
    Expr* d2 = enorm(subst_expr(diff, ip, two));
    expr_free(zero); expr_free(one); expr_free(two);
    Expr* A  = enorm(expr_new_function(expr_new_symbol(SYM_Plus),
        (Expr*[]){ expr_copy(d1),
                   expr_new_function(expr_new_symbol(SYM_Times),
                       (Expr*[]){ expr_new_integer(-1), expr_copy(B) }, 2) }, 2));

    /* Need A, B free of ip and of the variable inside v (constants), A numeric != 0,
     * and true linearity (d2 == 2A + B). */
    bool ok = true;
    Expr* twoAB = enorm(expr_new_function(expr_new_symbol(SYM_Plus),
        (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Times),
                       (Expr*[]){ expr_new_integer(2), expr_copy(A) }, 2),
                   expr_copy(B) }, 2));
    if (!expr_eq(d2, twoAB)) ok = false;
    expr_free(twoAB); expr_free(d1); expr_free(d2);
    /* A,B must be numeric constants: N gives a plain number and they must not
     * mention the ip's inner variable(s). */
    if (ok) {
        double dummy; (void)dummy;
        Expr* nA = eval_and_free(expr_new_function(expr_new_symbol(SYM_N), (Expr*[]){ expr_copy(A) }, 1));
        Expr* nB = eval_and_free(expr_new_function(expr_new_symbol(SYM_N), (Expr*[]){ expr_copy(B) }, 1));
        bool aok = (nA->type == EXPR_REAL || nA->type == EXPR_INTEGER);
        bool bok = (nB->type == EXPR_REAL || nB->type == EXPR_INTEGER);
        double av = nA->type == EXPR_REAL ? nA->data.real : (nA->type == EXPR_INTEGER ? (double)nA->data.integer : 0.0);
        if (!aok || !bok || av == 0.0) ok = false;
        expr_free(nA); expr_free(nB);
    }
    if (!ok) { expr_free(diff); expr_free(A); expr_free(B); return NULL; }

    /* c = -B/A ; flip the relation when A < 0. */
    Expr* c = enorm(expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ expr_new_integer(-1), expr_copy(B),
                   expr_new_function(expr_new_symbol(SYM_Power),
                       (Expr*[]){ expr_copy(A), expr_new_integer(-1) }, 2) }, 3));
    Expr* nA2 = eval_and_free(expr_new_function(expr_new_symbol(SYM_N), (Expr*[]){ expr_copy(A) }, 1));
    bool neg = (nA2->type == EXPR_REAL && nA2->data.real < 0.0)
            || (nA2->type == EXPR_INTEGER && nA2->data.integer < 0);
    expr_free(nA2);

    const char* r = h;
    if (neg) {
        if (r == SYM_Less) r = SYM_Greater; else if (r == SYM_Greater) r = SYM_Less;
        else if (r == SYM_LessEqual) r = SYM_GreaterEqual; else if (r == SYM_GreaterEqual) r = SYM_LessEqual;
    }

    Expr* out = ip_defining(iphead, v, r, c);
    if (out) *changed = true;
    expr_free(diff); expr_free(A); expr_free(B); expr_free(c);
    return out;
}

/* ------------------------------------------------------------------ *
 *  Preprocessing 2b: general single-integer-part polynomial relation  *
 * ------------------------------------------------------------------ */

/* Translate an integer-domain solution `s` for the fresh variable `kname` (a
 * boolean combination of `kname == <int>` atoms, True, and False) into the
 * corresponding statement in the original variable, replacing every atom by the
 * defining relation for ip[v] == <int>.  Sets *ok = false (and returns NULL) on
 * any node that is not of this restricted shape, so a solution the translator
 * cannot justify makes the whole rewrite decline. */
/* The relation as seen from the other operand: a REL b  <=>  b flip(REL) a. */
static const char* flip_rel(const char* rel) {
    if (rel == SYM_Less)         return SYM_Greater;
    if (rel == SYM_Greater)      return SYM_Less;
    if (rel == SYM_LessEqual)    return SYM_GreaterEqual;
    if (rel == SYM_GreaterEqual) return SYM_LessEqual;
    return rel;   /* Equal / Unequal are symmetric */
}

static Expr* translate_ksol(const Expr* s, const char* kname,
                            const char* iphead, const Expr* v, bool* ok) {
    if (!*ok) return NULL;
    if (!s) { *ok = false; return NULL; }
    if (s->type == EXPR_SYMBOL) {
        if (s->data.symbol.name == SYM_True || s->data.symbol.name == SYM_False)
            return expr_copy((Expr*)s);
        *ok = false; return NULL;                 /* a bare `k` etc. -> decline */
    }
    if (s->type != EXPR_FUNCTION) { *ok = false; return NULL; }
    const char* h = s->data.function.head->type == EXPR_SYMBOL
                  ? s->data.function.head->data.symbol.name : NULL;

    if (h == SYM_And || h == SYM_Or || h == SYM_Not) {
        size_t n = s->data.function.arg_count;
        Expr** args = (n > 0) ? malloc(n * sizeof(Expr*)) : NULL;
        for (size_t i = 0; i < n; i++)
            args[i] = translate_ksol(s->data.function.args[i], kname, iphead, v, ok);
        if (!*ok) {
            for (size_t i = 0; i < n; i++) if (args[i]) expr_free(args[i]);
            free(args);
            return NULL;
        }
        Expr* out = expr_new_function(expr_copy(s->data.function.head), args, n);
        free(args);
        return out;
    }

    /* An atomic  k REL <int>  (either operand order): translate REL, as seen
     * from k, into the ip's defining relation on its inner value v. */
    bool is_rel = (h == SYM_Equal || h == SYM_Unequal || h == SYM_Less
                   || h == SYM_LessEqual || h == SYM_Greater || h == SYM_GreaterEqual);
    if (is_rel && s->data.function.arg_count == 2) {
        const Expr* A = s->data.function.args[0];
        const Expr* B = s->data.function.args[1];
        const char* rel = NULL; const Expr* cst = NULL;
        if (A->type == EXPR_SYMBOL && A->data.symbol.name == kname)      { rel = h;            cst = B; }
        else if (B->type == EXPR_SYMBOL && B->data.symbol.name == kname) { rel = flip_rel(h);  cst = A; }
        if (rel && cst && is_integer_expr(cst)) {
            Expr* d = ip_defining(iphead, v, rel, cst);   /* e.g. Floor[v]<=n <=> v<n+1 */
            if (d) return d;
        }
    }
    *ok = false; return NULL;
}

/* Handle a relational leaf that is a polynomial (any degree) in a SINGLE
 * Floor/Ceiling/Round node, which the linear try_ip_isolate could not isolate.
 * Since the node is integer-valued, substitute it by a fresh integer variable k
 * and solve the relation over the Integers.  When that yields a bounded set --
 * a clean boolean combination of `k == <int>` -- expand each value back to its
 * defining x-interval; otherwise decline (unbounded integer sets come back
 * unevaluated from integer Reduce, and an unjustifiable rewrite is never
 * emitted).  Returns the rewritten statement, or NULL to leave as-is. */
static Expr* try_ip_general(const Expr* rel, const Expr* x, bool* changed) {
    if (!rel || rel->type != EXPR_FUNCTION || rel->data.function.arg_count != 2) return NULL;
    const char* h = rel->data.function.head->type == EXPR_SYMBOL
                  ? rel->data.function.head->data.symbol.name : NULL;
    if (h != SYM_Equal && h != SYM_Unequal && h != SYM_Less && h != SYM_LessEqual
        && h != SYM_Greater && h != SYM_GreaterEqual) return NULL;

    const char* iphead = NULL; const Expr* ip = NULL;
    find_single_ip(rel, &iphead, &ip);
    if (!ip || !iphead) return NULL;
    const Expr* v = ip->data.function.args[0];

    /* Fresh integer variable; bail if (absurdly) it already occurs. */
    Expr* ksym = expr_new_symbol("$IPk");
    const char* kname = ksym->data.symbol.name;
    if (contains_symbol(rel, kname)) { expr_free(ksym); return NULL; }

    Expr* rel_k = subst_expr(rel, ip, ksym);           /* copies ksym */
    /* The relation must be pure in k: if x survives outside the ip node, this is
     * a mixed Floor[x]/x statement the substitution cannot linearise -> decline. */
    if (contains_x(rel_k, x)) { expr_free(rel_k); expr_free(ksym); return NULL; }

    Expr* sols = eval_and_free(expr_new_function(expr_new_symbol(SYM_Reduce),
        (Expr*[]){ rel_k, ksym, expr_new_symbol(SYM_Integers) }, 3));  /* adopts rel_k, ksym */

    bool ok = true;
    Expr* out = translate_ksol(sols, kname, iphead, v, &ok);
    expr_free(sols);
    if (!ok || !out) { if (out) expr_free(out); return NULL; }
    *changed = true;
    return out;
}

/* Walk the logical tree; rewrite integer-part relational leaves. */
static Expr* rewrite_ip_leaves(const Expr* e, const Expr* x, bool* changed) {
    if (!e || e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    const char* h = e->data.function.head->type == EXPR_SYMBOL
                  ? e->data.function.head->data.symbol.name : NULL;
    bool logical = (h == SYM_And || h == SYM_Or || h == SYM_Not
                    || h == SYM_Implies || h == SYM_Xor || h == SYM_Inequality);
    bool relational = (h == SYM_Equal || h == SYM_Unequal || h == SYM_Less
                       || h == SYM_LessEqual || h == SYM_Greater || h == SYM_GreaterEqual);

    if (relational) {
        Expr* iso = try_ip_isolate(e, changed);      /* linear fast path */
        if (iso) return iso;
        iso = try_ip_general(e, x, changed);         /* polynomial-in-one-ip path */
        if (iso) return iso;
        return expr_copy((Expr*)e);
    }
    if (logical) {
        size_t n = e->data.function.arg_count;
        Expr** args = (n > 0) ? malloc(n * sizeof(Expr*)) : NULL;
        for (size_t i = 0; i < n; i++) args[i] = rewrite_ip_leaves(e->data.function.args[i], x, changed);
        Expr* out = expr_new_function(expr_copy(e->data.function.head), args, n);
        free(args);
        return out;
    }
    return expr_copy((Expr*)e);
}

/* ------------------------------------------------------------------ *
 *  Preprocessing 3: Abs sign-splitting                                *
 * ------------------------------------------------------------------ */

/* Locate an innermost Abs[u] (u itself Abs-free); returns the Abs node
 * (borrowed) or NULL. */
static const Expr* find_inner_abs(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        const Expr* r = find_inner_abs(e->data.function.args[i]);
        if (r) return r;
    }
    if (is_head(e, SYM_Abs) && e->data.function.arg_count == 1
        && is_absfree(e->data.function.args[0]))
        return e;
    return NULL;
}

static Expr* eliminate_abs(const Expr* e, bool* changed) {
    const Expr* absnode = find_inner_abs(e);
    if (!absnode) return expr_copy((Expr*)e);
    *changed = true;
    const Expr* u = absnode->data.function.args[0];

    /* pos:  u >= 0  &&  e[Abs[u] -> u] */
    Expr* e_pos = subst_expr(e, absnode, u);
    Expr* negu = enorm(expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ expr_new_integer(-1), expr_copy((Expr*)u) }, 2));
    Expr* e_neg = subst_expr(e, absnode, negu);
    expr_free(negu);

    Expr* pos = expr_new_function(expr_new_symbol(SYM_And), (Expr*[]){
        rel2(SYM_GreaterEqual, expr_copy((Expr*)u), expr_new_integer(0)), e_pos }, 2);
    Expr* neg = expr_new_function(expr_new_symbol(SYM_And), (Expr*[]){
        rel2(SYM_Less, expr_copy((Expr*)u), expr_new_integer(0)), e_neg }, 2);
    Expr* both = expr_new_function(expr_new_symbol(SYM_Or), (Expr*[]){ pos, neg }, 2);

    /* recurse to clear the remaining Abs nodes */
    Expr* out = eliminate_abs(both, changed);
    expr_free(both);
    return out;
}

/* ------------------------------------------------------------------ *
 *  Preprocessing 4: Min/Max case-splitting                            *
 * ------------------------------------------------------------------ */

/* True iff `e` contains no Max/Min node anywhere. */
static bool is_minmaxfree(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return true;
    if (is_head(e, SYM_Max) || is_head(e, SYM_Min)) return false;
    if (!is_minmaxfree(e->data.function.head)) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!is_minmaxfree(e->data.function.args[i])) return false;
    return true;
}

/* Locate an innermost Max/Min node (all of whose arguments are Max/Min-free);
 * returns it (borrowed) or NULL. */
static const Expr* find_inner_minmax(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        const Expr* r = find_inner_minmax(e->data.function.args[i]);
        if (r) return r;
    }
    if ((is_head(e, SYM_Max) || is_head(e, SYM_Min))
        && e->data.function.arg_count >= 1
        && is_minmaxfree(e->data.function.head))
        return e;
    return NULL;
}

/* Eliminate one innermost Max/Min[a1,...,an] by branching on which argument is
 * the extremum: Max -> Or_i[ And[ a_i>=a_j (j!=i), e|_{node->a_i} ] ] (Min uses
 * a_i<=a_j).  Overlapping ties are harmless (both branches then substitute an
 * equal value), and every branch that asserts a_i as the extremum substitutes a
 * value that genuinely is one there, so no spurious solutions are introduced.
 * Recurses to clear nested Max/Min. */
static Expr* eliminate_minmax(const Expr* e, bool* changed) {
    const Expr* node = find_inner_minmax(e);
    if (!node) return expr_copy((Expr*)e);
    *changed = true;
    bool is_max = is_head(node, SYM_Max);
    const char* order = is_max ? SYM_GreaterEqual : SYM_LessEqual;
    size_t n = node->data.function.arg_count;

    Expr** branches = malloc(n * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) {
        const Expr* ai = node->data.function.args[i];
        Expr** conj = malloc(n * sizeof(Expr*));   /* up to n-1 constraints + stmt */
        size_t nc = 0;
        for (size_t j = 0; j < n; j++) {
            if (j == i) continue;
            conj[nc++] = rel2(order, expr_copy((Expr*)ai),
                              expr_copy((Expr*)node->data.function.args[j]));
        }
        conj[nc++] = subst_expr(e, node, ai);      /* e with this Max/Min -> a_i */
        branches[i] = expr_new_function(expr_new_symbol(SYM_And), conj, nc);
        free(conj);
    }
    Expr* both = expr_new_function(expr_new_symbol(SYM_Or), branches, n);
    free(branches);

    /* recurse to clear the remaining Max/Min nodes */
    Expr* out = eliminate_minmax(both, changed);
    expr_free(both);
    return out;
}

/* ------------------------------------------------------------------ *
 *  Preprocessing 5: general piecewise case-splitting                  *
 * ------------------------------------------------------------------ */

/* Heads this pass splits into (guard, value) clauses.  `Piecewise` is the
 * general container; the rest are fixed-shape piecewise builtins that reduce to
 * a clause list.  Abs and Max/Min have their own dedicated passes above. */
static bool is_piecewise_head(const char* h) {
    return h == SYM_Piecewise || h == SYM_Sign || h == SYM_UnitStep
        || h == SYM_Ramp || h == SYM_Clip || h == SYM_HeavisideTheta
        || h == SYM_Boole || h == SYM_UnitBox || h == SYM_IntegerPart
        || h == SYM_FractionalPart;
}

static bool node_is_piecewise(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && is_piecewise_head(e->data.function.head->data.symbol.name);
}

/* Innermost splittable piecewise node (deepest-first via arg recursion). */
static const Expr* find_inner_piecewise(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return NULL;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        const Expr* r = find_inner_piecewise(e->data.function.args[i]);
        if (r) return r;
    }
    if (node_is_piecewise(e)) return e;
    return NULL;
}

static void pw_push(Expr*** conds, Expr*** vals, int* n, int* cap, Expr* c, Expr* v) {
    if (*n == *cap) {
        *cap = *cap ? *cap * 2 : 4;
        *conds = realloc(*conds, (size_t)*cap * sizeof(Expr*));
        *vals  = realloc(*vals,  (size_t)*cap * sizeof(Expr*));
    }
    (*conds)[*n] = c; (*vals)[*n] = v; (*n)++;
}

/* u REL k (k an integer literal), unevaluated (copies u). */
static Expr* pw_rel_int(const char* rel, const Expr* u, long k) {
    return rel2(rel, expr_copy((Expr*)u), expr_new_integer(k));
}
/* u REL b (b an expression), unevaluated (copies both). */
static Expr* pw_rel(const char* rel, const Expr* u, const Expr* b) {
    return rel2(rel, expr_copy((Expr*)u), expr_copy((Expr*)b));
}
/* head[u] (Floor/Ceiling), unevaluated (copies u). */
static Expr* pw_call1(const char* head, const Expr* u) {
    return expr_new_function(expr_new_symbol(head), (Expr*[]){ expr_copy((Expr*)u) }, 1);
}
/* u - head[u]  (for FractionalPart), unevaluated (copies u). */
static Expr* pw_u_minus_call(const Expr* u, const char* head) {
    return expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){ expr_copy((Expr*)u),
        expr_new_function(expr_new_symbol(SYM_Times),
            (Expr*[]){ expr_new_integer(-1), pw_call1(head, u) }, 2) }, 2);
}
/* +-1/2 as a Rational, evaluated. */
static Expr* pw_half(int sign) {
    return enorm(expr_new_function(expr_new_symbol(SYM_Rational),
        (Expr*[]){ expr_new_integer(sign), expr_new_integer(2) }, 2));
}

/* Decompose a piecewise node into first-match (guard, value) clauses plus an
 * optional default value.  Returns false (allocating nothing) on a shape the
 * splitter does not model, so the node is left untouched and the diagram's
 * soundness gate can decline it.  On success the caller owns every returned
 * Expr and both arrays. */
static bool piecewise_clauses(const Expr* node, Expr*** conds_out, Expr*** vals_out,
                              int* nclauses, Expr** default_out, bool* has_default) {
    const char* h = node->data.function.head->data.symbol.name;
    size_t argc = node->data.function.arg_count;
    Expr** conds = NULL; Expr** vals = NULL; int n = 0, cap = 0;
    *default_out = NULL; *has_default = false;
    const Expr* u = argc >= 1 ? node->data.function.args[0] : NULL;

    if (h == SYM_Piecewise) {
        if (argc < 1 || !is_head(u, SYM_List)) return false;
        for (size_t i = 0; i < u->data.function.arg_count; i++) {
            const Expr* pr = u->data.function.args[i];
            if (!is_head(pr, SYM_List) || pr->data.function.arg_count != 2) {
                for (int t = 0; t < n; t++) { expr_free(conds[t]); expr_free(vals[t]); }
                free(conds); free(vals);
                return false;
            }
            pw_push(&conds, &vals, &n, &cap,
                    expr_copy(pr->data.function.args[1]),   /* condition */
                    expr_copy(pr->data.function.args[0]));  /* value     */
        }
        *default_out = (argc >= 2) ? expr_copy(node->data.function.args[1]) : expr_new_integer(0);
        *has_default = true;
    } else if (h == SYM_Sign && argc == 1) {
        pw_push(&conds, &vals, &n, &cap, pw_rel_int(SYM_Less, u, 0),    expr_new_integer(-1));
        pw_push(&conds, &vals, &n, &cap, pw_rel_int(SYM_Greater, u, 0), expr_new_integer(1));
        *default_out = expr_new_integer(0); *has_default = true;      /* u == 0 -> 0 */
    } else if (h == SYM_UnitStep && argc == 1) {
        pw_push(&conds, &vals, &n, &cap, pw_rel_int(SYM_Less, u, 0), expr_new_integer(0));
        *default_out = expr_new_integer(1); *has_default = true;      /* u >= 0 -> 1 */
    } else if (h == SYM_Ramp && argc == 1) {
        pw_push(&conds, &vals, &n, &cap, pw_rel_int(SYM_Less, u, 0), expr_new_integer(0));
        *default_out = expr_copy((Expr*)u); *has_default = true;      /* u >= 0 -> u */
    } else if (h == SYM_Boole && argc == 1) {
        pw_push(&conds, &vals, &n, &cap, expr_copy((Expr*)u), expr_new_integer(1));
        *default_out = expr_new_integer(0); *has_default = true;      /* !cond -> 0 */
    } else if (h == SYM_HeavisideTheta && argc == 1) {
        pw_push(&conds, &vals, &n, &cap, pw_rel_int(SYM_Less, u, 0),    expr_new_integer(0));
        pw_push(&conds, &vals, &n, &cap, pw_rel_int(SYM_Greater, u, 0), expr_new_integer(1));
        *has_default = false;                                        /* u == 0 excluded */
    } else if (h == SYM_UnitBox && argc == 1) {
        Expr* hn = pw_half(-1); Expr* hp = pw_half(1);
        pw_push(&conds, &vals, &n, &cap, pw_rel(SYM_Less, u, hn),    expr_new_integer(0));
        pw_push(&conds, &vals, &n, &cap, pw_rel(SYM_Greater, u, hp), expr_new_integer(0));
        expr_free(hn); expr_free(hp);
        *default_out = expr_new_integer(1); *has_default = true;     /* |u|<=1/2 -> 1 */
    } else if (h == SYM_Clip) {
        Expr *lo, *hi, *vlo, *vhi;
        if (argc == 1) { lo = expr_new_integer(-1); hi = expr_new_integer(1);
                         vlo = expr_new_integer(-1); vhi = expr_new_integer(1); }
        else if (argc >= 2 && is_head(node->data.function.args[1], SYM_List)
                 && node->data.function.args[1]->data.function.arg_count == 2) {
            const Expr* bnd = node->data.function.args[1];
            lo = expr_copy(bnd->data.function.args[0]);
            hi = expr_copy(bnd->data.function.args[1]);
            if (argc >= 3 && is_head(node->data.function.args[2], SYM_List)
                && node->data.function.args[2]->data.function.arg_count == 2) {
                const Expr* vb = node->data.function.args[2];
                vlo = expr_copy(vb->data.function.args[0]);
                vhi = expr_copy(vb->data.function.args[1]);
            } else { vlo = expr_copy(lo); vhi = expr_copy(hi); }
        } else return false;
        pw_push(&conds, &vals, &n, &cap, pw_rel(SYM_Less, u, lo),    vlo);   /* u<lo -> vlo */
        pw_push(&conds, &vals, &n, &cap, pw_rel(SYM_Greater, u, hi), vhi);   /* u>hi -> vhi */
        expr_free(lo); expr_free(hi);
        *default_out = expr_copy((Expr*)u); *has_default = true;             /* else -> u */
    } else if (h == SYM_IntegerPart && argc == 1) {
        pw_push(&conds, &vals, &n, &cap, pw_rel_int(SYM_Less, u, 0), pw_call1(SYM_Ceiling, u));
        *default_out = pw_call1(SYM_Floor, u); *has_default = true;   /* u>=0 -> Floor */
    } else if (h == SYM_FractionalPart && argc == 1) {
        pw_push(&conds, &vals, &n, &cap, pw_rel_int(SYM_Less, u, 0), pw_u_minus_call(u, SYM_Ceiling));
        *default_out = pw_u_minus_call(u, SYM_Floor); *has_default = true;
    } else {
        return false;
    }

    if (n == 0 && !*has_default) { free(conds); free(vals); return false; }
    *conds_out = conds; *vals_out = vals; *nclauses = n;
    return true;
}

/* Eliminate one innermost piecewise node by first-match case-split:
 * Or_i[ And[ Not[g_0],...,Not[g_{i-1}], g_i, e|_{node->v_i} ] ], plus a default
 * branch And[ Not[g_0],...,Not[g_{n-1}], e|_{node->def} ] when the head has a
 * default.  A head with no default (HeavisideTheta) leaves its uncovered region
 * (u==0) out of every branch, i.e. excluded from the solution set -- exactly the
 * points where the function is undefined.  Recurses to clear nesting. */
static Expr* eliminate_piecewise(const Expr* e, bool* changed) {
    const Expr* node = find_inner_piecewise(e);
    if (!node) return expr_copy((Expr*)e);

    Expr** conds = NULL; Expr** vals = NULL; int nc = 0;
    Expr* defval = NULL; bool hasdef = false;
    if (!piecewise_clauses(node, &conds, &vals, &nc, &defval, &hasdef))
        return expr_copy((Expr*)e);              /* unsupported shape: leave as-is */
    *changed = true;

    int nbranch = nc + (hasdef ? 1 : 0);
    Expr** branches = malloc((size_t)nbranch * sizeof(Expr*));
    for (int i = 0; i < nc; i++) {
        Expr** conj = malloc((size_t)(i + 2) * sizeof(Expr*));
        int m = 0;
        for (int j = 0; j < i; j++)
            conj[m++] = expr_new_function(expr_new_symbol(SYM_Not),
                          (Expr*[]){ expr_copy(conds[j]) }, 1);
        conj[m++] = expr_copy(conds[i]);
        conj[m++] = subst_expr(e, node, vals[i]);
        branches[i] = expr_new_function(expr_new_symbol(SYM_And), conj, m);
        free(conj);
    }
    if (hasdef) {
        Expr** conj = malloc((size_t)(nc + 1) * sizeof(Expr*));
        int m = 0;
        for (int j = 0; j < nc; j++)
            conj[m++] = expr_new_function(expr_new_symbol(SYM_Not),
                          (Expr*[]){ expr_copy(conds[j]) }, 1);
        conj[m++] = subst_expr(e, node, defval);
        branches[nc] = expr_new_function(expr_new_symbol(SYM_And), conj, m);
        free(conj);
    }
    Expr* both = expr_new_function(expr_new_symbol(SYM_Or), branches, nbranch);
    free(branches);
    for (int i = 0; i < nc; i++) { expr_free(conds[i]); expr_free(vals[i]); }
    free(conds); free(vals); if (defval) expr_free(defval);

    Expr* out = eliminate_piecewise(both, changed);
    expr_free(both);
    return out;
}

/* ------------------------------------------------------------------ *
 *  Driver                                                             *
 * ------------------------------------------------------------------ */

/* One pass of the domain- and variable-count-agnostic SELECTOR splits (Abs,
 * Min/Max, piecewise).  Shared by the univariate real-function driver and the
 * multivariate piecewise driver. */
static Expr* apply_selector_splits(const Expr* e, bool* changed) {
    Expr* a = eliminate_abs(e, changed);
    Expr* b = eliminate_minmax(a, changed);
    expr_free(a);
    Expr* c = eliminate_piecewise(b, changed);
    expr_free(b);
    return c;
}

/* ------------------------------------------------------------------ *
 *  Preprocessing 5: radical rationalization                          *
 *                                                                    *
 *  Rewrite a relation carrying a real square root Sqrt[u] (internally *
 *  Power[u, 1/2]) into an EQUIVALENT radical-free boolean combination *
 *  in the SAME variables, by isolating one radical and squaring under *
 *  exact sign guards.  Runs AFTER Abs/selector splitting (so radicands *
 *  are polynomial) and BEFORE atom canonicalisation (Together/Numerator *
 *  would mangle a surviving Sqrt), letting FM/CAD -- which accept only *
 *  polynomial atoms -- solve statements like Sqrt[Abs[x]]+Abs[y]<1.    *
 *                                                                    *
 *  Soundness over coverage: a relation whose radical cannot be        *
 *  isolated with a constant coefficient (or that has too many coupled  *
 *  radicals) is emitted UNCHANGED, so the downstream engine declines   *
 *  rather than risk a wrong answer.                                    */

#define RAT_RADICAL_DEPTH_CAP 12

/* Square-root radical node Power[u, 1/2]? */
static bool is_sqrt_radical(const Expr* e) {
    if (!is_head(e, SYM_Power) || e->data.function.arg_count != 2) return false;
    const Expr* ex = e->data.function.args[1];
    if (!is_head(ex, SYM_Rational) || ex->data.function.arg_count != 2) return false;
    const Expr* pn = ex->data.function.args[0];
    const Expr* pd = ex->data.function.args[1];
    return pn->type == EXPR_INTEGER && pn->data.integer == 1
        && pd->type == EXPR_INTEGER && pd->data.integer == 2;
}

/* Gather distinct (by structural equality) Sqrt-radical subterms of `e` whose
 * radicand mentions a reduce variable. */
static void collect_sqrt_radicals(const Expr* e, Expr** vars, int nv,
                                  const Expr*** arr, int* n, int* cap) {
    if (!e || e->type != EXPR_FUNCTION) return;
    if (is_sqrt_radical(e) && contains_any_var(e->data.function.args[0], vars, nv)) {
        bool seen = false;
        for (int i = 0; i < *n && !seen; i++) if (expr_eq((Expr*)(*arr)[i], (Expr*)e)) seen = true;
        if (!seen) {
            if (*n == *cap) { *cap = *cap ? *cap * 2 : 4; *arr = realloc(*arr, (size_t)*cap * sizeof(Expr*)); }
            (*arr)[(*n)++] = e;
        }
        /* still recurse: the radicand may hold nested radicals */
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        collect_sqrt_radicals(e->data.function.args[i], vars, nv, arr, n, cap);
}

/* The relation as seen after logical negation (De Morgan at a leaf). */
static const char* negate_rel(const char* h) {
    if (h == SYM_Less)         return SYM_GreaterEqual;
    if (h == SYM_LessEqual)    return SYM_Greater;
    if (h == SYM_Greater)      return SYM_LessEqual;
    if (h == SYM_GreaterEqual) return SYM_Less;
    if (h == SYM_Equal)        return SYM_Unequal;
    if (h == SYM_Unequal)      return SYM_Equal;
    return h;
}

static bool is_rel_head(const char* h) {
    return h == SYM_Less || h == SYM_LessEqual || h == SYM_Greater
        || h == SYM_GreaterEqual || h == SYM_Equal || h == SYM_Unequal;
}

/* a - b, evaluated. */
static Expr* mk_sub(const Expr* a, const Expr* b) {
    return enorm(expr_new_function(expr_new_symbol(SYM_Plus),
        (Expr*[]){ expr_copy((Expr*)a),
                   expr_new_function(expr_new_symbol(SYM_Times),
                       (Expr*[]){ expr_new_integer(-1), expr_copy((Expr*)b) }, 2) }, 2));
}
/* An (unevaluated) relation  a REL b, on fresh copies. */
static Expr* rel_copy(const char* h, const Expr* a, const Expr* b) {
    return rel2(h, expr_copy((Expr*)a), expr_copy((Expr*)b));
}
/* An (unevaluated) relation  a REL 0. */
static Expr* rel0(const char* h, const Expr* a) {
    return rel2(h, expr_copy((Expr*)a), expr_new_integer(0));
}
/* c^2, unevaluated (Power[c,2]). */
static Expr* mk_sq(const Expr* c) {
    return expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ expr_copy((Expr*)c), expr_new_integer(2) }, 2);
}
static Expr* mk_and2(Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(SYM_And), (Expr*[]){ a, b }, 2);
}
static Expr* mk_and3(Expr* a, Expr* b, Expr* c) {
    return expr_new_function(expr_new_symbol(SYM_And), (Expr*[]){ a, b, c }, 3);
}
static Expr* mk_or2(Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(SYM_Or), (Expr*[]){ a, b }, 2);
}

/* Exact radical-free equivalent of  Sqrt[u] REL c  over the Reals (u the
 * radicand, c the isolated other side).  Every row anchors u>=0; the != row is
 * derived directly (NOT as Not[==]): both == and != are False where u<0. */
static Expr* apply_radical_table(const char* rel, const Expr* u, const Expr* c) {
    Expr* dom = rel0(SYM_GreaterEqual, u);                 /* u >= 0 */
    if (rel == SYM_Less)                                   /* u>=0 && c>0  && u<c^2  */
        return mk_and3(dom, rel0(SYM_Greater, c),      rel_copy(SYM_Less, u, mk_sq(c)));
    if (rel == SYM_LessEqual)                             /* u>=0 && c>=0 && u<=c^2 */
        return mk_and3(dom, rel0(SYM_GreaterEqual, c), rel_copy(SYM_LessEqual, u, mk_sq(c)));
    if (rel == SYM_Greater)                              /* u>=0 && (c<0 || u>c^2) */
        return mk_and2(dom, mk_or2(rel0(SYM_Less, c), rel_copy(SYM_Greater, u, mk_sq(c))));
    if (rel == SYM_GreaterEqual)                         /* u>=0 && (c<0 || u>=c^2) */
        return mk_and2(dom, mk_or2(rel0(SYM_Less, c), rel_copy(SYM_GreaterEqual, u, mk_sq(c))));
    if (rel == SYM_Equal)                                /* u>=0 && c>=0 && u==c^2 */
        return mk_and3(dom, rel0(SYM_GreaterEqual, c), rel_copy(SYM_Equal, u, mk_sq(c)));
    /* SYM_Unequal:  u>=0 && (c<0 || u!=c^2)  */
    return mk_and2(dom, mk_or2(rel0(SYM_Less, c), rel_copy(SYM_Unequal, u, mk_sq(c))));
}

/* N[e] to a real; true iff e is a numeric constant (free of variables). */
static bool numeric_const(const Expr* e, double* out) {
    Expr* v = eval_and_free(expr_new_function(expr_new_symbol(SYM_N),
        (Expr*[]){ expr_copy((Expr*)e) }, 1));
    bool ok = false;
    if (v->type == EXPR_REAL)         { *out = v->data.real;          ok = true; }
    else if (v->type == EXPR_INTEGER) { *out = (double)v->data.integer; ok = true; }
    expr_free(v);
    return ok;
}

static Expr* rationalize_tree(const Expr* e, bool neg, Expr** vars, int nv,
                              int depth, bool* changed);

/* Rationalize a single relational leaf  L rel R  (rel already Not-folded). */
static Expr* rationalize_relation(const Expr* L, const Expr* R, const char* rel,
                                  Expr** vars, int nv, int depth, bool* changed) {
    Expr* diff = mk_sub(L, R);
    const Expr** rads = NULL; int nr = 0, rcap = 0;
    collect_sqrt_radicals(diff, vars, nv, &rads, &nr, &rcap);
    if (nr == 0 || depth >= RAT_RADICAL_DEPTH_CAP) {
        expr_free(diff); free(rads);
        return rel_copy(rel, L, R);                       /* nothing to do / cap: decline */
    }
    for (int i = 0; i < nr; i++) {
        const Expr* rho = rads[i];
        const Expr* u   = rho->data.function.args[0];     /* radicand */
        Expr* zero = expr_new_integer(0), *one = expr_new_integer(1), *two = expr_new_integer(2);
        Expr* B  = enorm(subst_expr(diff, rho, zero));    /* diff|_{rho=0}         */
        Expr* d1 = enorm(subst_expr(diff, rho, one));
        Expr* d2 = enorm(subst_expr(diff, rho, two));
        expr_free(zero); expr_free(one); expr_free(two);
        Expr* A = mk_sub(d1, B);                           /* diff|_{rho=1} - B      */
        /* linearity in rho: diff|_{rho=2} == 2A + B */
        Expr* twoAB = enorm(expr_new_function(expr_new_symbol(SYM_Plus),
            (Expr*[]){ expr_new_function(expr_new_symbol(SYM_Times),
                           (Expr*[]){ expr_new_integer(2), expr_copy(A) }, 2),
                       expr_copy(B) }, 2));
        bool lin = expr_eq(d2, twoAB);
        expr_free(twoAB); expr_free(d1); expr_free(d2);
        double av = 0.0;
        bool acon = numeric_const(A, &av);                 /* A a nonzero constant?  */
        if (!lin || !acon || av == 0.0) { expr_free(A); expr_free(B); continue; }
        /* c = -B/A ; flip the relation when A < 0. */
        Expr* c = enorm(expr_new_function(expr_new_symbol(SYM_Times),
            (Expr*[]){ expr_new_integer(-1), expr_copy(B),
                       expr_new_function(expr_new_symbol(SYM_Power),
                           (Expr*[]){ expr_copy(A), expr_new_integer(-1) }, 2) }, 3));
        const char* r2 = (av < 0.0) ? flip_rel(rel) : rel;
        Expr* table = apply_radical_table(r2, u, c);
        expr_free(A); expr_free(B); expr_free(c);
        expr_free(diff); free(rads);
        *changed = true;
        Expr* out = rationalize_tree(table, false, vars, nv, depth + 1, changed);
        expr_free(table);
        return out;
    }
    expr_free(diff); free(rads);
    return rel_copy(rel, L, R);                            /* no constant-coeff radical */
}

/* NNF walk carrying polarity `neg`, so no Not ever survives around a rationalized
 * u>=0 guard (which De Morgan would wrongly re-open onto the u<0 region). */
static Expr* rationalize_tree(const Expr* e, bool neg, Expr** vars, int nv,
                              int depth, bool* changed) {
    if (!e) return NULL;
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        size_t n = e->data.function.arg_count;
        if (is_rel_head(h) && n == 2) {
            const char* rel = neg ? negate_rel(h) : h;
            return rationalize_relation(e->data.function.args[0], e->data.function.args[1],
                                        rel, vars, nv, depth, changed);
        }
        if ((h == SYM_And || h == SYM_Or) && n >= 1) {
            const char* out = (h == SYM_And) == !neg ? SYM_And : SYM_Or;  /* De Morgan */
            Expr** parts = malloc(n * sizeof(Expr*));
            for (size_t i = 0; i < n; i++)
                parts[i] = rationalize_tree(e->data.function.args[i], neg, vars, nv, depth, changed);
            Expr* r = expr_new_function(expr_new_symbol(out), parts, n);
            free(parts);
            return r;
        }
        if (h == SYM_Not && n == 1)
            return rationalize_tree(e->data.function.args[0], !neg, vars, nv, depth, changed);
        if (h == SYM_Implies && n == 2) {                  /* a => b  ==  !a || b */
            Expr* na = expr_new_function(expr_new_symbol(SYM_Not),
                (Expr*[]){ expr_copy(e->data.function.args[0]) }, 1);
            Expr* orx = mk_or2(na, expr_copy(e->data.function.args[1]));
            Expr* r = rationalize_tree(orx, neg, vars, nv, depth, changed);
            expr_free(orx);
            return r;
        }
        if (h == SYM_Inequality && n >= 3 && (n % 2) == 1) { /* chained a op b op c ... */
            size_t nrel = (n - 1) / 2;
            Expr** parts = malloc(nrel * sizeof(Expr*));
            for (size_t j = 0; j < nrel; j++) {
                const Expr* a  = e->data.function.args[2 * j];
                const Expr* op = e->data.function.args[2 * j + 1];
                const Expr* b  = e->data.function.args[2 * j + 2];
                const char* oh = op->type == EXPR_SYMBOL ? op->data.symbol.name : SYM_Less;
                Expr* leaf = rel_copy(oh, a, b);
                parts[j] = rationalize_tree(leaf, neg, vars, nv, depth, changed);
                expr_free(leaf);
            }
            Expr* r = expr_new_function(expr_new_symbol(neg ? SYM_Or : SYM_And), parts, nrel);
            free(parts);
            return r;
        }
    }
    /* Non-boolean / unhandled node: it carries no radical relation we rewrite, so
     * copying it (Not-wrapped under negation) is sound. */
    if (neg) return expr_new_function(expr_new_symbol(SYM_Not), (Expr*[]){ expr_copy((Expr*)e) }, 1);
    return expr_copy((Expr*)e);
}

/* One pass: rationalize every square-root-bearing relational leaf. */
static Expr* rationalize_radical_leaves(const Expr* e, Expr** vars, int nv, bool* changed) {
    return rationalize_tree(e, false, vars, nv, 0, changed);
}

/* Rewrite Mod->Floor, Abs sign-splits, Min/Max case-splits and integer-part
 * relations away, so the sign-diagram engines see only polynomial atoms.  The
 * four transforms are cyclically dependent -- an Abs or Min/Max split can EXPOSE
 * a fresh Floor relational leaf for the integer-part pass, and an integer-part
 * rewrite of Floor[Max[...]] can expose a Min/Max -- so a single pass in any
 * fixed order is incomplete.  Iterate to a fixpoint (bounded; each transform
 * strictly removes an occurrence of its target construct without reintroducing
 * another's). */
Expr* reduce_realfn_preprocess(const Expr* e, const Expr* x, bool* changed) {
    bool any = false;
    Expr* cur = expr_copy((Expr*)e);
    for (int iter = 0; iter < 8; iter++) {
        bool ch = false;
        Expr* a = subst_mod(cur, &ch);              /* Mod -> Floor           */
        Expr* b = apply_selector_splits(a, &ch);    /* Abs / Min-Max / Piecew.*/
        expr_free(a);
        Expr* d = rewrite_ip_leaves(b, x, &ch);     /* integer-part relations */
        expr_free(b);
        expr_free(cur);
        cur = d;
        if (!ch) break;
        any = true;
    }
    if (changed) *changed = any;
    if (!any) { expr_free(cur); return NULL; }
    return cur;
}

/* Multivariate (any nv) piecewise preprocessing: the domain-agnostic selector
 * splits only (Abs, Min/Max, Piecewise/Sign/UnitStep/...).  The integer-part
 * machinery (Mod->Floor and Floor/Ceiling/Round isolation) is univariate, so a
 * residual integer-part atom is left for the CAD engine to decline soundly. */
Expr* reduce_piecewise_preprocess(const Expr* e, Expr** vars, int nv, bool* changed) {
    bool any = false;
    Expr* cur = expr_copy((Expr*)e);
    for (int iter = 0; iter < 8; iter++) {
        bool ch = false;
        Expr* nxt = apply_selector_splits(cur, &ch);       /* Abs / Min-Max / Piecew. */
        Expr* rad = rationalize_radical_leaves(nxt, vars, nv, &ch);  /* Sqrt[u] -> u<c^2 */
        expr_free(nxt);
        expr_free(cur);
        cur = rad;
        if (!ch) break;
        any = true;
    }
    if (changed) *changed = any;
    if (!any) { expr_free(cur); return NULL; }
    return cur;
}
