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
        || h == SYM_ArcCosh || h == SYM_ArcSech || h == SYM_Sqrt) {
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
        if (rel == SYM_Equal) {
            if (!c_int) { out = expr_new_symbol(SYM_False); }
            else out = expr_new_function(expr_new_symbol(SYM_Inequality),
                (Expr*[]){ plus_int(c, -1), expr_new_symbol(SYM_Less),
                           expr_copy((Expr*)v), expr_new_symbol(SYM_LessEqual),
                           expr_copy((Expr*)c) }, 5);
        }
        /* other Ceiling relations: decline (not needed by the corpus). */
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

/* Walk the logical tree; rewrite integer-part relational leaves. */
static Expr* rewrite_ip_leaves(const Expr* e, bool* changed) {
    if (!e || e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    const char* h = e->data.function.head->type == EXPR_SYMBOL
                  ? e->data.function.head->data.symbol.name : NULL;
    bool logical = (h == SYM_And || h == SYM_Or || h == SYM_Not
                    || h == SYM_Implies || h == SYM_Xor || h == SYM_Inequality);
    bool relational = (h == SYM_Equal || h == SYM_Unequal || h == SYM_Less
                       || h == SYM_LessEqual || h == SYM_Greater || h == SYM_GreaterEqual);

    if (relational) {
        Expr* iso = try_ip_isolate(e, changed);
        if (iso) return iso;
        return expr_copy((Expr*)e);
    }
    if (logical) {
        size_t n = e->data.function.arg_count;
        Expr** args = (n > 0) ? malloc(n * sizeof(Expr*)) : NULL;
        for (size_t i = 0; i < n; i++) args[i] = rewrite_ip_leaves(e->data.function.args[i], changed);
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
 *  Driver                                                             *
 * ------------------------------------------------------------------ */

Expr* reduce_realfn_preprocess(const Expr* e, const Expr* x, bool* changed) {
    (void)x;
    bool ch = false;
    Expr* a = subst_mod(e, &ch);            /* 1. Mod -> Floor        */
    Expr* b = rewrite_ip_leaves(a, &ch);    /* 2. integer-part isolate */
    expr_free(a);
    Expr* c = eliminate_abs(b, &ch);        /* 3. Abs sign-split       */
    expr_free(b);
    if (changed) *changed = ch;
    if (!ch) { expr_free(c); return NULL; }
    return c;
}
