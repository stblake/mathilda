/*
 * reduce_trigregion.c -- bounded periodic equation + inequality conjunction.
 * See reduce_trigregion.h for the contract and method.
 */
#include "reduce_trigregion.h"

#include <stdbool.h>
#include <stdlib.h>
#include <math.h>

#include "eval.h"
#include "expr.h"
#include "reduce_opts.h"
#include "sym_names.h"
#include "zero_test.h"

/* Guard against a mis-parsed / pathological interval becoming a runaway loop.
 * A well-posed bounded region selects a handful of members per family. */
#define TRIGREGION_MAX_MEMBERS 100000

/* ------------------------------------------------------------------ *
 *  Small builders / accessors                                         *
 * ------------------------------------------------------------------ */

static Expr* mk_fn1(const char* head_sym, Expr* a) {
    return expr_new_function(expr_new_symbol(head_sym), (Expr*[]){ a }, 1);
}
static Expr* mk_fn2(const char* head_sym, Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(head_sym), (Expr*[]){ a, b }, 2);
}

static bool head_is(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == name;
}

static bool is_true_sym(const Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_True;
}
static bool is_false_sym(const Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol.name == SYM_False;
}

/* Does `e` contain the symbol `var` anywhere? (interned-name compare) */
static bool tr_contains_var(const Expr* e, const Expr* var) {
    if (!e || !var || var->type != EXPR_SYMBOL) return false;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == var->data.symbol.name;
    if (e->type != EXPR_FUNCTION) return false;
    if (tr_contains_var(e->data.function.head, var)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (tr_contains_var(e->data.function.args[i], var)) return true;
    return false;
}

/* ------------------------------------------------------------------ *
 *  Reality + numeric helpers                                          *
 * ------------------------------------------------------------------ */

/* Im[e] provably zero? */
static bool tr_im_is_zero(const Expr* e) {
    Expr* im = eval_and_free(mk_fn1(SYM_Im, expr_copy((Expr*)e)));
    bool z = (zero_test_decide(im) == ZERO_TEST_TRUE);
    expr_free(im);
    return z;
}

/* Evaluate N[e]; store a finite machine value in *out.  false if N[e] is not a
 * plain finite real (a symbol, Infinity, complex, ...). */
static bool tr_to_double(const Expr* e, double* out) {
    Expr* n = eval_and_free(mk_fn1(SYM_N, expr_copy((Expr*)e)));
    bool ok = true;
    if (n->type == EXPR_REAL)         *out = n->data.real;
    else if (n->type == EXPR_INTEGER) *out = (double)n->data.integer;
    else ok = false;
    expr_free(n);
    if (ok && !isfinite(*out)) ok = false;
    return ok;
}

/* stmt /. var -> val   (owned result).  `stmt` and `var` are borrowed; `val` is
 * ADOPTED (freed here) so callers pass a freshly built value directly and a
 * value they must retain as an explicit expr_copy. */
static Expr* tr_subst(const Expr* stmt, const Expr* var, Expr* val) {
    Expr* rule = mk_fn2(SYM_Rule, expr_copy((Expr*)var), val);   /* adopts val */
    return eval_and_free(mk_fn2(SYM_ReplaceAll, expr_copy((Expr*)stmt), rule));
}

/* stmt /. var -> cand  evaluates to exactly True? (`cand` borrowed.) */
static bool tr_verify(const Expr* stmt, const Expr* var, const Expr* cand) {
    Expr* r = tr_subst(stmt, var, expr_copy((Expr*)cand));
    bool ok = is_true_sym(r);
    expr_free(r);
    return ok;
}

/* ------------------------------------------------------------------ *
 *  Feasible-interval extraction (over the real parameter k)           *
 * ------------------------------------------------------------------ */

typedef struct { double lo, hi; } Ivl;
typedef struct { Ivl* a; size_t n, cap; } IvlVec;

static void ivl_push(IvlVec* v, double lo, double hi) {
    if (v->n == v->cap) {
        v->cap = v->cap ? v->cap * 2 : 4;
        v->a = (Ivl*)realloc(v->a, v->cap * sizeof(Ivl));
    }
    v->a[v->n].lo = lo;
    v->a[v->n].hi = hi;
    v->n++;
}

/* Tighten (*has_lo,*lo,*has_hi,*hi) with one relation atom over `k`.
 * Returns false if the atom cannot be interpreted as a bound on `k`. */
static bool tr_atom_bound(const Expr* atom, const Expr* k,
                          bool* has_lo, double* lo,
                          bool* has_hi, double* hi) {
    /* Inequality[b0, o0, k, o1, b1] with k in the middle. */
    if (head_is(atom, SYM_Inequality) && atom->data.function.arg_count == 5) {
        if (!expr_eq(atom->data.function.args[2], k)) return false;
        double b0, b1;
        if (!tr_to_double(atom->data.function.args[0], &b0)) return false;
        if (!tr_to_double(atom->data.function.args[4], &b1)) return false;
        if (!*has_lo || b0 > *lo) { *has_lo = true; *lo = b0; }
        if (!*has_hi || b1 < *hi) { *has_hi = true; *hi = b1; }
        return true;
    }
    if (!(atom && atom->type == EXPR_FUNCTION
          && atom->data.function.arg_count == 2)) return false;
    const char* h = atom->data.function.head->type == EXPR_SYMBOL
                        ? atom->data.function.head->data.symbol.name : NULL;
    if (!h) return false;
    const Expr* L = atom->data.function.args[0];
    const Expr* R = atom->data.function.args[1];
    bool kL = expr_eq(L, k), kR = expr_eq(R, k);
    if (kL == kR) return false;   /* k must appear on exactly one side */
    double b;
    if (kL) {                     /* k  op  R  */
        if (!tr_to_double(R, &b)) return false;
        if (h == SYM_Less || h == SYM_LessEqual) {
            if (!*has_hi || b < *hi) { *has_hi = true; *hi = b; }
        } else if (h == SYM_Greater || h == SYM_GreaterEqual) {
            if (!*has_lo || b > *lo) { *has_lo = true; *lo = b; }
        } else if (h == SYM_Equal) {
            *has_lo = *has_hi = true; *lo = *hi = b;
        } else return false;
    } else {                      /* L  op  k  */
        if (!tr_to_double(L, &b)) return false;
        if (h == SYM_Less || h == SYM_LessEqual) {
            if (!*has_lo || b > *lo) { *has_lo = true; *lo = b; }
        } else if (h == SYM_Greater || h == SYM_GreaterEqual) {
            if (!*has_hi || b < *hi) { *has_hi = true; *hi = b; }
        } else if (h == SYM_Equal) {
            *has_lo = *has_hi = true; *lo = *hi = b;
        } else return false;
    }
    return true;
}

/* A conjunctive region (a single atom or And of atoms) -> one bounded [lo,hi].
 * Returns false if any atom is uninterpretable or a side is left unbounded. */
static bool tr_conj_interval(const Expr* region, const Expr* k, Ivl* out) {
    bool has_lo = false, has_hi = false;
    double lo = 0, hi = 0;
    if (head_is(region, SYM_And)) {
        for (size_t i = 0; i < region->data.function.arg_count; i++)
            if (!tr_atom_bound(region->data.function.args[i], k,
                               &has_lo, &lo, &has_hi, &hi))
                return false;
    } else {
        if (!tr_atom_bound(region, k, &has_lo, &lo, &has_hi, &hi))
            return false;
    }
    if (!has_lo || !has_hi) return false;   /* unbounded */
    out->lo = lo; out->hi = hi;
    return true;
}

/* Parse `regionC` (Reduce[..., k, Reals]) into bounded intervals in `out`.
 * Returns false to DECLINE (True/unbounded/unrecognised); true otherwise, with
 * `out` possibly empty (a False region contributes no members). */
static bool tr_region_intervals(const Expr* regionC, const Expr* k, IvlVec* out) {
    if (is_true_sym(regionC))  return false;   /* all reals -> unbounded */
    if (is_false_sym(regionC)) return true;    /* empty */
    if (head_is(regionC, SYM_Or)) {
        for (size_t i = 0; i < regionC->data.function.arg_count; i++)
            if (!tr_region_intervals(regionC->data.function.args[i], k, out))
                return false;
        return true;
    }
    Ivl iv;
    if (!tr_conj_interval(regionC, k, &iv)) return false;
    if (iv.lo <= iv.hi) ivl_push(out, iv.lo, iv.hi);
    return true;
}

/* ------------------------------------------------------------------ *
 *  Candidate collection (dedup + numeric sort)                        *
 * ------------------------------------------------------------------ */

typedef struct { Expr* e; double key; } Cand;
typedef struct { Cand* a; size_t n, cap; } CandVec;

static void cand_add(CandVec* v, Expr* e /* adopted */) {
    for (size_t i = 0; i < v->n; i++)
        if (expr_eq(v->a[i].e, e)) { expr_free(e); return; }   /* dedup */
    double key;
    if (!tr_to_double(e, &key)) key = 0.0;   /* unorderable -> stable-ish */
    if (v->n == v->cap) {
        v->cap = v->cap ? v->cap * 2 : 8;
        v->a = (Cand*)realloc(v->a, v->cap * sizeof(Cand));
    }
    v->a[v->n].e = e;
    v->a[v->n].key = key;
    v->n++;
}

static int cand_cmp(const void* pa, const void* pb) {
    double a = ((const Cand*)pa)->key, b = ((const Cand*)pb)->key;
    return (a < b) ? -1 : (a > b) ? 1 : 0;
}

/* ------------------------------------------------------------------ *
 *  Driver                                                             *
 * ------------------------------------------------------------------ */

Expr* reduce_periodic_region(const Expr* peq, const Expr* region_stmt,
                             const Expr* var, const ReduceOpts* opts) {
    if (!peq || !region_stmt || !var) return NULL;

    /* Solve the periodic equation over the Reals into its families. */
    Expr* base[3];
    base[0] = expr_copy((Expr*)peq);
    base[1] = expr_copy((Expr*)var);
    base[2] = expr_new_symbol(SYM_Reals);
    Expr* sols = eval_and_free(reduce_opts_build_solve(base, 3, opts));
    if (!head_is(sols, SYM_List)) { expr_free(sols); return NULL; }

    /* The full statement, for the exact re-check of every candidate. */
    Expr* full = mk_fn2(SYM_And, expr_copy((Expr*)peq),
                                 expr_copy((Expr*)region_stmt));
    /* Fresh plain symbol standing in for the family parameter C[k] (an indexed
     * C[k] is not a valid Reduce variable). */
    Expr* ksym = expr_new_symbol("$ReduceRegionK");

    CandVec cands = {0};
    bool declined = false;

    for (size_t i = 0; i < sols->data.function.arg_count && !declined; i++) {
        Expr* row = sols->data.function.args[i];
        if (!head_is(row, SYM_List) || row->data.function.arg_count != 1) {
            declined = true; break;
        }
        Expr* rule = row->data.function.args[0];
        if (!head_is(rule, SYM_Rule) || rule->data.function.arg_count != 2
            || !expr_eq(rule->data.function.args[0], var)) {
            declined = true; break;
        }
        Expr* val = rule->data.function.args[1];

        /* Peel ConditionalExpression[value, Element[param, Integers]]. */
        const Expr* value = val;
        const Expr* param = NULL;
        if (head_is(val, SYM_ConditionalExpression)
            && val->data.function.arg_count == 2) {
            value = val->data.function.args[0];
            Expr* cond = val->data.function.args[1];
            if (head_is(cond, SYM_Element) && cond->data.function.arg_count == 2
                && cond->data.function.args[1]->type == EXPR_SYMBOL
                && cond->data.function.args[1]->data.symbol.name == SYM_Integers) {
                param = cond->data.function.args[0];
            }
        }

        if (!param) {
            /* Isolated solution: verify directly against the full statement. */
            Expr* c = eval_and_free(expr_copy((Expr*)value));
            if (tr_verify(full, var, c)) cand_add(&cands, c);
            else expr_free(c);
            continue;
        }

        /* Family value = a + p*param.  period p = Coefficient[value, param];
         * offset a = value /. param -> 0. */
        Expr* p = eval_and_free(mk_fn2(SYM_Coefficient,
                        expr_copy((Expr*)value), expr_copy((Expr*)param)));
        Expr* a = tr_subst(value, param, expr_new_integer(0));
        bool period_real = tr_im_is_zero(p);
        bool offset_real = tr_im_is_zero(a);
        expr_free(p);

        if (!period_real) {
            /* Non-real period (e.g. Sinh/Exp's 2 I Pi C[k]): the only candidate
             * real member is the offset a (the k=0 member).  Verify it against
             * the full statement (which enforces the region). */
            if (offset_real) {
                Expr* c = eval_and_free(a);   /* adopts a */
                if (tr_verify(full, var, c)) cand_add(&cands, c);
                else expr_free(c);
            } else {
                expr_free(a);
            }
            continue;
        }
        expr_free(a);
        if (!offset_real) continue;   /* real period, complex offset: no members */

        /* Real period: valK = value /. param -> ksym (a plain-symbol parameter),
         * then bound k by reducing the induced constraint over the reals. */
        Expr* valK = tr_subst(value, param, expr_copy(ksym));
        Expr* constraintC = tr_subst(region_stmt, var, expr_copy(valK));
        Expr* regionC = eval_and_free(expr_new_function(
            expr_new_symbol(SYM_Reduce),
            (Expr*[]){ constraintC, expr_copy(ksym), expr_new_symbol(SYM_Reals) }, 3));

        IvlVec iv = {0};
        bool ok = tr_region_intervals(regionC, ksym, &iv);
        expr_free(regionC);
        if (!ok) { free(iv.a); expr_free(valK); declined = true; break; }

        for (size_t j = 0; j < iv.n && !declined; j++) {
            double clo = iv.a[j].lo, chi = iv.a[j].hi;
            /* Widen by one on each side so a borderline family member is never
             * missed; the exact re-check discards any that do not belong. */
            long kmin = (long)floor(clo) - 1;
            long kmax = (long)ceil(chi) + 1;
            if (kmax < kmin || (kmax - kmin) > TRIGREGION_MAX_MEMBERS) {
                declined = true; break;
            }
            for (long kk = kmin; kk <= kmax; kk++) {
                Expr* cand = tr_subst(valK, ksym, expr_new_integer((int64_t)kk));
                if (tr_verify(full, var, cand)) cand_add(&cands, cand);
                else expr_free(cand);
            }
        }
        free(iv.a);
        expr_free(valK);
    }

    expr_free(ksym);
    expr_free(full);
    expr_free(sols);

    if (declined) {
        for (size_t i = 0; i < cands.n; i++) expr_free(cands.a[i].e);
        free(cands.a);
        return NULL;
    }
    if (cands.n == 0) { free(cands.a); return expr_new_symbol(SYM_False); }

    qsort(cands.a, cands.n, sizeof(Cand), cand_cmp);

    Expr* out;
    if (cands.n == 1) {
        out = mk_fn2(SYM_Equal, expr_copy((Expr*)var), cands.a[0].e);
    } else {
        Expr** terms = (Expr**)malloc(cands.n * sizeof(Expr*));
        for (size_t i = 0; i < cands.n; i++)
            terms[i] = mk_fn2(SYM_Equal, expr_copy((Expr*)var), cands.a[i].e);
        out = expr_new_function(expr_new_symbol(SYM_Or), terms, cands.n);
        free(terms);
    }
    free(cands.a);   /* the Expr*s were adopted by the Equal[] nodes */
    return eval_and_free(out);
}

/* ================================================================== *
 *  Trig / hyperbolic INEQUALITY over a bounded region                 *
 *  (Sin[x] > 1/2 && 0 < x < 2 Pi  ->  Pi/6 < x < 5 Pi/6)              *
 * ================================================================== */

static bool is_pole_free_trig(const char* h) {
    return h == SYM_Sin || h == SYM_Cos || h == SYM_Sinh
        || h == SYM_Cosh || h == SYM_Tanh || h == SYM_Sech;
}
static bool is_pole_bearing_trig(const char* h) {
    return h == SYM_Tan || h == SYM_Cot || h == SYM_Sec
        || h == SYM_Csc || h == SYM_Coth || h == SYM_Csch;
}

/* Does `e` contain any trig/hyperbolic head over a `var`-bearing argument? */
static bool tr_has_trig_over_var(const Expr* e, const Expr* var) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL && e->data.function.arg_count == 1) {
        const char* hn = h->data.symbol.name;
        if ((is_pole_free_trig(hn) || is_pole_bearing_trig(hn))
            && tr_contains_var(e->data.function.args[0], var))
            return true;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (tr_has_trig_over_var(e->data.function.args[i], var)) return true;
    return false;
}

/* Any pole-bearing trig head over `var`? (Tan/Cot/Sec/Csc/Coth/Csch.) */
static bool tr_has_pole_bearing(const Expr* e, const Expr* var) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type == EXPR_SYMBOL && e->data.function.arg_count == 1
        && is_pole_bearing_trig(h->data.symbol.name)
        && tr_contains_var(e->data.function.args[0], var))
        return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (tr_has_pole_bearing(e->data.function.args[i], var)) return true;
    return false;
}

/* A binary relation atom (Less/…/Equal/Unequal) — its residual arg0-arg1. */
static bool is_rel_atom(const Expr* a) {
    if (!a || a->type != EXPR_FUNCTION || a->data.function.arg_count != 2
        || a->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = a->data.function.head->data.symbol.name;
    return h == SYM_Less || h == SYM_LessEqual || h == SYM_Greater
        || h == SYM_GreaterEqual || h == SYM_Equal || h == SYM_Unequal;
}

static Expr* rel_residual(const Expr* atom) {
    return eval_and_free(mk_fn2(SYM_Plus,
        expr_copy(atom->data.function.args[0]),
        expr_new_function(expr_new_symbol(SYM_Times),
            (Expr*[]){ expr_new_integer(-1),
                       expr_copy(atom->data.function.args[1]) }, 2)));
}

/* ---- exact bounded-interval parsing of a Reduce[…, var, Reals] region --- */

typedef struct { Expr* lo; Expr* hi; } EIvl;   /* owned */
typedef struct { EIvl* a; size_t n, cap; } EIvlVec;

static void eivl_push(EIvlVec* v, Expr* lo, Expr* hi) {
    if (v->n == v->cap) {
        v->cap = v->cap ? v->cap * 2 : 4;
        v->a = (EIvl*)realloc(v->a, v->cap * sizeof(EIvl));
    }
    v->a[v->n].lo = lo; v->a[v->n].hi = hi; v->n++;
}

/* Read the (borrowed) lower/upper bound Exprs a single atom places on `var`. */
static bool tr_atom_exact_bound(const Expr* atom, const Expr* var,
                                const Expr** lo, const Expr** hi) {
    if (head_is(atom, SYM_Inequality) && atom->data.function.arg_count == 5) {
        if (!expr_eq(atom->data.function.args[2], var)) return false;
        *lo = atom->data.function.args[0];
        *hi = atom->data.function.args[4];
        return true;
    }
    if (!is_rel_atom(atom)) return false;
    const char* h = atom->data.function.head->data.symbol.name;
    const Expr* L = atom->data.function.args[0];
    const Expr* R = atom->data.function.args[1];
    bool kL = expr_eq(L, var), kR = expr_eq(R, var);
    if (kL == kR) return false;
    const Expr* b = kL ? R : L;
    /* Orientation: with var on the left, `<`/`<=` bound above and `>`/`>=`
     * below; with var on the right the roles flip. */
    bool upper = kL ? (h == SYM_Less || h == SYM_LessEqual)
                    : (h == SYM_Greater || h == SYM_GreaterEqual);
    bool lower = kL ? (h == SYM_Greater || h == SYM_GreaterEqual)
                    : (h == SYM_Less || h == SYM_LessEqual);
    if (h == SYM_Equal) { *lo = b; *hi = b; return true; }
    if (upper) { *hi = b; return true; }
    if (lower) { *lo = b; return true; }
    return false;   /* Unequal or unhandled */
}

/* Keep the tighter of two lower bounds (larger value); NULL means "unset". */
static const Expr* tighter_lo(const Expr* cur, const Expr* cand) {
    if (!cand) return cur;
    if (!cur) return cand;
    double dc, dn;
    if (tr_to_double(cand, &dn) && tr_to_double(cur, &dc) && dn > dc) return cand;
    return cur;
}
static const Expr* tighter_hi(const Expr* cur, const Expr* cand) {
    if (!cand) return cur;
    if (!cur) return cand;
    double dc, dn;
    if (tr_to_double(cand, &dn) && tr_to_double(cur, &dc) && dn < dc) return cand;
    return cur;
}

/* Tighten (*lo,*hi) with a conjunctive region (a lone atom or And of atoms). */
static bool tr_conj_exact(const Expr* region, const Expr* var,
                          const Expr** lo, const Expr** hi) {
    bool is_and = head_is(region, SYM_And);
    size_t n = is_and ? region->data.function.arg_count : 1;
    for (size_t i = 0; i < n; i++) {
        const Expr* atom = is_and ? region->data.function.args[i] : region;
        const Expr *alo = NULL, *ahi = NULL;
        if (!tr_atom_exact_bound(atom, var, &alo, &ahi)) return false;
        *lo = tighter_lo(*lo, alo);
        *hi = tighter_hi(*hi, ahi);
    }
    return (*lo != NULL && *hi != NULL);
}

/* Parse a Reduce region into bounded exact intervals.  false => decline. */
static bool tr_region_exact(const Expr* region, const Expr* var, EIvlVec* out) {
    if (is_true_sym(region))  return false;   /* unbounded */
    if (is_false_sym(region)) return true;    /* empty */
    if (head_is(region, SYM_Or)) {
        for (size_t i = 0; i < region->data.function.arg_count; i++)
            if (!tr_region_exact(region->data.function.args[i], var, out))
                return false;
        return true;
    }
    const Expr *lo = NULL, *hi = NULL;
    if (!tr_conj_exact(region, var, &lo, &hi)) return false;
    double dlo, dhi;
    if (!tr_to_double(lo, &dlo) || !tr_to_double(hi, &dhi)) return false;
    eivl_push(out, expr_copy((Expr*)lo), expr_copy((Expr*)hi));
    return true;
}

/* Collect the real zeros of `residual` strictly inside (lo,hi) into `out`.
 * false => Solve declined (caller declines the whole problem). */
static bool tr_collect_zeros(const Expr* residual, const Expr* var,
                             double lo, double hi, const ReduceOpts* opts,
                             CandVec* out) {
    Expr* base[3];
    base[0] = mk_fn2(SYM_Equal, expr_copy((Expr*)residual), expr_new_integer(0));
    base[1] = expr_copy((Expr*)var);
    base[2] = expr_new_symbol(SYM_Reals);
    Expr* sols = eval_and_free(reduce_opts_build_solve(base, 3, opts));
    if (!head_is(sols, SYM_List)) { expr_free(sols); return false; }

    bool ok = true;
    for (size_t i = 0; i < sols->data.function.arg_count && ok; i++) {
        Expr* row = sols->data.function.args[i];
        if (!head_is(row, SYM_List) || row->data.function.arg_count != 1) continue;
        Expr* rule = row->data.function.args[0];
        if (!head_is(rule, SYM_Rule) || rule->data.function.arg_count != 2
            || !expr_eq(rule->data.function.args[0], var)) continue;
        Expr* val = rule->data.function.args[1];
        const Expr* value = val;
        const Expr* param = NULL;
        if (head_is(val, SYM_ConditionalExpression)
            && val->data.function.arg_count == 2) {
            value = val->data.function.args[0];
            Expr* cond = val->data.function.args[1];
            if (head_is(cond, SYM_Element) && cond->data.function.arg_count == 2
                && cond->data.function.args[1]->type == EXPR_SYMBOL
                && cond->data.function.args[1]->data.symbol.name == SYM_Integers)
                param = cond->data.function.args[0];
        }
        if (!param) {                       /* isolated root */
            double d;
            if (tr_to_double(value, &d) && d > lo && d < hi)
                cand_add(out, eval_and_free(expr_copy((Expr*)value)));
            continue;
        }
        Expr* p = eval_and_free(mk_fn2(SYM_Coefficient,
                        expr_copy((Expr*)value), expr_copy((Expr*)param)));
        Expr* a = tr_subst(value, param, expr_new_integer(0));
        /* N[·] returns a plain real only for a real value (a complex period or
         * offset -- e.g. Sinh/Exp's 2 I Pi -- yields a Complex head), so the
         * numeric test doubles as the reality test without a zero_test call. */
        double dp, da;
        bool prd = tr_to_double(p, &dp) && dp != 0.0;   /* real nonzero period */
        bool ard = tr_to_double(a, &da);                /* real offset         */
        expr_free(p);
        if (!prd) {                          /* non-real period: k=0 member */
            if (ard && da > lo && da < hi) cand_add(out, a);
            else expr_free(a);
            continue;
        }
        expr_free(a);
        if (!ard) continue;
        /* real period dp, offset da: enumerate integers with member in (lo,hi). */
        double t0 = (lo - da) / dp, t1 = (hi - da) / dp;
        long kmin = (long)floor(t0 < t1 ? t0 : t1) - 1;
        long kmax = (long)ceil(t0 < t1 ? t1 : t0) + 1;
        if (kmax < kmin || (kmax - kmin) > TRIGREGION_MAX_MEMBERS) { ok = false; break; }
        for (long kk = kmin; kk <= kmax; kk++) {
            Expr* m = tr_subst(value, param, expr_new_integer((int64_t)kk));
            double dm;
            if (tr_to_double(m, &dm) && dm > lo && dm < hi) cand_add(out, m);
            else expr_free(m);
        }
    }
    expr_free(sols);
    return ok;
}

/* Statement true at the exact point `pt`? */
static bool tr_stmt_true_at(const Expr* stmt, const Expr* var, const Expr* pt) {
    return tr_verify(stmt, var, pt);
}
/* Statement true at the machine point `x`? */
static bool tr_stmt_true_at_num(const Expr* stmt, const Expr* var, double x) {
    Expr* px = expr_new_real(x);
    bool ok = tr_verify(stmt, var, px);
    expr_free(px);
    return ok;
}

typedef struct { Expr* lo; bool lo_cl; Expr* hi; bool hi_cl; double key; } OutIvl;

Expr* reduce_trig_ineq_region(const Expr* conj, const Expr* var,
                              const ReduceOpts* opts) {
    if (!head_is(conj, SYM_And) || conj->data.function.arg_count < 2) return NULL;
    size_t na = conj->data.function.arg_count;

    /* Partition atoms into trig (over var) and non-trig; require at least one
     * trig atom, and reject any pole-bearing head over var. */
    const Expr** trig = (const Expr**)malloc(na * sizeof(Expr*));
    Expr** nontrig = (Expr**)malloc(na * sizeof(Expr*));
    size_t ntrig = 0, nnt = 0;
    bool bad = false;
    for (size_t i = 0; i < na; i++) {
        const Expr* atom = conj->data.function.args[i];
        if (tr_has_pole_bearing(atom, var)) { bad = true; break; }
        if (is_rel_atom(atom) && tr_has_trig_over_var(atom, var))
            trig[ntrig++] = atom;
        else
            nontrig[nnt++] = expr_copy((Expr*)atom);
    }
    if (bad || ntrig == 0 || nnt == 0) {
        for (size_t i = 0; i < nnt; i++) expr_free(nontrig[i]);
        free(trig); free(nontrig);
        return NULL;
    }

    /* Bounded feasible intervals from the non-trig atoms. */
    Expr* ntconj = (nnt == 1) ? nontrig[0]
                 : expr_new_function(expr_new_symbol(SYM_And),
                                     nontrig, nnt);   /* adopts nontrig[] */
    free(nontrig);
    Expr* region = eval_and_free(expr_new_function(expr_new_symbol(SYM_Reduce),
        (Expr*[]){ ntconj, expr_copy((Expr*)var), expr_new_symbol(SYM_Reals) }, 3));
    EIvlVec box = {0};
    bool ok = tr_region_exact(region, var, &box);
    expr_free(region);
    if (!ok) {
        for (size_t i = 0; i < box.n; i++) { expr_free(box.a[i].lo); expr_free(box.a[i].hi); }
        free(box.a); free(trig);
        return NULL;
    }

    OutIvl* outs = NULL; size_t nouts = 0, capouts = 0;
    bool declined = false;

    for (size_t bi = 0; bi < box.n && !declined; bi++) {
        double lo, hi;
        if (!tr_to_double(box.a[bi].lo, &lo) || !tr_to_double(box.a[bi].hi, &hi)
            || !(lo < hi)) { declined = true; break; }

        /* Interior critical points: zeros of every trig atom in (lo,hi). */
        CandVec zeros = {0};
        for (size_t t = 0; t < ntrig && !declined; t++) {
            Expr* resid = rel_residual(trig[t]);
            if (!tr_collect_zeros(resid, var, lo, hi, opts, &zeros)) declined = true;
            expr_free(resid);
        }
        if (declined) { for (size_t i=0;i<zeros.n;i++) expr_free(zeros.a[i].e); free(zeros.a); break; }
        qsort(zeros.a, zeros.n, sizeof(Cand), cand_cmp);

        /* Nodes: lo, sorted distinct interior zeros, hi. */
        size_t M = zeros.n + 2;
        Expr** pts = (Expr**)malloc(M * sizeof(Expr*));
        double* keys = (double*)malloc(M * sizeof(double));
        size_t m = 0;
        pts[m] = expr_copy(box.a[bi].lo); keys[m] = lo; m++;
        for (size_t i = 0; i < zeros.n; i++) {
            if (m > 1 && zeros.a[i].key - keys[m-1] < 1e-9) { expr_free(zeros.a[i].e); continue; }
            pts[m] = zeros.a[i].e; keys[m] = zeros.a[i].key; m++;   /* adopt */
        }
        pts[m] = expr_copy(box.a[bi].hi); keys[m] = hi; m++;
        free(zeros.a);   /* Expr*s adopted into pts[] or freed above */

        /* Walk cells point(0),open(0,1),point(1),…,point(m-1); merge true runs. */
        bool in_run = false;
        Expr* rlo = NULL; bool rlo_cl = false; double rlo_key = 0;
        bool last_pt = false; size_t last_i = 0;
        size_t ncells = 2 * m - 1;
        for (size_t c = 0; c < ncells; c++) {
            bool is_pt = (c % 2 == 0);
            size_t i = c / 2;
            bool t;
            if (is_pt) t = tr_stmt_true_at(conj, var, pts[i]);
            else       t = tr_stmt_true_at_num(conj, var, (keys[i] + keys[i+1]) / 2.0);
            if (t) {
                if (!in_run) {
                    in_run = true;
                    rlo = is_pt ? pts[i] : pts[i];
                    rlo_cl = is_pt;
                    rlo_key = keys[i];
                }
                last_pt = is_pt; last_i = i;
            } else if (in_run) {
                Expr* rhi = last_pt ? pts[last_i] : pts[last_i + 1];
                bool rhi_cl = last_pt;
                if (nouts == capouts) {
                    capouts = capouts ? capouts * 2 : 8;
                    outs = (OutIvl*)realloc(outs, capouts * sizeof(OutIvl));
                }
                outs[nouts].lo = expr_copy(rlo); outs[nouts].lo_cl = rlo_cl;
                outs[nouts].hi = expr_copy(rhi); outs[nouts].hi_cl = rhi_cl;
                outs[nouts].key = rlo_key; nouts++;
                in_run = false;
            }
        }
        if (in_run) {
            Expr* rhi = last_pt ? pts[last_i] : pts[last_i + 1];
            bool rhi_cl = last_pt;
            if (nouts == capouts) {
                capouts = capouts ? capouts * 2 : 8;
                outs = (OutIvl*)realloc(outs, capouts * sizeof(OutIvl));
            }
            outs[nouts].lo = expr_copy(rlo); outs[nouts].lo_cl = rlo_cl;
            outs[nouts].hi = expr_copy(rhi); outs[nouts].hi_cl = rhi_cl;
            outs[nouts].key = rlo_key; nouts++;
        }
        for (size_t i = 0; i < m; i++) expr_free(pts[i]);
        free(pts); free(keys);
    }

    for (size_t i = 0; i < box.n; i++) { expr_free(box.a[i].lo); expr_free(box.a[i].hi); }
    free(box.a); free(trig);

    if (declined) {
        for (size_t i = 0; i < nouts; i++) { expr_free(outs[i].lo); expr_free(outs[i].hi); }
        free(outs);
        return NULL;
    }
    if (nouts == 0) { free(outs); return expr_new_symbol(SYM_False); }

    /* Sort by left endpoint and render each piece. */
    for (size_t i = 1; i < nouts; i++) {          /* small n: insertion sort */
        OutIvl v = outs[i]; size_t j = i;
        while (j > 0 && outs[j-1].key > v.key) { outs[j] = outs[j-1]; j--; }
        outs[j] = v;
    }
    Expr** terms = (Expr**)malloc(nouts * sizeof(Expr*));
    for (size_t i = 0; i < nouts; i++) {
        double klo, khi;
        bool point = tr_to_double(outs[i].lo, &klo) && tr_to_double(outs[i].hi, &khi)
                     && outs[i].lo_cl && outs[i].hi_cl && (khi - klo) < 1e-9;
        if (point) {
            terms[i] = mk_fn2(SYM_Equal, expr_copy((Expr*)var), outs[i].lo);
            expr_free(outs[i].hi);
        } else {
            terms[i] = expr_new_function(expr_new_symbol(SYM_Inequality),
                (Expr*[]){ outs[i].lo,
                           expr_new_symbol(outs[i].lo_cl ? SYM_LessEqual : SYM_Less),
                           expr_copy((Expr*)var),
                           expr_new_symbol(outs[i].hi_cl ? SYM_LessEqual : SYM_Less),
                           outs[i].hi }, 5);
        }
    }
    free(outs);
    Expr* result = (nouts == 1) ? terms[0]
                 : expr_new_function(expr_new_symbol(SYM_Or), terms, nouts);
    free(terms);
    return eval_and_free(result);
}
