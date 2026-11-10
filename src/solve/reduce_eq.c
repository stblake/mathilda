/*
 * reduce_eq.c
 *
 * Complete univariate equation engine for `Reduce` (docs/notes/REDUCE_PLAN.md, Phase 1).
 * See reduce_eq.h for the algorithm.  Everything numeric is done by reusing the
 * existing Expr-level builtins (Exponent, Coefficient, Expand, Solve) rather
 * than reimplementing polynomial arithmetic.
 */
#include "reduce_eq.h"

#include "eval.h"
#include "sym_names.h"
#include "zero_test.h"

#include <stdlib.h>

/* ------------------------------------------------------------------ *
 *  Expr-level polynomial helpers (via evaluate)                       *
 * ------------------------------------------------------------------ */

static bool is_head(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == name;
}

/* Degree of `poly` in `var` via Exponent[poly, var]; -1 if not a plain int. */
static int poly_degree(const Expr* poly, const Expr* var) {
    Expr* e = eval_and_free(expr_new_function(expr_new_symbol(SYM_Exponent),
        (Expr*[]){ expr_copy((Expr*)poly), expr_copy((Expr*)var) }, 2));
    int d = (e->type == EXPR_INTEGER) ? (int)e->data.integer : -1;
    expr_free(e);
    return d;
}

/* Coefficient[poly, var, deg]. */
static Expr* poly_coeff(const Expr* poly, const Expr* var, int deg) {
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Coefficient),
        (Expr*[]){ expr_copy((Expr*)poly), expr_copy((Expr*)var),
                   expr_new_integer(deg) }, 3));
}

/* Expand[poly - lc*var^deg] -- `poly` with its leading term removed. */
static Expr* poly_drop_leading(const Expr* poly, const Expr* var,
                               const Expr* lc, int deg) {
    Expr* pw = expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ expr_copy((Expr*)var), expr_new_integer(deg) }, 2);
    Expr* term = expr_new_function(expr_new_symbol(SYM_Times),
        (Expr*[]){ expr_new_integer(-1), expr_copy((Expr*)lc), pw }, 3);
    Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus),
        (Expr*[]){ expr_copy((Expr*)poly), term }, 2);
    return eval_and_free(expr_new_function(expr_new_symbol(SYM_Expand),
        (Expr*[]){ sum }, 1));
}

static bool is_zero(const Expr* e) {
    return e && e->type == EXPR_INTEGER && e->data.integer == 0;
}

/* A nonzero numeric literal: the leading coefficient is definitely nonzero, so
 * no `lc != 0 / lc == 0` case split is needed. */
static bool is_nonzero_numeric(const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: return e->data.integer != 0;
        case EXPR_REAL:    return e->data.real != 0.0;
        case EXPR_BIGINT:  return mpz_sgn(e->data.bigint) != 0;
#ifdef USE_MPFR
        case EXPR_MPFR:    return mpfr_sgn(e->data.mpfr) != 0;
#endif
        case EXPR_FUNCTION:
            /* A Rational is always nonzero (a zero would be Integer 0). */
            return is_head(e, SYM_Rational);
        default: return false;
    }
}

/* Collect the distinct parameter symbols of `e` -- leaves that are NOT one of the
 * reduce vars (a condition is a relation among the free parameters).  Owned. */
static void collect_params(const Expr* e, Expr** vars, int nv, Expr*** out, int* n) {
    if (!e) return;
    if (e->type == EXPR_SYMBOL) {
        for (int i = 0; i < nv; i++) if (expr_eq(e, vars[i])) return;    /* a reduce var */
        for (int i = 0; i < *n; i++) if (expr_eq((*out)[i], e)) return;  /* dup */
        *out = realloc(*out, (size_t)(*n + 1) * sizeof(Expr*));
        (*out)[(*n)++] = expr_copy((Expr*)e);
        return;
    }
    if (e->type == EXPR_FUNCTION)
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            collect_params(e->data.function.args[i], vars, nv, out, n);
}

/* When `poly` is linear in a single parameter `a` with a nonzero numeric
 * coefficient, return the solved display  a REL (-const/coeff)  (owned); else
 * NULL.  Mirrors reduce_sys.c's sys_norm_condition so a leading-coefficient
 * condition prints minimally: `-1 + a != 0` becomes `a != 1`, `-1 + 2 a == 0`
 * becomes `a == 1/2`.  A nonlinear or multi-parameter condition is left verbatim. */
static Expr* cond_solved_display(const char* rel_head, const Expr* poly,
                                 Expr** vars, int nv) {
    Expr** params = NULL; int nparams = 0;
    collect_params(poly, vars, nv, &params, &nparams);
    Expr* result = NULL;
    if (nparams == 1 && poly_degree(poly, params[0]) == 1) {
        Expr* a = params[0];
        Expr* coeff = poly_coeff(poly, a, 1);
        if (is_nonzero_numeric(coeff)) {
            Expr* c0 = poly_coeff(poly, a, 0);              /* constant term */
            Expr* val = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                (Expr*[]){ expr_new_integer(-1), c0,
                           expr_new_function(expr_new_symbol(SYM_Power),
                               (Expr*[]){ expr_copy(coeff), expr_new_integer(-1) }, 2) }, 3));
            result = expr_new_function(expr_new_symbol(rel_head),
                (Expr*[]){ expr_copy(a), val }, 2);
        }
        expr_free(coeff);
    }
    for (int i = 0; i < nparams; i++) expr_free(params[i]);
    free(params);
    return result;
}

/* Build a bare relational atom  poly REL 0  (rel is Equal or Unequal).  When the
 * condition is single-parameter-linear, attach a minimal solved-form display. */
static RForm* condition_form(const char* rel_head, const Expr* poly,
                             Expr** vars, int nv, bool* ok) {
    Expr* rel = expr_new_function(expr_new_symbol(rel_head),
        (Expr*[]){ expr_copy((Expr*)poly), expr_new_integer(0) }, 2);
    RAtom a = reduce_atom_canonicalize(rel, vars, nv, ok);
    expr_free(rel);
    if (*ok) {
        Expr* disp = cond_solved_display(rel_head, poly, vars, nv);
        if (disp) { if (a.display) expr_free(a.display); a.display = disp; }
    }
    return rform_single_atom(a);
}

/* ------------------------------------------------------------------ *
 *  Generic solutions via Solve                                        *
 * ------------------------------------------------------------------ */

/* Solve[poly == 0, var] over Complexes, reformatted as an OR of solved atoms
 * (x == r_k).  Sets *ok=false on any unexpected shape (declined, conditional,
 * multi-rule, wrong variable). */
static RForm* solve_generic(const Expr* poly, const Expr* var,
                            Expr** vars, int nv, bool* ok,
                            const ReduceOpts* opts) {
    Expr* eqn = expr_new_function(expr_new_symbol(SYM_Equal),
        (Expr*[]){ expr_copy((Expr*)poly), expr_new_integer(0) }, 2);
    Expr* base[2] = { eqn, expr_copy((Expr*)var) };
    Expr* call = reduce_opts_build_solve(base, 2, opts);
    Expr* sols = eval_and_free(call);

    if (!is_head(sols, SYM_List)) { *ok = false; expr_free(sols); return rform_false(); }

    RForm* out = rform_false();
    size_t nsol = sols->data.function.arg_count;
    for (size_t i = 0; i < nsol; i++) {
        Expr* row = sols->data.function.args[i];
        if (!is_head(row, SYM_List)) { *ok = false; break; }
        if (row->data.function.arg_count == 0) {   /* {{}} -> all values */
            rform_free(out);
            out = rform_true();
            break;
        }
        if (row->data.function.arg_count != 1) { *ok = false; break; }
        Expr* rule = row->data.function.args[0];
        if (!(is_head(rule, SYM_Rule) || is_head(rule, SYM_RuleDelayed))
            || rule->data.function.arg_count != 2) { *ok = false; break; }
        Expr* lhs = rule->data.function.args[0];
        Expr* val = rule->data.function.args[1];
        if (!expr_eq(lhs, var)) { *ok = false; break; }
        if (is_head(val, SYM_ConditionalExpression)) { *ok = false; break; }
        out = rform_or(out, rform_single_atom(ratom_solved(var, val, vars, nv)));
    }
    expr_free(sols);
    if (!*ok) { rform_free(out); return rform_false(); }
    return out;
}

/* ------------------------------------------------------------------ *
 *  The leading-coefficient-vanishing recursion                        *
 * ------------------------------------------------------------------ */

static RForm* solve_case(const Expr* poly, const Expr* var,
                         Expr** vars, int nv, bool* ok,
                         const ReduceOpts* opts) {
    if (!*ok) return rform_false();

    /* An identically-zero polynomial is `0 == 0`, true for every value of the
     * remaining variable.  This is the base case of a leading-coefficient split
     * whose lower-order part vanished: `Reduce[a x == 0, x]` drops the `a x`
     * term in its `a == 0` branch, leaving 0 -> True, so the branch is `a == 0`
     * (not a decline).  Guard it here because Exponent[0, x] is -Infinity, which
     * poly_degree cannot classify. */
    if (is_zero(poly)) return rform_true();

    int deg = poly_degree(poly, var);
    if (deg < 0) { *ok = false; return rform_false(); }

    if (deg == 0) {
        /* No x left: the equation is the parameter condition poly == 0. */
        if (is_zero(poly)) return rform_true();
        return condition_form(SYM_Equal, poly, vars, nv, ok);
    }

    Expr* lc = poly_coeff(poly, var, deg);
    RForm* gen = solve_generic(poly, var, vars, nv, ok, opts);
    if (!*ok) { expr_free(lc); rform_free(gen); return rform_false(); }

    if (is_nonzero_numeric(lc)) {          /* terminal: no case split needed */
        expr_free(lc);
        return gen;
    }

    /* Branch A: lc != 0  &&  (generic solutions of degree deg). */
    RForm* A = rform_and(condition_form(SYM_Unequal, lc, vars, nv, ok), gen);

    /* Branch B: lc == 0  &&  solve(poly with leading term dropped). */
    Expr* reduced = poly_drop_leading(poly, var, lc, deg);
    RForm* B = rform_and(condition_form(SYM_Equal, lc, vars, nv, ok),
                         solve_case(reduced, var, vars, nv, ok, opts));
    expr_free(reduced);
    expr_free(lc);

    return rform_or(A, B);
}

RForm* reduce_eq_univariate(const Expr* poly, const Expr* var,
                            Expr** vars, int nv, bool* ok,
                            const ReduceOpts* opts) {
    *ok = true;
    return solve_case(poly, var, vars, nv, ok, opts);
}

/* ------------------------------------------------------------------ *
 *  Transcendental single equation -> logical formula                  *
 * ------------------------------------------------------------------ */

/* ------------------------------------------------------------------ *
 *  Reality of a solution value over Reals                             *
 * ------------------------------------------------------------------ */

static Expr* rc_fn2(const char* h, Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a, b }, 2);
}
static Expr* rc_fn3(const char* h, Expr* a, Expr* b, Expr* c) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a, b, c }, 3);
}
/* And[a, b] (both adopted; either may be NULL meaning "no condition"). */
static Expr* rc_and(Expr* a, Expr* b) {
    if (!a) return b;
    if (!b) return a;
    return rc_fn2(SYM_And, a, b);
}

/* True iff `e` has a leaf symbol that is not a numeric constant (Pi, E, ...):
 * a free parameter such as `a`.  Head symbols (ArcSin, Log) are not leaves. */
static bool has_free_symbol(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) {
        Expr* q = eval_and_free(expr_new_function(expr_new_symbol(SYM_NumericQ),
            (Expr*[]){ expr_copy((Expr*)e) }, 1));
        bool free_sym = !(q->type == EXPR_SYMBOL && q->data.symbol.name == SYM_True);
        expr_free(q);
        return free_sym;
    }
    if (e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_FUNCTION
        && has_free_symbol(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_free_symbol(e->data.function.args[i])) return true;
    return false;
}

/* Unevaluated EXACT condition (owned) on the real parameters of `e` under which
 * the value `e` is real, or NULL when no exact condition is known.  Over Reals
 * every parameter is real (Mathematica's Reduce convention), so the condition
 * is the conjunction of the domain restrictions of the real-valued branches of
 * the heads in `e`: Log[u] needs u > 0, ArcSin/ArcCos[u] need -1 <= u <= 1, a
 * fractional power needs a non-negative base, ArcCosh[u] needs u >= 1, ...;
 * heads that are real on the whole real line (Sin, Exp, ArcTan, ArcSinh, ...)
 * contribute only their argument's condition.  A complex literal or an unknown
 * head yields NULL: the caller then falls back to a provable decision or
 * declines -- it never guesses. */
static Expr* real_value_cond_rec(const Expr* e) {
    if (!e) return NULL;
    switch (e->type) {
        case EXPR_INTEGER: case EXPR_REAL: case EXPR_BIGINT:
            return expr_new_symbol(SYM_True);
        case EXPR_SYMBOL: {
            const char* n = e->data.symbol.name;
            if (n == SYM_I || n == SYM_ComplexInfinity || n == SYM_Indeterminate)
                return NULL;
            return expr_new_symbol(SYM_True);         /* real parameter / constant */
        }
        case EXPR_FUNCTION: break;
        default: return NULL;
    }
    const Expr* hd = e->data.function.head;
    if (hd->type != EXPR_SYMBOL) return NULL;
    const char* h = hd->data.symbol.name;
    size_t n = e->data.function.arg_count;
    if (h == SYM_Rational) return expr_new_symbol(SYM_True);
    if (h == SYM_Complex) return NULL;

    if (h == SYM_Plus || h == SYM_Times) {
        Expr* acc = NULL;
        for (size_t i = 0; i < n; i++) {
            Expr* c = real_value_cond_rec(e->data.function.args[i]);
            if (!c) { expr_free(acc); return NULL; }
            acc = rc_and(acc, c);
        }
        return acc ? acc : expr_new_symbol(SYM_True);
    }

    if (h == SYM_Power && n == 2) {
        const Expr* b = e->data.function.args[0];
        const Expr* x = e->data.function.args[1];
        if (x->type == EXPR_INTEGER) {
            Expr* cb = real_value_cond_rec(b);
            if (!cb) return NULL;
            if (x->data.integer < 0)
                cb = rc_and(cb, rc_fn2(SYM_Unequal, expr_copy((Expr*)b),
                                       expr_new_integer(0)));
            return cb;
        }
        if (is_head(x, SYM_Rational) && x->data.function.arg_count == 2) {
            /* Principal fractional power: real iff the base is >= 0 (> 0 for a
             * negative exponent); (-8)^(1/3) is NOT real. */
            Expr* cb = real_value_cond_rec(b);
            if (!cb) return NULL;
            const Expr* num = x->data.function.args[0];
            bool neg = (num->type == EXPR_INTEGER && num->data.integer < 0);
            return rc_and(cb, rc_fn2(neg ? SYM_Greater : SYM_GreaterEqual,
                                     expr_copy((Expr*)b), expr_new_integer(0)));
        }
        /* b^x with a positive constant base (E, 2, Pi, ...): real iff x is. */
        if (!has_free_symbol(b)) {
            Expr* pos = eval_and_free(rc_fn2(SYM_Greater, expr_copy((Expr*)b),
                                             expr_new_integer(0)));
            bool bpos = (pos->type == EXPR_SYMBOL && pos->data.symbol.name == SYM_True);
            expr_free(pos);
            if (bpos) return real_value_cond_rec(x);
        }
        return NULL;
    }

    if (n == 2 && h == SYM_ArcTan) {
        Expr* c0 = real_value_cond_rec(e->data.function.args[0]);
        if (!c0) return NULL;
        Expr* c1 = real_value_cond_rec(e->data.function.args[1]);
        if (!c1) { expr_free(c0); return NULL; }
        return rc_and(c0, c1);
    }
    if (n == 2 && h == SYM_ProductLog
        && e->data.function.args[0]->type == EXPR_INTEGER) {
        /* Branch k of W at u: W_0 is real iff u >= -1/E; W_{-1} iff
         * -1/E <= u < 0; every other branch is never real. */
        int64_t k = e->data.function.args[0]->data.integer;
        const Expr* u = e->data.function.args[1];
        if (k != 0 && k != -1) return expr_new_symbol(SYM_False);
        Expr* cu = real_value_cond_rec(u);
        if (!cu) return NULL;
        Expr* lo = rc_fn2(SYM_GreaterEqual, expr_copy((Expr*)u),
                          rc_fn2(SYM_Times, expr_new_integer(-1),
                                 rc_fn2(SYM_Power, expr_new_symbol(SYM_E), expr_new_integer(-1))));
        Expr* dom = k == 0 ? lo
                  : rc_fn2(SYM_And, lo, rc_fn2(SYM_Less, expr_copy((Expr*)u), expr_new_integer(0)));
        return rc_and(cu, dom);
    }
    if (n != 1) return NULL;

    const Expr* u = e->data.function.args[0];
    Expr* cu = real_value_cond_rec(u);
    if (!cu) return NULL;
    Expr* uc = expr_copy((Expr*)u);
    Expr* dom = NULL;                                  /* NULL: whole real line */
    if (h == SYM_Sin || h == SYM_Cos || h == SYM_Tan || h == SYM_Cot
        || h == SYM_Sec || h == SYM_Csc || h == SYM_Exp || h == SYM_Sinh
        || h == SYM_Cosh || h == SYM_Tanh || h == SYM_ArcTan || h == SYM_ArcCot
        || h == SYM_ArcSinh || h == SYM_Abs) {
        expr_free(uc);
    } else if (h == SYM_Log) {
        dom = rc_fn2(SYM_Greater, uc, expr_new_integer(0));
    } else if (h == SYM_Sqrt) {
        dom = rc_fn2(SYM_GreaterEqual, uc, expr_new_integer(0));
    } else if (h == SYM_ArcSin || h == SYM_ArcCos) {
        dom = rc_fn3(SYM_LessEqual, expr_new_integer(-1), uc, expr_new_integer(1));
    } else if (h == SYM_ArcCosh) {
        dom = rc_fn2(SYM_GreaterEqual, uc, expr_new_integer(1));
    } else if (h == SYM_ArcTanh) {
        dom = rc_fn3(SYM_Less, expr_new_integer(-1), uc, expr_new_integer(1));
    } else if (h == SYM_ArcSec || h == SYM_ArcCsc) {
        dom = rc_fn2(SYM_Or,
                     rc_fn2(SYM_LessEqual, expr_copy(uc), expr_new_integer(-1)),
                     rc_fn2(SYM_GreaterEqual, uc, expr_new_integer(1)));
    } else if (h == SYM_ArcCoth) {
        dom = rc_fn2(SYM_Or,
                     rc_fn2(SYM_Less, expr_copy(uc), expr_new_integer(-1)),
                     rc_fn2(SYM_Greater, uc, expr_new_integer(1)));
    } else if (h == SYM_ArcSech) {
        dom = rc_fn2(SYM_And,
                     rc_fn2(SYM_Greater, expr_copy(uc), expr_new_integer(0)),
                     rc_fn2(SYM_LessEqual, uc, expr_new_integer(1)));
    } else if (h == SYM_ArcCsch) {
        dom = rc_fn2(SYM_Unequal, uc, expr_new_integer(0));
    } else if (h == SYM_ProductLog) {
        /* Principal branch W0: real iff u >= -1/E. */
        dom = rc_fn2(SYM_GreaterEqual, uc,
                     rc_fn2(SYM_Times, expr_new_integer(-1),
                            rc_fn2(SYM_Power, expr_new_symbol(SYM_E),
                                   expr_new_integer(-1))));
    } else {
        expr_free(uc); expr_free(cu);
        return NULL;
    }
    return rc_and(cu, dom);
}

/* Realness of a solution value over Reals: 1 real, 0 never real, 2 real
 * exactly when the owned, evaluated *cond_out holds (a condition on the real
 * parameters), -1 undecidable (the caller must decline).
 *
 * A value carrying free parameters is decided by its exact reality condition
 * (real_value_cond_rec), NOT by the numeric zero test: Im[ArcSin[a]] is
 * nonzero for a generically sampled `a`, and trusting that verdict is
 * precisely what made Reduce[Sin[x] == a, x, Reals] drop both families and
 * return False.  A parameter-free value is decided by the zero test, with the
 * reality condition (which then evaluates to True/False) tried first. */
int reduce_value_reality(const Expr* v, Expr** cond_out) {
    *cond_out = NULL;
    bool params = has_free_symbol(v);
    Expr* rc = real_value_cond_rec(v);
    if (rc) {
        rc = eval_and_free(rc);
        if (rc->type == EXPR_SYMBOL && rc->data.symbol.name == SYM_True) {
            expr_free(rc); return 1;
        }
        if (rc->type == EXPR_SYMBOL && rc->data.symbol.name == SYM_False) {
            expr_free(rc); return 0;
        }
        if (params) {
            /* Restate the condition in its solved form over the parameters
             * (ArcSin[-1 - a] gives -1 <= -1 - a <= 1, i.e. -2 <= a <= 0) when
             * Reduce decides it; otherwise keep it as built. */
            Expr* red = eval_and_free(expr_new_function(expr_new_symbol(SYM_Reduce),
                (Expr*[]){ expr_copy(rc), expr_new_symbol(SYM_Reals) }, 2));
            if (is_head(red, SYM_Reduce)) {
                expr_free(red);
            } else {
                expr_free(rc);
                rc = red;
                if (rc->type == EXPR_SYMBOL && rc->data.symbol.name == SYM_True) {
                    expr_free(rc); return 1;
                }
                if (rc->type == EXPR_SYMBOL && rc->data.symbol.name == SYM_False) {
                    expr_free(rc); return 0;
                }
            }
            *cond_out = rc;
            return 2;
        }
        expr_free(rc);                 /* numeric but unevaluated: zero-test below */
    } else if (params && is_head(v, SYM_Plus)) {
        /* c + r with a parameter-free non-real constant c and an
         * UNCONDITIONALLY real rest r (Sinh's offset I Pi - ArcSinh[a]) is never
         * real.  A rest that is only conditionally real is not decided here: a
         * complex rest could cancel Im[c]. */
        size_t n = v->data.function.arg_count;
        Expr** cs = (Expr**)malloc(n * sizeof(Expr*));
        Expr** rs = (Expr**)malloc(n * sizeof(Expr*));
        size_t nc = 0, nr = 0;
        for (size_t i = 0; i < n; i++) {
            const Expr* t = v->data.function.args[i];
            if (has_free_symbol(t)) rs[nr++] = expr_copy((Expr*)t);
            else                    cs[nc++] = expr_copy((Expr*)t);
        }
        Expr* csum = expr_new_function(expr_new_symbol(SYM_Plus), cs, nc);
        Expr* rsum = expr_new_function(expr_new_symbol(SYM_Plus), rs, nr);
        free(cs); free(rs);
        Expr* rrc = real_value_cond_rec(rsum);
        bool rest_real = false;
        if (rrc) {
            rrc = eval_and_free(rrc);
            rest_real = (rrc->type == EXPR_SYMBOL && rrc->data.symbol.name == SYM_True);
            expr_free(rrc);
        }
        int verdict = -1;
        if (rest_real && nc > 0) {
            Expr* imc = eval_and_free(expr_new_function(expr_new_symbol("Im"),
                (Expr*[]){ csum }, 1));
            csum = NULL;
            if (zero_test_decide(imc) == ZERO_TEST_FALSE) verdict = 0;
            expr_free(imc);
        }
        expr_free(csum); expr_free(rsum);
        if (verdict == 0) return 0;
    }
    Expr* im = eval_and_free(expr_new_function(expr_new_symbol("Im"),
        (Expr*[]){ expr_copy((Expr*)v) }, 1));
    ZeroTestResult z = zero_test_decide(im);
    expr_free(im);
    if (z == ZERO_TEST_TRUE) return 1;
    if (z == ZERO_TEST_FALSE && !params) return 0;
    return -1;
}

/* A Solve side condition restated for REAL parameters (owned, evaluated; NULL
 * when it reduces to True): Im[u] -> 0 and Re[u] -> u wherever u is provably
 * real.  Log's branch strip -Pi < Im[a] <= Pi is vacuous for a real `a`. */
static Expr* realify_rec(const Expr* c) {
    if (c->type != EXPR_FUNCTION) return expr_copy((Expr*)c);
    if ((is_head(c, SYM_Im) || is_head(c, SYM_Re)) && c->data.function.arg_count == 1) {
        Expr* rc = NULL;
        if (reduce_value_reality(c->data.function.args[0], &rc) == 1)
            return is_head(c, SYM_Im) ? expr_new_integer(0)
                                      : expr_copy(c->data.function.args[0]);
        expr_free(rc);
        return expr_copy((Expr*)c);
    }
    size_t n = c->data.function.arg_count;
    Expr** a = (Expr**)malloc((n ? n : 1) * sizeof(Expr*));
    for (size_t i = 0; i < n; i++) a[i] = realify_rec(c->data.function.args[i]);
    Expr* out = expr_new_function(expr_copy(c->data.function.head), a, n);
    free(a);
    return out;
}
static Expr* realify_cond(const Expr* c) {
    if (!c) return NULL;
    Expr* r = eval_and_free(realify_rec(c));
    if (r->type == EXPR_SYMBOL && r->data.symbol.name == SYM_True) {
        expr_free(r); return NULL;
    }
    return r;
}

/* A `{lhs -> ConditionalExpression[value, cond]}` row (`value` and `cond`
 * adopted; a NULL cond gives the plain `{lhs -> value}` row). */
static Expr* cond_row(const Expr* lhs, Expr* value, Expr* cond) {
    Expr* v = cond ? rc_fn2(SYM_ConditionalExpression, value, cond) : value;
    Expr* r = rc_fn2(SYM_Rule, expr_copy((Expr*)lhs), v);
    return expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ r }, 1);
}

/* Over Reals, collapse each periodic solution family to its real members.
 * A ConditionalExpression value a + b*C[k] whose period b AND offset a are both
 * real keeps the whole family (e.g. a trig ArcSin[c] + 2 Pi C[k] with real c);
 * a real period with a NON-real offset (e.g. ArcSin[2] + 2 Pi C[k], where
 * ArcSin[2] is complex) has every member non-real and is dropped, so
 * Reduce[Sin[x]==2,x,Reals] -> False.  A non-real period (imaginary, e.g.
 * Sinh/Exp's 2 I Pi C[k]) has a single candidate real member at C[k]=0 (kept if
 * real, else dropped): Reduce[Sinh[y]==0,y,Reals] -> y==0 and
 * Reduce[Exp[y]==3,y,Reals] -> y==Log[3].
 *
 * A PARAMETRIC offset (ArcSin[a], Log[a]) is real exactly when its reality
 * condition on the real parameters holds, and that condition is conjoined onto
 * the row: Reduce[Sin[x]==a,x,Reals] -> -1 <= a <= 1 && C[1] in Integers &&
 * ..., Reduce[Exp[x]==a,x,Reals] -> a > 0 && x == Log[a].  Returns NULL
 * (decline) when the reality of some value can be neither proved nor
 * conditioned.  Consumes `sols`. */
Expr* reduce_collapse_reals_families(Expr* sols) {
    if (!is_head(sols, SYM_List)) return sols;
    size_t n = sols->data.function.arg_count;
    Expr** kept = (Expr**)malloc((n ? n : 1) * sizeof(Expr*));
    size_t nk = 0;
    bool declined = false;
    for (size_t i = 0; i < n && !declined; i++) {
        Expr* row = sols->data.function.args[i];
        if (!is_head(row, SYM_List) || row->data.function.arg_count != 1) {
            kept[nk++] = expr_copy(row);   /* {{}} (all reals) or unusual: keep */
            continue;
        }
        Expr* rule = row->data.function.args[0];
        if (!is_head(rule, SYM_Rule) || rule->data.function.arg_count != 2) {
            kept[nk++] = expr_copy(row);
            continue;
        }
        Expr* lhs = rule->data.function.args[0];
        Expr* val = rule->data.function.args[1];
        const Expr* value = val;
        const Expr* param = NULL;
        const Expr* row_cond = NULL;      /* the row's own (non-family) condition */
        if (is_head(val, SYM_ConditionalExpression)
            && val->data.function.arg_count == 2) {
            value = val->data.function.args[0];
            Expr* cond = val->data.function.args[1];
            if (is_head(cond, SYM_Element) && cond->data.function.arg_count == 2)
                param = cond->data.function.args[0];
            else
                row_cond = cond;
        }
        if (param) {
            Expr* coeff = eval_and_free(expr_new_function(
                expr_new_symbol("Coefficient"),
                (Expr*[]){ expr_copy((Expr*)value), expr_copy((Expr*)param) }, 2));
            /* The family's offset: its member at param = 0. */
            Expr* v0 = eval_and_free(expr_new_function(
                expr_new_symbol("ReplaceAll"),
                (Expr*[]){ expr_copy((Expr*)value),
                           expr_new_function(expr_new_symbol(SYM_Rule),
                               (Expr*[]){ expr_copy((Expr*)param),
                                          expr_new_integer(0) }, 2) }, 2));
            Expr* pc = NULL;
            int period = reduce_value_reality(coeff, &pc);
            expr_free(coeff);
            if (period == 2 || period == -1) {   /* parametric / unknown period */
                expr_free(pc); expr_free(v0);
                declined = true; break;
            }
            Expr* oc = NULL;
            int off = reduce_value_reality(v0, &oc);
            if (off == -1) { expr_free(v0); declined = true; break; }
            if (period == 1) {
                /* Real period: the whole family is real iff the offset is (a
                 * real period cannot cancel a complex offset). */
                expr_free(v0);
                if (off == 1)
                    kept[nk++] = expr_copy(row);
                else if (off == 2)
                    kept[nk++] = cond_row(lhs, expr_copy((Expr*)value),
                        rc_and(oc, expr_copy(val->data.function.args[1])));
            } else {                                   /* non-real period */
                if (off == 1)      kept[nk++] = cond_row(lhs, v0, NULL);
                else if (off == 2) kept[nk++] = cond_row(lhs, v0, oc);
                else               expr_free(v0);
            }
        } else {
            Expr* vc = NULL;
            int r = reduce_value_reality(value, &vc);
            if (r == -1) { declined = true; break; }
            Expr* extra = realify_cond(row_cond);
            if (extra && extra->type == EXPR_SYMBOL
                && extra->data.symbol.name == SYM_False) {
                expr_free(extra); expr_free(vc);    /* side condition fails */
                continue;
            }
            if (r == 1)
                kept[nk++] = cond_row(lhs, expr_copy((Expr*)value), extra);
            else if (r == 2)
                kept[nk++] = cond_row(lhs, expr_copy((Expr*)value), rc_and(vc, extra));
            else
                expr_free(extra);              /* no real member -> drop */
        }
    }
    Expr* out = NULL;
    if (declined) {
        for (size_t i = 0; i < nk; i++) expr_free(kept[i]);
    } else {
        out = expr_new_function(expr_new_symbol(SYM_List), kept, nk);
    }
    free(kept);
    expr_free(sols);
    return out;
}

Expr* reduce_eq_transcendental(const Expr* poly, const Expr* var,
                               const Expr* dom, const ReduceOpts* opts) {
    /* Re-enter Solve[poly == 0, var (, dom)]; base args are adopted by the
     * Solve call.  An explicit `dom` (Reals/Complexes) is forwarded so the
     * inverse-function solver produces domain-appropriate solutions and
     * conditions; a NULL dom leaves Solve at its Complexes default. */
    Expr* eqn = expr_new_function(expr_new_symbol(SYM_Equal),
        (Expr*[]){ expr_copy((Expr*)poly), expr_new_integer(0) }, 2);
    Expr* base[3];
    int nb = 0;
    base[nb++] = eqn;
    base[nb++] = expr_copy((Expr*)var);
    if (dom) base[nb++] = expr_copy((Expr*)dom);
    Expr* sols = eval_and_free(reduce_opts_build_solve(base, nb, opts));

    /* Solve declined (unevaluated / non-list): let the caller fall back. */
    if (!is_head(sols, SYM_List)) { expr_free(sols); return NULL; }

    /* Over Reals, keep only the real members of each periodic family. */
    if (dom && dom->type == EXPR_SYMBOL && dom->data.symbol.name == SYM_Reals) {
        sols = reduce_collapse_reals_families(sols);
        if (!sols) return NULL;          /* a value of undecidable reality */
    }

    size_t nsol = sols->data.function.arg_count;
    Expr** terms = (Expr**)malloc((nsol ? nsol : 1) * sizeof(Expr*));
    size_t nt = 0;
    bool ok = true;
    Expr* all_values = NULL;   /* set to True when a {{}} row means "any var" */

    for (size_t i = 0; i < nsol; i++) {
        Expr* row = sols->data.function.args[i];
        if (!is_head(row, SYM_List)) { ok = false; break; }
        if (row->data.function.arg_count == 0) {           /* {{}} -> True */
            all_values = expr_new_symbol(SYM_True);
            break;
        }
        if (row->data.function.arg_count != 1) { ok = false; break; }
        Expr* rule = row->data.function.args[0];
        if (!(is_head(rule, SYM_Rule) || is_head(rule, SYM_RuleDelayed))
            || rule->data.function.arg_count != 2) { ok = false; break; }
        Expr* lhs = rule->data.function.args[0];
        Expr* val = rule->data.function.args[1];
        if (!expr_eq(lhs, var)) { ok = false; break; }

        Expr* term;
        if (is_head(val, SYM_ConditionalExpression)
            && val->data.function.arg_count == 2) {
            /* ConditionalExpression[value, cond] -> cond && (var == value). */
            Expr* value = val->data.function.args[0];
            Expr* cond  = val->data.function.args[1];
            Expr* eq = expr_new_function(expr_new_symbol(SYM_Equal),
                (Expr*[]){ expr_copy((Expr*)var), expr_copy(value) }, 2);
            term = expr_new_function(expr_new_symbol(SYM_And),
                (Expr*[]){ expr_copy(cond), eq }, 2);
        } else {
            term = expr_new_function(expr_new_symbol(SYM_Equal),
                (Expr*[]){ expr_copy((Expr*)var), expr_copy(val) }, 2);
        }
        terms[nt++] = term;
    }
    expr_free(sols);

    if (!ok) {
        for (size_t i = 0; i < nt; i++) expr_free(terms[i]);
        free(terms);
        if (all_values) expr_free(all_values);
        return NULL;
    }
    if (all_values) {
        for (size_t i = 0; i < nt; i++) expr_free(terms[i]);
        free(terms);
        return all_values;
    }
    if (nt == 0) { free(terms); return expr_new_symbol(SYM_False); }  /* {} */

    Expr* out;
    if (nt == 1) { out = terms[0]; free(terms); }
    else { out = expr_new_function(expr_new_symbol(SYM_Or), terms, nt); free(terms); }
    /* Let the evaluator flatten / tidy the And/Or combination. */
    return eval_and_free(out);
}
