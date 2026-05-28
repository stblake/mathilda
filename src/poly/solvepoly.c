/*
 * solvepoly.c
 *
 * Single-variable polynomial-equality solver -- the first specialist
 * dispatched from `Solve` (src/solve.c).  Also reachable directly as
 * the context-qualified builtin `Solve`SolvePolynomialEquality`.
 *
 * Algorithm: plans/SOLVE_PLAN.md.  Memory contract: every helper returns a
 * freshly-owned Expr* (or Expr**) -- inputs are borrowed and deep-
 * copied wherever they appear in the output.  On NULL return from the
 * top-level dispatcher, `res` (the original builtin argument) is left
 * untouched so the evaluator preserves the unevaluated form.
 *
 * expr_new_function memcpys the args[] array (src/expr.c:88-93), so
 * passing a stack-allocated compound literal as args is safe -- the
 * pointers it holds are transferred to the new Expr.
 */

#include "solvepoly.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "arithmetic.h"
#include "attr.h"
#include "eval.h"
#include "expr.h"
#include "internal.h"
#include "poly.h"
#include "sym_names.h"
#include "symtab.h"
#include "zero_test.h"

/* ------------------------------------------------------------------ *
 *  Expression-building shorthand.                                     *
 *                                                                    *
 *  The mk_* helpers take ownership of their pointer arguments; they  *
 *  do not deep-copy.  The caller is expected to pass freshly built   *
 *  or expr_copy()'d Expr*'s.                                          *
 * ------------------------------------------------------------------ */

static Expr* mk_int(int64_t v) { return expr_new_integer(v); }
static Expr* mk_sym(const char* name) { return expr_new_symbol(name); }

static Expr* mk_fn1(const char* head, Expr* a) {
    return expr_new_function(mk_sym(head), (Expr*[]){ a }, 1);
}
static Expr* mk_fn2(const char* head, Expr* a, Expr* b) {
    return expr_new_function(mk_sym(head), (Expr*[]){ a, b }, 2);
}
static Expr* mk_fn3(const char* head, Expr* a, Expr* b, Expr* c) {
    return expr_new_function(mk_sym(head), (Expr*[]){ a, b, c }, 3);
}
static Expr* mk_fn4(const char* head, Expr* a, Expr* b, Expr* c, Expr* d) {
    return expr_new_function(mk_sym(head), (Expr*[]){ a, b, c, d }, 4);
}
static Expr* mk_fn5(const char* head, Expr* a, Expr* b, Expr* c, Expr* d,
                    Expr* e) {
    return expr_new_function(mk_sym(head), (Expr*[]){ a, b, c, d, e }, 5);
}

static Expr* mk_pow(Expr* base, Expr* exp) { return mk_fn2("Power", base, exp); }
static Expr* mk_inv(Expr* e) { return mk_pow(e, mk_int(-1)); }
static Expr* mk_neg(Expr* e) { return mk_fn2("Times", mk_int(-1), e); }
static Expr* mk_rat(int64_t p, int64_t q) {
    return mk_fn2("Rational", mk_int(p), mk_int(q));
}
static Expr* mk_sqrt(Expr* e) { return mk_pow(e, mk_rat(1, 2)); }
static Expr* mk_rule(Expr* lhs, Expr* rhs) {
    return mk_fn2("Rule", lhs, rhs);
}

/* List[a, b, ...] taking ownership of each.  `args` may be NULL when n == 0. */
static Expr* mk_list(Expr** args, size_t n) {
    return expr_new_function(mk_sym("List"), args, n);
}

/* True iff `e` is structurally the integer 0. */
static bool is_int_zero(const Expr* e) {
    return e && e->type == EXPR_INTEGER && e->data.integer == 0;
}

/* True iff `e` is a *concrete* numeric zero (Integer 0, BigInt 0,
 * Real 0.0, or Rational[0, _]).  Used by the extraneous-root filter:
 * we only drop a candidate when we can *prove* the denominator is
 * zero at it -- a symbolic residue like a/b means "undetermined" and
 * the candidate is kept. */
static bool is_definite_zero(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) return e->data.integer == 0;
    if (e->type == EXPR_BIGINT)  return mpz_sgn(e->data.bigint) == 0;
    if (e->type == EXPR_REAL)    return e->data.real == 0.0;
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Rational
        && e->data.function.arg_count == 2) {
        return is_definite_zero(e->data.function.args[0]);
    }
    return false;
}

/* Inspect a concrete-numeric expression and return its sign:
 *    1 / -1 / 0 for positive / negative / zero,
 *    INT_MIN if the sign cannot be determined symbolically.
 * Mathilda has no Sign[] builtin yet, so we walk the tree directly.
 * Handles Integer, BigInt, Real, and Rational[n, d] (canonical form
 * with d > 0).  Anything else returns INT_MIN. */
static int try_sign(Expr* e) {
    if (!e) return INT_MIN;
    if (e->type == EXPR_INTEGER) {
        if (e->data.integer > 0) return 1;
        if (e->data.integer < 0) return -1;
        return 0;
    }
    if (e->type == EXPR_BIGINT) {
        int s = mpz_sgn(e->data.bigint);
        return s == 0 ? 0 : (s > 0 ? 1 : -1);
    }
    if (e->type == EXPR_REAL) {
        if (e->data.real > 0.0) return 1;
        if (e->data.real < 0.0) return -1;
        return 0;
    }
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Rational
        && e->data.function.arg_count == 2) {
        /* Mathilda's canonical Rational has positive denominator,
         * so the sign matches the numerator. */
        return try_sign(e->data.function.args[0]);
    }
    return INT_MIN;
}

/* ------------------------------------------------------------------ *
 *  Solution-list accumulator.                                         *
 * ------------------------------------------------------------------ */

typedef struct {
    Expr** vals;
    size_t count;
    size_t capacity;
} SolList;

static void sl_init(SolList* sl) {
    sl->vals = NULL;
    sl->count = 0;
    sl->capacity = 0;
}

static void sl_push(SolList* sl, Expr* owned_val) {
    if (sl->count == sl->capacity) {
        size_t nc = sl->capacity ? sl->capacity * 2 : 8;
        sl->vals = (Expr**)realloc(sl->vals, sizeof(Expr*) * nc);
        sl->capacity = nc;
    }
    sl->vals[sl->count++] = owned_val;
}

static void sl_push_with_mult(SolList* sl, Expr* owned_val, int64_t mult) {
    if (mult <= 0) { expr_free(owned_val); return; }
    sl_push(sl, owned_val);
    for (int64_t i = 1; i < mult; i++) {
        sl_push(sl, expr_copy(sl->vals[sl->count - 1]));
    }
}

static void sl_extend_with_mult(SolList* sl, Expr** owned_arr, size_t n,
                                int64_t mult) {
    for (size_t i = 0; i < n; i++) sl_push_with_mult(sl, owned_arr[i], mult);
}

static void sl_free(SolList* sl) {
    for (size_t i = 0; i < sl->count; i++) expr_free(sl->vals[i]);
    free(sl->vals);
    sl->vals = NULL;
    sl->count = sl->capacity = 0;
}

/* ------------------------------------------------------------------ *
 *  Per-degree solvers.                                                *
 * ------------------------------------------------------------------ */

/* a*x + b == 0  →  x = -b/a */
static Expr* solve_linear(Expr* a, Expr* b) {
    return eval_and_free(mk_fn2("Times", mk_neg(expr_copy(b)),
                                          mk_inv(expr_copy(a))));
}

/* a*x^2 + b*x + c == 0.  Discriminant-aware in Reals. */
static Expr** solve_quadratic(Expr* a, Expr* b, Expr* c,
                              bool reals_only, size_t* out_n) {
    /* D = b^2 - 4 a c */
    Expr* b_sq    = eval_and_free(mk_pow(expr_copy(b), mk_int(2)));
    Expr* four_ac = eval_and_free(mk_fn3("Times", mk_int(4),
                                                  expr_copy(a),
                                                  expr_copy(c)));
    Expr* D = eval_and_free(mk_fn2("Plus", b_sq, mk_neg(four_ac)));

    int s = INT_MIN;
    if (reals_only) s = try_sign(D);

    if (reals_only && s == -1) {
        expr_free(D);
        *out_n = 0;
        return NULL;
    }

    Expr* sqrtD  = eval_and_free(mk_sqrt(D));
    Expr* inv_2a = eval_and_free(mk_inv(mk_fn2("Times", mk_int(2),
                                                        expr_copy(a))));
    Expr* neg_b  = mk_neg(expr_copy(b));

    /* r_minus = (-b - sqrtD) * inv_2a */
    Expr* num_m = mk_fn2("Plus", expr_copy(neg_b), mk_neg(expr_copy(sqrtD)));
    Expr* r_minus = eval_and_free(mk_fn2("Times", num_m, expr_copy(inv_2a)));

    /* r_plus = (-b + sqrtD) * inv_2a */
    Expr* num_p = mk_fn2("Plus", neg_b, sqrtD);
    Expr* r_plus = eval_and_free(mk_fn2("Times", num_p, inv_2a));

    /* When Δ = 0 we deliberately emit *both* r_minus and r_plus rather
     * than collapsing to a single rule: Mathematica's convention is to
     * preserve algebraic multiplicity in the solution list, so
     * `Solve[(x-1)^2 == 0, x, Reals]` returns `{{x -> 1}, {x -> 1}}`
     * just like the Complexes case.  r_minus and r_plus evaluate to
     * the same value here, which is the intended duplication. */
    Expr** out = (Expr**)malloc(sizeof(Expr*) * 2);
    out[0] = r_minus;
    out[1] = r_plus;
    *out_n = 2;
    return out;
}

/* a*x^n + b == 0  →  r = (-b/a)^(1/n);  Complexes: r * Exp[2 Pi I k / n]. */
static Expr** solve_binomial(Expr* a, Expr* b, int64_t n,
                             bool reals_only, size_t* out_n) {
    if (n <= 0) { *out_n = 0; return NULL; }

    Expr* base = eval_and_free(mk_fn2("Times", mk_neg(expr_copy(b)),
                                               mk_inv(expr_copy(a))));
    Expr* r = eval_and_free(mk_pow(expr_copy(base), mk_rat(1, n)));

    if (reals_only) {
        int s = try_sign(base);
        if (n % 2 == 1) {
            /* Odd-degree binomial: the unique real root is
             *   base^(1/n)         if base ≥ 0,
             *   0                  if base == 0,
             *   -((-base)^(1/n))   if base < 0.
             * The third branch matters because Power's principal-root
             * convention canonicalises `(-c)^(1/n)` (c > 0) to a
             * complex value `c^(1/n) * (-1)^(1/n)`, not to the real
             * `-c^(1/n)`.  For unresolved sign we keep `r` -- the
             * caller is responsible for treating it as the symbolic
             * positive branch in that case. */
            Expr* root;
            if (s == 0) {
                expr_free(r);
                root = mk_int(0);
            } else if (s == 1) {
                root = r;
            } else if (s == -1) {
                expr_free(r);
                Expr* neg_base = eval_and_free(mk_neg(expr_copy(base)));
                Expr* pos_root = eval_and_free(mk_pow(neg_base, mk_rat(1, n)));
                root = eval_and_free(mk_neg(pos_root));
            } else {
                root = r;
            }
            expr_free(base);
            Expr** out = (Expr**)malloc(sizeof(Expr*) * 1);
            out[0] = root;
            *out_n = 1;
            return out;
        }
        if (s == 0) {
            expr_free(base); expr_free(r);
            Expr** out = (Expr**)malloc(sizeof(Expr*) * 1);
            out[0] = mk_int(0);
            *out_n = 1;
            return out;
        }
        if (s == 1) {
            expr_free(base);
            Expr** out = (Expr**)malloc(sizeof(Expr*) * 2);
            out[0] = eval_and_free(mk_neg(expr_copy(r)));
            out[1] = r;
            *out_n = 2;
            return out;
        }
        if (s == -1) {
            expr_free(base); expr_free(r);
            *out_n = 0;
            return NULL;
        }
        /* unresolved sign: fall through */
    }

    expr_free(base);
    Expr** out = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    for (int64_t k = 0; k < n; k++) {
        if (k == 0) {
            out[0] = expr_copy(r);
            continue;
        }
        /* The k-th n-th root of unity is Exp[2 k Pi I / n] = (-1)^(2k/n).
         * Routing through Power[-1, 2k/n] rather than Exp[...] avoids
         * Euler's-formula expansion into Cos + I Sin: Power's existing
         * (-1)^(p/q) canonicalisation produces the clean form
         * `sign * (-1)^(r/q)`, and Times merges those Power[-1, _] factors
         * with the principal root r = (-b/a)^(1/n) when r itself is a
         * Power[-1, 1/n].  Hence Solve[x^5 + 1 == 0, x] emits roots as
         * (-1)^(1/5), (-1)^(3/5), -1, -(-1)^(2/5), -(-1)^(4/5) (modulo
         * canonical ordering) instead of trigonometric sums. */
        Expr* exp_arg = (2 * k) % n == 0
            ? mk_int((2 * k) / n)
            : eval_and_free(mk_rat(2 * k, n));
        Expr* expk = eval_and_free(mk_pow(mk_int(-1), exp_arg));
        out[k] = eval_and_free(mk_fn2("Times", expr_copy(r), expk));
    }
    expr_free(r);
    *out_n = (size_t)n;
    return out;
}

/* a*x^(2n) + b*x^n + c == 0.  Substitute u = x^n, then for each
 * inner-quadratic root u_i recover n binomial roots of x^n − u_i. */
static Expr** solve_nquadratic(Expr* a, Expr* b, Expr* c, int64_t n,
                               bool reals_only, size_t* out_n) {
    size_t un = 0;
    Expr** us = solve_quadratic(a, b, c, reals_only, &un);
    size_t cap = un * (size_t)n;
    Expr** out = cap ? (Expr**)malloc(sizeof(Expr*) * cap) : NULL;
    size_t cnt = 0;
    for (size_t i = 0; i < un; i++) {
        Expr* one = mk_int(1);
        Expr* neg_u = mk_neg(expr_copy(us[i]));
        size_t kc = 0;
        Expr** kr = solve_binomial(one, neg_u, n, reals_only, &kc);
        expr_free(one); expr_free(neg_u);
        for (size_t k = 0; k < kc; k++) out[cnt++] = kr[k];
        free(kr);
    }
    for (size_t i = 0; i < un; i++) expr_free(us[i]);
    free(us);
    *out_n = cnt;
    return out;
}

/* Cardano: a*x^3 + b*x^2 + c*x + d == 0.
 *
 *   Δ₀ = b² − 3 a c
 *   Δ₁ = 2 b³ − 9 a b c + 27 a² d
 *   C  = ((Δ₁ + Sqrt[Δ₁² − 4 Δ₀³]) / 2)^(1/3)
 *   ζ  = Exp[2 π I / 3]
 *
 *   x_k = -1/(3a) (b + ζ^k C + Δ₀ / (ζ^k C))    for k = 0, 1, 2.
 *
 * Reals filter uses the cubic discriminant
 *   Δ = 18 a b c d − 4 b³ d + b² c² − 4 a c³ − 27 a² d²
 * (Δ < 0 ⇒ exactly one real root; we keep x_0 only.) */
static Expr** solve_cubic_radical(Expr* a, Expr* b, Expr* c, Expr* d,
                                  bool reals_only, size_t* out_n) {
    /* Δ₀ */
    Expr* D0 = eval_and_free(mk_fn2("Plus",
        eval_and_free(mk_pow(expr_copy(b), mk_int(2))),
        mk_neg(eval_and_free(mk_fn3("Times", mk_int(3),
                                              expr_copy(a),
                                              expr_copy(c))))));

    /* Δ₁ */
    Expr* term1 = eval_and_free(mk_fn2("Times", mk_int(2),
        eval_and_free(mk_pow(expr_copy(b), mk_int(3)))));
    Expr* term2 = eval_and_free(mk_fn4("Times", mk_int(9),
                                                expr_copy(a),
                                                expr_copy(b),
                                                expr_copy(c)));
    Expr* term3 = eval_and_free(mk_fn3("Times", mk_int(27),
        eval_and_free(mk_pow(expr_copy(a), mk_int(2))),
        expr_copy(d)));
    Expr* D1 = eval_and_free(mk_fn3("Plus", term1, mk_neg(term2), term3));

    /* inner = Δ₁² − 4 Δ₀³ */
    Expr* d1_sq = eval_and_free(mk_pow(expr_copy(D1), mk_int(2)));
    Expr* four_d0_3 = eval_and_free(mk_fn2("Times", mk_int(4),
        eval_and_free(mk_pow(expr_copy(D0), mk_int(3)))));
    Expr* inner = eval_and_free(mk_fn2("Plus", d1_sq, mk_neg(four_d0_3)));
    Expr* sqrt_inner = eval_and_free(mk_sqrt(inner));

    /* C_arg = (Δ₁ + sqrt_inner) / 2 */
    Expr* C_arg = eval_and_free(mk_fn2("Times",
        mk_fn2("Plus", expr_copy(D1), expr_copy(sqrt_inner)),
        mk_rat(1, 2)));

    /* If C_arg == 0, swap to the (Δ₁ - sqrt_inner)/2 branch to avoid the
     * triple-root degeneracy.  If still zero, the triple-root case
     * applies: x = -b/(3a) thrice. */
    bool triple_root = false;
    Expr* C_val;
    if (is_int_zero(C_arg)) {
        Expr* alt = eval_and_free(mk_fn2("Times",
            mk_fn2("Plus", expr_copy(D1), mk_neg(expr_copy(sqrt_inner))),
            mk_rat(1, 2)));
        expr_free(C_arg);
        C_arg = alt;
        if (is_int_zero(C_arg)) triple_root = true;
    }
    C_val = triple_root
        ? mk_int(0)
        : eval_and_free(mk_pow(expr_copy(C_arg), mk_rat(1, 3)));
    expr_free(C_arg);

    /* ζ = Exp[2 π I / 3] */
    Expr* zeta = eval_and_free(mk_fn1("Exp",
        mk_fn3("Times", mk_rat(2, 3), mk_sym("Pi"), mk_sym("I"))));

    /* inv_3a = 1 / (3 a) */
    Expr* inv_3a = eval_and_free(mk_inv(
        mk_fn2("Times", mk_int(3), expr_copy(a))));

    Expr** out_arr = (Expr**)malloc(sizeof(Expr*) * 3);
    for (int k = 0; k < 3; k++) {
        Expr* zk = (k == 0)
            ? mk_int(1)
            : eval_and_free(mk_pow(expr_copy(zeta), mk_int(k)));

        Expr* zkC;
        Expr* D0_term;
        if (triple_root) {
            zkC = mk_int(0);
            D0_term = mk_int(0);
            expr_free(zk);
        } else {
            zkC = eval_and_free(mk_fn2("Times", zk, expr_copy(C_val)));
            if (is_int_zero(D0)) {
                D0_term = mk_int(0);
            } else {
                D0_term = eval_and_free(mk_fn2("Times",
                    expr_copy(D0), mk_inv(expr_copy(zkC))));
            }
        }

        Expr* sum_in = mk_fn3("Plus", expr_copy(b), zkC, D0_term);
        Expr* xk = eval_and_free(mk_fn3("Times", mk_int(-1),
                                                  expr_copy(inv_3a),
                                                  sum_in));
        out_arr[k] = xk;
    }

    expr_free(zeta);
    expr_free(C_val);
    expr_free(D0);
    expr_free(D1);
    expr_free(sqrt_inner);
    expr_free(inv_3a);

    if (reals_only) {
        /* cubic discriminant Δ = 18 abcd − 4 b³d + b²c² − 4 a c³ − 27 a²d² */
        Expr* t1 = eval_and_free(mk_fn5("Times", mk_int(18), expr_copy(a),
                                                  expr_copy(b), expr_copy(c),
                                                  expr_copy(d)));
        Expr* t2 = eval_and_free(mk_fn3("Times", mk_int(4),
            eval_and_free(mk_pow(expr_copy(b), mk_int(3))),
            expr_copy(d)));
        Expr* t3 = eval_and_free(mk_fn2("Times",
            eval_and_free(mk_pow(expr_copy(b), mk_int(2))),
            eval_and_free(mk_pow(expr_copy(c), mk_int(2)))));
        Expr* t4 = eval_and_free(mk_fn3("Times", mk_int(4),
            expr_copy(a),
            eval_and_free(mk_pow(expr_copy(c), mk_int(3)))));
        Expr* t5 = eval_and_free(mk_fn3("Times", mk_int(27),
            eval_and_free(mk_pow(expr_copy(a), mk_int(2))),
            eval_and_free(mk_pow(expr_copy(d), mk_int(2)))));
        Expr* disc = eval_and_free(mk_fn5("Plus",
            t1, mk_neg(t2), t3, mk_neg(t4), mk_neg(t5)));
        int s = try_sign(disc);
        expr_free(disc);
        if (s == -1) {
            /* Only x_0 is real; the other two are complex conjugates. */
            expr_free(out_arr[1]);
            expr_free(out_arr[2]);
            Expr** out = (Expr**)malloc(sizeof(Expr*) * 1);
            out[0] = out_arr[0];
            free(out_arr);
            *out_n = 1;
            return out;
        }
        /* Δ ≥ 0 or unresolved: keep all three branches.  In casus
         * irreducibilis (Δ > 0) Cardano's formula still yields the
         * three real roots, expressed via complex intermediaries. */
    }

    *out_n = 3;
    return out_arr;
}

/* ------------------------------------------------------------------ *
 *  Root[] object construction.                                        *
 * ------------------------------------------------------------------ */

/* Deep-copy `e`, replacing every occurrence of the interned symbol
 * `vname` with Slot[1].  Mirrors substitute_bvar_with_slot in
 * src/root.c. */
static Expr* sub_var_with_slot(Expr* e, const char* vname) {
    if (!e) return NULL;
    if (e->type == EXPR_SYMBOL && e->data.symbol == vname) {
        return mk_fn1("Slot", mk_int(1));
    }
    if (e->type != EXPR_FUNCTION) return expr_copy(e);
    size_t n = e->data.function.arg_count;
    Expr* head = sub_var_with_slot(e->data.function.head, vname);
    Expr** args = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) {
        args[i] = sub_var_with_slot(e->data.function.args[i], vname);
    }
    Expr* out = expr_new_function(head, args, n);
    free(args);
    return out;
}

/* Root[Function[poly_in_slot1], k]. */
static Expr* root_object(Expr* poly_in_var, Expr* var, int k) {
    Expr* poly_in_slot = sub_var_with_slot(poly_in_var, var->data.symbol);
    Expr* fn = mk_fn1("Function", poly_in_slot);
    return mk_fn2("Root", fn, mk_int(k));
}

/* ------------------------------------------------------------------ *
 *  Polynomial-shape detectors.                                        *
 * ------------------------------------------------------------------ */

/* gcd helper for the exponent-GCD substitution (mirrors Maxima's
 * solventhp).  Inputs may be negative; result is non-negative. */
static int64_t gcd_i64(int64_t a, int64_t b) {
    if (a < 0) a = -a;
    if (b < 0) b = -b;
    while (b != 0) { int64_t t = a % b; a = b; b = t; }
    return a;
}

/* Compute g = gcd of every var-exponent k > 0 at which `poly` has a
 * nonzero coefficient.  Returns true iff g >= 2 *and* the substituted
 * polynomial would still have degree d/g >= 2 (so this isn't the
 * binomial special case that solve_binomial already handles directly)
 * *and* at least two distinct nonzero exponents > 0 exist (otherwise
 * the polynomial is one term + maybe a constant, again a binomial).
 * Mirrors Maxima's solventhp. */
static bool find_exponent_gcd(Expr* poly, Expr* var, int d,
                              int64_t* g_out) {
    int64_t g = 0;
    int n_nonzero_pos = 0;
    for (int k = 1; k <= d; k++) {
        Expr* ck = get_coeff(poly, var, k);
        if (!ck) continue;
        bool zero = is_int_zero(ck);
        expr_free(ck);
        if (zero) continue;
        n_nonzero_pos++;
        g = (g == 0) ? (int64_t)k : gcd_i64(g, (int64_t)k);
        if (g == 1) return false;
    }
    if (g < 2) return false;
    if (n_nonzero_pos < 2) return false;
    if ((int64_t)d / g < 2) return false;
    *g_out = g;
    return true;
}

/* Forward declaration -- referenced by solve_gcd_substituted. */
Expr* solvepoly_solve_polynomial_equality(Expr* equation, Expr* var,
                                          Expr* dom,
                                          const SolvePolyOpts* opts);

/* For an inner-poly root `u_val`, push the g binomial roots of
 * x^g - u_val into `sl`.  When u_val == 0 the only root is x = 0,
 * pushed g times (the original poly had a factor of x^g at that root). */
static void push_binomial_x_roots(Expr* u_val, int64_t g, bool reals_only,
                                  int64_t mult, SolList* sl) {
    if (is_int_zero(u_val)) {
        for (int64_t i = 0; i < g; i++) {
            sl_push_with_mult(sl, mk_int(0), mult);
        }
        return;
    }
    Expr* one = mk_int(1);
    Expr* neg_u = mk_neg(expr_copy(u_val));
    size_t outc = 0;
    Expr** rs = solve_binomial(one, neg_u, g, reals_only, &outc);
    expr_free(one); expr_free(neg_u);
    sl_extend_with_mult(sl, rs, outc, mult);
    free(rs);
}

/* GCD-of-exponents substitution.  Builds u_poly = sum_j c_{j g} u^j
 * (reusing `var` as the dummy u-symbol), solves it via solvepoly
 * recursively, then for each u-root u_i pushes the binomial roots of
 * x^g - u_i into `sl`.  Multiplicity is carried both through duplicate
 * u-root entries (collated by solvepoly's multiplicity machinery) and
 * through the outer `mult` parameter from handle_factor.  Mirrors
 * Maxima's solventh. */
static bool solve_gcd_substituted(Expr* poly, Expr* var, int d, int64_t g,
                                  bool reals_only,
                                  const SolvePolyOpts* opts, Expr* dom,
                                  int64_t mult, SolList* sl) {
    int u_deg = d / (int)g;
    Expr** u_terms = (Expr**)malloc(sizeof(Expr*) * (u_deg + 1));
    size_t cnt = 0;
    for (int j = 0; j <= u_deg; j++) {
        Expr* cj = get_coeff(poly, var, j * (int)g);
        if (!cj) {
            for (size_t i = 0; i < cnt; i++) expr_free(u_terms[i]);
            free(u_terms);
            return false;
        }
        if (is_int_zero(cj)) { expr_free(cj); continue; }
        Expr* term;
        if (j == 0) {
            term = cj;
        } else if (j == 1) {
            term = mk_fn2("Times", cj, expr_copy(var));
        } else {
            term = mk_fn2("Times", cj,
                          mk_pow(expr_copy(var), mk_int(j)));
        }
        u_terms[cnt++] = term;
    }
    Expr* u_poly;
    if (cnt == 0) {
        u_poly = mk_int(0);
        free(u_terms);
    } else if (cnt == 1) {
        u_poly = u_terms[0];
        free(u_terms);
    } else {
        u_poly = expr_new_function(mk_sym("Plus"), u_terms, cnt);
    }
    u_poly = eval_and_free(u_poly);

    Expr* u_eq = eval_and_free(mk_fn2("Equal", u_poly, mk_int(0)));
    Expr* u_solutions = solvepoly_solve_polynomial_equality(
        u_eq, var, dom, opts);
    expr_free(u_eq);
    if (!u_solutions) return false;

    if (u_solutions->type != EXPR_FUNCTION
        || u_solutions->data.function.head->type != EXPR_SYMBOL
        || u_solutions->data.function.head->data.symbol != SYM_List) {
        expr_free(u_solutions);
        return false;
    }
    for (size_t i = 0; i < u_solutions->data.function.arg_count; i++) {
        Expr* sol = u_solutions->data.function.args[i];
        if (sol->type != EXPR_FUNCTION
            || sol->data.function.head->type != EXPR_SYMBOL
            || sol->data.function.head->data.symbol != SYM_List
            || sol->data.function.arg_count != 1) continue;
        Expr* rule = sol->data.function.args[0];
        if (rule->type != EXPR_FUNCTION
            || rule->data.function.head->type != EXPR_SYMBOL
            || rule->data.function.head->data.symbol != SYM_Rule
            || rule->data.function.arg_count != 2) continue;
        Expr* u_val = rule->data.function.args[1];
        push_binomial_x_roots(u_val, g, reals_only, mult, sl);
    }
    expr_free(u_solutions);
    return true;
}

/* True iff `e` contains the var symbol anywhere in its tree.  Used
 * by the polydecomp remainder check (every g_i must be free of var). */
static bool expr_contains_var(const Expr* e, const Expr* var) {
    if (!e || !var) return false;
    if (e->type == EXPR_SYMBOL && var->type == EXPR_SYMBOL
        && e->data.symbol == var->data.symbol) return true;
    if (e->type != EXPR_FUNCTION) return false;
    if (expr_contains_var(e->data.function.head, var)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (expr_contains_var(e->data.function.args[i], var)) return true;
    }
    return false;
}

/* Polynomial decomposition (Maxima `polydecomp` / `solve-by-
 * decomposition`), m = 2 inner-degree case.  When a monic polynomial
 * of even degree n >= 4 admits a decomposition f(x) = g(h(x)) with
 * h(x) = x^2 + h_1 x (h_0 absorbed into g), this routine recovers
 * h_1 directly from the second-highest coefficient (h^k binomial
 * expansion: [x^{n-1}] h^k = k h_1) and then uses repeated
 * polynomial division by h to read off the outer coefficients of
 * g(y) = sum g_i y^i.  Decomposition succeeds iff every remainder is
 * free of var and the top-degree coefficient of g is exactly 1
 * (matching the monic input).  On success g(y) and h(x) - y_i are
 * solved recursively through the polynomial specialist.
 *
 * This is the smallest decomposition class that *isn't* already
 * caught by Factor or the GCD-of-exponents substitution: it cracks
 * irreducible polynomials whose decomposition's inner shape has both
 * x^2 and x terms, e.g. the canonical example
 *   f(x) = (x^2 + x)^3 + 2 (x^2 + x) - 5
 *        = x^6 + 3 x^5 + 3 x^4 + x^3 + 2 x^2 + 2 x - 5
 * which factors as 6 monolithic Root[] objects without this pass and
 * as `(g of y = x^2 + x)(then quadratic per y_i)` with it. */
static bool solve_polydecomp_m2(Expr* poly, Expr* var, int d,
                                const SolvePolyOpts* opts, Expr* dom,
                                int64_t mult, SolList* sl) {
    if (d < 4 || (d % 2) != 0) return false;
    int k = d / 2;
    if (k < 2) return false;

    /* Require monic lead coefficient (extending to non-monic is just
     * a uniform scaling we don't need for the common cases). */
    Expr* a_lead = get_coeff(poly, var, d);
    bool monic = a_lead && a_lead->type == EXPR_INTEGER
                 && a_lead->data.integer == 1;
    if (a_lead) expr_free(a_lead);
    if (!monic) return false;

    /* h_1 = coeff_{n-1}(f) / k.  (Binomial expansion of h^k.) */
    Expr* c_nm1 = get_coeff(poly, var, d - 1);
    if (!c_nm1) return false;
    Expr* h1 = eval_and_free(mk_fn2("Times", c_nm1,
                                    mk_inv(mk_int(k))));
    Expr* h = eval_and_free(mk_fn2("Plus",
        mk_pow(expr_copy(var), mk_int(2)),
        mk_fn2("Times", h1, expr_copy(var))));

    /* Iteratively divide by h to read off g_0, g_1, ..., g_k. */
    Expr** g_coeffs = (Expr**)malloc(sizeof(Expr*) * (size_t)(k + 1));
    for (int i = 0; i <= k; i++) g_coeffs[i] = NULL;
    Expr* current = expr_copy(poly);
    bool ok = true;
    for (int i = 0; i <= k && ok; i++) {
        Expr* r = eval_and_free(internal_polynomialremainder(
            (Expr*[]){ expr_copy(current), expr_copy(h),
                       expr_copy(var) }, 3));
        if (!r || expr_contains_var(r, var)) {
            if (r) expr_free(r);
            ok = false;
            break;
        }
        g_coeffs[i] = r;
        if (i == k) break;
        Expr* q = eval_and_free(internal_polynomialquotient(
            (Expr*[]){ expr_copy(current), expr_copy(h),
                       expr_copy(var) }, 3));
        expr_free(current);
        current = q;
    }
    expr_free(current);
    if (!ok) {
        for (int i = 0; i <= k; i++) if (g_coeffs[i]) expr_free(g_coeffs[i]);
        free(g_coeffs);
        expr_free(h);
        return false;
    }

    /* Outer must be monic too (g_k == 1) for our normalisation. */
    if (!g_coeffs[k]
        || g_coeffs[k]->type != EXPR_INTEGER
        || g_coeffs[k]->data.integer != 1) {
        for (int i = 0; i <= k; i++) expr_free(g_coeffs[i]);
        free(g_coeffs);
        expr_free(h);
        return false;
    }

    /* Build g(y) using var as a temporary y-symbol. */
    Expr** g_terms = (Expr**)malloc(sizeof(Expr*) * (size_t)(k + 1));
    size_t cnt = 0;
    for (int i = 0; i <= k; i++) {
        if (is_int_zero(g_coeffs[i])) continue;
        Expr* term;
        if (i == 0) {
            term = expr_copy(g_coeffs[i]);
        } else if (i == 1) {
            term = mk_fn2("Times", expr_copy(g_coeffs[i]),
                          expr_copy(var));
        } else {
            term = mk_fn2("Times", expr_copy(g_coeffs[i]),
                          mk_pow(expr_copy(var), mk_int(i)));
        }
        g_terms[cnt++] = term;
    }
    for (int i = 0; i <= k; i++) expr_free(g_coeffs[i]);
    free(g_coeffs);
    Expr* g_poly;
    if (cnt == 1) {
        g_poly = g_terms[0];
        free(g_terms);
    } else {
        g_poly = expr_new_function(mk_sym("Plus"), g_terms, cnt);
    }
    g_poly = eval_and_free(g_poly);

    /* Solve g(y) == 0. */
    Expr* g_eq = eval_and_free(mk_fn2("Equal", g_poly, mk_int(0)));
    Expr* y_solutions = solvepoly_solve_polynomial_equality(
        g_eq, var, dom, opts);
    expr_free(g_eq);
    if (!y_solutions
        || y_solutions->type != EXPR_FUNCTION
        || y_solutions->data.function.head->type != EXPR_SYMBOL
        || y_solutions->data.function.head->data.symbol != SYM_List) {
        if (y_solutions) expr_free(y_solutions);
        expr_free(h);
        return false;
    }

    /* For each y_i, solve h(x) - y_i == 0 and push roots. */
    for (size_t i = 0; i < y_solutions->data.function.arg_count; i++) {
        Expr* sol = y_solutions->data.function.args[i];
        if (sol->type != EXPR_FUNCTION
            || sol->data.function.head->type != EXPR_SYMBOL
            || sol->data.function.head->data.symbol != SYM_List
            || sol->data.function.arg_count != 1) continue;
        Expr* rule = sol->data.function.args[0];
        if (rule->type != EXPR_FUNCTION
            || rule->data.function.head->type != EXPR_SYMBOL
            || rule->data.function.head->data.symbol != SYM_Rule
            || rule->data.function.arg_count != 2) continue;
        Expr* y_val = rule->data.function.args[1];
        Expr* h_eq = eval_and_free(mk_fn2("Equal",
            mk_fn2("Plus", expr_copy(h),
                   mk_neg(expr_copy(y_val))),
            mk_int(0)));
        Expr* x_solutions = solvepoly_solve_polynomial_equality(
            h_eq, var, dom, opts);
        expr_free(h_eq);
        if (!x_solutions) continue;
        if (x_solutions->type == EXPR_FUNCTION
            && x_solutions->data.function.head->type == EXPR_SYMBOL
            && x_solutions->data.function.head->data.symbol == SYM_List) {
            for (size_t j = 0; j < x_solutions->data.function.arg_count; j++) {
                Expr* xsol = x_solutions->data.function.args[j];
                if (xsol->type != EXPR_FUNCTION
                    || xsol->data.function.head->type != EXPR_SYMBOL
                    || xsol->data.function.head->data.symbol != SYM_List
                    || xsol->data.function.arg_count != 1) continue;
                Expr* xrule = xsol->data.function.args[0];
                if (xrule->type != EXPR_FUNCTION
                    || xrule->data.function.head->type != EXPR_SYMBOL
                    || xrule->data.function.head->data.symbol != SYM_Rule
                    || xrule->data.function.arg_count != 2) continue;
                Expr* x_val = xrule->data.function.args[1];
                sl_push_with_mult(sl, expr_copy(x_val), mult);
            }
        }
        expr_free(x_solutions);
    }
    expr_free(y_solutions);
    expr_free(h);
    return true;
}

/* Detect a*var^d + b == 0 (only the top and constant coefficients
 * nonzero).  On true, returns freshly-owned a, b and the degree n=d. */
static bool is_binomial(Expr* poly, Expr* var, int d,
                        Expr** a_out, Expr** b_out, int64_t* n_out) {
    if (d < 2) return false;
    Expr* a = get_coeff(poly, var, d);
    if (!a || is_int_zero(a)) { expr_free(a); return false; }
    Expr* b = get_coeff(poly, var, 0);
    if (!b) { expr_free(a); return false; }
    for (int k = 1; k < d; k++) {
        Expr* ck = get_coeff(poly, var, k);
        bool zero = is_int_zero(ck);
        expr_free(ck);
        if (!zero) { expr_free(a); expr_free(b); return false; }
    }
    *a_out = a; *b_out = b; *n_out = d;
    return true;
}

/* Detect a*var^(2n) + b*var^n + c == 0 (n ≥ 2). */
static bool is_nquadratic(Expr* poly, Expr* var, int d,
                          Expr** a_out, Expr** b_out, Expr** c_out,
                          int64_t* n_out) {
    if (d < 4 || (d % 2) != 0) return false;
    int n = d / 2;
    Expr* a = get_coeff(poly, var, d);
    if (!a || is_int_zero(a)) { expr_free(a); return false; }
    Expr* b = get_coeff(poly, var, n);
    if (!b || is_int_zero(b)) { expr_free(a); expr_free(b); return false; }
    Expr* c = get_coeff(poly, var, 0);
    if (!c) { expr_free(a); expr_free(b); return false; }
    for (int k = 1; k <= d - 1; k++) {
        if (k == n) continue;
        Expr* ck = get_coeff(poly, var, k);
        bool zero = is_int_zero(ck);
        expr_free(ck);
        if (!zero) { expr_free(a); expr_free(b); expr_free(c); return false; }
    }
    *a_out = a; *b_out = b; *c_out = c; *n_out = n;
    return true;
}

/* ------------------------------------------------------------------ *
 *  Per-factor dispatcher.                                             *
 * ------------------------------------------------------------------ */

static bool handle_factor(Expr* g, Expr* var, bool reals_only,
                          const SolvePolyOpts* opts, Expr* dom,
                          int64_t mult, SolList* sl) {
    int dg = get_degree_poly(g, var);
    if (dg <= 0) return true;

    if (dg == 1) {
        Expr* a = get_coeff(g, var, 1);
        Expr* b = get_coeff(g, var, 0);
        if (!a || !b) { expr_free(a); expr_free(b); return false; }
        Expr* root = solve_linear(a, b);
        expr_free(a); expr_free(b);
        sl_push_with_mult(sl, root, mult);
        return true;
    }
    if (dg == 2) {
        Expr* a = get_coeff(g, var, 2);
        Expr* b = get_coeff(g, var, 1);
        Expr* c = get_coeff(g, var, 0);
        if (!a || !b || !c) { expr_free(a); expr_free(b); expr_free(c); return false; }
        /* If the constant term is zero, the polynomial has var as a
         * factor: a*x^2 + b*x = x*(a*x + b).  Split into the trivial
         * root x = 0 and the linear root x = -b/a directly -- otherwise
         * the quadratic formula produces a spurious Sqrt[b^2]
         * discriminant that does not simplify for symbolic b. */
        if (is_int_zero(c)) {
            expr_free(c);
            sl_push_with_mult(sl, mk_int(0), mult);
            Expr* root = solve_linear(a, b);
            expr_free(a); expr_free(b);
            sl_push_with_mult(sl, root, mult);
            return true;
        }
        size_t n = 0;
        Expr** rs = solve_quadratic(a, b, c, reals_only, &n);
        expr_free(a); expr_free(b); expr_free(c);
        sl_extend_with_mult(sl, rs, n, mult);
        free(rs);
        return true;
    }

    { Expr *a, *b; int64_t n;
      if (is_binomial(g, var, dg, &a, &b, &n)) {
          size_t outc = 0;
          Expr** rs = solve_binomial(a, b, n, reals_only, &outc);
          expr_free(a); expr_free(b);
          sl_extend_with_mult(sl, rs, outc, mult);
          free(rs);
          return true;
      } }

    { Expr *a, *b, *c; int64_t n;
      if (is_nquadratic(g, var, dg, &a, &b, &c, &n)) {
          size_t outc = 0;
          Expr** rs = solve_nquadratic(a, b, c, n, reals_only, &outc);
          expr_free(a); expr_free(b); expr_free(c);
          sl_extend_with_mult(sl, rs, outc, mult);
          free(rs);
          return true;
      } }

    /* GCD-of-exponents substitution (Maxima solventh).  Activates
     * when every nonzero coefficient lies at an exponent that is a
     * multiple of some g >= 2, generalising the binomial / nquadratic
     * fast paths to any sparsely-supported polynomial. */
    { int64_t gc = 0;
      if (find_exponent_gcd(g, var, dg, &gc)) {
          if (solve_gcd_substituted(g, var, dg, gc, reals_only, opts,
                                    dom, mult, sl)) {
              return true;
          }
          /* If the inner recursive solve failed we fall through to the
           * cubic/Root[] paths below rather than aborting. */
      } }

    if (dg == 3 && opts->cubics_radical) {
        Expr* a = get_coeff(g, var, 3);
        Expr* b = get_coeff(g, var, 2);
        Expr* c = get_coeff(g, var, 1);
        Expr* d = get_coeff(g, var, 0);
        if (!a || !b || !c || !d) {
            expr_free(a); expr_free(b); expr_free(c); expr_free(d);
            return false;
        }
        size_t outc = 0;
        Expr** rs = solve_cubic_radical(a, b, c, d, reals_only, &outc);
        expr_free(a); expr_free(b); expr_free(c); expr_free(d);
        sl_extend_with_mult(sl, rs, outc, mult);
        free(rs);
        return true;
    }

    /* Polynomial decomposition last-ditch (Maxima solve-by-
     * decomposition, m=2 inner-degree case).  Only fires on degree
     * >= 6 -- smaller polynomials are either already cracked by the
     * fast paths above or are not interesting to decompose. */
    if (dg >= 6
        && solve_polydecomp_m2(g, var, dg, opts, dom, mult, sl)) {
        return true;
    }

    /* Default: dg Root[] objects per unit of multiplicity. */
    for (int64_t m = 0; m < mult; m++) {
        for (int k = 1; k <= dg; k++) {
            sl_push(sl, root_object(g, var, k));
        }
    }
    return true;
}

/* ------------------------------------------------------------------ *
 *  Factored-product walker.                                           *
 * ------------------------------------------------------------------ */

static bool walk_product(Expr* prod, Expr* var, bool reals_only,
                         const SolvePolyOpts* opts, Expr* dom,
                         int64_t outer_mult, SolList* sl);

static bool walk_one_factor(Expr* f, Expr* var, bool reals_only,
                            const SolvePolyOpts* opts, Expr* dom,
                            int64_t outer_mult, SolList* sl) {
    if (f->type == EXPR_FUNCTION
        && f->data.function.head->type == EXPR_SYMBOL
        && f->data.function.head->data.symbol == SYM_Power
        && f->data.function.arg_count == 2
        && f->data.function.args[1]->type == EXPR_INTEGER
        && f->data.function.args[1]->data.integer > 0) {
        int64_t e = f->data.function.args[1]->data.integer;
        return walk_one_factor(f->data.function.args[0], var, reals_only,
                               opts, dom, outer_mult * e, sl);
    }
    if (f->type == EXPR_FUNCTION
        && f->data.function.head->type == EXPR_SYMBOL
        && f->data.function.head->data.symbol == SYM_Times) {
        return walk_product(f, var, reals_only, opts, dom, outer_mult, sl);
    }
    return handle_factor(f, var, reals_only, opts, dom, outer_mult, sl);
}

static bool walk_product(Expr* prod, Expr* var, bool reals_only,
                         const SolvePolyOpts* opts, Expr* dom,
                         int64_t outer_mult, SolList* sl) {
    if (prod->type != EXPR_FUNCTION
        || prod->data.function.head->type != EXPR_SYMBOL
        || prod->data.function.head->data.symbol != SYM_Times) {
        return walk_one_factor(prod, var, reals_only, opts, dom,
                               outer_mult, sl);
    }
    for (size_t i = 0; i < prod->data.function.arg_count; i++) {
        if (!walk_one_factor(prod->data.function.args[i], var,
                             reals_only, opts, dom, outer_mult, sl)) {
            return false;
        }
    }
    return true;
}

/* ------------------------------------------------------------------ *
 *  Solution ordering.                                                 *
 *                                                                    *
 *  Mathematica's `Solve` orders its roots in a way that interleaves   *
 *  Power[-1, p/q] and Times[-1, Power[-1, p/q]] forms by the value of *
 *  the (-1)-exponent modulo 1.  For x^5 + 1 == 0 this yields          *
 *      -1, (-1)^(1/5), -(-1)^(2/5), (-1)^(3/5), -(-1)^(4/5)           *
 *  (canonical exponents 0/5, 1/5, 2/5, 3/5, 4/5).  Mathilda's         *
 *  canonical Sort groups all Power[..] before all Times[-1, Power[..]] *
 *  and therefore does not reproduce this ordering -- we apply a       *
 *  domain-specific comparator here.                                   *
 *                                                                    *
 *  Rule:                                                              *
 *    1. Concrete real numbers (Integer/Rational/Real/BigInt) sort     *
 *       first, by numerical value.                                    *
 *    2. Other expressions are decomposed as v = mag * (-1)^exp by     *
 *       absorbing every Integer(-1) and Power[-1, e] factor into the  *
 *       exponent.  The remaining factors form `mag`.  Comparison is   *
 *       (mag, exp mod 1, exp descending), with `expr_compare` as the  *
 *       final tiebreak.                                               *
 * ------------------------------------------------------------------ */

/* True for Integer/Rational/Real/BigInt.  Complex[_, _] and other     *
 * compound numerics are deliberately excluded so they go through the  *
 * (mag, exp)-decomposition path. */
static bool is_concrete_real(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL
        || e->type == EXPR_BIGINT) {
        return true;
    }
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Rational
        && e->data.function.arg_count == 2) {
        return true;
    }
    return false;
}

/* Decompose `v` as v = mag * (-1)^exp.  Walks v's factors:              *
 *   Integer(-1)             contributes exp += 1.                       *
 *   Power[-1, e]            contributes exp += e.                       *
 *   Everything else         is multiplied into mag.                     *
 *                                                                       *
 * On return, *mag_out and *exp_out are freshly allocated and owned by   *
 * the caller. */
static void decompose_root(const Expr* v, Expr** mag_out, Expr** exp_out) {
    Expr* const* factors;
    size_t fc;
    Expr* single_array[1] = { (Expr*)v };
    if (v->type == EXPR_FUNCTION
        && v->data.function.head->type == EXPR_SYMBOL
        && v->data.function.head->data.symbol == SYM_Times) {
        factors = v->data.function.args;
        fc = v->data.function.arg_count;
    } else {
        factors = single_array;
        fc = 1;
    }

    Expr* exp = mk_int(0);
    Expr** mag_args = (Expr**)malloc(sizeof(Expr*) * (fc ? fc : 1));
    size_t mag_count = 0;

    for (size_t i = 0; i < fc; i++) {
        const Expr* f = factors[i];
        if (f->type == EXPR_INTEGER && f->data.integer == -1) {
            exp = eval_and_free(mk_fn2("Plus", exp, mk_int(1)));
            continue;
        }
        if (f->type == EXPR_FUNCTION
            && f->data.function.head->type == EXPR_SYMBOL
            && f->data.function.head->data.symbol == SYM_Power
            && f->data.function.arg_count == 2
            && f->data.function.args[0]->type == EXPR_INTEGER
            && f->data.function.args[0]->data.integer == -1) {
            exp = eval_and_free(mk_fn2("Plus", exp,
                                       expr_copy(f->data.function.args[1])));
            continue;
        }
        mag_args[mag_count++] = expr_copy((Expr*)f);
    }

    Expr* mag;
    if (mag_count == 0) {
        mag = mk_int(1);
    } else if (mag_count == 1) {
        mag = mag_args[0];
    } else {
        mag = expr_new_function(mk_sym("Times"), mag_args, mag_count);
        mag = eval_and_free(mag);
    }
    free(mag_args);

    *mag_out = mag;
    *exp_out = exp;
}

/* Return `exp mod 1` as an Integer or Rational in [0, 1).  Falls back  *
 * to evaluator-driven `exp - Floor[exp]` for shapes we don't recognise. */
static Expr* exp_mod_one(const Expr* exp) {
    if (exp->type == EXPR_INTEGER || exp->type == EXPR_BIGINT) {
        return mk_int(0);
    }
    int64_t p, q;
    if (is_rational((Expr*)exp, &p, &q)) {
        int64_t r = p % q;
        if (r < 0) r += q;
        if (r == 0) return mk_int(0);
        return eval_and_free(mk_rat(r, q));
    }
    Expr* floored = eval_and_free(mk_fn1("Floor", expr_copy((Expr*)exp)));
    return eval_and_free(mk_fn2("Plus", expr_copy((Expr*)exp), mk_neg(floored)));
}

typedef struct {
    Expr* val;        /* borrowed -- points into sl.vals[] */
    Expr* mag;        /* owned */
    Expr* exp;        /* owned, the raw exponent */
    Expr* exp_mod;    /* owned, exp mod 1 */
    bool  concrete;   /* true => sort by numeric value of val alone */
} RootKey;

static int root_key_qsort_cmp(const void* pa, const void* pb) {
    const RootKey* ka = (const RootKey*)pa;
    const RootKey* kb = (const RootKey*)pb;

    if (ka->concrete && kb->concrete) {
        int c = expr_compare(ka->val, kb->val);
        if (c != 0) return c;
        return 0;
    }
    if (ka->concrete) return -1;
    if (kb->concrete) return 1;

    int c = expr_compare(ka->mag, kb->mag);
    if (c != 0) return c;
    c = expr_compare(ka->exp_mod, kb->exp_mod);
    if (c != 0) return c;
    /* Raw exponent descending: (-1)^1 (== -1) before (-1)^0 (== 1)
     * when both collapse to exp mod 1 == 0. */
    c = expr_compare(kb->exp, ka->exp);
    if (c != 0) return c;

    return expr_compare(ka->val, kb->val);
}

/* Reorder sl->vals[] in place so that successive roots appear in the    *
 * canonical Solve order described in the comment above.  Stable across  *
 * equal keys: when two roots share every key, the `expr_compare` final  *
 * tiebreak in `root_key_qsort_cmp` keeps duplicates adjacent so that    *
 * multiplicity-pushed copies stay grouped. */
static void sl_sort_canonical(SolList* sl) {
    if (sl->count <= 1) return;

    RootKey* keys = (RootKey*)malloc(sizeof(RootKey) * sl->count);
    for (size_t i = 0; i < sl->count; i++) {
        keys[i].val = sl->vals[i];
        keys[i].concrete = is_concrete_real(sl->vals[i]);
        if (keys[i].concrete) {
            keys[i].mag = NULL;
            keys[i].exp = NULL;
            keys[i].exp_mod = NULL;
        } else {
            decompose_root(sl->vals[i], &keys[i].mag, &keys[i].exp);
            keys[i].exp_mod = exp_mod_one(keys[i].exp);
        }
    }

    qsort(keys, sl->count, sizeof(RootKey), root_key_qsort_cmp);

    for (size_t i = 0; i < sl->count; i++) {
        sl->vals[i] = keys[i].val;
        expr_free(keys[i].mag);
        expr_free(keys[i].exp);
        expr_free(keys[i].exp_mod);
    }
    free(keys);
}

/* ------------------------------------------------------------------ *
 *  Public entry.                                                      *
 * ------------------------------------------------------------------ */

Expr* solvepoly_solve_polynomial_equality(Expr* equation,
                                          Expr* var,
                                          Expr* dom,
                                          const SolvePolyOpts* opts) {
    if (!equation || !var) return NULL;
    if (var->type != EXPR_SYMBOL) return NULL;
    if (equation->type != EXPR_FUNCTION
        || equation->data.function.head->type != EXPR_SYMBOL
        || equation->data.function.head->data.symbol != SYM_Equal
        || equation->data.function.arg_count != 2) {
        return NULL;
    }

    /* Domain handling.  `Integers` is implemented as "solve in Reals
     * then drop every candidate that is not a concrete Integer / BigInt"
     * -- a coarse but sound strategy: it never invents a non-integer
     * solution, and the Reals filter (discriminant tests, real-branch
     * binomial selection) already prunes everything that is provably
     * complex.  The slow path (FactorSquareFree → Factor → linear
     * roots) handles rational-root polynomials such as
     * `x^3 - 6 x^2 + 11 x - 6` correctly because each irreducible
     * factor is degree 1 and emits an integer rule.  Irreducible
     * factors of degree ≥ 3 (default Cubics/Quartics → False) become
     * held Root[] objects which never satisfy the integer test, so we
     * silently drop them -- the user can opt into radical output via
     * `Cubics -> True` to recover those integer solutions that
     * collapse to an Integer after Cardano's substitution. */
    bool reals_only = false;
    bool integers_only = false;
    if (dom && dom->type == EXPR_SYMBOL) {
        const char* dsym = dom->data.symbol;
        if (dsym == SYM_Reals) {
            reals_only = true;
        } else if (dsym == SYM_Integers) {
            reals_only = true;
            integers_only = true;
        } else if (dsym != SYM_Complexes) {
            /* Rationals / Algebraics / Primes / Booleans etc. are not
             * yet wired up; leave the call unevaluated. */
            return NULL;
        }
    }

    Expr* lhs = equation->data.function.args[0];
    Expr* rhs = equation->data.function.args[1];

    /* Canonicalise the equation before polynomial dispatch.
     *
     *     Together both sides:   lhs → N1/D1,  rhs → N2/D2
     *     Cross-multiply:        poly = N1*D2 - N2*D1
     *
     * This brings rational inputs like `a/x + b == 0` into polynomial
     * form (here: a + b*x == 0).  When both sides are already
     * polynomial in `var`, Together leaves them untouched and the
     * denominators are 1, so the formula degenerates to the previous
     * `lhs - rhs` -- existing tests are unaffected.
     *
     * Cross-multiplying may *introduce* extraneous roots at points
     * where the original denominator vanished; the filter at the end
     * of this function drops any candidate `var = r` for which
     * D1*D2 evaluates to a concrete zero. */
    Expr* tl = eval_and_free(internal_together(
        (Expr*[]){ expr_copy(lhs) }, 1));
    Expr* tr = eval_and_free(internal_together(
        (Expr*[]){ expr_copy(rhs) }, 1));
    Expr* n1 = eval_and_free(internal_numerator(
        (Expr*[]){ expr_copy(tl) }, 1));
    Expr* d1 = eval_and_free(internal_denominator(
        (Expr*[]){ expr_copy(tl) }, 1));
    Expr* n2 = eval_and_free(internal_numerator(
        (Expr*[]){ expr_copy(tr) }, 1));
    Expr* d2 = eval_and_free(internal_denominator(
        (Expr*[]){ expr_copy(tr) }, 1));
    expr_free(tl); expr_free(tr);

    Expr* poly = eval_and_free(mk_fn2("Plus",
        mk_fn2("Times", expr_copy(n1), expr_copy(d2)),
        mk_neg(mk_fn2("Times", expr_copy(n2), expr_copy(d1)))));

    if (!is_polynomial(poly, &var, 1)) {
        expr_free(poly);
        expr_free(n1); expr_free(d1);
        expr_free(n2); expr_free(d2);
        return NULL;
    }

    /* Build the denominator product only when at least one side
     * actually carried a `var`-dependent denominator -- otherwise
     * cross-multiplication did not introduce any extraneous roots
     * and the filter would be wasted work. */
    Expr* denom_prod = NULL;
    if (contains_any_symbol_from(d1, var)
        || contains_any_symbol_from(d2, var)) {
        denom_prod = eval_and_free(mk_fn2("Times",
            expr_copy(d1), expr_copy(d2)));
    }
    expr_free(n1); expr_free(d1);
    expr_free(n2); expr_free(d2);

    /* Collect terms wrt `var` so the polynomial appears in standard
     * form before the fast-path classifier inspects it.  This is the
     * "collect wrt x" step in the canonicalisation pipeline:
     * Together-then-cross-multiply can leave the polynomial as e.g.
     * `1 - 2*(x-1)`; Collect rewrites it as `3 - 2*x`. */
    Expr* collected = internal_collect(
        (Expr*[]){ expr_copy(poly), expr_copy(var) }, 2);
    expr_free(poly);
    poly = collected;

    /* get_degree_poly only recognises var-bare Power nodes, so factored
     * forms like (x-1)^2 report degree 0.  Expand once up front so the
     * fast-path classifier sees the true degree.  Any subsequent
     * FactorSquareFree / Factor in the slow path will re-detect the
     * factorisation and recover multiplicities. */
    Expr* expanded = internal_expand((Expr*[]){ expr_copy(poly) }, 1);
    expr_free(poly);
    poly = expanded;

    int d = get_degree_poly(poly, var);

    /* Hidden-zero coefficient pass.  After Collect+Expand the structural
     * degree can overstate the true degree when a leading coefficient is
     * a non-trivial zero -- e.g. Sqrt[5 + 2 Sqrt[6]] - Sqrt[3] - Sqrt[2],
     * which standard rational normalisation cannot collapse but
     * PossibleZeroQ recognises through Stage-2 numeric verification.
     * Walk from k = d downward, testing each coefficient with
     * `zero_test_decide` until we hit one we cannot prove zero -- that
     * index is the true polynomial degree.  When every coefficient
     * (including the constant) collapses to zero, the equation is a
     * tautology.  Otherwise rebuild the polynomial from the surviving
     * coefficients so the downstream fast-path classifier sees the
     * correct shape. */
    int true_d = -1;
    for (int k = d; k >= 0; k--) {
        Expr* ck = get_coeff(poly, var, k);
        if (!ck) { true_d = k; break; }
        if (is_int_zero(ck)) { expr_free(ck); continue; }
        bool zero = (zero_test_decide(ck) == ZERO_TEST_TRUE);
        expr_free(ck);
        if (!zero) { true_d = k; break; }
    }

    if (true_d < 0) {
        /* Every coefficient -- including the constant -- collapses to
         * zero: the equation is a tautology. */
        expr_free(poly);
        Expr** outer = (Expr**)malloc(sizeof(Expr*) * 1);
        outer[0] = mk_list(NULL, 0);
        return mk_list(outer, 1);
    }

    if (true_d < d) {
        /* Rebuild the polynomial from c_0 ... c_{true_d}.  Any hidden-
         * zero coefficient at or below true_d is forced to a structural
         * 0 so the polynomial-shape detectors below operate on a clean
         * canonical form. */
        Expr** terms = (Expr**)malloc(sizeof(Expr*) * (size_t)(true_d + 2));
        size_t cnt = 0;
        for (int k = 0; k <= true_d; k++) {
            Expr* ck = get_coeff(poly, var, k);
            if (!ck) continue;
            bool drop = is_int_zero(ck)
                || (k < true_d && zero_test_decide(ck) == ZERO_TEST_TRUE);
            if (drop) { expr_free(ck); continue; }
            Expr* term;
            if (k == 0) {
                term = ck;
            } else if (k == 1) {
                term = eval_and_free(mk_fn2("Times", ck, expr_copy(var)));
            } else {
                term = eval_and_free(mk_fn2("Times", ck,
                    eval_and_free(mk_pow(expr_copy(var), mk_int(k)))));
            }
            terms[cnt++] = term;
        }
        Expr* new_poly;
        if (cnt == 0) {
            new_poly = mk_int(0);
        } else if (cnt == 1) {
            new_poly = terms[0];
        } else {
            new_poly = eval_and_free(
                expr_new_function(mk_sym("Plus"), terms, cnt));
            terms = NULL;
        }
        free(terms);
        expr_free(poly);
        poly = new_poly;
        d = true_d;
    }

    if (d == 0) {
        bool zero = is_int_zero(poly);
        expr_free(poly);
        if (zero) {
            /* tautology: { {} } */
            Expr** outer = (Expr**)malloc(sizeof(Expr*) * 1);
            outer[0] = mk_list(NULL, 0);
            return mk_list(outer, 1);
        }
        return mk_list(NULL, 0);
    }

    SolList sl; sl_init(&sl);

    bool used_fast = false;
    if (d == 1 || d == 2) {
        if (!handle_factor(poly, var, reals_only, opts, dom, 1, &sl)) {
            sl_free(&sl); expr_free(poly); return NULL;
        }
        used_fast = true;
    } else {
        Expr *a, *b; int64_t n;
        if (is_binomial(poly, var, d, &a, &b, &n)) {
            expr_free(a); expr_free(b);
            if (!handle_factor(poly, var, reals_only, opts, dom, 1, &sl)) {
                sl_free(&sl); expr_free(poly); return NULL;
            }
            used_fast = true;
        } else {
            Expr *na, *nb, *nc; int64_t nn;
            if (is_nquadratic(poly, var, d, &na, &nb, &nc, &nn)) {
                expr_free(na); expr_free(nb); expr_free(nc);
                if (!handle_factor(poly, var, reals_only, opts, dom,
                                   1, &sl)) {
                    sl_free(&sl); expr_free(poly); return NULL;
                }
                used_fast = true;
            } else {
                /* Exponent-GCD fast path -- catches polynomials like
                 *   x^9 + x^6 + x^3 + 1  (gcd 3)
                 * that aren't binomial / nquadratic but still collapse
                 * cleanly under u = x^g. */
                int64_t gc = 0;
                if (find_exponent_gcd(poly, var, d, &gc)) {
                    if (!handle_factor(poly, var, reals_only, opts, dom,
                                       1, &sl)) {
                        sl_free(&sl); expr_free(poly); return NULL;
                    }
                    used_fast = true;
                }
            }
        }
    }

    if (!used_fast) {
        Expr* collected = internal_collect((Expr*[]){ expr_copy(poly),
                                                       expr_copy(var) }, 2);
        Expr* sqfree    = internal_factorsquarefree((Expr*[]){ collected }, 1);

        /* Decompose sqfree as a list of (factor, mult). */
        size_t cap = 1;
        if (sqfree->type == EXPR_FUNCTION
            && sqfree->data.function.head->type == EXPR_SYMBOL
            && sqfree->data.function.head->data.symbol == SYM_Times) {
            cap = sqfree->data.function.arg_count;
        }
        Expr** sf_factors = (Expr**)malloc(sizeof(Expr*) * cap);
        int64_t* sf_mults = (int64_t*)malloc(sizeof(int64_t) * cap);
        size_t sf_count = 0;

        if (sqfree->type == EXPR_FUNCTION
            && sqfree->data.function.head->type == EXPR_SYMBOL
            && sqfree->data.function.head->data.symbol == SYM_Times) {
            for (size_t i = 0; i < sqfree->data.function.arg_count; i++) {
                Expr* fi = sqfree->data.function.args[i];
                if (fi->type == EXPR_FUNCTION
                    && fi->data.function.head->type == EXPR_SYMBOL
                    && fi->data.function.head->data.symbol == SYM_Power
                    && fi->data.function.arg_count == 2
                    && fi->data.function.args[1]->type == EXPR_INTEGER
                    && fi->data.function.args[1]->data.integer > 0) {
                    sf_factors[sf_count] = expr_copy(fi->data.function.args[0]);
                    sf_mults[sf_count] = fi->data.function.args[1]->data.integer;
                } else {
                    sf_factors[sf_count] = expr_copy(fi);
                    sf_mults[sf_count] = 1;
                }
                sf_count++;
            }
        } else {
            if (sqfree->type == EXPR_FUNCTION
                && sqfree->data.function.head->type == EXPR_SYMBOL
                && sqfree->data.function.head->data.symbol == SYM_Power
                && sqfree->data.function.arg_count == 2
                && sqfree->data.function.args[1]->type == EXPR_INTEGER
                && sqfree->data.function.args[1]->data.integer > 0) {
                sf_factors[0] = expr_copy(sqfree->data.function.args[0]);
                sf_mults[0] = sqfree->data.function.args[1]->data.integer;
            } else {
                sf_factors[0] = expr_copy(sqfree);
                sf_mults[0] = 1;
            }
            sf_count = 1;
        }
        expr_free(sqfree);

        bool ok = true;
        for (size_t i = 0; i < sf_count && ok; i++) {
            if (get_degree_poly(sf_factors[i], var) <= 0) continue;
            Expr* factored = internal_factor(
                (Expr*[]){ expr_copy(sf_factors[i]) }, 1);
            ok = walk_product(factored, var, reals_only, opts, dom,
                              sf_mults[i], &sl);
            expr_free(factored);
        }
        for (size_t i = 0; i < sf_count; i++) expr_free(sf_factors[i]);
        free(sf_factors);
        free(sf_mults);

        if (!ok) {
            sl_free(&sl); expr_free(poly); return NULL;
        }
    }

    expr_free(poly);

    /* Extraneous-root filter.  Cross-multiplying by D1*D2 may have
     * introduced solutions at the poles of the original equation; we
     * drop any candidate `var = r` for which D1*D2 evaluates to a
     * concrete (provable) zero.  Symbolic / undetermined values are
     * kept -- consistent with Mathematica's "verify when you can"
     * default. */
    if (denom_prod) {
        size_t kept = 0;
        for (size_t i = 0; i < sl.count; i++) {
            Expr* rule = mk_rule(expr_copy(var), expr_copy(sl.vals[i]));
            Expr* sub = eval_and_free(internal_replace_all(
                (Expr*[]){ expr_copy(denom_prod), rule }, 2));
            bool extraneous = is_definite_zero(sub);
            expr_free(sub);
            if (extraneous) {
                expr_free(sl.vals[i]);
            } else {
                sl.vals[kept++] = sl.vals[i];
            }
        }
        sl.count = kept;
        expr_free(denom_prod);
    }

    /* Integers-domain filter.  After the standard solve + extraneous-
     * root pass, drop every candidate whose value is not a *provably*
     * concrete integer (EXPR_INTEGER or EXPR_BIGINT).  We never trust
     * a `Rational[p, q]`, `Power[_, Rational[_, _]]`, or symbolic
     * residue here -- the only safe answers in `Integers` are values
     * the system has already collapsed to integer form.  Mathematica's
     * `Solve[..., Integers]` exhibits the same "drop everything that
     * is not provably an integer" semantics. */
    if (integers_only) {
        size_t kept = 0;
        for (size_t i = 0; i < sl.count; i++) {
            Expr* v = sl.vals[i];
            if (v->type == EXPR_INTEGER || v->type == EXPR_BIGINT) {
                sl.vals[kept++] = v;
            } else {
                expr_free(v);
            }
        }
        sl.count = kept;
    }

    /* Re-order solutions into Mathematica's canonical Solve order
     * (concrete reals first by value, then by exponent-of-(-1) mod 1).
     * Applied here, *after* the extraneous-root filter, so the surviving
     * roots are exactly what the user sees. */
    sl_sort_canonical(&sl);

    /* Wrap each solution value into List[Rule[var, val]]. */
    Expr** outer = sl.count ? (Expr**)malloc(sizeof(Expr*) * sl.count) : NULL;
    for (size_t i = 0; i < sl.count; i++) {
        Expr* rule = mk_rule(expr_copy(var), sl.vals[i]);
        outer[i] = mk_list((Expr*[]){ rule }, 1);
    }
    free(sl.vals);
    return mk_list(outer, sl.count);
}

/* ------------------------------------------------------------------ *
 *  Builtin entry & init.                                              *
 * ------------------------------------------------------------------ */

Expr* builtin_solve_polynomial_equality(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 3) return NULL;
    Expr* equation = res->data.function.args[0];
    Expr* var      = res->data.function.args[1];
    Expr* dom      = (argc >= 3) ? res->data.function.args[2] : NULL;

    SolvePolyOpts opts = { false, false };
    return solvepoly_solve_polynomial_equality(equation, var, dom, &opts);
    /* evaluator frees res on non-NULL return */
}

void solvepoly_init(void) {
    symtab_add_builtin("Solve`SolvePolynomialEquality",
                       builtin_solve_polynomial_equality);
    SymbolDef* def = symtab_get_def("Solve`SolvePolynomialEquality");
    if (def) def->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Solve`SolvePolynomialEquality",
        "Solve`SolvePolynomialEquality[lhs == rhs, var]\n"
        "Solve`SolvePolynomialEquality[lhs == rhs, var, dom]\n"
        "\tThe polynomial-equality specialist used by Solve.  Returns\n"
        "\tthe solutions of the given polynomial equation in `var` as a\n"
        "\tList of singleton-rule Lists, e.g. {{x -> 2}, {x -> 3}}.\n"
        "\tdom = Complexes (default) returns every complex root.\n"
        "\tdom = Reals returns only real roots (per-degree discriminant\n"
        "\tand sign-aware branch selection).\n"
        "\tdom = Integers solves over the reals and then keeps only the\n"
        "\troots that are provably concrete integers; non-integer real,\n"
        "\trational, irrational, and Root[] outputs are dropped.\n"
        "\tCubic and quartic factors are emitted as held Root[] objects\n"
        "\tby default; the Cubics and Quartics options on the parent\n"
        "\tSolve switch to closed-form radical output where supported.");
}
