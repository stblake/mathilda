/* radicals.c
 *
 * ToRadicals: convert held Root[Function[poly], k] objects into closed-form
 * radical expressions.  See radicals.h for the public contract.
 *
 * Algorithm (per Root node):
 *   1. Extract the polynomial body from Root[Function[..], k].  Both the
 *      Slot[1] form `Function[expr]` and the bound-variable form
 *      `Function[t, expr]` are accepted.
 *   2. Substitute Slot[1] (or t) with a fresh symbol `x$` so the existing
 *      get_coeff / get_degree_poly polynomial machinery operates on a
 *      standard univariate polynomial.
 *   3. Dispatch on degree d:
 *        d == 1 : linear, x = -c0/c1
 *        d == 2 : quadratic formula
 *        d == 3 : Cardano
 *        d == 4 : Ferrari (depressed quartic + resolvent cubic)
 *        d >= 5 : binomial fast-path a x^n + b only; otherwise leave the
 *                 Root untouched.
 *      Each path produces ALL d radical roots as a freshly-owned Expr**.
 *   4. Select the k-th root in Mathilda's canonical Root ordering by
 *      computing N[Root[poly, k]] at machine precision (via
 *      root_numericalize) and picking the radical root whose numeric
 *      value lies closest in the complex plane.
 *      When numeric evaluation is unavailable (the polynomial has
 *      parametric coefficients), fall back to the natural per-formula
 *      order with k - 1 as the index.
 *
 * Threading: the top-level walker is a structural recurrence that
 * reconstructs every EXPR_FUNCTION node it visits, so a Root buried
 * inside List, Equal, Less, And, Or, ... is processed identically.
 *
 * Memory: every internal helper returns a freshly-owned Expr*; inputs
 * are borrowed and deep-copied wherever they appear in the output.  The
 * exact-arithmetic core (eval_and_free, Plus/Times/Power normalisation)
 * does the bookkeeping for intermediate trees.
 */

#include "radicals.h"

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "attr.h"
#include "common.h"
#include "eval.h"
#include "expr.h"
#include "numeric.h"
#include "poly/poly.h"
#include "root_numeric.h"
#include "sym_names.h"
#include "symtab.h"

/* ------------------------------------------------------------------ *
 *  Expression-building shorthand.                                     *
 *                                                                    *
 *  Each helper takes ownership of its pointer arguments; the caller   *
 *  is expected to pass freshly built or expr_copy()'d Expr*'s.         *
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

static Expr* mk_pow(Expr* base, Expr* exp) { return mk_fn2("Power", base, exp); }
static Expr* mk_inv(Expr* e) { return mk_pow(e, mk_int(-1)); }
static Expr* mk_neg(Expr* e) { return mk_fn2("Times", mk_int(-1), e); }
static Expr* mk_rat(int64_t p, int64_t q) {
    return mk_fn2("Rational", mk_int(p), mk_int(q));
}
static Expr* mk_sqrt(Expr* e) { return mk_pow(e, mk_rat(1, 2)); }

static bool is_int_zero(const Expr* e) {
    return e && e->type == EXPR_INTEGER && e->data.integer == 0;
}

/* ------------------------------------------------------------------ *
 *  Polynomial extraction.                                             *
 * ------------------------------------------------------------------ */

/* Deep-copy `e`, replacing every Slot[1] occurrence with a deep copy of
 * `xvar` (a fresh dummy symbol). */
static Expr* slot1_to_var(const Expr* e, const Expr* xvar) {
    if (!e) return NULL;
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Slot
        && e->data.function.arg_count == 1
        && e->data.function.args[0]->type == EXPR_INTEGER
        && e->data.function.args[0]->data.integer == 1) {
        return expr_copy((Expr*)xvar);
    }
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t n = e->data.function.arg_count;
    Expr** args = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) {
        args[i] = slot1_to_var(e->data.function.args[i], xvar);
    }
    Expr* head = slot1_to_var(e->data.function.head, xvar);
    Expr* out = expr_new_function(head, args, n);
    free(args);
    return out;
}

/* Deep-copy `e`, replacing every symbol whose interned name pointer
 * equals `from_name` with a deep copy of `to`. */
static Expr* subst_sym(const Expr* e, const char* from_name, const Expr* to) {
    if (!e) return NULL;
    if (e->type == EXPR_SYMBOL && e->data.symbol.name == from_name) {
        return expr_copy((Expr*)to);
    }
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);
    size_t n = e->data.function.arg_count;
    Expr** args = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) {
        args[i] = subst_sym(e->data.function.args[i], from_name, to);
    }
    Expr* head = subst_sym(e->data.function.head, from_name, to);
    Expr* out = expr_new_function(head, args, n);
    free(args);
    return out;
}

/* Extract the polynomial body from Root[Function[..], k] as a fresh
 * Expr* in terms of `xvar`.  Accepts both Slot[1] (1-arg Function) and
 * bound-variable (2-arg Function) forms.  Returns NULL on malformed
 * input. */
static Expr* extract_root_poly(const Expr* root_expr, const Expr* xvar) {
    if (root_expr->type != EXPR_FUNCTION
        || !head_is(root_expr, SYM_Root)
        || root_expr->data.function.arg_count != 2) return NULL;
    Expr* fn = root_expr->data.function.args[0];
    if (fn->type != EXPR_FUNCTION || !head_is(fn, SYM_Function)) return NULL;

    Expr* raw;
    if (fn->data.function.arg_count == 1) {
        raw = slot1_to_var(fn->data.function.args[0], xvar);
    } else if (fn->data.function.arg_count == 2
               && fn->data.function.args[0]->type == EXPR_SYMBOL) {
        const char* bvar = fn->data.function.args[0]->data.symbol.name;
        raw = subst_sym(fn->data.function.args[1], bvar, xvar);
    } else {
        return NULL;
    }
    return eval_and_free(raw);
}

/* ------------------------------------------------------------------ *
 *  Closed-form radical solvers.                                       *
 *                                                                    *
 *  Each returns a fresh Expr** of length d carrying the d radical     *
 *  roots, in the natural order produced by the formula (NOT the       *
 *  canonical Root ordering — selection is done numerically by the     *
 *  caller).  Inputs are borrowed; outputs are owned by the caller.    *
 * ------------------------------------------------------------------ */

static Expr* radical_linear(Expr* a, Expr* b) {
    return eval_and_free(mk_fn2("Times", mk_neg(expr_copy(b)),
                                          mk_inv(expr_copy(a))));
}

/* a x^2 + b x + c == 0 */
static Expr** radical_quadratic(Expr* a, Expr* b, Expr* c, size_t* out_n) {
    Expr* b_sq    = eval_and_free(mk_pow(expr_copy(b), mk_int(2)));
    Expr* four_ac = eval_and_free(mk_fn3("Times", mk_int(4),
                                                  expr_copy(a),
                                                  expr_copy(c)));
    Expr* D       = eval_and_free(mk_fn2("Plus", b_sq, mk_neg(four_ac)));
    Expr* sqrtD   = eval_and_free(mk_sqrt(D));
    Expr* inv_2a  = eval_and_free(mk_inv(mk_fn2("Times", mk_int(2),
                                                          expr_copy(a))));
    Expr* neg_b   = mk_neg(expr_copy(b));

    Expr* r_minus = eval_and_free(mk_fn2("Times",
        mk_fn2("Plus", expr_copy(neg_b), mk_neg(expr_copy(sqrtD))),
        expr_copy(inv_2a)));
    Expr* r_plus  = eval_and_free(mk_fn2("Times",
        mk_fn2("Plus", neg_b, sqrtD), inv_2a));

    Expr** out = (Expr**)malloc(sizeof(Expr*) * 2);
    out[0] = r_minus;
    out[1] = r_plus;
    *out_n = 2;
    return out;
}

/* Cardano: a x^3 + b x^2 + c x + d == 0.  Mirrors the formula used by
 * solvepoly's solve_cubic_radical, kept local so radicals.c does not
 * depend on solvepoly's private internals. */
static Expr** radical_cubic(Expr* a, Expr* b, Expr* c, Expr* d, size_t* out_n) {
    /* Δ₀ = b² − 3 a c */
    Expr* D0 = eval_and_free(mk_fn2("Plus",
        eval_and_free(mk_pow(expr_copy(b), mk_int(2))),
        mk_neg(eval_and_free(mk_fn3("Times", mk_int(3),
                                              expr_copy(a),
                                              expr_copy(c))))));

    /* Δ₁ = 2 b³ − 9 a b c + 27 a² d */
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

    /* C_arg = (Δ₁ + sqrt_inner) / 2; switch to the other branch when
     * the first is structurally zero (avoids the triple-root degeneracy
     * collapsing the cube root). */
    Expr* C_arg = eval_and_free(mk_fn2("Times",
        mk_fn2("Plus", expr_copy(D1), expr_copy(sqrt_inner)),
        mk_rat(1, 2)));
    bool triple_root = false;
    if (is_int_zero(C_arg)) {
        Expr* alt = eval_and_free(mk_fn2("Times",
            mk_fn2("Plus", expr_copy(D1), mk_neg(expr_copy(sqrt_inner))),
            mk_rat(1, 2)));
        expr_free(C_arg);
        C_arg = alt;
        if (is_int_zero(C_arg)) triple_root = true;
    }
    Expr* C_val = triple_root
        ? mk_int(0)
        : eval_and_free(mk_pow(expr_copy(C_arg), mk_rat(1, 3)));
    expr_free(C_arg);

    /* ζ = Exp[2 π I / 3] */
    Expr* zeta = eval_and_free(mk_fn1("Exp",
        mk_fn3("Times", mk_rat(2, 3), mk_sym("Pi"), mk_sym("I"))));

    /* inv_3a = 1 / (3 a) */
    Expr* inv_3a = eval_and_free(mk_inv(
        mk_fn2("Times", mk_int(3), expr_copy(a))));

    Expr** out = (Expr**)malloc(sizeof(Expr*) * 3);
    for (int k = 0; k < 3; k++) {
        Expr* zk = (k == 0)
            ? mk_int(1)
            : eval_and_free(mk_pow(expr_copy(zeta), mk_int(k)));
        Expr* zkC, *D0_term;
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
        out[k] = eval_and_free(mk_fn3("Times", mk_int(-1),
                                                expr_copy(inv_3a),
                                                sum_in));
    }
    expr_free(zeta);
    expr_free(C_val);
    expr_free(D0);
    expr_free(D1);
    expr_free(sqrt_inner);
    expr_free(inv_3a);

    *out_n = 3;
    return out;
}

/* Ferrari: a x^4 + b x^3 + c x^2 + d x + e == 0.
 *
 * Substitute x = y - b/(4a) to get the depressed quartic
 *     y^4 + p y^2 + q y + r == 0
 * with
 *     p = (8 a c − 3 b²) / (8 a²)
 *     q = (b³ − 4 a b c + 8 a² d) / (8 a³)
 *     r = (−3 b⁴ + 256 a³ e − 64 a² b d + 16 a b² c) / (256 a⁴)
 *
 * When q == 0 the depressed quartic is biquadratic — solve as a
 * quadratic in y² and take the four square roots directly.
 *
 * Otherwise let t be any root of the resolvent cubic
 *     t³ + 2 p t² + (p² − 4 r) t − q² == 0
 * (we pick t = first Cardano root).  Then
 *     y² − √t · y + (p/2 + t/2 + q/(2 √t)) == 0
 *     y² + √t · y + (p/2 + t/2 − q/(2 √t)) == 0
 * each contribute two roots y_i; the four x_i = y_i − b/(4a) follow. */
static Expr** radical_quartic(Expr* a, Expr* b, Expr* c, Expr* d, Expr* e,
                              size_t* out_n) {
    /* p = (8 a c − 3 b²) / (8 a²) */
    Expr* eight_ac = eval_and_free(mk_fn3("Times", mk_int(8),
                                                    expr_copy(a),
                                                    expr_copy(c)));
    Expr* three_b2 = eval_and_free(mk_fn2("Times", mk_int(3),
        eval_and_free(mk_pow(expr_copy(b), mk_int(2)))));
    Expr* num_p = eval_and_free(mk_fn2("Plus", eight_ac, mk_neg(three_b2)));
    Expr* den_p = eval_and_free(mk_fn2("Times", mk_int(8),
        eval_and_free(mk_pow(expr_copy(a), mk_int(2)))));
    Expr* p = eval_and_free(mk_fn2("Times", num_p, mk_inv(den_p)));

    /* q = (b³ − 4 a b c + 8 a² d) / (8 a³) */
    Expr* b3 = eval_and_free(mk_pow(expr_copy(b), mk_int(3)));
    Expr* four_abc = eval_and_free(mk_fn4("Times", mk_int(4),
                                                    expr_copy(a),
                                                    expr_copy(b),
                                                    expr_copy(c)));
    Expr* eight_a2d = eval_and_free(mk_fn3("Times", mk_int(8),
        eval_and_free(mk_pow(expr_copy(a), mk_int(2))), expr_copy(d)));
    Expr* num_q = eval_and_free(mk_fn3("Plus", b3, mk_neg(four_abc), eight_a2d));
    Expr* den_q = eval_and_free(mk_fn2("Times", mk_int(8),
        eval_and_free(mk_pow(expr_copy(a), mk_int(3)))));
    Expr* q = eval_and_free(mk_fn2("Times", num_q, mk_inv(den_q)));

    /* r = (−3 b⁴ + 256 a³ e − 64 a² b d + 16 a b² c) / (256 a⁴) */
    Expr* b4 = eval_and_free(mk_pow(expr_copy(b), mk_int(4)));
    Expr* neg_three_b4 = eval_and_free(mk_fn2("Times", mk_int(-3), b4));
    Expr* a3 = eval_and_free(mk_pow(expr_copy(a), mk_int(3)));
    Expr* term_e = eval_and_free(mk_fn3("Times", mk_int(256), a3, expr_copy(e)));
    Expr* a2 = eval_and_free(mk_pow(expr_copy(a), mk_int(2)));
    Expr* term_d = eval_and_free(mk_fn4("Times", mk_int(-64),
                                                  a2,
                                                  expr_copy(b),
                                                  expr_copy(d)));
    Expr* b2 = eval_and_free(mk_pow(expr_copy(b), mk_int(2)));
    Expr* term_c = eval_and_free(mk_fn4("Times", mk_int(16),
                                                  expr_copy(a),
                                                  b2,
                                                  expr_copy(c)));
    Expr* num_r = eval_and_free(mk_fn4("Plus", neg_three_b4, term_e,
                                                term_d, term_c));
    Expr* a4 = eval_and_free(mk_pow(expr_copy(a), mk_int(4)));
    Expr* den_r = eval_and_free(mk_fn2("Times", mk_int(256), a4));
    Expr* r = eval_and_free(mk_fn2("Times", num_r, mk_inv(den_r)));

    /* shift = -b/(4a) (the constant added to y to recover x). */
    Expr* shift = eval_and_free(mk_neg(mk_fn2("Times", expr_copy(b),
        mk_inv(mk_fn2("Times", mk_int(4), expr_copy(a))))));

    Expr** out = (Expr**)malloc(sizeof(Expr*) * 4);

    if (is_int_zero(q)) {
        /* Biquadratic branch: y² satisfies y² = (-p ± √(p² − 4r)) / 2 */
        Expr* p_sq    = eval_and_free(mk_pow(expr_copy(p), mk_int(2)));
        Expr* four_r  = eval_and_free(mk_fn2("Times", mk_int(4), expr_copy(r)));
        Expr* disc    = eval_and_free(mk_fn2("Plus", p_sq, mk_neg(four_r)));
        Expr* sqd     = eval_and_free(mk_sqrt(disc));
        Expr* u_minus = eval_and_free(mk_fn2("Times",
            mk_fn2("Plus", mk_neg(expr_copy(p)), mk_neg(expr_copy(sqd))),
            mk_rat(1, 2)));
        Expr* u_plus  = eval_and_free(mk_fn2("Times",
            mk_fn2("Plus", mk_neg(expr_copy(p)), sqd), mk_rat(1, 2)));
        Expr* y_a = eval_and_free(mk_sqrt(u_minus));
        Expr* y_b = eval_and_free(mk_neg(expr_copy(y_a)));
        Expr* y_c = eval_and_free(mk_sqrt(u_plus));
        Expr* y_d = eval_and_free(mk_neg(expr_copy(y_c)));
        out[0] = eval_and_free(mk_fn2("Plus", expr_copy(shift), y_a));
        out[1] = eval_and_free(mk_fn2("Plus", expr_copy(shift), y_b));
        out[2] = eval_and_free(mk_fn2("Plus", expr_copy(shift), y_c));
        out[3] = eval_and_free(mk_fn2("Plus", expr_copy(shift), y_d));
    } else {
        /* Resolvent cubic in t: t³ + 2 p t² + (p² − 4 r) t − q² = 0. */
        Expr* p2_e   = eval_and_free(mk_pow(expr_copy(p), mk_int(2)));
        Expr* fourR  = eval_and_free(mk_fn2("Times", mk_int(4), expr_copy(r)));
        Expr* Cc     = eval_and_free(mk_fn2("Plus", p2_e, mk_neg(fourR)));
        Expr* q_sq   = eval_and_free(mk_pow(expr_copy(q), mk_int(2)));
        Expr* Dc     = mk_neg(q_sq);
        Expr* A_c    = mk_int(1);
        Expr* B_c    = eval_and_free(mk_fn2("Times", mk_int(2), expr_copy(p)));
        size_t cn = 0;
        Expr** trs = radical_cubic(A_c, B_c, Cc, Dc, &cn);
        expr_free(A_c); expr_free(B_c); expr_free(Cc); expr_free(Dc);

        /* Pick the natural-first Cardano branch. */
        Expr* t = expr_copy(trs[0]);
        for (size_t i = 0; i < cn; i++) expr_free(trs[i]);
        free(trs);

        Expr* sqrt_t          = eval_and_free(mk_sqrt(expr_copy(t)));
        Expr* inv_2_sqrt_t    = eval_and_free(mk_inv(
            mk_fn2("Times", mk_int(2), expr_copy(sqrt_t))));
        Expr* q_over_2sqrt_t  = eval_and_free(mk_fn2("Times",
            expr_copy(q), expr_copy(inv_2_sqrt_t)));
        Expr* half_p_t        = eval_and_free(mk_fn2("Plus",
            mk_fn2("Times", mk_rat(1, 2), expr_copy(p)),
            mk_fn2("Times", mk_rat(1, 2), expr_copy(t))));

        /* Branch a: y² − √t y + (p/2 + t/2 + q/(2 √t)) == 0.
         * radical_quadratic borrows its a/b/c arguments, so each must be
         * held in a local for the trailing expr_free pass. */
        Expr* a_a = mk_int(1);
        Expr* b_a = eval_and_free(mk_neg(expr_copy(sqrt_t)));
        Expr* c_a = eval_and_free(mk_fn2("Plus",
            expr_copy(half_p_t), expr_copy(q_over_2sqrt_t)));
        size_t qn = 0;
        Expr** ys_a = radical_quadratic(a_a, b_a, c_a, &qn);
        expr_free(a_a); expr_free(b_a); expr_free(c_a);

        /* Branch b: y² + √t y + (p/2 + t/2 − q/(2 √t)) == 0 */
        Expr* a_b = mk_int(1);
        Expr* b_b = expr_copy(sqrt_t);
        Expr* c_b = eval_and_free(mk_fn2("Plus",
            expr_copy(half_p_t), mk_neg(expr_copy(q_over_2sqrt_t))));
        Expr** ys_b = radical_quadratic(a_b, b_b, c_b, &qn);
        expr_free(a_b); expr_free(b_b); expr_free(c_b);

        expr_free(half_p_t);
        expr_free(q_over_2sqrt_t);
        expr_free(inv_2_sqrt_t);
        expr_free(sqrt_t);
        expr_free(t);

        out[0] = eval_and_free(mk_fn2("Plus", expr_copy(shift), ys_a[0]));
        out[1] = eval_and_free(mk_fn2("Plus", expr_copy(shift), ys_a[1]));
        out[2] = eval_and_free(mk_fn2("Plus", expr_copy(shift), ys_b[0]));
        out[3] = eval_and_free(mk_fn2("Plus", expr_copy(shift), ys_b[1]));
        free(ys_a);
        free(ys_b);
    }

    expr_free(shift);
    expr_free(p);
    expr_free(q);
    expr_free(r);
    *out_n = 4;
    return out;
}

/* a x^n + b == 0  ⇒  x = (-b/a)^(1/n) · (-1)^(2k/n) for k = 0 .. n-1.
 * Routing through Power[-1, 2k/n] (rather than Exp[2 π I k / n]) lets
 * Power's canonical simplification produce clean radical expressions
 * like `(-1)^(4/5) 2^(1/5)`. */
static Expr** radical_binomial(Expr* a, Expr* b, int64_t n, size_t* out_n) {
    Expr* base = eval_and_free(mk_fn2("Times", mk_neg(expr_copy(b)),
                                               mk_inv(expr_copy(a))));
    Expr* r = eval_and_free(mk_pow(expr_copy(base), mk_rat(1, n)));
    expr_free(base);
    Expr** out = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    for (int64_t k = 0; k < n; k++) {
        if (k == 0) {
            out[0] = expr_copy(r);
            continue;
        }
        Expr* exp_arg = (2 * k) % n == 0
            ? mk_int((2 * k) / n)
            : mk_rat(2 * k, n);
        Expr* expk = eval_and_free(mk_pow(mk_int(-1), exp_arg));
        out[k] = eval_and_free(mk_fn2("Times", expr_copy(r), expk));
    }
    expr_free(r);
    *out_n = (size_t)n;
    return out;
}

/* True iff `poly` is a x^d + b in `var`, i.e. only the top and constant
 * coefficients are nonzero.  Outputs are freshly owned on success. */
static bool is_binomial_shape(Expr* poly, Expr* var, int d,
                              Expr** a_out, Expr** b_out, int64_t* n_out) {
    if (d < 2) return false;
    Expr* aa = get_coeff(poly, var, d);
    if (!aa || is_int_zero(aa)) { expr_free(aa); return false; }
    Expr* bb = get_coeff(poly, var, 0);
    if (!bb) { expr_free(aa); return false; }
    for (int k = 1; k < d; k++) {
        Expr* ck = get_coeff(poly, var, k);
        bool zero = is_int_zero(ck);
        expr_free(ck);
        if (!zero) { expr_free(aa); expr_free(bb); return false; }
    }
    *a_out = aa;
    *b_out = bb;
    *n_out = d;
    return true;
}

/* ------------------------------------------------------------------ *
 *  k-th root selection.                                               *
 *                                                                    *
 *  Each formula above emits its d roots in the order natural to the   *
 *  derivation (Cardano: x_0, x_1 = ζ-rotated, x_2 = ζ²-rotated;        *
 *  Ferrari: two y-branches × two quadratic roots; binomial: k =       *
 *  0 .. n-1).  None of these match Mathilda's canonical Root[..., k]  *
 *  ordering (real-first ascending, then complex by Re asc / |Im| asc /
 *  negative-Im first).  We pick the k-th radical by numerically       *
 *  evaluating Root[poly, k] via the existing root_numericalize        *
 *  pipeline (machine precision is enough to distinguish well-          *
 *  separated roots) and choosing the radical root with the smallest   *
 *  Euclidean distance in the complex plane.                            *
 * ------------------------------------------------------------------ */

/* Extract (re, im) doubles from a numericalised Expr.  Accepts Integer,
 * Real, Rational (after one evaluate), and Complex[Integer|Real, ...].
 * Returns false for symbolic residues. */
static bool to_complex_double(const Expr* e, double* re_out, double* im_out) {
    if (!e) return false;
    double re = 0.0, im = 0.0;
    if (e->type == EXPR_INTEGER) {
        re = (double)e->data.integer;
    } else if (e->type == EXPR_REAL) {
        re = e->data.real;
    } else if (e->type == EXPR_FUNCTION
               && head_is(e, SYM_Complex)
               && e->data.function.arg_count == 2) {
        Expr* re_e = e->data.function.args[0];
        Expr* im_e = e->data.function.args[1];
        if      (re_e->type == EXPR_INTEGER) re = (double)re_e->data.integer;
        else if (re_e->type == EXPR_REAL)    re = re_e->data.real;
        else return false;
        if      (im_e->type == EXPR_INTEGER) im = (double)im_e->data.integer;
        else if (im_e->type == EXPR_REAL)    im = im_e->data.real;
        else return false;
    } else {
        return false;
    }
    *re_out = re;
    *im_out = im;
    return true;
}

/* N[e] at machine precision; project onto (re, im) doubles. */
static bool numericalize_to_complex(const Expr* e, double* re_out, double* im_out) {
    NumericSpec spec = numeric_machine_spec();
    Expr* num = numericalize(e, spec);
    if (!num) return false;
    Expr* evald = evaluate(num);
    expr_free(num);
    bool ok = to_complex_double(evald, re_out, im_out);
    expr_free(evald);
    return ok;
}

/* Pick the index of the radical root in roots[0..n-1] closest to the
 * complex target (target_re, target_im).  Returns -1 if no radical
 * could be numerically evaluated. */
static int pick_closest_root(Expr** roots, size_t n,
                             double target_re, double target_im) {
    int best = -1;
    double best_dist = 0.0;
    for (size_t i = 0; i < n; i++) {
        double re = 0.0, im = 0.0;
        if (!numericalize_to_complex(roots[i], &re, &im)) continue;
        double dre = re - target_re;
        double dim = im - target_im;
        double dist = dre * dre + dim * dim;
        if (best < 0 || dist < best_dist) {
            best_dist = dist;
            best = (int)i;
        }
    }
    return best;
}

/* ------------------------------------------------------------------ *
 *  Per-Root conversion.                                               *
 * ------------------------------------------------------------------ */

/* Solve poly[xvar] == 0 for radical roots; on success returns a freshly
 * owned Expr** of length *nroots_out and a positive int degree.  On
 * unsupported shapes (degree >= 5 non-binomial) returns NULL.  Caller
 * still owns `poly`, `xvar`, and the returned root array on success. */
static Expr** solve_radical_roots(Expr* poly, Expr* xvar, int d,
                                  size_t* nroots_out) {
    if (d == 1) {
        Expr* a = get_coeff(poly, xvar, 1);
        Expr* b = get_coeff(poly, xvar, 0);
        if (!a || !b) { expr_free(a); expr_free(b); return NULL; }
        Expr** out = (Expr**)malloc(sizeof(Expr*));
        out[0] = radical_linear(a, b);
        *nroots_out = 1;
        expr_free(a); expr_free(b);
        return out;
    }
    if (d == 2) {
        Expr* a = get_coeff(poly, xvar, 2);
        Expr* b = get_coeff(poly, xvar, 1);
        Expr* c = get_coeff(poly, xvar, 0);
        if (!a || !b || !c) {
            expr_free(a); expr_free(b); expr_free(c);
            return NULL;
        }
        Expr** out = radical_quadratic(a, b, c, nroots_out);
        expr_free(a); expr_free(b); expr_free(c);
        return out;
    }
    if (d == 3) {
        Expr* a = get_coeff(poly, xvar, 3);
        Expr* b = get_coeff(poly, xvar, 2);
        Expr* c = get_coeff(poly, xvar, 1);
        Expr* dd = get_coeff(poly, xvar, 0);
        if (!a || !b || !c || !dd) {
            expr_free(a); expr_free(b); expr_free(c); expr_free(dd);
            return NULL;
        }
        Expr** out = radical_cubic(a, b, c, dd, nroots_out);
        expr_free(a); expr_free(b); expr_free(c); expr_free(dd);
        return out;
    }
    if (d == 4) {
        Expr* a = get_coeff(poly, xvar, 4);
        Expr* b = get_coeff(poly, xvar, 3);
        Expr* c = get_coeff(poly, xvar, 2);
        Expr* dd = get_coeff(poly, xvar, 1);
        Expr* ee = get_coeff(poly, xvar, 0);
        if (!a || !b || !c || !dd || !ee) {
            expr_free(a); expr_free(b); expr_free(c);
            expr_free(dd); expr_free(ee);
            return NULL;
        }
        Expr** out = radical_quartic(a, b, c, dd, ee, nroots_out);
        expr_free(a); expr_free(b); expr_free(c); expr_free(dd); expr_free(ee);
        return out;
    }
    /* d >= 5: binomial fast path only. */
    {
        Expr *aa, *bb;
        int64_t nn;
        if (is_binomial_shape(poly, xvar, d, &aa, &bb, &nn)) {
            Expr** out = radical_binomial(aa, bb, nn, nroots_out);
            expr_free(aa); expr_free(bb);
            return out;
        }
    }
    return NULL;
}

/* Try to convert one Root[Function[..], k] node into its radical form.
 * Returns NULL when the input is not a convertible Root (caller must
 * leave the node unchanged). */
static Expr* convert_root_to_radical(const Expr* root_expr) {
    if (root_expr->type != EXPR_FUNCTION
        || !head_is(root_expr, SYM_Root)
        || root_expr->data.function.arg_count != 2) return NULL;

    Expr* k_expr = root_expr->data.function.args[1];
    if (k_expr->type != EXPR_INTEGER || k_expr->data.integer < 1) return NULL;
    int64_t k = k_expr->data.integer;

    /* Fresh dummy variable in a private context so it cannot clash with
     * any user-bound symbol. */
    Expr* xvar = expr_new_symbol("ToRadicals`Private`x$");
    Expr* poly = extract_root_poly(root_expr, xvar);
    if (!poly) { expr_free(xvar); return NULL; }

    Expr* vars[1] = { xvar };
    if (!is_polynomial(poly, vars, 1)) {
        expr_free(poly); expr_free(xvar);
        return NULL;
    }
    int d = get_degree_poly(poly, xvar);
    if (d < 1 || (int64_t)d < k) {
        expr_free(poly); expr_free(xvar);
        return NULL;
    }

    size_t nroots = 0;
    Expr** roots = solve_radical_roots(poly, xvar, d, &nroots);
    expr_free(poly);
    expr_free(xvar);
    if (!roots || nroots == 0) {
        free(roots);
        return NULL;
    }

    /* Single root (d == 1): k must already be 1. */
    if (nroots == 1) {
        Expr* sel = roots[0];
        free(roots);
        return sel;
    }

    /* Numeric matching against the canonical Root[..., k] value. */
    int idx = -1;
    NumericSpec spec = numeric_machine_spec();
    Expr* target = root_numericalize(root_expr, spec);
    if (target) {
        double tre = 0.0, tim = 0.0;
        if (to_complex_double(target, &tre, &tim)) {
            idx = pick_closest_root(roots, nroots, tre, tim);
        }
        expr_free(target);
    }
    if (idx < 0) {
        /* Parametric / symbolic input — root_numericalize cannot
         * decide.  Fall back to natural per-formula index. */
        if ((size_t)(k - 1) < nroots) idx = (int)(k - 1);
    }

    if (idx < 0 || (size_t)idx >= nroots) {
        for (size_t i = 0; i < nroots; i++) expr_free(roots[i]);
        free(roots);
        return NULL;
    }
    Expr* sel = roots[idx];
    for (size_t i = 0; i < nroots; i++) {
        if ((int)i != idx) expr_free(roots[i]);
    }
    free(roots);
    return sel;
}

/* ------------------------------------------------------------------ *
 *  Tree walker.                                                       *
 * ------------------------------------------------------------------ */

/* Recursively reconstruct `e`, replacing every Root[Function[..], k]
 * node with its radical form when one can be produced.  Atoms and Root
 * nodes that fail to convert are deep-copied verbatim.  This is what
 * gives ToRadicals its threading semantics: a Root inside List, Equal,
 * Less, And, Or, Rule, ... is processed identically to one at the top
 * level. */
static Expr* radicals_walk(const Expr* e) {
    if (!e) return NULL;
    if (e->type != EXPR_FUNCTION) return expr_copy((Expr*)e);

    if (head_is(e, SYM_Root)) {
        Expr* converted = convert_root_to_radical(e);
        if (converted) return converted;
        return expr_copy((Expr*)e);
    }

    size_t n = e->data.function.arg_count;
    Expr* head = radicals_walk(e->data.function.head);
    Expr** args = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++) {
        args[i] = radicals_walk(e->data.function.args[i]);
    }
    Expr* out = expr_new_function(head, args, n);
    free(args);
    return out;
}

Expr* builtin_to_radicals(Expr* res) {
    if (res->type != EXPR_FUNCTION
        || res->data.function.arg_count != 1) return NULL;
    Expr* walked = radicals_walk(res->data.function.args[0]);
    return eval_and_free(walked);
}

void radicals_init(void) {
    symtab_add_builtin("ToRadicals", builtin_to_radicals);
    symtab_get_def("ToRadicals")->attributes |= ATTR_PROTECTED;
}
