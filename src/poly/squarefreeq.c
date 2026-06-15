/*
 * squarefreeq.c -- SquareFreeQ[expr] / SquareFreeQ[expr, vars] / opts.
 *
 * Always returns True or False on a structurally valid call.
 * Wrong arg count emits `SquareFreeQ::argb` to stderr and returns NULL.
 * Malformed options emit `SquareFreeQ::nonopt` and return NULL.
 *
 * Algorithms:
 *   - Integer n: factor with FactorInteger, check every prime exponent <= 1.
 *               Special-case 0 -> False, +/-1 -> True.
 *   - Rational p/q: square-free iff both numerator and denominator are.
 *   - Gaussian integer a + b I (with GaussianIntegers -> True or auto-detected
 *     from a Complex[Integer, Integer] input): factor N(z) = a^2 + b^2 over Z
 *     and dispatch by the rational prime's residue mod 4 (see sqfree_gaussian).
 *   - Polynomial in `vars`: for every var x_i that the polynomial has positive
 *     degree in, compute PolynomialGCD(p, dp/dx_i); the polynomial is
 *     square-free iff every such gcd is independent of x_i (degree 0 in x_i).
 *   - Anything else (Real, symbolic): False -- per Mathematica's "expr is not
 *     manifestly square free" semantics.
 *
 * The Modulus option is parsed but only Modulus -> 0 is honoured; non-zero
 * values emit `SquareFreeQ::modnotimpl` and return the call unevaluated until
 * a real mod-p polynomial sqfree test is wired in.
 */

#include "squarefreeq.h"
#include "facpoly.h"
#include "poly.h"
#include "eval.h"
#include "expand.h"
#include "symtab.h"
#include "attr.h"
#include "arithmetic.h"
#include "internal.h"
#include "sym_names.h"
#include "print.h"
#include "expr.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <gmp.h>

/* ===================================================================== */
/* Forward declarations                                                  */
/* ===================================================================== */

static bool sqfree_dispatch(Expr* e, Expr** vars, size_t v_count,
                            bool gaussian);
static bool sqfree_integer_mpz(const mpz_t n);
static bool sqfree_rational(int64_t n, int64_t d);
static bool sqfree_rational_mpz(const mpz_t n, const mpz_t d);
static bool sqfree_gaussian(const mpz_t a, const mpz_t b);
static bool sqfree_polynomial(Expr* p, Expr** vars, size_t v_count);
static Expr* squarefreeq_deriv(Expr* p, Expr* var);

/* ===================================================================== */
/* Small utilities                                                       */
/* ===================================================================== */

static bool is_sym_eq(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol, name) == 0;
}

static bool is_rule_head(const Expr* e) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           (e->data.function.head->data.symbol == SYM_Rule ||
            e->data.function.head->data.symbol == SYM_RuleDelayed) &&
           e->data.function.arg_count == 2;
}

/* True iff `e` is a "variable-like" leaf usable as a polynomial generator:
 * a bare Symbol, OR a List of bare Symbols.  Anything else (numbers,
 * function calls, etc.) is rejected by the vars parser. */
static bool is_var_spec(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return true;
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_List) {
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            if (e->data.function.args[i]->type != EXPR_SYMBOL) return false;
        }
        return true;
    }
    return false;
}

/* True if `e` is Complex[a, b] with both parts integer-like and b != 0.
 * Used by the GaussianIntegers -> Automatic mode to flip into Z[i] for
 * Gaussian-integer inputs without the user spelling out the option. */
static bool is_complex_integer(const Expr* e, mpz_t a_out, mpz_t b_out) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    if (e->data.function.head->data.symbol != SYM_Complex) return false;
    if (e->data.function.arg_count != 2) return false;
    Expr* re = e->data.function.args[0];
    Expr* im = e->data.function.args[1];
    if (!(re->type == EXPR_INTEGER || re->type == EXPR_BIGINT)) return false;
    if (!(im->type == EXPR_INTEGER || im->type == EXPR_BIGINT)) return false;
    expr_to_mpz(re, a_out);
    expr_to_mpz(im, b_out);
    return true;
}

/* ===================================================================== */
/* Diagnostics                                                           */
/* ===================================================================== */

static Expr* sqfree_emit_argb(size_t argc) {
    fprintf(stderr,
            "SquareFreeQ::argb: SquareFreeQ called with %zu argument%s; "
            "between 1 and 3 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

static Expr* sqfree_emit_nonopt(Expr* bad, size_t pos, Expr* res) {
    char* bad_str = expr_to_string(bad);
    char* call_str = expr_to_string(res);
    fprintf(stderr,
            "SquareFreeQ::nonopt: Options expected (instead of %s) beyond "
            "position %zu in %s. An option must be a rule or a list of "
            "rules.\n",
            bad_str ? bad_str : "?", pos, call_str ? call_str : "?");
    free(bad_str);
    free(call_str);
    return NULL;
}

/* `Modulus -> n` with n other than 0 is parsed but not yet implemented.
 * Emit a Mathematica-style diagnostic and signal "leave unevaluated" by
 * returning NULL so the surface call stays visible to the user. */
static Expr* sqfree_emit_modnotimpl(Expr* val) {
    char* val_str = expr_to_string(val);
    fprintf(stderr,
            "SquareFreeQ::modnotimpl: Modulus -> %s is not yet supported; "
            "only Modulus -> 0 (the default integer ring) is currently "
            "implemented.\n",
            val_str ? val_str : "?");
    free(val_str);
    return NULL;
}

/* ===================================================================== */
/* Integer / rational square-free tests                                  */
/* ===================================================================== */

/* True iff |n| has no rational prime factor of multiplicity >= 2.
 * Handles 0 (False) and units +/-1 (True) without calling FactorInteger. */
static bool sqfree_integer_mpz(const mpz_t n) {
    if (mpz_sgn(n) == 0) return false;
    mpz_t a;
    mpz_init_set(a, n);
    if (mpz_sgn(a) < 0) mpz_neg(a, a);
    if (mpz_cmp_ui(a, 1) == 0) { mpz_clear(a); return true; }

    /* Call FactorInteger via the public internal entry; it returns
     * {{p, e}, ...} for the rational-prime factorisation of |n|. */
    Expr* n_expr = expr_bigint_normalize(expr_new_bigint_from_mpz(a));
    mpz_clear(a);
    Expr** args = malloc(sizeof(Expr*) * 1);
    args[0] = n_expr;
    Expr* fi = internal_factorinteger(args, 1);
    free(args);
    if (!fi) return false;

    bool ok = true;
    if (fi->type == EXPR_FUNCTION &&
        fi->data.function.head->type == EXPR_SYMBOL &&
        fi->data.function.head->data.symbol == SYM_List) {
        for (size_t i = 0; i < fi->data.function.arg_count; i++) {
            Expr* pair = fi->data.function.args[i];
            if (pair->type != EXPR_FUNCTION ||
                pair->data.function.arg_count != 2) { ok = false; break; }
            Expr* exp_e = pair->data.function.args[1];
            int64_t exp_val = 0;
            if (exp_e->type == EXPR_INTEGER) exp_val = exp_e->data.integer;
            else if (exp_e->type == EXPR_BIGINT) {
                /* No 64-bit prime can have a fits-in-int exponent >= 2^63;
                 * treat unrepresentably-large exponents as "not square-free". */
                ok = false; break;
            }
            if (exp_val >= 2) { ok = false; break; }
        }
    } else {
        ok = false;
    }
    expr_free(fi);
    return ok;
}

static bool sqfree_rational(int64_t n, int64_t d) {
    mpz_t mn, md;
    mpz_init_set_si(mn, n);
    mpz_init_set_si(md, d);
    bool r = sqfree_rational_mpz(mn, md);
    mpz_clear(mn);
    mpz_clear(md);
    return r;
}

static bool sqfree_rational_mpz(const mpz_t n, const mpz_t d) {
    /* p/q is square-free iff both p and q are square-free integers
     * (the implicit Rational normaliser already guarantees gcd(p,q)=1). */
    return sqfree_integer_mpz(n) && sqfree_integer_mpz(d);
}

/* ===================================================================== */
/* Gaussian integer square-free test                                     */
/* ===================================================================== */

/* Cornacchia: given an odd rational prime p with p ≡ 1 (mod 4), find a
 * representation p = c^2 + d^2 with c, d >= 0.  Returns true on success.
 *
 * Method: brute search c in [1, floor(sqrt(p))], test whether p - c^2 is
 * a perfect square.  Fine for any prime that fits in a few thousand bits
 * because the loop runs O(sqrt(p)) operations, all O(1) mpz_sqrt /
 * mpz_perfect_square_p calls.  For modest-size primes this is well under
 * a millisecond; for very large primes Cornacchia proper (via Tonelli)
 * would scale better but is overkill for the SquareFreeQ use case. */
static bool cornacchia_p1_mod4(const mpz_t p, mpz_t c_out, mpz_t d_out) {
    mpz_t c, c2, rem, root;
    mpz_inits(c, c2, rem, root, NULL);
    mpz_sqrt(root, p);  /* root = floor(sqrt(p)) */
    bool found = false;
    for (mpz_set_ui(c, 1); mpz_cmp(c, root) <= 0; mpz_add_ui(c, c, 1)) {
        mpz_mul(c2, c, c);
        mpz_sub(rem, p, c2);
        if (mpz_perfect_square_p(rem)) {
            mpz_sqrt(d_out, rem);
            mpz_set(c_out, c);
            found = true;
            break;
        }
    }
    mpz_clears(c, c2, rem, root, NULL);
    return found;
}

/* Divide z = a + b I by w = c + d I over Z[i] *only when the result is
 * itself a Gaussian integer*.  On exact divisibility writes the quotient
 * to (a, b) in-place and returns true; otherwise leaves (a, b) untouched
 * and returns false.
 *
 * (a + bi)/(c + di) = ((ac + bd) + (bc - ad) i) / (c^2 + d^2),
 * exact iff (c^2 + d^2) divides both ac + bd and bc - ad. */
static bool gaussian_divide(mpz_t a, mpz_t b, const mpz_t c, const mpz_t d) {
    mpz_t norm, ac, bd, bc, ad, num_re, num_im, q_re, q_im, rem;
    mpz_inits(norm, ac, bd, bc, ad, num_re, num_im, q_re, q_im, rem, NULL);

    mpz_mul(ac, a, c);
    mpz_mul(bd, b, d);
    mpz_mul(bc, b, c);
    mpz_mul(ad, a, d);
    mpz_mul(norm, c, c);
    mpz_addmul(norm, d, d);  /* norm = c^2 + d^2 */

    mpz_add(num_re, ac, bd);
    mpz_sub(num_im, bc, ad);

    mpz_fdiv_qr(q_re, rem, num_re, norm);
    bool ok = (mpz_sgn(rem) == 0);
    if (ok) {
        mpz_fdiv_qr(q_im, rem, num_im, norm);
        ok = (mpz_sgn(rem) == 0);
    }
    if (ok) {
        mpz_set(a, q_re);
        mpz_set(b, q_im);
    }
    mpz_clears(norm, ac, bd, bc, ad, num_re, num_im, q_re, q_im, rem, NULL);
    return ok;
}

/* Multiplicity of the Gaussian prime (c + d I) in (a + b I).  Returns
 * the largest k with (c + d I)^k | (a + b I); modifies (a, b) to the
 * fully stripped quotient.  Caller owns a, b. */
static int gaussian_strip_multiplicity(mpz_t a, mpz_t b,
                                       const mpz_t c, const mpz_t d) {
    int k = 0;
    while (gaussian_divide(a, b, c, d)) k++;
    return k;
}

/* Decide square-freeness of z = a + b I in Z[i].  Returns False for
 * z == 0; True for units (norm == 1). */
static bool sqfree_gaussian(const mpz_t a_in, const mpz_t b_in) {
    /* Norm 0 ↔ z == 0; treat zero as non-square-free for parity with
     * SquareFreeQ[0] -> False over Z. */
    mpz_t norm;
    mpz_init(norm);
    mpz_mul(norm, a_in, a_in);
    mpz_addmul(norm, b_in, b_in);
    if (mpz_sgn(norm) == 0) { mpz_clear(norm); return false; }

    /* Norm 1 ↔ z is a unit (+/-1, +/-i); units are square-free. */
    if (mpz_cmp_ui(norm, 1) == 0) { mpz_clear(norm); return true; }

    /* Factor N(z) over Z; for each rational prime p | N(z) classify by
     * its residue mod 4 to recover the Gaussian-prime structure above p. */
    Expr* norm_e = expr_bigint_normalize(expr_new_bigint_from_mpz(norm));
    mpz_clear(norm);
    Expr** ar = malloc(sizeof(Expr*) * 1);
    ar[0] = norm_e;
    Expr* fi = internal_factorinteger(ar, 1);
    free(ar);
    if (!fi || fi->type != EXPR_FUNCTION ||
        fi->data.function.head->type != EXPR_SYMBOL ||
        fi->data.function.head->data.symbol != SYM_List) {
        if (fi) expr_free(fi);
        return false;
    }

    mpz_t a, b;
    mpz_init_set(a, a_in);
    mpz_init_set(b, b_in);
    bool ok = true;

    for (size_t i = 0; ok && i < fi->data.function.arg_count; i++) {
        Expr* pair = fi->data.function.args[i];
        if (pair->type != EXPR_FUNCTION || pair->data.function.arg_count != 2) {
            ok = false; break;
        }
        Expr* p_e = pair->data.function.args[0];
        Expr* e_e = pair->data.function.args[1];
        mpz_t p;
        mpz_init(p);
        expr_to_mpz(p_e, p);
        int64_t e_in_norm = (e_e->type == EXPR_INTEGER) ? e_e->data.integer : 0;
        if (e_e->type != EXPR_INTEGER || e_in_norm <= 0) {
            mpz_clear(p);
            ok = false;
            break;
        }

        unsigned long pmod4 = mpz_fdiv_ui(p, 4);

        if (mpz_cmp_ui(p, 2) == 0) {
            /* 2 = -i (1+i)^2; the Gaussian prime above 2 is (1+i).
             * Its multiplicity in z equals the multiplicity of 2 in N(z). */
            if (e_in_norm >= 2) ok = false;
        } else if (pmod4 == 3) {
            /* p ≡ 3 (mod 4): p itself is a Gaussian prime, N(p) = p^2,
             * so each occurrence of p in z contributes p^2 to N(z). */
            if (e_in_norm >= 4) ok = false;
            else if (e_in_norm == 2) {
                /* Confirm p | z (it must, since p^2 | N(z)). The multiplicity
                 * is 1, which is square-free. */
                /* No further check needed -- e_in_norm == 2 ↔ mult(p in z) == 1. */
            }
        } else {
            /* p ≡ 1 (mod 4): p = pi * conj(pi).  Multiplicity of p in N(z) =
             * mult(pi in z) + mult(conj(pi) in z).  Square-free iff both
             * are ≤ 1; the only ambiguous case is e_in_norm == 2 (could be
             * pi^2, conj(pi)^2, or pi*conj(pi)). */
            if (e_in_norm >= 3) {
                ok = false;
            } else if (e_in_norm == 2) {
                mpz_t c, d;
                mpz_inits(c, d, NULL);
                if (!cornacchia_p1_mod4(p, c, d)) {
                    /* Should never happen for a prime p ≡ 1 (mod 4).
                     * Defensive: bail out as "not provably square-free". */
                    mpz_clears(c, d, NULL);
                    ok = false;
                } else {
                    mpz_t a2, b2;
                    mpz_init_set(a2, a);
                    mpz_init_set(b2, b);
                    int kpi = gaussian_strip_multiplicity(a2, b2, c, d);
                    if (kpi >= 2) {
                        ok = false;
                    } else {
                        /* mult(conj(pi)) = e_in_norm - kpi; must also be <= 1. */
                        int kcj = (int)e_in_norm - kpi;
                        if (kcj >= 2) ok = false;
                    }
                    mpz_clears(a2, b2, NULL);
                    mpz_clears(c, d, NULL);
                }
            }
            /* e_in_norm == 1: exactly one of pi, conj(pi) divides z, mult 1.
             * Square-free. */
        }
        mpz_clear(p);
    }

    mpz_clear(a);
    mpz_clear(b);
    expr_free(fi);
    return ok;
}

/* ===================================================================== */
/* Polynomial square-free test                                           */
/* ===================================================================== */

/* Take the partial derivative dp/dvar.  Built by routing through the
 * `D` builtin (defined in calculus/deriv.c via DownValues + builtin); we
 * don't reach into the polynomial substrate ourselves because SquareFreeQ
 * accepts arbitrary polynomial expressions and the public D[] already
 * implements the same chain/product/quotient rules. */
static Expr* squarefreeq_deriv(Expr* p, Expr* var) {
    Expr** args = malloc(sizeof(Expr*) * 2);
    args[0] = expr_copy(p);
    args[1] = expr_copy(var);
    Expr* call = expr_new_function(expr_new_symbol(SYM_D), args, 2);
    free(args);
    Expr* d = eval_and_free(call);
    return d;
}

/* Square-free criterion: p is square-free over K[vars] iff for every var
 * x_i in `vars` that p actually depends on, PolynomialGCD(p, dp/dx_i) has
 * degree 0 in x_i.  Proof: if g^2 | p with g non-constant, choose any x_i
 * appearing in g; then g | gcd(p, p_{x_i}) and the gcd has positive degree
 * in x_i.  Converse: vanishing of every such gcd-in-x_i (i.e. the gcd is
 * x_i-independent for all i) forces p and p_{x_i} to share no common
 * factor that depends on any variable, i.e. p is square-free. */
static bool sqfree_polynomial(Expr* p, Expr** vars, size_t v_count) {
    if (is_zero_poly(p)) return false;  /* 0 is not square-free */

    /* Constant in all vars: defer to numeric handling.  This is reached
     * when SquareFreeQ[4] is dispatched through the polynomial branch
     * with an empty var set; never via the top-level entry. */
    bool any_var_present = false;
    for (size_t i = 0; i < v_count; i++) {
        if (get_degree_poly(p, vars[i]) > 0) { any_var_present = true; break; }
    }
    if (!any_var_present) {
        /* Pure constant: numeric branch handles it.  Anything that reaches
         * here from the polynomial dispatcher is a polynomial-shaped expr
         * with no varying terms (e.g. user passed wrong vars list); treat
         * as square-free wrt those vars -- vacuously true. */
        return true;
    }

    for (size_t i = 0; i < v_count; i++) {
        Expr* x = vars[i];
        if (get_degree_poly(p, x) <= 0) continue;  /* p constant in x_i */

        Expr* dpdx_raw = squarefreeq_deriv(p, x);
        Expr* dpdx = expr_expand(dpdx_raw);
        expr_free(dpdx_raw);

        /* dp/dx_i can vanish identically over Z (e.g. in characteristic
         * p, x^p has zero derivative).  Over Z that means p actually has
         * no x_i-dependence -- skip, since the gcd would be p itself and
         * always positive-degree, masquerading as not-square-free. */
        if (is_zero_poly(dpdx)) { expr_free(dpdx); continue; }

        Expr* g = poly_gcd_internal(p, dpdx, vars, v_count);
        expr_free(dpdx);

        int deg = get_degree_poly(g, x);
        expr_free(g);
        if (deg > 0) return false;
    }
    return true;
}

/* ===================================================================== */
/* Top-level dispatcher                                                  */
/* ===================================================================== */

/* Drive the appropriate test for `e` given the resolved `vars` (possibly
 * empty) and `gaussian` mode (already resolved from Automatic).  Returns
 * a True/False decision; never reports "unknown" -- ambiguous inputs go
 * to False per the *Q always-bool contract. */
static bool sqfree_dispatch(Expr* e, Expr** vars, size_t v_count,
                            bool gaussian) {
    /* Plain integers (machine and big). */
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) {
        mpz_t n;
        expr_to_mpz(e, n);
        bool r;
        if (gaussian) {
            mpz_t zero;
            mpz_init_set_ui(zero, 0);
            r = sqfree_gaussian(n, zero);
            mpz_clear(zero);
        } else {
            r = sqfree_integer_mpz(n);
        }
        mpz_clear(n);
        return r;
    }

    /* Rationals: Mathilda stores them as Rational[n, d]. */
    {
        int64_t n, d;
        if (is_rational(e, &n, &d)) {
            if (gaussian) {
                /* Gaussian rationals are not integers; outside the domain. */
                return false;
            }
            return sqfree_rational(n, d);
        }
    }

    /* Complex with integer real and integer imag: dispatch to Z[i] only
     * when Gaussian mode is enabled (Automatic resolved to True, or the
     * user explicitly passed GaussianIntegers -> True).  With Gaussian
     * mode off, a Complex-integer literal is not an element of Z and
     * falls through to the final False branch -- mirroring Mathematica's
     * SquareFreeQ[3 + 2 I, GaussianIntegers -> False] -> False. */
    {
        mpz_t a, b;
        mpz_inits(a, b, NULL);
        if (gaussian && is_complex_integer(e, a, b)) {
            bool r = sqfree_gaussian(a, b);
            mpz_clears(a, b, NULL);
            return r;
        }
        mpz_clears(a, b, NULL);
        if (e->type == EXPR_FUNCTION &&
            e->data.function.head->type == EXPR_SYMBOL &&
            e->data.function.head->data.symbol == SYM_Complex) {
            /* Gaussian-mode-off Complex literal: not square-free in Z. */
            return false;
        }
    }

    /* Reals are never square-free in Mathematica's sense (not manifestly so). */
    if (e->type == EXPR_REAL) return false;

    /* Polynomials in `vars` (auto-collected when v_count == 0). */
    bool free_vars = false;
    Expr** use_vars = vars;
    size_t use_count = v_count;
    if (use_count == 0) {
        size_t cap = 16;
        use_vars = malloc(sizeof(Expr*) * cap);
        use_count = 0;
        collect_variables(e, &use_vars, &use_count, &cap);
        /* Mathematica's `Variables` only counts bare symbols; non-symbol
         * "variables" emitted by collect_variables (e.g. Sin[x], Sqrt[2])
         * signal "not a manifest polynomial".  Drop them so the auto-detect
         * path agrees with Mathematica's `SquareFreeQ[Sin[x]] -> False`. */
        size_t kept = 0;
        for (size_t i = 0; i < use_count; i++) {
            if (use_vars[i]->type == EXPR_SYMBOL) {
                use_vars[kept++] = use_vars[i];
            } else {
                expr_free(use_vars[i]);
            }
        }
        if (kept != use_count) {
            /* Owned-Expr cells past `kept` already freed above; shrink count. */
            use_count = kept;
        }
        free_vars = true;
    }

    /* No vars and not numeric: treat as a non-polynomial symbolic input.
     * Per Mathematica, return False ("not manifestly square free"). */
    if (use_count == 0) {
        if (free_vars) free(use_vars);
        return false;
    }

    /* Must be a polynomial in the chosen vars. */
    if (!is_polynomial(e, use_vars, use_count)) {
        if (free_vars) {
            for (size_t i = 0; i < use_count; i++) expr_free(use_vars[i]);
            free(use_vars);
        }
        return false;
    }

    /* Expand and test. */
    Expr* expanded = expr_expand(e);
    bool r = sqfree_polynomial(expanded, use_vars, use_count);
    expr_free(expanded);

    if (free_vars) {
        for (size_t i = 0; i < use_count; i++) expr_free(use_vars[i]);
        free(use_vars);
    }
    return r;
}

/* ===================================================================== */
/* Entry point                                                           */
/* ===================================================================== */

Expr* builtin_squarefreeq(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc == 0) return sqfree_emit_argb(0);

    Expr* expr = res->data.function.args[0];

    /* Argument layout:
     *   args[0] = expr (required, position 1)
     *   args[1] = optional vars (Symbol or List of Symbols) OR first option
     *   args[2..] = options (Rule or RuleDelayed)
     *
     * Resolution: if args[1] is a Rule, it is the first option (vars
     * defaults to auto-collected).  Otherwise args[1] is consumed as the
     * vars slot; if the value isn't a valid vars-spec we silently fall
     * back to auto-detection -- matching Mathematica's lenient slot
     * consumption, which only surfaces nonopt at the FIRST non-Rule
     * found at position >= 3.
     *
     * For the nonopt diagnostic we mirror Mathematica's phrasing:
     * "beyond position 1" with the LAST offending non-Rule reported
     * (Mathematica reports the trailing bad arg, not the first). */
    Expr** vars = NULL;
    size_t var_count = 0;
    bool free_vars = false;
    size_t opt_start = 1;

    if (argc >= 2) {
        Expr* a1 = res->data.function.args[1];
        if (is_rule_head(a1)) {
            opt_start = 1;
        } else {
            /* Permissive vars consumption: valid spec populates `vars`;
             * a bogus value (e.g. integer 2) is accepted as a no-op slot
             * and triggers fall-through to auto-detection.  Either way
             * options now start at position 2. */
            if (is_var_spec(a1)) {
                if (a1->type == EXPR_SYMBOL) {
                    vars = malloc(sizeof(Expr*) * 1);
                    vars[0] = a1;
                    var_count = 1;
                } else {
                    var_count = a1->data.function.arg_count;
                    if (var_count > 0) {
                        vars = malloc(sizeof(Expr*) * var_count);
                        for (size_t i = 0; i < var_count; i++) vars[i] = a1->data.function.args[i];
                    }
                }
                free_vars = true;
            }
            opt_start = 2;
        }
    }

    /* Parse options.  Recognised: GaussianIntegers -> True|False|Automatic,
     * Modulus -> integer (only 0 is meaningfully honoured for now).
     * Anything else (or unknown option name) marks `last_bad`; we report
     * the LAST one after the full scan to match Mathematica's diagnostic. */
    int gaussian_setting = 0;  /* 0 = Automatic (default), 1 = True, -1 = False */
    Expr* last_bad = NULL;
    for (size_t i = opt_start; i < argc; i++) {
        Expr* opt = res->data.function.args[i];
        if (!is_rule_head(opt) || opt->data.function.args[0]->type != EXPR_SYMBOL) {
            last_bad = opt;
            continue;
        }
        const char* name = opt->data.function.args[0]->data.symbol;
        Expr* val = opt->data.function.args[1];
        if (name == SYM_GaussianIntegers || strcmp(name, "GaussianIntegers") == 0) {
            if (is_sym_eq(val, "True")) gaussian_setting = 1;
            else if (is_sym_eq(val, "False")) gaussian_setting = -1;
            else if (is_sym_eq(val, "Automatic")) gaussian_setting = 0;
            else last_bad = opt;
        } else if (name == SYM_Modulus || strcmp(name, "Modulus") == 0) {
            /* Only Modulus -> 0 is implemented; any other value (non-zero
             * integer, non-integer, or symbolic) emits modnotimpl and the
             * call returns unevaluated.  Mathematica accepts arbitrary
             * integer moduli; that is a deferred extension. */
            bool is_zero_modulus = false;
            if (val->type == EXPR_INTEGER && val->data.integer == 0) {
                is_zero_modulus = true;
            } else if (val->type == EXPR_BIGINT && mpz_sgn(val->data.bigint) == 0) {
                is_zero_modulus = true;
            }
            if (!is_zero_modulus) {
                if (free_vars) free(vars);
                return sqfree_emit_modnotimpl(val);
            }
        } else {
            last_bad = opt;
        }
    }
    if (last_bad != NULL) {
        if (free_vars) free(vars);
        return sqfree_emit_nonopt(last_bad, 1, res);
    }

    /* Resolve GaussianIntegers -> Automatic: enable Z[i] mode iff the
     * input is a Gaussian-integer Complex literal. */
    bool gaussian;
    if (gaussian_setting == 1) gaussian = true;
    else if (gaussian_setting == -1) gaussian = false;
    else {
        mpz_t a, b;
        mpz_inits(a, b, NULL);
        gaussian = is_complex_integer(expr, a, b);
        mpz_clears(a, b, NULL);
    }

    bool result = sqfree_dispatch(expr, vars, var_count, gaussian);
    if (free_vars) free(vars);
    return expr_new_symbol(result ? "True" : "False");
}

/* ===================================================================== */
/* Init                                                                  */
/* ===================================================================== */

void squarefreeq_init(void) {
    symtab_add_builtin("SquareFreeQ", builtin_squarefreeq);
    symtab_get_def("SquareFreeQ")->attributes |= ATTR_PROTECTED;
    /* Docstring lives in info.c alongside the other *Q predicates. */
}
