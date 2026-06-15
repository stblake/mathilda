
#include "arithmetic.h"
#include "eval.h"
#include "sym_names.h"
#include "internal.h"
#include "print.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <gmp.h>

int g_arith_warnings_muted = 0;

int64_t gcd(int64_t a, int64_t b) {
    a = llabs(a);
    b = llabs(b);
    while (b) {
        a %= b;
        int64_t tmp = a;
        a = b;
        b = tmp;
    }
    return a;
}

int64_t lcm(int64_t a, int64_t b) {
    if (a == 0 || b == 0) return 0;
    a = llabs(a);
    b = llabs(b);
    return (a / gcd(a, b)) * b;
}

Expr* make_rational(int64_t n, int64_t d) {
    if (d == 0) return NULL; // Error
    if (n == 0) return expr_new_integer(0);
    
    int64_t common = gcd(n, d);
    n /= common;
    d /= common;

    if (d < 0) {
        n = -n;
        d = -d;
    }

    if (d == 1) return expr_new_integer(n);

    Expr* args[2];
    args[0] = expr_new_integer(n);
    args[1] = expr_new_integer(d);
    return expr_new_function(expr_new_symbol(SYM_Rational), args, 2);
}

Expr* builtin_rational(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* n_expr = res->data.function.args[0];
    Expr* d_expr = res->data.function.args[1];
    
    if (n_expr->type == EXPR_INTEGER && d_expr->type == EXPR_INTEGER) {
        int64_t n = n_expr->data.integer;
        int64_t d = d_expr->data.integer;
        if (d == 0) {
            if (!arith_warnings_muted())
                fprintf(stderr, "Power::infy: Infinite expression 1/0 encountered.\n");
            if (n == 0) {
                if (!arith_warnings_muted())
                    fprintf(stderr,
                        "Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.\n");
                return expr_new_symbol(SYM_Indeterminate);
            }
            return expr_new_symbol(SYM_ComplexInfinity);
        }
        
        Expr* r = make_rational(n, d);
        if (r && r->type == EXPR_FUNCTION && r->data.function.head->type == EXPR_SYMBOL && r->data.function.head->data.symbol == SYM_Rational) {
            Expr* rn = r->data.function.args[0];
            Expr* rd = r->data.function.args[1];
            if (rn->type == EXPR_INTEGER && rd->type == EXPR_INTEGER && rn->data.integer == n && rd->data.integer == d) {
                // No simplification happened
                expr_free(r);
                return NULL;
            }
        }
        return r;
    }
    return NULL;
}

bool is_rational(const Expr* e, int64_t* n, int64_t* d) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) {
        if (n) *n = e->data.integer;
        if (d) *d = 1;
        return true;
    }
    /* Function nodes built during evaluator transitions can transiently
     * carry a NULL head pointer; guard before dereferencing.  Without
     * this check, a Function with a NULL head reaching this fast-path
     * crashes the dispatcher (observed on Linux under nested-radical
     * Simplify; macOS heap layout makes the deref usually land in
     * mapped memory). */
    if (e->type == EXPR_FUNCTION && e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Rational) {
        if (e->data.function.arg_count == 2 &&
            e->data.function.args[0] &&
            e->data.function.args[1] &&
            e->data.function.args[0]->type == EXPR_INTEGER &&
            e->data.function.args[1]->type == EXPR_INTEGER) {
            if (n) *n = e->data.function.args[0]->data.integer;
            if (d) *d = e->data.function.args[1]->data.integer;
            return true;
        }
    }
    return false;
}

/* Bigint-aware version: matches Integer, BigInt, or Rational[X, Y] where
 * both components are integer-like (Integer or BigInt).  Used in numeric
 * predicates so the Times/Plus folders recognise rationals whose
 * components have overflowed int64 — without this, a BigInt times a
 * Rational[1, BigInt] is left unsimplified and can defeat exact
 * polynomial division. */
bool is_rational_like(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) return true;
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Rational &&
        e->data.function.arg_count == 2 &&
        expr_is_integer_like(e->data.function.args[0]) &&
        expr_is_integer_like(e->data.function.args[1])) {
        return true;
    }
    return false;
}

bool is_complex(Expr* e, Expr** re, Expr** im) {
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Complex &&
        e->data.function.arg_count == 2) {
        if (re) *re = e->data.function.args[0];
        if (im) *im = e->data.function.args[1];
        return true;
    }
    return false;
}

Expr* make_complex(Expr* re, Expr* im) {
    Expr* args[2] = { re, im };
    return expr_new_function(expr_new_symbol(SYM_Complex), args, 2);
}

Expr* builtin_subtract(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;

    Expr* a = res->data.function.args[0];
    Expr* b = res->data.function.args[1];

    Expr* minus_one = expr_new_integer(-1);
    Expr* mb_args[2] = { minus_one, expr_copy(b) };
    Expr* minus_b = expr_new_function(expr_new_symbol(SYM_Times), mb_args, 2);

    Expr* p_args[2] = { expr_copy(a), minus_b };
    return expr_new_function(expr_new_symbol(SYM_Plus), p_args, 2);
}

Expr* builtin_complex(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;

    Expr* re = res->data.function.args[0];
    Expr* im = res->data.function.args[1];

    if (im->type == EXPR_INTEGER && im->data.integer == 0) {
        return expr_copy(re);
    }
    if (im->type == EXPR_REAL && im->data.real == 0.0) {
        if (re->type == EXPR_INTEGER) {
            return expr_new_real((double)re->data.integer);
        }
        return expr_copy(re);
    }

    return NULL;
}

Expr* builtin_divide(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    
    Expr* num = res->data.function.args[0];
    Expr* den = res->data.function.args[1];

    if (num->type == EXPR_REAL || den->type == EXPR_REAL) {
        double vnum = (num->type == EXPR_REAL) ? num->data.real : (num->type == EXPR_INTEGER) ? (double)num->data.integer : (num->type == EXPR_BIGINT) ? mpz_get_d(num->data.bigint) : 0.0;
        double vden = (den->type == EXPR_REAL) ? den->data.real : (den->type == EXPR_INTEGER) ? (double)den->data.integer : (den->type == EXPR_BIGINT) ? mpz_get_d(den->data.bigint) : 0.0;
        if (vden == 0.0) {
            if (!arith_warnings_muted())
                fprintf(stderr, "Power::infy: Infinite expression 1/0 encountered.\n");
            return expr_new_symbol(SYM_ComplexInfinity);
        }
        return expr_new_real(vnum / vden);
    }

    int64_t n1, d1, n2, d2;
    if (is_rational(num, &n1, &d1) && is_rational(den, &n2, &d2)) {
        if (n2 == 0) {
            /* x / 0 with rational/integer x: 0/0 -> Indeterminate (handled in
             * Times when 0 multiplies ComplexInfinity); otherwise emit the
             * Power::infy message and yield ComplexInfinity. */
            if (!arith_warnings_muted())
                fprintf(stderr, "Power::infy: Infinite expression 1/0 encountered.\n");
            if (n1 == 0) {
                if (!arith_warnings_muted())
                    fprintf(stderr,
                        "Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.\n");
                return expr_new_symbol(SYM_Indeterminate);
            }
            return expr_new_symbol(SYM_ComplexInfinity);
        }
        Expr* r = make_rational(n1 * d2, d1 * n2);
        if (r) return r;
    }

    Expr* minus_one = expr_new_integer(-1);
    Expr* p_args[2] = { expr_copy(den), minus_one };
    Expr* power = expr_new_function(expr_new_symbol(SYM_Power), p_args, 2);
    
    Expr* t_args[2] = { expr_copy(num), power };
    return expr_new_function(expr_new_symbol(SYM_Times), t_args, 2);
}

bool is_infinity_sym(Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol == SYM_Infinity;
}

bool is_complex_infinity_sym(Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol == SYM_ComplexInfinity;
}

bool is_indeterminate_sym(Expr* e) {
    return e && e->type == EXPR_SYMBOL && e->data.symbol == SYM_Indeterminate;
}

int expr_numeric_sign(Expr* e) {
    if (!e) return 0;
    if (e->type == EXPR_INTEGER) {
        if (e->data.integer > 0) return 1;
        if (e->data.integer < 0) return -1;
        return 0;
    }
    if (e->type == EXPR_REAL) {
        if (e->data.real > 0.0) return 1;
        if (e->data.real < 0.0) return -1;
        return 0;
    }
    if (e->type == EXPR_BIGINT) return mpz_sgn(e->data.bigint);
    int64_t n, d;
    if (is_rational(e, &n, &d)) {
        // d is conventionally positive in Mathilda Rational[n, d].
        if (n > 0) return (d > 0) ? 1 : -1;
        if (n < 0) return (d > 0) ? -1 : 1;
        return 0;
    }
    return 0;
}

bool expr_is_superficially_negative(Expr* e) {
    if (!e) return false;
    int s = expr_numeric_sign(e);
    if (s < 0) return true;
    if (s > 0) return false;
    /* s == 0: either a zero numeric (not negative) or a non-numeric --
     * fall through to Complex and Times shape checks. */
    Expr* re; Expr* im;
    if (is_complex(e, &re, &im)) {
        int rs = expr_numeric_sign(re);
        if (rs < 0) return true;
        if (rs > 0) return false;
        /* Real part is zero: use the imaginary part's sign. Catches the
         * "pure imaginary with negative coefficient" case, e.g. -2 I =
         * Complex[0, -2]. */
        return expr_numeric_sign(im) < 0;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Times &&
        e->data.function.arg_count > 0) {
        /* Leading factor carries the syntactic sign, per Times canonical
         * ordering (numerics sort first). */
        return expr_is_superficially_negative(e->data.function.args[0]);
    }
    return false;
}

bool is_neg_infinity_form(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    if (e->data.function.head->data.symbol != SYM_Times) return false;
    if (e->data.function.arg_count != 2) return false;
    Expr* a = e->data.function.args[0];
    Expr* b = e->data.function.args[1];
    if (!is_infinity_sym(b)) return false;
    return expr_numeric_sign(a) < 0;
}
