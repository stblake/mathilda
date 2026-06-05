/*
 * piecewise.c
 * 
 * Implementation of Mathilda piecewise mathematical functions:
 * Floor, Ceiling, Round, IntegerPart, FractionalPart.
 * Suitable for numeric manipulation, automatically threading over lists.
 */

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include "piecewise.h"
#include "symtab.h"
#include "eval.h"
#include "arithmetic.h"
#include "complex.h"
#include "numeric.h"
#include "sym_names.h"

void piecewise_init(void) {
    symtab_add_builtin("Floor", builtin_floor);
    symtab_add_builtin("Ceiling", builtin_ceiling);
    symtab_add_builtin("Round", builtin_round);
    symtab_add_builtin("IntegerPart", builtin_integerpart);
    symtab_add_builtin("FractionalPart", builtin_fractionalpart);

    const char* funcs[] = {"Floor", "Ceiling", "Round", "IntegerPart", "FractionalPart", NULL};
    for (int i = 0; funcs[i] != NULL; i++) {
        symtab_get_def(funcs[i])->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
    }
}

enum { OP_FLOOR, OP_CEILING, OP_ROUND, OP_INTPART, OP_FRACPART };

static const char* op_head_name(int op) {
    switch (op) {
        case OP_FLOOR:   return "Floor";
        case OP_CEILING: return "Ceiling";
        case OP_ROUND:   return "Round";
        default:         return NULL;
    }
}

/* Return OP_FLOOR/OP_CEILING/OP_ROUND when `e` is a 1-arg call to one of
 * those heads, or -1 otherwise. Used to detect nested integer-valued
 * forms so the outer rounding becomes a no-op. */
static int classify_int_valued_head(Expr* e) {
    if (e->type != EXPR_FUNCTION) return -1;
    if (e->data.function.arg_count != 1) return -1;
    Expr* h = e->data.function.head;
    if (!h || h->type != EXPR_SYMBOL) return -1;
    if (h->data.symbol == SYM_Floor)   return OP_FLOOR;
    if (h->data.symbol == SYM_Ceiling) return OP_CEILING;
    if (h->data.symbol == SYM_Round)   return OP_ROUND;
    return -1;
}

/*
 * round_half_even:
 * Rounds numbers of the form x.5 toward the nearest even integer.
 */
static double round_half_even(double x) {
    double f = floor(x);
    double r = x - f;
    if (r < 0.5) return f;
    if (r > 0.5) return f + 1.0;
    // Exactly 0.5
    if (fmod(fabs(f), 2.0) == 0.0) return f;
    return f + 1.0;
}

static bool is_infinity(Expr* e) {
    return e->type == EXPR_SYMBOL && (e->data.symbol == SYM_Infinity || e->data.symbol == SYM_ComplexInfinity);
}

static bool is_minus_infinity(Expr* e) {
    if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 2 && 
        e->data.function.head->type == EXPR_SYMBOL && e->data.function.head->data.symbol == SYM_Times) {
        Expr* a1 = e->data.function.args[0];
        Expr* a2 = e->data.function.args[1];
        if (a1->type == EXPR_INTEGER && a1->data.integer == -1 && is_infinity(a2)) return true;
        if (a2->type == EXPR_INTEGER && a2->data.integer == -1 && is_infinity(a1)) return true;
    }
    return false;
}

#ifdef USE_MPFR
/* Resolve Floor/Ceiling/Round/IntegerPart of an *exact* real numeric
 * quantity (e.g. 10000000 3^(2/3), 25000000000000000000 Pi) that the leaf
 * branches cannot handle. Defined below; forward-declared so do_piecewise_1
 * can fall back to it. */
static Expr* do_piecewise_numeric_exact(Expr* x, int op);
#endif

/*
 * do_piecewise_1:
 * Applies the piecewise operation to a single numeric argument.
 */
static Expr* do_piecewise_1(Expr* x, int op) {
    if (is_infinity(x) || is_minus_infinity(x)) {
        if (op == OP_FRACPART) return expr_new_integer(0);
        return expr_copy(x);
    }

    if (x->type == EXPR_INTEGER) {
        if (op == OP_FRACPART) return expr_new_integer(0);
        return expr_copy(x);
    }

    /* GMP BigInt: integer-valued for all five ops. */
    if (x->type == EXPR_BIGINT) {
        if (op == OP_FRACPART) return expr_new_integer(0);
        return expr_copy(x);
    }
    
    if (x->type == EXPR_REAL) {
        double v = x->data.real;
        double res = 0.0;
        if (op == OP_FLOOR) res = floor(v);
        else if (op == OP_CEILING) res = ceil(v);
        else if (op == OP_ROUND) res = round_half_even(v);
        else if (op == OP_INTPART) res = trunc(v);
        else if (op == OP_FRACPART) return expr_new_real(v - trunc(v));

        return expr_new_integer((int64_t)res);
    }

#ifdef USE_MPFR
    if (x->type == EXPR_MPFR) {
        mpfr_prec_t prec = mpfr_get_prec(x->data.mpfr);
        if (op == OP_FRACPART) {
            /* FractionalPart preserves precision. Result = x - trunc(x). */
            mpfr_t tmp;
            mpfr_init2(tmp, prec);
            mpfr_trunc(tmp, x->data.mpfr);
            Expr* result = expr_new_mpfr_bits(prec);
            mpfr_sub(result->data.mpfr, x->data.mpfr, tmp, MPFR_RNDN);
            mpfr_clear(tmp);
            return result;
        }
        /* Floor/Ceiling/Round/IntegerPart all produce an exact integer.
         * Route through mpz_t so arbitrarily large MPFR values do not
         * silently truncate when they overflow int64_t. */
        mpfr_t r;
        mpfr_init2(r, prec);
        if (op == OP_FLOOR)        mpfr_floor(r, x->data.mpfr);
        else if (op == OP_CEILING) mpfr_ceil(r, x->data.mpfr);
        else if (op == OP_ROUND)   mpfr_rint(r, x->data.mpfr, MPFR_RNDN);
        else /* OP_INTPART */      mpfr_trunc(r, x->data.mpfr);
        mpz_t z;
        mpz_init(z);
        mpfr_get_z(z, r, MPFR_RNDN);
        Expr* result = expr_bigint_normalize(expr_new_bigint_from_mpz(z));
        mpz_clear(z);
        mpfr_clear(r);
        return result;
    }
#endif
    
    int64_t n, d;
    if (is_rational(x, &n, &d)) {
        double v = (double)n / d;
        double res = 0.0;
        if (op == OP_FLOOR) res = floor(v);
        else if (op == OP_CEILING) res = ceil(v);
        else if (op == OP_ROUND) res = round_half_even(v);
        else if (op == OP_INTPART) res = trunc(v);
        else if (op == OP_FRACPART) {
            int64_t int_part = (int64_t)trunc(v);
            return make_rational(n - int_part * d, d);
        }
        return expr_new_integer((int64_t)res);
    }

    /* Rational[BigInt, _] / Rational[_, BigInt] / Rational[BigInt, BigInt] —
     * the int64-only is_rational() refuses to extract these, so do exact
     * mpz arithmetic. */
    if (x->type == EXPR_FUNCTION
        && x->data.function.head->type == EXPR_SYMBOL
        && x->data.function.head->data.symbol == SYM_Rational
        && x->data.function.arg_count == 2
        && expr_is_integer_like(x->data.function.args[0])
        && expr_is_integer_like(x->data.function.args[1])) {
        mpz_t num, den;
        expr_to_mpz(x->data.function.args[0], num);
        expr_to_mpz(x->data.function.args[1], den);
        if (mpz_sgn(den) == 0) {
            mpz_clears(num, den, NULL);
            return NULL;
        }
        /* Canonical Rationals carry a positive denominator, but tolerate
         * a structural Rational[*, -k] reaching here from a builder. */
        if (mpz_sgn(den) < 0) {
            mpz_neg(num, num);
            mpz_neg(den, den);
        }
        if (op == OP_FRACPART) {
            mpz_t int_num, frac_num;
            mpz_inits(int_num, frac_num, NULL);
            mpz_tdiv_q(int_num, num, den);      /* int_num = trunc(num/den) */
            mpz_mul(int_num, int_num, den);
            mpz_sub(frac_num, num, int_num);    /* frac_num = num - trunc*den */
            Expr* num_e = expr_bigint_normalize(expr_new_bigint_from_mpz(frac_num));
            Expr* den_e = expr_bigint_normalize(expr_new_bigint_from_mpz(den));
            Expr* args[2] = { num_e, den_e };
            mpz_clears(num, den, int_num, frac_num, NULL);
            return eval_and_free(expr_new_function(expr_new_symbol("Rational"), args, 2));
        }
        mpz_t q;
        mpz_init(q);
        if (op == OP_FLOOR)        mpz_fdiv_q(q, num, den);
        else if (op == OP_CEILING) mpz_cdiv_q(q, num, den);
        else if (op == OP_INTPART) mpz_tdiv_q(q, num, den);
        else /* OP_ROUND, banker's */ {
            /* Round half to even. floor((2*num + den) / (2*den)) gives
             * round-half-up; correct it when the residue is exactly half. */
            mpz_t two_num, two_den, two_num_plus_den, rem, half_den;
            mpz_inits(two_num, two_den, two_num_plus_den, rem, half_den, NULL);
            mpz_mul_ui(two_num, num, 2);
            mpz_mul_ui(two_den, den, 2);
            mpz_add(two_num_plus_den, two_num, den);
            mpz_fdiv_qr(q, rem, two_num_plus_den, two_den);
            /* Tie case: rem == 0 means the input was exactly q + 1/2; force
             * q to even. */
            if (mpz_sgn(rem) == 0 && mpz_odd_p(q)) {
                /* Round toward the nearest even by subtracting 1 when q is odd
                 * and positive (or adding 1 when q is odd and negative). The
                 * floor-based formula systematically biases away from zero on
                 * the half-step, so the adjustment is -1 for even-but-odd-q
                 * always: this brings q toward the even neighbour. */
                mpz_sub_ui(q, q, 1);
            }
            mpz_clears(two_num, two_den, two_num_plus_den, rem, half_den, NULL);
        }
        Expr* out = expr_bigint_normalize(expr_new_bigint_from_mpz(q));
        mpz_clears(num, den, q, NULL);
        return out;
    }
    
    Expr *re, *im;
    if (is_complex(x, &re, &im)) {
        Expr* re_res = do_piecewise_1(re, op);
        Expr* im_res = do_piecewise_1(im, op);
        if (re_res && im_res) {
            Expr* combined = make_complex(re_res, im_res);
            return combined;
        }
        if (re_res) expr_free(re_res);
        if (im_res) expr_free(im_res);
    }

    /* Symbolic simplifications for Floor/Ceiling/Round only.
     * Any concrete numeric input has already been resolved above. */
    if (op == OP_FLOOR || op == OP_CEILING || op == OP_ROUND) {
        /* Idempotency / composition:
         *   Floor[Floor[y]]    -> Floor[y]
         *   Ceiling[Ceiling[y]]-> Ceiling[y]
         *   Round[Round[y]]    -> Round[y]
         *   Floor[Ceiling[y]]  -> Ceiling[y]
         *   Floor[Round[y]]    -> Round[y]
         *   Ceiling[Floor[y]]  -> Floor[y]
         *   Ceiling[Round[y]]  -> Round[y]
         *   Round[Floor[y]]    -> Floor[y]
         *   Round[Ceiling[y]]  -> Ceiling[y]
         * In every case the inner result is an integer, so the outer
         * rounding is a no-op and the inner expression is the answer. */
        if (classify_int_valued_head(x) >= 0) {
            return expr_copy(x);
        }

        /* Sign extraction:
         *   Floor[-y]   -> -Ceiling[y]
         *   Ceiling[-y] -> -Floor[y]
         *   Round[-y]   -> -Round[y]
         * Triggered by any superficially negative argument (Times[-c, ...],
         * Complex[neg, ...] that wasn't fully resolved above, etc.). */
        if (expr_is_superficially_negative(x)) {
            int swapped_op = (op == OP_FLOOR) ? OP_CEILING
                          : (op == OP_CEILING) ? OP_FLOOR
                          : OP_ROUND;
            Expr* neg_args[2] = { expr_new_integer(-1), expr_copy(x) };
            Expr* pos = eval_and_free(expr_new_function(expr_new_symbol("Times"), neg_args, 2));
            Expr* inner_args[1] = { pos };
            Expr* inner = expr_new_function(expr_new_symbol(op_head_name(swapped_op)), inner_args, 1);
            Expr* outer_args[2] = { expr_new_integer(-1), inner };
            return expr_new_function(expr_new_symbol("Times"), outer_args, 2);
        }
    }

#ifdef USE_MPFR
    /* Exact real numeric argument that no leaf branch resolved (e.g.
     * Round[10000000 3^(2/3)], Floor[25000000000000000000 Pi]). Certify
     * the result by high-precision numeric evaluation. */
    if (op == OP_FRACPART) {
        /* FractionalPart[x] = x - IntegerPart[x], kept exact (matches
         * Mathematica: FractionalPart[10000000 3^(2/3)] stays symbolic as
         * 10000000 3^(2/3) - 20800838). */
        Expr* ip = do_piecewise_numeric_exact(x, OP_INTPART);
        if (ip) {
            Expr* neg_args[2] = { expr_new_integer(-1), ip };
            Expr* neg_ip = expr_new_function(expr_new_symbol("Times"), neg_args, 2);
            Expr* sum_args[2] = { expr_copy(x), neg_ip };
            return eval_and_free(expr_new_function(expr_new_symbol("Plus"), sum_args, 2));
        }
        return NULL;
    }
    {
        Expr* r = do_piecewise_numeric_exact(x, op);
        if (r) return r;
    }
#endif

    return NULL;
}

#ifdef USE_MPFR
/*
 * do_piecewise_numeric_exact:
 *   Floor/Ceiling/Round/IntegerPart of an exact real numeric quantity the
 *   leaf branches above could not resolve. We numericalize `x` to MPFR at
 *   increasing precision and apply the operation, accepting the result only
 *   once two successive precisions agree on the integer — an interval-style
 *   certification that guards against a value sitting arbitrarily close to a
 *   rounding boundary. Returns NULL (leaving the call unevaluated, never
 *   wrong) when `x` is not a pure real number or convergence is not reached
 *   within the precision cap.
 */
static Expr* do_piecewise_numeric_exact(Expr* x, int op) {
    const long max_prec = 1L << 16;   /* ~19,700 decimal digits */
    Expr* prev = NULL;                /* previous iteration's integer result */
    for (long prec = 256; prec <= max_prec; prec *= 2) {
        NumericSpec spec;
        spec.mode = NUMERIC_MODE_MPFR;
        spec.bits = prec;
        Expr* approx = numericalize(x, spec);
        /* Must collapse to a pure real number; anything else (still
         * symbolic, or Complex[...]) means `x` is not a real numeric and we
         * have no business rounding it. */
        if (!approx || (approx->type != EXPR_MPFR && approx->type != EXPR_REAL
                        && approx->type != EXPR_INTEGER
                        && approx->type != EXPR_BIGINT)) {
            if (approx) expr_free(approx);
            break;
        }
        Expr* cur = do_piecewise_1(approx, op);   /* MPFR/real -> exact int */
        expr_free(approx);
        if (!cur) break;
        if (prev && expr_eq(prev, cur)) {
            expr_free(prev);
            return cur;
        }
        if (prev) expr_free(prev);
        prev = cur;
    }
    if (prev) expr_free(prev);
    return NULL;
}
#endif

static Expr* make_divide(Expr* num, Expr* den) {
    Expr* pow_args[2] = { expr_copy(den), expr_new_integer(-1) };
    Expr* inv = expr_new_function(expr_new_symbol("Power"), pow_args, 2);
    Expr* times_args[2] = { expr_copy(num), inv };
    return expr_new_function(expr_new_symbol("Times"), times_args, 2);
}

static Expr* do_piecewise(Expr* res, int op, const char* name, bool allow_2_args) {
    if (res->type != EXPR_FUNCTION) return NULL;
    
    if (res->data.function.arg_count == 1) {
        return do_piecewise_1(res->data.function.args[0], op);
    } else if (allow_2_args && res->data.function.arg_count == 2) {
        Expr* x = res->data.function.args[0];
        Expr* a = res->data.function.args[1];
        
        Expr* div = make_divide(x, a);
        Expr* func_args[1] = { div }; // div ownership passes to function
        Expr* func = expr_new_function(expr_new_symbol(name), func_args, 1);
        
        Expr* times_args[2] = { expr_copy(a), func };
        return expr_new_function(expr_new_symbol("Times"), times_args, 2);
    }
    
    return NULL;
}

Expr* builtin_floor(Expr* res) { return do_piecewise(res, OP_FLOOR, "Floor", true); }
Expr* builtin_ceiling(Expr* res) { return do_piecewise(res, OP_CEILING, "Ceiling", true); }
Expr* builtin_round(Expr* res) { return do_piecewise(res, OP_ROUND, "Round", true); }
Expr* builtin_integerpart(Expr* res) { return do_piecewise(res, OP_INTPART, "IntegerPart", false); }
Expr* builtin_fractionalpart(Expr* res) { return do_piecewise(res, OP_FRACPART, "FractionalPart", false); }
