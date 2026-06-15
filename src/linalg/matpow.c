/* MatrixPower[m, e] and MatrixPower[m, e, v].
 *
 * Implements integer matrix exponentiation via square-and-multiply.
 * Negative exponents are handled by composing with Inverse[m] (which lives
 * in matinv.c).  Fractional exponents are flagged as unsupported.
 */

#include "linalg.h"
#include "eval.h"
#include "print.h"
#include "sym_names.h"
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

/* A.B for two square n x n matrices.  Thin wrapper around dot2 that
 * discards the (always-false here) error-printed channel. */
static Expr* mat_dot(Expr* a, Expr* b) {
    bool err = false;
    return dot2(a, b, &err);
}

Expr* builtin_matrixpower(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 3) return NULL;

    Expr* m = res->data.function.args[0];
    Expr* exp_arg = res->data.function.args[1];

    /* Validate matrix: must be rank 2, square, non-empty */
    int64_t dims[64];
    int rank = get_tensor_dims(m, dims);
    if (rank != 2 || dims[0] != dims[1] || dims[0] == 0) {
        if (rank >= 0) {
            char* m_str = expr_to_string(m);
            fprintf(stderr, "MatrixPower::matsq: Argument %s at position 1 is not a non-empty square matrix.\n", m_str);
            free(m_str);
        }
        return NULL;
    }
    int n = (int)dims[0];

    /* Validate exponent: must be an integer */
    bool is_int = (exp_arg->type == EXPR_INTEGER);
    bool is_bigint = (exp_arg->type == EXPR_BIGINT);
    bool is_rational = false;
    bool is_real = (exp_arg->type == EXPR_REAL);

    if (exp_arg->type == EXPR_FUNCTION && exp_arg->data.function.head->type == EXPR_SYMBOL
        && exp_arg->data.function.head->data.symbol == SYM_Rational) {
        is_rational = true;
    }

    /* Fractional powers: warn and return unevaluated */
    if (is_rational || is_real) {
        char* exp_str = expr_to_string(exp_arg);
        fprintf(stderr, "MatrixPower::fract: Fractional matrix powers are not currently supported. Exponent %s is not an integer.\n", exp_str);
        free(exp_str);
        return NULL;
    }

    /* Symbolic or non-numeric exponent: return unevaluated */
    if (!is_int && !is_bigint) return NULL;

    /* Get the exponent value; for bigint, check it fits in int64_t */
    int64_t exp_val = 0;
    if (is_int) {
        exp_val = exp_arg->data.integer;
    } else {
        /* EXPR_BIGINT: check if it fits in int64_t range */
        if (mpz_fits_slong_p(exp_arg->data.bigint)) {
            exp_val = mpz_get_si(exp_arg->data.bigint);
        } else {
            /* Exponent too large to compute */
            return NULL;
        }
    }

    /* If argc == 3, validate vector argument */
    Expr* vec = NULL;
    if (argc == 3) {
        vec = res->data.function.args[2];
        int64_t vdims[64];
        int vrank = get_tensor_dims(vec, vdims);
        if (vrank != 1 || vdims[0] != n) {
            char* v_str = expr_to_string(vec);
            fprintf(stderr, "MatrixPower::vecsh: Vector %s has incompatible length for matrix of size %d.\n", v_str, n);
            free(v_str);
            return NULL;
        }
    }

    /* For negative exponents, compute inverse first */
    Expr* base = NULL;
    int64_t abs_exp = exp_val;
    if (exp_val < 0) {
        abs_exp = -exp_val;
        /* Compute Inverse[m] */
        Expr* inv_call = expr_new_function(expr_new_symbol(SYM_Inverse),
            (Expr*[]){expr_copy(m)}, 1);
        Expr* inv_result = evaluate(inv_call);
        expr_free(inv_call);

        /* Check if Inverse returned unevaluated (singular matrix) */
        if (inv_result->type == EXPR_FUNCTION && inv_result->data.function.head->type == EXPR_SYMBOL
            && inv_result->data.function.head->data.symbol == SYM_Inverse) {
            expr_free(inv_result);
            return NULL; /* Singular: Inverse already printed warning */
        }
        base = inv_result;
    } else {
        base = expr_copy(m);
    }

    /* Compute matrix power via binary exponentiation */
    Expr* result = NULL;

    if (abs_exp == 0) {
        /* M^0 = IdentityMatrix[n] */
        expr_free(base);
        Expr* id_args[1];
        id_args[0] = expr_new_integer(n);
        Expr* id_call = expr_new_function(expr_new_symbol(SYM_IdentityMatrix), id_args, 1);
        result = evaluate(id_call);
        expr_free(id_call);
    } else {
        /* Binary exponentiation: square-and-multiply */
        result = NULL;
        Expr* sq = base; /* current square, takes ownership of base */

        int64_t e = abs_exp;
        while (e > 0) {
            if (e & 1) {
                if (result == NULL) {
                    result = expr_copy(sq);
                } else {
                    Expr* new_result = mat_dot(result, sq);
                    if (!new_result) {
                        /* Should not happen for valid square matrices */
                        expr_free(result);
                        expr_free(sq);
                        return NULL;
                    }
                    Expr* evaluated = evaluate(new_result);
                    expr_free(new_result);
                    expr_free(result);
                    result = evaluated;
                }
            }
            e >>= 1;
            if (e > 0) {
                Expr* new_sq = mat_dot(sq, sq);
                if (!new_sq) {
                    if (result) expr_free(result);
                    expr_free(sq);
                    return NULL;
                }
                Expr* evaluated = evaluate(new_sq);
                expr_free(new_sq);
                expr_free(sq);
                sq = evaluated;
            }
        }
        expr_free(sq);
    }

    /* If vector argument provided, compute result . v */
    if (vec && result) {
        Expr* dotted = mat_dot(result, vec);
        if (dotted) {
            Expr* evaluated = evaluate(dotted);
            expr_free(dotted);
            expr_free(result);
            result = evaluated;
        }
        /* If dot fails, just return matrix power without applying to vector */
    }

    /* Note: evaluator frees res after we return non-NULL */
    return result;
}
