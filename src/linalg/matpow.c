/* MatrixPower[m, e] and MatrixPower[m, e, v].
 *
 * Implements integer matrix exponentiation via square-and-multiply.
 * Negative exponents are handled by composing with Inverse[m] (which lives
 * in matinv.c).  Fractional exponents are flagged as unsupported.
 */

#include "linalg.h"
#include "ndlinalg.h"
#include "ndarray.h"
#include "common.h"
#include "pack.h"
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

    /* Validate matrix: must be rank 2, square, non-empty.
     *
     * A BUFFER carries its own shape. This head used to open with
     * linalg_delist_and_reeval, which threw the buffer away and left every
     * mat_dot below contracting boxed Expr nodes: MatrixPower[A300, 4] on
     * random Reals took 14.7 s, where the same three products written A.A.A.A
     * reach dgemm and cost 2.6 ms each (plan item C.2). Nothing else in this
     * function reads an element -- the whole algorithm is square-and-multiply
     * over dot2, and dot2 has had the BLAS path since the Phase 1 routing work
     * -- so keeping the buffer needs only a shape probe that understands one.
     *
     * An INTEGER matrix still materialises, upstream: MatrixPower is not on
     * pack.c's INT64_OK list, so the transparency gate hands this function a
     * plain List and the exact path runs. That is the intended split -- the
     * exact fourth power of an integer matrix is an integer matrix, and a
     * float64 buffer cannot hold one past 2^53. */
    int64_t dims[64];
    int rank;
    if (is_ndarray(m)) {
        rank = m->data.ndarray.rank;
        for (int i = 0; i < rank && i < 64; i++) dims[i] = m->data.ndarray.dims[i];
    } else {
        rank = get_tensor_dims(m, dims);
    }
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
        && exp_arg->data.function.head->data.symbol.name == SYM_Rational) {
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
        int vrank;
        if (is_ndarray(vec)) {
            vrank = vec->data.ndarray.rank;
            for (int i = 0; i < vrank && i < 64; i++) vdims[i] = vec->data.ndarray.dims[i];
        } else {
            vrank = get_tensor_dims(vec, vdims);
        }
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
            && inv_result->data.function.head->data.symbol.name == SYM_Inverse) {
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
        /* M^0 IS THE IDENTITY, AT M's OWN EXACTNESS. Mathematica's
         * MatrixPower[{{2., 0.}, {0., 2.}}, 0] is {{1., 0.}, {0., 1.}} and the
         * exact matrix's is {{1, 0}, {0, 1}}. Going through IdentityMatrix[n]
         * gave the exact one either way, which for a machine matrix is both
         * wrong and unpackable -- see common.h on machine-real contagion. */
        expr_free(base);
        if (common_has_machine_real(m)) {
            int64_t dims[2] = { n, n };
            void* raw = NULL;
            result = ndbuild_open(2, dims, NDT_FLOAT64, &raw);
            if (result) {
                double* b = (double*)raw;
                for (int64_t z = 0; z < (int64_t)n * n; z++) b[z] = 0.0;
                for (int i = 0; i < n; i++) b[(int64_t)i * n + i] = 1.0;
            } else {
                /* Under the packing threshold, or packing off: the same matrix,
                 * boxed. */
                Expr** rows = malloc(sizeof(Expr*) * (size_t)n);
                for (int i = 0; i < n; i++) {
                    Expr** cells = malloc(sizeof(Expr*) * (size_t)n);
                    for (int j = 0; j < n; j++)
                        cells[j] = expr_new_real(i == j ? 1.0 : 0.0);
                    rows[i] = expr_new_function(expr_new_symbol(SYM_List), cells, (size_t)n);
                    free(cells);
                }
                result = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)n);
                free(rows);
            }
        } else {
            Expr* id_args[1];
            id_args[0] = expr_new_integer(n);
            Expr* id_call = expr_new_function(expr_new_symbol(SYM_IdentityMatrix), id_args, 1);
            result = evaluate(id_call);
            expr_free(id_call);
        }
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
