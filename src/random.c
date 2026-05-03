/*
 * random.c - Random number generation builtins for PicoCAS
 *
 * RandomInteger[{imin, imax}]  - pseudorandom integer in [imin, imax]
 * RandomInteger[imax]          - pseudorandom integer in [0, imax]
 * RandomInteger[]              - pseudorandomly gives 0 or 1
 * RandomInteger[range, n]      - list of n pseudorandom integers
 * RandomInteger[range, {n1, n2, ...}] - n1*n2*... array
 *
 * RandomReal[]                 - pseudorandom real in [0, 1)
 * RandomReal[xmax]             - pseudorandom real in [0, xmax)
 * RandomReal[{xmin, xmax}]    - pseudorandom real in [xmin, xmax)
 * RandomReal[range, n]         - list of n pseudorandom reals
 * RandomReal[range, {n1, n2, ...}] - n1*n2*... array
 *
 * RandomComplex[]              - pseudorandom complex in rectangle [0,1]x[0,1]
 * RandomComplex[zmax]          - pseudorandom complex in rectangle [0,zmax]
 * RandomComplex[{zmin, zmax}]  - pseudorandom complex in rectangle [zmin,zmax]
 * RandomComplex[range, n]      - list of n pseudorandom complex numbers
 * RandomComplex[range, {n1, n2, ...}] - n1*n2*... array
 *
 * RandomChoice[list]           - pseudorandom choice from list
 * RandomChoice[list, n]        - list of n pseudorandom choices
 * RandomChoice[list, {n1, n2, ...}] - n1*n2*... array of choices
 * RandomChoice[wts->elems]     - weighted pseudorandom choice
 * RandomChoice[wts->elems, n]  - list of n weighted choices
 * RandomChoice[wts->elems, {n1,...}] - array of weighted choices
 *
 * SeedRandom[n]                - seeds the RNG with integer n
 * SeedRandom[]                 - reseeds from system entropy
 *
 * Uses GMP's random state (gmp_randstate_t) for bignum support.
 */

#include "random.h"
#include "arithmetic.h"
#include "symtab.h"
#include "eval.h"
#include "attr.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <gmp.h>

/* Global GMP random state */
static gmp_randstate_t g_rand_state;
static int g_rand_initialized = 0;

/* Ensure RNG is initialized */
static void ensure_rand_init(void) {
    if (!g_rand_initialized) {
        gmp_randinit_mt(g_rand_state);
        gmp_randseed_ui(g_rand_state, (unsigned long)time(NULL) ^ (unsigned long)clock());
        g_rand_initialized = 1;
    }
}

/*
 * Generate a single random integer in [imin, imax] (inclusive).
 * imin and imax are given as mpz_t values.
 * Returns a new Expr* (EXPR_INTEGER or EXPR_BIGINT, normalized).
 */
static Expr* random_integer_range(const mpz_t imin, const mpz_t imax) {
    ensure_rand_init();

    /* range_size = imax - imin + 1 */
    mpz_t range_size, result;
    mpz_init(range_size);
    mpz_init(result);

    mpz_sub(range_size, imax, imin);
    mpz_add_ui(range_size, range_size, 1);

    /* If range_size <= 0, the range is empty or invalid */
    if (mpz_sgn(range_size) <= 0) {
        mpz_clear(range_size);
        mpz_clear(result);
        return NULL;
    }

    /* result = random in [0, range_size) */
    mpz_urandomm(result, g_rand_state, range_size);

    /* result = result + imin */
    mpz_add(result, result, imin);

    Expr* e = expr_bigint_normalize(expr_new_bigint_from_mpz(result));
    mpz_clear(range_size);
    mpz_clear(result);
    return e;
}

/*
 * Parse a range argument into imin, imax (as mpz_t).
 * Supports:
 *   - Integer or BigInt atom: range is [0, val]
 *   - List[imin, imax]: range is [imin, imax]
 * Returns true on success, false if the argument doesn't match.
 */
static bool parse_range(Expr* arg, mpz_t imin, mpz_t imax) {
    if (arg->type == EXPR_INTEGER) {
        mpz_set_si(imin, 0);
        mpz_set_si(imax, arg->data.integer);
        return true;
    }
    if (arg->type == EXPR_BIGINT) {
        mpz_set_si(imin, 0);
        mpz_set(imax, arg->data.bigint);
        return true;
    }
    /* Check for List[imin, imax] */
    if (arg->type == EXPR_FUNCTION &&
        arg->data.function.head->type == EXPR_SYMBOL &&
        arg->data.function.head->data.symbol == SYM_List &&
        arg->data.function.arg_count == 2) {

        Expr* lo = arg->data.function.args[0];
        Expr* hi = arg->data.function.args[1];

        if (!expr_is_integer_like(lo) || !expr_is_integer_like(hi))
            return false;

        expr_to_mpz(lo, imin);
        expr_to_mpz(hi, imax);
        return true;
    }
    return false;
}

/*
 * Build a multi-dimensional array of random integers.
 * dims is a list expression {n1, n2, ...}, dim_idx is current dimension.
 */
static Expr* random_array(const mpz_t imin, const mpz_t imax,
                          Expr* dims, size_t dim_idx) {
    if (dim_idx >= dims->data.function.arg_count)
        return NULL;

    Expr* dim_expr = dims->data.function.args[dim_idx];
    if (dim_expr->type != EXPR_INTEGER || dim_expr->data.integer < 0)
        return NULL;

    size_t n = (size_t)dim_expr->data.integer;
    Expr** elems = malloc(sizeof(Expr*) * n);
    if (!elems) return NULL;

    bool is_last = (dim_idx == dims->data.function.arg_count - 1);

    for (size_t i = 0; i < n; i++) {
        if (is_last) {
            elems[i] = random_integer_range(imin, imax);
            if (!elems[i]) {
                for (size_t j = 0; j < i; j++) expr_free(elems[j]);
                free(elems);
                return NULL;
            }
        } else {
            elems[i] = random_array(imin, imax, dims, dim_idx + 1);
            if (!elems[i]) {
                for (size_t j = 0; j < i; j++) expr_free(elems[j]);
                free(elems);
                return NULL;
            }
        }
    }

    Expr* list = expr_new_function(expr_new_symbol("List"), elems, n);
    free(elems);
    return list;
}

/*
 * builtin_randominteger - implements RandomInteger[...]
 *
 * Forms:
 *   RandomInteger[]              -> 0 or 1
 *   RandomInteger[imax]          -> random in [0, imax]
 *   RandomInteger[{imin, imax}]  -> random in [imin, imax]
 *   RandomInteger[range, n]      -> flat list of n values
 *   RandomInteger[range, {n1, n2, ...}] -> nested array
 */
Expr* builtin_randominteger(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    /* RandomInteger[] -> 0 or 1 */
    if (argc == 0) {
        mpz_t lo, hi;
        mpz_init_set_si(lo, 0);
        mpz_init_set_si(hi, 1);
        Expr* result = random_integer_range(lo, hi);
        mpz_clear(lo);
        mpz_clear(hi);
        return result;
    }

    /* RandomInteger[range] or RandomInteger[range, dims] */
    if (argc == 1 || argc == 2) {
        Expr* range_arg = res->data.function.args[0];

        mpz_t imin, imax;
        mpz_init(imin);
        mpz_init(imax);

        if (!parse_range(range_arg, imin, imax)) {
            mpz_clear(imin);
            mpz_clear(imax);
            return NULL; /* Can't evaluate: symbolic args etc. */
        }

        if (argc == 1) {
            /* Single random integer */
            Expr* result = random_integer_range(imin, imax);
            mpz_clear(imin);
            mpz_clear(imax);
            if (!result) return NULL;
            return result;
        }

        /* argc == 2: RandomInteger[range, n] or RandomInteger[range, {n1, n2, ...}] */
        Expr* dim_arg = res->data.function.args[1];

        if (dim_arg->type == EXPR_INTEGER && dim_arg->data.integer >= 0) {
            /* Flat list of n values */
            size_t n = (size_t)dim_arg->data.integer;
            Expr** elems = malloc(sizeof(Expr*) * n);
            if (!elems) {
                mpz_clear(imin);
                mpz_clear(imax);
                return NULL;
            }
            for (size_t i = 0; i < n; i++) {
                elems[i] = random_integer_range(imin, imax);
                if (!elems[i]) {
                    for (size_t j = 0; j < i; j++) expr_free(elems[j]);
                    free(elems);
                    mpz_clear(imin);
                    mpz_clear(imax);
                    return NULL;
                }
            }
            Expr* list = expr_new_function(expr_new_symbol("List"), elems, n);
            free(elems);
            mpz_clear(imin);
            mpz_clear(imax);
            return list;
        }

        if (dim_arg->type == EXPR_FUNCTION &&
            dim_arg->data.function.head->type == EXPR_SYMBOL &&
            dim_arg->data.function.head->data.symbol == SYM_List) {
            /* Multi-dimensional array */
            /* Validate all dimensions are non-negative integers */
            for (size_t i = 0; i < dim_arg->data.function.arg_count; i++) {
                Expr* d = dim_arg->data.function.args[i];
                if (d->type != EXPR_INTEGER || d->data.integer < 0) {
                    mpz_clear(imin);
                    mpz_clear(imax);
                    return NULL;
                }
            }
            Expr* result = random_array(imin, imax, dim_arg, 0);
            mpz_clear(imin);
            mpz_clear(imax);
            if (!result) return NULL;
            return result;
        }

        mpz_clear(imin);
        mpz_clear(imax);
        return NULL; /* Unrecognized second argument */
    }

    return NULL; /* Too many arguments */
}

/*
 * Convert a numeric expression to a double value.
 * Supports EXPR_INTEGER, EXPR_REAL, EXPR_BIGINT, and Rational[n,d].
 * Returns true on success, false if the expression is not numeric.
 */
static bool expr_to_real(Expr* e, double* out) {
    if (e->type == EXPR_INTEGER) {
        *out = (double)e->data.integer;
        return true;
    }
    if (e->type == EXPR_REAL) {
        *out = e->data.real;
        return true;
    }
    if (e->type == EXPR_BIGINT) {
        *out = mpz_get_d(e->data.bigint);
        return true;
    }
    int64_t n, d;
    if (is_rational(e, &n, &d)) {
        *out = (double)n / (double)d;
        return true;
    }
    return false;
}

/*
 * Generate a single uniform random real in [0, 1).
 * Uses 53 bits of randomness from GMP for full double precision.
 */
static double random_uniform_01(void) {
    ensure_rand_init();
    mpz_t big;
    mpz_init(big);
    /* 2^53 = 9007199254740992 — enough bits for full double mantissa */
    mpz_t modulus;
    mpz_init(modulus);
    mpz_ui_pow_ui(modulus, 2, 53);
    mpz_urandomm(big, g_rand_state, modulus);
    double val = mpz_get_d(big) / 9007199254740992.0;
    mpz_clear(big);
    mpz_clear(modulus);
    return val;
}

/*
 * Generate a single random real in [xmin, xmax).
 */
static Expr* random_real_range(double xmin, double xmax) {
    double u = random_uniform_01();
    double val = xmin + u * (xmax - xmin);
    return expr_new_real(val);
}

/*
 * Parse a real-valued range argument into xmin, xmax.
 * Supports:
 *   - Numeric atom: range is [0, val]
 *   - List[xmin, xmax]: range is [xmin, xmax]
 * Returns true on success, false if the argument doesn't match.
 */
static bool parse_real_range(Expr* arg, double* xmin, double* xmax) {
    double val;
    if (expr_to_real(arg, &val)) {
        *xmin = 0.0;
        *xmax = val;
        return true;
    }
    /* Check for List[xmin, xmax] */
    if (arg->type == EXPR_FUNCTION &&
        arg->data.function.head->type == EXPR_SYMBOL &&
        arg->data.function.head->data.symbol == SYM_List &&
        arg->data.function.arg_count == 2) {

        Expr* lo = arg->data.function.args[0];
        Expr* hi = arg->data.function.args[1];

        if (!expr_to_real(lo, xmin) || !expr_to_real(hi, xmax))
            return false;

        return true;
    }
    return false;
}

/*
 * Build a multi-dimensional array of random reals.
 * dims is a list expression {n1, n2, ...}, dim_idx is current dimension.
 */
static Expr* random_real_array(double xmin, double xmax,
                               Expr* dims, size_t dim_idx) {
    if (dim_idx >= dims->data.function.arg_count)
        return NULL;

    Expr* dim_expr = dims->data.function.args[dim_idx];
    if (dim_expr->type != EXPR_INTEGER || dim_expr->data.integer < 0)
        return NULL;

    size_t n = (size_t)dim_expr->data.integer;
    Expr** elems = malloc(sizeof(Expr*) * n);
    if (!elems) return NULL;

    bool is_last = (dim_idx == dims->data.function.arg_count - 1);

    for (size_t i = 0; i < n; i++) {
        if (is_last) {
            elems[i] = random_real_range(xmin, xmax);
        } else {
            elems[i] = random_real_array(xmin, xmax, dims, dim_idx + 1);
            if (!elems[i]) {
                for (size_t j = 0; j < i; j++) expr_free(elems[j]);
                free(elems);
                return NULL;
            }
        }
    }

    Expr* list = expr_new_function(expr_new_symbol("List"), elems, n);
    free(elems);
    return list;
}

/*
 * builtin_randomreal - implements RandomReal[...]
 *
 * Forms:
 *   RandomReal[]                   -> random real in [0, 1)
 *   RandomReal[xmax]              -> random real in [0, xmax)
 *   RandomReal[{xmin, xmax}]     -> random real in [xmin, xmax)
 *   RandomReal[range, n]          -> flat list of n values
 *   RandomReal[range, {n1, n2, ...}] -> nested array
 */
Expr* builtin_randomreal(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    /* RandomReal[] -> random real in [0, 1) */
    if (argc == 0) {
        return random_real_range(0.0, 1.0);
    }

    /* RandomReal[range] or RandomReal[range, dims] */
    if (argc == 1 || argc == 2) {
        Expr* range_arg = res->data.function.args[0];

        double xmin, xmax;
        if (!parse_real_range(range_arg, &xmin, &xmax)) {
            return NULL; /* Can't evaluate: symbolic args etc. */
        }

        if (argc == 1) {
            return random_real_range(xmin, xmax);
        }

        /* argc == 2: RandomReal[range, n] or RandomReal[range, {n1, n2, ...}] */
        Expr* dim_arg = res->data.function.args[1];

        if (dim_arg->type == EXPR_INTEGER && dim_arg->data.integer >= 0) {
            /* Flat list of n values */
            size_t n = (size_t)dim_arg->data.integer;
            Expr** elems = malloc(sizeof(Expr*) * n);
            if (!elems) return NULL;
            for (size_t i = 0; i < n; i++) {
                elems[i] = random_real_range(xmin, xmax);
            }
            Expr* list = expr_new_function(expr_new_symbol("List"), elems, n);
            free(elems);
            return list;
        }

        if (dim_arg->type == EXPR_FUNCTION &&
            dim_arg->data.function.head->type == EXPR_SYMBOL &&
            dim_arg->data.function.head->data.symbol == SYM_List) {
            /* Multi-dimensional array */
            for (size_t i = 0; i < dim_arg->data.function.arg_count; i++) {
                Expr* d = dim_arg->data.function.args[i];
                if (d->type != EXPR_INTEGER || d->data.integer < 0) {
                    return NULL;
                }
            }
            Expr* result = random_real_array(xmin, xmax, dim_arg, 0);
            if (!result) return NULL;
            return result;
        }

        return NULL; /* Unrecognized second argument */
    }

    return NULL; /* Too many arguments */
}

/*
 * builtin_seedrandom - implements SeedRandom[n] and SeedRandom[]
 */
Expr* builtin_seedrandom(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    if (argc == 0) {
        /* Re-seed from system entropy */
        ensure_rand_init();
        gmp_randseed_ui(g_rand_state, (unsigned long)time(NULL) ^ (unsigned long)clock());
        return expr_new_symbol("Null");
    }

    if (argc == 1) {
        Expr* arg = res->data.function.args[0];
        if (arg->type == EXPR_INTEGER) {
            ensure_rand_init();
            gmp_randseed_ui(g_rand_state, (unsigned long)arg->data.integer);
            return expr_new_symbol("Null");
        }
        if (arg->type == EXPR_BIGINT) {
            ensure_rand_init();
            gmp_randseed(g_rand_state, arg->data.bigint);
            return expr_new_symbol("Null");
        }
        return NULL; /* Non-integer seed */
    }

    return NULL; /* Too many arguments */
}

/*
 * Generate a single random complex number with real part in [re_min, re_max)
 * and imaginary part in [im_min, im_max).
 * Returns a new Complex[re, im] expression.
 */
static Expr* random_complex_range(double re_min, double re_max,
                                   double im_min, double im_max) {
    double u1 = random_uniform_01();
    double u2 = random_uniform_01();
    double re = re_min + u1 * (re_max - re_min);
    double im = im_min + u2 * (im_max - im_min);
    return make_complex(expr_new_real(re), expr_new_real(im));
}

/*
 * Parse a complex-valued range argument into corner coordinates.
 * A complex range is a rectangle in the complex plane defined by two corners.
 * Supports:
 *   - Complex atom: range is [0+0i, val]
 *   - Real/Integer atom: range is [0+0i, val+0i]
 *   - List[zmin, zmax]: range is the rectangle with corners zmin, zmax
 * Returns true on success, false if the argument doesn't match.
 */
static bool parse_complex_range(Expr* arg,
                                double* re_min, double* re_max,
                                double* im_min, double* im_max) {
    /* Try as a single value (upper corner, lower corner is origin) */
    Expr *re_part, *im_part;
    double val;

    if (is_complex(arg, &re_part, &im_part)) {
        double re_val, im_val;
        if (!expr_to_real(re_part, &re_val) || !expr_to_real(im_part, &im_val))
            return false;
        *re_min = 0.0; *re_max = re_val;
        *im_min = 0.0; *im_max = im_val;
        return true;
    }
    if (expr_to_real(arg, &val)) {
        *re_min = 0.0; *re_max = val;
        *im_min = 0.0; *im_max = 0.0;
        return true;
    }

    /* Check for List[zmin, zmax] */
    if (arg->type == EXPR_FUNCTION &&
        arg->data.function.head->type == EXPR_SYMBOL &&
        arg->data.function.head->data.symbol == SYM_List &&
        arg->data.function.arg_count == 2) {

        Expr* lo = arg->data.function.args[0];
        Expr* hi = arg->data.function.args[1];

        double lo_re = 0.0, lo_im = 0.0, hi_re = 0.0, hi_im = 0.0;

        /* Parse lo */
        Expr *lo_re_part, *lo_im_part;
        if (is_complex(lo, &lo_re_part, &lo_im_part)) {
            if (!expr_to_real(lo_re_part, &lo_re) || !expr_to_real(lo_im_part, &lo_im))
                return false;
        } else if (expr_to_real(lo, &lo_re)) {
            lo_im = 0.0;
        } else {
            return false;
        }

        /* Parse hi */
        Expr *hi_re_part, *hi_im_part;
        if (is_complex(hi, &hi_re_part, &hi_im_part)) {
            if (!expr_to_real(hi_re_part, &hi_re) || !expr_to_real(hi_im_part, &hi_im))
                return false;
        } else if (expr_to_real(hi, &hi_re)) {
            hi_im = 0.0;
        } else {
            return false;
        }

        *re_min = lo_re; *re_max = hi_re;
        *im_min = lo_im; *im_max = hi_im;
        return true;
    }
    return false;
}

/*
 * Build a multi-dimensional array of random complex numbers.
 */
static Expr* random_complex_array(double re_min, double re_max,
                                   double im_min, double im_max,
                                   Expr* dims, size_t dim_idx) {
    if (dim_idx >= dims->data.function.arg_count)
        return NULL;

    Expr* dim_expr = dims->data.function.args[dim_idx];
    if (dim_expr->type != EXPR_INTEGER || dim_expr->data.integer < 0)
        return NULL;

    size_t n = (size_t)dim_expr->data.integer;
    Expr** elems = malloc(sizeof(Expr*) * n);
    if (!elems) return NULL;

    bool is_last = (dim_idx == dims->data.function.arg_count - 1);

    for (size_t i = 0; i < n; i++) {
        if (is_last) {
            elems[i] = random_complex_range(re_min, re_max, im_min, im_max);
        } else {
            elems[i] = random_complex_array(re_min, re_max, im_min, im_max,
                                            dims, dim_idx + 1);
            if (!elems[i]) {
                for (size_t j = 0; j < i; j++) expr_free(elems[j]);
                free(elems);
                return NULL;
            }
        }
    }

    Expr* list = expr_new_function(expr_new_symbol("List"), elems, n);
    free(elems);
    return list;
}

/*
 * builtin_randomcomplex - implements RandomComplex[...]
 *
 * Forms:
 *   RandomComplex[]                       -> random complex in [0,1]+[0,1]i
 *   RandomComplex[zmax]                   -> random complex in rectangle [0, zmax]
 *   RandomComplex[{zmin, zmax}]           -> random complex in rectangle [zmin, zmax]
 *   RandomComplex[range, n]               -> flat list of n values
 *   RandomComplex[range, {n1, n2, ...}]   -> nested array
 */
Expr* builtin_randomcomplex(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    /* RandomComplex[] -> random complex in [0,1]+[0,1]i */
    if (argc == 0) {
        return random_complex_range(0.0, 1.0, 0.0, 1.0);
    }

    /* RandomComplex[range] or RandomComplex[range, dims] */
    if (argc == 1 || argc == 2) {
        Expr* range_arg = res->data.function.args[0];

        double re_min, re_max, im_min, im_max;
        if (!parse_complex_range(range_arg, &re_min, &re_max, &im_min, &im_max)) {
            return NULL; /* Can't evaluate: symbolic args etc. */
        }

        if (argc == 1) {
            return random_complex_range(re_min, re_max, im_min, im_max);
        }

        /* argc == 2 */
        Expr* dim_arg = res->data.function.args[1];

        if (dim_arg->type == EXPR_INTEGER && dim_arg->data.integer >= 0) {
            /* Flat list of n values */
            size_t n = (size_t)dim_arg->data.integer;
            Expr** elems = malloc(sizeof(Expr*) * n);
            if (!elems) return NULL;
            for (size_t i = 0; i < n; i++) {
                elems[i] = random_complex_range(re_min, re_max, im_min, im_max);
            }
            Expr* list = expr_new_function(expr_new_symbol("List"), elems, n);
            free(elems);
            return list;
        }

        if (dim_arg->type == EXPR_FUNCTION &&
            dim_arg->data.function.head->type == EXPR_SYMBOL &&
            dim_arg->data.function.head->data.symbol == SYM_List) {
            /* Multi-dimensional array */
            for (size_t i = 0; i < dim_arg->data.function.arg_count; i++) {
                Expr* d = dim_arg->data.function.args[i];
                if (d->type != EXPR_INTEGER || d->data.integer < 0) {
                    return NULL;
                }
            }
            Expr* result = random_complex_array(re_min, re_max, im_min, im_max,
                                                dim_arg, 0);
            if (!result) return NULL;
            return result;
        }

        return NULL; /* Unrecognized second argument */
    }

    return NULL; /* Too many arguments */
}

/*
 * Pick a random index in [0, n) using uniform distribution.
 */
static size_t random_index(size_t n) {
    ensure_rand_init();
    mpz_t modulus, result;
    mpz_init(modulus);
    mpz_init(result);
    mpz_set_ui(modulus, (unsigned long)n);
    mpz_urandomm(result, g_rand_state, modulus);
    size_t idx = (size_t)mpz_get_ui(result);
    mpz_clear(modulus);
    mpz_clear(result);
    return idx;
}

/*
 * Pick a weighted random index given an array of cumulative weights.
 * cum_weights[i] is the sum of weights[0..i].
 * total is cum_weights[count-1].
 * Uses binary search for efficiency.
 */
static size_t weighted_random_index(double* cum_weights, size_t count, double total) {
    double u = random_uniform_01() * total;
    /* Binary search for the first index where cum_weights[i] > u */
    size_t lo = 0, hi = count;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        if (cum_weights[mid] <= u)
            lo = mid + 1;
        else
            hi = mid;
    }
    /* Clamp to valid range */
    if (lo >= count) lo = count - 1;
    return lo;
}

/*
 * Build a multi-dimensional array of random choices from a list.
 * Uses uniform selection (each element equally likely).
 */
static Expr* random_choice_array(Expr** choices, size_t choice_count,
                                 Expr* dims, size_t dim_idx) {
    if (dim_idx >= dims->data.function.arg_count)
        return NULL;

    Expr* dim_expr = dims->data.function.args[dim_idx];
    if (dim_expr->type != EXPR_INTEGER || dim_expr->data.integer < 0)
        return NULL;

    size_t n = (size_t)dim_expr->data.integer;
    Expr** elems = malloc(sizeof(Expr*) * n);
    if (!elems) return NULL;

    bool is_last = (dim_idx == dims->data.function.arg_count - 1);

    for (size_t i = 0; i < n; i++) {
        if (is_last) {
            size_t idx = random_index(choice_count);
            elems[i] = expr_copy(choices[idx]);
        } else {
            elems[i] = random_choice_array(choices, choice_count,
                                           dims, dim_idx + 1);
            if (!elems[i]) {
                for (size_t j = 0; j < i; j++) expr_free(elems[j]);
                free(elems);
                return NULL;
            }
        }
    }

    Expr* list = expr_new_function(expr_new_symbol("List"), elems, n);
    free(elems);
    return list;
}

/*
 * Build a multi-dimensional array of weighted random choices.
 */
static Expr* weighted_choice_array(Expr** choices, size_t choice_count,
                                   double* cum_weights, double total,
                                   Expr* dims, size_t dim_idx) {
    if (dim_idx >= dims->data.function.arg_count)
        return NULL;

    Expr* dim_expr = dims->data.function.args[dim_idx];
    if (dim_expr->type != EXPR_INTEGER || dim_expr->data.integer < 0)
        return NULL;

    size_t n = (size_t)dim_expr->data.integer;
    Expr** elems = malloc(sizeof(Expr*) * n);
    if (!elems) return NULL;

    bool is_last = (dim_idx == dims->data.function.arg_count - 1);

    for (size_t i = 0; i < n; i++) {
        if (is_last) {
            size_t idx = weighted_random_index(cum_weights, choice_count, total);
            elems[i] = expr_copy(choices[idx]);
        } else {
            elems[i] = weighted_choice_array(choices, choice_count,
                                             cum_weights, total,
                                             dims, dim_idx + 1);
            if (!elems[i]) {
                for (size_t j = 0; j < i; j++) expr_free(elems[j]);
                free(elems);
                return NULL;
            }
        }
    }

    Expr* list = expr_new_function(expr_new_symbol("List"), elems, n);
    free(elems);
    return list;
}

/*
 * Check if an expression is Rule[wlist, elist] (i.e., wlist -> elist).
 * If so, set *wlist and *elist to the left and right sides.
 */
static bool is_rule(Expr* e, Expr** wlist, Expr** elist) {
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_Rule &&
        e->data.function.arg_count == 2) {
        *wlist = e->data.function.args[0];
        *elist = e->data.function.args[1];
        return true;
    }
    return false;
}

/*
 * Check if an expression is a List with at least one element.
 * If so, set *args and *count.
 */
static bool is_nonempty_list(Expr* e, Expr*** args, size_t* count) {
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_List &&
        e->data.function.arg_count > 0) {
        *args = e->data.function.args;
        *count = e->data.function.arg_count;
        return true;
    }
    return false;
}

/*
 * builtin_randomchoice - implements RandomChoice[...]
 *
 * Forms:
 *   RandomChoice[{e1, e2, ...}]                     -> uniform random choice
 *   RandomChoice[list, n]                            -> list of n uniform choices
 *   RandomChoice[list, {n1, n2, ...}]                -> array of uniform choices
 *   RandomChoice[{w1,w2,...}->{e1,e2,...}]            -> weighted random choice
 *   RandomChoice[wlist->elist, n]                    -> list of n weighted choices
 *   RandomChoice[wlist->elist, {n1, n2, ...}]        -> array of weighted choices
 */
Expr* builtin_randomchoice(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    if (argc < 1 || argc > 2) return NULL;

    Expr* first_arg = res->data.function.args[0];
    Expr *wlist_expr, *elist_expr;

    /* Check if first arg is Rule[wlist, elist] (weighted form) */
    if (is_rule(first_arg, &wlist_expr, &elist_expr)) {
        /* Weighted random choice */
        Expr** weights_args;
        size_t weights_count;
        Expr** elems_args;
        size_t elems_count;

        if (!is_nonempty_list(wlist_expr, &weights_args, &weights_count))
            return NULL;
        if (!is_nonempty_list(elist_expr, &elems_args, &elems_count))
            return NULL;
        if (weights_count != elems_count)
            return NULL;

        /* Convert weights to doubles and compute cumulative weights */
        double* cum_weights = malloc(sizeof(double) * weights_count);
        if (!cum_weights) return NULL;

        double total = 0.0;
        bool valid = true;
        for (size_t i = 0; i < weights_count; i++) {
            double w;
            if (!expr_to_real(weights_args[i], &w) || w < 0.0) {
                valid = false;
                break;
            }
            total += w;
            cum_weights[i] = total;
        }

        if (!valid || total <= 0.0) {
            free(cum_weights);
            return NULL;
        }

        if (argc == 1) {
            /* Single weighted choice */
            size_t idx = weighted_random_index(cum_weights, elems_count, total);
            free(cum_weights);
            return expr_copy(elems_args[idx]);
        }

        /* argc == 2 */
        Expr* dim_arg = res->data.function.args[1];

        if (dim_arg->type == EXPR_INTEGER && dim_arg->data.integer >= 0) {
            size_t n = (size_t)dim_arg->data.integer;
            Expr** result_elems = malloc(sizeof(Expr*) * n);
            if (!result_elems) { free(cum_weights); return NULL; }
            for (size_t i = 0; i < n; i++) {
                size_t idx = weighted_random_index(cum_weights, elems_count, total);
                result_elems[i] = expr_copy(elems_args[idx]);
            }
            free(cum_weights);
            Expr* list = expr_new_function(expr_new_symbol("List"), result_elems, n);
            free(result_elems);
            return list;
        }

        if (dim_arg->type == EXPR_FUNCTION &&
            dim_arg->data.function.head->type == EXPR_SYMBOL &&
            dim_arg->data.function.head->data.symbol == SYM_List) {
            for (size_t i = 0; i < dim_arg->data.function.arg_count; i++) {
                Expr* d = dim_arg->data.function.args[i];
                if (d->type != EXPR_INTEGER || d->data.integer < 0) {
                    free(cum_weights);
                    return NULL;
                }
            }
            Expr* result = weighted_choice_array(elems_args, elems_count,
                                                 cum_weights, total,
                                                 dim_arg, 0);
            free(cum_weights);
            if (!result) return NULL;
            return result;
        }

        free(cum_weights);
        return NULL;
    }

    /* Uniform random choice: first arg must be a non-empty list */
    Expr** choices;
    size_t choice_count;
    if (!is_nonempty_list(first_arg, &choices, &choice_count))
        return NULL;

    if (argc == 1) {
        /* Single uniform choice */
        size_t idx = random_index(choice_count);
        return expr_copy(choices[idx]);
    }

    /* argc == 2 */
    Expr* dim_arg = res->data.function.args[1];

    if (dim_arg->type == EXPR_INTEGER && dim_arg->data.integer >= 0) {
        size_t n = (size_t)dim_arg->data.integer;
        Expr** elems = malloc(sizeof(Expr*) * n);
        if (!elems) return NULL;
        for (size_t i = 0; i < n; i++) {
            size_t idx = random_index(choice_count);
            elems[i] = expr_copy(choices[idx]);
        }
        Expr* list = expr_new_function(expr_new_symbol("List"), elems, n);
        free(elems);
        return list;
    }

    if (dim_arg->type == EXPR_FUNCTION &&
        dim_arg->data.function.head->type == EXPR_SYMBOL &&
        dim_arg->data.function.head->data.symbol == SYM_List) {
        for (size_t i = 0; i < dim_arg->data.function.arg_count; i++) {
            Expr* d = dim_arg->data.function.args[i];
            if (d->type != EXPR_INTEGER || d->data.integer < 0) {
                return NULL;
            }
        }
        Expr* result = random_choice_array(choices, choice_count, dim_arg, 0);
        if (!result) return NULL;
        return result;
    }

    return NULL;
}

/*
 * Fisher-Yates partial shuffle: select n elements from indices [0, total)
 * without replacement. Returns an array of n selected indices.
 * Caller must free the returned array.
 */
static size_t* fisher_yates_sample(size_t total, size_t n) {
    size_t* indices = malloc(sizeof(size_t) * total);
    if (!indices) return NULL;
    for (size_t i = 0; i < total; i++) indices[i] = i;

    /* Partial shuffle: only need first n elements */
    for (size_t i = 0; i < n; i++) {
        size_t remaining = total - i;
        size_t j = i + random_index(remaining);
        /* Swap indices[i] and indices[j] */
        size_t tmp = indices[i];
        indices[i] = indices[j];
        indices[j] = tmp;
    }

    /* Return only the first n indices */
    size_t* result = malloc(sizeof(size_t) * n);
    if (!result) { free(indices); return NULL; }
    for (size_t i = 0; i < n; i++) result[i] = indices[i];
    free(indices);
    return result;
}

/*
 * Weighted sampling without replacement.
 * Selects n indices from [0, count) with probabilities proportional to weights.
 * Uses the "remove and renormalize" approach.
 * Returns an array of n selected indices. Caller must free.
 */
static size_t* weighted_sample_without_replacement(double* weights, size_t count, size_t n) {
    /* Make a mutable copy of weights */
    double* w = malloc(sizeof(double) * count);
    if (!w) return NULL;
    for (size_t i = 0; i < count; i++) w[i] = weights[i];

    size_t* result = malloc(sizeof(size_t) * n);
    if (!result) { free(w); return NULL; }

    for (size_t k = 0; k < n; k++) {
        /* Compute total remaining weight */
        double total = 0.0;
        for (size_t i = 0; i < count; i++) total += w[i];
        if (total <= 0.0) { free(w); free(result); return NULL; }

        /* Build cumulative weights for active elements */
        double u = random_uniform_01() * total;
        double cum = 0.0;
        size_t chosen = count - 1; /* fallback */
        for (size_t i = 0; i < count; i++) {
            if (w[i] <= 0.0) continue;
            cum += w[i];
            if (cum > u) {
                chosen = i;
                break;
            }
        }

        result[k] = chosen;
        w[chosen] = 0.0; /* Remove from future selection */
    }

    free(w);
    return result;
}

/*
 * Check if an expression is UpTo[n] and extract n.
 * Returns true if it matches, with *val set to n.
 */
static bool is_upto(Expr* e, int64_t* val) {
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_UpTo &&
        e->data.function.arg_count == 1 &&
        e->data.function.args[0]->type == EXPR_INTEGER) {
        *val = e->data.function.args[0]->data.integer;
        return true;
    }
    return false;
}

/*
 * builtin_randomsample - implements RandomSample[...]
 *
 * Forms:
 *   RandomSample[{e1, e2, ...}]                        -> random permutation
 *   RandomSample[{e1, e2, ...}, n]                     -> sample of n elements (no replacement)
 *   RandomSample[{e1, e2, ...}, UpTo[n]]               -> sample of min(n, len) elements
 *   RandomSample[{w1,...}->{e1,...}]                    -> weighted random permutation
 *   RandomSample[{w1,...}->{e1,...}, n]                 -> weighted sample of n elements
 *   RandomSample[{w1,...}->{e1,...}, UpTo[n]]           -> weighted sample of min(n, len)
 */
Expr* builtin_randomsample(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;

    if (argc < 1 || argc > 2) return NULL;

    Expr* first_arg = res->data.function.args[0];
    Expr *wlist_expr, *elist_expr;

    /* Determine the sample size */
    size_t sample_n;
    bool have_n = false;

    /* Check if first arg is Rule[wlist, elist] (weighted form) */
    if (is_rule(first_arg, &wlist_expr, &elist_expr)) {
        /* Weighted random sample */
        Expr** weights_args;
        size_t weights_count;
        Expr** elems_args;
        size_t elems_count;

        if (!is_nonempty_list(wlist_expr, &weights_args, &weights_count))
            return NULL;
        if (!is_nonempty_list(elist_expr, &elems_args, &elems_count))
            return NULL;
        if (weights_count != elems_count)
            return NULL;

        /* Parse weights to doubles */
        double* weights = malloc(sizeof(double) * weights_count);
        if (!weights) return NULL;

        bool valid = true;
        for (size_t i = 0; i < weights_count; i++) {
            if (!expr_to_real(weights_args[i], &weights[i]) || weights[i] < 0.0) {
                valid = false;
                break;
            }
        }
        if (!valid) { free(weights); return NULL; }

        if (argc == 2) {
            Expr* n_arg = res->data.function.args[1];
            int64_t upto_val;
            if (n_arg->type == EXPR_INTEGER && n_arg->data.integer >= 0) {
                sample_n = (size_t)n_arg->data.integer;
                if (sample_n > elems_count) {
                    free(weights);
                    return NULL; /* Cannot sample more than available */
                }
                have_n = true;
            } else if (is_upto(n_arg, &upto_val) && upto_val >= 0) {
                sample_n = (size_t)upto_val;
                if (sample_n > elems_count) sample_n = elems_count;
                have_n = true;
            } else {
                free(weights);
                return NULL;
            }
        } else {
            /* Full permutation */
            sample_n = elems_count;
            have_n = true;
        }

        size_t* selected = weighted_sample_without_replacement(weights, elems_count, sample_n);
        free(weights);
        if (!selected) return NULL;

        Expr** result_elems = malloc(sizeof(Expr*) * sample_n);
        if (!result_elems) { free(selected); return NULL; }
        for (size_t i = 0; i < sample_n; i++) {
            result_elems[i] = expr_copy(elems_args[selected[i]]);
        }
        free(selected);

        Expr* list = expr_new_function(expr_new_symbol("List"), result_elems, sample_n);
        free(result_elems);
        return list;
    }

    /* Uniform random sample: first arg must be a non-empty list */
    Expr** choices;
    size_t choice_count;
    if (!is_nonempty_list(first_arg, &choices, &choice_count))
        return NULL;

    if (argc == 2) {
        Expr* n_arg = res->data.function.args[1];
        int64_t upto_val;
        if (n_arg->type == EXPR_INTEGER && n_arg->data.integer >= 0) {
            sample_n = (size_t)n_arg->data.integer;
            if (sample_n > choice_count) return NULL; /* Can't sample more than available */
            have_n = true;
        } else if (is_upto(n_arg, &upto_val) && upto_val >= 0) {
            sample_n = (size_t)upto_val;
            if (sample_n > choice_count) sample_n = choice_count;
            have_n = true;
        } else {
            return NULL;
        }
    } else {
        /* Full permutation */
        sample_n = choice_count;
        have_n = true;
    }

    (void)have_n;

    size_t* selected = fisher_yates_sample(choice_count, sample_n);
    if (!selected) return NULL;

    Expr** result_elems = malloc(sizeof(Expr*) * sample_n);
    if (!result_elems) { free(selected); return NULL; }
    for (size_t i = 0; i < sample_n; i++) {
        result_elems[i] = expr_copy(choices[selected[i]]);
    }
    free(selected);

    Expr* list = expr_new_function(expr_new_symbol("List"), result_elems, sample_n);
    free(result_elems);
    return list;
}

void random_init(void) {
    symtab_add_builtin("RandomInteger", builtin_randominteger);
    symtab_get_def("RandomInteger")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("RandomReal", builtin_randomreal);
    symtab_get_def("RandomReal")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("SeedRandom", builtin_seedrandom);
    symtab_get_def("SeedRandom")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("RandomComplex", builtin_randomcomplex);
    symtab_get_def("RandomComplex")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("RandomChoice", builtin_randomchoice);
    symtab_get_def("RandomChoice")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("RandomSample", builtin_randomsample);
    symtab_get_def("RandomSample")->attributes |= ATTR_PROTECTED;
}
