/*
 * bitlength.c
 *
 * BitLength[n] -- the number of binary bits needed to represent the integer n.
 *
 *   n > 0 : Floor[Log[2, n]] + 1   (= mpz_sizeinbase(n, 2)).
 *   n = 0 : 0.
 *   n < 0 : BitLength[BitNot[n]], and BitNot[n] = -n - 1 in two's complement,
 *           so BitLength[-1] = 0, BitLength[-2] = 1, BitLength[-2^k] = k.
 *
 * BitLength is Listable (threads over lists automatically). The packed /
 * NDArray fast path is the int64 kernel ndk_BitLength_ii in src/ndinteger.c and
 * the Compile[] lowering is OP_BLEN_I in src/compile/. All arithmetic here is
 * done in GMP, so machine integers and arbitrary-precision bignums are handled
 * uniformly.
 */

#include "bitwise.h"
#include "print.h"        /* expr_to_string */

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

/* `BitLength::argx: BitLength called with N arguments; 1 argument is
 * expected.`  argx never fires at N == 1, so the count is always plural. */
static Expr* bitlen_emit_argx(size_t argc) {
    fprintf(stderr,
            "BitLength::argx: BitLength called with %zu argument%s; "
            "1 argument is expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

/* `BitLength::int: Integer expected at position <pos> in <call>.`  Matches the
 * IntegerLength::int surface diagnostic. */
static Expr* bitlen_emit_int(size_t pos, Expr* res) {
    char* call_str = expr_to_string(res);
    fprintf(stderr,
            "BitLength::int: Integer expected at position %zu in %s.\n",
            pos, call_str ? call_str : "?");
    free(call_str);
    return NULL;
}

Expr* builtin_bitlength(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) return bitlen_emit_argx(argc);

    Expr* n_expr = res->data.function.args[0];
    if (!expr_is_integer_like(n_expr)) {
        /* Symbolic n flows through silently; a concrete non-integer type
         * (Real, Rational, Complex, ...) triggers the ::int diagnostic and
         * leaves the call unevaluated. */
        if (expr_is_numeric_like(n_expr)) return bitlen_emit_int(1, res);
        return NULL;
    }

    /* m is the magnitude whose base-2 size is the bit length: n itself for
     * n >= 0, and BitNot[n] = -n - 1 for n < 0 (which is >= 0 there).
     * mpz_sizeinbase is exact for base 2 (a power of two). */
    mpz_t n, m;
    expr_to_mpz(n_expr, n);
    mpz_init(m);
    if (mpz_sgn(n) >= 0) {
        mpz_set(m, n);
    } else {
        mpz_add_ui(m, n, 1);   /* n + 1  (<= 0) */
        mpz_neg(m, m);         /* -(n + 1) = -n - 1 = BitNot[n] */
    }
    size_t len = (mpz_sgn(m) == 0) ? 0 : mpz_sizeinbase(m, 2);
    mpz_clear(n);
    mpz_clear(m);

    /* `len` is bounded by the mpz bit-size and lands well within EXPR_INTEGER
     * range on every real target. */
    return expr_new_integer((int64_t)len);
}
