/*
 * stringrotate.c - StringRotateLeft / StringRotateRight.
 *
 *   StringRotateLeft[s, n]   cycles the characters of s n places to the left
 *   StringRotateLeft[s]      one place
 *   StringRotateRight[s, n]  n places to the right (StringRotateLeft[s, -n])
 *
 * n may be any integer (negative rotates the other way; |n| > length wraps).
 * Positions count characters (Unicode code points). The result is two memcpys
 * of the subject's bytes.
 */

#include "picostrings.h"
#include "strings_internal.h"
#include "utf8.h"

#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <gmp.h>

static Expr* rotate(Expr* res, const char* head, int left) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;
    Expr* s = res->data.function.args[0];
    if (s->type != EXPR_STRING)
        return strings_arg_message(head, "string", "String", 1, res);
    Expr* nx = (argc == 2) ? res->data.function.args[1] : NULL;
    if (nx && !expr_is_integer_like(nx))
        return strings_arg_message(head, "int", "Integer", 2, res);

    const char* str = s->data.string;
    size_t len = strlen(str);
    size_t clen = u8_length(str, len);
    if (clen == 0) return expr_new_string("");

    /* Left rotation by k in [0, clen): the result starts at character k.
     * n is reduced mod the length in GMP, so a BigInt shift (10^20) and
     * n = -2^63 (whose negation overflows int64) both rotate correctly. */
    int64_t m = (int64_t)clen;
    int64_t k;
    if (!nx) {
        k = left ? 1 % m : (m - 1 % m) % m;
    } else {
        mpz_t z;
        mpz_init(z);
        expr_to_mpz(nx, z);
        if (!left) mpz_neg(z, z);
        k = (int64_t)mpz_fdiv_ui(z, (unsigned long)m);   /* in [0, m) */
        mpz_clear(z);
    }
    size_t cut = u8_offset(str, len, (size_t)k);

    char* buf = malloc(len + 1);
    if (!buf) return NULL;
    memcpy(buf, str + cut, len - cut);
    memcpy(buf + (len - cut), str, cut);
    buf[len] = '\0';
    return expr_new_string_adopt(buf);
}

Expr* builtin_stringrotateleft(Expr* res)  { return rotate(res, "StringRotateLeft", 1); }
Expr* builtin_stringrotateright(Expr* res) { return rotate(res, "StringRotateRight", 0); }
