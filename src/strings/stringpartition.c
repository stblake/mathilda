/*
 * stringpartition.c - StringPartition builtin for Mathilda
 *
 * StringPartition["string", n]        - non-overlapping length-n blocks
 * StringPartition["string", n, d]     - length-n blocks starting every d chars
 * StringPartition["string", UpTo[n]]  - length-<=n blocks; final may be shorter
 * StringPartition[{s1, s2, ...}, spec] - threads over a list of strings
 *
 * Strings are treated as raw byte arrays (consistent with StringTake/StringDrop
 * across this subsystem); no UTF-8 codepoint decoding is performed.
 */

#include "picostrings.h"
#include "strings_internal.h"
#include "common.h"
#include "eval.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>

/*
 * sp_block:
 * Builds a new EXPR_STRING from the byte range str[start .. start+n) (0-based,
 * caller guarantees start + n <= strlen(str)). Returns NULL on allocation
 * failure.
 */
static Expr* sp_block(const char* str, int64_t start, int64_t n) {
    char* buf = malloc((size_t)n + 1);
    if (!buf) return NULL;
    memcpy(buf, str + start, (size_t)n);
    buf[n] = '\0';
    Expr* result = expr_new_string(buf);
    free(buf);
    return result;
}

/*
 * builtin_stringpartition:
 * See file header for the accepted forms. Returns NULL (unevaluated) for
 * non-string / non-integer arguments and for non-positive n or d.
 */
Expr* builtin_stringpartition(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;

    if (argc < 2 || argc > 3)
        return builtin_arg_error("StringPartition", argc, 2, 3);

    Expr* arg0 = args[0];

    /* StringPartition[{s1, s2, ...}, spec, ...] - thread over list of strings */
    if (arg0->type == EXPR_FUNCTION &&
        arg0->data.function.head->type == EXPR_SYMBOL &&
        arg0->data.function.head->data.symbol.name == SYM_List) {
        size_t m = arg0->data.function.arg_count;
        Expr** results = malloc(sizeof(Expr*) * (m ? m : 1));
        if (!results) return NULL;
        for (size_t i = 0; i < m; i++) {
            Expr** inner = malloc(sizeof(Expr*) * argc);
            if (!inner) {
                for (size_t j = 0; j < i; j++) expr_free(results[j]);
                free(results);
                return NULL;
            }
            inner[0] = expr_copy(arg0->data.function.args[i]);
            for (size_t k = 1; k < argc; k++) inner[k] = expr_copy(args[k]);
            Expr* call = expr_new_function(expr_new_symbol(SYM_StringPartition),
                                           inner, argc);
            free(inner);
            results[i] = evaluate(call);
            expr_free(call);
        }
        Expr* out = expr_new_function(expr_new_symbol(SYM_List), results, m);
        free(results);
        return out;
    }

    if (arg0->type != EXPR_STRING) return NULL;

    const char* str = arg0->data.string;
    int64_t len = (int64_t)strlen(str);

    /* Parse block length n and whether a short final block is allowed. */
    int64_t n;
    bool allow_short;
    if (args[1]->type == EXPR_INTEGER) {
        n = args[1]->data.integer;
        allow_short = false;
    } else if (strings_is_upto(args[1], &n)) {
        allow_short = true;
    } else {
        return NULL;
    }

    /* Offset d defaults to n; an explicit third argument overrides it. */
    int64_t d;
    if (argc == 3) {
        if (args[2]->type != EXPR_INTEGER) return NULL;
        d = args[2]->data.integer;
    } else {
        d = n;
    }

    if (n <= 0 || d <= 0) return NULL;

    /* Collect blocks. Upper bound on count: one per offset step, plus a
     * possible short tail. */
    int64_t max_blocks = len / d + 2;
    Expr** blocks = malloc(sizeof(Expr*) * (size_t)max_blocks);
    if (!blocks) return NULL;

    size_t count = 0;
    for (int64_t start = 0; start < len; start += d) {
        int64_t end = start + n;
        Expr* block;
        if (end <= len) {
            block = sp_block(str, start, n);
        } else if (allow_short) {
            block = sp_block(str, start, len - start);
        } else {
            break;
        }
        if (!block) {
            for (size_t j = 0; j < count; j++) expr_free(blocks[j]);
            free(blocks);
            return NULL;
        }
        blocks[count++] = block;
        if (end > len) break;  /* short tail emitted; nothing further fits */
    }

    Expr* result = expr_new_function(expr_new_symbol(SYM_List), blocks, count);
    free(blocks);
    return result;
}
