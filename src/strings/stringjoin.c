/*
 * stringjoin.c - StringJoin builtin for Mathilda
 *
 * StringJoin["s1", "s2", ...] - concatenates strings together
 */

#include "picostrings.h"
#include "strings_internal.h"
#include "sym_names.h"
#include "linkedlist.h"   /* a LinkedList of strings */
#include <stdlib.h>
#include <string.h>

/* A borrowed piece of the result: pointer + byte length (measured once). */
typedef struct { const char* s; size_t n; } JoinPiece;
typedef struct { JoinPiece* v; size_t n, cap; size_t total; } JoinVec;

/*
 * Recursively collect the strings of `e`, flattening List wrappers, into jv
 * (borrowed pointers). Returns 1 on success, 0 when a non-string, non-list
 * leaf is found (*bad receives it), -1 on OOM.
 */
static int collect_strings(Expr* e, JoinVec* jv, Expr** bad) {
    if (e->type == EXPR_STRING) {
        if (jv->n == jv->cap) {
            size_t nc = jv->cap ? 2 * jv->cap : 16;
            JoinPiece* nv = realloc(jv->v, nc * sizeof(JoinPiece));
            if (!nv) return -1;
            jv->v = nv;
            jv->cap = nc;
        }
        size_t n = strlen(e->data.string);
        jv->v[jv->n].s = e->data.string;
        jv->v[jv->n].n = n;
        jv->n++;
        jv->total += n;
        return 1;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_List) {
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            int rc = collect_strings(e->data.function.args[i], jv, bad);
            if (rc != 1) return rc;
        }
        return 1;
    }
    /* A LinkedList flattens like a List, its strings borrowed from the boxed
     * slots (no code runs during the join, so the slots stay put). A typed
     * list holds numbers, never strings. */
    const LinkedList* ll = expr_linkedlist(e);
    if (ll && (ll_kind(ll) == LL_BOXED || ll_length(ll) == 0)) {
        for (int64_t i = 0, n = ll_length(ll); i < n; i++) {
            int rc = collect_strings((Expr*)ll_peek_boxed(ll, i), jv, bad);
            if (rc != 1) return rc;
        }
        return 1;
    }
    *bad = e;
    return 0;
}

/*
 * builtin_stringjoin:
 * StringJoin[s1, s2, ...] concatenates all string arguments into a single
 * string. Any List wrappers in the arguments are recursively flattened so
 * that all enclosed strings participate in the join. One pass measures every
 * piece, a second copies it into a single buffer the result adopts.
 *
 * StringJoin[] with no arguments returns the empty string "".
 * A non-string, non-list leaf emits StringJoin::string (naming the top-level
 * argument that contains it) and leaves the call unevaluated.
 */
Expr* builtin_stringjoin(Expr* res) {
    if (res->type != EXPR_FUNCTION)
        return NULL;

    size_t argc = res->data.function.arg_count;
    if (argc == 0)
        return expr_new_string("");

    JoinVec jv = { NULL, 0, 0, 0 };
    for (size_t i = 0; i < argc; i++) {
        Expr* bad = NULL;
        int rc = collect_strings(res->data.function.args[i], &jv, &bad);
        if (rc != 1) {
            free(jv.v);
            if (rc == 0)
                return strings_arg_message("StringJoin", "string", "String", i + 1, res);
            return NULL;
        }
    }

    char* buf = malloc(jv.total + 1);
    if (!buf) { free(jv.v); return NULL; }
    size_t off = 0;
    for (size_t i = 0; i < jv.n; i++) {
        memcpy(buf + off, jv.v[i].s, jv.v[i].n);
        off += jv.v[i].n;
    }
    buf[off] = '\0';
    free(jv.v);
    return expr_new_string_adopt(buf);
}
