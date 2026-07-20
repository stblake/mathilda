/*
 * stringextract.c - StringExtract[subject, spec_1, ..., spec_k]
 *
 * Splits a string into blocks and selects blocks by position. Each spec_j is a
 * level: either a bare position (n, -n, {n1,...}, n1;;n2, All) or a rule
 * `sep -> pos`. Level j splits its input on separator sep_j, then selects with
 * pos_j; if more levels remain it recurses into each selected block.
 *
 * The split is delegated to StringSplit (built here as StringSplit[str, sep] and
 * evaluated), so the entire WL string-pattern engine, whitespace-run collapsing,
 * and empty-end trimming are reused verbatim. This makes the documented
 * equivalences exact:
 *     StringExtract[s, patt -> All]  ==  StringSplit[s, patt]
 *     StringExtract[s, {p1, p2, ...}] == Part[StringSplit[s], {p1, p2, ...}]
 *
 * Bare positions get depth-default separators: the last (lowest) level splits on
 * whitespace, the level above on a single "\n", the next on "\n\n", and so on.
 * An out-of-range single index yields Missing["PartAbsent", n] (matching Wolfram
 * and the Missing["KeyAbsent", ...] convention in part.c), not the usual
 * NULL-unevaluated result.
 *
 * Registered by regex_init() (regex_init.c) alongside its StringSplit dependency.
 * Docstring lives in info.c. Attribute: Protected.
 */

#include "picostrings.h"
#include "sym_names.h"
#include "eval.h"
#include "part.h"
#include "common.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* True if e is a List[...] expression. */
static int is_list_expr(const Expr* e) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_List;
}

/* Missing["PartAbsent", n]. */
static Expr* make_missing(int64_t n) {
    Expr* margs[2];
    margs[0] = expr_new_string("PartAbsent");
    margs[1] = expr_new_integer(n);
    return expr_new_function(expr_new_symbol(SYM_Missing), margs, 2);
}

/*
 * Split `str` (an EXPR_STRING) into a List of block strings using separator
 * `sep` (a borrowed pattern Expr: a Whitespace symbol, a "\n"* string, or any
 * string pattern). Returns a fresh List, or NULL if StringSplit did not yield a
 * List (e.g. the regex engine is unavailable) -> the caller bails unevaluated.
 */
static Expr* do_split(Expr* str, Expr* sep) {
    Expr* args[2];
    args[0] = expr_copy(str);
    args[1] = expr_copy(sep);
    Expr* call = expr_new_function(expr_new_symbol(SYM_StringSplit), args, 2);
    Expr* blocks = eval_and_free(call);
    if (is_list_expr(blocks)) return blocks;
    expr_free(blocks);
    return NULL;
}

/*
 * Resolve a single integer position (1-based; negative counts from the end)
 * against `blocks` (a List). Returns a fresh copy of the selected block, or
 * Missing["PartAbsent", n] if out of range.
 */
static Expr* pos_single(Expr* blocks, int64_t n) {
    int64_t N = (int64_t)blocks->data.function.arg_count;
    int64_t k = (n < 0) ? N + n + 1 : n;
    if (k < 1 || k > N) return make_missing(n);
    return expr_copy(blocks->data.function.args[k - 1]);
}

/*
 * Apply position spec `pos` to `blocks`. On success returns a fresh Expr and
 * sets *is_multi (1 if the result is a List of blocks, 0 if a single
 * block/Missing). Returns NULL for an unsupported spec (caller bails).
 */
static Expr* apply_position(Expr* blocks, Expr* pos, int* is_multi) {
    /* Single integer: nth / -nth. */
    if (pos->type == EXPR_INTEGER) {
        *is_multi = 0;
        return pos_single(blocks, pos->data.integer);
    }
    /* All -> every block (== StringSplit). */
    if (pos->type == EXPR_SYMBOL && pos->data.symbol.name == SYM_All) {
        *is_multi = 1;
        return expr_copy(blocks);
    }
    /* {n1, n2, ...} -> collection of blocks. */
    if (is_list_expr(pos)) {
        *is_multi = 1;
        size_t m = pos->data.function.arg_count;
        Expr** out = malloc(sizeof(Expr*) * (m ? m : 1));
        for (size_t i = 0; i < m; i++) {
            Expr* pi = pos->data.function.args[i];
            if (pi->type == EXPR_INTEGER) {
                out[i] = pos_single(blocks, pi->data.integer);
            } else {
                /* Span/All/etc. nested in the list: defer to Part semantics. */
                Expr* v = expr_part(blocks, &pi, 1);
                out[i] = v ? v : expr_copy(pi);
            }
        }
        Expr* r = expr_new_function(expr_new_symbol(SYM_List), out, m);
        free(out);
        return r;
    }
    /* n1;;n2 (Span) -> range of blocks (Part clips to valid range). */
    if (head_is(pos, SYM_Span)) {
        *is_multi = 1;
        return expr_part(blocks, &pos, 1);   /* NULL -> bail */
    }
    return NULL;
}

/* One level's separator + position, borrowed or owned (see builtin). */
typedef struct { Expr* sep; Expr* pos; } LevelSpec;

/*
 * Recursively extract from a single string `str` at level `level` (0-based) of
 * `k` levels. Returns a fresh Expr (string, Missing, or nested List), or NULL to
 * bail (leave the whole StringExtract unevaluated).
 */
static Expr* extract(Expr* str, LevelSpec* specs, int k, int level) {
    Expr* blocks = do_split(str, specs[level].sep);
    if (!blocks) return NULL;

    int is_multi = 0;
    Expr* sel = apply_position(blocks, specs[level].pos, &is_multi);
    expr_free(blocks);
    if (!sel) return NULL;

    if (level == k - 1) return sel;   /* last (lowest) level: done */

    /* Descend into the next level. */
    if (!is_multi) {
        if (head_is(sel, SYM_Missing)) return sel;   /* absent block propagates */
        Expr* r = extract(sel, specs, k, level + 1);
        expr_free(sel);
        return r;
    }

    /* sel is a List: map the next level over its string elements. */
    size_t m = sel->data.function.arg_count;
    Expr** out = malloc(sizeof(Expr*) * (m ? m : 1));
    for (size_t i = 0; i < m; i++) {
        Expr* e = sel->data.function.args[i];
        if (e->type == EXPR_STRING) {
            out[i] = extract(e, specs, k, level + 1);
            if (!out[i]) {                       /* bail: free partials */
                for (size_t j = 0; j < i; j++) expr_free(out[j]);
                free(out);
                expr_free(sel);
                return NULL;
            }
        } else {
            out[i] = expr_copy(e);               /* Missing / non-string kept */
        }
    }
    Expr* r = expr_new_function(expr_new_symbol(SYM_List), out, m);
    free(out);
    expr_free(sel);
    return r;
}

Expr* builtin_stringextract(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** a = res->data.function.args;

    if (argc < 2) return builtin_arg_error("StringExtract", argc, 2, SIZE_MAX);

    Expr* subject = a[0];
    int k = (int)(argc - 1);   /* number of levels */

    /* Build the per-level separator/position table. Explicit `sep -> pos` rules
     * carry their own (borrowed) separator; bare positions get depth-default
     * separators (whitespace at the lowest level, growing runs of "\n" above). */
    LevelSpec* specs = malloc(sizeof(LevelSpec) * (size_t)k);
    Expr** owned = malloc(sizeof(Expr*) * (size_t)k);   /* owned default seps */
    int n_owned = 0;
    for (int j = 0; j < k; j++) {
        Expr* arg = a[1 + j];
        if ((head_is(arg, SYM_Rule) || head_is(arg, SYM_RuleDelayed)) &&
            arg->data.function.arg_count == 2) {
            specs[j].sep = arg->data.function.args[0];   /* borrowed */
            specs[j].pos = arg->data.function.args[1];   /* borrowed */
        } else {
            specs[j].pos = arg;                          /* borrowed */
            int nl = (k - 1) - j;                        /* newlines for this level */
            Expr* sep;
            if (nl == 0) {
                sep = expr_new_symbol(SYM_Whitespace);
            } else {
                char* buf = malloc((size_t)nl + 1);
                for (int t = 0; t < nl; t++) buf[t] = '\n';
                buf[nl] = '\0';
                sep = expr_new_string(buf);
                free(buf);
            }
            specs[j].sep = sep;
            owned[n_owned++] = sep;
        }
    }

    /* Apply. Thread over a list subject; a single string extracts directly. */
    Expr* result;
    if (is_list_expr(subject)) {
        size_t m = subject->data.function.arg_count;
        Expr** out = malloc(sizeof(Expr*) * (m ? m : 1));
        int ok = 1;
        size_t i;
        for (i = 0; i < m; i++) {
            Expr* si = subject->data.function.args[i];
            if (si->type == EXPR_STRING) {
                out[i] = extract(si, specs, k, 0);
                if (!out[i]) ok = 0;
            } else {
                out[i] = NULL;
                ok = 0;
            }
            if (!ok) { for (size_t j = 0; j < i; j++) expr_free(out[j]); break; }
        }
        result = ok ? expr_new_function(expr_new_symbol(SYM_List), out, m) : NULL;
        free(out);
    } else if (subject->type == EXPR_STRING) {
        result = extract(subject, specs, k, 0);
    } else {
        result = NULL;   /* non-string subject: leave unevaluated */
    }

    for (int i = 0; i < n_owned; i++) expr_free(owned[i]);
    free(owned);
    free(specs);
    return result;
}
