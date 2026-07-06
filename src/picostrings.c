/*
 * strings.c - String manipulation builtins for Mathilda
 *
 * StringLength["string"]        - gives the number of characters in a string
 * Characters["string"]          - gives a list of the individual characters in a string
 * StringJoin["s1", "s2", ...]   - concatenates strings together
 * StringPart["string", spec]    - extracts characters by index, list, or Span
 * StringTake["string", spec]    - extracts substrings by count, range, or step
 */

#include "picostrings.h"
#include "symtab.h"
#include "attr.h"
#include "eval.h"
#include "sym_names.h"
#include <string.h>
#include <stdio.h>

/*
 * builtin_stringlength:
 * StringLength["string"] returns the number of characters in the string
 * as an integer. Returns NULL (unevaluated) for non-string arguments.
 */
Expr* builtin_stringlength(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1)
        return NULL;

    Expr* arg = res->data.function.args[0];
    if (arg->type != EXPR_STRING)
        return NULL;

    int64_t len = (int64_t)strlen(arg->data.string);
    Expr* result = expr_new_integer(len);

    return result;
}

/*
 * builtin_characters:
 * Characters["string"] returns a List of single-character strings.
 * Each character in the input string becomes a length-1 string element.
 * Returns NULL (unevaluated) for non-string arguments.
 */
Expr* builtin_characters(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1)
        return NULL;

    Expr* arg = res->data.function.args[0];
    if (arg->type != EXPR_STRING)
        return NULL;

    const char* str = arg->data.string;
    size_t len = strlen(str);

    Expr** chars = malloc(sizeof(Expr*) * len);
    if (!chars) {
    
        return expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    }

    char buf[2];
    buf[1] = '\0';
    for (size_t i = 0; i < len; i++) {
        buf[0] = str[i];
        chars[i] = expr_new_string(buf);
    }

    Expr* result = expr_new_function(expr_new_symbol(SYM_List), chars, len);
    free(chars);

    return result;
}

/*
 * collect_strings:
 * Recursively collects all strings from an expression, flattening any List
 * wrappers encountered. Appends string pointers (borrowed, not copied) to
 * the output array. Returns true if all leaves are strings, false otherwise.
 */
static bool collect_strings(Expr* e, const char*** strs, size_t* count, size_t* cap) {
    if (e->type == EXPR_STRING) {
        if (*count >= *cap) {
            *cap *= 2;
            *strs = realloc(*strs, sizeof(const char*) * (*cap));
        }
        (*strs)[(*count)++] = e->data.string;
        return true;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_List) {
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            if (!collect_strings(e->data.function.args[i], strs, count, cap))
                return false;
        }
        return true;
    }
    return false;
}

/*
 * builtin_stringjoin:
 * StringJoin[s1, s2, ...] concatenates all string arguments into a single
 * string. Any List wrappers in the arguments are recursively flattened so
 * that all enclosed strings participate in the join.
 *
 * StringJoin[] with no arguments returns the empty string "".
 * Returns NULL (unevaluated) if any non-string, non-list leaf is found.
 */
Expr* builtin_stringjoin(Expr* res) {
    if (res->type != EXPR_FUNCTION)
        return NULL;

    size_t argc = res->data.function.arg_count;

    /* StringJoin[] -> "" */
    if (argc == 0)
        return expr_new_string("");

    /* Collect all strings, flattening lists */
    size_t cap = 16, count = 0;
    const char** strs = malloc(sizeof(const char*) * cap);

    for (size_t i = 0; i < argc; i++) {
        if (!collect_strings(res->data.function.args[i], &strs, &count, &cap)) {
            free(strs);
            return NULL;
        }
    }

    /* Compute total length */
    size_t total_len = 0;
    for (size_t i = 0; i < count; i++) {
        total_len += strlen(strs[i]);
    }

    /* Build concatenated string */
    char* buf = malloc(total_len + 1);
    buf[0] = '\0';
    size_t offset = 0;
    for (size_t i = 0; i < count; i++) {
        size_t slen = strlen(strs[i]);
        memcpy(buf + offset, strs[i], slen);
        offset += slen;
    }
    buf[total_len] = '\0';

    free(strs);

    Expr* result = expr_new_string(buf);
    free(buf);
    return result;
}

/*
 * stringpart_single:
 * Helper to extract a single character from a string by 1-based index.
 * Negative indices count from the end. Returns a new single-character
 * string Expr*, or NULL if the index is out of bounds or not an integer.
 */
static Expr* stringpart_single(const char* str, int64_t len, Expr* idx) {
    if (idx->type != EXPR_INTEGER) return NULL;
    int64_t k = idx->data.integer;
    if (k < 0) k = len + k + 1;
    if (k < 1 || k > len) return NULL;
    char buf[2] = { str[k - 1], '\0' };
    return expr_new_string(buf);
}

/*
 * builtin_stringpart:
 * StringPart["string", n]        - gives the nth character as a string
 * StringPart["string", {n1,n2,...}] - gives a list of characters
 * StringPart["string", m;;n]     - gives characters m through n
 * StringPart["string", m;;n;;s]  - gives characters m through n in steps of s
 * StringPart[{s1,s2,...}, spec]   - gives the list of results for each si
 *
 * Negative indices count from the end. Returns NULL (unevaluated) for
 * invalid arguments.
 */
Expr* builtin_stringpart(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)
        return NULL;

    Expr* arg0 = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];

    /* StringPart[{s1, s2, ...}, spec] - map over list of strings */
    if (arg0->type == EXPR_FUNCTION &&
        arg0->data.function.head->type == EXPR_SYMBOL &&
        arg0->data.function.head->data.symbol == SYM_List) {
        size_t n = arg0->data.function.arg_count;
        Expr** results = malloc(sizeof(Expr*) * n);
        if (!results) return NULL;
        for (size_t i = 0; i < n; i++) {
            /* Build StringPart[si, spec] and evaluate */
            Expr** inner_args = malloc(sizeof(Expr*) * 2);
            inner_args[0] = expr_copy(arg0->data.function.args[i]);
            inner_args[1] = expr_copy(spec);
            Expr* call = expr_new_function(expr_new_symbol(SYM_StringPart),
                                           inner_args, 2);
            free(inner_args);
            Expr* r = evaluate(call);
            expr_free(call);
            results[i] = r;
        }
        Expr* result = expr_new_function(expr_new_symbol(SYM_List), results, n);
        free(results);
        return result;
    }

    if (arg0->type != EXPR_STRING) return NULL;

    const char* str = arg0->data.string;
    int64_t len = (int64_t)strlen(str);

    /* StringPart["string", n] - single integer index */
    if (spec->type == EXPR_INTEGER) {
        return stringpart_single(str, len, spec);
    }

    /* StringPart["string", {n1, n2, ...}] - list of indices */
    if (spec->type == EXPR_FUNCTION &&
        spec->data.function.head->type == EXPR_SYMBOL &&
        spec->data.function.head->data.symbol == SYM_List) {
        size_t n = spec->data.function.arg_count;
        Expr** results = malloc(sizeof(Expr*) * n);
        if (!results) return NULL;
        for (size_t i = 0; i < n; i++) {
            Expr* ch = stringpart_single(str, len, spec->data.function.args[i]);
            if (!ch) {
                for (size_t j = 0; j < i; j++) expr_free(results[j]);
                free(results);
                return NULL;
            }
            results[i] = ch;
        }
        Expr* result = expr_new_function(expr_new_symbol(SYM_List), results, n);
        free(results);
        return result;
    }

    /* StringPart["string", m;;n] or StringPart["string", m;;n;;s] - Span */
    if (spec->type == EXPR_FUNCTION &&
        spec->data.function.head->type == EXPR_SYMBOL &&
        spec->data.function.head->data.symbol == SYM_Span) {
        int64_t start = 1, end = len, step = 1;
        size_t span_argc = spec->data.function.arg_count;

        if (span_argc >= 1) {
            Expr* a1 = spec->data.function.args[0];
            if (a1->type == EXPR_INTEGER) {
                start = a1->data.integer;
                if (start < 0) start = len + start + 1;
            } else if (a1->type == EXPR_SYMBOL &&
                       a1->data.symbol == SYM_All) {
                start = 1;
            } else return NULL;
        }
        if (span_argc >= 2) {
            Expr* a2 = spec->data.function.args[1];
            if (a2->type == EXPR_INTEGER) {
                end = a2->data.integer;
                if (end < 0) end = len + end + 1;
            } else if (a2->type == EXPR_SYMBOL &&
                       a2->data.symbol == SYM_All) {
                end = len;
            } else return NULL;
        }
        if (span_argc >= 3) {
            Expr* a3 = spec->data.function.args[2];
            if (a3->type == EXPR_INTEGER) step = a3->data.integer;
            else return NULL;
            if (step == 0) return NULL;
        }

        /* Calculate number of elements */
        int64_t count = 0;
        if (step > 0) {
            if (start >= 1 && end <= len && start <= end)
                count = (end - start) / step + 1;
        } else {
            if (start <= len && end >= 1 && start >= end)
                count = (start - end) / (-step) + 1;
        }
        if (count < 0) count = 0;

        Expr** results = NULL;
        if (count > 0) {
            results = malloc(sizeof(Expr*) * (size_t)count);
            if (!results) return NULL;
        }

        int64_t current = start;
        for (int64_t i = 0; i < count; i++) {
            char buf[2] = { str[current - 1], '\0' };
            results[i] = expr_new_string(buf);
            current += step;
        }

        Expr* result = expr_new_function(expr_new_symbol(SYM_List),
                                         results, (size_t)count);
        free(results);
        return result;
    }

    return NULL;
}

/*
 * is_upto:
 * Returns true if expr is UpTo[n] where n is a positive integer,
 * and writes the integer value to *out.
 */
static bool is_upto(Expr* e, int64_t* out) {
    if (e->type == EXPR_FUNCTION &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol == SYM_UpTo &&
        e->data.function.arg_count == 1 &&
        e->data.function.args[0]->type == EXPR_INTEGER) {
        *out = e->data.function.args[0]->data.integer;
        return true;
    }
    return false;
}

/*
 * stringtake_substring:
 * Helper to extract a substring from str[start-1..end-1] (1-based indices,
 * already resolved to positive). Returns a new EXPR_STRING.
 */
static Expr* stringtake_substring(const char* str, int64_t start, int64_t end) {
    size_t slen = (size_t)(end - start + 1);
    char* buf = malloc(slen + 1);
    if (!buf) return NULL;
    memcpy(buf, str + start - 1, slen);
    buf[slen] = '\0';
    Expr* result = expr_new_string(buf);
    free(buf);
    return result;
}

/*
 * builtin_stringtake:
 * StringTake["string", n]       - first n characters
 * StringTake["string", -n]      - last n characters
 * StringTake["string", UpTo[n]] - first n chars, or as many as available
 * StringTake["string", {n}]     - nth character as a string
 * StringTake["string", {m,n}]   - characters m through n
 * StringTake["string", {m,n,s}] - characters m through n in steps of s
 * StringTake[{s1,s2,...}, spec]  - maps over list of strings
 *
 * Returns NULL (unevaluated) for invalid arguments.
 */
Expr* builtin_stringtake(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)
        return NULL;

    Expr* arg0 = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];

    /* StringTake[{s1, s2, ...}, spec] - map over list of strings */
    if (arg0->type == EXPR_FUNCTION &&
        arg0->data.function.head->type == EXPR_SYMBOL &&
        arg0->data.function.head->data.symbol == SYM_List) {
        size_t n = arg0->data.function.arg_count;
        Expr** results = malloc(sizeof(Expr*) * n);
        if (!results) return NULL;
        for (size_t i = 0; i < n; i++) {
            Expr** inner_args = malloc(sizeof(Expr*) * 2);
            inner_args[0] = expr_copy(arg0->data.function.args[i]);
            inner_args[1] = expr_copy(spec);
            Expr* call = expr_new_function(expr_new_symbol(SYM_StringTake),
                                           inner_args, 2);
            free(inner_args);
            Expr* r = evaluate(call);
            expr_free(call);
            results[i] = r;
        }
        Expr* result = expr_new_function(expr_new_symbol(SYM_List), results, n);
        free(results);
        return result;
    }

    if (arg0->type != EXPR_STRING) return NULL;

    const char* str = arg0->data.string;
    int64_t len = (int64_t)strlen(str);

    /* StringTake["string", n] - first n chars (n>0) or last |n| chars (n<0) */
    if (spec->type == EXPR_INTEGER) {
        int64_t n = spec->data.integer;
        if (n == 0) return expr_new_string("");
        if (n > 0) {
            if (n > len) return NULL;
            return stringtake_substring(str, 1, n);
        } else {
            int64_t take = -n;
            if (take > len) return NULL;
            return stringtake_substring(str, len - take + 1, len);
        }
    }

    /* StringTake["string", UpTo[n]] */
    int64_t upto_n;
    if (is_upto(spec, &upto_n)) {
        if (upto_n <= 0) return NULL;
        int64_t take = upto_n < len ? upto_n : len;
        if (take == 0) return expr_new_string("");
        return stringtake_substring(str, 1, take);
    }

    /* StringTake["string", {spec}] - list spec */
    if (spec->type == EXPR_FUNCTION &&
        spec->data.function.head->type == EXPR_SYMBOL &&
        spec->data.function.head->data.symbol == SYM_List) {
        size_t spec_argc = spec->data.function.arg_count;

        /* StringTake["string", {n}] - single character */
        if (spec_argc == 1) {
            Expr* a1 = spec->data.function.args[0];
            if (a1->type != EXPR_INTEGER) return NULL;
            int64_t k = a1->data.integer;
            if (k < 0) k = len + k + 1;
            if (k < 1 || k > len) return NULL;
            char buf[2] = { str[k - 1], '\0' };
            return expr_new_string(buf);
        }

        /* StringTake["string", {m, n}] - range */
        if (spec_argc == 2) {
            Expr* a1 = spec->data.function.args[0];
            Expr* a2 = spec->data.function.args[1];
            if (a1->type != EXPR_INTEGER || a2->type != EXPR_INTEGER)
                return NULL;
            int64_t m = a1->data.integer;
            int64_t n = a2->data.integer;
            if (m < 0) m = len + m + 1;
            if (n < 0) n = len + n + 1;
            if (m < 1 || m > len || n < 1 || n > len || m > n) return NULL;
            return stringtake_substring(str, m, n);
        }

        /* StringTake["string", {m, n, s}] - range with step */
        if (spec_argc == 3) {
            Expr* a1 = spec->data.function.args[0];
            Expr* a2 = spec->data.function.args[1];
            Expr* a3 = spec->data.function.args[2];
            if (a1->type != EXPR_INTEGER || a2->type != EXPR_INTEGER ||
                a3->type != EXPR_INTEGER)
                return NULL;
            int64_t m = a1->data.integer;
            int64_t n = a2->data.integer;
            int64_t s = a3->data.integer;
            if (s == 0) return NULL;
            if (m < 0) m = len + m + 1;
            if (n < 0) n = len + n + 1;

            /* Calculate count */
            int64_t count = 0;
            if (s > 0) {
                if (m >= 1 && n <= len && m <= n)
                    count = (n - m) / s + 1;
            } else {
                if (m <= len && n >= 1 && m >= n)
                    count = (m - n) / (-s) + 1;
            }
            if (count <= 0) return expr_new_string("");

            char* buf = malloc((size_t)count + 1);
            if (!buf) return NULL;
            int64_t current = m;
            for (int64_t i = 0; i < count; i++) {
                buf[i] = str[current - 1];
                current += s;
            }
            buf[count] = '\0';
            Expr* result = expr_new_string(buf);
            free(buf);
            return result;
        }
    }

    return NULL;
}

/*
 * sd_emit_argrx:
 * Emit `StringDrop::argrx: StringDrop called with N arguments; 2 arguments
 * are expected.` to stderr and return NULL so the call is left unevaluated.
 * Matches Mathematica's diagnostic for a wrong argument count.
 */
static Expr* sd_emit_argrx(size_t argc) {
    fprintf(stderr,
            "StringDrop::argrx: StringDrop called with %zu argument%s; "
            "2 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

/*
 * stringdrop_build_kept:
 * Rebuild a string from str keeping only the positions i (0-based) for which
 * keep[i] is true, in order. Returns a new EXPR_STRING (possibly ""). The
 * complementary "keep" mask is how StringDrop is expressed: the sequence spec
 * marks characters to drop, and everything else survives.
 */
static Expr* stringdrop_build_kept(const char* str, int64_t len, const bool* keep) {
    char* buf = malloc((size_t)len + 1);
    if (!buf) return NULL;
    int64_t j = 0;
    for (int64_t i = 0; i < len; i++) {
        if (keep[i]) buf[j++] = str[i];
    }
    buf[j] = '\0';
    Expr* result = expr_new_string(buf);
    free(buf);
    return result;
}

/*
 * builtin_stringdrop:
 * StringDrop["string", n]       - drop the first n characters
 * StringDrop["string", -n]      - drop the last n characters
 * StringDrop["string", UpTo[n]] - drop n characters, or as many as available
 * StringDrop["string", {n}]     - drop the nth character
 * StringDrop["string", {m, n}]  - drop characters m through n
 * StringDrop["string", {m, n, s}] - drop characters m through n in steps of s
 * StringDrop[{s1, s2, ...}, spec] - maps over the list of strings
 *
 * StringDrop is the complement of StringTake: the same standard Wolfram
 * Language sequence specification selects the characters to remove, and the
 * remaining characters are concatenated in order. Negative indices count from
 * the end. A call whose argument count is not two emits StringDrop::argrx and
 * is left unevaluated; an out-of-range position or a symbolic argument leaves
 * the call unevaluated (NULL) so it flows through the evaluator unchanged.
 */
Expr* builtin_stringdrop(Expr* res) {
    if (res->type != EXPR_FUNCTION)
        return NULL;

    size_t argc = res->data.function.arg_count;
    if (argc != 2)
        return sd_emit_argrx(argc);

    Expr* arg0 = res->data.function.args[0];
    Expr* spec = res->data.function.args[1];

    /* StringDrop[{s1, s2, ...}, spec] - map over list of strings */
    if (arg0->type == EXPR_FUNCTION &&
        arg0->data.function.head->type == EXPR_SYMBOL &&
        arg0->data.function.head->data.symbol == SYM_List) {
        size_t n = arg0->data.function.arg_count;
        Expr** results = malloc(sizeof(Expr*) * n);
        if (!results) return NULL;
        for (size_t i = 0; i < n; i++) {
            Expr** inner_args = malloc(sizeof(Expr*) * 2);
            inner_args[0] = expr_copy(arg0->data.function.args[i]);
            inner_args[1] = expr_copy(spec);
            Expr* call = expr_new_function(expr_new_symbol(SYM_StringDrop),
                                           inner_args, 2);
            free(inner_args);
            Expr* r = evaluate(call);
            expr_free(call);
            results[i] = r;
        }
        Expr* result = expr_new_function(expr_new_symbol(SYM_List), results, n);
        free(results);
        return result;
    }

    if (arg0->type != EXPR_STRING) return NULL;

    const char* str = arg0->data.string;
    int64_t len = (int64_t)strlen(str);

    /* keep[i] tracks whether position i (0-based) survives the drop. */
    bool* keep = malloc(sizeof(bool) * (len > 0 ? (size_t)len : 1));
    if (!keep) return NULL;
    for (int64_t i = 0; i < len; i++) keep[i] = true;

    /* StringDrop["string", n] - drop first n (n>0) or last |n| (n<0) chars */
    if (spec->type == EXPR_INTEGER) {
        int64_t n = spec->data.integer;
        if (n == 0) {
            /* drop nothing */
        } else if (n > 0) {
            if (n > len) { free(keep); return NULL; }
            for (int64_t i = 0; i < n; i++) keep[i] = false;
        } else {
            int64_t take = -n;
            if (take > len) { free(keep); return NULL; }
            for (int64_t i = len - take; i < len; i++) keep[i] = false;
        }
        Expr* result = stringdrop_build_kept(str, len, keep);
        free(keep);
        return result;
    }

    /* StringDrop["string", UpTo[n]] - drop n, or as many as are available */
    int64_t upto_n;
    if (is_upto(spec, &upto_n)) {
        if (upto_n <= 0) { free(keep); return NULL; }
        int64_t take = upto_n < len ? upto_n : len;
        for (int64_t i = 0; i < take; i++) keep[i] = false;
        Expr* result = stringdrop_build_kept(str, len, keep);
        free(keep);
        return result;
    }

    /* StringDrop["string", {spec}] - list spec */
    if (spec->type == EXPR_FUNCTION &&
        spec->data.function.head->type == EXPR_SYMBOL &&
        spec->data.function.head->data.symbol == SYM_List) {
        size_t spec_argc = spec->data.function.arg_count;

        /* StringDrop["string", {n}] - drop the single character n */
        if (spec_argc == 1) {
            Expr* a1 = spec->data.function.args[0];
            if (a1->type != EXPR_INTEGER) { free(keep); return NULL; }
            int64_t k = a1->data.integer;
            if (k < 0) k = len + k + 1;
            if (k < 1 || k > len) { free(keep); return NULL; }
            keep[k - 1] = false;
            Expr* result = stringdrop_build_kept(str, len, keep);
            free(keep);
            return result;
        }

        /* StringDrop["string", {m, n}] - drop the range m through n */
        if (spec_argc == 2) {
            Expr* a1 = spec->data.function.args[0];
            Expr* a2 = spec->data.function.args[1];
            if (a1->type != EXPR_INTEGER || a2->type != EXPR_INTEGER) {
                free(keep); return NULL;
            }
            int64_t m = a1->data.integer;
            int64_t n = a2->data.integer;
            if (m < 0) m = len + m + 1;
            if (n < 0) n = len + n + 1;
            /* A decreasing range (m > n) is empty and drops nothing. */
            if (m <= n) {
                if (m < 1 || m > len || n < 1 || n > len) {
                    free(keep); return NULL;
                }
                for (int64_t i = m; i <= n; i++) keep[i - 1] = false;
            }
            Expr* result = stringdrop_build_kept(str, len, keep);
            free(keep);
            return result;
        }

        /* StringDrop["string", {m, n, s}] - drop the range m..n in steps of s */
        if (spec_argc == 3) {
            Expr* a1 = spec->data.function.args[0];
            Expr* a2 = spec->data.function.args[1];
            Expr* a3 = spec->data.function.args[2];
            if (a1->type != EXPR_INTEGER || a2->type != EXPR_INTEGER ||
                a3->type != EXPR_INTEGER) {
                free(keep); return NULL;
            }
            int64_t m = a1->data.integer;
            int64_t n = a2->data.integer;
            int64_t s = a3->data.integer;
            if (s == 0) { free(keep); return NULL; }
            if (m < 0) m = len + m + 1;
            if (n < 0) n = len + n + 1;

            for (int64_t p = m; s > 0 ? p <= n : p >= n; p += s) {
                if (p < 1 || p > len) { free(keep); return NULL; }
                keep[p - 1] = false;
            }
            Expr* result = stringdrop_build_kept(str, len, keep);
            free(keep);
            return result;
        }
    }

    free(keep);
    return NULL;
}

/*
 * sr_emit_argx:
 * Emit `StringReverse::argx: StringReverse called with N arguments; 1
 * argument is expected.` to stderr and return NULL so the call is left
 * unevaluated. Matches Mathematica's diagnostic for a wrong argument count.
 */
static Expr* sr_emit_argx(size_t argc) {
    fprintf(stderr,
            "StringReverse::argx: StringReverse called with %zu argument%s; "
            "1 argument is expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

/*
 * builtin_stringreverse:
 * StringReverse["string"] reverses the order of the characters in "string",
 * returning a new string. StringReverse is Listable, so the evaluator threads
 * it element-wise over a list argument before this builtin is reached; each
 * threaded call therefore sees a single string.
 *
 * A call with a number of arguments other than one emits StringReverse::argx
 * and is left unevaluated. A single non-string argument also leaves the call
 * unevaluated (NULL) so symbolic arguments flow through unchanged.
 */
Expr* builtin_stringreverse(Expr* res) {
    if (res->type != EXPR_FUNCTION)
        return NULL;

    size_t argc = res->data.function.arg_count;
    if (argc != 1)
        return sr_emit_argx(argc);

    Expr* arg = res->data.function.args[0];
    if (arg->type != EXPR_STRING)
        return NULL;

    const char* str = arg->data.string;
    size_t len = strlen(str);

    char* buf = malloc(len + 1);
    if (!buf)
        return NULL;

    for (size_t i = 0; i < len; i++)
        buf[i] = str[len - 1 - i];
    buf[len] = '\0';

    Expr* result = expr_new_string(buf);
    free(buf);
    return result;
}

/*
 * strings_init:
 * Registers string builtins with appropriate attributes and docstrings.
 */
void strings_init(void) {
    symtab_add_builtin("StringLength", builtin_stringlength);
    symtab_get_def("StringLength")->attributes |= ATTR_LISTABLE | ATTR_PROTECTED;
    symtab_set_docstring("StringLength",
        "StringLength[\"string\"]\n"
        "\tGives the number of characters in a string.");

    symtab_add_builtin("Characters", builtin_characters);
    symtab_get_def("Characters")->attributes |= ATTR_LISTABLE | ATTR_PROTECTED;
    symtab_set_docstring("Characters",
        "Characters[\"string\"]\n"
        "\tGives a list of the characters in a string.\n"
        "\tEach character is given as a length-1 string.");

    symtab_add_builtin("StringJoin", builtin_stringjoin);
    symtab_get_def("StringJoin")->attributes |= ATTR_FLAT | ATTR_ONEIDENTITY | ATTR_PROTECTED;
    symtab_set_docstring("StringJoin",
        "StringJoin[\"s1\", \"s2\", ...]\n"
        "\tConcatenates strings together.\n"
        "\tStringJoin[{\"s1\", \"s2\", ...}] flattens all lists.\n"
        "\tThe infix form is \"s1\" <> \"s2\" <> ...");

    symtab_add_builtin("StringPart", builtin_stringpart);
    symtab_get_def("StringPart")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("StringPart",
        "StringPart[\"string\", n]\n"
        "\tGives the nth character in \"string\".\n"
        "StringPart[\"string\", {n1, n2, ...}]\n"
        "\tGives a list of the ni-th characters.\n"
        "StringPart[\"string\", m;;n;;s]\n"
        "\tGives characters m through n in steps of s.\n"
        "StringPart[{s1, s2, ...}, spec]\n"
        "\tGives the list of results for each si.\n\n"
        "\tNegative indices count from the end.");

    symtab_add_builtin("StringTake", builtin_stringtake);
    symtab_get_def("StringTake")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("StringTake",
        "StringTake[\"string\", n]\n"
        "\tGives a string containing the first n characters.\n"
        "StringTake[\"string\", -n]\n"
        "\tGives the last n characters.\n"
        "StringTake[\"string\", {n}]\n"
        "\tGives the nth character.\n"
        "StringTake[\"string\", {m, n}]\n"
        "\tGives characters m through n.\n"
        "StringTake[\"string\", {m, n, s}]\n"
        "\tGives characters m through n in steps of s.\n"
        "StringTake[\"string\", UpTo[n]]\n"
        "\tGives n characters, or as many as are available.\n"
        "StringTake[{s1, s2, ...}, spec]\n"
        "\tGives the list of results for each si.");

    symtab_add_builtin("StringDrop", builtin_stringdrop);
    symtab_get_def("StringDrop")->attributes |= ATTR_PROTECTED;
    /* Docstring lives in info.c (info_init). */

    symtab_add_builtin("StringReverse", builtin_stringreverse);
    symtab_get_def("StringReverse")->attributes |= ATTR_LISTABLE | ATTR_PROTECTED;
    /* Docstring lives in info.c (info_init). */
}
