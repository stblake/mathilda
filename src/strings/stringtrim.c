/*
 * stringtrim.c - StringTrim[...], trims matching substrings from both ends.
 *
 *   StringTrim["string"]           trims whitespace runs from start and end
 *   StringTrim["string", patt]     trims substrings matching patt from both ends
 *   StringTrim[{s1, s2, ...}, ...] threads over a list of subject strings
 *
 * The trim pattern is translated to PCRE by the shared string-pattern engine
 * (string_pattern.c: wl_pattern_to_regex), so it accepts literal strings,
 * RegularExpression["re"], the character-class heads (Whitespace, ...),
 * StringExpression (~~), Alternatives (|), Repeated (..), Except, and so on.
 * The default pattern is Whitespace (a run of spaces / tabs / newlines).
 *
 * Anchoring: rather than the whole-string \A...\z wrap the other regex builtins
 * use, StringTrim needs the pattern anchored at just one end.
 *   - The front is stripped by matching `(?:src)` starting at the current
 *     position and accepting only a match that begins exactly there.  PCRE's
 *     `\A`/`^` anchor to absolute offset 0, so a start-relative anchor is done
 *     by the ov[0] == start check, not by wrapping.
 *   - The back is stripped by matching `(?:src)\z` while passing a truncated
 *     subject length (the current end) to regex_match, so `\z` refers to the
 *     current end rather than the absolute string end.
 * Each end is stripped repeatedly to a fixed point (so StringTrim["xxabcxx",
 * "x"] -> "abc"); a zero-width match ends the loop.
 */

#include "picostrings.h"
#include "regex_common.h"
#include "sym_names.h"
#include "common.h"

#include <stdlib.h>
#include <string.h>

/*
 * Compile a one-end-anchored trim program from raw PCRE source `src`.
 * `anchor_end == 0` builds the front matcher `(?:src)`; `anchor_end != 0`
 * builds the back matcher `(?:src)\z`.  Returns an owned RegexProgram (free
 * with regex_free), or NULL on a compile error / when regex support is absent.
 */
static RegexProgram* st_compile(const char* src, int anchor_end) {
    size_t n = strlen(src);
    size_t extra = anchor_end ? 6 : 4;    /* "(?:" + ")" (+ "\z") + NUL slack */
    char* pat = malloc(n + extra + 1);
    if (!pat) return NULL;

    size_t j = 0;
    memcpy(pat + j, "(?:", 3); j += 3;
    memcpy(pat + j, src, n);   j += n;
    pat[j++] = ')';
    if (anchor_end) { pat[j++] = '\\'; pat[j++] = 'z'; }
    pat[j] = '\0';

    char err[128];
    RegexProgram* p = regex_compile(pat, err, sizeof err);
    free(pat);
    return p;
}

/* Trim one string subject, returning a fresh EXPR_STRING. */
static Expr* st_scalar(const char* subj, RegexProgram* front, RegexProgram* back) {
    size_t len = strlen(subj);
    size_t start = 0, end = len;
    size_t ov[2];

    /* Strip the leading run: repeated matches anchored at the current start. */
    while (start < end) {
        if (regex_match(front, subj, end, start, ov, 1) != 1) break;
        if (ov[0] != start || ov[1] <= start) break;   /* not at start / zero-width */
        start = ov[1];
    }

    /* Strip the trailing run: `\z`-anchored matches over subj[0..end). */
    while (end > start) {
        if (regex_match(back, subj, end, start, ov, 1) != 1) break;
        if (ov[0] >= end) break;                        /* zero-width at end */
        end = ov[0];
    }

    size_t n = end - start;
    char* buf = malloc(n + 1);
    if (!buf) return expr_new_string("");
    memcpy(buf, subj + start, n);
    buf[n] = '\0';
    Expr* r = expr_new_string(buf);
    free(buf);
    return r;
}

Expr* builtin_stringtrim(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** a = res->data.function.args;

    if (argc < 1 || argc > 2) return builtin_arg_error("StringTrim", argc, 1, 2);

    Expr* subject = a[0];

    /* Trim pattern (default Whitespace). */
    int own_patt = 0;
    Expr* patt;
    if (argc >= 2) {
        patt = a[1];
    } else {
        patt = expr_new_symbol(SYM_Whitespace);
        own_patt = 1;
    }

    char* src = wl_pattern_to_regex(patt, NULL);
    if (own_patt) expr_free(patt);
    if (!src) return NULL;                              /* unsupported pattern */

    RegexProgram* front = st_compile(src, 0);
    RegexProgram* back  = st_compile(src, 1);
    free(src);
    if (!front || !back) {                              /* compile error / no regex */
        regex_free(front);
        regex_free(back);
        return NULL;
    }

    Expr* result;
    if (subject->type == EXPR_FUNCTION &&
        subject->data.function.head->type == EXPR_SYMBOL &&
        subject->data.function.head->data.symbol.name == SYM_List) {
        size_t m = subject->data.function.arg_count;
        Expr** out = malloc(sizeof(Expr*) * (m ? m : 1));
        for (size_t i = 0; i < m; i++) {
            Expr* si = subject->data.function.args[i];
            out[i] = (si->type == EXPR_STRING)
                         ? st_scalar(si->data.string, front, back)
                         : expr_copy(si);
        }
        result = expr_new_function(expr_new_symbol(SYM_List), out, m);
        free(out);
    } else if (subject->type == EXPR_STRING) {
        result = st_scalar(subject->data.string, front, back);
    } else {
        result = NULL;                                  /* non-string subject */
    }

    regex_free(front);
    regex_free(back);
    return result;
}
