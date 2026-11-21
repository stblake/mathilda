/*
 * stringtrim.c - StringTrim[...], trims matching substrings from both ends.
 *
 *   StringTrim["string"]           trims whitespace runs from start and end
 *   StringTrim["string", patt]     trims a substring matching patt from each end
 *   StringTrim[{s1, s2, ...}, ...] threads over a list of subject strings
 *   IgnoreCase -> True             case-insensitive pattern
 *
 * As in Mathematica, ONE match is removed at each end -- StringTrim is
 * StringReplace[s, {StartOfString ~~ patt -> "", patt ~~ EndOfString -> ""}]:
 * StringTrim["xxabcxx", "x"] is "xabcx" (use "x".. to strip runs). The back
 * match must start at or after the end of the front match, so StringTrim["xxx",
 * "xx"] is "x".
 *
 * Both ends go through the shared string-pattern engine (regex_common.c), so
 * the full pattern language works, Condition and PatternTest included:
 *   - front: the pattern compiled REGEX_ANCHORED, one match at offset 0 -- no
 *     scan of the subject at all;
 *   - back: the pattern compiled REGEX_ENDANCHORED, searched from the end of the
 *     front match; PCRE rejects a candidate as soon as it cannot reach the end.
 */

#include "picostrings.h"
#include "strings_internal.h"
#include "regex_common.h"
#include "sym_names.h"
#include "symtab.h"
#include "common.h"
#include "message.h"
#include "print.h"

#include <stdlib.h>
#include <string.h>

/* Trim one string subject, returning a fresh EXPR_STRING. */
static Expr* st_scalar(const char* subj, RegexRule* front, RegexRule* back) {
    size_t len = strlen(subj);
    size_t start = 0, end = len;
    size_t ov[2];

    if (regex_match(front->prog, subj, len, 0, ov, 1) == 1)
        start = ov[1];
    if (start < len && regex_match(back->prog, subj, len, start, ov, 1) == 1 &&
        ov[1] == len)
        end = ov[0];

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

    /* Trailing IgnoreCase option (seeded from Options[StringTrim]). */
    int caseless = 0;
    Expr* defs = symtab_get_options("StringTrim");   /* borrowed */
    if (defs && defs->type == EXPR_FUNCTION)
        for (size_t i = 0; i < defs->data.function.arg_count; i++)
            regex_match_opt(defs->data.function.args[i], SYM_IgnoreCase, &caseless, 0);
    size_t pargc = argc;
    while (pargc >= 2 && regex_match_opt(a[pargc - 1], SYM_IgnoreCase, &caseless, 0))
        pargc--;

    if (pargc < 1 || pargc > 2) return builtin_arg_error("StringTrim", argc, 1, 2);

    Expr* subject = a[0];
    int is_list = strings_is_seq(subject);      /* a List or a LinkedList */
    int ok = is_list ? strings_seq_all_strings(subject) : subject->type == EXPR_STRING;
    if (!ok) {
        char* shown = expr_to_string_outputform(res);
        mth_message("StringTrim", "strse",
                    "String or list of strings expected at position 1 in %s.",
                    shown ? shown : "StringTrim");
        free(shown);
        return NULL;
    }

    /* Trim pattern (default Whitespace). */
    Expr* patt = (pargc == 2) ? expr_copy(a[1]) : expr_new_symbol(SYM_Whitespace);
    /* A list of patterns is a set of alternatives (one rule, not several). */
    if (patt->type == EXPR_FUNCTION && patt->data.function.head->type == EXPR_SYMBOL &&
        patt->data.function.head->data.symbol.name == SYM_List) {
        patt = expr_unshare(patt);
        expr_free(patt->data.function.head);
        patt->data.function.head = expr_new_symbol(SYM_Alternatives);
    } else if (patt->type == EXPR_LINKEDLIST) {    /* the same, from the slots */
        Expr* alt = ll_call_with_elements(expr_new_symbol(SYM_Alternatives), patt->data.llist);
        expr_free(patt);
        patt = alt;
    }

    RegexRule* front = NULL;
    RegexRule* back  = NULL;
    unsigned cf = caseless ? REGEX_CASELESS : 0u;
    int nf = regex_rules_build_opts(patt, 0, REGEX_ANCHORED | cf, &front, "StringTrim");
    int nb = (nf == 1)
           ? regex_rules_build_opts(patt, 0, REGEX_ENDANCHORED | cf, &back, "StringTrim")
           : -1;
    Expr* result = NULL;
    if (nf == 1 && nb == 1) {
        if (is_list) {
            LLElems sv;
            ll_elems_open(subject, &sv);
            size_t m = sv.n;
            Expr** out = malloc(sizeof(Expr*) * (m ? m : 1));
            if (out) {
                for (size_t i = 0; i < m; i++)
                    out[i] = st_scalar(sv.v[i]->data.string, front, back);
                result = ll_seq_like(subject, out, m);
                free(out);
            }
            ll_elems_close(&sv);
        } else {
            result = st_scalar(subject->data.string, front, back);
        }
    }
    /* An unsupported pattern leaves the call unevaluated. */
    if (nf >= 0) regex_rules_free(front, nf);
    if (nb >= 0) regex_rules_free(back, nb);
    expr_free(patt);       /* after the rules: they borrow into it */
    return result;
}
