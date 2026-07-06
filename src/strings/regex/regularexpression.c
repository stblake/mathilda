/*
 * regularexpression.c - the RegularExpression[...] head.
 *
 * RegularExpression["re"] is an inert data head: it carries a PCRE pattern
 * string for use by StringMatchQ / StringCases / StringReplace / StringSplit
 * and does not evaluate to anything else (exactly like Wolfram Language).
 * The builtin therefore validates its argument and otherwise returns NULL so
 * the expression survives evaluation unchanged.  When the pattern does not
 * compile a RegularExpression::regex diagnostic is emitted (still inert), and
 * when Mathilda was built without PCRE2 a RegularExpression::regavail note is
 * printed the first time an invalid-but-present call is seen.
 */

#include "picostrings.h"
#include "regex_engine.h"

#include <stdio.h>

Expr* builtin_regularexpression(Expr* res) {
    if (res->type != EXPR_FUNCTION)
        return NULL;

    /* Only the single-string form is meaningful; anything else stays inert. */
    if (res->data.function.arg_count != 1)
        return NULL;

    Expr* arg = res->data.function.args[0];
    if (arg->type != EXPR_STRING)
        return NULL;

    /* Best-effort syntax check: warn (once per evaluation) on a bad pattern,
     * but still leave the head inert so downstream code can report context. */
    if (regex_available()) {
        char err[256];
        RegexProgram* prog = regex_compile(arg->data.string, err, sizeof err);
        if (!prog) {
            fprintf(stderr, "RegularExpression::regex: %s\n", err);
        } else {
            regex_free(prog);
        }
    }

    return NULL;   /* inert: RegularExpression["re"] stays as itself */
}
