/*
 * stringcontainsq.c - the substring-predicate family:
 *
 *   StringContainsQ["string", patt]   True if SOME substring matches patt
 *   StringFreeQ    ["string", patt]   True if NO   substring matches patt
 *   StringStartsQ  ["string", patt]   True if a PREFIX of "string" matches patt
 *   StringEndsQ    ["string", patt]   True if a SUFFIX of "string" matches patt
 *
 * All four ask one question - "does any rule match anywhere?" - so they share a
 * single core (sq_dispatch) and differ only in an SqKind tag.  The Wolfram
 * Language identities they implement are exactly the ones the four faces are
 * derived from:
 *
 *   StringFreeQ  [s, p] == !StringContainsQ[s, p]
 *   StringStartsQ[s, p] ==  StringContainsQ[s, StartOfString ~~ p]
 *   StringEndsQ  [s, p] ==  StringContainsQ[s, p ~~ EndOfString]
 *
 * Anchoring.  Rather than thread a new anchor mode through the shared builder
 * regex_rules_build_ex() (which eight call sites depend on), StringStartsQ /
 * StringEndsQ anchor at the *pattern* level: they synthesise a temporary
 * StringExpression[StartOfString, patt] (resp. StringExpression[patt,
 * EndOfString]) and hand that to the unmodified builder unanchored.  This is
 * safe because the translator's group_join() wraps every child in `(?:...)`, so
 * an alternation-bearing child such as RegularExpression["a|b"] becomes
 * `(?:(?:\A)(?:a|b))` rather than the broken `\Aa|b`.  A top-level List pattern
 * lands in the translator's nested-list branch and becomes alternatives, which
 * for a boolean any-match test is equivalent to the usual one-rule-per-element
 * split.  The wrapper Expr must outlive the rule set: RegexRule.lhs borrows
 * into it.
 *
 * Matching uses one regex_match() per rule with an early exit on the first hit,
 * not the regex_scan() enumerator the counting/extraction builtins use - a
 * predicate never needs to know where or how often the pattern occurred.
 *
 * Threading.  Like the rest of the regex string family these are NOT Listable;
 * each hand-threads over a list of subject strings so the pattern argument is
 * never threaded.  They do deviate from StringCount / StringMatchQ in one way:
 * a non-string subject (or any non-string element of a subject list) leaves the
 * whole call unevaluated rather than passing the offending element through.  A
 * boolean predicate answering {False, 7} would be wrong; unevaluated is what
 * the Wolfram Language does.
 *
 * Options: IgnoreCase -> True | False (default False, overridable via
 * SetOptions).  Operator form: with a single positional argument the call
 * returns Function[head[#1, patt, opts...]], so
 * StringFreeQ["ac", IgnoreCase -> True]["BACCD"] is False.
 *
 * Offsets and lengths are byte-based, matching the rest of the string family.
 */

#include "picostrings.h"
#include "regex_common.h"
#include "sym_names.h"
#include "symtab.h"
#include "common.h"

#include <stdlib.h>
#include <string.h>

/* Which of the four faces is being evaluated. */
typedef enum { SQ_CONTAINS, SQ_FREE, SQ_STARTS, SQ_ENDS } SqKind;

/* Does any rule match anywhere in `subj`?  Short-circuits on the first hit. */
static int sq_match(const char* subj, RegexRule* rules, int nr) {
    size_t len = strlen(subj);
    size_t ov[2];
    for (int i = 0; i < nr; i++)
        if (regex_match(rules[i].prog, subj, len, 0, ov, 1) == 1)
            return 1;
    return 0;
}

/* Turn a match/no-match into the answer for this face (StringFreeQ inverts). */
static Expr* sq_answer(int matched, SqKind kind) {
    int yes = (kind == SQ_FREE) ? !matched : matched;
    return expr_new_symbol(yes ? SYM_True : SYM_False);
}

/*
 * Apply the rule set to the subject: a single string gives a single Boolean, a
 * list of strings gives the list of Booleans.  Anything else (a non-string, or
 * a list holding one) returns NULL so the call is left unevaluated.
 *
 * `nr == 0` is a legitimate input meaning "no patterns", so an empty pattern
 * list falls out as no-match without a special case.
 */
static Expr* sq_build_result(Expr* subject, RegexRule* rules, int nr, SqKind kind) {
    if (subject->type == EXPR_STRING)
        return sq_answer(sq_match(subject->data.string, rules, nr), kind);

    if (!head_is(subject, SYM_List)) return NULL;

    size_t n = subject->data.function.arg_count;
    Expr** elems = subject->data.function.args;
    for (size_t i = 0; i < n; i++)
        if (elems[i]->type != EXPR_STRING) return NULL;

    Expr** out = malloc(sizeof(Expr*) * (n ? n : 1));
    if (!out) return NULL;
    for (size_t i = 0; i < n; i++)
        out[i] = sq_answer(sq_match(elems[i]->data.string, rules, nr), kind);

    Expr* result = expr_new_function(expr_new_symbol(SYM_List), out, n);
    free(out);                              /* expr_new_function copies the array */
    return result;
}

/*
 * The one-end-anchored pattern for StringStartsQ / StringEndsQ, as a fresh
 * StringExpression the caller owns.  NULL for the unanchored faces.
 */
static Expr* sq_anchor_pattern(Expr* patt, SqKind kind) {
    Expr* args[2];
    if (kind == SQ_STARTS) {
        args[0] = expr_new_symbol(SYM_StartOfString);
        args[1] = expr_copy(patt);
    } else if (kind == SQ_ENDS) {
        args[0] = expr_copy(patt);
        args[1] = expr_new_symbol(SYM_EndOfString);
    } else {
        return NULL;
    }
    return expr_new_function(expr_new_symbol(SYM_StringExpression), args, 2);
}

/* The interned head symbol for a face, needed to rebuild the operator form. */
static const char* sq_head_symbol(SqKind kind) {
    switch (kind) {
        case SQ_FREE:   return SYM_StringFreeQ;
        case SQ_STARTS: return SYM_StringStartsQ;
        case SQ_ENDS:   return SYM_StringEndsQ;
        default:        return SYM_StringContainsQ;
    }
}

/*
 * Operator (curried) form: head[patt, opts...] becomes the pure function
 * Function[head[#1, patt, opts...]], mirroring the Cases / StringReplacePart
 * precedent.  The trailing options are carried through so a curried
 * IgnoreCase -> True still applies when the operator is finally applied.
 */
static Expr* sq_operator_form(Expr** a, size_t argc, SqKind kind) {
    Expr** inner = malloc(sizeof(Expr*) * (argc + 1));
    if (!inner) return NULL;

    Expr* slot_args[1] = { expr_new_integer(1) };
    inner[0] = expr_new_function(expr_new_symbol(SYM_Slot), slot_args, 1);
    for (size_t i = 0; i < argc; i++) inner[i + 1] = expr_copy(a[i]);

    Expr* body = expr_new_function(expr_new_symbol(sq_head_symbol(kind)),
                                   inner, argc + 1);
    free(inner);

    Expr* func_args[1] = { body };
    return expr_new_function(expr_new_symbol(SYM_Function), func_args, 1);
}

/* Shared evaluation core; `head` is the printable name used in diagnostics. */
static Expr* sq_dispatch(Expr* res, SqKind kind, const char* head) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** a = res->data.function.args;

    if (argc == 0) return builtin_arg_error(head, 0, 1, 2);

    /* Seed IgnoreCase from the registered defaults so SetOptions takes effect;
     * an explicit trailing option below overrides it. */
    int caseless = 0;
    Expr* defs = symtab_get_options(head);          /* borrowed */
    if (defs && defs->type == EXPR_FUNCTION) {
        for (size_t i = 0; i < defs->data.function.arg_count; i++) {
            int v;
            if (regex_match_opt(defs->data.function.args[i], SYM_IgnoreCase, &v, 0))
                caseless = v;
        }
    }

    /* Strip trailing IgnoreCase options, leaving the positional arguments.
     * Unlike StringCount / StringSplit this strips down to a single positional
     * argument, because head[patt, IgnoreCase -> True] is the operator form
     * rather than a subject/pattern pair. */
    size_t pargc = argc;
    while (pargc >= 1) {
        int v;
        if (regex_match_opt(a[pargc - 1], SYM_IgnoreCase, &v, 0)) { caseless = v; pargc--; }
        else break;
    }
    if (pargc == 0 || pargc > 2) return builtin_arg_error(head, argc, 1, 2);

    if (pargc == 1) return sq_operator_form(a, argc, kind);

    Expr* subject = a[0];
    Expr* patt = a[1];

    /* An empty pattern list matches nothing.  regex_rules_build_ex reports that
     * as -1, indistinguishable from a genuine failure, so answer it directly:
     * sq_build_result with an empty rule set is exactly "no match". */
    if (head_is(patt, SYM_List) && patt->data.function.arg_count == 0)
        return sq_build_result(subject, NULL, 0, kind);

    /* StringStartsQ / StringEndsQ anchor by wrapping the pattern; the wrapper
     * has to outlive the rule set because RegexRule.lhs borrows into it. */
    Expr* wrapper = sq_anchor_pattern(patt, kind);
    RegexRule* rules;
    int nr = regex_rules_build_ex(wrapper ? wrapper : patt, /*anchored=*/0,
                                  caseless, &rules, head);
    if (nr < 0) {
        if (wrapper) expr_free(wrapper);
        return NULL;
    }

    Expr* result = sq_build_result(subject, rules, nr, kind);
    regex_rules_free(rules, nr);
    if (wrapper) expr_free(wrapper);
    return result;
}

Expr* builtin_stringcontainsq(Expr* res) {
    return sq_dispatch(res, SQ_CONTAINS, "StringContainsQ");
}

Expr* builtin_stringfreeq(Expr* res) {
    return sq_dispatch(res, SQ_FREE, "StringFreeQ");
}

Expr* builtin_stringstartsq(Expr* res) {
    return sq_dispatch(res, SQ_STARTS, "StringStartsQ");
}

Expr* builtin_stringendsq(Expr* res) {
    return sq_dispatch(res, SQ_ENDS, "StringEndsQ");
}
