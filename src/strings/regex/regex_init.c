/*
 * regex_init.c - registration of the regular-expression string builtins.
 *
 * Registers RegularExpression (an inert pattern head) and the regex-aware
 * string functions.  None are Listable: each hand-threads over a list of
 * subject strings so the pattern/rule argument is never threaded (matching the
 * StringReplacePart / Replace precedent).  Docstrings live in info.c.
 *
 * Called from core_init() (core.c) immediately after strings_init().
 */

#include "picostrings.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include "expr.h"
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

/* Register a symbol as an inert protected head (no builtin, no rules). Used for
 * the string-pattern class heads that the WL-pattern translator recognises. */
static void register_inert(const char* name, const char* doc) {
    symtab_get_def(name)->attributes |= ATTR_PROTECTED;
    if (doc) symtab_set_docstring(name, doc);
}

/* Options[name] = {IgnoreCase -> False, Overlaps -> overlaps_default}. Each
 * builtin seeds its option state from these defaults (so SetOptions works) and
 * lets explicit trailing options in the call override them. */
static void register_scan_options(const char* name, const char* overlaps_default) {
    Expr* ic = expr_new_function(expr_new_symbol(SYM_Rule),
        (Expr*[]){ expr_new_symbol(SYM_IgnoreCase), expr_new_symbol(SYM_False) }, 2);
    Expr* ov = expr_new_function(expr_new_symbol(SYM_Rule),
        (Expr*[]){ expr_new_symbol(SYM_Overlaps), expr_new_symbol(overlaps_default) }, 2);
    Expr* opts = expr_new_function(expr_new_symbol(SYM_List),
        (Expr*[]){ ic, ov }, 2);
    symtab_set_options(name, opts);
}

/* Options[name] = {IgnoreCase -> False}.  The substring predicates take no
 * Overlaps option (a Boolean answer does not depend on the overlap policy), so
 * they get this one-option sibling of register_scan_options. */
/* StringExpression[p1, p2, ...]: adjacent literal strings are one literal
 * ("a" ~~ "b" is "ab", "a" ~~ "b" ~~ x ~~ "c" is "ab" ~~ x ~~ "c"), and a
 * lone string is itself. The match is unchanged: consecutive literals match
 * consecutively either way. Anything else stays as it is. */
static Expr* builtin_stringexpression(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t n = res->data.function.arg_count;
    Expr** a = res->data.function.args;
    bool adjacent = false;
    for (size_t i = 0; i + 1 < n; i++)
        if (a[i]->type == EXPR_STRING && a[i + 1]->type == EXPR_STRING) { adjacent = true; break; }
    if (!adjacent) {
        if (n == 1 && a[0]->type == EXPR_STRING) return expr_copy(a[0]);
        return NULL;
    }
    Expr** out = malloc(sizeof(Expr*) * n);
    if (!out) return NULL;
    size_t m = 0;
    for (size_t i = 0; i < n; ) {
        if (a[i]->type != EXPR_STRING) { out[m++] = expr_copy(a[i]); i++; continue; }
        size_t j = i, total = 0;
        while (j < n && a[j]->type == EXPR_STRING) total += strlen(a[j++]->data.string);
        char* buf = malloc(total + 1);
        if (!buf) { for (size_t k = 0; k < m; k++) expr_free(out[k]); free(out); return NULL; }
        size_t o = 0;
        for (size_t k = i; k < j; k++) {
            size_t l = strlen(a[k]->data.string);
            memcpy(buf + o, a[k]->data.string, l); o += l;
        }
        buf[o] = '\0';
        out[m++] = expr_new_string_adopt(buf);
        i = j;
    }
    Expr* r = m == 1 ? out[0]
                     : expr_new_function(expr_copy(res->data.function.head), out, m);
    free(out);
    return r;
}

static void register_ignorecase_option(const char* name) {
    Expr* ic = expr_new_function(expr_new_symbol(SYM_Rule),
        (Expr*[]){ expr_new_symbol(SYM_IgnoreCase), expr_new_symbol(SYM_False) }, 2);
    Expr* opts = expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ ic }, 1);
    symtab_set_options(name, opts);
}

/* Options[StringReplace]: IgnoreCase -> False, MetaCharacters -> None. */
static void register_replace_options(const char* name) {
    Expr* ic = expr_new_function(expr_new_symbol(SYM_Rule),
        (Expr*[]){ expr_new_symbol(SYM_IgnoreCase), expr_new_symbol(SYM_False) }, 2);
    Expr* mc = expr_new_function(expr_new_symbol(SYM_Rule),
        (Expr*[]){ expr_new_symbol("MetaCharacters"), expr_new_symbol(SYM_None) }, 2);
    symtab_set_options(name, expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){ ic, mc }, 2));
}

void regex_init(void) {
    symtab_add_builtin("RegularExpression", builtin_regularexpression);
    symtab_get_def("RegularExpression")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("StringMatchQ", builtin_stringmatchq);
    symtab_get_def("StringMatchQ")->attributes |= ATTR_PROTECTED;
    /* Options[StringMatchQ]: IgnoreCase -> False, SpellingCorrection -> False
     * (accepted; no spelling correction is ever applied). */
    symtab_set_options("StringMatchQ", expr_new_function(expr_new_symbol(SYM_List), (Expr*[]){
        expr_new_function(expr_new_symbol(SYM_Rule),
            (Expr*[]){ expr_new_symbol(SYM_IgnoreCase), expr_new_symbol(SYM_False) }, 2),
        expr_new_function(expr_new_symbol(SYM_Rule),
            (Expr*[]){ expr_new_symbol(SYM_SpellingCorrection), expr_new_symbol(SYM_False) }, 2) }, 2));

    symtab_add_builtin("StringCases", builtin_stringcases);
    symtab_get_def("StringCases")->attributes |= ATTR_PROTECTED;
    register_scan_options("StringCases", SYM_False);

    symtab_add_builtin("StringCount", builtin_stringcount);
    symtab_get_def("StringCount")->attributes |= ATTR_PROTECTED;
    register_scan_options("StringCount", SYM_False);

    symtab_add_builtin("StringReplace", builtin_stringreplace);
    symtab_get_def("StringReplace")->attributes |= ATTR_PROTECTED;
    register_replace_options("StringReplace");

    symtab_add_builtin("StringDelete", builtin_stringdelete);
    symtab_get_def("StringDelete")->attributes |= ATTR_PROTECTED;
    register_ignorecase_option("StringDelete");
    symtab_set_docstring("StringDelete",
        "StringDelete[\"string\", patt]\n"
        "\tDeletes every substring of \"string\" that matches patt.\n"
        "StringDelete[{s1, s2, ...}, patt]\n"
        "\tGives the list of results for each of the si.\n"
        "StringDelete[patt]\n"
        "\tRepresents an operator form that can be applied to a string.\n\n"
        "\tEquivalent to StringReplace[\"string\", patt -> \"\"].\n"
        "\tOptions: IgnoreCase -> True.");

    symtab_add_builtin("StringReplaceList", builtin_stringreplacelist);
    symtab_get_def("StringReplaceList")->attributes |= ATTR_PROTECTED;
    register_ignorecase_option("StringReplaceList");
    symtab_set_docstring("StringReplaceList",
        "StringReplaceList[\"string\", rules]\n"
        "\tGives the list of strings obtained by replacing ONE match of the\n"
        "\trules at a time, for every match start of every rule, in order of\n"
        "\tposition.\n"
        "StringReplaceList[\"string\", rules, n]\n"
        "\tGives at most n results.\n"
        "StringReplaceList[{s1, s2, ...}, rules]\n"
        "\tGives the list of results for each of the si.\n\n"
        "\tOptions: IgnoreCase -> True.");

    symtab_add_builtin("StringSplit", builtin_stringsplit);
    symtab_get_def("StringSplit")->attributes |= ATTR_PROTECTED;
    register_ignorecase_option("StringSplit");
    symtab_set_docstring("StringSplit",
        "StringSplit[\"string\"]\n"
        "\tSplits \"string\" at runs of whitespace.\n"
        "StringSplit[\"string\", patt]\n"
        "\tSplits at delimiters matching the string pattern patt.\n"
        "StringSplit[\"string\", {p1, p2, ...}]\n"
        "\tSplits at any of the pi.\n"
        "StringSplit[\"string\", patt -> val]\n"
        "\tInserts val at the position of each delimiter.\n"
        "StringSplit[\"string\", patt, n]\n"
        "\tSplits into at most n substrings.\n"
        "StringSplit[{s1, s2, ...}, patt]\n"
        "\tGives the list of results for each of the si.\n\n"
        "\tEmpty substrings between adjacent interior delimiters are kept; those\n"
        "\tat the start or end are dropped unless All is given as the third\n"
        "\targument. \"\" splits at every character. Option: IgnoreCase -> True.");

    symtab_add_builtin("StringTrim", builtin_stringtrim);
    symtab_get_def("StringTrim")->attributes |= ATTR_PROTECTED;
    register_ignorecase_option("StringTrim");

    symtab_add_builtin("StringExtract", builtin_stringextract);
    symtab_get_def("StringExtract")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("StringPosition", builtin_stringposition);
    symtab_get_def("StringPosition")->attributes |= ATTR_PROTECTED;
    /* StringPosition defaults to Overlaps -> True, unlike StringCases and
     * StringCount, which default to False. This matches Wolfram Language. */
    register_scan_options("StringPosition", SYM_True);

    /* The substring predicates.  StringFreeQ is the negation of
     * StringContainsQ; StringStartsQ / StringEndsQ are the one-end-anchored
     * variants.  All four share one core in stringcontainsq.c. */
    symtab_add_builtin("StringContainsQ", builtin_stringcontainsq);
    symtab_get_def("StringContainsQ")->attributes |= ATTR_PROTECTED;
    register_ignorecase_option("StringContainsQ");

    symtab_add_builtin("StringFreeQ", builtin_stringfreeq);
    symtab_get_def("StringFreeQ")->attributes |= ATTR_PROTECTED;
    register_ignorecase_option("StringFreeQ");

    symtab_add_builtin("StringStartsQ", builtin_stringstartsq);
    symtab_get_def("StringStartsQ")->attributes |= ATTR_PROTECTED;
    register_ignorecase_option("StringStartsQ");

    symtab_add_builtin("StringEndsQ", builtin_stringendsq);
    symtab_get_def("StringEndsQ")->attributes |= ATTR_PROTECTED;
    register_ignorecase_option("StringEndsQ");

    /* Inert protected heads understood by the string-pattern translator
     * (string_pattern.c). StringExpression (the ~~ operator's head) is Flat so
     * a ~~ b ~~ c collapses to StringExpression[a, b, c]. */
    register_inert("Whitespace",
        "Whitespace\n\tA string pattern matching a run of whitespace characters.");
    register_inert("WhitespaceCharacter",
        "WhitespaceCharacter\n\tA string pattern matching a single whitespace character.");
    register_inert("LetterCharacter",
        "LetterCharacter\n\tA string pattern matching a single letter.");
    register_inert("DigitCharacter",
        "DigitCharacter\n\tA string pattern matching a single digit.");
    register_inert("WordCharacter",
        "WordCharacter\n\tA string pattern matching a single letter or digit.");
    register_inert("NumberString",
        "NumberString\n\tA string pattern matching a signed integer or decimal number.");
    register_inert("StartOfString",
        "StartOfString\n\tA zero-width string pattern matching the start of the whole string.");
    register_inert("EndOfString",
        "EndOfString\n\tA zero-width string pattern matching the end of the whole string.");
    register_inert("HexadecimalCharacter",
        "HexadecimalCharacter\n\tA string pattern matching a hexadecimal digit 0-9, a-f, A-F.");
    register_inert("PunctuationCharacter",
        "PunctuationCharacter\n\tA string pattern matching a punctuation character.");
    register_inert("StartOfLine",
        "StartOfLine\n\tA zero-width string pattern matching the start of a line.");
    register_inert("EndOfLine",
        "EndOfLine\n\tA zero-width string pattern matching the end of a line.");
    register_inert("WordBoundary",
        "WordBoundary\n\tA zero-width string pattern matching a boundary between word and non-word characters.");
    register_inert("IgnoreCase", NULL);
    register_inert("Overlaps", NULL);

    symtab_add_builtin("StringExpression", builtin_stringexpression);
    symtab_get_def("StringExpression")->attributes |=
        ATTR_FLAT | ATTR_ONEIDENTITY | ATTR_PROTECTED;
    symtab_set_docstring("StringExpression",
        "StringExpression[p1, p2, ...] or p1 ~~ p2 ~~ ...\n"
        "\tRepresents a sequence of string patterns to be matched consecutively.");
}
