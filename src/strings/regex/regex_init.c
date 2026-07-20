/*
 * regex_init.c - registration of the regular-expression string builtins.
 *
 * Registers RegularExpression (an inert pattern head) and the four regex-aware
 * string functions.  None are Listable: each hand-threads over a list of
 * subject strings so the pattern/rule argument is never threaded (matching the
 * StringReplacePart / Replace precedent).  Docstrings live in info.c.
 *
 * Called from core_init() (core.c) immediately after strings_init().
 */

#include "picostrings.h"
#include "symtab.h"
#include "attr.h"

/* Register a symbol as an inert protected head (no builtin, no rules). Used for
 * the string-pattern class heads that the WL-pattern translator recognises. */
static void register_inert(const char* name, const char* doc) {
    symtab_get_def(name)->attributes |= ATTR_PROTECTED;
    if (doc) symtab_set_docstring(name, doc);
}

void regex_init(void) {
    symtab_add_builtin("RegularExpression", builtin_regularexpression);
    symtab_get_def("RegularExpression")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("StringMatchQ", builtin_stringmatchq);
    symtab_get_def("StringMatchQ")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("StringCases", builtin_stringcases);
    symtab_get_def("StringCases")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("StringReplace", builtin_stringreplace);
    symtab_get_def("StringReplace")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("StringSplit", builtin_stringsplit);
    symtab_get_def("StringSplit")->attributes |= ATTR_PROTECTED;
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

    symtab_add_builtin("StringExtract", builtin_stringextract);
    symtab_get_def("StringExtract")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("StringPosition", builtin_stringposition);
    symtab_get_def("StringPosition")->attributes |= ATTR_PROTECTED;

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
    register_inert("IgnoreCase", NULL);
    register_inert("Overlaps", NULL);

    symtab_get_def("StringExpression")->attributes |=
        ATTR_FLAT | ATTR_ONEIDENTITY | ATTR_PROTECTED;
    symtab_set_docstring("StringExpression",
        "StringExpression[p1, p2, ...] or p1 ~~ p2 ~~ ...\n"
        "\tRepresents a sequence of string patterns to be matched consecutively.");
}
