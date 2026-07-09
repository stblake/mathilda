/*
 * strings_init.c - registration of the string manipulation builtins.
 *
 * Each builtin is implemented in its own translation unit under src/strings/;
 * this file wires them into the symbol table with their attributes and
 * docstrings. Some docstrings live in info.c (info_init) instead, noted below.
 */

#include "picostrings.h"
#include "symtab.h"
#include "attr.h"

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

    symtab_add_builtin("StringInsert", builtin_stringinsert);
    symtab_get_def("StringInsert")->attributes |= ATTR_PROTECTED;
    /* Docstring lives in info.c (info_init). */

    symtab_add_builtin("StringReplacePart", builtin_stringreplacepart);
    symtab_get_def("StringReplacePart")->attributes |= ATTR_PROTECTED;
    /* Docstring lives in info.c (info_init). */
}
