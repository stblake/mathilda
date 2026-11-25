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

    symtab_add_builtin("StringPartition", builtin_stringpartition);
    symtab_get_def("StringPartition")->attributes |= ATTR_PROTECTED;
    /* Docstring lives in info.c (info_init). */

    symtab_add_builtin("StringRepeat", builtin_stringrepeat);
    symtab_get_def("StringRepeat")->attributes |= ATTR_PROTECTED;
    /* Docstring lives in info.c (info_init). */

    /* StringPadLeft/StringPadRight thread over a list of strings internally
     * (the 1-arg form pads to the longest element), so they must NOT be
     * ATTR_LISTABLE - that would thread before the handler sees the list. */
    symtab_add_builtin("StringPadLeft", builtin_stringpadleft);
    symtab_get_def("StringPadLeft")->attributes |= ATTR_PROTECTED;
    /* Docstring lives in info.c (info_init). */

    symtab_add_builtin("StringPadRight", builtin_stringpadright);
    symtab_get_def("StringPadRight")->attributes |= ATTR_PROTECTED;
    /* Docstring lives in info.c (info_init). */

    /* StringRiffle takes a (possibly nested) list as its first argument and must
     * inspect the whole structure, so it must NOT be ATTR_LISTABLE. */
    symtab_add_builtin("StringRiffle", builtin_stringriffle);
    symtab_get_def("StringRiffle")->attributes |= ATTR_PROTECTED;
    /* Docstring lives in info.c (info_init). */

    symtab_add_builtin("StringTakeDrop", builtin_stringtakedrop);
    symtab_get_def("StringTakeDrop")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("StringTakeDrop",
        "StringTakeDrop[\"string\", spec]\n"
        "\tGives {StringTake[\"string\", spec], StringDrop[\"string\", spec]}.\n"
        "StringTakeDrop[{s1, s2, ...}, spec]\n"
        "\tGives one pair for each of the si.");

    /* Rotation (stringrotate.c). */
    symtab_add_builtin("StringRotateLeft", builtin_stringrotateleft);
    symtab_get_def("StringRotateLeft")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("StringRotateLeft",
        "StringRotateLeft[\"string\", n]\n"
        "\tCycles the characters of \"string\" n places to the left.\n"
        "StringRotateLeft[\"string\"]\n"
        "\tCycles one place to the left.");
    symtab_add_builtin("StringRotateRight", builtin_stringrotateright);
    symtab_get_def("StringRotateRight")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("StringRotateRight",
        "StringRotateRight[\"string\", n]\n"
        "\tCycles the characters of \"string\" n places to the right.\n"
        "StringRotateRight[\"string\"]\n"
        "\tCycles one place to the right.");

    /* Letter case and character-class predicates (stringcase.c). */
    symtab_add_builtin("ToUpperCase", builtin_touppercase);
    symtab_get_def("ToUpperCase")->attributes |= ATTR_LISTABLE | ATTR_PROTECTED;
    symtab_set_docstring("ToUpperCase",
        "ToUpperCase[\"string\"]\n"
        "\tGives \"string\" with every letter converted to uppercase (full\n"
        "\tUnicode case mapping: ToUpperCase[\"\\[SZ]\"] is \"SS\").");
    symtab_add_builtin("ToLowerCase", builtin_tolowercase);
    symtab_get_def("ToLowerCase")->attributes |= ATTR_LISTABLE | ATTR_PROTECTED;
    symtab_set_docstring("ToLowerCase",
        "ToLowerCase[\"string\"]\n"
        "\tGives \"string\" with every letter converted to lowercase.");
    symtab_add_builtin("Capitalize", builtin_capitalize);
    symtab_get_def("Capitalize")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Capitalize",
        "Capitalize[\"string\"]\n"
        "\tGives \"string\" with its first character in uppercase.\n"
        "Capitalize[{s1, s2, ...}]\n"
        "\tCapitalizes each string.\n"
        "Capitalize[s, \"AllWords\"]\n"
        "\tCapitalizes the first letter of every word.\n"
        "Capitalize[s, \"TitleCase\"]\n"
        "\tCapitalizes every word except short function words not at either end.");
    symtab_add_builtin("Decapitalize", builtin_decapitalize);
    symtab_get_def("Decapitalize")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Decapitalize",
        "Decapitalize[\"string\"]\n"
        "\tGives \"string\" with its first character in lowercase.\n"
        "Decapitalize[{s1, s2, ...}]\n"
        "\tDecapitalizes each string.\n"
        "Decapitalize[s, \"AllWords\"]\n"
        "\tDecapitalizes the first letter of every word.");
    symtab_add_builtin("LetterQ", builtin_letterq);
    symtab_get_def("LetterQ")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("LetterQ",
        "LetterQ[\"string\"]\n"
        "\tGives True if every character of \"string\" is a letter (any\n"
        "\tscript), and False otherwise. LetterQ[\"\"] is True.");
    symtab_add_builtin("DigitQ", builtin_digitq);
    symtab_get_def("DigitQ")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DigitQ",
        "DigitQ[\"string\"]\n"
        "\tGives True if every character of \"string\" is a digit 0-9, and\n"
        "\tFalse otherwise. DigitQ[\"\"] is True.");
    symtab_add_builtin("UpperCaseQ", builtin_uppercaseq);
    symtab_get_def("UpperCaseQ")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("UpperCaseQ",
        "UpperCaseQ[\"string\"]\n"
        "\tGives True if every character of \"string\" is an uppercase letter.");
    symtab_add_builtin("LowerCaseQ", builtin_lowercaseq);
    symtab_get_def("LowerCaseQ")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("LowerCaseQ",
        "LowerCaseQ[\"string\"]\n"
        "\tGives True if every character of \"string\" is a lowercase letter.");

    /* Character codes and ranges (charcode.c). */
    symtab_add_builtin("ToCharacterCode", builtin_tocharactercode);
    symtab_get_def("ToCharacterCode")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("ToCharacterCode",
        "ToCharacterCode[\"string\"]\n"
        "\tGives the list of the Unicode code points of the characters of\n"
        "\t\"string\" (a packed integer array).\n"
        "ToCharacterCode[{s1, s2, ...}]\n"
        "\tGives one list of codes per string.\n"
        "ToCharacterCode[s, \"enc\"]\n"
        "\tGives the codes in encoding \"enc\": \"UTF8\" (the bytes),\n"
        "\t\"ASCII\" or \"ISOLatin1\" (None for a character outside the\n"
        "\tencoding), or \"Unicode\".");
    symtab_add_builtin("FromCharacterCode", builtin_fromcharactercode);
    symtab_get_def("FromCharacterCode")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("FromCharacterCode",
        "FromCharacterCode[n]\n"
        "\tGives the one-character string with Unicode code point n.\n"
        "FromCharacterCode[{n1, n2, ...}]\n"
        "\tGives the string of the characters with those codes.\n"
        "FromCharacterCode[{{n11, ...}, {n21, ...}, ...}]\n"
        "\tGives one string per sublist.\n"
        "FromCharacterCode[codes, \"enc\"]\n"
        "\tDecodes codes in encoding \"enc\": \"UTF8\" (a byte sequence),\n"
        "\t\"ISOLatin1\", \"ASCII\" or \"Unicode\".");
    symtab_add_builtin("CharacterRange", builtin_characterrange);
    symtab_get_def("CharacterRange")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("CharacterRange",
        "CharacterRange[\"c1\", \"c2\"]\n"
        "\tGives the list of characters with codes from that of c1 to that of c2.\n"
        "CharacterRange[n1, n2]\n"
        "\tGives the characters with codes n1 through n2.");
    symtab_add_builtin("Alphabet", builtin_alphabet);
    symtab_get_def("Alphabet")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Alphabet",
        "Alphabet[]\n"
        "\tGives the list of the lowercase letters a through z.\n"
        "Alphabet[\"lang\"]\n"
        "\tGives the lowercase letters of the alphabet of \"English\", \"Greek\" or \"Russian\".");
}
