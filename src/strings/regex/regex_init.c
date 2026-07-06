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
}
