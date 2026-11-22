/*
 * bitwise.c
 *
 * Registration hub for the bit-level integer builtins (src/bitwise/). Each
 * builtin's implementation lives in its own file; this hub only registers the
 * head and assigns its attributes. Docstrings live in src/info.c (info_init),
 * per the project convention.
 */

#include "bitwise.h"
#include "symtab.h"
#include "attr.h"

void bitwise_init(void) {
    symtab_add_builtin("BitLength", builtin_bitlength);
    symtab_get_def("BitLength")->attributes |=
        (ATTR_PROTECTED | ATTR_LISTABLE);

    symtab_add_builtin("BitXor", builtin_bitxor);
    symtab_get_def("BitXor")->attributes |=
        (ATTR_FLAT | ATTR_LISTABLE | ATTR_ONEIDENTITY | ATTR_ORDERLESS |
         ATTR_PROTECTED);

    symtab_add_builtin("BitAnd", builtin_bitand);
    symtab_get_def("BitAnd")->attributes |=
        (ATTR_FLAT | ATTR_LISTABLE | ATTR_ONEIDENTITY | ATTR_ORDERLESS |
         ATTR_PROTECTED);

    symtab_add_builtin("BitOr", builtin_bitor);
    symtab_get_def("BitOr")->attributes |=
        (ATTR_FLAT | ATTR_LISTABLE | ATTR_ONEIDENTITY | ATTR_ORDERLESS |
         ATTR_PROTECTED);

    symtab_add_builtin("BitNot", builtin_bitnot);
    symtab_get_def("BitNot")->attributes |= (ATTR_LISTABLE | ATTR_PROTECTED);

    symtab_add_builtin("BitSet", builtin_bitset);
    symtab_get_def("BitSet")->attributes |= (ATTR_LISTABLE | ATTR_PROTECTED);

    symtab_add_builtin("BitClear", builtin_bitclear);
    symtab_get_def("BitClear")->attributes |= (ATTR_LISTABLE | ATTR_PROTECTED);

    symtab_add_builtin("BitGet", builtin_bitget);
    symtab_get_def("BitGet")->attributes |= (ATTR_LISTABLE | ATTR_PROTECTED);

    symtab_add_builtin("BitShiftLeft", builtin_bitshiftleft);
    symtab_get_def("BitShiftLeft")->attributes |= (ATTR_LISTABLE | ATTR_PROTECTED);

    symtab_add_builtin("BitShiftRight", builtin_bitshiftright);
    symtab_get_def("BitShiftRight")->attributes |= (ATTR_LISTABLE | ATTR_PROTECTED);
}
