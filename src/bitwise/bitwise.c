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
}
