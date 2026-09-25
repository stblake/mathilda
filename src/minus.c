/*
 * minus.c -- Minus[x], the functional form of unary negation.
 *
 * The parser already reads `-x` as Times[-1, x], so Minus only has to exist as
 * a callable head: Minus[x] rewrites to Times[-1, x] and lets Times do the
 * arithmetic (numbers, Rationals, Complex, Plus distribution, packed arrays).
 * This is what makes SortBy[list, Minus] and KeySortBy[a, Minus] work.
 *
 * Attributes match Mathematica: Listable, NumericFunction, Protected.  Any
 * other argument count emits Minus::argx and stays unevaluated (WL prints
 * Minus[x, y] as x − y but does not evaluate it).
 */
#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "common.h"
#include "sym_names.h"

Expr* builtin_minus(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) return builtin_arg_error("Minus", argc, 1, 1);
    Expr* targs[2] = { expr_new_integer(-1), expr_copy(res->data.function.args[0]) };
    return expr_new_function(expr_new_symbol(SYM_Times), targs, 2);
}

void minus_init(void) {
    symtab_add_builtin("Minus", builtin_minus);
    symtab_get_def("Minus")->attributes |= ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED;
    symtab_set_docstring("Minus",
        "Minus[x] is the arithmetic negation of x, equivalent to -x (Times[-1, x]).\n"
        "\tListable; Minus[x, y] (any count other than one) is left unevaluated with Minus::argx.");
}
