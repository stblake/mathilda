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

/* Internal`SyntacticNegativeQ[e] -- does `e` carry a minus sign in the LEADING
 * position of its printed form? A *syntactic* question, not a sign test:
 *
 *   -3, -2.5, -1/2, -x, -2 x y, -1 - I     True
 *   3, x, x - y, -a - b, 1 - 2 x           False
 *
 * Note the Plus rows. `-a - b` prints with a leading minus but is False, because
 * the head is Plus and the predicate does not look inside it -- which is exactly
 * why callers that want "prints with a leading minus" write
 * `SyntacticNegativeQ[e] || (Head[e] === Plus && SyntacticNegativeQ[First[e]])`.
 * Keeping that asymmetry is the point: it is what SymPy's
 * could_extract_minus_sign is built on, and the ported integrator reproduces
 * that function's term ordering from it.
 *
 * Times defers to its first argument, since the numeric coefficient sorts first
 * in canonical order. A negative real part decides a Complex, its imaginary part
 * breaking the tie when the real part is zero. */
static bool syntactic_negative(const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: return e->data.integer < 0;
        case EXPR_REAL:    return e->data.real < 0.0;
        case EXPR_BIGINT:  return mpz_sgn(e->data.bigint) < 0;
        case EXPR_FUNCTION: break;
        default: return false;
    }
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    size_t argc  = e->data.function.arg_count;
    if (h == SYM_Rational && argc == 2)
        return syntactic_negative(e->data.function.args[0]);   /* denominator > 0 */
    if (h == SYM_Complex && argc == 2) {
        const Expr* re = e->data.function.args[0];
        if (syntactic_negative(re)) return true;
        bool re_zero = (re->type == EXPR_INTEGER && re->data.integer == 0) ||
                       (re->type == EXPR_REAL && re->data.real == 0.0);
        return re_zero && syntactic_negative(e->data.function.args[1]);
    }
    if (h == SYM_Times && argc >= 1)
        return syntactic_negative(e->data.function.args[0]);
    /* A literal Minus[x] head needs no clause: the parser reads -x as
     * Times[-1, x], and builtin_minus rewrites any explicit Minus[x] to the
     * same, so nothing reaches here still wearing that head. */
    return false;
}

static Expr* builtin_syntacticnegativeq(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) return builtin_arg_error("Internal`SyntacticNegativeQ", argc, 1, 1);
    return expr_new_symbol(syntactic_negative(res->data.function.args[0])
                           ? SYM_True : SYM_False);
}

void minus_init(void) {
    symtab_add_builtin("Minus", builtin_minus);
    symtab_get_def("Minus")->attributes |= ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED;
    symtab_set_docstring("Minus",
        "Minus[x] is the arithmetic negation of x, equivalent to -x (Times[-1, x]).\n"
        "\tListable; Minus[x, y] (any count other than one) is left unevaluated with Minus::argx.");

    symtab_add_builtin("Internal`SyntacticNegativeQ", builtin_syntacticnegativeq);
    symtab_get_def("Internal`SyntacticNegativeQ")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Internal`SyntacticNegativeQ",
        "Internal`SyntacticNegativeQ[expr] gives True if expr carries a minus sign "
        "in the leading position of its printed form: a negative number, a "
        "Rational or Complex with a negative leading part, or a Times whose first "
        "factor is one of those. A Plus is never syntactically negative, so "
        "-a - b gives False.");
}
