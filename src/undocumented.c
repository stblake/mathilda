/*
 * undocumented.c -- Indexed and a set of the reference system's widely used
 * Internal` utilities (GitHub issue #127).
 *
 *   Indexed[expr, i] / Indexed[expr, {i, j, ...}]
 *       a symbolic indexed component; evaluates to the part when expr is a
 *       list deep enough for the index (Indexed[{a, b, c}, 2] is b), and
 *       otherwise stays as an inert symbolic component (Indexed[x, 1]).
 *   Internal`ExceptionFreeQ[expr]
 *       False when expr contains an exceptional value anywhere (heads
 *       included): Indeterminate, an infinity (DirectedInfinity,
 *       ComplexInfinity, Infinity), Overflow[], Underflow[], or a non-finite
 *       machine real; True otherwise.
 *   Internal`PositiveIntegerQ / Internal`NonNegativeIntegerQ
 *       syntactic integer tests (no evaluation of symbolic arguments).
 *   Internal`RealValuedNumberQ[x]
 *       True for an Integer, Rational or Real (a number with no imaginary
 *       part); False for Complex and non-numbers.
 *   Internal`RealValuedNumericQ[x]
 *       True for a NumericQ quantity whose numerical value is real (Pi,
 *       Sqrt[2], Log[3]); False for I, Sqrt[-2] and anything non-numeric.
 *   Internal`PartitionRagged[list, {n1, n2, ...}]
 *       splits list into consecutive sublists of lengths n1, n2, ...; the
 *       lengths must be non-negative integers summing to Length[list], else
 *       the call stays unevaluated.
 *
 * Internal`InheritedBlock lives with Block in modular.c.  All heads here are
 * structural or predicates (no numeric kernel applies).
 */
#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "common.h"
#include "eval.h"
#include "numeric.h"
#include "sym_names.h"
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

static bool is_head(const Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head &&
           e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == name;
}

static Expr* bool_expr(bool b) { return expr_new_symbol(b ? SYM_True : SYM_False); }

/* ---------------------------------------------------------------- Indexed */

/* The machine-integer index at `e`, or false. */
static bool index_value(const Expr* e, int64_t* out) {
    if (e->type != EXPR_INTEGER) return false;
    *out = e->data.integer;
    return true;
}

/* Descend one level: list[[i]] for 1 <= |i| <= Length (negative from the end),
 * or NULL when list is not a List or i is out of range. Borrowed. */
static const Expr* list_part(const Expr* list, int64_t i) {
    if (!is_head(list, SYM_List)) return NULL;
    int64_t n = (int64_t)list->data.function.arg_count;
    if (i < 0) i += n + 1;
    if (i < 1 || i > n) return NULL;
    return list->data.function.args[i - 1];
}

static Expr* builtin_indexed(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 2) return builtin_arg_error("Indexed", argc, 2, 2);
    const Expr* cur = res->data.function.args[0];
    const Expr* idx = res->data.function.args[1];
    if (!is_head(cur, SYM_List)) return NULL;     /* symbolic: stays inert */
    if (is_head(idx, SYM_List)) {
        size_t k = idx->data.function.arg_count;
        if (k == 0) return NULL;
        for (size_t j = 0; j < k; j++) {
            int64_t i;
            if (!index_value(idx->data.function.args[j], &i)) return NULL;
            cur = list_part(cur, i);
            if (!cur) return NULL;
        }
        return expr_copy((Expr*)cur);
    }
    int64_t i;
    if (!index_value(idx, &i)) return NULL;
    cur = list_part(cur, i);
    return cur ? expr_copy((Expr*)cur) : NULL;
}

/* -------------------------------------------------------- ExceptionFreeQ */

static bool exception_free(const Expr* e) {
    switch (e->type) {
    case EXPR_REAL:
        return isfinite(e->data.real) != 0;
    case EXPR_SYMBOL: {
        const char* s = e->data.symbol.name;
        return s != SYM_Indeterminate && s != SYM_ComplexInfinity && s != SYM_Infinity;
    }
    case EXPR_FUNCTION: {
        const Expr* h = e->data.function.head;
        if (h->type == EXPR_SYMBOL) {
            const char* s = h->data.symbol.name;
            if (s == SYM_DirectedInfinity || s == SYM_Overflow || s == SYM_Underflow)
                return false;
        }
        if (!exception_free(h)) return false;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (!exception_free(e->data.function.args[i])) return false;
        return true;
    }
    default:
        return true;
    }
}

static Expr* builtin_exceptionfreeq(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) return builtin_arg_error("Internal`ExceptionFreeQ", argc, 1, 1);
    return bool_expr(exception_free(res->data.function.args[0]));
}

/* ------------------------------------------------------------ predicates */

/* Sign of an integer-like expression: -1, 0, 1. */
static int int_sign(const Expr* e) {
    if (e->type == EXPR_INTEGER)
        return (e->data.integer > 0) - (e->data.integer < 0);
    return mpz_sgn(e->data.bigint);
}

static Expr* builtin_positiveintegerq(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) return builtin_arg_error("Internal`PositiveIntegerQ", argc, 1, 1);
    const Expr* a = res->data.function.args[0];
    return bool_expr(expr_is_integer_like(a) && int_sign(a) > 0);
}

static Expr* builtin_nonnegativeintegerq(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) return builtin_arg_error("Internal`NonNegativeIntegerQ", argc, 1, 1);
    const Expr* a = res->data.function.args[0];
    return bool_expr(expr_is_integer_like(a) && int_sign(a) >= 0);
}

/* A number with no imaginary part: Integer, Rational, machine or MPFR Real. */
static bool real_number(const Expr* a) {
    if (a->type == EXPR_INTEGER || a->type == EXPR_BIGINT || a->type == EXPR_REAL)
        return true;
#ifdef USE_MPFR
    if (a->type == EXPR_MPFR) return true;
#endif
    return is_head(a, SYM_Rational);
}

static Expr* builtin_realvaluednumberq(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) return builtin_arg_error("Internal`RealValuedNumberQ", argc, 1, 1);
    return bool_expr(real_number(res->data.function.args[0]));
}

static Expr* builtin_realvaluednumericq(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) return builtin_arg_error("Internal`RealValuedNumericQ", argc, 1, 1);
    const Expr* a = res->data.function.args[0];
    if (real_number(a)) return bool_expr(true);
    if (is_head(a, SYM_Complex) || !expr_is_numeric_quantity(a)) return bool_expr(false);
    /* A numeric quantity: real iff its numerical value has no imaginary part. */
    Expr* n = expr_new_function(expr_new_symbol(SYM_N), (Expr*[]){ expr_copy((Expr*)a) }, 1);
    Expr* v = evaluate(n);
    expr_free(n);
    bool real = real_number(v) && exception_free(v);
    expr_free(v);
    return bool_expr(real);
}

/* -------------------------------------------------------- PartitionRagged */

static Expr* builtin_partitionragged(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 2) return builtin_arg_error("Internal`PartitionRagged", argc, 2, 2);
    Expr* list = res->data.function.args[0];
    const Expr* lens = res->data.function.args[1];
    if (!is_head(list, SYM_List) || !is_head(lens, SYM_List)) return NULL;
    size_t n = list->data.function.arg_count, k = lens->data.function.arg_count;
    size_t total = 0;
    for (size_t j = 0; j < k; j++) {
        const Expr* l = lens->data.function.args[j];
        if (l->type != EXPR_INTEGER || l->data.integer < 0) return NULL;
        if ((uint64_t)l->data.integer > (uint64_t)(n - total)) return NULL;
        total += (size_t)l->data.integer;
    }
    if (total != n) return NULL;
    Expr** parts = malloc((k ? k : 1) * sizeof *parts);
    if (!parts) return NULL;
    size_t pos = 0;
    for (size_t j = 0; j < k; j++) {
        size_t m = (size_t)lens->data.function.args[j]->data.integer;
        Expr** elems = malloc((m ? m : 1) * sizeof *elems);
        if (!elems) {
            for (size_t q = 0; q < j; q++) expr_free(parts[q]);
            free(parts);
            return NULL;
        }
        for (size_t q = 0; q < m; q++) elems[q] = expr_copy(list->data.function.args[pos + q]);
        pos += m;
        parts[j] = expr_new_function(expr_new_symbol(SYM_List), elems, m);
        free(elems);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), parts, k);
    free(parts);
    return out;
}

/* ------------------------------------------------------------------ init */

static void reg(const char* name, BuiltinFunc f, uint32_t attrs, const char* doc) {
    symtab_add_builtin(name, f);
    symtab_get_def(name)->attributes |= attrs;
    symtab_set_docstring(name, doc);
}

void undocumented_init(void) {
    reg("Indexed", builtin_indexed, ATTR_NHOLDREST | ATTR_PROTECTED,
        "Indexed[expr, i] represents the i-th component of expr; Indexed[expr, {i, j, ...}]\n"
        "\tthe component with indices i, j, .... When expr is a list deep enough for the\n"
        "\tindices it evaluates to Part[expr, i, j, ...]; otherwise it stays symbolic.");
    reg("Internal`ExceptionFreeQ", builtin_exceptionfreeq, ATTR_PROTECTED,
        "Internal`ExceptionFreeQ[expr] gives False if expr contains Indeterminate, an\n"
        "\tinfinity (DirectedInfinity, ComplexInfinity), Overflow[], Underflow[] or a\n"
        "\tnon-finite machine real anywhere (heads included), and True otherwise.");
    reg("Internal`PositiveIntegerQ", builtin_positiveintegerq, ATTR_PROTECTED,
        "Internal`PositiveIntegerQ[expr] gives True if expr is an integer greater than 0, False otherwise.");
    reg("Internal`NonNegativeIntegerQ", builtin_nonnegativeintegerq, ATTR_PROTECTED,
        "Internal`NonNegativeIntegerQ[expr] gives True if expr is an integer >= 0, False otherwise.");
    reg("Internal`RealValuedNumberQ", builtin_realvaluednumberq, ATTR_PROTECTED,
        "Internal`RealValuedNumberQ[expr] gives True if expr is an Integer, Rational or Real\n"
        "\tnumber, and False otherwise (Complex numbers included).");
    reg("Internal`RealValuedNumericQ", builtin_realvaluednumericq, ATTR_PROTECTED,
        "Internal`RealValuedNumericQ[expr] gives True if expr is a numeric quantity (NumericQ)\n"
        "\twhose numerical value is real, e.g. Pi or Sqrt[2]; False otherwise (I, Sqrt[-2], x).");
    reg("Internal`PartitionRagged", builtin_partitionragged, ATTR_PROTECTED,
        "Internal`PartitionRagged[list, {n1, n2, ...}] splits list into consecutive sublists of\n"
        "\tlengths n1, n2, ...; the ni must be non-negative integers summing to Length[list].");
}
