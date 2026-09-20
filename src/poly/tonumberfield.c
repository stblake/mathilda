/* tonumberfield.c — ToNumberField builtin (argument-form dispatch).
 *
 * See tonumberfield.h. All field/primitive-element computation is delegated to
 * flint_qqbar.c; this file only classifies the argument forms and assembles the
 * result. Each produced AlgebraicNumber is re-canonicalised by the evaluator.
 */

#include "tonumberfield.h"

#include "flint_qqbar.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"

#include <stdlib.h>
#include <string.h>

static int is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION && e->data.function.head &&
           e->data.function.head->type == EXPR_SYMBOL &&
           strcmp(e->data.function.head->data.symbol.name, "List") == 0;
}

static int is_symbol_named(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL &&
           strcmp(e->data.symbol.name, name) == 0;
}

/* ToNumberField[{a1..ak}] / [{a1..ak}, All]: common field via primitive element. */
static Expr* common_form(const Expr* list, int smallest) {
    size_t n = list->data.function.arg_count;
    if (n == 0) return NULL;
    const Expr** as = malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) as[i] = list->data.function.args[i];
    Expr* r = flint_qqbar_to_number_field_common((const Expr* const*)as, n, smallest);
    free(as);
    return r;
}

/* ToNumberField[{a1..ak}, theta]: express each ai in Q(theta). */
static Expr* list_in_field(const Expr* list, const Expr* theta) {
    size_t n = list->data.function.arg_count;
    Expr** items = malloc(sizeof(Expr*) * (n ? n : 1));
    size_t built = 0;
    int good = 1;
    for (size_t i = 0; i < n; i++) {
        items[i] = flint_qqbar_to_number_field(list->data.function.args[i], theta);
        if (!items[i]) { good = 0; break; }   /* some ai not in Q(theta) */
        built++;
    }
    Expr* r = NULL;
    if (good) r = expr_new_function(expr_new_symbol(SYM_List), items, n);
    else for (size_t i = 0; i < built; i++) expr_free(items[i]);
    free(items);
    return r;
}

Expr* builtin_tonumberfield(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t n = res->data.function.arg_count;
    Expr** a = res->data.function.args;

    if (n == 1) {
        if (is_list(a[0])) return common_form(a[0], 0);       /* Automatic */
        return flint_qqbar_to_number_field_self(a[0]);        /* explicit AlgebraicNumber */
    }
    if (n == 2) {
        if (is_list(a[0])) {
            if (is_symbol_named(a[1], "Automatic")) return common_form(a[0], 0);
            if (is_symbol_named(a[1], "All"))       return common_form(a[0], 1);
            return list_in_field(a[0], a[1]);                 /* each ai in Q(theta) */
        }
        return flint_qqbar_to_number_field(a[0], a[1]);       /* a in Q(theta) */
    }
    return NULL;
}

void tonumberfield_init(void) {
    symtab_add_builtin("ToNumberField", builtin_tonumberfield);
    /* Not Listable: the {a1,..} form defines a common field over the whole list,
     * so threading element-wise would destroy its meaning. */
    symtab_get_def("ToNumberField")->attributes |= ATTR_PROTECTED;
}
