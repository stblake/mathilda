/* The machine-scalar gate shared by the special-function builtins.
 * See sf_scalar.h for the contract. */
#include "sf_scalar.h"

#include <gmp.h>
#include "sym_names.h"
#include "numeric.h"

bool sf_machine_leaf(const Expr* e, double* out, bool* inexact) {
    switch (e->type) {
        case EXPR_REAL:
            *out = e->data.real;
            *inexact = true;
            return true;
        case EXPR_INTEGER:
            *out = (double)e->data.integer;
            *inexact = false;
            return true;
        case EXPR_BIGINT:
            *out = numeric_mpz_to_double(e->data.bigint);
            *inexact = false;
            return true;
        case EXPR_FUNCTION: {
            /* Rational[p, q] with integer parts. */
            const Expr* h = e->data.function.head;
            if (h->type != EXPR_SYMBOL || h->data.symbol.name != SYM_Rational ||
                e->data.function.arg_count != 2)
                return false;
            const Expr* p = e->data.function.args[0];
            const Expr* q = e->data.function.args[1];
            /* one rounding of p/q (numeric_machine_value), not the
             * double-rounded -- or inf/inf = NaN -- quotient of two doubles */
            if (!expr_is_integer_like(p) || !expr_is_integer_like(q)) return false;
            if (!numeric_machine_value(e, out)) return false;
            *inexact = false;
            return true;
        }
        default:
            return false;   /* MPFR, symbol, string, ... */
    }
}

bool sf_machine_scalar_args(const Expr* res, SFMachineArgs* out) {
    if (res->type != EXPR_FUNCTION) return false;
    size_t n = res->data.function.arg_count;
    if (n == 0 || n > SF_MACHINE_MAXARGS) return false;
    bool any_inexact = false;
    out->n = (int)n;
    out->any_complex = false;
    for (size_t i = 0; i < n; i++) {
        const Expr* a = res->data.function.args[i];
        bool ix = false;
        if (sf_machine_leaf(a, &out->re[i], &ix)) {
            out->im[i] = 0.0;
            out->is_complex[i] = false;
        } else if (a->type == EXPR_FUNCTION &&
                   a->data.function.head->type == EXPR_SYMBOL &&
                   a->data.function.head->data.symbol.name == SYM_Complex &&
                   a->data.function.arg_count == 2) {
            bool ix_re = false, ix_im = false;
            if (!sf_machine_leaf(a->data.function.args[0], &out->re[i], &ix_re) ||
                !sf_machine_leaf(a->data.function.args[1], &out->im[i], &ix_im))
                return false;
            ix = ix_re || ix_im;
            out->is_complex[i] = true;
            out->any_complex = true;
        } else {
            return false;
        }
        out->is_exact[i] = !ix;
        any_inexact |= ix;
    }
    return any_inexact;
}

Expr* sf_machine_complex_result(double re, double im) {
    /* One shared `Complex` head node, handed out by refcount. expr_new_symbol
     * would allocate a node and re-hash the name through intern_symbol on every
     * call, which is a measurable share of a ~0.2 us machine-complex call.
     * Symbol nodes are immutable apart from the benign cached-def pointer, so
     * sharing one is safe (the same idiom as expr_copy everywhere else). */
    static Expr* complex_head = NULL;
    if (!complex_head) complex_head = expr_new_symbol(SYM_Complex);
    Expr* args[2] = { expr_new_real(re), expr_new_real(im) };
    return expr_new_function(expr_copy(complex_head), args, 2);
}
