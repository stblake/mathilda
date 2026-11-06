/* Mathilda -- ZetaZero[k], the k-th nontrivial zero of the Riemann zeta
 * function: 1/2 + I t_k, t_k the k-th positive zero of the Hardy function
 * RiemannSiegelZ (on the critical line; every zero computed so far lies on
 * it). ZetaZero[-k] is the conjugate, the k-th zero below the real axis.
 *
 * The index is the whole argument, so an exact ZetaZero[k] stays symbolic and
 * N gives its value: the head is NHoldFirst (N must not turn k into a Real,
 * which is not an index), and N reaches the value through
 * numeric_register_held_value. ZetaZero[x] for a number x that is not a
 * nonzero integer (0, 2.5, 1.) emits ZetaZero::intnz and stays unevaluated.
 *
 * Numerics: flint_num_zeta_zero (arb's acb_dirichlet_hardy_z_zero), which
 * certifies the index -- the zeros are separated by Gram / Rosser blocks and
 * counted with Turing's method, so the k-th zero is the k-th and not merely a
 * zero near the Gram-point estimate -- and refines the zero at the requested
 * precision. A machine request is computed at 53 bits and lowered to a
 * machine Complex. Without FLINT the call stays unevaluated under N.
 *
 * Surfaces: no machine-number argument exists (the index is an exact
 * integer, and a Real index is an error), so there is no element-wise
 * NDArray / packed / Compile[] kernel to provide; the audits' EXEMPT entries
 * record this.
 *
 * Memory contract: builtins take ownership of res but never free it; they
 * return NULL (unevaluated) or an owned result.
 */
#include "zetazero.h"

#include <stdlib.h>

#include "attr.h"
#include "expr.h"
#include "flint_num_bridge.h"
#include "message.h"
#include "numeric.h"
#include "print.h"
#include "sym_names.h"
#include "symtab.h"

/* A number that is not a symbolic expression: an Integer, BigInt, Real,
 * MPFR, Rational or Complex. */
static bool zz_is_number(const Expr* e) {
    switch (e->type) {
        case EXPR_INTEGER: case EXPR_BIGINT: case EXPR_REAL:
#ifdef USE_MPFR
        case EXPR_MPFR:
#endif
            return true;
        case EXPR_FUNCTION: {
            const Expr* h = e->data.function.head;
            return h->type == EXPR_SYMBOL &&
                   (h->data.symbol.name == SYM_Rational || h->data.symbol.name == SYM_Complex);
        }
        default:
            return false;
    }
}

static bool zz_is_nonzero_integer(const Expr* e) {
    return (e->type == EXPR_INTEGER && e->data.integer != 0) || e->type == EXPR_BIGINT;
}

Expr* builtin_zetazero(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    const Expr* k = res->data.function.args[0];
    if (zz_is_number(k) && !zz_is_nonzero_integer(k)) {
        char* s = expr_to_string((Expr*)k);
        mth_message("ZetaZero", "intnz",
                    "Nonzero integer expected at position 1 in ZetaZero[%s].",
                    s ? s : "");
        free(s);
    }
    return NULL;   /* exact index: symbolic; N gives the value */
}

/* N[ZetaZero[k], p]. */
static Expr* zz_numeric(const Expr* call, NumericSpec spec) {
    if (call->type != EXPR_FUNCTION || call->data.function.arg_count != 1) return NULL;
    const Expr* k = call->data.function.args[0];
    if (!zz_is_nonzero_integer(k)) return NULL;
    bool machine = !numeric_spec_is_mpfr(spec);
    long bits = machine ? 53 : spec.bits;
    if (bits < 53) bits = 53;
    Expr* r = flint_num_zeta_zero(k, bits);
    if (r && machine) r = numeric_lower_mpfr_to_machine(r);
    return r;
}

void zetazero_init(void) {
    symtab_add_builtin("ZetaZero", builtin_zetazero);
    symtab_get_def("ZetaZero")->attributes |=
        (ATTR_LISTABLE | ATTR_NHOLDFIRST | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    numeric_register_held_value("ZetaZero", zz_numeric);
    symtab_set_docstring("ZetaZero",
        "ZetaZero[k]\n"
        "\tThe k-th nontrivial zero 1/2 + I t of the Riemann zeta function above\n"
        "\tthe real axis (k > 0), or the |k|-th below it (k < 0). Symbolic for an\n"
        "\texact k; N gives the value at any precision, with the index certified\n"
        "\t(Gram/Rosser isolation and Turing's method).");
}
