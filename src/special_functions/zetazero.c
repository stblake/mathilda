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

#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "attr.h"
#include "eval.h"
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

/* True when Im[z] > t, decided by the comparison builtins (z the certified
 * zero, t a real number); *ok cleared when undecided. */
static bool zz_im_above(const Expr* z, const Expr* t, bool* ok) {
    Expr* im = expr_new_function(expr_new_symbol(SYM_Im), (Expr*[]){ expr_copy((Expr*)z) }, 1);
    Expr* r = eval_and_free(expr_new_function(expr_new_symbol(SYM_Greater),
                  (Expr*[]){ im, expr_copy((Expr*)t) }, 2));
    bool above = r && r->type == EXPR_SYMBOL && r->data.symbol.name == SYM_True;
    bool below = r && r->type == EXPR_SYMBOL && r->data.symbol.name == SYM_False;
    if (!above && !below) *ok = false;
    expr_free(r);
    return above;
}

/* N[ZetaZero[k, t]]: the k-th zero above height t. The index of the first
 * zero above t starts at the Riemann-von Mangoldt count N(t) ~ theta(t)/pi + 1
 * and is corrected against the certified zeros, so the answer is the
 * certified ZetaZero[j + k - 1] (ZetaZero#1 of the 2026-11 audit). */
static Expr* zz_numeric_above(const Expr* k, const Expr* t, long bits) {
    if (k->type != EXPR_INTEGER || k->data.integer < 1) return NULL;
    Expr* td = numericalize((Expr*)t, numeric_machine_spec());
    bool real = td && td->type == EXPR_REAL;
    double tv = real ? td->data.real : 0.0;
    expr_free(td);
    if (!real || !(tv < 1.0e15)) return NULL;
    int64_t j = 1;
    if (tv > 14.0) {
        double th = 0.5 * tv * log(tv / (2.0 * 3.14159265358979323846)) - 0.5 * tv
                  - 3.14159265358979323846 / 8.0 + 1.0 / (48.0 * tv);
        double nt = th / 3.14159265358979323846 + 1.0;
        j = (int64_t)floor(nt);
        if (j < 1) j = 1;
    }
    /* j: the first index whose zero lies above t */
    for (int it = 0; it < 2000; it++) {
        Expr* jz = expr_new_integer(j);
        Expr* z = flint_num_zeta_zero(jz, bits);
        expr_free(jz);
        if (!z) return NULL;
        bool ok = true;
        bool above = zz_im_above(z, t, &ok);
        expr_free(z);
        if (!ok) return NULL;
        if (above) {
            if (j == 1) break;
            Expr* pz = expr_new_integer(j - 1);
            Expr* z2 = flint_num_zeta_zero(pz, bits);
            expr_free(pz);
            if (!z2) return NULL;
            bool prev_above = zz_im_above(z2, t, &ok);
            expr_free(z2);
            if (!ok) return NULL;
            if (!prev_above) break;
            j--;
        } else {
            j++;
        }
        if (it == 1999) return NULL;
    }
    Expr* idx = expr_new_integer(j + k->data.integer - 1);
    Expr* r = flint_num_zeta_zero(idx, bits);
    expr_free(idx);
    return r;
}

/* N[ZetaZero[k], p]. */
static Expr* zz_numeric(const Expr* call, NumericSpec spec) {
    if (call->type == EXPR_FUNCTION && call->data.function.arg_count == 2) {
        bool machine = !numeric_spec_is_mpfr(spec);
        long bits = machine ? 53 : spec.bits;
        if (bits < 53) bits = 53;
        Expr* r = zz_numeric_above(call->data.function.args[0], call->data.function.args[1], bits);
        if (r && machine) r = numeric_lower_mpfr_to_machine(r);
        return r;
    }
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
        "ZetaZero[k]\nZetaZero[k, t]\n"
        "\tThe k-th nontrivial zero 1/2 + I t of the Riemann zeta function above\n"
        "\tthe real axis (k > 0), or the |k|-th below it (k < 0). Symbolic for an\n"
        "\texact k; N gives the value at any precision, with the index certified\n"
        "\t(Gram/Rosser isolation and Turing's method). ZetaZero[k, t] is the k-th\n"
        "\tzero above height t.");
}
