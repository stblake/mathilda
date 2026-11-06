/* Mathilda -- the zeros of the Airy functions:
 *
 *   AiryAiZero[k]        a_k,  the k-th zero of AiryAi
 *   AiryAiPrimeZero[k]   a'_k, the k-th zero of AiryAiPrime
 *   AiryBiZero[k]        b_k,  the k-th zero of AiryBi
 *   AiryBiPrimeZero[k]   b'_k, the k-th zero of AiryBiPrime
 *
 * All four sequences are real and negative, decreasing in k = 1, 2, ...
 * (DLMF 9.9). The index is the whole argument, so an exact call stays
 * symbolic and N gives the value: the heads are NHoldFirst and N reaches them
 * through numeric_register_held_value. A non-integer or non-positive index
 * stays unevaluated.
 *
 * Numerics: flint_num_airy_zero (arb's arb_hypgeom_airy_zero -- the
 * asymptotic expansion of the zero in (3 pi (4k - 1)/8)^(2/3), DLMF 9.9.6-9.9.9,
 * refined by Newton with a certified enclosure), at the requested precision;
 * a machine request at 53 bits, lowered to a machine Real. Without FLINT the
 * call stays unevaluated under N.
 *
 * Surfaces: the argument is an exact index (a Real index is not evaluated),
 * so there is no machine-number argument to map an NDArray / packed / Compile[]
 * kernel over; the audits' EXEMPT entries record this.
 *
 * Memory contract: builtins take ownership of res but never free it; they
 * return NULL (unevaluated) or an owned result.
 */
#include "airyzero.h"

#include "attr.h"
#include "expr.h"
#include "flint_num_bridge.h"
#include "numeric.h"
#include "sym_names.h"
#include "sym_intern.h"
#include "symtab.h"

/* The index -> arb selector: 0 Ai, 1 Ai', 2 Bi, 3 Bi'. */
static const char* const AZ_NAMES[4] = {
    "AiryAiZero", "AiryAiPrimeZero", "AiryBiZero", "AiryBiPrimeZero"
};

static int az_which(const Expr* call) {
    if (call->type != EXPR_FUNCTION || call->data.function.head->type != EXPR_SYMBOL)
        return -1;
    const char* h = call->data.function.head->data.symbol.name;
    for (int i = 0; i < 4; i++)
        if (h == intern_symbol(AZ_NAMES[i])) return i;
    return -1;
}

/* Exact calls stay symbolic (N supplies the value); nothing else evaluates. */
Expr* builtin_airyzero(Expr* res) {
    (void)res;
    return NULL;
}

/* N[AiryAiZero[k], p] and its three siblings. */
static Expr* az_numeric(const Expr* call, NumericSpec spec) {
    if (call->type != EXPR_FUNCTION || call->data.function.arg_count != 1) return NULL;
    int which = az_which(call);
    if (which < 0) return NULL;
    const Expr* k = call->data.function.args[0];
    if (!((k->type == EXPR_INTEGER && k->data.integer > 0) ||
          (k->type == EXPR_BIGINT && mpz_sgn(k->data.bigint) > 0)))
        return NULL;
    bool machine = !numeric_spec_is_mpfr(spec);
    long bits = machine ? 53 : spec.bits;
    if (bits < 53) bits = 53;
    Expr* r = flint_num_airy_zero(which, k, bits);
    if (r && machine) r = numeric_lower_mpfr_to_machine(r);
    return r;
}

void airyzero_init(void) {
    static const char* const DOCS[4] = {
        "AiryAiZero[k]\n"
        "\tThe k-th zero of AiryAi[x] (k = 1, 2, ...), a negative real.\n"
        "\tSymbolic for an exact k; N gives the value at any precision.",
        "AiryAiPrimeZero[k]\n"
        "\tThe k-th zero of AiryAiPrime[x] (k = 1, 2, ...), a negative real.\n"
        "\tSymbolic for an exact k; N gives the value at any precision.",
        "AiryBiZero[k]\n"
        "\tThe k-th zero of AiryBi[x] (k = 1, 2, ...), a negative real.\n"
        "\tSymbolic for an exact k; N gives the value at any precision.",
        "AiryBiPrimeZero[k]\n"
        "\tThe k-th zero of AiryBiPrime[x] (k = 1, 2, ...), a negative real.\n"
        "\tSymbolic for an exact k; N gives the value at any precision.",
    };
    for (int i = 0; i < 4; i++) {
        symtab_add_builtin(AZ_NAMES[i], builtin_airyzero);
        symtab_get_def(AZ_NAMES[i])->attributes |=
            (ATTR_LISTABLE | ATTR_NHOLDFIRST | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
        numeric_register_held_value(AZ_NAMES[i], az_numeric);
        symtab_set_docstring(AZ_NAMES[i], DOCS[i]);
    }
}
