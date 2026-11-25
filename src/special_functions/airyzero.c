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

#include <math.h>
#include <stdint.h>

#include "attr.h"
#include "eval.h"
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

#ifdef USE_MPFR
/* The asymptotic zero (DLMF 9.9.6-9.9.9) for an index past arb's 62-bit
 * range: -T(t) for Ai, Bi and -U(t) for Ai', Bi', t = 3 Pi (4k - c)/8 with
 * c = 1 (Ai, Bi'), 3 (Ai', Bi),
 *   T(t) = t^(2/3) (1 + 5/48 t^-2 - 5/36 t^-4 + 77125/82944 t^-6 - ...),
 *   U(t) = t^(2/3) (1 - 7/48 t^-2 + 35/288 t^-4 - 181223/207360 t^-6 + ...).
 * At k >= 2^62 the omitted t^-8 term is below 2^-500 relative, so the value
 * is exact to any precision up to that (AiryAiPrimeZero#1, AiryBiZero#1 of
 * the 2026-11 audit). NULL past that precision. */
static Expr* az_asymptotic(int which, const Expr* k, long bits) {
    if (bits > 480) return NULL;
    mpfr_prec_t wp = (mpfr_prec_t)bits + 32;
    mpfr_t t, u, v, w;
    mpfr_inits2(wp, t, u, v, w, (mpfr_ptr)0);
    mpz_t kz;
    expr_to_mpz(k, kz);
    mpfr_set_z(t, kz, MPFR_RNDN);
    mpz_clear(kz);
    mpfr_mul_ui(t, t, 4, MPFR_RNDN);
    bool c1 = (which == 0 || which == 3);
    mpfr_sub_ui(t, t, c1 ? 1 : 3, MPFR_RNDN);
    mpfr_const_pi(u, MPFR_RNDN);
    mpfr_mul(t, t, u, MPFR_RNDN);
    mpfr_mul_ui(t, t, 3, MPFR_RNDN);
    mpfr_div_ui(t, t, 8, MPFR_RNDN);              /* t */
    mpfr_sqr(u, t, MPFR_RNDN);
    mpfr_ui_div(u, 1, u, MPFR_RNDN);              /* u = t^-2 */
    bool T = (which == 0 || which == 2);
    /* Horner in u: 1 + a1 u + a2 u^2 + a3 u^3 */
    mpfr_set_si(v, T ? 77125 : -181223, MPFR_RNDN);
    mpfr_div_ui(v, v, T ? 82944 : 207360, MPFR_RNDN);
    mpfr_mul(v, v, u, MPFR_RNDN);
    mpfr_set_si(w, T ? -5 : 35, MPFR_RNDN);
    mpfr_div_ui(w, w, T ? 36 : 288, MPFR_RNDN);
    mpfr_add(v, v, w, MPFR_RNDN);
    mpfr_mul(v, v, u, MPFR_RNDN);
    mpfr_set_si(w, T ? 5 : -7, MPFR_RNDN);
    mpfr_div_ui(w, w, 48, MPFR_RNDN);
    mpfr_add(v, v, w, MPFR_RNDN);
    mpfr_mul(v, v, u, MPFR_RNDN);
    mpfr_add_ui(v, v, 1, MPFR_RNDN);              /* the series */
    mpfr_set_ui(w, 2, MPFR_RNDN);
    mpfr_div_ui(w, w, 3, MPFR_RNDN);
    mpfr_pow(t, t, w, MPFR_RNDN);                  /* t^(2/3) */
    mpfr_mul(t, t, v, MPFR_RNDN);
    mpfr_neg(t, t, MPFR_RNDN);
    Expr* r = expr_new_mpfr_bits((mpfr_prec_t)bits);
    mpfr_set(r->data.mpfr, t, MPFR_RNDN);
    mpfr_clears(t, u, v, w, (mpfr_ptr)0);
    return r;
}
#endif

/* The k-th zero at `bits`, the asymptotic form past arb's index range. */
static Expr* az_zero(int which, const Expr* k, long bits) {
    Expr* r = flint_num_airy_zero(which, k, bits);
#ifdef USE_MPFR
    /* arb declines past a 62-bit index; from k = 2^40 the omitted term is
     * below 2^-330, so the asymptotic form is exact up to 300 bits */
    if (!r && (k->type == EXPR_BIGINT ||
               (k->type == EXPR_INTEGER && k->data.integer >= (INT64_C(1) << 40) && bits <= 300)))
        r = az_asymptotic(which, k, bits);
#endif
    return r;
}

/* True when the zero z lies below x0 (both real numbers); *ok cleared when
 * the comparison does not decide. */
static bool az_below(const Expr* z, const Expr* x0, bool* ok) {
    Expr* r = eval_and_free(expr_new_function(expr_new_symbol(SYM_Less),
                  (Expr*[]){ expr_copy((Expr*)z), expr_copy((Expr*)x0) }, 2));
    bool lt = r && r->type == EXPR_SYMBOL && r->data.symbol.name == SYM_True;
    bool ge = r && r->type == EXPR_SYMBOL && r->data.symbol.name == SYM_False;
    if (!lt && !ge) *ok = false;
    expr_free(r);
    return lt;
}

/* The k-th zero less than x0 (AiryBiPrimeZero#1 of the 2026-11 audit): the
 * index j of the first zero below x0 starts at the inverse of the leading
 * asymptotic, (4j - c) 3 Pi/8 ~ (-x0)^(3/2), and is corrected against the
 * computed zeros. */
static Expr* az_below_x0(int which, const Expr* k, const Expr* x0, long bits) {
    if (k->type != EXPR_INTEGER || k->data.integer < 1) return NULL;
    Expr* xd = numericalize((Expr*)x0, numeric_machine_spec());
    bool real = xd && xd->type == EXPR_REAL;
    double xv = real ? xd->data.real : 0.0;
    expr_free(xd);
    if (!real || !(xv > -1.0e12)) return NULL;
    double c = (which == 0 || which == 3) ? 1.0 : 3.0;
    double j = 1.0;
    if (xv < 0.0) j = floor((pow(-xv, 1.5) * 8.0 / (3.0 * 3.14159265358979323846) + c) / 4.0);
    if (j < 1.0) j = 1.0;
    for (int it = 0; it < 64; it++) {
        Expr* je = expr_new_integer((int64_t)j);
        Expr* z = az_zero(which, je, bits);
        expr_free(je);
        if (!z) return NULL;
        bool ok = true;
        bool below = az_below(z, x0, &ok);
        expr_free(z);
        if (!ok) return NULL;
        if (below) {
            if (j == 1.0) break;
            Expr* pe = expr_new_integer((int64_t)j - 1);
            Expr* zp = az_zero(which, pe, bits);
            expr_free(pe);
            if (!zp) return NULL;
            bool prev_below = az_below(zp, x0, &ok);
            expr_free(zp);
            if (!ok) return NULL;
            if (!prev_below) break;
            j -= 1.0;
        } else {
            j += 1.0;
        }
        if (it == 63) return NULL;
    }
    Expr* idx = expr_new_integer((int64_t)j + k->data.integer - 1);
    Expr* r = az_zero(which, idx, bits);
    expr_free(idx);
    return r;
}

/* N[AiryAiZero[k], p] and its three siblings; N[AiryAiZero[k, x0], p] the
 * k-th zero less than x0. */
static Expr* az_numeric(const Expr* call, NumericSpec spec) {
    if (call->type != EXPR_FUNCTION) return NULL;
    size_t na = call->data.function.arg_count;
    if (na != 1 && na != 2) return NULL;
    int which = az_which(call);
    if (which < 0) return NULL;
    const Expr* k = call->data.function.args[0];
    if (!((k->type == EXPR_INTEGER && k->data.integer > 0) ||
          (k->type == EXPR_BIGINT && mpz_sgn(k->data.bigint) > 0)))
        return NULL;
    bool machine = !numeric_spec_is_mpfr(spec);
    long bits = machine ? 53 : spec.bits;
    if (bits < 53) bits = 53;
    Expr* r = na == 2 ? az_below_x0(which, k, call->data.function.args[1], bits)
                      : az_zero(which, k, bits);
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
