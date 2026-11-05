/* Mathilda -- the Kelvin functions ber, bei, ker, kei.
 *
 *   KelvinBer[nu, z] + I KelvinBei[nu, z] = BesselJ[nu, z Exp[3 Pi I/4]]
 *   KelvinKer[nu, z] + I KelvinKei[nu, z] = Exp[-nu Pi I/2] BesselK[nu, z Exp[Pi I/4]]
 * for z > 0 (DLMF 10.61.1-2), continued analytically with the only branch cut
 * along the negative real axis. The one-argument forms are order 0 and
 * evaluate to it: KelvinBer[z] -> KelvinBer[0, z].
 *
 * Exact reductions:
 *   head[-n, z] = (-1)^n head[n, z]          integer n > 0 (DLMF 10.61.6)
 *   KelvinBer[0, 0] = 1, KelvinBer[nu, 0] = 0 (nu > 0), KelvinBei[nu, 0] = 0 (nu >= 0)
 *   KelvinKer[0, 0] = Infinity, KelvinKei[0, 0] = -Pi/4,
 *   KelvinKer[1, 0] = KelvinKei[1, 0] = ComplexInfinity
 *
 * Numerics: machine real arguments through sf_machine_kelvin (the shared complex
 * Bessel kernels at the rotated argument); arbitrary precision and complex
 * arguments through arb (flint_num_kelvin). Both shared with NDArray / Compile[].
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "kelvin.h"

#include "symtab.h"
#include "attr.h"
#include "eval.h"
#include "common.h"
#include "flint_num_bridge.h"
#include "sf_machine.h"
#include "sf_dispatch2.h"

#include <gmp.h>
#include <stdbool.h>

static const char* const KELVIN_NAME[4] = { "KelvinBer", "KelvinBei", "KelvinKer", "KelvinKei" };

static Expr* mk_int(int64_t v) { return expr_new_integer(v); }
static Expr* mk_fn2(const char* name, Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(name), (Expr*[]){ a, b }, 2);
}
static Expr* kv_eval(Expr* call) {
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

static Expr* kelvin(Expr* res, int which) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc == 1)                                /* head[z] -> head[0, z] */
        return mk_fn2(KELVIN_NAME[which], mk_int(0), expr_copy(res->data.function.args[0]));
    if (argc != 2) return NULL;
    Expr* nu = res->data.function.args[0];
    Expr* z  = res->data.function.args[1];

    /* ker_0 has a logarithmic singularity at 0: Infinity at a machine 0. too
     * (the kernels decline there, and arb's non-finite ball would leave it
     * unevaluated). */
    if (which == 2 && z->type == EXPR_REAL && z->data.real == 0.0) {
        double nuv;
        if (common_machine_real_value(nu, &nuv) && nuv == 0.0)
            return expr_new_symbol("Infinity");
    }
    bool handled;
    Expr* v = sf2_numeric(nu, z, which, sf_machine_kelvin, flint_num_kelvin, &handled);
    if (handled) return v;

    if (sf2_exact_integer(nu)) {
        mpz_t n;
        mpz_init(n);
        expr_to_mpz(nu, n);
        if (mpz_sgn(n) < 0) {                     /* head[-n, z] = (-1)^n head[n, z] */
            bool odd = mpz_odd_p(n);
            mpz_neg(n, n);
            Expr* pos = mk_fn2(KELVIN_NAME[which],
                               expr_bigint_normalize(expr_new_bigint_from_mpz(n)),
                               expr_copy(z));
            mpz_clear(n);
            return kv_eval(odd ? mk_fn2("Times", mk_int(-1), pos) : pos);
        }
        mpz_clear(n);
    }
    if (sf2_exact_zero(z)) {
        double nuv;
        if (!common_machine_real_value(nu, &nuv)) return NULL;
        bool n0 = sf2_exact_zero(nu);
        bool n1 = nu->type == EXPR_INTEGER && nu->data.integer == 1;
        switch (which) {
        case 0: if (n0) return mk_int(1);
                if (nuv > 0.0) return mk_int(0);
                return NULL;
        case 1: return nuv >= 0.0 ? mk_int(0) : NULL;
        case 2: if (n0) return expr_new_symbol("Infinity");
                if (n1) return expr_new_symbol("ComplexInfinity");
                return NULL;
        default:
                if (n0) return kv_eval(mk_fn2("Times",
                              mk_fn2("Rational", mk_int(-1), mk_int(4)), expr_new_symbol("Pi")));
                if (n1) return expr_new_symbol("ComplexInfinity");
                return NULL;
        }
    }
    return NULL;
}

Expr* builtin_kelvinber(Expr* res) { return kelvin(res, 0); }
Expr* builtin_kelvinbei(Expr* res) { return kelvin(res, 1); }
Expr* builtin_kelvinker(Expr* res) { return kelvin(res, 2); }
Expr* builtin_kelvinkei(Expr* res) { return kelvin(res, 3); }

void kelvin_init(void) {
    /* Literal names: the audit tools (tools/numeric_coverage.py and friends)
     * read the registry out of symtab_add_builtin("...") call sites. */
    symtab_add_builtin("KelvinBer", builtin_kelvinber);
    symtab_add_builtin("KelvinBei", builtin_kelvinbei);
    symtab_add_builtin("KelvinKer", builtin_kelvinker);
    symtab_add_builtin("KelvinKei", builtin_kelvinkei);
    static const char* const DOC[4] = {
        "KelvinBer[z] or KelvinBer[nu, z]\n"
        "The Kelvin function ber_nu(z): KelvinBer[nu, x] + I KelvinBei[nu, x] =\n"
        "BesselJ[nu, x Exp[3 Pi I/4]] for x > 0, continued analytically (cut along the\n"
        "negative real axis). KelvinBer[z] = KelvinBer[0, z]; KelvinBer[-n, z] =\n"
        "(-1)^n KelvinBer[n, z]. Machine, arbitrary-precision and complex arguments;\n"
        "fast paths over packed arrays and in Compile[].",
        "KelvinBei[z] or KelvinBei[nu, z]\n"
        "The Kelvin function bei_nu(z): KelvinBer[nu, x] + I KelvinBei[nu, x] =\n"
        "BesselJ[nu, x Exp[3 Pi I/4]] for x > 0, continued analytically (cut along the\n"
        "negative real axis). KelvinBei[z] = KelvinBei[0, z]; KelvinBei[-n, z] =\n"
        "(-1)^n KelvinBei[n, z]. Machine, arbitrary-precision and complex arguments;\n"
        "fast paths over packed arrays and in Compile[].",
        "KelvinKer[z] or KelvinKer[nu, z]\n"
        "The Kelvin function ker_nu(z): KelvinKer[nu, x] + I KelvinKei[nu, x] =\n"
        "Exp[-nu Pi I/2] BesselK[nu, x Exp[Pi I/4]] for x > 0, continued analytically\n"
        "(cut along the negative real axis). KelvinKer[z] = KelvinKer[0, z];\n"
        "KelvinKer[-n, z] = (-1)^n KelvinKer[n, z]. Machine, arbitrary-precision and\n"
        "complex arguments; fast paths over packed arrays and in Compile[].",
        "KelvinKei[z] or KelvinKei[nu, z]\n"
        "The Kelvin function kei_nu(z): KelvinKer[nu, x] + I KelvinKei[nu, x] =\n"
        "Exp[-nu Pi I/2] BesselK[nu, x Exp[Pi I/4]] for x > 0, continued analytically\n"
        "(cut along the negative real axis). KelvinKei[z] = KelvinKei[0, z];\n"
        "KelvinKei[-n, z] = (-1)^n KelvinKei[n, z]; KelvinKei[0, 0] = -Pi/4. Machine,\n"
        "arbitrary-precision and complex arguments; fast paths over packed arrays and\n"
        "in Compile[]." };
    for (int i = 0; i < 4; i++) {
        symtab_get_def(KELVIN_NAME[i])->attributes |=
            (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
        symtab_set_docstring(KELVIN_NAME[i], DOC[i]);
    }
}
