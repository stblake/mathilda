/* Mathilda -- the Anger function J_nu(z) and the Weber function E_nu(z).
 *
 *   AngerJ[nu, z]   J_nu(z) = (1/Pi) Int_0^Pi Cos[nu t - z Sin[t]] dt
 *   WeberE[nu, z]   E_nu(z) = (1/Pi) Int_0^Pi Sin[nu t - z Sin[t]] dt
 *
 * Exact reductions (each checked against the integral definition at 30 digits
 * in tests/test_struvel_anger_kelvin.c):
 *   AngerJ[n, z]   = BesselJ[n, z]                       integer n
 *   AngerJ[nu, 0]  = Sinc[nu Pi]
 *   WeberE[n, z]   = (1/Pi) Sum_{k=0}^{Floor[(n-1)/2]} Gamma[k+1/2] (z/2)^(n-2k-1)
 *                    / Gamma[n+1/2-k] - StruveH[n, z]      n >= 0 (DLMF 11.10.15)
 *   WeberE[-n, z]  = (-1)^n WeberE[n, z]                   (t -> Pi - t in the integral)
 *   WeberE[nu, 0]  = nu Pi Sinc[nu Pi/2]^2 / 2              (= (1 - Cos[nu Pi])/(nu Pi))
 *
 * Numerics: machine real arguments through sf_machine_anger_weber (the series
 * for small |z|, J_nu / -Y_nu plus a Gauss-Legendre Laplace integral beyond);
 * arbitrary precision and complex arguments through two regularized 1F2's in
 * arb (flint_num_anger_weber). Both shared with NDArray / Compile[].
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "angerweber.h"

#include "symtab.h"
#include "attr.h"
#include "eval.h"
#include "common.h"
#include "numeric.h"             /* expr_is_numeric_quantity */
#include "flint_num_bridge.h"
#include "sf_machine.h"
#include "sf_dispatch2.h"

#include <gmp.h>
#include <stdbool.h>
#include <stdlib.h>

static Expr* mk_int(int64_t v) { return expr_new_integer(v); }
static Expr* mk_fn2(const char* name, Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(name), (Expr*[]){ a, b }, 2);
}
static Expr* aw_eval(Expr* call) {
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

static Expr* aw_arb(const Expr* nu, const Expr* z, int which) {
    return flint_num_anger_weber(nu, z, which);
}

/* E_n(z), n >= 0, as the DLMF 11.10.15 polynomial minus StruveH[n, z]. The
 * coefficient of z^(n-2k-1)/Pi is Gamma(k+1/2)/(Gamma(n-k+1/2) 2^(n-2k-1))
 * = 2 / Prod_{j=k}^{n-k-1} (2j+1), built from the top k down so each product
 * extends the previous by (2k+1)(2n-2k-1). */
static Expr* weber_integer(const mpz_t n_mpz, const Expr* z) {
    if (!mpz_fits_slong_p(n_mpz)) return NULL;
    long n = mpz_get_si(n_mpz);
    Expr* struve = mk_fn2("StruveH", expr_bigint_normalize(expr_new_bigint_from_mpz(n_mpz)), expr_copy((Expr*)z));
    if (n == 0) return aw_eval(mk_fn2("Times", mk_int(-1), struve));
    long kmax = (n - 1) / 2;
    Expr** t = malloc(sizeof(Expr*) * (size_t)(kmax + 2));
    if (!t) { expr_free(struve); return NULL; }
    mpz_t P;
    mpz_init_set_ui(P, 1);
    for (long j = kmax; j <= n - kmax - 1; j++) mpz_mul_ui(P, P, (unsigned long)(2 * j + 1));
    for (long k = kmax; k >= 0; k--) {
        if (k < kmax) {
            mpz_mul_ui(P, P, (unsigned long)(2 * k + 1));
            mpz_mul_ui(P, P, (unsigned long)(2 * n - 2 * k - 1));
        }
        t[k] = expr_new_function(expr_new_symbol("Times"), (Expr*[]){
            mk_int(2),
            mk_fn2("Power", expr_bigint_normalize(expr_new_bigint_from_mpz(P)), mk_int(-1)),
            mk_fn2("Power", expr_new_symbol("Pi"), mk_int(-1)),
            mk_fn2("Power", expr_copy((Expr*)z), mk_int(n - 2 * k - 1)) }, 4);
    }
    mpz_clear(P);
    t[kmax + 1] = mk_fn2("Times", mk_int(-1), struve);
    Expr* sum = expr_new_function(expr_new_symbol("Plus"), t, (size_t)(kmax + 2));
    free(t);
    return aw_eval(sum);
}

static Expr* anger_weber(Expr* res, int which) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* nu = res->data.function.args[0];
    Expr* z  = res->data.function.args[1];

    bool handled;
    Expr* v = sf2_numeric(nu, z, which, sf_machine_anger_weber, aw_arb, &handled);
    if (handled) return v;

    if (sf2_exact_integer(nu)) {
        if (which == 0)                       /* AngerJ[n, z] = BesselJ[n, z] */
            return aw_eval(mk_fn2("BesselJ", expr_copy(nu), expr_copy(z)));
        mpz_t n;
        mpz_init(n);
        expr_to_mpz(nu, n);
        bool neg = mpz_sgn(n) < 0, odd = mpz_odd_p(n);
        mpz_abs(n, n);
        Expr* e = weber_integer(n, z);
        mpz_clear(n);
        if (e && neg && odd) e = aw_eval(mk_fn2("Times", mk_int(-1), e));
        return e;
    }
    if (sf2_exact_zero(z)) {
        /* An exact numeric order (nonzero; 0 is an integer, handled above):
         * the elementary values Sin[nu Pi]/(nu Pi) and (1 - Cos[nu Pi])/(nu Pi),
         * which evaluate (AngerJ[1/2, 0] = 2/Pi); Sinc stays for a symbolic nu. */
        if (expr_is_numeric_quantity(nu)) {
            Expr* inv = mk_fn2("Power", mk_fn2("Times", expr_copy(nu), expr_new_symbol("Pi")),
                               mk_int(-1));
            Expr* npi = mk_fn2("Times", expr_copy(nu), expr_new_symbol("Pi"));
            Expr* top = which == 0
                ? expr_new_function(expr_new_symbol("Sin"), (Expr*[]){ npi }, 1)
                : mk_fn2("Plus", mk_int(1), mk_fn2("Times", mk_int(-1),
                      expr_new_function(expr_new_symbol("Cos"), (Expr*[]){ npi }, 1)));
            return aw_eval(mk_fn2("Times", top, inv));
        }
        Expr* npi = mk_fn2("Times", expr_copy(nu), expr_new_symbol("Pi"));
        if (which == 0)                       /* AngerJ[nu, 0] = Sinc[nu Pi] */
            return aw_eval(expr_new_function(expr_new_symbol("Sinc"), (Expr*[]){ npi }, 1));
        /* WeberE[nu, 0] = nu Pi Sinc[nu Pi/2]^2 / 2 */
        Expr* half = mk_fn2("Times", mk_fn2("Rational", mk_int(1), mk_int(2)), expr_copy(npi));
        Expr* sinc = expr_new_function(expr_new_symbol("Sinc"), (Expr*[]){ half }, 1);
        return aw_eval(expr_new_function(expr_new_symbol("Times"), (Expr*[]){
            mk_fn2("Rational", mk_int(1), mk_int(2)), npi,
            mk_fn2("Power", sinc, mk_int(2)) }, 3));
    }
    return NULL;
}

Expr* builtin_angerj(Expr* res) { return anger_weber(res, 0); }
Expr* builtin_webere(Expr* res) { return anger_weber(res, 1); }

void angerweber_init(void) {
    symtab_add_builtin("AngerJ", builtin_angerj);
    symtab_get_def("AngerJ")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    symtab_set_docstring("AngerJ",
        "AngerJ[nu, z]\n"
        "The Anger function J_nu(z) = (1/Pi) Integrate[Cos[nu t - z Sin[t]], {t, 0, Pi}].\n"
        "An integer order gives BesselJ[n, z]; AngerJ[nu, 0] = Sinc[nu Pi].\n"
        "Machine, arbitrary-precision and complex arguments; fast paths over packed\n"
        "arrays and in Compile[]. D[AngerJ[nu, z], z] = (AngerJ[nu-1, z] - AngerJ[nu+1, z])/2.");
    symtab_add_builtin("WeberE", builtin_webere);
    symtab_get_def("WeberE")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    symtab_set_docstring("WeberE",
        "WeberE[nu, z]\n"
        "The Weber function E_nu(z) = (1/Pi) Integrate[Sin[nu t - z Sin[t]], {t, 0, Pi}].\n"
        "An integer order reduces to a polynomial minus StruveH[n, z] (WeberE[0, z] =\n"
        "-StruveH[0, z], WeberE[-n, z] = (-1)^n WeberE[n, z]); WeberE[nu, 0] =\n"
        "nu Pi Sinc[nu Pi/2]^2/2. Machine, arbitrary-precision and complex arguments;\n"
        "fast paths over packed arrays and in Compile[].\n"
        "D[WeberE[nu, z], z] = (WeberE[nu-1, z] - WeberE[nu+1, z])/2.");
}
