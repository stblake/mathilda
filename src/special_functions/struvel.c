/* Mathilda -- the modified Struve function L_nu(z).
 *
 *   StruveL[nu, z]   L_nu(z) = (z/2)^(nu+1) (2/(Sqrt[Pi] Gamma[nu+3/2]))
 *                              1F2(1; 3/2, nu+3/2; z^2/4)        (DLMF 11.2.2)
 *
 * Evaluation (mirrors StruveH, struveh.c):
 *   StruveL[nu, 0] = 0             exact, for a real nu > -1
 *   machine precision, real nu, z  the shared double kernel sf_machine_struve_l
 *                                  (bit-for-bit with NDArray / Compile[])
 *   anything else numeric          arb's regularized 1F2 (flint_num_struvel),
 *                                  demoted to machine numbers for machine input
 *   exact half-odd-integer nu      elementary closed form in Sinh / Cosh
 *   everything else                stays symbolic
 *
 * D[StruveL[nu, z], z] lives in calculus/deriv.c with StruveH's rule.
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "struvel.h"

#include "symtab.h"
#include "attr.h"
#include "eval.h"
#include "common.h"
#include "arithmetic.h"         /* is_rational_like */
#include "flint_num_bridge.h"   /* flint_num_struvel */
#include "sf_machine.h"         /* sf_machine_struve_l */
#include "sf_dispatch2.h"

#include <stdbool.h>
#include <stdlib.h>

static Expr* mk_int(int64_t v) { return expr_new_integer(v); }
static Expr* mk_fn1(const char* name, Expr* a) {
    return expr_new_function(expr_new_symbol(name), (Expr*[]){ a }, 1);
}
static Expr* mk_fn2(const char* name, Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(name), (Expr*[]){ a, b }, 2);
}
static Expr* mk_rat(int64_t p, int64_t q) { return mk_fn2("Rational", mk_int(p), mk_int(q)); }
static Expr* sl_eval(Expr* call) {
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

static bool sl_machine(int which, double nu, double z, double* out) {
    (void)which;
    return sf_machine_struve_l(nu, z, out);
}
static Expr* sl_arb(const Expr* nu, const Expr* z, int which) {
    (void)which;
    return flint_num_struvel(nu, z);
}

/* ---- half-odd-integer order: elementary closed forms ------------------- */

/* a_k(n + 1/2) = (n + k)! / (2^k k! (n - k)!)  (DLMF 10.49.1), unevaluated. */
static Expr* sl_hankel_a(int64_t n, int64_t k) {
    return expr_new_function(expr_new_symbol("Times"), (Expr*[]){
        mk_fn1("Factorial", mk_int(n + k)),
        mk_fn2("Power", mk_int(2), mk_int(-k)),
        mk_fn2("Power", mk_fn1("Factorial", mk_int(k)), mk_int(-1)),
        mk_fn2("Power", mk_fn1("Factorial", mk_int(n - k)), mk_int(-1)) }, 4);
}

/* Sum_{k = par, par + 2, ... <= n} a_k(n + 1/2) z^-k (no sign alternation: the
 * modified functions take e^(+-z) where the ordinary ones take e^(+-iz)). */
static Expr* sl_hankel_sum(int64_t n, int par, const Expr* z) {
    int64_t cnt = n >= par ? (n - par) / 2 + 1 : 0;
    if (cnt == 0) return mk_int(0);
    Expr** t = malloc(sizeof(Expr*) * (size_t)cnt);
    if (!t) return NULL;
    for (int64_t i = 0; i < cnt; i++) {
        int64_t k = par + 2 * i;
        t[i] = mk_fn2("Times", sl_hankel_a(n, k),
                      mk_fn2("Power", expr_copy((Expr*)z), mk_int(-k)));
    }
    Expr* s = expr_new_function(expr_new_symbol("Plus"), t, (size_t)cnt);
    free(t);
    return s;
}

/* The modified spherical-Bessel forms (DLMF 10.49.8 / 10.49.9 times
 * Sqrt(2z/pi)), with E / O the even / odd Hankel sums:
 *   I_{n+1/2}(z)    = Sqrt(2/(pi z)) (E Sinh z - O Cosh z)  n even
 *                                    (E Cosh z - O Sinh z)  n odd
 *   I_{-(n+1/2)}(z) = the same with Sinh and Cosh exchanged. */
static Expr* sl_sph_bessel_i(int64_t n, bool negative_order, const Expr* z) {
    Expr* E = sl_hankel_sum(n, 0, z);
    Expr* O = sl_hankel_sum(n, 1, z);
    if (!E || !O) { expr_free(E); expr_free(O); return NULL; }
    bool sinh_first = ((n & 1) == 0) != negative_order;
    Expr* body = mk_fn2("Plus",
        mk_fn2("Times", E, mk_fn1(sinh_first ? "Sinh" : "Cosh", expr_copy((Expr*)z))),
        expr_new_function(expr_new_symbol("Times"), (Expr*[]){ mk_int(-1), O,
            mk_fn1(sinh_first ? "Cosh" : "Sinh", expr_copy((Expr*)z)) }, 3));
    return expr_new_function(expr_new_symbol("Times"), (Expr*[]){
        mk_fn2("Power", mk_int(2), mk_rat(1, 2)),
        mk_fn2("Power", expr_new_symbol("Pi"), mk_rat(-1, 2)),
        mk_fn2("Power", expr_copy((Expr*)z), mk_rat(-1, 2)),
        body }, 4);
}

/* A half-odd-integer order nu = k/2 (k odd), n >= 0 (DLMF 11.4.4, and 11.6.2,
 * whose asymptotic series terminates for this order):
 *   L_{-(n+1/2)}(z) = I_{n+1/2}(z)
 *   L_{n+1/2}(z)    = I_{-(n+1/2)}(z)
 *                     - (1/pi) Sum_{j=0}^n (-1)^j Gamma(j+1/2) (z/2)^(n-2j-1/2) / (n-j)! */
static Expr* sl_half_integer(const Expr* nu, const Expr* z) {
    if (nu->type != EXPR_FUNCTION || !is_rational_like(nu)) return NULL;
    const Expr* num = nu->data.function.args[0];
    const Expr* den = nu->data.function.args[1];
    if (num->type != EXPR_INTEGER || den->type != EXPR_INTEGER || den->data.integer != 2)
        return NULL;
    int64_t k = num->data.integer;
    int64_t n = k < 0 ? (-k - 1) / 2 : (k - 1) / 2;
    Expr* bes = sl_sph_bessel_i(n, k > 0, z);
    if (!bes) return NULL;
    if (k < 0) return sl_eval(bes);
    Expr** t = malloc(sizeof(Expr*) * (size_t)(n + 1));
    if (!t) { expr_free(bes); return NULL; }
    for (int64_t j = 0; j <= n; j++) {
        Expr* e = mk_rat(2 * (n - 2 * j) - 1, 2);
        t[j] = expr_new_function(expr_new_symbol("Times"), (Expr*[]){
            mk_int((j & 1) ? 1 : -1),        /* the leading minus folded in */
            mk_fn1("Gamma", mk_rat(2 * j + 1, 2)),
            mk_fn2("Power", mk_int(2), mk_fn2("Times", mk_int(-1), expr_copy(e))),
            mk_fn2("Power", expr_copy((Expr*)z), e),
            mk_fn2("Power", mk_fn1("Factorial", mk_int(n - j)), mk_int(-1)) }, 5);
    }
    Expr* poly = expr_new_function(expr_new_symbol("Plus"), t, (size_t)(n + 1));
    free(t);
    Expr* part = mk_fn2("Times", mk_fn2("Power", expr_new_symbol("Pi"), mk_int(-1)), poly);
    return sl_eval(mk_fn2("Plus", part, bes));
}

Expr* builtin_struvel(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* nu = res->data.function.args[0];
    Expr* z  = res->data.function.args[1];

    /* L_nu(0) = 0 for a real nu > -1 (exact). */
    if (sf2_exact_zero(z)) {
        double nuv;
        if (common_machine_real_value(nu, &nuv) && nuv > -1.0) return mk_int(0);
    }
    bool handled;
    Expr* v = sf2_numeric(nu, z, 0, sl_machine, sl_arb, &handled);
    if (handled) return v;
    return sl_half_integer(nu, z);
}

void struvel_init(void) {
    symtab_add_builtin("StruveL", builtin_struvel);
    symtab_get_def("StruveL")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    symtab_set_docstring("StruveL",
        "StruveL[nu, z]\n"
        "The modified Struve function L_nu(z), the particular solution of\n"
        "z^2 y'' + z y' - (z^2 + nu^2) y = 4 (z/2)^(nu+1)/(Sqrt[Pi] Gamma[nu+1/2]);\n"
        "L_nu(z) = (z/2)^(nu+1) Sum_k (z/2)^(2k)/(Gamma[k+3/2] Gamma[k+nu+3/2]).\n"
        "StruveL[nu, 0] = 0 for nu > -1. A half-odd-integer order reduces to\n"
        "Sinh, Cosh and powers. Machine, arbitrary-precision and complex arguments;\n"
        "fast paths over packed arrays and in Compile[]. D[StruveL[nu, z], z] =\n"
        "((z/2)^nu/(Sqrt[Pi] Gamma[nu+3/2]) + StruveL[nu-1, z] + StruveL[nu+1, z])/2.");
}
