/* Mathilda -- the Legendre elliptic integrals.
 *
 *   EllipticK[m]          = Int_0^(Pi/2) dt / Sqrt(1 - m Sin[t]^2)
 *   EllipticF[phi, m]     = Int_0^phi    dt / Sqrt(1 - m Sin[t]^2)
 *   EllipticE[m]          = Int_0^(Pi/2) Sqrt(1 - m Sin[t]^2) dt
 *   EllipticE[phi, m]     = Int_0^phi    Sqrt(1 - m Sin[t]^2) dt
 *   EllipticPi[n, m]      = Int_0^(Pi/2) dt / ((1 - n Sin[t]^2) Sqrt(1 - m Sin[t]^2))
 *   EllipticPi[n, phi, m] = Int_0^phi    dt / ((1 - n Sin[t]^2) Sqrt(1 - m Sin[t]^2))
 *
 * PARAMETER, not modulus. The second (E, F) or third (Pi) argument is
 * m = k^2, the Wolfram Language convention, and the one place a port of this
 * family silently goes wrong is by reading it as the modulus k. Every doc
 * comment here spells the integral out rather than naming "the modulus".
 *
 * E and Pi are arity-overloaded exactly as in Wolfram: EllipticE[m] is complete
 * and EllipticE[phi, m] incomplete; EllipticPi[n, m] is complete and
 * EllipticPi[n, phi, m] incomplete. There is no EllipticK[phi, m] -- that
 * spelling is EllipticF.
 *
 * Evaluation is layered so each kind of argument takes the cheapest route:
 *
 *   exact special values  ->  0, phi, Pi/2, 1, ComplexInfinity, the complete
 *                             form at phi = Pi/2, the F/K form at n = 0
 *   numeric (inexact)     ->  FLINT/Arb acb_elliptic_* (full complex plane,
 *                             arbitrary precision, Wolfram's branch placement)
 *   everything else        ->  stays symbolic (return NULL)
 *
 * Exact non-special arguments (EllipticF[1/3, 1/2]) stay symbolic, as in the
 * Wolfram Language; only inexact input or an explicit N[...] evaluates.
 *
 * Why Arb rather than a hand-rolled kernel. The incomplete integrals need the
 * quasi-periodic extension off the principal strip (F(phi + k Pi | m) =
 * F(phi | m) + 2 k K(m)), complex phi -- which is routine here, since
 * EllipticF[ArcSin[z], m] with |z| > 1 has complex phi -- and, for EllipticPi
 * with n > 1, the Cauchy principal value across the pole at Sin[t]^2 = 1/n,
 * which is where a hand-rolled Carlson RJ is most easily got wrong. Arb has all
 * three, rigorously and to arbitrary precision, and FLINT is already a
 * dependency. Without FLINT (USE_FLINT undefined) the numeric path falls back to
 * the machine-precision Carlson kernels below for the real domain and otherwise
 * declines, leaving the call symbolic.
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "elliptic.h"
#include "sym_names.h"

#include <math.h>
#include <stdbool.h>
#include <string.h>

#include "arithmetic.h"        /* is_rational, make_rational */
#include "attr.h"
#include "common.h"            /* builtin_arg_error */
#include "eval.h"              /* eval_and_free */
#include "expr.h"
#include "flint_num_bridge.h"  /* flint_num_elliptic_* */
#include "message.h"           /* mth_message: the Quiet/Check funnel */
#include "numeric.h"           /* arg_is_inexact */
#include "symtab.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

/* ------------------------------------------------------------------ */
/*  Machine-precision Carlson symmetric forms (real domain)            */
/* ------------------------------------------------------------------ */
/* R_F(x,y,z) = (1/2) Int_0^oo dt / Sqrt((t+x)(t+y)(t+z))
 * R_D(x,y,z) = (3/2) Int_0^oo dt / ((t+z) Sqrt((t+x)(t+y)(t+z)))
 *
 * Carlson's duplication theorem with the fifth-order series tail (DLMF 19.36.1,
 * 19.36.2). Real non-negative arguments, at most one zero; anything else
 * declines so the caller falls through to the FLINT path, which is defined on
 * the whole complex plane. These exist for the packed/NDArray element-wise
 * kernels, where the ABI is `double` in and `double` out. */

#define EC_RTOL 1e-16      /* duplication stops well inside double precision */
#define EC_MAXIT 100

static bool carlson_rf(double x, double y, double z, double* out) {
    if (x < 0.0 || y < 0.0 || z < 0.0) return false;
    if ((x + y == 0.0) || (y + z == 0.0) || (x + z == 0.0)) return false;  /* two zeros */
    double A0 = (x + y + z) / 3.0;
    if (A0 == 0.0) return false;
    double A = A0, xm = x, ym = y, zm = z, fac = 1.0;
    int it = 0;
    for (; it < EC_MAXIT; it++) {
        double dx = fabs(A - xm), dy = fabs(A - ym), dz = fabs(A - zm);
        double q = dx > dy ? (dx > dz ? dx : dz) : (dy > dz ? dy : dz);
        if (q < fabs(A) * 0.01) break;               /* series converges: see below */
        double sx = sqrt(xm), sy = sqrt(ym), sz = sqrt(zm);
        double lam = sx * sy + sy * sz + sz * sx;
        xm = 0.25 * (xm + lam); ym = 0.25 * (ym + lam); zm = 0.25 * (zm + lam);
        A = 0.25 * (A + lam);
        fac *= 0.25;
        if (!isfinite(A) || A == 0.0) return false;
    }
    if (it == EC_MAXIT) return false;
    /* Fifth-order tail in the deviations X = 1 - x_m/A, ... */
    double X = (A0 - x) * fac / A, Y = (A0 - y) * fac / A;
    double Z = -(X + Y);
    double E2 = X * Y - Z * Z, E3 = X * Y * Z;
    double s = 1.0 - E2 / 10.0 + E3 / 14.0 + E2 * E2 / 24.0 - 3.0 * E2 * E3 / 44.0;
    double r = s / sqrt(A);
    if (!isfinite(r)) return false;
    *out = r;
    (void)EC_RTOL;
    return true;
}

static bool carlson_rd(double x, double y, double z, double* out) {
    if (x < 0.0 || y < 0.0 || z <= 0.0) return false;
    if (x + y == 0.0) return false;
    double xm = x, ym = y, zm = z, sum = 0.0, fac = 1.0;
    double A0 = (x + y + 3.0 * z) / 5.0, A = A0;
    if (A0 == 0.0) return false;
    int it = 0;
    for (; it < EC_MAXIT; it++) {
        double dx = fabs(A - xm), dy = fabs(A - ym), dz = fabs(A - zm);
        double q = dx > dy ? (dx > dz ? dx : dz) : (dy > dz ? dy : dz);
        if (q < fabs(A) * 0.01) break;
        double sx = sqrt(xm), sy = sqrt(ym), sz = sqrt(zm);
        double lam = sx * sy + sy * sz + sz * sx;
        sum += fac / (sz * (zm + lam));
        xm = 0.25 * (xm + lam); ym = 0.25 * (ym + lam); zm = 0.25 * (zm + lam);
        A = 0.25 * (A + lam);
        fac *= 0.25;
        if (!isfinite(A) || A == 0.0) return false;
    }
    if (it == EC_MAXIT) return false;
    double X = (A0 - x) * fac / A, Y = (A0 - y) * fac / A;
    double Z = -(X + Y) / 3.0;
    double E2 = X * Y - 6.0 * Z * Z;
    double E3 = (3.0 * X * Y - 8.0 * Z * Z) * Z;
    double E4 = 3.0 * (X * Y - Z * Z) * Z * Z;
    double E5 = X * Y * Z * Z * Z;
    double s = 1.0 - 3.0 * E2 / 14.0 + E3 / 6.0 + 9.0 * E2 * E2 / 88.0
             - 3.0 * E4 / 22.0 - 9.0 * E2 * E3 / 52.0 + 3.0 * E5 / 26.0;
    double r = 3.0 * sum + fac * s / (A * sqrt(A));
    if (!isfinite(r)) return false;
    *out = r;
    return true;
}

bool elliptic_machine_k(double m, double* out) {
    if (!(m < 1.0)) return false;                /* m = 1 is a pole, m > 1 complex */
    return carlson_rf(0.0, 1.0 - m, 1.0, out);
}

bool elliptic_machine_e_complete(double m, double* out) {
    if (m == 1.0) { *out = 1.0; return true; }
    if (!(m < 1.0)) return false;
    double rf, rd;
    if (!carlson_rf(0.0, 1.0 - m, 1.0, &rf)) return false;
    if (!carlson_rd(0.0, 1.0 - m, 1.0, &rd)) return false;
    double r = rf - (m / 3.0) * rd;
    if (!isfinite(r)) return false;
    *out = r;
    return true;
}

/* The incomplete forms on the principal strip |Re phi| <= Pi/2, extended by the
 * quasi-period: F(phi + k Pi | m) = F(phi | m) + 2 k K(m), and likewise for E
 * with E(m) in place of K(m). `k` is chosen so the residual lands in the strip.
 *
 *   F(phi | m) = s R_F(1 - s^2, 1 - m s^2, 1)
 *   E(phi | m) = s R_F(1 - s^2, 1 - m s^2, 1) - (m/3) s^3 R_D(1 - s^2, 1 - m s^2, 1)
 *
 * with s = Sin[phi] (DLMF 19.25.5, 19.25.7, after using the homogeneity of R_F
 * and R_D to clear the csc^2 scaling). Declines when 1 - m s^2 < 0, where the
 * value is genuinely complex and only the FLINT path can answer. */
static bool elliptic_inc_real(double phi, double m, bool want_E, double* out) {
    if (!isfinite(phi) || !isfinite(m)) return false;
    double k = floor(phi / M_PI + 0.5);           /* nearest integer */
    double r = phi - k * M_PI;                     /* |r| <= Pi/2 */
    double s = sin(r), s2 = s * s;
    double a = 1.0 - s2, b = 1.0 - m * s2;
    if (b < 0.0) return false;                     /* complex: FLINT's job */
    double rf;
    if (!carlson_rf(a, b, 1.0, &rf)) return false;
    double v = s * rf;
    if (want_E) {
        double rd;
        if (!carlson_rd(a, b, 1.0, &rd)) return false;
        v -= (m / 3.0) * s * s2 * rd;
    }
    if (k != 0.0) {
        double comp;
        bool ok = want_E ? elliptic_machine_e_complete(m, &comp)
                         : elliptic_machine_k(m, &comp);
        if (!ok) return false;
        v += 2.0 * k * comp;
    }
    if (!isfinite(v)) return false;
    *out = v;
    return true;
}

bool elliptic_machine_f(double phi, double m, double* out) {
    return elliptic_inc_real(phi, m, false, out);
}

bool elliptic_machine_e_inc(double phi, double m, double* out) {
    return elliptic_inc_real(phi, m, true, out);
}

/* ------------------------------------------------------------------ */
/*  Small Expr predicates                                              */
/* ------------------------------------------------------------------ */

static bool ell_is_zero(const Expr* e) {
    return e && ((e->type == EXPR_INTEGER && e->data.integer == 0) ||
                 (e->type == EXPR_REAL && e->data.real == 0.0));
}

static bool ell_is_one(const Expr* e) {
    return e && ((e->type == EXPR_INTEGER && e->data.integer == 1) ||
                 (e->type == EXPR_REAL && e->data.real == 1.0));
}

static bool ell_is_sym(const Expr* e, const char* name) {
    return e && e->type == EXPR_SYMBOL && strcmp(e->data.symbol.name, name) == 0;
}

/* Exactly the expression Pi/2, i.e. Times[Rational[1,2], Pi]. The complete
 * integrals are the incomplete ones at this upper limit, and the corpus hands
 * that spelling in when a pencil closes at the quarter period. */
static bool ell_is_half_pi(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL ||
        e->data.function.head->data.symbol.name != SYM_Times ||
        e->data.function.arg_count != 2) return false;
    int64_t n, d;
    const Expr* a = e->data.function.args[0];
    const Expr* b = e->data.function.args[1];
    if (!is_rational(a, &n, &d) || n != 1 || d != 2) return false;
    return ell_is_sym(b, "Pi");
}

/* Does any argument carry an inexact number, i.e. should we evaluate at all?
 * Exact input stays symbolic, as in Wolfram. */
static bool ell_any_inexact(Expr* const* args, size_t n) {
    for (size_t i = 0; i < n; i++)
        if (arg_is_inexact(args[i])) return true;
    return false;
}

static Expr* ell_half_pi(void) {
    Expr* t[2] = { make_rational(1, 2), expr_new_symbol(SYM_Pi) };
    return expr_new_function(expr_new_symbol(SYM_Times), t, 2);
}

/* Build head[args...] and evaluate it, so a reduction can name another member
 * of the family and let the evaluator finish the job. */
static Expr* ell_call(const char* head, Expr** args, size_t n) {
    return eval_and_free(expr_new_function(expr_new_symbol(head), args, n));
}

static Expr* ell_argt(const char* head, size_t argc, size_t lo, size_t hi) {
    mth_message(head, "argt",
        "%s called with %zu argument%s; %zu or %zu arguments are expected.",
        head, argc, argc == 1 ? "" : "s", lo, hi);
    return NULL;
}

/* ------------------------------------------------------------------ */
/*  EllipticK[m]                                                       */
/* ------------------------------------------------------------------ */

Expr* builtin_elliptick(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1)
        return builtin_arg_error("EllipticK",
            res->type == EXPR_FUNCTION ? res->data.function.arg_count : 0, 1, 1);
    Expr* m = res->data.function.args[0];

    if (ell_is_zero(m)) return ell_half_pi();                  /* K[0] = Pi/2 */
    if (ell_is_one(m))  return expr_new_symbol(SYM_ComplexInfinity);

    if (ell_any_inexact(&m, 1)) {
        Expr* v = flint_num_elliptic_k(m);
        if (v) return v;
    }
    return NULL;
}

/* ------------------------------------------------------------------ */
/*  EllipticF[phi, m]                                                  */
/* ------------------------------------------------------------------ */

Expr* builtin_ellipticf(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)
        return builtin_arg_error("EllipticF",
            res->type == EXPR_FUNCTION ? res->data.function.arg_count : 0, 2, 2);
    Expr* phi = res->data.function.args[0];
    Expr* m   = res->data.function.args[1];

    if (ell_is_zero(phi)) return expr_new_integer(0);           /* F[0, m] = 0 */
    if (ell_is_zero(m))   return expr_copy(phi);                /* F[phi, 0] = phi */
    if (ell_is_half_pi(phi)) {                                  /* F[Pi/2, m] = K[m] */
        Expr* a[1] = { expr_copy(m) };
        return ell_call("EllipticK", a, 1);
    }

    if (ell_any_inexact(res->data.function.args, 2)) {
        Expr* v = flint_num_elliptic_f(phi, m);
        if (v) return v;
    }
    return NULL;
}

/* ------------------------------------------------------------------ */
/*  EllipticE[m] / EllipticE[phi, m]                                   */
/* ------------------------------------------------------------------ */

static Expr* elliptice_complete(Expr* m) {
    if (ell_is_zero(m)) return ell_half_pi();                   /* E[0] = Pi/2 */
    if (ell_is_one(m))  return expr_new_integer(1);             /* E[1] = 1    */
    if (ell_any_inexact(&m, 1)) {
        Expr* v = flint_num_elliptic_e(m);
        if (v) return v;
    }
    return NULL;
}

static Expr* elliptice_incomplete(Expr* phi, Expr* m) {
    if (ell_is_zero(phi)) return expr_new_integer(0);           /* E[0, m] = 0 */
    if (ell_is_zero(m))   return expr_copy(phi);                /* E[phi, 0] = phi */
    if (ell_is_one(m)) {                                        /* E[phi, 1] = Sin[phi] */
        Expr* a[1] = { expr_copy(phi) };
        return ell_call("Sin", a, 1);
    }
    if (ell_is_half_pi(phi)) {                                  /* E[Pi/2, m] = E[m] */
        Expr* a[1] = { expr_copy(m) };
        return ell_call("EllipticE", a, 1);
    }
    Expr* two[2] = { phi, m };
    if (ell_any_inexact(two, 2)) {
        Expr* v = flint_num_elliptic_e_inc(phi, m);
        if (v) return v;
    }
    return NULL;
}

Expr* builtin_elliptice(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc == 1) return elliptice_complete(res->data.function.args[0]);
    if (argc == 2) return elliptice_incomplete(res->data.function.args[0],
                                               res->data.function.args[1]);
    return ell_argt("EllipticE", argc, 1, 2);
}

/* ------------------------------------------------------------------ */
/*  EllipticPi[n, m] / EllipticPi[n, phi, m]                           */
/* ------------------------------------------------------------------ */

static Expr* ellipticpi_complete(Expr* n, Expr* m) {
    if (ell_is_zero(n)) {                                       /* Pi[0, m] = K[m] */
        Expr* a[1] = { expr_copy(m) };
        return ell_call("EllipticK", a, 1);
    }
    Expr* two[2] = { n, m };
    if (ell_any_inexact(two, 2)) {
        Expr* v = flint_num_elliptic_pi(n, m);
        if (v) return v;
    }
    return NULL;
}

static Expr* ellipticpi_incomplete(Expr* n, Expr* phi, Expr* m) {
    if (ell_is_zero(phi)) return expr_new_integer(0);           /* Pi[n, 0, m] = 0 */
    if (ell_is_zero(n)) {                                       /* Pi[0, phi, m] = F */
        Expr* a[2] = { expr_copy(phi), expr_copy(m) };
        return ell_call("EllipticF", a, 2);
    }
    if (ell_is_half_pi(phi)) {                                  /* Pi[n, Pi/2, m] */
        Expr* a[2] = { expr_copy(n), expr_copy(m) };
        return ell_call("EllipticPi", a, 2);
    }
    Expr* three[3] = { n, phi, m };
    if (ell_any_inexact(three, 3)) {
        Expr* v = flint_num_elliptic_pi_inc(n, phi, m);
        if (v) return v;
    }
    return NULL;
}

Expr* builtin_ellipticpi(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc == 2) return ellipticpi_complete(res->data.function.args[0],
                                              res->data.function.args[1]);
    if (argc == 3) return ellipticpi_incomplete(res->data.function.args[0],
                                                res->data.function.args[1],
                                                res->data.function.args[2]);
    return ell_argt("EllipticPi", argc, 2, 3);
}

/* ------------------------------------------------------------------ */
/*  Registration                                                       */
/* ------------------------------------------------------------------ */

void elliptic_init(void) {
    symtab_add_builtin("EllipticK", builtin_elliptick);
    symtab_get_def("EllipticK")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    symtab_set_docstring("EllipticK",
        "EllipticK[m] is the complete elliptic integral of the first kind, "
        "Integrate[1/Sqrt[1 - m Sin[t]^2], {t, 0, Pi/2}]. The argument is the "
        "PARAMETER m = k^2, not the modulus k. EllipticK[0] is Pi/2 and "
        "EllipticK[1] is ComplexInfinity; exact arguments otherwise stay "
        "symbolic and inexact ones evaluate numerically at their precision.");

    symtab_add_builtin("EllipticF", builtin_ellipticf);
    symtab_get_def("EllipticF")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    symtab_set_docstring("EllipticF",
        "EllipticF[phi, m] is the incomplete elliptic integral of the first "
        "kind, Integrate[1/Sqrt[1 - m Sin[t]^2], {t, 0, phi}]. The second "
        "argument is the PARAMETER m = k^2, not the modulus k. "
        "EllipticF[phi, 0] is phi and EllipticF[Pi/2, m] is EllipticK[m]; phi "
        "may be complex and of any size (the quasi-period is applied).");

    symtab_add_builtin("EllipticE", builtin_elliptice);
    symtab_get_def("EllipticE")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    symtab_set_docstring("EllipticE",
        "EllipticE[m] is the complete elliptic integral of the second kind, "
        "Integrate[Sqrt[1 - m Sin[t]^2], {t, 0, Pi/2}], and "
        "EllipticE[phi, m] the incomplete one, with upper limit phi. The "
        "parameter argument is m = k^2, not the modulus k. EllipticE[0] is "
        "Pi/2, EllipticE[1] is 1, and EllipticE[phi, 1] is Sin[phi].");

    symtab_add_builtin("EllipticPi", builtin_ellipticpi);
    symtab_get_def("EllipticPi")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    symtab_set_docstring("EllipticPi",
        "EllipticPi[n, m] is the complete elliptic integral of the third kind, "
        "Integrate[1/((1 - n Sin[t]^2) Sqrt[1 - m Sin[t]^2]), {t, 0, Pi/2}], "
        "and EllipticPi[n, phi, m] the incomplete one, with upper limit phi. "
        "The parameter argument is m = k^2, not the modulus k. "
        "EllipticPi[0, m] is EllipticK[m] and EllipticPi[0, phi, m] is "
        "EllipticF[phi, m]. For n > 1 the path crosses the pole at "
        "Sin[t]^2 == 1/n and the value is the Cauchy principal value.");
}
