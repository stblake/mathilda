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
 * past the pole at Sin[t]^2 = 1/n, a value that is genuinely COMPLEX rather than
 * a real Cauchy principal value: Pi[3/2 | 1/2] is
 * -0.456720313453 - 2.72069904635 I (mpmath, and Arb agrees). An earlier version
 * of this comment claimed a real principal value there, which is why the machine
 * kernels were thought impossible to write; what they actually have to do is
 * DECLINE when 1 - n (complete) or 1 - n Sin[phi]^2 (incomplete) is <= 0, the
 * same contract K/E/F already follow outside their real domains. Arb has all of
 * it, rigorously and to arbitrary precision, and FLINT is already a dependency.
 * Without FLINT (USE_FLINT undefined) the numeric path falls back to the
 * machine-precision Carlson kernels below for the real domain and otherwise
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
#include "ndarray.h"           /* is_ndarray, ndarray_delist_and_reeval */
#include "numeric.h"           /* arg_is_inexact */
#include "symtab.h"

#ifdef USE_MPFR
#include <mpfr.h>
#endif

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

/* Where the duplication stops. The tail above is fifth order in the deviations
 * X, Y, Z, so its truncation error goes as (q/A)^6 -- and that is what fixes
 * these constants, which are Numerical Recipes' (rf 0.0025, rd 0.0015) and
 * deliver ~2e-16.
 *
 * They were 0.01 for both, with a `#define EC_RTOL 1e-16` alongside that was
 * never read (`(void)EC_RTOL`) and a comment claiming "well inside double
 * precision". (q/A)^6 at 0.01 is 1e-12, not 1e-16, and the delivered accuracy
 * was measured at 61 ulp (K), 106 ulp (E) and 70 ulp (F) against mpmath -- two
 * digits short of the scalar path's 1.5 ulp on the same inputs, which is a
 * SURFACE SKEW: the same expression answered differently packed and unpacked.
 * Tightening to 0.0025/0.0015 costs exactly one more duplication (the deviations
 * fall ~4x per step, and log4(0.01/0.0025) = 1). */
#define EC_ERRTOL_RF 0.0025
#define EC_ERRTOL_RD 0.0015
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
        if (q < fabs(A) * EC_ERRTOL_RF) break;        /* tail is accurate: see EC_ERRTOL_* */
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
    return true;
}

/* R_F and R_D together, in ONE duplication loop.
 *
 * Every caller that wants R_D wants R_F at the same arguments -- E(m) is
 * R_F - (m/3) R_D and E(phi|m) likewise -- and the two recurrences walk the
 * IDENTICAL x_m, y_m, z_m sequence: same lambda, same quarter-update. Only the
 * weighted mean differs ((x+y+z)/3 against (x+y+3z)/5) and R_D's running sum,
 * whose `sz` the shared step has already computed. Running them separately pays
 * for the three square roots per step twice; fusing halves the sqrt count of
 * every E call, which is what buys back the extra duplication that
 * EC_ERRTOL_RF/RD cost (E over 10^6 elements: 174 ns/elt before the tolerance
 * fix, 210 after it, and this brings it back under).
 *
 * z > 0 is required (R_D's integrand has a t+z in the denominator), which costs
 * nothing: every call site here passes z = 1. */
static bool carlson_rf_rd(double x, double y, double z, double* out_rf, double* out_rd) {
    if (x < 0.0 || y < 0.0 || z <= 0.0) return false;
    if (x + y == 0.0) return false;
    double Af0 = (x + y + z) / 3.0;
    double Ad0 = (x + y + 3.0 * z) / 5.0;
    if (Af0 == 0.0 || Ad0 == 0.0) return false;
    double Af = Af0, Ad = Ad0;
    double xm = x, ym = y, zm = z, sum = 0.0, fac = 1.0;
    int it = 0;
    for (; it < EC_MAXIT; it++) {
        double dx = fabs(Ad - xm), dy = fabs(Ad - ym), dz = fabs(Ad - zm);
        double qd = dx > dy ? (dx > dz ? dx : dz) : (dy > dz ? dy : dz);
        dx = fabs(Af - xm); dy = fabs(Af - ym); dz = fabs(Af - zm);
        double qf = dx > dy ? (dx > dz ? dx : dz) : (dy > dz ? dy : dz);
        /* Both tails must be in range before either may stop; R_D's is the
         * tighter demand, so this normally costs R_F nothing but accuracy. */
        if (qd < fabs(Ad) * EC_ERRTOL_RD && qf < fabs(Af) * EC_ERRTOL_RF) break;
        double sx = sqrt(xm), sy = sqrt(ym), sz = sqrt(zm);
        double lam = sx * sy + sy * sz + sz * sx;
        sum += fac / (sz * (zm + lam));
        xm = 0.25 * (xm + lam); ym = 0.25 * (ym + lam); zm = 0.25 * (zm + lam);
        Af = 0.25 * (Af + lam);
        Ad = 0.25 * (Ad + lam);
        fac *= 0.25;
        if (!isfinite(Af) || Af == 0.0 || !isfinite(Ad) || Ad == 0.0) return false;
    }
    if (it == EC_MAXIT) return false;

    /* R_F's fifth-order tail (DLMF 19.36.1). */
    double Xf = (Af0 - x) * fac / Af, Yf = (Af0 - y) * fac / Af;
    double Zf = -(Xf + Yf);
    double F2 = Xf * Yf - Zf * Zf, F3 = Xf * Yf * Zf;
    double sf = 1.0 - F2 / 10.0 + F3 / 14.0 + F2 * F2 / 24.0 - 3.0 * F2 * F3 / 44.0;
    double rf = sf / sqrt(Af);

    /* R_D's (DLMF 19.36.2). */
    double X = (Ad0 - x) * fac / Ad, Y = (Ad0 - y) * fac / Ad;
    double Z = -(X + Y) / 3.0;
    double E2 = X * Y - 6.0 * Z * Z;
    double E3 = (3.0 * X * Y - 8.0 * Z * Z) * Z;
    double E4 = 3.0 * (X * Y - Z * Z) * Z * Z;
    double E5 = X * Y * Z * Z * Z;
    double sd = 1.0 - 3.0 * E2 / 14.0 + E3 / 6.0 + 9.0 * E2 * E2 / 88.0
              - 3.0 * E4 / 22.0 - 9.0 * E2 * E3 / 52.0 + 3.0 * E5 / 26.0;
    double rd = 3.0 * sum + fac * sd / (Ad * sqrt(Ad));

    if (!isfinite(rf) || !isfinite(rd)) return false;
    *out_rf = rf;
    *out_rd = rd;
    return true;
}

/* R_C(x,y) = R_F(x,y,y), the degenerate symmetric form. It exists here only
 * because R_J needs it once per duplication step. Duplication with the
 * fourth-order tail (DLMF 19.36.3). y < 0 is reached through Carlson's
 * transformation, which is what makes the running sum in R_J legitimate. */
#define EC_ERRTOL_RC 0.0012

static bool carlson_rc(double x, double y, double* out) {
    if (x < 0.0 || y == 0.0) return false;
    double xt, yt, w;
    if (y > 0.0) {
        xt = x; yt = y; w = 1.0;
    } else {
        xt = x - y; yt = -y;
        if (xt <= 0.0) return false;
        w = sqrt(x) / sqrt(xt);
    }
    double ave, s;
    int it = 0;
    do {
        if (it++ == EC_MAXIT) return false;
        double lam = 2.0 * sqrt(xt) * sqrt(yt) + yt;
        xt = 0.25 * (xt + lam);
        yt = 0.25 * (yt + lam);
        ave = (xt + yt + yt) / 3.0;
        if (!isfinite(ave) || ave == 0.0) return false;
        s = (yt - ave) / ave;
    } while (fabs(s) > EC_ERRTOL_RC);
    double r = w * (1.0 + s * s * (0.3 + s * (1.0 / 7.0 + s * (0.375 + s * (9.0 / 22.0)))))
             / sqrt(ave);
    if (!isfinite(r)) return false;
    *out = r;
    return true;
}

/* R_J(x,y,z,p) = (3/2) Int_0^oo dt / ((t+p) Sqrt((t+x)(t+y)(t+z))).
 *
 * This is the kernel EllipticPi needs and did not have, which left the third
 * kind with NO machine path at all: ~40 us per element against SciPy's 316 ns
 * for the same Carlson composition.
 *
 * p <= 0 DECLINES, and that is a real mathematical boundary rather than
 * caution. p is 1 - n for the complete form and 1 - n Sin[phi]^2 for the
 * incomplete one, so p <= 0 is exactly where Mathematica's value stops being
 * real: Pi[3/2 | 1/2] is -0.456720313453 - 2.72069904635 I, not a real Cauchy
 * principal value. A `double` out parameter cannot carry that, so the honest
 * machine answer is to decline and let Arb place the branch.
 *
 * Duplication with the fifth-order tail, NR's arrangement (Carlson 1979). */
#define EC_ERRTOL_RJ 0.0015

static bool carlson_rj(double x, double y, double z, double p, double* out) {
    if (x < 0.0 || y < 0.0 || z < 0.0) return false;
    if (p <= 0.0) return false;                       /* genuinely complex: Arb's job */
    if (x + y == 0.0 || y + z == 0.0 || x + z == 0.0) return false;  /* two zeros */

    double xt = x, yt = y, zt = z, pt = p;
    double sum = 0.0, fac = 1.0, ave = 0.0;
    double dx = 0.0, dy = 0.0, dz = 0.0, dp = 0.0;
    int it = 0;
    do {
        if (it++ == EC_MAXIT) return false;
        double sx = sqrt(xt), sy = sqrt(yt), sz = sqrt(zt);
        double lam   = sx * (sy + sz) + sy * sz;
        double alpha = pt * (sx + sy + sz) + sx * sy * sz;
        alpha *= alpha;
        double beta  = pt * (pt + lam) * (pt + lam);
        double rc;
        if (!carlson_rc(alpha, beta, &rc)) return false;
        sum += fac * rc;
        fac *= 0.25;
        xt = 0.25 * (xt + lam); yt = 0.25 * (yt + lam);
        zt = 0.25 * (zt + lam); pt = 0.25 * (pt + lam);
        ave = 0.2 * (xt + yt + zt + pt + pt);
        if (!isfinite(ave) || ave == 0.0) return false;
        dx = (ave - xt) / ave; dy = (ave - yt) / ave;
        dz = (ave - zt) / ave; dp = (ave - pt) / ave;
    } while (fabs(dx) > EC_ERRTOL_RJ || fabs(dy) > EC_ERRTOL_RJ ||
             fabs(dz) > EC_ERRTOL_RJ || fabs(dp) > EC_ERRTOL_RJ);

    const double C1 = 3.0 / 14.0, C2 = 1.0 / 3.0, C3 = 3.0 / 22.0, C4 = 3.0 / 26.0;
    const double C5 = 0.75 * C3, C6 = 1.5 * C4, C7 = 0.5 * C2, C8 = C3 + C3;
    double ea = dx * (dy + dz) + dy * dz;
    double eb = dx * dy * dz;
    double ec = dp * dp;
    double ed = ea - 3.0 * ec;
    double ee = eb + 2.0 * dp * (ea - ec);
    double r = 3.0 * sum
             + fac * (1.0 + ed * (-C1 + C5 * ed - C6 * ee)
                          + eb * (C7 + dp * (-C8 + dp * C4))
                          + dp * ea * (C2 - dp * C3) - C2 * dp * ec)
               / (ave * sqrt(ave));
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
    if (!carlson_rf_rd(0.0, 1.0 - m, 1.0, &rf, &rd)) return false;
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
    /* a is Cos[r]^2, and it must be computed AS Cos[r]^2. Spelling it 1 - s2
     * loses everything near r = Pi/2: the subtraction carries an absolute error
     * of ~1e-16 while the true Cos[r]^2 is ~1e-16 itself, so R_F's first
     * argument arrives with ~100% relative error. Measured against the 30-digit
     * path at the same double, worst case over a Pi/2 approach ladder with
     * m = 0.99: EllipticF was wrong by 2.7e-08 relative (1.2e8 ulp, an
     * EIGHT-digit loss) at r = Pi/2 - 1e-8, and EllipticE by 9.8e-10 -- against
     * 6.7e-16 at a generic amplitude. cos(r) is correctly rounded, so squaring
     * it is accurate to 1 ulp wherever the answer is.
     *
     * fmax guards the one case squaring cannot: at r for which cos(r) underflows
     * to exactly 0, a must be +0 and not a rounding-negative value, or carlson_rf
     * declines at exactly the amplitude where F is simply K. */
    double c = cos(r);
    double a = c * c;
    if (a < 0.0) a = 0.0;
    double b = 1.0 - m * s2;
    if (b < 0.0) return false;                     /* complex: FLINT's job */
    double rf, rd;
    if (want_E) {
        if (!carlson_rf_rd(a, b, 1.0, &rf, &rd)) return false;
    } else {
        if (!carlson_rf(a, b, 1.0, &rf)) return false;
    }
    double v = s * rf;
    if (want_E) v -= (m / 3.0) * s * s2 * rd;
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

/* Pi(n | m) = R_F(0, 1-m, 1) + (n/3) R_J(0, 1-m, 1, 1-n).
 * Checked against mpmath at <= 2 ulp over n, m in (0,1) and <= 3.3 ulp with
 * n in (-6, 1). */
bool elliptic_machine_pi(double n, double m, double* out) {
    if (!isfinite(n) || !isfinite(m)) return false;
    if (!(m < 1.0)) return false;           /* m >= 1: complex (or a pole) */
    if (!(n < 1.0)) return false;           /* n >= 1: complex, see carlson_rj */
    double rf, rj;
    if (!carlson_rf(0.0, 1.0 - m, 1.0, &rf)) return false;
    if (!carlson_rj(0.0, 1.0 - m, 1.0, 1.0 - n, &rj)) return false;
    double r = rf + (n / 3.0) * rj;
    if (!isfinite(r)) return false;
    *out = r;
    return true;
}

/* Pi(n; phi | m) = s R_F(c^2, 1 - m s^2, 1) + (n/3) s^3 R_J(c^2, 1 - m s^2, 1, 1 - n s^2),
 * s = Sin[phi], c = Cos[phi]; checked against mpmath at <= 3.8 ulp over 400
 * random (n, phi, m) with n in (-6, 0.98), m in (-3, 0.98), phi in (0.01, 1.5).
 *
 * Extended off the principal strip by the quasi-period
 * Pi(n; phi + k Pi | m) = Pi(n; phi | m) + 2 k Pi(n | m), verified exactly
 * against mpmath at k = -1, 1, 2. Note the complete form appears in that
 * extension, so a shift declines whenever the complete form does.
 *
 * c^2 is Cos[phi]^2 computed as such, never 1 - Sin[phi]^2 -- see the comment
 * in elliptic_inc_real for the eight digits that spelling costs near Pi/2. */
bool elliptic_machine_pi_inc(double n, double phi, double m, double* out) {
    if (!isfinite(n) || !isfinite(phi) || !isfinite(m)) return false;
    double k = floor(phi / M_PI + 0.5);
    double r = phi - k * M_PI;                     /* |r| <= Pi/2 */
    double sn = sin(r), s2 = sn * sn;
    double c = cos(r);
    double a = c * c;
    if (a < 0.0) a = 0.0;
    double b = 1.0 - m * s2;
    double pp = 1.0 - n * s2;
    if (b < 0.0 || pp <= 0.0) return false;        /* complex: Arb's job */
    double rf, rj;
    if (!carlson_rf(a, b, 1.0, &rf)) return false;
    if (!carlson_rj(a, b, 1.0, pp, &rj)) return false;
    double v = sn * rf + (n / 3.0) * sn * s2 * rj;
    if (k != 0.0) {
        double comp;
        if (!elliptic_machine_pi(n, m, &comp)) return false;
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

/* Is `e` the number 0 (resp. 1), in any inexact spelling as well as the exact
 * one? EXPR_MPFR has to be in here: without it EllipticK[SetPrecision[1, 30]]
 * missed the pole, handed Arb an exact 1, got a non-finite ball back and
 * declined -- so the answer was an unevaluated EllipticK[1.0] where
 * EllipticK[1.] gives ComplexInfinity. The two spellings of the same number
 * must not disagree about whether they are at a pole. */
static bool ell_is_zero(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) return e->data.integer == 0;
    if (e->type == EXPR_REAL)    return e->data.real == 0.0;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR)    return mpfr_cmp_ui(e->data.mpfr, 0) == 0;
#endif
    return false;
}

static bool ell_is_one(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) return e->data.integer == 1;
    if (e->type == EXPR_REAL)    return e->data.real == 1.0;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR)    return mpfr_cmp_ui(e->data.mpfr, 1) == 0;
#endif
    return false;
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

/* True only when phi PROVABLY lies in the principal strip |phi| <= Pi/2.
 *
 * E(phi | 1) = Int_0^phi |Cos[t]| dt, which is Sin[phi] on that strip and
 * Sin[phi - k Pi] + 2k outside it. The rule used to be applied unconditionally,
 * so EllipticE[2, 1] answered Sin[2] = 0.909297 where the true value is
 * 2 - Sin[2] = 1.090703 -- a wrong answer, and one that jumped 0.18 against its
 * own neighbour (N[EllipticE[2, 1 - 10^-18], 20] is 1.0907025731743183...).
 *
 * Decided by asking the evaluator rather than by matching a pattern, so an
 * exact Rational, a Real, an MPFR and a symbol carrying assumptions all get an
 * answer; a symbolic phi is undecidable and answers false, which leaves
 * EllipticE[phi, 1] symbolic instead of wrong. Costs one evaluation, on a branch
 * that only fires when m is exactly 1. */
static bool ell_in_principal_strip(const Expr* phi) {
    Expr* ab[1] = { expr_copy((Expr*)phi) };
    Expr* le[2] = { ell_call(SYM_Abs, ab, 1), ell_half_pi() };
    Expr* q = ell_call(SYM_LessEqual, le, 2);
    bool in = q && q->type == EXPR_SYMBOL && q->data.symbol.name == SYM_True;
    expr_free(q);
    return in;
}

static Expr* ell_argt(const char* head, size_t argc, size_t lo, size_t hi) {
    mth_message(head, "argt",
        "%s called with %zu argument%s; %zu or %zu arguments are expected.",
        head, argc, argc == 1 ? "" : "s", lo, hi);
    return NULL;
}

/* A MACHINE-precision double for `e`, or false.
 *
 * EXPR_REAL and the exact numbers convert; EXPR_MPFR deliberately does NOT. A
 * 53-bit MPFR still carries a DECLARED precision, and SetPrecision[m, 16] must
 * not be silently answered by a double kernel; anything above 53 bits obviously
 * must not either. Complex and symbolic decline. */
static bool ell_double(const Expr* e, double* out) {
    if (!e) return false;
    int64_t n, d;
    switch (e->type) {
        case EXPR_REAL:    *out = e->data.real;        return true;
        case EXPR_INTEGER: *out = (double)e->data.integer; return true;
        case EXPR_BIGINT:  *out = mpz_get_d(e->data.bigint);
                           return isfinite(*out);
        default:
            if (is_rational(e, &n, &d) && d != 0) {
                *out = (double)n / (double)d;
                return isfinite(*out);
            }
            return false;
    }
}

/* Does this argument list justify a MACHINE answer, and what are the doubles?
 *
 * At least one argument must be an EXPR_REAL: with none, the call is either
 * wholly exact (and must stay symbolic) or carries an MPFR (and wants arbitrary
 * precision). Every argument must then convert, exact ones included -- the
 * result is machine precision either way, which is the same rule
 * pick_out_bits applies on the Arb side.
 *
 * This is what makes the elliptic scalars MACHINE NUMBERS. Every inexact scalar
 * used to go through acb + MPFR and come back as an EXPR_MPFR at 53 bits, so
 * Precision[EllipticK[0.5]] was 15.9546 where Zeta and Gamma give
 * MachinePrecision, MachineNumberQ was False, and -- the part that actually
 * costs something -- a 53-bit MPFR leaf does not pack, so every
 * elliptic-produced list was off the buffer for all its consumers. The house
 * convention is zeta.c:603-607: Real in, Real out. */
static bool ell_machine_args(Expr* const* a, size_t n, double* out) {
    bool any_real = false;
    for (size_t i = 0; i < n; i++) {
        if (a[i] && a[i]->type == EXPR_REAL) any_real = true;
        if (!ell_double(a[i], &out[i])) return false;
    }
    return any_real;
}

/* Small builders, so the closed forms below read like the identities they are.
 * Each goes through the evaluator (ell_call), which is this file's idiom and
 * what canonicalises the result. */
static Expr* ell_pow(Expr* b, Expr* e)   { Expr* a[2] = { b, e }; return ell_call(SYM_Power, a, 2); }
static Expr* ell_mul(Expr* x, Expr* y)   { Expr* a[2] = { x, y }; return ell_call(SYM_Times, a, 2); }
static Expr* ell_add(Expr* x, Expr* y)   { Expr* a[2] = { x, y }; return ell_call(SYM_Plus,  a, 2); }
static Expr* ell_sqrt(Expr* x)           { Expr* a[1] = { x };    return ell_call(SYM_Sqrt,  a, 1); }
static Expr* ell_gamma(Expr* x)          { Expr* a[1] = { x };    return ell_call(SYM_Gamma, a, 1); }
static Expr* ell_pi_sym(void)            { return expr_new_symbol(SYM_Pi); }

/* The two singular values of K that have a closed form, both checked against
 * mpmath to 30 digits:
 *   K(-1)  = Gamma[1/4]^2 / (4 Sqrt[2 Pi])          = 1.31102877714605990523...
 *   K(1/2) = 8 Pi^(3/2) / Gamma[-1/4]^2             = 1.85407467730137191843...
 * These are the lemniscatic cases; Mathematica gives both, and E has no
 * corresponding closed form at either point (it leaves them alone, so we do). */
static Expr* ell_k_closed_minus_one(void) {
    Expr* num = ell_pow(ell_gamma(make_rational(1, 4)), expr_new_integer(2));
    Expr* den = ell_mul(expr_new_integer(4), ell_sqrt(ell_mul(expr_new_integer(2), ell_pi_sym())));
    return ell_mul(num, ell_pow(den, expr_new_integer(-1)));
}

static Expr* ell_k_closed_half(void) {
    Expr* num = ell_mul(expr_new_integer(8), ell_pow(ell_pi_sym(), make_rational(3, 2)));
    Expr* den = ell_pow(ell_gamma(make_rational(-1, 4)), expr_new_integer(2));
    return ell_mul(num, ell_pow(den, expr_new_integer(-1)));
}

/* Pi(n | 0) = Pi / (2 Sqrt[1 - n]) -- the m = 0 integrand is 1/(1 - n Sin[t]^2),
 * whose quarter-period integral is elementary. */
static Expr* ell_pi_closed_m_zero(const Expr* n) {
    Expr* root = ell_sqrt(ell_add(expr_new_integer(1),
                                  ell_mul(expr_new_integer(-1), expr_copy((Expr*)n))));
    return ell_mul(ell_pi_sym(), ell_pow(ell_mul(expr_new_integer(2), root),
                                         expr_new_integer(-1)));
}

/* Pi(n; phi | 0) = ArcTanh[Sqrt[n-1] Tan[phi]] / Sqrt[n-1], the indefinite form
 * of the same elementary integrand. Verified numerically in-engine. */
static Expr* ell_pi_closed_m_zero_inc(const Expr* n, const Expr* phi) {
    Expr* root = ell_sqrt(ell_add(expr_new_integer(-1), expr_copy((Expr*)n)));
    Expr* tn[1] = { expr_copy((Expr*)phi) };
    Expr* at[1] = { ell_mul(expr_copy(root), ell_call(SYM_Tan, tn, 1)) };
    Expr* num = ell_call(SYM_ArcTanh, at, 1);
    return ell_mul(num, ell_pow(root, expr_new_integer(-1)));
}

/* Infinity in any direction, i.e. DirectedInfinity[...] or the Infinity /
 * ComplexInfinity symbols. The elliptic limits at infinity do not care which
 * direction, only that the argument is unbounded. */
static bool ell_is_infinite(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL)
        return e->data.symbol.name == SYM_Infinity ||
               e->data.symbol.name == SYM_ComplexInfinity;
    if (e->type == EXPR_FUNCTION && e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_DirectedInfinity) return true;
    return false;
}

/* f[-phi, ...] -> -f[phi, ...]. All three incomplete integrals are ODD in the
 * amplitude, because each integrand is EVEN in t. Uses the same superficial
 * negativity test the trig heads use (odd_fold, src/trig.c), so -x, -2x and
 * -x/3 fold while -x - y does not -- and the folded argument's leading
 * coefficient is then positive, so this cannot recur. `phi_at` is the amplitude
 * slot: 0 for EllipticF/EllipticE, 1 for EllipticPi. */
static Expr* ell_odd_in_phi(const char* head, Expr* const* args, size_t n, size_t phi_at) {
    if (!expr_is_superficially_negative(args[phi_at])) return NULL;
    Expr* a[3];
    for (size_t i = 0; i < n; i++) a[i] = expr_copy(args[i]);
    Expr* neg[2] = { expr_new_integer(-1), a[phi_at] };
    a[phi_at] = ell_call(SYM_Times, neg, 2);
    return ell_mul(expr_new_integer(-1), ell_call(head, a, n));
}

/* ------------------------------------------------------------------ */
/*  EllipticK[m]                                                       */
/* ------------------------------------------------------------------ */

Expr* builtin_elliptick(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1)
        return builtin_arg_error("EllipticK",
            res->type == EXPR_FUNCTION ? res->data.function.arg_count : 0, 1, 1);
    Expr* m = res->data.function.args[0];

    /* The pole is checked BEFORE the numeric path: Arb answers a non-finite ball
     * there and the bridge turns that into NULL, which would leave
     * EllipticK[1.] unevaluated rather than ComplexInfinity. */
    if (ell_is_one(m))  return expr_new_symbol(SYM_ComplexInfinity);
    if (ell_is_zero(m) && !arg_is_inexact(m)) return ell_half_pi();   /* K[0] = Pi/2 */
    /* K(m) ~ (Pi - I Log[16 m]) / (2 Sqrt[m]) for large m, so the magnitude
     * vanishes: K[Infinity] = 0. */
    if (ell_is_infinite(m)) return expr_new_integer(0);
    { int64_t kn, kd;
      if (is_rational(m, &kn, &kd) && kn == 1 && kd == 2) return ell_k_closed_half();
      if (m->type == EXPR_INTEGER && m->data.integer == -1) return ell_k_closed_minus_one(); }

    if (ell_any_inexact(&m, 1)) {
        double mv[1], r;
        if (ell_machine_args(&m, 1, mv) && elliptic_machine_k(mv[0], &r))
            return expr_new_real(r);
        Expr* v = flint_num_elliptic_k(m);
        if (v) return v;
    }
    if (ell_is_zero(m)) return ell_half_pi();
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

/* The numeric path runs BEFORE the remaining exact reductions, so that an
 * inexact argument gets an inexact answer. The rest of the system is consistent
 * about this -- Sin[0.] is 0., Gamma[1.] is 1., Cos[0.] is 1. -- and these
 * reductions were not: EllipticK[0.] returned the exact Pi/2 (head Times), and
 * EllipticPi[0., 1/2] returned the SYMBOLIC EllipticK[1/2], so a numeric call
 * came back un-numeric. Arb answers all of those correctly at the argument's own
 * precision, so the fix is an ordering, not a conversion.
 *
 * The exact reductions still sit underneath as the fallback, and they are still
 * reached when the numeric path declines -- which is exactly what happens for a
 * spelling Arb cannot read, such as the exact Pi/2 upper limit (scalar_to_arb
 * takes numbers, not a Times). So EllipticE[Pi/2, 0.5] still reduces to
 * EllipticE[0.5] and then evaluates. */
    if (ell_any_inexact(res->data.function.args, 2)) {
        double a2[2], r;
        if (ell_machine_args(res->data.function.args, 2, a2)
            && elliptic_machine_f(a2[0], a2[1], &r))
            return expr_new_real(r);
        Expr* v = flint_num_elliptic_f(phi, m);
        if (v) return v;
    }

    if (ell_is_zero(phi)) return expr_new_integer(0);           /* F[0, m] = 0 */
    if (ell_is_zero(m))   return expr_copy(phi);                /* F[phi, 0] = phi */
    if (ell_is_infinite(m)) return expr_new_integer(0);         /* F[phi, Infinity] = 0 */
    if (ell_is_half_pi(phi)) {                                  /* F[Pi/2, m] = K[m] */
        Expr* a[1] = { expr_copy(m) };
        return ell_call("EllipticK", a, 1);
    }
    { Expr* o = ell_odd_in_phi("EllipticF", res->data.function.args, 2, 0);
      if (o) return o; }
    return NULL;
}

/* ------------------------------------------------------------------ */
/*  EllipticE[m] / EllipticE[phi, m]                                   */
/* ------------------------------------------------------------------ */

static Expr* elliptice_complete(Expr* m) {
    if (ell_any_inexact(&m, 1)) {               /* see the note in builtin_ellipticf */
        double mv[1], r;
        if (ell_machine_args(&m, 1, mv) && elliptic_machine_e_complete(mv[0], &r))
            return expr_new_real(r);
        Expr* v = flint_num_elliptic_e(m);
        if (v) return v;
    }
    if (ell_is_zero(m)) return ell_half_pi();                   /* E[0] = Pi/2 */
    if (ell_is_one(m))  return expr_new_integer(1);             /* E[1] = 1    */
    /* E(m) ~ I Sqrt[m] for large m, so the magnitude diverges. */
    if (ell_is_infinite(m)) return expr_new_symbol(SYM_ComplexInfinity);
    return NULL;
}

static Expr* elliptice_incomplete(Expr* phi, Expr* m) {
    Expr* two[2] = { phi, m };
    if (ell_any_inexact(two, 2)) {              /* see the note in builtin_ellipticf */
        double a2[2], r;
        if (ell_machine_args(two, 2, a2) && elliptic_machine_e_inc(a2[0], a2[1], &r))
            return expr_new_real(r);
        Expr* v = flint_num_elliptic_e_inc(phi, m);
        if (v) return v;
    }

    if (ell_is_zero(phi)) return expr_new_integer(0);           /* E[0, m] = 0 */
    if (ell_is_zero(m))   return expr_copy(phi);                /* E[phi, 0] = phi */
    if (ell_is_one(m) && ell_in_principal_strip(phi)) {  /* E[phi, 1] = Sin[phi] there */
        Expr* a[1] = { expr_copy(phi) };
        return ell_call("Sin", a, 1);
    }
    if (ell_is_half_pi(phi)) {                                  /* E[Pi/2, m] = E[m] */
        Expr* a[1] = { expr_copy(m) };
        return ell_call("EllipticE", a, 1);
    }
    if (ell_is_infinite(m) || ell_is_infinite(phi))
        return expr_new_symbol(SYM_ComplexInfinity);
    { Expr* o = ell_odd_in_phi("EllipticE", two, 2, 0);
      if (o) return o; }
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
    /* n = 1 is a pole of the COMPLETE form only: the integrand carries a
     * 1/Cos[t]^2 and the upper limit is Pi/2. Checked before the numeric path
     * for the same reason as EllipticK's -- Arb returns a non-finite ball and
     * the bridge turns that into NULL, which would leave EllipticPi[1., 1/2]
     * unevaluated. (Verified divergent: Pi[1-e|1/2] grows as 1/Sqrt[e] --
     * 21.5, 221, 2221, 22214 at e = 10^-2, 10^-4, 10^-6, 10^-8.) The INCOMPLETE
     * form has no such pole: Pi[1, 1, 1/2] is 1.73199154202. */
    if (ell_is_one(n)) return expr_new_symbol(SYM_ComplexInfinity);

    Expr* two[2] = { n, m };
    if (ell_any_inexact(two, 2)) {              /* see the note in builtin_ellipticf */
        double a2[2], r;
        if (ell_machine_args(two, 2, a2) && elliptic_machine_pi(a2[0], a2[1], &r))
            return expr_new_real(r);
        Expr* v = flint_num_elliptic_pi(n, m);
        if (v) return v;
    }
    if (ell_is_zero(n)) {                                       /* Pi[0, m] = K[m] */
        Expr* a[1] = { expr_copy(m) };
        return ell_call("EllipticK", a, 1);
    }
    if (ell_is_zero(m)) return ell_pi_closed_m_zero(n);         /* Pi[n, 0] */
    if (ell_is_infinite(n) || ell_is_infinite(m)) return expr_new_integer(0);
    return NULL;
}

static Expr* ellipticpi_incomplete(Expr* n, Expr* phi, Expr* m) {
    Expr* three[3] = { n, phi, m };
    if (ell_any_inexact(three, 3)) {            /* see the note in builtin_ellipticf */
        double a3[3], r;
        if (ell_machine_args(three, 3, a3)
            && elliptic_machine_pi_inc(a3[0], a3[1], a3[2], &r))
            return expr_new_real(r);
        Expr* v = flint_num_elliptic_pi_inc(n, phi, m);
        if (v) return v;
    }

    if (ell_is_zero(phi)) return expr_new_integer(0);           /* Pi[n, 0, m] = 0 */
    if (ell_is_zero(n)) {                                       /* Pi[0, phi, m] = F */
        Expr* a[2] = { expr_copy(phi), expr_copy(m) };
        return ell_call("EllipticF", a, 2);
    }
    if (ell_is_half_pi(phi)) {                                  /* Pi[n, Pi/2, m] */
        Expr* a[2] = { expr_copy(n), expr_copy(m) };
        return ell_call("EllipticPi", a, 2);
    }
    if (ell_is_zero(m)) return ell_pi_closed_m_zero_inc(n, phi); /* Pi[n, phi, 0] */
    if (ell_is_infinite(n) || ell_is_infinite(m)) return expr_new_integer(0);
    { Expr* o = ell_odd_in_phi("EllipticPi", three, 3, 1);
      if (o) return o; }
    return NULL;
}

Expr* builtin_ellipticpi(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr* v;
    if (argc == 2)      v = ellipticpi_complete(res->data.function.args[0],
                                                res->data.function.args[1]);
    else if (argc == 3) v = ellipticpi_incomplete(res->data.function.args[0],
                                                  res->data.function.args[1],
                                                  res->data.function.args[2]);
    else return ell_argt("EllipticPi", argc, 2, 3);
    if (v) return v;

    /* A VISIBLE NDArray that no kernel consumed reaches this builtin untouched
     * -- where, before this, it was simply left unevaluated in all three
     * argument positions. Per CLAUDE.md an unevaluated visible NDArray is a
     * wrong answer, not a slow one, so delist and re-evaluate: the List path
     * answers correctly, element by element.
     *
     * Still load-bearing now that the R_J kernels exist, and MORE so. The
     * element-wise ND layer tops out at arity 2, so the three-argument form has
     * no kernel to dispatch to -- yet packed_aware is a property of the SYMBOL,
     * not of one arity, so registering the two-argument kernel also stopped the
     * transparency gate materialising packed Lists here. This fallback is what
     * keeps that correct, which is why it had to land before REG_B. */
    for (size_t i = 0; i < argc; i++)
        if (is_ndarray(res->data.function.args[i]))
            return ndarray_delist_and_reeval(res);
    return NULL;
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
        "Pi/2 and EllipticE[1] is 1. EllipticE[phi, 1] is Sin[phi] only for "
        "|phi| <= Pi/2, since E(phi|1) is the integral of Abs[Cos[t]]: "
        "EllipticE[2, 1] is 2 - Sin[2], not Sin[2].");

    symtab_add_builtin("EllipticPi", builtin_ellipticpi);
    symtab_get_def("EllipticPi")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    symtab_set_docstring("EllipticPi",
        "EllipticPi[n, m] is the complete elliptic integral of the third kind, "
        "Integrate[1/((1 - n Sin[t]^2) Sqrt[1 - m Sin[t]^2]), {t, 0, Pi/2}], "
        "and EllipticPi[n, phi, m] the incomplete one, with upper limit phi. "
        "The parameter argument is m = k^2, not the modulus k. "
        "EllipticPi[0, m] is EllipticK[m] and EllipticPi[0, phi, m] is "
        "EllipticF[phi, m]; EllipticPi[1, m] is ComplexInfinity. Where the path "
        "crosses the pole at Sin[t]^2 == 1/n the value is complex, not a real "
        "principal value: EllipticPi[3/2, 1/2] is "
        "-0.4567203134529101 - 2.720699046351328 I.");
}
