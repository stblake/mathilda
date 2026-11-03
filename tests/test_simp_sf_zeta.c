#include "simp_sf_test_util.h"
#include "symtab.h"
#include "core.h"

/*
 * Simplify, special-function engine, zeta module (src/simp/simp_sf_zeta.c):
 * PolyGamma, HarmonicNumber, HurwitzZeta / Zeta[s, a], LerchPhi, the
 * Riemann functional equation, PolyLog.
 *
 * Derived identities (docs/design/simp_special_functions.md, sections 4-5),
 * each also checked numerically at two real and one complex point:
 *
 *   psi^(n)(z+1) = psi^(n)(z) + (-1)^n n! z^(-n-1)     (n-th derivative of
 *       psi(z+1) = psi(z) + 1/z, the log-derivative of Gamma(z+1) = z Gamma(z))
 *   psi^(n)(1-w) = (-1)^n (psi^(n)(w) + Pi D^n Cot[Pi w])   (from reflection)
 *   psi(m z) = Log m + (1/m) Sum_j psi(z + j/m);
 *   psi^(n)(m z) = m^(-n-1) Sum_j psi^(n)(z + j/m)  (n >= 1)  (multiplication)
 *   H(n+1) = H(n) + 1/(n+1);  H(n) = psi(n+1) + EulerGamma;
 *   H(n+1, r) = H(n, r) + (n+1)^-r
 *   zeta_H(s, a+1) = zeta_H(s, a) - a^-s        (first term of the series)
 *   Zeta(s, a+1) = Zeta(s, a) - (a^2)^(-s/2)    (Zeta[s, a] sums ((k+a)^2)^(-s/2))
 *   Phi(z, s, a) = a^-s + z Phi(z, s, a+1)      (Lerch series)
 *   zeta(1-w) = 2^(1-w) Pi^-w Cos[Pi w/2] Gamma(w) zeta(w)   (Riemann)
 *   Li2(u) + Li2(1-u) = Pi^2/6 - Log u Log(1-u)       (Euler reflection)
 *   Li2(u) + Li2(u/(u-1)) = -Log(1-u)^2/2             (Landen; u < 1)
 *   Li_n(u) + Li_n(-u) = 2^(1-n) Li_n(u^2)            (series, any order)
 *   Li2(u) + Li2(1/u) = -Pi^2/6 - Log(-u)^2/2   (u < 0; fails on (0, Inf))
 */

static char g_buf[8192];
#define FMT(...) (snprintf(g_buf, sizeof g_buf, __VA_ARGS__), g_buf)

static const char* Z_ARGS[] = { "x", "x + 1/3", "2 x", "a + b", "x/2" };
#define N_Z (sizeof Z_ARGS / sizeof Z_ARGS[0])

/* negative-real points for the inversion (u < 0) */
static const char* NEG_PTS[] = { "{x -> -17/10}", "{x -> -3/10}", "{x -> -5}" };
/* -1 < z < 1/2: every shape u(z) the Landen cases use stays below 1 there
 * (z + 1/10 < 3/5, z^2 < 1, ...), as the identity requires */
static const char* LT1_PTS[] = { "{z -> 2/5}", "{z -> -1/2}", "{z -> 3/10}" };

/* ---- PolyGamma ------------------------------------------------------------ */

/* Sum_{j=0}^{k-1} (-1)^n n! (z + j)^(-n-1), as source text */
static const char* pg_sum(int n, const char* z, int k) {
    static char sb[2048];
    sb[0] = 0;
    for (int j = 0; j < k; j++) {
        char t[256];
        snprintf(t, sizeof t, "%s(-1)^%d %d! (%s + %d)^(-%d)", j ? " + " : "", n, n, z, j, n + 1);
        strcat(sb, t);
    }
    return sb;
}

static void test_polygamma_recurrence(void) {
    sf_family("PolyGamma");
    for (int n = 0; n <= 3; n++)
        for (size_t i = 0; i < N_Z; i++)
            for (int k = 1; k <= 3; k++) {
                char sum[2048];
                snprintf(sum, sizeof sum, "%s", pg_sum(n, Z_ARGS[i], k));
                sf_zero(FMT("PolyGamma[%d, %s + %d] - PolyGamma[%d, %s] - (%s)", n, Z_ARGS[i], k, n, Z_ARGS[i], sum));
            }
    sf_to("PolyGamma[0, x + 1] - PolyGamma[0, x]", "1/x");
    sf_to("PolyGamma[2, x + 1] - PolyGamma[2, x]", "2/x^3");
    sf_to("PolyGamma[1, x + 2] - PolyGamma[1, x] - 1/x^2 - 1/(x + 1)^2", "-2/x^2 - 2/(1 + x)^2");
    sf_zero("PolyGamma[x + 1] - PolyGamma[x] - 1/x");
    sf_zero("PolyGamma[0, x] - PolyGamma[0, x - 1] - 1/(x - 1)");
    sf_zero("D[LogGamma[x + 1], x] - PolyGamma[0, x] - 1/x");
    sf_identity("PolyGamma[3, x + 1] - PolyGamma[3, x] + 6/x^4");
}

static void test_polygamma_reflection(void) {
    sf_family("PolyGamma");
    for (size_t i = 0; i < N_Z; i++) {
        sf_zero(FMT("PolyGamma[0, 1 - (%s)] - PolyGamma[0, %s] - Pi Cot[Pi (%s)]", Z_ARGS[i], Z_ARGS[i], Z_ARGS[i]));
        sf_zero(FMT("PolyGamma[1, 1 - (%s)] + PolyGamma[1, %s] - Pi^2 Csc[Pi (%s)]^2", Z_ARGS[i], Z_ARGS[i], Z_ARGS[i]));
        sf_zero(FMT("PolyGamma[2, 1 - (%s)] - PolyGamma[2, %s] - 2 Pi^3 Cot[Pi (%s)] Csc[Pi (%s)]^2",
                    Z_ARGS[i], Z_ARGS[i], Z_ARGS[i], Z_ARGS[i]));
    }
    sf_to("PolyGamma[0, 1 - x] - PolyGamma[0, x]", "Pi Cot[Pi x]");
    sf_to("PolyGamma[1, 1 - x] + PolyGamma[1, x]", "Pi^2 Csc[Pi x]^2");
    sf_zero("PolyGamma[0, -x] - PolyGamma[0, x + 1] - Pi Cot[Pi (x + 1)]");
    sf_identity("PolyGamma[0, 1 - x] - PolyGamma[0, x] - Pi Cot[Pi x]");
}

static void test_polygamma_multiplication(void) {
    sf_family("PolyGamma");
    static const char* Z[] = { "x", "x + 1/3", "a" };
    for (size_t i = 0; i < 3; i++) {
        sf_zero(FMT("PolyGamma[0, 2 (%s)] - Log[2] - (PolyGamma[0, %s] + PolyGamma[0, %s + 1/2])/2", Z[i], Z[i], Z[i]));
        sf_zero(FMT("PolyGamma[1, 2 (%s)] - (PolyGamma[1, %s] + PolyGamma[1, %s + 1/2])/4", Z[i], Z[i], Z[i]));
        sf_zero(FMT("PolyGamma[0, 3 (%s)] - Log[3] - (PolyGamma[0, %s] + PolyGamma[0, %s + 1/3] + PolyGamma[0, %s + 2/3])/3",
                    Z[i], Z[i], Z[i], Z[i]));
    }
    sf_identity("PolyGamma[0, 2 x] - Log[2] - (PolyGamma[0, x] + PolyGamma[0, x + 1/2])/2");
}

/* ---- HarmonicNumber -------------------------------------------------------- */

static void test_harmonic(void) {
    sf_family("HarmonicNumber");
    static const char* N[] = { "n", "2 n", "n + 1/2", "a + b", "n/3", "m" };
    for (size_t i = 0; i < 6; i++) {
        for (int k = 1; k <= 4; k++) {
            /* H(N + k) - H(N) = Sum_{j=1..k} 1/(N + j) */
            char sum[1024] = "";
            for (int j = 1; j <= k; j++) {
                char t[128];
                snprintf(t, sizeof t, "%s1/(%s + %d)", j > 1 ? " + " : "", N[i], j);
                strcat(sum, t);
            }
            sf_zero(FMT("HarmonicNumber[%s + %d] - HarmonicNumber[%s] - (%s)", N[i], k, N[i], sum));
        }
        for (int k = 1; k <= 2; k++) {
            /* H(N) - H(N - k) = Sum_{j=0..k-1} 1/(N - j) */
            char sum[1024] = "";
            for (int j = 0; j < k; j++) {
                char t[128];
                snprintf(t, sizeof t, "%s1/(%s - %d)", j ? " + " : "", N[i], j);
                strcat(sum, t);
            }
            sf_zero(FMT("HarmonicNumber[%s] - HarmonicNumber[%s - %d] - (%s)", N[i], N[i], k, sum));
        }
        /* H(N) = psi(N + 1) + EulerGamma */
        sf_zero(FMT("HarmonicNumber[%s] - PolyGamma[0, %s + 1] - EulerGamma", N[i], N[i]));
    }
    static const char* NN[] = { "n", "m" };
    for (size_t i = 0; i < 2; i++)
        for (int r = 2; r <= 4; r++)
            for (int k = 1; k <= 3; k++) {
                char sum[1024] = "";
                for (int j = 1; j <= k; j++) {
                    char t[128];
                    snprintf(t, sizeof t, "%s(%s + %d)^(-%d)", j > 1 ? " + " : "", NN[i], j, r);
                    strcat(sum, t);
                }
                sf_zero(FMT("HarmonicNumber[%s + %d, %d] - HarmonicNumber[%s, %d] - (%s)", NN[i], k, r, NN[i], r, sum));
            }
    sf_to("HarmonicNumber[n + 1] - HarmonicNumber[n]", "1/(1 + n)");
    sf_to("HarmonicNumber[n] - PolyGamma[0, n + 1]", "EulerGamma");
    sf_zero("HarmonicNumber[n, s] - HarmonicNumber[n - 1, s] - n^(-s)");
    sf_unchanged("HarmonicNumber[n] + HarmonicNumber[m]");
    sf_identity("HarmonicNumber[n] - PolyGamma[0, n + 1] - EulerGamma");
}

/* ---- Hurwitz zeta, Zeta[s, a], LerchPhi ----------------------------------- */

static const char* A_SH[] = { "a", "a + 1/2", "2 a", "a + b", "b", "a/3" };
#define N_A (sizeof A_SH / sizeof A_SH[0])

static void test_hurwitz(void) {
    sf_family("HurwitzZeta");
    for (size_t i = 0; i < N_A; i++) {
        for (int k = 1; k <= 4; k++) {
            char sum[1024] = "";
            for (int j = 0; j < k; j++) {
                char t[128];
                snprintf(t, sizeof t, "%s(%s + %d)^(-s)", j ? " + " : "", A_SH[i], j);
                strcat(sum, t);
            }
            sf_zero(FMT("HurwitzZeta[s, %s + %d] - HurwitzZeta[s, %s] + (%s)", A_SH[i], k, A_SH[i], sum));
        }
        for (int k = 1; k <= 3; k++) {
            char sum[1024] = "";
            for (int j = 1; j <= k; j++) {
                char t[128];
                snprintf(t, sizeof t, "%s(%s - %d)^(-s)", j > 1 ? " + " : "", A_SH[i], j);
                strcat(sum, t);
            }
            sf_zero(FMT("HurwitzZeta[s, %s] - HurwitzZeta[s, %s - %d] + (%s)", A_SH[i], A_SH[i], k, sum));
        }
    }
    for (size_t i = 0; i < 4; i++) {
        sf_zero(FMT("HurwitzZeta[3, %s + 1] - HurwitzZeta[3, %s] + (%s)^-3", A_SH[i], A_SH[i], A_SH[i]));
        sf_zero(FMT("HurwitzZeta[5/2, %s + 1] - HurwitzZeta[5/2, %s] + (%s)^(-5/2)", A_SH[i], A_SH[i], A_SH[i]));
    }
    sf_to("HurwitzZeta[s, a + 1] - HurwitzZeta[s, a]", "-1/a^s");
    sf_unchanged("HurwitzZeta[s, a] + HurwitzZeta[t, a + 1]");
    sf_identity("HurwitzZeta[s, a + 1] - HurwitzZeta[s, a] + a^(-s)");
}

static void test_zeta2(void) {
    sf_family("Zeta[s,a]");
    /* one-step identities with the convention's (a^2)^(-s/2) term; a - 1 < 0
     * at the second real point exercises the Re a < 0 branch */
    static const char* S_SH[] = { "s", "3", "5/2", "2 s" };
    for (size_t i = 0; i < N_A; i++)
        for (size_t j = 0; j < 4; j++) {
            sf_zero(FMT("Zeta[%s, %s + 1] - Zeta[%s, %s] + ((%s)^2)^(-(%s)/2)", S_SH[j], A_SH[i], S_SH[j], A_SH[i], A_SH[i], S_SH[j]));
            sf_zero(FMT("Zeta[%s, %s] - Zeta[%s, %s - 1] + ((%s - 1)^2)^(-(%s)/2)", S_SH[j], A_SH[i], S_SH[j], A_SH[i], A_SH[i], S_SH[j]));
        }
    sf_to("Zeta[s, a + 1] - Zeta[s, a]", "-1/(a^2)^(1/2 s)");
    sf_unchanged("Zeta[s, a] + Zeta[s, b]");
    sf_identity("Zeta[s, a] - Zeta[s, a - 1] + ((a - 1)^2)^(-s/2)");
}

static void test_lerch(void) {
    sf_family("LerchPhi");
    static const char* ZZ[] = { "z", "z/2" };
    for (size_t i = 0; i < 5; i++)
        for (size_t zi = 0; zi < 2; zi++) {
            for (int k = 1; k <= 3; k++) {
                /* Phi(Z, s, A+k) = Z^-k Phi(Z, s, A) - Sum_{j<k} Z^(j-k) (A+j)^-s */
                char sum[1024] = "";
                for (int j = 0; j < k; j++) {
                    char t[160];
                    snprintf(t, sizeof t, "%s(%s)^(%d) (%s + %d)^(-s)", j ? " + " : "", ZZ[zi], j - k, A_SH[i], j);
                    strcat(sum, t);
                }
                sf_zero(FMT("LerchPhi[%s, s, %s + %d] - (%s)^(-%d) LerchPhi[%s, s, %s] + (%s)",
                            ZZ[zi], A_SH[i], k, ZZ[zi], k, ZZ[zi], A_SH[i], sum));
            }
            /* down one step: Phi(Z, s, A - 1) = (A-1)^-s + Z Phi(Z, s, A) */
            sf_zero(FMT("LerchPhi[%s, s, %s - 1] - (%s - 1)^(-s) - (%s) LerchPhi[%s, s, %s]",
                        ZZ[zi], A_SH[i], A_SH[i], ZZ[zi], ZZ[zi], A_SH[i]));
        }
    for (size_t i = 0; i < 6; i++)
        sf_zero(FMT("z LerchPhi[z, 2, %s + 1] - LerchPhi[z, 2, %s] + (%s)^-2", A_SH[i], A_SH[i], A_SH[i]));
    sf_to("z LerchPhi[z, s, a + 1] - LerchPhi[z, s, a]", "-1/a^s");
    sf_zero("LerchPhi[z, 3, a + 1] z - LerchPhi[z, 3, a] + a^-3");
    sf_unchanged("LerchPhi[z, s, a] + LerchPhi[w, s, a + 1]");
    sf_identity("z LerchPhi[z, s, a + 1] - LerchPhi[z, s, a] + a^(-s)");
}

/* ---- Riemann functional equation ----------------------------------------- */

static void test_riemann(void) {
    sf_family("Zeta (Riemann)");
    /* shapes stay off the integers at every sample point (s = 9/4, 11/3, ...):
     * at an integer w, Sin[Pi w/2] Gamma[1 - w] is 0 times a pole */
    static const char* W[] = { "s", "s + 1", "2 s", "s/2", "s + 1/5", "a + b",
                               "3 s + 1/7", "s - 1/6", "s + 2", "nu", "2 nu", "nu + 1/2" };
    for (size_t i = 0; i < 12; i++) {
        const char* w = W[i];
        sf_zero(FMT("Zeta[1 - (%s)] - 2^(1 - (%s)) Pi^(-(%s)) Cos[Pi (%s)/2] Gamma[%s] Zeta[%s]", w, w, w, w, w, w));
        sf_zero(FMT("Zeta[%s] - 2^(%s) Pi^((%s) - 1) Sin[Pi (%s)/2] Gamma[1 - (%s)] Zeta[1 - (%s)]", w, w, w, w, w, w));
        sf_zero(FMT("x (Zeta[1 - (%s)] - 2^(1 - (%s)) Pi^(-(%s)) Cos[Pi (%s)/2] Gamma[%s] Zeta[%s])", w, w, w, w, w, w));
        sf_zero(FMT("(Zeta[1 - (%s)] - 2^(1 - (%s)) Pi^(-(%s)) Cos[Pi (%s)/2] Gamma[%s] Zeta[%s]) Zeta[%s]", w, w, w, w, w, w, w));
    }
    sf_zero("Zeta[-s] - 2^(-s) Pi^(-s - 1) Sin[-Pi s/2] Gamma[1 + s] Zeta[1 + s]");
    sf_unchanged("Zeta[s] + Zeta[t]");
    sf_identity("Zeta[s] - 2^s Pi^(s - 1) Sin[Pi s/2] Gamma[1 - s] Zeta[1 - s]");
}

/* ---- PolyLog --------------------------------------------------------------- */

static void test_polylog(void) {
    sf_family("PolyLog");
    static const char* U[] = { "z", "z/2", "z^2", "z + 1/10", "2 z/3", "z - 1/10", "z/3", "3 z/4" };
    for (size_t i = 0; i < 8; i++) {
        sf_zero(FMT("PolyLog[2, %s] + PolyLog[2, 1 - (%s)] - Pi^2/6 + Log[%s] Log[1 - (%s)]", U[i], U[i], U[i], U[i]));
        sf_zero_assuming(FMT("PolyLog[2, %s] + PolyLog[2, (%s)/((%s) - 1)] + Log[1 - (%s)]^2/2", U[i], U[i], U[i], U[i]),
                         "-1 < z < 1/2", LT1_PTS, 3);
    }
    for (int n = 2; n <= 5; n++)
        for (size_t i = 0; i < 5; i++)
            sf_zero(FMT("PolyLog[%d, (%s)^2] - 2^%d (PolyLog[%d, %s] + PolyLog[%d, -(%s)])",
                        n, U[i], n - 1, n, U[i], n, U[i]));
    static const char* NEG[] = { "x", "2 x", "x/3", "x - 1", "3 x - 2" };
    for (size_t i = 0; i < 5; i++)
        sf_zero_assuming(FMT("PolyLog[2, %s] + PolyLog[2, 1/(%s)] + Pi^2/6 + Log[-(%s)]^2/2", NEG[i], NEG[i], NEG[i]),
                         "x < 0", NEG_PTS, 3);
    static const char* C[] = { "13/20", "1/3", "1/4", "2/5", "1/10" };
    for (size_t i = 0; i < 5; i++)
        sf_zero(FMT("PolyLog[2, %s] + PolyLog[2, 1 - %s] - Pi^2/6 + Log[%s] Log[1 - %s]", C[i], C[i], C[i], C[i]));
    sf_zero("PolyLog[s, z^2] - 2^(s - 1) (PolyLog[s, z] + PolyLog[s, -z])");
    /* the reference's Simplify takes the content out (v0.6147) */
    sf_to("PolyLog[2, z] + PolyLog[2, 1 - z]", "1/6 (Pi^2 - 6 Log[1 - z] Log[z])");
    sf_to_assuming("PolyLog[2, z] + PolyLog[2, z/(z - 1)]", "z < 1", "-1/2 Log[1 - z]^2", LT1_PTS, 3);
    sf_unchanged("PolyLog[2, z] + PolyLog[2, z/(z - 1)]");    /* Landen fails on (1, Inf) */
    sf_unchanged("PolyLog[2, z] + PolyLog[2, 1/z]");          /* inversion fails on (0, Inf) */
    sf_unchanged("PolyLog[2, z] + PolyLog[3, 1 - z]");        /* different orders */
    sf_identity("PolyLog[2, z] + PolyLog[2, 1 - z] - Pi^2/6 + Log[z] Log[1 - z]");
}

/* ---- soundness ------------------------------------------------------------ */

static void test_unchanged(void) {
    sf_family("PolyGamma");
    sf_unchanged("PolyGamma[0, x] + PolyGamma[0, y]");
    sf_unchanged("PolyGamma[0, x] + PolyGamma[1, x + 1]");    /* different orders */
    sf_unchanged("PolyGamma[n, x + 1] - PolyGamma[n, x]");    /* symbolic order: no recurrence */
    sf_unchanged("PolyGamma[1/2, x + 1] - PolyGamma[1/2, x]");/* non-integer order */
    sf_unchanged("HarmonicNumber[n] + HarmonicNumber[m]");
    sf_unchanged("HurwitzZeta[s, a] + HurwitzZeta[t, a + 1]");
    sf_unchanged("LerchPhi[z, s, a] + LerchPhi[w, s, a + 1]");
    sf_sound("PolyGamma[0, x + 2] + PolyGamma[0, 1 - x]");
    sf_sound("HarmonicNumber[n + 3] - 2 HarmonicNumber[n]");
    sf_sound("HurwitzZeta[s, a + 2] + HurwitzZeta[s, a]");
}

static void test_fast(void) {
    sf_family("PolyGamma");
    sf_fast("PolyGamma[0, x + 30] - PolyGamma[0, x]", 5.0);
    sf_fast("Sum[HarmonicNumber[n + k], {k, 0, 8}]", 5.0);
    sf_fast("PolyLog[2, z] + PolyLog[2, 1 - z] + PolyLog[2, z/(z - 1)]", 5.0);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_polygamma_recurrence);
    TEST(test_polygamma_reflection);
    TEST(test_polygamma_multiplication);
    TEST(test_harmonic);
    TEST(test_hurwitz);
    TEST(test_zeta2);
    TEST(test_lerch);
    TEST(test_riemann);
    TEST(test_polylog);
    TEST(test_unchanged);
    TEST(test_fast);
    return sf_finish("simp_sf_zeta", 1600);
}
