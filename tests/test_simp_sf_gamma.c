#include "simp_sf_test_util.h"
#include "symtab.h"
#include "core.h"

/*
 * Simplify, special-function engine, Gamma module (src/simp/simp_sf_gamma.c).
 *
 * Every identity is DERIVED (docs/design/simp_special_functions.md, section
 * Gamma), and every case is checked numerically at two real and one complex
 * point (simp_sf_test_util.h):
 *
 *   Gamma(z+1) = z Gamma(z)                (Euler integral, by parts)
 *     => Gamma(z+k) = (z)_k Gamma(z), (z)_k = z (z+1) ... (z+k-1)
 *   Gamma(a+1, z) = a Gamma(a, z) + z^a e^-z    (incomplete, by parts)
 *   LogGamma(z+1) = LogGamma(z) + Log(z)   (on the slit plane)
 *   G(z+1) = Gamma(z) G(z)                 (Barnes G functional equation)
 *   Gamma(w) Gamma(1-w) = Pi / Sin(Pi w)   (Euler reflection)
 *     => Gamma(x+k) Gamma(1-x-k) = (-1)^k Pi / Sin(Pi x)
 *   Gamma(2z) = 2^(2z-1) Pi^(-1/2) Gamma(z) Gamma(z+1/2)   (Legendre)
 *   Gamma(m z) = (2 Pi)^((1-m)/2) m^(m z-1/2) Prod_j Gamma(z + j/m)  (Gauss)
 *   (a)_n = Gamma(a+n)/Gamma(a);  B(a,b) = Gamma(a) Gamma(b)/Gamma(a+b);
 *   C(n,k) = Gamma(n+1)/(Gamma(k+1) Gamma(n-k+1));  n! = Gamma(n+1)
 *     => (a)_{n+k} = (a)_n (a+n)_k,   B(a+k,b) = (a)_k/(a+b)_k B(a,b),
 *        C(n+1,k) = C(n,k) + C(n,k-1) (Pascal), C(n,k) = C(n,n-k),
 *        k C(n,k) = n C(n-1,k-1), C(n,k) C(k,j) = C(n,j) C(n-j,k-j).
 *
 * The generated blocks instantiate one derived identity over many shifts
 * and argument shapes; each instance is a separate case.
 */

static char g_buf[8192];
#define FMT(...) (snprintf(g_buf, sizeof g_buf, __VA_ARGS__), g_buf)

/* argument shapes for a symbolic Gamma argument */
static const char* S_ARGS[] = { "x", "2 x", "x + 1/3", "a + b", "x/2", "a - 1/4" };
#define N_S_ARGS (sizeof S_ARGS / sizeof S_ARGS[0])

/* ---- Gamma: recurrence ------------------------------------------------ */

static void test_gamma_recurrence_generated(void) {
    sf_family("Gamma");
    /* up: Gamma[s + k] - Pochhammer[s, k] Gamma[s] (Pochhammer of integer k
     * expands to the product, so the right side is the derived polynomial) */
    for (size_t i = 0; i < N_S_ARGS; i++)
        for (int k = 1; k <= 4; k++)
            sf_zero(FMT("Gamma[%s + %d] - Pochhammer[%s, %d] Gamma[%s]", S_ARGS[i], k, S_ARGS[i], k, S_ARGS[i]));
    /* down: Gamma[s - k] Pochhammer[s - k, k] - Gamma[s] */
    for (size_t i = 0; i < 4; i++)
        for (int k = 1; k <= 3; k++)
            sf_zero(FMT("Gamma[%s - %d] Pochhammer[%s - %d, %d] - Gamma[%s]",
                        S_ARGS[i], k, S_ARGS[i], k, k, S_ARGS[i]));
}

static void test_gamma_recurrence_forms(void) {
    sf_family("Gamma");
    sf_zero("Gamma[x + 1] - x Gamma[x]");
    sf_to("Gamma[x + 2]/Gamma[x]", "x (1 + x)");
    sf_to("Gamma[x + 5]/Gamma[x + 2]", "(2 + x) (3 + x) (4 + x)");
    sf_to("Gamma[x + 1] + Gamma[x]", "(1 + x) Gamma[x]");
    sf_to("Gamma[x + 1/2]/Gamma[x - 1/2]", "-1/2 + x");
    sf_to("Gamma[n + 1]/Gamma[n]", "n");
    /* (a span of 500 is the sf_fast guard; numerically Gamma[x + 500] at 30
     * digits leaves ~1e-14 absolute on a value ~2.5e5) */
    sf_to("Gamma[x + 50]/Gamma[x + 48]", "(48 + x) (49 + x)");
    sf_zero("Gamma[x + 3] - (x + 2) (x + 1) x Gamma[x]");
    sf_zero("x Gamma[x] - Gamma[x + 1] + Gamma[y] - Gamma[y]");
    sf_zero("Gamma[x + 1]^2 - x^2 Gamma[x]^2");
    sf_zero("Gamma[x + 2] Gamma[y] - x (x + 1) Gamma[x] Gamma[y]");
    sf_zero("Gamma[x + 1]/Gamma[x] - Gamma[y + 1]/Gamma[y] - x + y");
    sf_zero("(Gamma[x + 2] - Gamma[x + 1])/Gamma[x + 1] - x");
    sf_identity("Gamma[x + 3] - (x + 2) (x + 1) x Gamma[x]");
}

/* ---- Gamma: reflection ----------------------------------------------------- */

static void test_gamma_reflection(void) {
    sf_family("Gamma");
    for (size_t i = 0; i < N_S_ARGS; i++)
        sf_zero(FMT("Gamma[%s] Gamma[1 - (%s)] - Pi/Sin[Pi (%s)]", S_ARGS[i], S_ARGS[i], S_ARGS[i]));
    /* shifted pairs: Gamma[x + k] Gamma[1 - x - k] = (-1)^k Pi / Sin[Pi x] */
    for (int k = -3; k <= 3; k++)
        sf_zero(FMT("Gamma[x + %d] Gamma[1 - x - %d] - (-1)^%d Pi/Sin[Pi x]", k, k, k < 0 ? -k : k));
    /* numeric classes */
    static const char* C[] = { "1/3", "1/4", "1/5", "1/6", "2/5", "3/7", "1/8" };
    for (size_t i = 0; i < sizeof C / sizeof C[0]; i++)
        sf_zero(FMT("Gamma[%s] Gamma[1 - %s] - Pi/Sin[Pi %s]", C[i], C[i], C[i]));
    sf_to("Gamma[x] Gamma[1 - x]", "Pi Csc[Pi x]");
    sf_to("Gamma[1/3] Gamma[2/3]", "(2 Pi)/Sqrt[3]");
    sf_to("Gamma[x + 1/2] Gamma[1/2 - x]", "Pi Sec[Pi x]");
    sf_zero("Gamma[x] Gamma[1 - x] Sin[Pi x] - Pi");
    sf_zero("Gamma[x + 1] Gamma[1 - x] - Pi x/Sin[Pi x]");
    sf_zero("Gamma[x] Gamma[-x] + Pi/(x Sin[Pi x])");
    sf_zero("Gamma[1/2 + x] Gamma[1/2 - x] Cos[Pi x] - Pi");
}

/* ---- Gamma: multiplication ------------------------------------------------- */

static void test_gamma_multiplication(void) {
    sf_family("Gamma");
    static const char* Z[] = { "x", "x + 1", "a", "x + 1/3", "a + b" };
    for (size_t i = 0; i < sizeof Z / sizeof Z[0]; i++)
        sf_zero(FMT("Gamma[2 (%s)] - 2^(2 (%s) - 1)/Sqrt[Pi] Gamma[%s] Gamma[%s + 1/2]", Z[i], Z[i], Z[i], Z[i]));
    sf_to("Gamma[2 x]/(Gamma[x] Gamma[x + 1/2])", "2^(-1 + 2 x)/Sqrt[Pi]");   /* 2 2^-1 absorbed (O4) */
    sf_to("Gamma[3 x]/(Gamma[x] Gamma[x + 1/3] Gamma[x + 2/3])", "3^(-1/2 + 3 x)/(2 Pi)");
    sf_zero("Gamma[3 x] - (2 Pi)^(-1) 3^(3 x - 1/2) Gamma[x] Gamma[x + 1/3] Gamma[x + 2/3]");
    sf_zero("Gamma[4 x] - (2 Pi)^(-3/2) 4^(4 x - 1/2) Gamma[x] Gamma[x + 1/4] Gamma[x + 1/2] Gamma[x + 3/4]");
    sf_identity("Gamma[2 x] - 2^(2 x - 1)/Sqrt[Pi] Gamma[x] Gamma[x + 1/2]");
}

/* ---- Pochhammer ---------------------------------------------------------- */

static const char* P_ARGS[] = { "a", "a + 1/2", "2 a", "a + b", "x", "x/3" };
#define N_P (sizeof P_ARGS / sizeof P_ARGS[0])

static void test_pochhammer(void) {
    sf_family("Pochhammer");
    for (size_t i = 0; i < N_P; i++)
        for (int k = 1; k <= 4; k++) {
            /* (a)_{n+k} = (a)_n (a+n)_k */
            sf_zero(FMT("Pochhammer[%s, n + %d] - Pochhammer[%s, n] Pochhammer[%s + n, %d]",
                        P_ARGS[i], k, P_ARGS[i], P_ARGS[i], k));
            /* (a+k)_n (a)_k = (a)_n (a+n)_k */
            sf_zero(FMT("Pochhammer[%s + %d, n] Pochhammer[%s, %d] - Pochhammer[%s, n] Pochhammer[%s + n, %d]",
                        P_ARGS[i], k, P_ARGS[i], k, P_ARGS[i], P_ARGS[i], k));
        }
    for (size_t i = 0; i < 3; i++)
        for (int k = 1; k <= 2; k++)
            sf_zero(FMT("Pochhammer[%s, n + %d]/Pochhammer[%s, n] - Pochhammer[%s + n, %d]",
                        P_ARGS[i], k, P_ARGS[i], P_ARGS[i], k));
    sf_to("Pochhammer[a, n + 1]/Pochhammer[a, n]", "a + n");
    sf_to("Pochhammer[a, n + 1] + Pochhammer[a, n]", "(1 + a + n) Pochhammer[a, n]");
    sf_zero("Pochhammer[a, n + 1] - (a + n) Pochhammer[a, n]");
    sf_zero("Pochhammer[a + 1, n] a - Pochhammer[a, n] (a + n)");
    sf_zero("Pochhammer[a, n] Gamma[a] - Gamma[a + n]");
    sf_zero("Pochhammer[1, n] - n!");
    sf_zero("Pochhammer[a, n + 2] - (a + n) (a + n + 1) Pochhammer[a, n]");
    sf_zero("Pochhammer[a, n] Pochhammer[a + n, m] - Pochhammer[a, m] Pochhammer[a + m, n]");
    sf_zero("Pochhammer[x, n] - Gamma[x + n]/Gamma[x]");
    sf_unchanged("Pochhammer[a, n]");
    sf_unchanged("Pochhammer[a, n] + Pochhammer[b, n + 1]");
    sf_identity("Pochhammer[a, n + 2] - (a + n) (a + n + 1) Pochhammer[a, n]");
}

/* ---- Beta ---------------------------------------------------------------- */

static void test_beta(void) {
    sf_family("Beta");
    static const char* PAIRS[][2] = { { "a", "b" }, { "x", "y" }, { "a + 1/2", "b" }, { "2 a", "b" } };
    for (size_t p = 0; p < 4; p++) {
        const char* A = PAIRS[p][0];
        const char* B = PAIRS[p][1];
        for (int k = 1; k <= 4; k++) {
            sf_zero(FMT("Beta[%s + %d, %s] - Pochhammer[%s, %d]/Pochhammer[%s + %s, %d] Beta[%s, %s]",
                        A, k, B, A, k, A, B, k, A, B));
            sf_zero(FMT("Beta[%s, %s + %d] - Pochhammer[%s, %d]/Pochhammer[%s + %s, %d] Beta[%s, %s]",
                        A, B, k, B, k, A, B, k, A, B));
            sf_zero(FMT("Beta[%s + %d, %s + %d] - Pochhammer[%s, %d] Pochhammer[%s, %d]/Pochhammer[%s + %s, %d] Beta[%s, %s]",
                        A, k, B, k, A, k, B, k, A, B, 2 * k, A, B));
        }
    }
    sf_to("Beta[a + 1, b]/Beta[a, b]", "a/(a + b)");
    sf_zero("Beta[a, b + 1] + Beta[a + 1, b] - Beta[a, b]");
    sf_zero("Beta[a, b] Gamma[a + b] - Gamma[a] Gamma[b]");
    sf_zero("Beta[a, 1 - a] - Pi/Sin[Pi a]");
    sf_zero("Beta[a, b] - Beta[b, a]");
    sf_zero("(a + b) Beta[a + 1, b] - a Beta[a, b]");
    sf_zero("Beta[x, x] - 2^(1 - 2 x) Sqrt[Pi] Gamma[x]/Gamma[x + 1/2]");
    sf_zero("Beta[a, b] Beta[a + b, c] - Beta[b, c] Beta[a, b + c]");
    sf_zero("Beta[a, 1] - 1/a");
    sf_unchanged("Beta[a, b]");
    sf_unchanged("Beta[a, b] + Beta[c, b]");
    sf_identity("Beta[a, b + 1] + Beta[a + 1, b] - Beta[a, b]");
}

/* ---- Binomial ------------------------------------------------------------ */

static void test_binomial(void) {
    sf_family("Binomial");
    static const char* PAIRS[][2] = { { "n", "k" }, { "m", "k" }, { "n + 1/2", "k" }, { "2 n", "k" }, { "n", "2 k" } };
    for (size_t p = 0; p < 5; p++) {
        const char* N = PAIRS[p][0];
        const char* K = PAIRS[p][1];
        for (int j = 1; j <= 3; j++) {
            /* C(N, K+j) = C(N, K) (N-K-j+1)_j / (K+1)_j */
            sf_zero(FMT("Binomial[%s, %s + %d] - Binomial[%s, %s] Pochhammer[%s - (%s) - %d + 1, %d]/Pochhammer[%s + 1, %d]",
                        N, K, j, N, K, N, K, j, j, K, j));
            /* C(N+j, K) = C(N, K) (N+1)_j / (N-K+1)_j */
            sf_zero(FMT("Binomial[%s + %d, %s] - Binomial[%s, %s] Pochhammer[%s + 1, %d]/Pochhammer[%s - (%s) + 1, %d]",
                        N, j, K, N, K, N, j, N, K, j));
        }
        sf_zero(FMT("Binomial[%s + 1, %s] - Binomial[%s, %s] - Binomial[%s, %s - 1]", N, K, N, K, N, K));
        sf_zero(FMT("Binomial[%s, %s] - Binomial[%s, %s - (%s)]", N, K, N, N, K));
        sf_zero(FMT("(%s) Binomial[%s, %s] - (%s) Binomial[%s - 1, %s - 1]", K, N, K, N, N, K));
        sf_zero(FMT("(%s - (%s)) Binomial[%s, %s] - (%s) Binomial[%s - 1, %s]", N, K, N, K, N, N, K));
    }
    sf_zero("Binomial[n + 2, k] - Binomial[n, k] - 2 Binomial[n, k - 1] - Binomial[n, k - 2]");
    sf_zero("Binomial[n, k] Binomial[k, m] - Binomial[n, m] Binomial[n - m, k - m]");
    sf_zero("Binomial[n, k] - Gamma[n + 1]/(Gamma[k + 1] Gamma[n - k + 1])");
    sf_to("Binomial[n, k + 1]/Binomial[n, k]", "(-k + n)/(1 + k)");
    sf_zero("Binomial[n + 1, k + 1] - (n + 1)/(k + 1) Binomial[n, k]");
    sf_zero("Binomial[n, k] - Binomial[n - 1, k] - Binomial[n - 1, k - 1]");
    sf_unchanged("Binomial[n, k]");
    sf_unchanged("Binomial[n, k] + Binomial[m, j]");
    sf_identity("Binomial[n + 1, k] - Binomial[n, k] - Binomial[n, k - 1]");
}

/* ---- Factorial ----------------------------------------------------------- */

static void test_factorial(void) {
    sf_family("Factorial");
    static const char* N[] = { "n", "2 n", "m", "a + b", "n + 1/2" };
    for (size_t i = 0; i < 5; i++) {
        for (int k = 1; k <= 4; k++)
            sf_zero(FMT("(%s + %d)! - Pochhammer[%s + 1, %d] (%s)!", N[i], k, N[i], k, N[i]));
        for (int k = 1; k <= 3; k++)
            sf_zero(FMT("(%s + %d)!/(%s + %d)! - Pochhammer[%s + %d, 2]", N[i], k + 2, N[i], k, N[i], k + 1));
        sf_zero(FMT("(%s)! - Gamma[%s + 1]", N[i], N[i]));
        sf_zero(FMT("(%s + 1)!/Gamma[%s + 1] - (%s + 1)", N[i], N[i], N[i]));
    }
    sf_to("(n + 1)!/n!", "1 + n");
    sf_zero("(n + 1)! + n! - (n + 2) n!");
    sf_zero("(n + 2)! - (n + 2) (n + 1) n!");
    sf_zero("n! (n + 1) - (n + 1)!");
    sf_zero("(2 n)! - 2^(2 n) n! Gamma[n + 1/2]/Sqrt[Pi]");
    sf_zero("(n - 1)! n - n!");
    sf_unchanged("n!");
    sf_unchanged("n! + m!");
    sf_identity("(2 n)! - 2^(2 n) n! Gamma[n + 1/2]/Sqrt[Pi]");
}

/* ---- LogGamma -------------------------------------------------------------- */

static void test_loggamma(void) {
    sf_family("LogGamma");
    static const char* Z[] = { "x", "x + 1/3", "a + b", "2 x", "x + 2", "a + 1/2" };
    for (size_t i = 0; i < 6; i++)
        for (int k = 1; k <= 3; k++) {
            /* LogGamma[z + k] - LogGamma[z] - Sum_{j<k} Log[z + j] */
            char sum[512] = "";
            for (int j = 0; j < k; j++) {
                char t[64];
                snprintf(t, sizeof t, "%sLog[%s + %d]", j ? " + " : "", Z[i], j);
                strcat(sum, t);
            }
            sf_zero(FMT("LogGamma[%s + %d] - LogGamma[%s] - (%s)", Z[i], k, Z[i], sum));
            /* down: LogGamma[z] - LogGamma[z - k] - Sum_{j=1..k} Log[z - j] */
            char dsum[512] = "";
            for (int j = 1; j <= k; j++) {
                char t[64];
                snprintf(t, sizeof t, "%sLog[%s - %d]", j > 1 ? " + " : "", Z[i], j);
                strcat(dsum, t);
            }
            sf_zero(FMT("LogGamma[%s] - LogGamma[%s - %d] - (%s)", Z[i], Z[i], k, dsum));
        }
    for (size_t i = 0; i < 6; i++) {
        sf_zero(FMT("LogGamma[%s + 2] - LogGamma[%s + 1] - Log[%s + 1]", Z[i], Z[i], Z[i]));
        sf_zero(FMT("LogGamma[%s + 3] - LogGamma[%s + 1] - Log[%s + 1] - Log[%s + 2]", Z[i], Z[i], Z[i], Z[i]));
    }
    sf_zero("LogGamma[x + 1] - LogGamma[x] - Log[x]");
    sf_to("LogGamma[x + 3] - LogGamma[x]", "Log[x] + Log[1 + x] + Log[2 + x]");
    sf_unchanged("LogGamma[x] + LogGamma[y]");
    sf_unchanged("LogGamma[x] + LogGamma[x + 1/2]");
    sf_identity("LogGamma[x + 1] - LogGamma[x] - Log[x]");
}

/* ---- BarnesG ---------------------------------------------------------------- */

static void test_barnesg(void) {
    sf_family("BarnesG");
    static const char* Z[] = { "x", "x + 1/3", "a", "2 x", "a + b", "x/2" };
    for (size_t i = 0; i < 6; i++) {
        for (int k = 1; k <= 4; k++) {
            /* G(z + k) = G(z) Prod_{j<k} Gamma(z + j) */
            char prod[512] = "";
            for (int j = 0; j < k; j++) {
                char t[64];
                snprintf(t, sizeof t, "%sGamma[%s + %d]", j ? " " : "", Z[i], j);
                strcat(prod, t);
            }
            sf_zero(FMT("BarnesG[%s + %d] - %s BarnesG[%s]", Z[i], k, prod, Z[i]));
        }
        for (int k = 1; k <= 3; k++) {
            /* G(z) = G(z - k) Prod_{j=1..k} Gamma(z - j) */
            char prod[512] = "";
            for (int j = 1; j <= k; j++) {
                char t[64];
                snprintf(t, sizeof t, "%sGamma[%s - %d]", j > 1 ? " " : "", Z[i], j);
                strcat(prod, t);
            }
            sf_zero(FMT("BarnesG[%s] - %s BarnesG[%s - %d]", Z[i], prod, Z[i], k));
        }
    }
    for (size_t i = 0; i < 6; i++)
        sf_zero(FMT("BarnesG[%s + 2] - (%s) Gamma[%s]^2 BarnesG[%s]", Z[i], Z[i], Z[i], Z[i]));
    sf_to("BarnesG[x + 2]/BarnesG[x]", "x Gamma[x]^2");
    sf_unchanged("BarnesG[x] + BarnesG[y]");
    sf_identity("BarnesG[x + 1] - Gamma[x] BarnesG[x]");
}

/* ---- incomplete Gamma ------------------------------------------------------ */

static void test_gamma_incomplete(void) {
    sf_family("GammaIncomplete");
    static const char* A[] = { "a", "a + 1/2", "2 a", "b" };
    static const char* Zs[] = { "z", "2 z", "w" };
    for (size_t i = 0; i < 4; i++)
        for (size_t zi = 0; zi < 3; zi++)
            for (int k = 1; k <= 4; k++) {
                /* Gamma[A+k, Z] = (A)_k Gamma[A, Z] + Sum_{j<k} (A+j+1)_{k-j-1} Z^(A+j) E^-Z.
                 * E^-Z is distributed over the sum: a single E^-Z (polynomial)
                 * factor sends Simplify's bottom-up pass into minutes of
                 * TrigReduce on that special-function-free subterm before the
                 * engine sees the whole (pre-existing; tracked in
                 * tasks/todo.md). */
                char sum[1024] = "";
                for (int j = 0; j < k; j++) {
                    char t[200];
                    snprintf(t, sizeof t, " - Pochhammer[%s + %d, %d] (%s)^(%s + %d) E^-(%s)",
                             A[i], j + 1, k - j - 1, Zs[zi], A[i], j, Zs[zi]);
                    strcat(sum, t);
                }
                sf_zero(FMT("Gamma[%s + %d, %s] - Pochhammer[%s, %d] Gamma[%s, %s]%s",
                            A[i], k, Zs[zi], A[i], k, A[i], Zs[zi], sum));
            }
    sf_zero("Gamma[a + 1, z] - a Gamma[a, z] - z^a E^-z");
    sf_to("Gamma[a + 2, z] - (a + 1) Gamma[a + 1, z]", "z^(1 + a)/E^z");
    sf_zero("Gamma[a, z] - (Gamma[a + 1, z] - z^a E^-z)/a");
    sf_unchanged("Gamma[a, z] + Gamma[a, w]");
    sf_identity("Gamma[a + 1, z] - a Gamma[a, z] - z^a E^-z");
}

/* ---- soundness ---------------------------------------------------------------- */

static void test_unchanged(void) {
    sf_family("Gamma");
    sf_unchanged("Gamma[x] Gamma[y]");
    sf_unchanged("Gamma[x] + Gamma[x + 1/3]");                /* classes differ */
    sf_unchanged("Gamma[x] Gamma[1/2 - x]");                 /* not a reflection pair */
    sf_unchanged("Gamma[x]");
    sf_sound("Gamma[x] + Gamma[2 x]");
    sf_sound("Gamma[x + 1] Gamma[y] + Gamma[x] Gamma[y + 1]");
    sf_sound("Beta[a, b] + Binomial[n, k]");
    sf_sound("Gamma[x] Gamma[1 - x] + Gamma[x + 1/3] Gamma[2/3 - x]");
}

static void test_fast(void) {
    sf_family("Gamma");
    sf_fast("Gamma[x + 40]/Gamma[x]", 5.0);
    sf_fast("Gamma[x + 500]/Gamma[x + 498]", 5.0);
    sf_fast("Sum[Gamma[x + k], {k, 0, 10}]", 5.0);
    sf_fast("Binomial[n + 3, k + 2] - Binomial[n, k]", 5.0);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_gamma_recurrence_generated);
    TEST(test_gamma_recurrence_forms);
    TEST(test_gamma_reflection);
    TEST(test_gamma_multiplication);
    TEST(test_pochhammer);
    TEST(test_beta);
    TEST(test_binomial);
    TEST(test_factorial);
    TEST(test_loggamma);
    TEST(test_barnesg);
    TEST(test_gamma_incomplete);
    TEST(test_unchanged);
    TEST(test_fast);
    return sf_finish("simp_sf_gamma", 1500);
}
