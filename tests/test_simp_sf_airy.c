#include "simp_sf_test_util.h"
#include "symtab.h"
#include "core.h"

/*
 * Simplify, special-function engine, Airy / Struve / Whittaker / Fibonacci
 * module (src/simp/simp_sf_airy.c).  Every identity is derived
 * (docs/design/simp_special_functions.md, section 9) and was checked
 * numerically at non-integer and complex parameters before it was encoded:
 *
 *   Ai Bi' - Ai' Bi = 1/Pi                     (Abel: w'' = z w; value at 0)
 *   H[v-1] + H[v+1] = (2v/z) H[v] + (z/2)^v/(Sqrt[Pi] Gamma[v + 3/2])
 *   z^2 H'' + z H' + (z^2 - v^2) H = 4 (z/2)^(v+1)/(Sqrt[Pi] Gamma[v + 1/2])
 *   W[k+1] + (2k - z) W[k] - (m - k + 1/2)(m + k - 1/2) W[k-1] = 0
 *   (2m + 2k + 1) M[k+1] + 2 (z - 2k) M[k] - (2m - 2k + 1) M[k-1] = 0
 *   F[n+2] = F[n+1] + F[n];  F[n+2, x] = x F[n+1, x] + F[n, x]  (Binet, all n)
 *     => F[n+3] = 2 F[n+1] + F[n],  F[n+4] = 3 F[n+1] + 2 F[n]   (likewise L)
 *   F[n+1] L[n] - F[n] L[n+1] = 2 Cos[Pi n]        (Binet, all n)
 */

static char g_buf[8192];
#define FMT(...) (snprintf(g_buf, sizeof g_buf, __VA_ARGS__), g_buf)

static const char* Z[] = { "x", "2 x", "x + 1/3", "a + b", "x/2", "-x", "x^2", "Sqrt[x]" };
static const char* N[] = { "n", "n + 1", "1/3", "4/3", "2/5", "n/2", "2 n", "a + b" };

/* ---- Airy ---------------------------------------------------------------- */

static void test_airy(void) {
    sf_family("Airy");
    for (size_t i = 0; i < 8; i++) {
        const char* z = Z[i];
        char W[512];
        snprintf(W, sizeof W, "(AiryAi[%s] AiryBiPrime[%s] - AiryAiPrime[%s] AiryBi[%s])", z, z, z, z);
        sf_zero(FMT("%s - 1/Pi", W));
        sf_zero(FMT("x (%s - 1/Pi)", W));
        sf_zero(FMT("(%s - 1/Pi) AiryAi[%s]", W, z));
        sf_zero(FMT("(%s - 1/Pi)^2", W));
        sf_zero(FMT("Pi %s - 1", W));
        sf_zero(FMT("%s AiryAi[%s] - AiryAi[%s]/Pi", W, z, z));
    }
    sf_zero("AiryAi[x] D[AiryBi[x], x] - D[AiryAi[x], x] AiryBi[x] - 1/Pi");
    sf_zero("D[AiryAi[x], {x, 2}] - x AiryAi[x]");
    sf_to("AiryAi[z] AiryBiPrime[z] - AiryAiPrime[z] AiryBi[z]", "1/Pi");
    sf_unchanged("AiryAi[x] AiryBiPrime[y] - AiryAiPrime[x] AiryBi[y]");
}

/* ---- StruveH --------------------------------------------------------------- */

static void test_struve(void) {
    sf_family("StruveH");
    static const char* ZS[] = { "z", "2 z", "w" };
    for (size_t i = 0; i < 8; i++)
        for (size_t j = 0; j < 3; j++)
            sf_zero(FMT("StruveH[(%s) - 1, %s] + StruveH[(%s) + 1, %s] - 2 (%s)/(%s) StruveH[%s, %s] - ((%s)/2)^(%s)/(Sqrt[Pi] Gamma[(%s) + 3/2])",
                        N[i], ZS[j], N[i], ZS[j], N[i], ZS[j], N[i], ZS[j], ZS[j], N[i], N[i]));
    for (size_t i = 0; i < 8; i++) {
        const char* n = N[i];
        sf_zero(FMT("z^2 D[StruveH[%s, z], {z, 2}] + z D[StruveH[%s, z], z] + (z^2 - (%s)^2) StruveH[%s, z] - 4 (z/2)^((%s) + 1)/(Sqrt[Pi] Gamma[(%s) + 1/2])",
                    n, n, n, n, n, n));
        sf_zero(FMT("x (StruveH[(%s) - 1, z] + StruveH[(%s) + 1, z] - 2 (%s)/z StruveH[%s, z] - (z/2)^(%s)/(Sqrt[Pi] Gamma[(%s) + 3/2]))",
                    n, n, n, n, n, n));
        sf_zero(FMT("(StruveH[(%s) - 1, z] + StruveH[(%s) + 1, z] - 2 (%s)/z StruveH[%s, z] - (z/2)^(%s)/(Sqrt[Pi] Gamma[(%s) + 3/2])) StruveH[%s, z]",
                    n, n, n, n, n, n, n));
    }
    sf_to("StruveH[n - 1, z] + StruveH[n + 1, z] - 2 n/z StruveH[n, z]", "z^n/(2^n Sqrt[Pi] Gamma[3/2 + n])");
    sf_unchanged("StruveH[n, z] + StruveH[n, w]");
}

/* ---- Whittaker ------------------------------------------------------------- */

static const char* K[] = { "k", "k + 1/2", "2 k", "k/3", "a + b", "1/3" };
static const char* M[] = { "m", "m + 1/3" };
static const char* ZW[] = { "z", "2 z" };

static void test_whittakerW(void) {
    sf_family("WhittakerW");
    for (size_t i = 0; i < 6; i++)
        for (size_t j = 0; j < 2; j++)
            for (size_t l = 0; l < 2; l++) {
                const char *k = K[i], *m = M[j], *z = ZW[l];
                char R[1024];
                snprintf(R, sizeof R,
                         "(WhittakerW[(%s) + 1, %s, %s] + (2 (%s) - (%s)) WhittakerW[%s, %s, %s] - ((%s) - (%s) + 1/2) ((%s) + (%s) - 1/2) WhittakerW[(%s) - 1, %s, %s])",
                         k, m, z, k, z, k, m, z, m, k, m, k, k, m, z);
                sf_zero(R);
                if (l == 0) sf_zero(FMT("x %s", R));
                if (l == 1) sf_zero(FMT("%s WhittakerW[%s, %s, %s]", R, k, m, z));
            }
    sf_zero("WhittakerW[k + 1, m, z] + (2 k - z) WhittakerW[k, m, z] - (m - k + 1/2) (m + k - 1/2) WhittakerW[k - 1, m, z]");
    sf_unchanged("WhittakerW[k, m, z] + WhittakerW[k, m, w]");
}

static void test_whittakerM(void) {
    sf_family("WhittakerM");
    for (size_t i = 0; i < 6; i++)
        for (size_t j = 0; j < 2; j++)
            for (size_t l = 0; l < 2; l++) {
                const char *k = K[i], *m = M[j], *z = ZW[l];
                char R[1024];
                snprintf(R, sizeof R,
                         "((2 (%s) + 2 (%s) + 1) WhittakerM[(%s) + 1, %s, %s] + 2 ((%s) - 2 (%s)) WhittakerM[%s, %s, %s] - (2 (%s) - 2 (%s) + 1) WhittakerM[(%s) - 1, %s, %s])",
                         m, k, k, m, z, z, k, k, m, z, m, k, k, m, z);
                sf_zero(R);
                if (l == 0) sf_zero(FMT("x %s", R));
                if (l == 1) sf_zero(FMT("%s WhittakerM[%s, %s, %s]", R, k, m, z));
            }
    sf_zero("(2 m + 2 k + 1) WhittakerM[k + 1, m, z] + 2 (z - 2 k) WhittakerM[k, m, z] - (2 m - 2 k + 1) WhittakerM[k - 1, m, z]");
    sf_unchanged("WhittakerM[k, m, z] + WhittakerM[k, m, w]");
}

/* ---- Fibonacci / LucasL ------------------------------------------------- */

static void test_fib_lucas(const char* H, bool pair) {
    sf_family(H);
    for (size_t i = 0; i < 8; i++) {
        const char* n = N[i];
        sf_zero(FMT("%s[(%s) + 2] - %s[(%s) + 1] - %s[%s]", H, n, H, n, H, n));
        sf_zero(FMT("%s[(%s) + 3] - 2 %s[(%s) + 1] - %s[%s]", H, n, H, n, H, n));
        sf_zero(FMT("%s[(%s) + 4] - 3 %s[(%s) + 1] - 2 %s[%s]", H, n, H, n, H, n));
        sf_zero(FMT("%s[(%s) + 2, x] - x %s[(%s) + 1, x] - %s[%s, x]", H, n, H, n, H, n));
        sf_zero(FMT("%s[(%s) + 2, 2 x] - 2 x %s[(%s) + 1, 2 x] - %s[%s, 2 x]", H, n, H, n, H, n));
        if (pair) sf_zero(FMT("Fibonacci[(%s) + 1] LucasL[%s] - Fibonacci[%s] LucasL[(%s) + 1] - 2 Cos[Pi (%s)]", n, n, n, n, n));
        else      sf_zero(FMT("x (%s[(%s) + 2] - %s[(%s) + 1] - %s[%s])", H, n, H, n, H, n));
    }
    sf_zero(FMT("%s[n + 5] - 5 %s[n + 1] - 3 %s[n]", H, H, H));
    sf_unchanged(FMT("%s[n] + %s[m]", H, H));
}

static void test_fibonacci(void) {
    test_fib_lucas("Fibonacci", true);
    sf_to("Fibonacci[n + 1] LucasL[n] - Fibonacci[n] LucasL[n + 1]", "2 Cos[n Pi]");
}
static void test_lucas(void) { test_fib_lucas("LucasL", false); }

int main(void) {
    symtab_init();
    core_init();
    TEST(test_airy);
    TEST(test_struve);
    TEST(test_whittakerW);
    TEST(test_whittakerM);
    TEST(test_fibonacci);
    TEST(test_lucas);
    return sf_finish("simp_sf_airy", 1000);
}
