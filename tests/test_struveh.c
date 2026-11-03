/* Tests for StruveH[nu, z], the Struve function H_nu(z).
 *
 * Covers: the exact origin value (StruveH[nu,0]=0 for real nu > -1) and
 * symbolic/exact passthrough; machine-precision real numerics across the
 * convergent-series (small/moderate z) and Y_nu+K_nu asymptotic (large z)
 * regimes; a closed-form oracle for half-integer order (H_{1/2}, H_{-1/2});
 * the intrinsically ill-conditioned transition band (looser tolerance, the
 * double floor); arbitrary-precision real values (vs mpmath, exact arguments)
 * and precision tracking; a genuinely-complex arbitrary-precision value; list
 * threading; argument-count diagnostics; attributes.
 *
 * Reference values are from mpmath 1.4 (mp.struveh). StruveH has no symbolic
 * .m rewrites, so -- unlike Bessel -- every path here is the C kernel the REPL
 * also uses; this harness loads core_init() only (no init.m), which is exactly
 * what the numeric paths need.
 */

#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---- numeric helpers ------------------------------------------------- */

static double eval_real(const char* input) {
    Expr* e = parse_expression(input);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    expr_free(e);
    ASSERT_MSG(r->type == EXPR_REAL, "%s: expected a Real result", input);
    double v = r->data.real;
    expr_free(r);
    return v;
}

static void assert_close(const char* input, double expected, double tol) {
    double v = eval_real(input);
    ASSERT_MSG(fabs(v - expected) <= tol,
               "%s: expected %.12g, got %.12g", input, expected, v);
}

/* ---- exact origin value and passthrough ----------------------------- */

void test_struveh_exact() {
    /* H_nu(0) = 0 for a real nu > -1 (exact). */
    assert_eval_eq("StruveH[0, 0]", "0", 0);
    assert_eval_eq("StruveH[1, 0]", "0", 0);
    assert_eval_eq("StruveH[3, 0]", "0", 0);
    assert_eval_eq("StruveH[1/2, 0]", "0", 0);
    assert_eval_eq("StruveH[-1/2, 0]", "0", 0);
    /* Exact non-zero / symbolic arguments stay symbolic (no auto-numericize). */
    assert_eval_eq("StruveH[0, 2]", "StruveH[0, 2]", 0);
    assert_eval_eq("StruveH[nu, x]", "StruveH[nu, x]", 0);
    assert_eval_eq("StruveH[1, z]", "StruveH[1, z]", 0);
}

/* ---- machine-precision real (vs mpmath) ----------------------------- */

void test_struveh_machine_real() {
    assert_close("StruveH[0, 0.5]", 0.30955591458375472,  1e-9);
    assert_close("StruveH[0, 1.0]", 0.56865662704828795,  1e-9);
    assert_close("StruveH[0, 1.5]", 0.73672346560439987,  1e-9);
    assert_close("StruveH[1, 1.0]", 0.19845733620194440,  1e-9);
    assert_close("StruveH[1, 1.3]", 0.32012314597642024,  1e-9);
    assert_close("StruveH[2, 3.0]", 0.74238666967748319,  1e-9);
    assert_close("StruveH[2, 5.5]", 1.62394289163940280,  1e-9);
}

/* ---- large-z asymptotic (Y_nu + K_nu) ------------------------------- */

void test_struveh_asymptotic() {
    assert_close("StruveH[0, 30.0]", -0.09609842155416211, 1e-9);
    assert_close("StruveH[1, 30.0]",  0.72175037834698999, 1e-9);
    assert_close("StruveH[0, 100.0]", -0.070878751689647343, 1e-9);
}

/* ---- transition band (the intrinsic double floor, ~1e-9 abs) -------- */
/* z ~ [17, 24] (small nu) is where both elementary series are
 * ill-conditioned in double; the kernel returns the best available value and
 * the surfaces still agree bit-for-bit. A looser tolerance documents the floor
 * (scipy's struve shares this). Arbitrary precision is exact there. */
void test_struveh_transition_band() {
    assert_close("StruveH[0, 16.0]", 0.13544931808186468, 1e-7);
    assert_close("StruveH[2, 17.0]", 3.75765108226997820, 1e-6);
}

/* ---- closed-form oracle for half-integer order ---------------------- */
/* H_{1/2}(x) = sqrt(2/(pi x)) (1 - cos x),  H_{-1/2}(x) = sqrt(2/(pi x)) sin x
 * (DLMF 11.4.5-6). For a half-integer order the K_nu asymptotic terminates, so
 * both regimes are near-exact -- the strongest oracle-free check of the core. */
void test_struveh_halfinteger_oracle() {
    for (int t = 1; t <= 30; t++) {
        double x = t * 0.5;
        char in[64];
        double s = sqrt(2.0 / (M_PI * x));

        snprintf(in, sizeof(in), "StruveH[1/2, %.4f]", x);
        double got = eval_real(in);
        double want = s * (1.0 - cos(x));
        ASSERT_MSG(fabs(got - want) <= 1e-9 * (1.0 + fabs(want)),
                   "StruveH[1/2, %.4f]: builtin %.14g vs closed form %.14g",
                   x, got, want);

        snprintf(in, sizeof(in), "StruveH[-1/2, %.4f]", x);
        got = eval_real(in);
        want = s * sin(x);
        ASSERT_MSG(fabs(got - want) <= 1e-9 * (1.0 + fabs(want)),
                   "StruveH[-1/2, %.4f]: builtin %.14g vs closed form %.14g",
                   x, got, want);
    }
}

/* ---- nu = -(n+1/2): the Gamma(nu+3/2) poles ------------------------- */
/* H_{-(n+1/2)}(x) = (-1)^n J_{n+1/2}(x) (DLMF 11.4.4): the 1F2 series loses its
 * first n terms to 1/Gamma(k+1-n) = 0 and K_nu vanishes identically. The machine
 * kernel used to decline at the pole and the unregularized fallback template
 * returned `0.0 HypergeometricPFQ[...]`. Oracle: the elementary J_{3/2}, J_{5/2};
 * spot values from mpmath 1.4 (mp.struveh). */
void test_struveh_negative_halfinteger() {
    for (int t = 1; t <= 40; t++) {
        double x = t * 0.5;
        char in[64];
        double s = sqrt(2.0 / (M_PI * x)), sn = sin(x), cs = cos(x);
        double j32 = s * (sn / x - cs);
        double j52 = s * ((3.0 / (x * x) - 1.0) * sn - 3.0 * cs / x);

        snprintf(in, sizeof(in), "StruveH[-3/2, %.4f]", x);
        double got = eval_real(in);
        ASSERT_MSG(fabs(got + j32) <= 1e-13 * (1.0 + fabs(j32)),
                   "StruveH[-3/2, %.4f]: %.16g vs -J_{3/2} %.16g", x, got, -j32);

        snprintf(in, sizeof(in), "StruveH[-5/2, %.4f]", x);
        got = eval_real(in);
        ASSERT_MSG(fabs(got - j52) <= 1e-13 * (1.0 + fabs(j52)),
                   "StruveH[-5/2, %.4f]: %.16g vs J_{5/2} %.16g", x, got, j52);
    }
    assert_close("StruveH[-3/2, 2.]", -0.491293778687162, 1e-14);
    assert_close("StruveH[-5/2, 2.]", 0.223924531468916, 1e-14);
    assert_close("StruveH[-7/2, 2.]", -0.0685175499851271, 1e-15);
    assert_close("StruveH[-3/2, 30.]", 0.0272679457111777, 1e-15);
    /* the Hankel series of a half-integer order is finite: summed whole, not
     * cut at its first growing term (that left an 8e-10 floor here) */
    assert_close("StruveH[-13/2, 20.]", -0.147368651190097, 1e-14);
    assert_close("StruveH[9/2, 20.]", 75.9102805710853, 1e-12);
    /* kernel declines (complex value) -> arb at 53 bits, demoted to machine */
    assert_close("Re[StruveH[-3/2, 2. + I]]", -0.646752436117891, 1e-14);
    assert_close("Im[StruveH[-3/2, 2. + I]]", -0.172080645894644, 1e-14);
    assert_close("Im[StruveH[-3/2, -2.]]", 0.491293778687162, 1e-14);
    assert_eval_eq("FreeQ[StruveH[-3/2, -2.], HypergeometricPFQ]", "True", 0);
    /* scalar == packed == NDArray == Compile, bit for bit */
    assert_eval_eq(
        "With[{zs = 0.37 N[Range[40]]}, "
        "Table[StruveH[-7/2, z], {z, zs}] === StruveH[-7/2, zs] === "
        "Normal[StruveH[-7/2, NDArray[zs]]] === "
        "Compile[{{v, _Real, 1}}, StruveH[-7/2, v]][zs] === "
        "Table[Compile[{{x, _Real}}, StruveH[-7/2, x]][z], {z, zs}]]", "True", 0);
}

/* ---- negative-argument parity (integer nu is real) ------------------ */
/* H_nu(-z) = (-1)^(nu+1) H_nu(z); for integer nu this is a real value. */
void test_struveh_negative_arg() {
    /* H_1 is even: H_1(-z) = H_1(z). */
    assert_close("StruveH[1, -1.3]", 0.32012314597642024, 1e-9);
    /* H_0 is odd: H_0(-z) = -H_0(z). */
    assert_close("StruveH[0, -1.5]", -0.73672346560439987, 1e-9);
}

/* ---- arbitrary-precision real (vs mpmath, exact arguments) ---------- */

void test_struveh_arbitrary_real() {
    /* 50 digits from an exact argument; last printed digit is unpinned by N. */
    assert_eval_startswith("N[StruveH[0, 3], 50]",
        "0.5743061488143983979756175433066279369885366610067");
    assert_eval_startswith("N[StruveH[1, 13/10], 50]",
        "0.3201231459764202395239403961772816052031112705381");
    assert_eval_startswith("N[StruveH[2, 7/2], 50]",
        "1.005657904832457848395587132374085637990379463552");
    /* Precision tracks the input. */
    assert_eval_startswith("StruveH[0, 3.0`45]",
        "0.57430614881439839797561754330662793698853666");
}

/* ---- genuinely-complex arbitrary precision -------------------------- */
/* H_1(1 + 2 I) = -0.486907231297804... + 1.203779527377793... I (mpmath). */
void test_struveh_complex() {
    assert_eval_startswith("N[StruveH[1, 1 + 2 I], 30]",
        "-0.4869072312978040009409313970");
    assert_eval_startswith("N[Im[StruveH[1, 1 + 2 I]], 30]",
        "1.2037795273777930258452114499");
}

/* ---- list threading (Listable) -------------------------------------- */

void test_struveh_listable() {
    assert_eval_eq("StruveH[1, {a, b}]", "{StruveH[1, a], StruveH[1, b]}", 0);
    assert_eval_eq("StruveH[1, {}]", "{}", 0);
    assert_close("StruveH[0, {0.5, 1.0, 1.5}][[1]]", 0.30955591458375472, 1e-9);
    assert_close("StruveH[0, {0.5, 1.0, 1.5}][[2]]", 0.56865662704828795, 1e-9);
    assert_close("StruveH[0, {0.5, 1.0, 1.5}][[3]]", 0.73672346560439987, 1e-9);
    /* Threading also over the order. */
    assert_eval_eq("StruveH[{0, 1}, x]", "{StruveH[0, x], StruveH[1, x]}", 0);
}

/* ---- argument count ------------------------------------------------- */

void test_struveh_argcount() {
    assert_eval_eq("StruveH[]", "StruveH[]", 0);
    assert_eval_eq("StruveH[1]", "StruveH[1]", 0);
    assert_eval_eq("StruveH[1, 2, 3]", "StruveH[1, 2, 3]", 0);
}

/* ---- attributes ----------------------------------------------------- */

void test_struveh_attributes() {
    SymbolDef* d = symtab_get_def("StruveH");
    ASSERT_MSG(d != NULL, "StruveH not registered");
    ASSERT_MSG((d->attributes & ATTR_LISTABLE) != 0, "StruveH not Listable");
    ASSERT_MSG((d->attributes & ATTR_NUMERICFUNCTION) != 0, "StruveH not NumericFunction");
    ASSERT_MSG((d->attributes & ATTR_PROTECTED) != 0, "StruveH not Protected");
}

/* ---- half-odd-integer order: elementary closed forms ---------------- */

void test_struveh_halfinteger_elementary() {
    /* DLMF 11.4.4/11.4.5 with the spherical-Bessel forms written out. */
    assert_eval_eq("StruveH[1/2, x]",
        "Sqrt[2]/(Sqrt[Pi] Sqrt[x]) - (Sqrt[2] Cos[x])/(Sqrt[Pi] Sqrt[x])", 0);
    assert_eval_eq("StruveH[-1/2, x]", "(Sqrt[2] Sin[x])/(Sqrt[Pi] Sqrt[x])", 0);
    assert_eval_eq("StruveH[1/2, 2]", "1/Sqrt[Pi] - Cos[2]/Sqrt[Pi]", 0);
    /* Values at 13/10 against the reference system's 25-digit StruveH. */
    assert_close("N[StruveH[5/2, x] /. x -> 13/10]", 0.038804081097075112823, 1e-15);
    assert_close("N[StruveH[7/2, x] /. x -> 13/10]", 0.0063957243807799023608, 1e-16);
    assert_close("N[StruveH[9/2, x] /. x -> 13/10]", 0.00083931298990262896181, 1e-17);
    assert_close("N[StruveH[-9/2, x] /. x -> 13/10]", 0.0025450232156187135632, 1e-16);
    assert_close("N[StruveH[-3/2, x] /. x -> 13/10]", -0.33149086451171976044, 1e-15);
    /* Other rational orders stay symbolic. */
    assert_eval_eq("StruveH[1/3, x]", "StruveH[1/3, x]", 0);
}

/* ---- derivatives (DLMF 11.4.23-26) ----------------------------------- */

void test_struve_derivative() {
    assert_eval_eq("D[StruveH[nu, x], x]",
        "1/2 (x^nu/(2^nu Sqrt[Pi] Gamma[3/2 + nu]) + StruveH[-1 + nu, x] - StruveH[1 + nu, x])", 0);
    assert_eval_eq("D[StruveL[nu, x], x]",
        "1/2 (x^nu/(2^nu Sqrt[Pi] Gamma[3/2 + nu]) + StruveL[-1 + nu, x] + StruveL[1 + nu, x])", 0);
    assert_eval_eq("D[StruveH[nu, x], nu]", "Derivative[1, 0][StruveH][nu, x]", 0);
    assert_eval_eq("D[StruveH[nu, y], x]", "0", 0);
    /* The derivative of the elementary H_{1/2} agrees with the rule's value. */
    assert_eval_eq("Abs[N[(D[StruveH[1/2, x], x] - (1/2 (StruveH[-1/2, x] - "
                   "StruveH[3/2, x] + (2^(-1/2) x^(1/2))/(Sqrt[Pi] Gamma[2])))) "
                   "/. x -> 13/10]] < 10^-14", "True", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_struveh_exact);
    TEST(test_struveh_machine_real);
    TEST(test_struveh_asymptotic);
    TEST(test_struveh_transition_band);
    TEST(test_struveh_halfinteger_oracle);
    TEST(test_struveh_negative_halfinteger);
    TEST(test_struveh_negative_arg);
    TEST(test_struveh_arbitrary_real);
    TEST(test_struveh_complex);
    TEST(test_struveh_listable);
    TEST(test_struveh_argcount);
    TEST(test_struveh_attributes);
    TEST(test_struveh_halfinteger_elementary);
    TEST(test_struve_derivative);

    printf("All StruveH tests passed.\n");
    return 0;
}
