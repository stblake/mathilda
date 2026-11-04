/* test_matfunc.c -- MatrixExp, MatrixLog, MatrixFunction and the non-integer
 * exponents of MatrixPower (src/linalg/matfunc*.c).
 *
 * Covers every precision class (exact, symbolic, machine real and complex,
 * MPFR real and complex), the List / packed / NDArray / SparseArray surfaces
 * (bit-identical machine answers), the 1 x 1 scalar case, the shape and
 * singularity messages, Compile[] lowering, the invariants that pin a
 * non-unique-looking answer (exp(log m) == m, exp(a + b) == exp(a) exp(b) for
 * commuting a and b, MatrixFunction[Exp, m] == MatrixExp[m], (m^(1/2))^2 == m),
 * and accuracy against the reference's 20-digit answers on hard cases (the
 * Moler--Van Loan examples, a large norm, a nonnormal / near-defective matrix,
 * a badly scaled one).
 *
 * Failures are counted rather than asserted: the CMake Release build passes
 * -DNDEBUG, which would no-op assert(). */

#include "test_utils.h"
#include "symtab.h"
#include "core.h"

static int failures = 0;

static void chk(const char* input, const char* expected) {
    struct Expr* p = parse_expression(input);
    if (!p) { fprintf(stderr, "FAIL(parse): %s\n", input); failures++; return; }
    struct Expr* e = evaluate(p);
    expr_free(p);
    char* s = expr_to_string(e);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  Expected: %s\n  Actual:   %s\n", input, expected, s);
        failures++;
    }
    free(s);
    expr_free(e);
}

static void chk_true(const char* input) { chk(input, "True"); }

/* The InputForm of an exact / symbolic answer (the reference's printed form). */
static void chk_if(const char* input, const char* expected) {
    char buf[4096];
    snprintf(buf, sizeof buf, "ToString[InputForm[%s]]", input);
    char exp[4096];
    snprintf(exp, sizeof exp, "\"%s\"", expected);
    chk(buf, exp);
}

/* Quiet[...] wrapper that also reports whether a message fired. */
static void chk_msg(const char* input, const char* expected_unevaluated) {
    char buf[4096];
    snprintf(buf, sizeof buf, "Check[%s, \"msg\"]", input);
    chk(buf, "\"msg\"");
    snprintf(buf, sizeof buf, "Quiet[%s]", input);
    chk(buf, expected_unevaluated);
}

static void setup(void) {
    /* relative normwise error, and an "is a numeric matrix" test */
    chk("relerr[a_, r_] := Max[Abs[Flatten[N[a] - N[r]]]] / Max[Abs[Flatten[N[r]]]]; Null", "Null");
    chk("m2 = {{1., 2.}, {3., 4.}}; c2 = {{1. + I, 2.}, {3., 4.}}; Null", "Null");
}

/* ---- MatrixExp --------------------------------------------------------- */

static void test_exp_exact(void) {
    chk_if("MatrixExp[{{1, 1}, {0, 1}}]", "{{E, E}, {0, E}}");
    chk_if("MatrixExp[{{0, 1}, {0, 0}}]", "{{1, 1}, {0, 1}}");
    chk_if("MatrixExp[{{0, 1}, {-1, 0}}]", "{{Cos[1], Sin[1]}, {-Sin[1], Cos[1]}}");
    chk_if("MatrixExp[{{1, -2}, {2, 1}}]", "{{E*Cos[2], -(E*Sin[2])}, {E*Sin[2], E*Cos[2]}}");
    chk_if("MatrixExp[{{a, 0}, {0, b}}]", "{{E^a, 0}, {0, E^b}}");
    chk_if("MatrixExp[{{1, 0, 0}, {1, 1, 0}, {0, 1, 1}}]", "{{E, 0, 0}, {E, E, 0}, {E/2, E, E}}");
    chk_if("MatrixExp[{{2, 1, 0}, {0, 2, 0}, {0, 0, 3}}]", "{{E^2, E^2, 0}, {0, E^2, 0}, {0, 0, E^3}}");
    chk_if("MatrixExp[{{2}}]", "{{E^2}}");
    chk_if("MatrixExp[{{0, 1}, {-1, 0}}, Parameter -> 2]", "{{Cos[2], Sin[2]}, {-Sin[2], Cos[2]}}");
    /* the closed form agrees with the reference's numeric value */
    chk_true("relerr[MatrixExp[{{1, 2}, {3, 4}}], {{51.96895619870503, 74.73656456700326}, "
             "{112.10484685050487, 164.0738030492099}}] < 1*^-14");
    chk_true("relerr[MatrixExp[{{a, b}, {c, d}}] /. {a -> 1, b -> 2, c -> 3, d -> 4}, "
             "MatrixExp[{{1., 2.}, {3., 4.}}]] < 1*^-13");
    chk_true("relerr[MatrixExp[{{1, 2}, {3, 4}}, {1, 1}], {126.70552076570829, 276.1786498997148}] < 1*^-14");
    /* an integer NDArray takes the exact path */
    chk_if("MatrixExp[NDArray[{{0, 1}, {0, 0}}, DataType -> \"int64\"]]", "{{1, 1}, {0, 1}}");
}

static void test_exp_machine(void) {
    chk_true("relerr[MatrixExp[m2], {{51.96895619870503, 74.73656456700326}, "
             "{112.10484685050487, 164.0738030492099}}] < 1*^-14");
    chk_true("relerr[MatrixExp[c2], {{44.04632415126094 + 24.33979511886785 I, "
             "69.06122774908384 + 23.418865710876865 I}, {103.59184162362575 + 35.1282985663153 I, "
             "159.3475986303251 + 24.93747981064122 I}}] < 1*^-14");
    /* 1 x 1: the scalar exp itself */
    chk("MatrixExp[{{2.}}] === {{Exp[2.]}}", "True");
    chk("MatrixExp[{{1. + 2. I}}] == {{Exp[1. + 2. I]}}", "True");
    /* real in, real out (packed) */
    chk("Head[MatrixExp[m2][[1, 1]]]", "Real");
    chk("NDArrayQ[MatrixExp[m2]]", "True");
    /* List == packed == NDArray, bit for bit */
    chk_true("MatrixExp[m2] === MatrixExp[ToPackedArray[m2]]");
    chk_true("Normal[MatrixExp[NDArray[m2]]] === MatrixExp[m2]");
    chk("Head[MatrixExp[NDArray[m2]]]", "NDArray");
    chk_true("MatrixExp[c2] === MatrixExp[ToPackedArray[c2]]");
    /* exp(0) = I, exp of a diagonal */
    chk("MatrixExp[{{0., 0.}, {0., 0.}}]", "{{1.0, 0.0}, {0.0, 1.0}}");
    /* MatrixExp[m, v] (the action) */
    chk_true("relerr[MatrixExp[m2, {1., 1.}], {126.70552076570829, 276.1786498997148}] < 1*^-14");
    chk_true("relerr[MatrixExp[m2, {1., 1.}], MatrixExp[m2] . {1., 1.}] < 1*^-14");
    chk_true("relerr[MatrixExp[c2, {1., 2. I}], MatrixExp[c2] . {1., 2. I}] < 1*^-14");
    chk_true("MatrixExp[m2, {1., 1.}] === Normal[MatrixExp[NDArray[m2], NDArray[{1., 1.}]]]");
    /* a large matrix: the Taylor action stays accurate */
    chk("SeedRandom[7]; big = RandomReal[{-1, 1}, {60, 60}] / 10; bv = RandomReal[{-1, 1}, 60]; Null", "Null");
    chk_true("relerr[MatrixExp[big, bv], MatrixExp[big] . bv] < 1*^-13");
}

static void test_exp_mpfr(void) {
    chk_true("relerr[MatrixExp[N[{{1, 2}, {3, 4}}, 30]], N[MatrixExp[{{1, 2}, {3, 4}}], 40]] < 1*^-28");
    chk("Precision[MatrixExp[N[{{1, 2}, {3, 4}}, 30]][[1, 1]]] >= 29", "True");
    chk_true("relerr[MatrixExp[N[{{1 + I, 2}, {3, 4}}, 30]], N[MatrixExp[{{1 + I, 2}, {3, 4}}], 40]] < 1*^-27");
    chk_true("relerr[MatrixExp[N[{{0, 1}, {-1, 0}}, 50]], N[{{Cos[1], Sin[1]}, {-Sin[1], Cos[1]}}, 60]] < 1*^-48");
    chk_true("relerr[MatrixExp[N[{{1, 2}, {3, 4}}, 30], N[{1, 1}, 30]], N[MatrixExp[{{1, 2}, {3, 4}}, {1, 1}], 40]] < 1*^-28");
}

/* The Moler--Van Loan examples and other hard cases, against the reference's
 * 20-digit answers. */
static void test_exp_hard(void) {
    chk_true("relerr[MatrixExp[{{-49., 24.}, {-64., 31.}}], {{-0.73575875814475307964, 0.55181909965809770062}, "
             "{-1.47151759908826053498, 1.10363824071557258908}}] < 1*^-12");
    chk_true("relerr[MatrixExp[N[{{-131, 19, 18}, {-390, 56, 54}, {-387, 57, 52}}]], "
             "{{-1.50964415879608968340, 0.36787943911028869916, 0.13533528117545906946}, "
             "{-5.63257079990259601498, 1.47151775850230841907, 0.40600584352637720837}, "
             "{-4.93493832609810712588, 1.10363831733086609747, 0.54134112676298990026}}] < 1*^-12");
    /* nonnormal, near-defective */
    chk_true("relerr[MatrixExp[{{1., 10.^4}, {0., 1. + 10.^-8}}], {{2.71828182845904523536, 27182.8184205045442296}, "
             "{0, 2.71828185564186365586}}] < 1*^-13");
    /* large norm */
    chk_true("relerr[MatrixExp[{{100., 15.}, {-20., 50.}}], {{3.0075875843774318083*^40, 1.0484816677413626239*^40}, "
             "{-1.397975556988483499*^40, -4.8735130809377693791*^39}}] < 1*^-12");
    /* badly scaled */
    chk_true("relerr[MatrixExp[{{1., 10.^-8}, {10.^8, 1.}}], {{4.19452804946532511362, 3.1945280494653251136*^-8}, "
             "{3.1945280494653251136*^8, 4.19452804946532511362}}] < 1*^-12");
    /* nilpotent: exact in floating point */
    chk("MatrixExp[{{0., 6., 0., 0.}, {0., 0., 6., 0.}, {0., 0., 0., 6.}, {0., 0., 0., 0.}}]",
        "{{1.0, 6.0, 18.0, 36.0}, {0.0, 1.0, 6.0, 18.0}, {0.0, 0.0, 1.0, 6.0}, {0.0, 0.0, 0.0, 1.0}}");
}

static void test_exp_messages(void) {
    chk_msg("MatrixExp[{{1, 2, 3}, {3, 4, 5}}]", "MatrixExp[{{1, 2, 3}, {3, 4, 5}}]");
    chk_msg("MatrixExp[{}]", "MatrixExp[{}]");
    chk_msg("MatrixExp[{1, 2}]", "MatrixExp[{1, 2}]");
    chk("MatrixExp[x]", "MatrixExp[x]");
    chk_msg("MatrixExp[m2, {1., 1., 1.}]", "MatrixExp[{{1.0, 2.0}, {3.0, 4.0}}, {1.0, 1.0, 1.0}]");
    chk_msg("MatrixExp[m2, {{1., 0.}, {0., 1.}}]", "MatrixExp[{{1.0, 2.0}, {3.0, 4.0}}, {{1.0, 0.0}, {0.0, 1.0}}]");
    chk("Attributes[MatrixExp]", "{Protected}");
}

/* ---- MatrixLog --------------------------------------------------------- */

static void test_log(void) {
    chk_if("MatrixLog[{{2, 1}, {1, 2}}]", "{{Log[3]/2, Log[3]/2}, {Log[3]/2, Log[3]/2}}");
    chk_if("MatrixLog[{{1, 1}, {0, 1}}]", "{{0, 1}, {0, 0}}");
    chk_if("MatrixLog[{{4, 0}, {0, 9}}]", "{{Log[4], 0}, {0, Log[9]}}");
    chk_if("MatrixLog[{{-1, 0}, {0, 2}}]", "{{I*Pi, 0}, {0, Log[2]}}");
    chk_if("MatrixLog[{{a, 0}, {0, b}}]", "{{Log[a], 0}, {0, Log[b]}}");
    chk_msg("MatrixLog[{{0, 1}, {0, 0}}]", "MatrixLog[{{0, 1}, {0, 0}}]");
    chk_msg("MatrixLog[{{1, 0}, {0, 0}}]", "MatrixLog[{{1, 0}, {0, 0}}]");
    chk_msg("MatrixLog[{{1., 0.}, {0., 0.}}]", "MatrixLog[{{1.0, 0.0}, {0.0, 0.0}}]");
    chk_msg("MatrixLog[{{1, 2, 3}, {3, 4, 5}}]", "MatrixLog[{{1, 2, 3}, {3, 4, 5}}]");
    /* machine */
    chk_true("relerr[MatrixLog[{{2., 1.}, {1., 2.}}], {{0.5493061443340548, 0.5493061443340548}, "
             "{0.5493061443340548, 0.5493061443340548}}] < 1*^-14");
    chk_true("relerr[MatrixLog[m2], {{-0.3504398139985535 + 2.391117954451724 I, 0.9293512057047016 - 1.0937621702091036 I}, "
             "{1.3940268085570522 - 1.6406432553136552 I, 1.0435869945584988 + 0.7504746991380686 I}}] < 1*^-14");
    chk("MatrixLog[{{-2.}}] == {{Log[-2.]}}", "True");
    chk("MatrixLog[{{1., 1.}, {0., 1.}}]", "{{0.0, 1.0}, {0.0, 0.0}}");
    /* a real matrix with a complex-conjugate pair has a REAL logarithm */
    chk("Head[MatrixLog[{{0., 1.}, {-1., 0.}}][[1, 2]]]", "Real");
    chk_true("relerr[MatrixLog[{{0., 1.}, {-1., 0.}}], {{0, Pi/2}, {-Pi/2, 0}}] < 1*^-14");
    /* a negative real eigenvalue gives the principal (complex) logarithm */
    chk_true("relerr[MatrixLog[{{-1., 0.}, {0., 2.}}], {{I Pi, 0}, {0, Log[2]}}] < 1*^-15");
    chk_true("relerr[MatrixLog[{{-49., 24.}, {-64., 31.}}], {{8.4996400321686482407 + 3.1415926535897932385 I, "
             "-4.2498200160843241204}, {11.332853376224864321, -5.6664266881124321605 + 3.1415926535897932385 I}}] < 1*^-12");
    chk_true("relerr[MatrixLog[N[{{4, 2, 0}, {1, 4, 1}, {0, 1, 4}}]], {{1.3170812395271424516, 0.53531766275749593227, "
             "-0.069213121592748167205}, {0.26765883137874796613, 1.2824746787307683680, 0.26765883137874796613}, "
             "{-0.034606560796374083603, 0.26765883137874796613, 1.3516878003235165352}}] < 1*^-13");
    chk_true("relerr[MatrixLog[{{1., 1000.}, {0., 2.}}], {{0, 693.14718055994530942}, {0, 0.69314718055994530942}}] < 1*^-14");
    /* the invariant: exp(log m) == m */
    chk_true("relerr[MatrixExp[MatrixLog[m2]], m2] < 1*^-13");
    chk_true("relerr[MatrixExp[MatrixLog[c2]], c2] < 1*^-13");
    chk_true("relerr[MatrixExp[MatrixLog[N[{{4, 2, 0}, {1, 4, 1}, {0, 1, 4}}]]], {{4, 2, 0}, {1, 4, 1}, {0, 1, 4}}] < 1*^-14");
    /* surfaces */
    chk_true("MatrixLog[m2] === MatrixLog[ToPackedArray[m2]]");
    chk_true("Normal[MatrixLog[NDArray[m2]]] === MatrixLog[m2]");
    /* MPFR keeps the precision */
    chk_true("relerr[MatrixLog[N[{{2, 1}, {1, 2}}, 25]], N[{{Log[3]/2, Log[3]/2}, {Log[3]/2, Log[3]/2}}, 40]] < 1*^-23");
    chk("Precision[MatrixLog[N[{{2, 1}, {1, 2}}, 25]][[1, 1]]] >= 24", "True");
    chk_true("relerr[MatrixLog[N[{{-1, 0}, {0, 2}}, 30]], N[{{I Pi, 0}, {0, Log[2]}}, 40]] < 1*^-28");
    chk_true("relerr[MatrixLog[N[{{1, 2}, {3, 4}}, 30]], MatrixLog[m2]] < 1*^-13");
    chk_true("relerr[MatrixExp[MatrixLog[N[{{1, 2}, {3, 4}}, 30]]], {{1, 2}, {3, 4}}] < 1*^-27");
    chk("Attributes[MatrixLog]", "{Protected}");
}

/* ---- MatrixFunction ---------------------------------------------------- */

static void test_function(void) {
    chk_if("MatrixFunction[Exp, {{1, 1}, {0, 1}}]", "{{E, E}, {0, E}}");
    chk_if("MatrixFunction[Sin, {{1, 1}, {0, 1}}]", "{{Sin[1], Cos[1]}, {0, Sin[1]}}");
    chk_if("MatrixFunction[Sqrt, {{4, 1}, {0, 4}}]", "{{2, 1/4}, {0, 2}}");
    chk_if("MatrixFunction[Sqrt, {{2, 1}, {1, 2}}]", "{{1/2 + Sqrt[3]/2, -1/2 + Sqrt[3]/2}, {-1/2 + Sqrt[3]/2, 1/2 + Sqrt[3]/2}}");
    chk_if("MatrixFunction[f, {{1, 1}, {0, 1}}]", "{{f[1], Derivative[1][f][1]}, {0, f[1]}}");
    chk_if("MatrixFunction[f, {{1, 0}, {0, 2}}]", "{{f[1], 0}, {0, f[2]}}");
    chk_if("MatrixFunction[f, {{1, 1, 0}, {0, 1, 1}, {0, 0, 1}}]", "{{f[1], Derivative[1][f][1], Derivative[2][f][1]/2}, {0, f[1], Derivative[1][f][1]}, {0, 0, f[1]}}");
    chk_if("MatrixFunction[Abs, {{-1, 0}, {0, 2}}]", "{{1, 0}, {0, 2}}");
    chk_true("Simplify[MatrixFunction[#^2 + 1 &, {{1, 2}, {3, 4}}] - ({{1, 2}, {3, 4}} . {{1, 2}, {3, 4}} + IdentityMatrix[2])] "
             "== {{0, 0}, {0, 0}}");
    /* machine */
    chk_true("relerr[MatrixFunction[Sin, m2], {{-0.4655814863137303, -0.14842445991317654}, "
             "{-0.22263668986976481, -0.6882181761834951}}] < 1*^-14");
    chk("Head[MatrixFunction[Sin, m2][[1, 1]]]", "Real");
    chk_true("relerr[MatrixFunction[Sin, {{1., 1.}, {0., 1.}}], {{Sin[1], Cos[1]}, {0, Sin[1]}}] < 1*^-15");
    chk_true("relerr[MatrixFunction[Sin[#] + Cos[2 #] &, {{1., 0.5}, {0.25, 2.}}], "
             "{{0.5556395238751893, -0.075194305636404}, {-0.03759715281820197, 0.4052509126023812}}] < 1*^-14");
    chk_true("relerr[MatrixFunction[Sin, {{1 + I, 2}, {3, 4.}}], {{-0.4864223289528402 + 0.6172649752412648 I, "
             "-0.116737828952654 - 0.2404478960490102 I}, {-0.17510674342898092 - 0.36067184407351544 I, "
             "-0.7817530204063263 + 0.31496204564407626 I}}] < 1*^-14");
    chk_true("relerr[MatrixFunction[Sin, N[{{1, 2, 0}, {0, 1, 3}, {0, 0, 2}}]], {{0.84147098480789650665, "
             "1.0806046117362794348, -2.8348551831021271719}, {0, 0.84147098480789650665, 0.20347932605335556623}, "
             "{0, 0, 0.90929742682568169540}}] < 1*^-14");
    chk_true("relerr[MatrixFunction[Cos, {{1., 1.}, {0., 1. + 10.^-6}}], {{0.54030230586813971740, "
             "-0.84147125495890919554}, {0, 0.54030146439688475849}}] < 1*^-12");
    /* an undefined f on a machine matrix: the eigen-decomposition form */
    chk_true("relerr[MatrixFunction[g, m2] /. g -> Exp, MatrixExp[m2]] < 1*^-13");
    /* the invariant MatrixFunction[Exp, m] == MatrixExp[m] (the same kernel) */
    chk_true("MatrixFunction[Exp, m2] === MatrixExp[m2]");
    chk_true("relerr[MatrixFunction[Exp[#] &, m2], MatrixExp[m2]] < 1*^-13");
    chk_true("relerr[MatrixFunction[Exp, c2], MatrixExp[c2]] < 1*^-15");
    /* surfaces */
    chk_true("MatrixFunction[Sin, m2] === MatrixFunction[Sin, ToPackedArray[m2]]");
    chk_true("Normal[MatrixFunction[Sin, NDArray[m2]]] === MatrixFunction[Sin, m2]");
    /* MPFR */
    chk_true("relerr[MatrixFunction[Sin, N[{{1, 2}, {3, 4}}, 25]], {{-0.46558148631373039292721425973, "
             "-0.14842445991317655963529454645}, {-0.22263668986976483945294181968, -0.68821817618349523238015607941}}] < 1*^-23");
    /* messages */
    chk_msg("MatrixFunction[Abs, {{1., 1.}, {0., 1.}}]", "MatrixFunction[Abs, {{1.0, 1.0}, {0.0, 1.0}}]");
    chk_msg("MatrixFunction[Log, {{1, 0}, {0, 0}}]", "MatrixFunction[Log, {{1, 0}, {0, 0}}]");
    chk_msg("MatrixFunction[Exp]", "MatrixFunction[Exp]");
    chk_msg("MatrixFunction[Exp, {{1, 2, 3}, {3, 4, 5}}]", "MatrixFunction[Exp, {{1, 2, 3}, {3, 4, 5}}]");
    chk("MatrixFunction[Exp, x]", "MatrixFunction[Exp, x]");
    chk("Attributes[MatrixFunction]", "{Protected}");
}

/* ---- MatrixPower, non-integer exponents -------------------------------- */

static void test_power(void) {
    chk_if("MatrixPower[{{4, 0}, {0, 9}}, 1/2]", "{{2, 0}, {0, 3}}");
    chk_if("MatrixPower[{{2, 1}, {1, 2}}, 1/2]", "{{1/2 + Sqrt[3]/2, -1/2 + Sqrt[3]/2}, {-1/2 + Sqrt[3]/2, 1/2 + Sqrt[3]/2}}");
    chk_if("MatrixPower[{{4, 1}, {0, 4}}, 1/2]", "{{2, 1/4}, {0, 2}}");
    chk_if("MatrixPower[{{2, 1}, {1, 2}}, 1/3]", "{{1/2 + 3^(1/3)/2, -1/2 + 3^(1/3)/2}, {-1/2 + 3^(1/3)/2, 1/2 + 3^(1/3)/2}}");
    chk_if("MatrixPower[{{2, 1}, {1, 2}}, s]", "{{1/2 + 3^s/2, -1/2 + 3^s/2}, {-1/2 + 3^s/2, 1/2 + 3^s/2}}");
    chk_if("MatrixPower[{{2, 1}, {0, 2}}, s]", "{{2^s, 2^(-1 + s)*s}, {0, 2^s}}");
    chk_if("MatrixPower[{{a, 0}, {0, b}}, s]", "{{a^s, 0}, {0, b^s}}");
    chk_if("MatrixPower[{{1, 0}, {0, 0}}, 1/2]", "{{1, 0}, {0, 0}}");
    chk_if("MatrixPower[{{2, 1}, {1, 2}}, 1/2, {1, 0}]", "{1/2 + Sqrt[3]/2, -1/2 + Sqrt[3]/2}");
    chk_if("MatrixPower[{{0, 1}, {0, 0}}, 1/2]", "MatrixPower[{{0, 1}, {0, 0}}, 1/2]");
    /* the integer path is untouched */
    chk_if("MatrixPower[{{1, 2}, {3, 4}}, 2]", "{{7, 10}, {15, 22}}");
    chk("MatrixPower[{{1., 2.}, {3., 4.}}, 2]", "{{7.0, 10.0}, {15.0, 22.0}}");
    /* machine */
    chk_true("relerr[MatrixPower[{{2., 1.}, {1., 2.}}, 0.5], {{1.3660254037844386, 0.3660254037844386}, "
             "{0.3660254037844386, 1.3660254037844386}}] < 1*^-15");
    chk_true("MatrixPower[{{2., 1.}, {1., 2.}}, 0.5] === MatrixPower[{{2., 1.}, {1., 2.}}, 1/2]");
    chk_true("relerr[MatrixPower[m2, 0.5], {{0.5536885671459112 + 0.46439416283907065 I, 0.8069607270132164 - 0.21242647876640203 I}, "
             "{1.2104410905198244 - 0.318639718149603 I, 1.7641296576657357 + 0.14575444468946763 I}}] < 1*^-14");
    chk_true("relerr[MatrixPower[m2, 2.5], {{15.980230875219622 + 0.06436195837746439 I, 23.290021665749872 - 0.029440904470137808 I}, "
             "{34.93503249862481 - 0.0441613567052066 I, 50.91526337384444 + 0.020200601672257346 I}}] < 1*^-13");
    chk_true("relerr[MatrixPower[m2, -0.5], {{0.10306395622800202 - 1.247428043827744 I, 0.15020820363930287 + 0.5706074022222716 I}, "
             "{0.22531230545895453 + 0.8559111033334074 I, 0.32837626168695677 - 0.39151694049433683 I}}] < 1*^-13");
    chk_true("relerr[MatrixPower[{{2., 1.}, {1., 2.}}, 0.5 + I], {{0.8938964326326693 + 0.7712623421114616 I, "
             "-0.10610356736733079 + 0.7712623421114616 I}, {-0.10610356736733079 + 0.7712623421114616 I, "
             "0.8938964326326693 + 0.7712623421114616 I}}] < 1*^-14");
    chk_true("relerr[MatrixPower[c2, 0.5], {{0.9603587567972733 + 0.7165310911013546 I, 0.6562116352644154 - 0.19112422247830346 I}, "
             "{0.9843174528966231 - 0.2866863337174551 I, 1.849114098454745 + 0.10173893975169165 I}}] < 1*^-14");
    chk_true("relerr[MatrixPower[N[{{4, 2, 0}, {1, 4, 1}, {0, 1, 4}}], 1/3], {{1.5636349452490148706, 0.27471014120740597439, "
             "-0.023766106719184604155}, {0.13735507060370298719, 1.5517518918894225685, 0.13735507060370298719}, "
             "{-0.011883053359592302078, 0.13735507060370298719, 1.5755179986086071727}}] < 1*^-14");
    chk("MatrixPower[{{4.}}, 0.5]", "{{2.0}}");
    chk("MatrixPower[{{1., 1.}, {0., 1.}}, 0.5]", "{{1.0, 0.5}, {0.0, 1.0}}");
    chk("Head[MatrixPower[{{0., 1.}, {-1., 0.}}, 0.5][[1, 1]]]", "Real");
    chk_true("relerr[MatrixPower[{{2., 1.}, {1., 2.}}, 0.5, {1., 0.}], {1.3660254037844386, 0.3660254037844386}] < 1*^-15");
    chk("MatrixPower[{{1., 0.}, {0., 0.}}, 0.5]", "{{1.0, 0.0}, {0.0, 0.0}}");
    chk_msg("MatrixPower[{{0., 1.}, {0., 0.}}, 0.5]", "MatrixPower[{{0.0, 1.0}, {0.0, 0.0}}, 0.5]");
    chk_msg("MatrixPower[{{1., 0.}, {0., 0.}}, -0.5]", "MatrixPower[{{1.0, 0.0}, {0.0, 0.0}}, -0.5]");
    /* an exact matrix with a machine exponent answers in machine numbers */
    chk_true("relerr[MatrixPower[{{2, 1}, {1, 2}}, 0.5], {{1.3660254037844386, 0.3660254037844386}, "
             "{0.3660254037844386, 1.3660254037844386}}] < 1*^-15");
    /* a symbolic exponent on a machine matrix */
    chk_true("relerr[MatrixPower[{{2., 1.}, {1., 2.}}, x] /. x -> 0.5, MatrixPower[{{2., 1.}, {1., 2.}}, 0.5]] < 1*^-14");
    /* the invariant (m^(1/2))^2 == m */
    chk_true("relerr[MatrixPower[m2, 0.5] . MatrixPower[m2, 0.5], m2] < 1*^-14");
    chk_true("relerr[MatrixPower[c2, 1/2] . MatrixPower[c2, 1/2], c2] < 1*^-14");
    chk_true("relerr[MatrixPower[m2, 1/3] . MatrixPower[m2, 1/3] . MatrixPower[m2, 1/3], m2] < 1*^-13");
    chk_true("Simplify[MatrixPower[{{2, 1}, {1, 2}}, 1/2] . MatrixPower[{{2, 1}, {1, 2}}, 1/2]] == {{2, 1}, {1, 2}}");
    /* surfaces */
    chk_true("MatrixPower[m2, 0.5] === MatrixPower[ToPackedArray[m2], 0.5]");
    chk_true("Normal[MatrixPower[NDArray[m2], 0.5]] === MatrixPower[m2, 0.5]");
    /* MPFR */
    chk_true("relerr[MatrixPower[N[{{2, 1}, {1, 2}}, 25], 1/2], N[MatrixPower[{{2, 1}, {1, 2}}, 1/2], 40]] < 1*^-23");
    chk("Precision[MatrixPower[N[{{2, 1}, {1, 2}}, 25], 1/2][[1, 1]]] >= 24", "True");
    chk_true("relerr[MatrixPower[N[{{1, 2}, {3, 4}}, 30], 1/3], MatrixPower[m2, 1/3]] < 1*^-13");
    chk_true("relerr[MatrixPower[N[{{1, 2}, {3, 4}}, 30], 1/2] . MatrixPower[N[{{1, 2}, {3, 4}}, 30], 1/2], {{1, 2}, {3, 4}}] < 1*^-27");
    chk_true("relerr[MatrixPower[N[{{1, 2}, {3, 4}}, 30], N[1/3, 30]], MatrixPower[m2, 1/3]] < 1*^-13");
}

/* ---- commuting exponentials, SparseArray, Compile ----------------------- */

static void test_invariants(void) {
    chk("ca = {{1., 2.}, {2., -1.}}; cb = 2.5 ca + 0.75 IdentityMatrix[2]; Null", "Null");
    chk_true("relerr[MatrixExp[ca + cb], MatrixExp[ca] . MatrixExp[cb]] < 1*^-13");
    chk_true("relerr[MatrixExp[2 big], MatrixExp[big] . MatrixExp[big]] < 1*^-13");
    chk_true("relerr[MatrixExp[-big] . MatrixExp[big], IdentityMatrix[60]] < 1*^-13");
}

static void test_sparse(void) {
    chk("MatrixExp[SparseArray[{{1, 2}, {3, 4}}]] == MatrixExp[{{1, 2}, {3, 4}}]", "True");
    chk_true("relerr[Normal[MatrixExp[SparseArray[m2]]], MatrixExp[m2]] < 1*^-15");
    chk_true("relerr[Normal[MatrixLog[SparseArray[{{2., 1.}, {1., 2.}}]]], MatrixLog[{{2., 1.}, {1., 2.}}]] < 1*^-15");
    chk_true("relerr[Normal[MatrixFunction[Sin, SparseArray[m2]]], MatrixFunction[Sin, m2]] < 1*^-15");
    chk_true("relerr[Normal[MatrixPower[SparseArray[{{2., 1.}, {1., 2.}}], 0.5]], MatrixPower[{{2., 1.}, {1., 2.}}, 0.5]] < 1*^-15");
}

static void test_compile(void) {
    chk("\"Compiled\" /. CompileDiagnostics[{{m, _Real, 2}}, MatrixExp[m]]", "True");
    chk("\"Compiled\" /. CompileDiagnostics[{{m, _Complex, 2}}, MatrixExp[m]]", "True");
    chk("\"Compiled\" /. CompileDiagnostics[{{m, _Real, 2}}, MatrixLog[m]]", "True");
    chk("\"Compiled\" /. CompileDiagnostics[{{m, _Complex, 2}}, MatrixLog[m]]", "True");
    chk("\"Compiled\" /. CompileDiagnostics[{{m, _Real, 2}, {s, _Real}}, MatrixPower[m, s]]", "True");
    chk("\"Compiled\" /. CompileDiagnostics[{{m, _Real, 2}}, MatrixPower[m, 0.5]]", "True");
    chk("\"Compiled\" /. CompileDiagnostics[{{m, _Real, 2}, {v, _Real, 1}}, MatrixExp[m, v]]", "True");
    chk("\"Compiled\" /. CompileDiagnostics[{{m, _Complex, 2}, {v, _Real, 1}}, MatrixExp[m, v]]", "True");
    /* compiled == interpreted, bit for bit */
    chk("cfe = Compile[{{m, _Real, 2}}, MatrixExp[m]]; cfl = Compile[{{m, _Real, 2}}, MatrixLog[m]]; "
        "cfp = Compile[{{m, _Real, 2}, {s, _Real}}, MatrixPower[m, s]]; "
        "cfv = Compile[{{m, _Real, 2}, {v, _Real, 1}}, MatrixExp[m, v]]; Null", "Null");
    chk_true("Normal[cfe[m2]] === MatrixExp[m2]");
    chk_true("Normal[cfl[{{2., 1.}, {1., 2.}}]] === MatrixLog[{{2., 1.}, {1., 2.}}]");
    chk_true("Normal[cfp[{{2., 1.}, {1., 2.}}, 0.5]] === MatrixPower[{{2., 1.}, {1., 2.}}, 0.5]");
    chk_true("Normal[cfv[m2, {1., 1.}]] === MatrixExp[m2, {1., 1.}]");
    chk_true("Normal[cfe[big]] === MatrixExp[big]");
    /* a real matrix whose log is complex: the compiled call falls back and
     * still gives the interpreter's answer */
    chk_true("Normal[Quiet[cfl[m2]]] === MatrixLog[m2]");
    chk_true("Normal[Quiet[cfp[m2, 0.5]]] === MatrixPower[m2, 0.5]");
    /* 1 x 1: the scalar functions */
    chk_true("Normal[cfe[{{2.}}]] === {{Exp[2.]}}");
}

int main(void) {
    symtab_init();
    core_init();
    setup();
    test_exp_exact();
    test_exp_machine();
    test_exp_mpfr();
    test_exp_hard();
    test_exp_messages();
    test_log();
    test_function();
    test_power();
    test_invariants();
    test_sparse();
    test_compile();
    if (failures) {
        fprintf(stderr, "\n%d matrix-function test(s) FAILED\n", failures);
        return 1;
    }
    printf("All matrix-function tests passed.\n");
    return 0;
}
