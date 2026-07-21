/* test_ndarray_functions.c — element-wise evaluation of elementary and special
 * functions over NDArray objects (see src/ndkernels.c, ndarray_map_unary /
 * ndarray_map_binary, and the eval.c dispatch hook).
 *
 * The core correctness property is "matches the List path": for machine data D,
 * f[NDArray[D]] must equal f[D] elementwise. We assert that via a tolerance
 * comparison in the language itself
 *     Max[Abs[Flatten[Normal[f[NDArray[D]]] - f[D]]]] < eps
 * which is robust to last-ULP differences between the C buffer loop and the
 * scalar builtin, and treats real and complex results uniformly (Abs is the
 * complex magnitude). Structural properties (result stays a packed NDArray,
 * dtype promotion on real->complex escape, faithful degrade at poles) are
 * asserted directly.
 *
 * Uses exit(1) on failure rather than assert(): the CMake Release build passes
 * -DNDEBUG, which no-ops assert(). */

#include "test_utils.h"
#include "symtab.h"
#include "core.h"

/* Evaluate `input`; fail unless it prints exactly `expected`. */
static void chk_eq(const char* input, const char* expected) {
    struct Expr* p = parse_expression(input);
    if (!p) { fprintf(stderr, "FAIL(parse): %s\n", input); exit(1); }
    struct Expr* e = evaluate(p);
    expr_free(p);
    char* s = expr_to_string(e);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  Expected: %s\n  Actual:   %s\n",
                input, expected, s);
        exit(1);
    }
    free(s);
    expr_free(e);
}

/* Assert f[NDArray[data]] matches f[data] elementwise within 1e-9. */
static void chk_matches(const char* fn, const char* data) {
    char buf[512];
    snprintf(buf, sizeof(buf),
             "Max[Abs[Flatten[Normal[%s[NDArray[%s]]] - %s[%s]]]] < 1/1000000000",
             fn, data, fn, data);
    chk_eq(buf, "True");
}

/* Like chk_matches but packs the NDArray at an explicit dtype (needed for
 * complex data, which the default float64 dtype cannot hold). */
static void chk_matches_dt(const char* fn, const char* data, const char* dtype) {
    /* 32-bit dtypes carry ~1e-7 precision, so loosen the bound for them. */
    const char* eps = strstr(dtype, "32") ? "1/100000" : "1/1000000000";
    char buf[640];
    snprintf(buf, sizeof(buf),
             "Max[Abs[Flatten[Normal[%s[NDArray[%s, DataType -> \"%s\"]]] - %s[%s]]]] < %s",
             fn, data, dtype, fn, data, eps);
    chk_eq(buf, "True");
}

/* Binary f[scalar, NDArray[data]] vs f[scalar, data]. */
static void chk_matches2(const char* fn, const char* scalar, const char* data) {
    char buf[640];
    snprintf(buf, sizeof(buf),
             "Max[Abs[Flatten[Normal[%s[%s, NDArray[%s]]] - %s[%s, %s]]]] < 1/1000000000",
             fn, scalar, data, fn, scalar, data);
    chk_eq(buf, "True");
}

static void test_elementary_real(void) {
    const char* d = "{0.3, 1.4, 2.7, -0.6}";
    const char* fns[] = {
        "Sin","Cos","Tan","Cot","Sec","Csc",
        "Sinh","Cosh","Tanh","Coth","Sech","Csch",
        "Exp","ArcTan","ArcCot","ArcSinh","ArcCsch",
        "Abs","Re","Im","Arg","Conjugate","Sign",
    };
    for (unsigned i = 0; i < sizeof(fns)/sizeof(fns[0]); i++)
        chk_matches(fns[i], d);
}

static void test_elementary_escape(void) {
    /* Functions that leave the real axis for some real inputs. */
    const char* d = "{0.3, 1.4, 2.7, -0.6}";
    const char* fns[] = {
        "Log","ArcSin","ArcCos","ArcSec","ArcCsc",
        "ArcCosh","ArcTanh","ArcCoth","ArcSech",
    };
    for (unsigned i = 0; i < sizeof(fns)/sizeof(fns[0]); i++)
        chk_matches(fns[i], d);
    /* Real->complex escape promotes the whole array to complex, matching WL. */
    chk_eq("DataType[Log[NDArray[{-1., 4.}]]]", "\"complex64\"");
    chk_eq("DataType[ArcSin[NDArray[{0.5, 2.}]]]", "\"complex64\"");
    /* A fully in-domain real input stays real (narrowed back). */
    chk_eq("DataType[Log[NDArray[{1., 4.}]]]", "\"float64\"");
}

static void test_structure_and_dtype(void) {
    /* Result stays a packed NDArray for the C-loop path. */
    chk_eq("Head[Sin[NDArray[{1., 2.}]]]", "NDArray");
    chk_eq("Head[Gamma[NDArray[{0.5, 1.}]]]", "NDArray");
    /* dtype (component width) is preserved. */
    chk_eq("DataType[Sin[NDArray[{1., 2.}, DataType -> \"float32\"]]]", "\"float32\"");
    /* Projections yield a real dtype even from complex input. */
    chk_eq("DataType[Abs[NDArray[{1.+2. I}, DataType -> \"complex64\"]]]", "\"float64\"");
    chk_eq("DataType[Re[NDArray[{1.+2. I}, DataType -> \"complex64\"]]]", "\"float64\"");
}

static void test_complex_and_float32_inputs(void) {
    chk_matches_dt("Sin", "{1.+2. I, 0.5-1. I, -0.3+0.7 I}", "complex64");
    chk_matches_dt("Exp", "{1.+2. I, 0.5-1. I}", "complex64");
    chk_matches_dt("Conjugate", "{1.+2. I, 3.-4. I}", "complex64");
    chk_matches_dt("Cos", "{1.+2. I}", "complex32");
    chk_matches_dt("Tan", "{0.3, 1.4, -0.6}", "float32");   /* float32 real path */
}

static void test_binary(void) {
    chk_matches2("Log", "10", "{1., 10., 100.}");     /* narrows to real */
    chk_matches2("Log", "2", "{-1., 4.}");            /* escapes to complex */
    chk_matches2("ArcTan", "1.", "{1., 0., -1.}");    /* two-arg ArcTan */
    /* array-first ordering: ArcTan[array, y]. */
    chk_eq("Max[Abs[Flatten[Normal[ArcTan[NDArray[{1., 0.}], 1.]] "
           "- ArcTan[{1., 0.}, 1.]]]] < 1/1000000000", "True");
}

static void test_special_cloop(void) {
    const char* d = "{0.5, 1.4, 2.6, 3.1}";
    chk_matches("Gamma", d);
    chk_matches("LogGamma", d);
    chk_matches("Erf", d);
    chk_matches("Erfc", d);
    chk_matches("Factorial", d);
    /* Factorial degrades at a negative-integer pole (non-finite). */
    chk_eq("Head[Factorial[NDArray[{1., 2.}]]]", "NDArray");
    chk_eq("Head[Factorial[NDArray[{-1., 2.}]]]", "List");
    chk_matches2("BesselJ", "2", "{1., 2., 3.}");
    chk_matches2("BesselY", "1", "{1., 2., 3.}");
    chk_matches2("Beta", "2.", "{1., 2., 3.}");
    /* Gamma stays packed off the poles... */
    chk_eq("Head[Gamma[NDArray[{0.5, 1.5}]]]", "NDArray");
    /* ...and degrades faithfully at a pole (ComplexInfinity), yielding a List. */
    chk_eq("Head[Gamma[NDArray[{0., 1.}]]]", "List");
    chk_eq("Gamma[NDArray[{0., 1., 2.}]]", "{ComplexInfinity, 1.0, 1.0}");
}

static void test_special_degrade(void) {
    /* Libc-free special functions still evaluate NDArray inputs correctly via
     * the List degrade path (result is a List, values match). */
    const char* d = "{0.4, 1.3, 2.6}";
    const char* uf[] = {
        "Erfi","ExpIntegralEi","SinIntegral","CosIntegral","FresnelS","FresnelC",
        "AiryAi","AiryBi","ProductLog","Zeta",
    };
    for (unsigned i = 0; i < sizeof(uf)/sizeof(uf[0]); i++)
        chk_matches(uf[i], d);
    chk_matches2("BesselI", "2", d);
    chk_matches2("BesselK", "2", d);
    chk_matches2("PolyLog", "2", d);
    chk_matches2("LegendreP", "2", d);
}

static void test_rounding(void) {
    /* Floor/Ceiling/Round/IntegerPart/FractionalPart over a real array match the
     * List path, incl. banker's rounding at the .5 boundary and floored sign. */
    const char* d = "{1.2, 2.5, 3.5, -1.7, -2.5}";
    const char* fns[] = {
        "Floor", "Ceiling", "Round", "IntegerPart", "FractionalPart",
    };
    for (unsigned i = 0; i < sizeof(fns)/sizeof(fns[0]); i++)
        chk_matches(fns[i], d);
    /* Round is half-to-even: 2.5->2, 3.5->4, -2.5->-2. */
    chk_eq("Round[NDArray[{2.5, 3.5, -2.5}]]", "NDArray[{2.0, 4.0, -2.0}]");
    /* Stays packed on real input; complex input degrades to the List path. */
    chk_eq("Head[Floor[NDArray[{1.5, 2.5}]]]", "NDArray");
    chk_eq("Head[Floor[NDArray[{1.5 + 0.5 I}, DataType -> \"complex64\"]]]", "List");
}

static void test_sqrt_and_rational_power(void) {
    /* Sqrt[NDArray] = NDArray^(1/2) now routes through the C power loop instead
     * of renaming back to Sqrt[NDArray]. */
    chk_eq("Head[Sqrt[NDArray[Range[4]]]]", "NDArray");
    chk_eq("Max[Abs[Flatten[Normal[Sqrt[NDArray[{1., 2., 4., 9.}]]] "
           "- Sqrt[{1., 2., 4., 9.}]]]] < 1/1000000000", "True");
    /* Rational exponent p/q, and real->complex escape (negative base). */
    chk_eq("Max[Abs[Flatten[Normal[NDArray[{1., 8., 27.}]^(1/3)] "
           "- {1., 8., 27.}^(1/3)]]] < 1/1000000000", "True");
    chk_eq("DataType[Sqrt[NDArray[{-1., 4.}]]]", "\"complex64\"");
}

static void test_rational_scalar_arithmetic(void) {
    /* A Rational scalar broadcasts over Plus/Times (1/2 NDArray = Times[1/2, ...])
     * instead of leaving the product symbolic. */
    chk_eq("Head[NDArray[Range[6]] / 2]", "NDArray");
    chk_eq("NDArray[{1., 2., 3.}] / 2", "NDArray[{0.5, 1.0, 1.5}]");
    chk_eq("3/4 * NDArray[{4., 8.}]", "NDArray[{3.0, 6.0}]");
    chk_eq("NDArray[{1., 2.}] + 1/2", "NDArray[{1.5, 2.5}]");
    /* A genuinely symbolic scalar still declines (stays unevaluated). */
    chk_eq("Head[x * NDArray[{1., 2.}]]", "Times");
}

static void test_mod_quotient(void) {
    /* Mod[array, n] / Quotient[array, n] match the List path (floored). */
    chk_eq("Mod[NDArray[{1., 2., 3., 4., 5.}], 3]", "NDArray[{1.0, 2.0, 0.0, 1.0, 2.0}]");
    chk_eq("Quotient[NDArray[{1., 2., 3., 4., 5.}], 3]", "NDArray[{0.0, 0.0, 1.0, 1.0, 1.0}]");
    /* Floored: negatives carry the divisor's sign. */
    chk_eq("Mod[NDArray[{-1., -2., -3.}], 3]", "NDArray[{2.0, 1.0, 0.0}]");
    /* Complex operand / zero divisor decline -> faithful List degrade. */
    chk_eq("Head[Mod[NDArray[{1. + 1. I}, DataType -> \"complex64\"], 3]]", "List");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_elementary_real);
    TEST(test_elementary_escape);
    TEST(test_structure_and_dtype);
    TEST(test_complex_and_float32_inputs);
    TEST(test_binary);
    TEST(test_rounding);
    TEST(test_sqrt_and_rational_power);
    TEST(test_rational_scalar_arithmetic);
    TEST(test_mod_quotient);
    TEST(test_special_cloop);
    TEST(test_special_degrade);

    printf("All NDArray function tests passed.\n");
    return 0;
}
