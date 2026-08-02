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
    /* Round is half-to-even: 2.5->2, 3.5->4, -2.5->-2 -- and the result is an
     * INT64 buffer of exact Integers, not float64.
     *
     * This expectation was NDArray[{2.0, 4.0, -2.0}] until 2026-07-30. That was
     * the old limitation showing through: with no narrowing kernel category the
     * only thing a real-closed kernel could do was keep the float64 dtype, so
     * the visible NDArray surface disagreed with BOTH the List path and
     * Mathematica, which give exact Integers -- while this very test's stated
     * purpose (line above) is that the two match. NDUnaryKernel.to_int now
     * writes NDT_INT64, and all three agree. */
    chk_eq("Round[NDArray[{2.5, 3.5, -2.5}]]", "NDArray[{2, 4, -2}]");
    chk_eq("DataType[Round[NDArray[{2.5, 3.5, -2.5}]]]", "\"int64\"");
    chk_eq("Head[Round[NDArray[{2.5}]][[1]]]", "Integer");
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
    /* Quotient NARROWS: real in, exact Integer out. This line asserted
     * "{0.0, 0.0, 1.0, 1.0, 1.0}" until 2026-08-01 and so enshrined the defect
     * it was meant to guard -- the List path has always answered
     * Quotient[{1., 2., 3., 4., 5.}, 3] with the exact {0, 0, 1, 1, 1}
     * (builtin_quotient ends in expr_new_integer on every branch), and so does
     * Mathematica, whose heads for that call are all Integer. Only the buffer
     * path disagreed, because NDKB_Quotient had no narrowing arm and the
     * real-closed branch wrote the double straight back into a Real buffer.
     * Mod is the contrast and is correct as it stands: its exactness FOLLOWS
     * the argument, so a Real array gives Real results. */
    chk_eq("Quotient[NDArray[{1., 2., 3., 4., 5.}], 3]", "NDArray[{0, 0, 1, 1, 1}]");
    /* Floored: negatives carry the divisor's sign. */
    chk_eq("Mod[NDArray[{-1., -2., -3.}], 3]", "NDArray[{2.0, 1.0, 0.0}]");
    chk_eq("Quotient[NDArray[{-1., -2., -3.}], 3]", "NDArray[{-1, -1, -1}]");
    /* Complex operand / zero divisor decline -> faithful List degrade. */
    chk_eq("Head[Mod[NDArray[{1. + 1. I}, DataType -> \"complex64\"], 3]]", "List");

    /* The int64 arms. An integer buffer keeps its dtype and its exact Integer
     * elements through both heads -- the property that lets pack.c carry
     * "Mod" and "Quotient" in INT64_OK, which is what keeps an integer pipeline
     * packed for the set operations downstream of it. */
    chk_eq("Mod[NDArray[{-7, -1, 0, 1, 7}, DataType -> \"int64\"], 3]",
           "NDArray[{2, 2, 0, 1, 1}]");
    chk_eq("Mod[NDArray[{-7, -1, 0, 1, 7}, DataType -> \"int64\"], -3]",
           "NDArray[{-1, -1, 0, -2, -2}]");
    chk_eq("Quotient[NDArray[{-7, -1, 0, 1, 7}, DataType -> \"int64\"], 3]",
           "NDArray[{-3, -1, 0, 0, 2}]");
    chk_eq("Quotient[NDArray[{-7, -1, 0, 1, 7}, DataType -> \"int64\"], -3]",
           "NDArray[{2, 0, 0, -1, -3}]");
    /* A REAL scalar over an integer buffer has no int64 answer; the kernel
     * declines and the List path supplies the Real elements. */
    chk_eq("Mod[NDArray[{1, 2, 3}, DataType -> \"int64\"], 2.]", "{1.0, 0.0, 1.0}");
    /* Zero divisor: the kernel declines and the List path answers, which for a
     * Listable head is a List of unevaluated Mod calls -- exactly what
     * Mod[{1, 2, 3}, 0] gives. The point is the agreement, not the shape. */
    chk_eq("Mod[NDArray[{1, 2, 3}, DataType -> \"int64\"], 0]",
           "{Mod[1, 0], Mod[2, 0], Mod[3, 0]}");

    /* An int64 buffer with a scalar that is not an exact Integer must NOT take
     * Quotient's REAL arm: that arm reads elements through ndt_get, which is
     * exact only to 2^53. 2^53 + 1 rounds to 2^53 as a double, so the real arm
     * would answer 18014398509481984 where the List path divides in mpq and
     * gives ...986. The kernel declines and the List path answers. */
    chk_eq("Quotient[NDArray[{9007199254740993}, DataType -> \"int64\"], 1/2]",
           "{18014398509481986}");
}

/* An int64 buffer through a kernel with no EXACT arm must degrade to the List
 * path, not compute in double and truncate back into the integer slot.
 *
 * Before 2026-08-01 every real_closed kernel did the latter on a VISIBLE
 * NDArray -- Sin gave {0, 0, 0}, Exp[{1,2,3}] gave {2, 7, 20}, and Erf, Tanh,
 * Cos, BesselJ and Log[b, .] likewise -- silently. The packed representation
 * was never exposed to it because src/eval.c's transparency gate materialises
 * an int64 buffer for any head that has not claimed packed_int64_ok; the
 * visible surface has no such gate, and nothing else stood between an integer
 * buffer and a double kernel. */
static void test_int64_ndarray_no_truncation(void) {
    chk_eq("Sin[NDArray[{1, 2, 3}, DataType -> \"int64\"]]",
           "{Sin[1], Sin[2], Sin[3]}");
    chk_eq("Exp[NDArray[{1, 2, 3}, DataType -> \"int64\"]]", "{E, E^2, E^3}");
    chk_eq("Tanh[NDArray[{1, 2, 3}, DataType -> \"int64\"]]",
           "{Tanh[1], Tanh[2], Tanh[3]}");
    chk_eq("Erf[NDArray[{1, 2, 3}, DataType -> \"int64\"]]",
           "{Erf[1], Erf[2], Erf[3]}");
    chk_eq("BesselJ[0, NDArray[{1, 2, 3}, DataType -> \"int64\"]]",
           "{BesselJ[0, 1], BesselJ[0, 2], BesselJ[0, 3]}");
    chk_eq("Log[2, NDArray[{1, 2, 3}, DataType -> \"int64\"]]",
           "{0, 1, Log[3]/Log[2]}");
    /* The kernels that DO have an exact arm keep their buffer fast path. */
    chk_eq("Floor[NDArray[{1, 2, 3}, DataType -> \"int64\"]]", "NDArray[{1, 2, 3}]");
    chk_eq("Abs[NDArray[{-1, 2, -3}, DataType -> \"int64\"]]", "NDArray[{1, 2, 3}]");
    chk_eq("Sign[NDArray[{-1, 0, 3}, DataType -> \"int64\"]]", "NDArray[{-1, 0, 1}]");
}

/* The pattern-matching family treats a VISIBLE NDArray as an atom: the matcher
 * walks data.function.args and an EXPR_NDARRAY has none, so each head searched
 * an expression with no elements and answered, confidently, that it found
 * nothing. Before 2026-08-01:
 *
 *     MemberQ[NDArray[Range[1., 300.]], 5.]  ->  False   (List: True)
 *     Count[NDArray[Range[300]], 5]          ->  0       (List: 1)
 *     Position[NDArray[Range[300]], 5]       ->  {}      (List: {{5}})
 *     FreeQ[NDArray[Range[300]], 5]          ->  True    (List: False)
 *
 * A wrong answer that looks like a right one, which is worse than declining.
 * The PACKED form was never affected: these heads are not on pack.c's AWARE
 * list, so the transparency gate materialises their arguments first. Each now
 * materialises a visible array itself.
 *
 * The lists are 300 elements deliberately -- above PACK_MIN_ELEMENTS, so the
 * packed comparison in each pair is genuinely packed rather than a plain List
 * that never exercised the path. */
static void test_patterns_visible_ndarray(void) {
    chk_eq("MemberQ[NDArray[Range[1., 300.]], 5.]", "True");
    chk_eq("MemberQ[NDArray[Range[1., 300.]], 999.]", "False");
    chk_eq("Count[NDArray[Range[300], DataType -> \"int64\"], 5]", "1");
    chk_eq("Position[NDArray[Range[300], DataType -> \"int64\"], 5]", "{{5}}");
    chk_eq("Cases[NDArray[Range[300], DataType -> \"int64\"], 5]", "{5}");
    chk_eq("FreeQ[NDArray[Range[300], DataType -> \"int64\"], 5]", "False");
    chk_eq("FirstCase[NDArray[Range[300], DataType -> \"int64\"], 5]", "5");
    chk_eq("Length[DeleteCases[NDArray[Range[300], DataType -> \"int64\"], 5]]",
           "299");
    /* and the packed form, which the gate has always covered, still agrees */
    chk_eq("MemberQ[ToNDArray[Range[1., 300.]], 5.]", "True");
    chk_eq("Count[ToNDArray[Range[300]], 5]", "1");
}

/* The set operations are a rank-1 int64 fast path, and it must not be reachable
 * only through the PACKED representation. `setop_i64` tested is_packed_list
 * until 2026-08-01, so a visible NDArray fell through to
 * ndarray_delist_and_reeval and ran the boxed hash path over a materialised
 * List: at 10^6 elements Union was 5.85 ms packed against 850 ms visible and
 * DeleteDuplicates 2.05 ms against 147 ms, on identical data. */
static void test_setops_visible_ndarray(void) {
    chk_eq("Union[NDArray[{3, 1, 2, 1}, DataType -> \"int64\"]]",
           "NDArray[{1, 2, 3}]");
    chk_eq("DeleteDuplicates[NDArray[{3, 1, 2, 1}, DataType -> \"int64\"]]",
           "NDArray[{3, 1, 2}]");
    chk_eq("Intersection[NDArray[{1, 2, 3}, DataType -> \"int64\"], "
           "NDArray[{2, 3, 4}, DataType -> \"int64\"]]", "NDArray[{2, 3}]");
    chk_eq("Complement[NDArray[{1, 2, 3}, DataType -> \"int64\"], "
           "NDArray[{2}, DataType -> \"int64\"]]", "NDArray[{1, 3}]");
    /* Same values, same answers, through the packed representation. */
    chk_eq("Union[ToNDArray[{3, 1, 2, 1}]]", "{1, 2, 3}");
    chk_eq("DeleteDuplicates[ToNDArray[{3, 1, 2, 1}]]", "{3, 1, 2}");
}


/* ---------------------------------------------------------------------------
 * NINTH ROUND: the 26 heads that declined a visible NDArray outright.
 *
 * The eighth round's audit classified these as ND-UNSUPPORTED -- they answered
 * on a plain List and on a packed List and left the call UNEVALUATED on a
 * visible NDArray. They now have buffer fast paths. Every assertion below is
 * the value the List path gives, because that is the whole contract: a
 * representation may never change a value (pack.h), and these tests exist to
 * hold the three surfaces to one answer.
 * ------------------------------------------------------------------------- */

/* Batch A -- ndarray_map_binary2. The dispatch required exactly ONE array
 * operand, so a two-argument kernel was unreachable when both arguments were
 * arrays, which is the ordinary way to call one. Fifteen registered kernels
 * were in that position; ArcTan and Beta are the two the sweep probes. */
static void test_binary_two_arrays(void) {
    chk_eq("Normal[ArcTan[NDArray[{1., 0., -1.}], NDArray[{1., 1., 1.}]]] === "
           "ArcTan[{1., 0., -1.}, {1., 1., 1.}]", "True");
    chk_eq("Head[ArcTan[NDArray[{1., 2.}], NDArray[{3., 4.}]]]", "NDArray");
    /* atan2, not the csqrt/clog form: bit-identical to the scalar builtin,
     * which the old kernel was not (68 of 400 elements differed by 1-2 ulp). */
    chk_eq("Normal[ArcTan[NDArray[{1., 2., 3.}], NDArray[{3., 2., 1.}]]] === "
           "ArcTan[{1., 2., 3.}, {3., 2., 1.}]", "True");
    chk_eq("Head[Beta[NDArray[{1., 2.}], NDArray[{3., 4.}]]]", "NDArray");
    /* Mismatched shapes must not answer: same Thread::tdlen as the List path. */
    chk_eq("Head[ArcTan[NDArray[{1., 2., 3.}], NDArray[{1., 2.}]]]", "ArcTan");
    /* The exact-integer arm carries over to two arrays. */
    chk_eq("Mod[NDArray[{7, 8, 9}, DataType -> \"int64\"], "
           "NDArray[{2, 3, 4}, DataType -> \"int64\"]]", "NDArray[{1, 2, 1}]");
}

/* Batch B -- the integer domain. int64 in, int64 out, through ci_*_i64, and
 * ABANDONING the whole array rather than wrapping. Each assertion is checked
 * against the scalar builtin's own convention: the sign of n is ignored by
 * EulerPhi / MoebiusMu / DivisorSigma, IntegerLength[0] is 0, MoebiusMu[0] and
 * DivisorSigma[k, 0] stay unevaluated. */
static void test_integer_kernels(void) {
    chk_eq("IntegerLength[NDArray[{0, 1, 9, 10, -345}, DataType -> \"int64\"]]",
           "NDArray[{0, 1, 1, 2, 3}]");
    chk_eq("EulerPhi[NDArray[{1, 2, 12, 97}, DataType -> \"int64\"]]",
           "NDArray[{1, 1, 4, 96}]");
    chk_eq("MoebiusMu[NDArray[{1, 2, 4, 6}, DataType -> \"int64\"]]",
           "NDArray[{1, -1, 0, 1}]");
    chk_eq("GCD[NDArray[{12, 0, -12}, DataType -> \"int64\"], 18]",
           "NDArray[{6, 18, 6}]");
    chk_eq("LCM[NDArray[{4, 0, -4}, DataType -> \"int64\"], 6]",
           "NDArray[{12, 0, 12}]");
    chk_eq("DivisorSigma[1, NDArray[{1, 6, -6}, DataType -> \"int64\"]]",
           "NDArray[{1, 12, 12}]");
    chk_eq("DivisorSigma[0, NDArray[{12}, DataType -> \"int64\"]]", "NDArray[{6}]");
    chk_eq("Normal[PowerMod[NDArray[{2, 3, 4}, DataType -> \"int64\"], 3, 1009]] "
           "=== PowerMod[{2, 3, 4}, 3, 1009]", "True");
    /* Prime over an array is a SIEVE, not one Meissel count per element. */
    chk_eq("Prime[NDArray[{1, 10, 100}, DataType -> \"int64\"]]",
           "NDArray[{2, 29, 541}]");
    chk_eq("IntegerDigits[NDArray[{0, 123, -123}, DataType -> \"int64\"]]",
           "{{0}, {1, 2, 3}, {1, 2, 3}}");
    /* Arguments the scalar builtin leaves unevaluated abandon the whole array,
     * so the List path reproduces it element for element. */
    chk_eq("Normal[MoebiusMu[NDArray[{1, 0}, DataType -> \"int64\"]]] === "
           "MoebiusMu[{1, 0}]", "True");
    chk_eq("Normal[DivisorSigma[1, NDArray[{6, 0}, DataType -> \"int64\"]]] === "
           "DivisorSigma[1, {6, 0}]", "True");
    /* A REAL buffer is not this domain: EulerPhi[3.] is not EulerPhi[3]. */
    chk_eq("Normal[EulerPhi[NDArray[{1., 2., 3.}]]] === EulerPhi[{1., 2., 3.}]",
           "True");
    /* Overflow abandons rather than wrapping -- GMP answers exactly. */
    chk_eq("Normal[LCM[NDArray[{1000000000000000000, 2}, DataType -> \"int64\"], "
           "999999999]] === LCM[{1000000000000000000, 2}, 999999999]", "True");
    /* Both operands arrays, through ndarray_map_binary2's integer arm. */
    chk_eq("GCD[NDArray[{12, 18}, DataType -> \"int64\"], "
           "NDArray[{18, 24}, DataType -> \"int64\"]]", "NDArray[{6, 6}]");
}

/* Batch C -- the sign predicates. The answer is a List of True/False, which no
 * dtype holds; the win is that the INPUT is read off the buffer. An
 * Indeterminate element is unordered and must stay symbolic, not answer False,
 * so it abandons the array. */
static void test_sign_predicates(void) {
    chk_eq("Positive[NDArray[{-1., 0., 1.}]]", "{False, False, True}");
    chk_eq("Negative[NDArray[{-1., 0., 1.}]]", "{True, False, False}");
    chk_eq("NonNegative[NDArray[{-1., 0., 1.}]]", "{False, True, True}");
    chk_eq("NonPositive[NDArray[{-1., 0., 1.}]]", "{True, True, False}");
    chk_eq("Positive[NDArray[{-1, 0, 1}, DataType -> \"int64\"]]",
           "{False, False, True}");
    chk_eq("Normal[Positive[ToNDArray[Table[1., {300}]]]] === "
           "Positive[Table[1., {300}]]", "True");
}

/* Batch D -- structural. */
static void test_structural_ninth(void) {
    chk_eq("Ratios[NDArray[{1., 2., 4., 8.}]]", "NDArray[{2.0, 2.0, 2.0}]");
    /* Ratios of INTEGERS are exact Rationals, which no buffer holds. */
    chk_eq("Ratios[NDArray[{1, 2, 3}, DataType -> \"int64\"]]", "{2, 3/2}");
    chk_eq("Append[NDArray[{1., 2.}], 3.]", "NDArray[{1.0, 2.0, 3.0}]");
    chk_eq("Prepend[NDArray[{1., 2.}], 0.]", "NDArray[{0.0, 1.0, 2.0}]");
    chk_eq("Append[NDArray[{1, 2}, DataType -> \"int64\"], 3]", "NDArray[{1, 2, 3}]");
    /* An EXACT append into a float64 buffer must NOT be coerced: the answer is
     * mixed exact/inexact and no buffer holds it. Repacking here gave 0. */
    chk_eq("Normal[Append[NDArray[{1., 2.}], 0]] === Append[{1., 2.}, 0]", "True");
    chk_eq("Catenate[{NDArray[{1., 2.}], NDArray[{3., 4.}]}]",
           "NDArray[{1.0, 2.0, 3.0, 4.0}]");
    chk_eq("Catenate[NDArray[{{1., 2.}, {3., 4.}}]]", "NDArray[{1.0, 2.0, 3.0, 4.0}]");
    chk_eq("TakeLargest[NDArray[{3., 1., 4., 1., 5.}], 3]", "NDArray[{5.0, 4.0, 3.0}]");
    chk_eq("TakeSmallest[NDArray[{3., 1., 4., 1., 5.}], 3]", "NDArray[{1.0, 1.0, 3.0}]");
    chk_eq("TakeLargest[NDArray[{3., 1., 4.}], 3]", "NDArray[{4.0, 3.0, 1.0}]");
    chk_eq("Counts[NDArray[{1, 2, 2, 3, 3, 3}, DataType -> \"int64\"]] === "
           "Counts[{1, 2, 2, 3, 3, 3}]", "True");
    chk_eq("Counts[NDArray[{1., 2., 2.}]] === Counts[{1., 2., 2.}]", "True");
    /* Inner[Times, a, b, Plus] IS a Dot; any other operator pair degrades
     * rather than falling out unevaluated. */
    chk_eq("Inner[Times, NDArray[{1., 2., 3.}], NDArray[{4., 5., 6.}], Plus]", "32.0");
    chk_eq("Inner[Plus, NDArray[{1., 2.}], NDArray[{3., 4.}], Times] === "
           "Inner[Plus, {1., 2.}, {3., 4.}, Times]", "True");
    /* Rank-2 Append takes a whole matching row. */
    chk_eq("Append[NDArray[{{1., 2.}, {3., 4.}}], NDArray[{5., 6.}]]",
           "NDArray[{{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}}]");
}

/* RandomSample / RandomChoice draw from the SAME generator sequence the List
 * path uses, so a seeded run agrees across surfaces. A producer that consumed
 * the generator differently would make SeedRandom[1]; RandomSample[v] give one
 * permutation for a packed v and another for the identical plain v. */
static void test_random_buffer_paths(void) {
    chk_eq("SeedRandom[42]; a = RandomSample[NDArray[{1., 2., 3., 4., 5.}]]; "
           "SeedRandom[42]; b = RandomSample[{1., 2., 3., 4., 5.}]; Normal[a] === b",
           "True");
    chk_eq("SeedRandom[7]; a = RandomSample[NDArray[{1., 2., 3., 4., 5.}], 3]; "
           "SeedRandom[7]; b = RandomSample[{1., 2., 3., 4., 5.}, 3]; Normal[a] === b",
           "True");
    chk_eq("SeedRandom[9]; a = RandomChoice[NDArray[{1., 2., 3.}], 8]; "
           "SeedRandom[9]; b = RandomChoice[{1., 2., 3.}, 8]; Normal[a] === b",
           "True");
    chk_eq("Head[RandomSample[NDArray[{1., 2., 3.}]]]", "NDArray");
}

/* Batch E -- the hypergeometric family, which rewrites to HypergeometricPFQ
 * before anything else sees it. The buffer path takes ONLY inexact parameters:
 * try_cancel / try_terminate / try_reduce answer with exact or closed forms and
 * fire on exact parameters, so those must keep reaching the scalar steps. */
static void test_hypergeometric_buffer(void) {
    chk_eq("Head[Hypergeometric1F1[0.5, 1.5, NDArray[{0.5, 1., 1.5}]]]", "NDArray");
    chk_eq("Head[Hypergeometric2F1[0.5, 0.5, 1.5, NDArray[{0.1, 0.2}]]]", "NDArray");
    chk_eq("Max[Abs[Flatten[Normal[Hypergeometric1F1[0.5, 1.5, NDArray[{0.5, 1.}]]] "
           "- Hypergeometric1F1[0.5, 1.5, {0.5, 1.}]]]] < 1/1000000000", "True");
    /* Exact parameters keep the closed form -- 2F1[1,1,2,z] is -Log[1-z]/z. */
    chk_eq("Normal[Hypergeometric2F1[1, 1, 2, NDArray[{0.5, 0.25}]]] === "
           "Hypergeometric2F1[1, 1, 2, {0.5, 0.25}]", "True");
    /* A terminating upper parameter is a polynomial, not the series. */
    chk_eq("Normal[Hypergeometric2F1[-2, 1, 3, NDArray[{0.5, 0.25}]]] === "
           "Hypergeometric2F1[-2, 1, 3, {0.5, 0.25}]", "True");
    /* An exact-integer buffer is evaluated in its own right by the List path. */
    chk_eq("Normal[Hypergeometric1F1[0.5, 1.5, NDArray[{1, 2}, "
           "DataType -> \"int64\"]]] === Hypergeometric1F1[0.5, 1.5, {1, 2}]",
           "True");
}

/* Flatten treated a visible NDArray nested in an ordinary List as an atom --
 * head_is is false for an EXPR_NDARRAY -- so Flatten[{ndA, ndB}] answered with
 * the list unchanged. It is a list of values by every other measure, and this
 * is what made the sweep read small-matrix Dot/Inverse/LinearSolve as
 * unsupported: the heads were right and the probe's checksum could not reduce. */
static void test_flatten_nested_ndarray(void) {
    chk_eq("Flatten[{NDArray[{1., 2.}], NDArray[{3., 4.}]}]",
           "{1.0, 2.0, 3.0, 4.0}");
    chk_eq("Total[Flatten[{NDArray[{1., 2.}], NDArray[{3., 4.}]}]]", "10.0");
    chk_eq("Flatten[{NDArray[{{1., 2.}, {3., 4.}}]}]", "{1.0, 2.0, 3.0, 4.0}");
    /* A level spec still bounds the descent. */
    chk_eq("Flatten[{NDArray[{{1., 2.}, {3., 4.}}]}, 1]",
           "{{1.0, 2.0}, {3.0, 4.0}}");
    /* Flatten of the array itself is still the reshape, and stays an NDArray. */
    chk_eq("Flatten[NDArray[{{1., 2.}, {3., 4.}}]]", "NDArray[{1.0, 2.0, 3.0, 4.0}]");
}

/* ---------------------------------------------------------------------------
 * ORDERING: a packed list compared against a scalar must not materialise, and
 * must give the same answer as if it had.
 *
 * expr_compare orders a packed list as the List it is -- correct -- but reached
 * that through ndarray_to_nested_list whenever only ONE operand was packed. The
 * evaluator sorts `GCD[1234, cv]` before dispatch (GCD is Orderless), so 200000
 * Expr nodes were built and discarded per call to settle an ordering that the
 * type tiers decide without reading an element: GCD 46.7 ms packed against
 * 18.1 ms visible, LCM 31.5 against 8.3, while Mod and Quotient -- the same
 * kernels WITHOUT ATTR_ORDERLESS -- showed no gap at all.
 *
 * These pin the answer, not the speed. Every tier, both directions.
 * ------------------------------------------------------------------------- */
static void test_packed_scalar_ordering(void) {
    /* A packed list sorts AFTER every atom tier: numbers, complex, strings. */
    chk_eq("OrderedQ[{5, Table[1., {300}]}]", "True");
    chk_eq("OrderedQ[{Table[1., {300}], 5}]", "False");
    chk_eq("OrderedQ[{-2, Table[1., {300}]}]", "True");
    chk_eq("OrderedQ[{7/3, Table[1., {300}]}]", "True");
    chk_eq("OrderedQ[{2.5, Table[1., {300}]}]", "True");
    chk_eq("OrderedQ[{3 + 4 I, Table[1., {300}]}]", "True");
    chk_eq("OrderedQ[{2.5 + 0.5 I, Table[1., {300}]}]", "True");
    chk_eq("OrderedQ[{\"abc\", Table[1., {300}]}]", "True");
    chk_eq("OrderedQ[{Table[1., {300}], \"abc\"}]", "False");
    /* ... and BEFORE a Symbol and a general expression, which still go through
     * the materialising path (step 3 reads elements). */
    chk_eq("OrderedQ[{Table[1., {300}], xzz}]", "True");
    chk_eq("OrderedQ[{Table[1., {300}], Sin[xzz]}]", "True");
    /* An integer buffer behaves identically. */
    chk_eq("OrderedQ[{5, Range[300]}]", "True");
    chk_eq("OrderedQ[{Range[300], 5}]", "False");
    /* Two packed lists still compare on their buffers. */
    chk_eq("OrderedQ[{Table[1., {300}], Table[2., {300}]}]", "True");
    chk_eq("OrderedQ[{Table[2., {300}], Table[1., {300}]}]", "False");
    /* Sort of a mixed bag: the full canonical order, in one assertion. */
    chk_eq("Map[Head, Sort[{Table[1., {300}], 5, \"abc\", xzz, 3 + 4 I, 2.5}]]",
           "{Real, Integer, Complex, String, List, Symbol}");
    /* Orderless canonicalisation through a real head still sees a full list. */
    chk_eq("Length[Plus[Table[1., {300}], 5]]", "300");
    /* And the values GCD/LCM produce are unaffected by the sort shortcut. */
    chk_eq("Total[GCD[Range[300], 12]] == Total[GCD[12, Range[300]]]", "True");
}

/* PACK_MIN_ELEMENTS is 4 -- a 2x2 matrix -- so any matrix packs and the LAPACK
 * path is reachable from an ordinary computed list. It was 250, which left a
 * 6x6 unpacked and Det[A6] at 102.8 ms against 0.189 ms. A LITERAL list is a
 * separate gate and still never packs: packing is opt-in per producer. */
static void test_pack_threshold(void) {
    chk_eq("NDArrayQ[Table[1.*i, {i, 4}]]", "True");
    chk_eq("NDArrayQ[Table[If[i == j, 4., 1.], {i, 2}, {j, 2}]]", "True");
    chk_eq("NDArrayQ[Table[1.*i, {i, 3}]]", "False");     /* under the threshold */
    chk_eq("NDArrayQ[Range[1., 6.]]", "True");
    /* A literal never packs, whatever its size -- no producer to opt in. */
    chk_eq("NDArrayQ[{1., 2., 3., 4., 5., 6.}]", "False");
    /* Packing an ordinary 2x2 must not change a single value. */
    chk_eq("Det[Table[If[i == j, 4., 1.], {i, 2}, {j, 2}]]", "15.0");
    chk_eq("Inverse[Table[If[i == j, 2., 0.], {i, 2}, {j, 2}]]",
           "{{0.5, 0.0}, {0.0, 0.5}}");
    /* Exactness is unchanged at the low threshold: an exact 2x2 stays exact. */
    chk_eq("Det[Table[If[i == j, 4, 1], {i, 2}, {j, 2}]]", "15");
    chk_eq("Total[Range[4]]", "10");
    chk_eq("Mean[Range[4]]", "5/2");
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
    TEST(test_int64_ndarray_no_truncation);
    TEST(test_patterns_visible_ndarray);
    TEST(test_setops_visible_ndarray);
    TEST(test_special_cloop);
    TEST(test_special_degrade);
    TEST(test_binary_two_arrays);
    TEST(test_integer_kernels);
    TEST(test_sign_predicates);
    TEST(test_structural_ninth);
    TEST(test_random_buffer_paths);
    TEST(test_hypergeometric_buffer);
    TEST(test_flatten_nested_ndarray);
    TEST(test_packed_scalar_ordering);
    TEST(test_pack_threshold);

    printf("All NDArray function tests passed.\n");
    return 0;
}
