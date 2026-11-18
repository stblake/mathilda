/* Mathilda — regression corpus for the 2026-11 numeric-tower review
 * (docs/design/numeric_tower_review.md): the int64 <-> BigInt boundary, the
 * machine <-> MPFR boundary, and fixed-precision MPFR edge cases.
 *
 * Each section is keyed by the review's finding id (A = crash / undefined
 * behaviour, B = exact-layer wrong answer, C = machine/software interplay,
 * D = arbitrary precision). Expected values are derived from first principles
 * -- exact identities, checked with SameQ wherever the printed form is not
 * itself the point -- never from another system's output.
 *
 * ASSERT_MSG (test_utils.h) aborts unconditionally, unlike assert(), which a
 * CMake Release build compiles out.
 */

#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

static int g_fail = 0;

/* Evaluate `input`; its standard printed form must equal `expected`. */
static void chk(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    ASSERT_MSG(e != NULL, "parse failed: %s", input);
    Expr* res = evaluate(e);
    char* s = expr_to_string(res);
    if (strcmp(s, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, s);
        g_fail++;
    }
    free(s);
    expr_free(e); expr_free(res);
}

/* `input` must evaluate to True. */
static void chk_true(const char* input) { chk(input, "True"); }

/* ---- A: crashes and undefined behaviour --------------------------------- */

static void test_a1_int64_min_rational_power(void) {
    /* (-a)^r = (-1)^r a^r for real a > 0 (principal branch); 2^63 = 2^62 * 2,
     * so Sqrt[2^63] = 2^31 Sqrt[2]. Previously recursed on -INT64_MIN and
     * segfaulted. */
    chk_true("Sqrt[-2^63] === 2147483648 I Sqrt[2]");
    chk_true("(-2^63)^(1/3) === 2097152 (-1)^(1/3)");      /* 2^(63/3) = 2^21 */
    chk_true("(-2^63)^(3/2) === -I 2^94 Sqrt[2]");          /* (-1)^(3/2) = -I */
    chk_true("Sqrt[-2^62] === 2147483648 I");               /* 2^31 I */
}

static void test_a2_quotient_huge_real(void) {
    /* floor of an integer-valued double outside int64 is converted exactly
     * (it was a UB cast giving -2^63). 1.5*^20 = 3 * 5^20 * 2^19 is exact. */
    chk_true("Quotient[1.*^300, 1] === Floor[1.*^300]");
    chk("Quotient[1.5*^20, 3]", "50000000000000000000");
    chk_true("Quotient[-1.*^300, 7.] === Floor[-1.*^300/7.]");
    chk("Quotient[7.5, 2]", "3");
    chk("Quotient[-7.5, 2]", "-4");
}

static void test_a5_power_past_gmp_range(void) {
    /* A result past GMP's representable range (INT_MAX limbs) aborted the
     * process inside mpz_pow_ui; now General::ovfl and Overflow[]. */
    chk("Quiet[3^(2^62)]", "Overflow[]");
    chk("Quiet[(2^1000/3)^(10^9)]", "Overflow[]");
    chk("Quiet[2^(2^70)]", "Overflow[]");
    /* Small cases still exact. */
    chk("2^100", "1267650600228229401496703205376");
    chk("(2/3)^5", "32/243");
}

static void test_a6_print_int64_min_exponent(void) {
    /* The printer flipped a negative exponent with int64 negation; -INT64_MIN
     * wraps, so x^(-2^63/3) printed 1/(1/(1/(...))) until the stack ran out. */
    chk("x^(-2^63/3)", "1/x^(9223372036854775808/3)");
    chk("x^(-2^63)", "1/x^9223372036854775808");
}

/* ---- B: exact-layer wrong answers ---------------------------------------- */

static void test_b1_divide_cross_product_overflow(void) {
    /* (1/3^39)/3^39 = 1/3^78; the int64 cross product wrapped negative. */
    chk_true("Divide[1/3^39, 3^39] === 1/3^78");
    chk_true("Divide[3^39, 1/3^39] === 3^78");
    chk_true("Divide[2^40/3, 3/2^40] === 2^80/9");
    chk("Divide[6, 4]", "3/2");
}

static void test_b2_gcd_lcm_int64_min(void) {
    chk("GCD[-2^63]", "9223372036854775808");
    chk("GCD[-2^63, 0]", "9223372036854775808");
    chk("GCD[-2^63, 6]", "2");
    chk("LCM[-2^63]", "9223372036854775808");
    chk_true("LCM[-2^63, 3] === 3 * 2^63");
    chk("GCD[-1/2]", "1/2");
}

static void test_b3_negative_denominator(void) {
    /* A canonical Rational has a positive denominator. 3^41 = 36472996377170786403. */
    chk("FullForm[(-3)^-41]", "Rational[-1, 36472996377170786403]");
    chk("FullForm[(-2^70)^-1]", "Rational[-1, 1180591620717411303424]");
    chk("FullForm[(-3)^-3]", "Rational[-1, 27]");
    chk_true("(-3)^-41 * (-3)^41 === 1");
}

static void test_b4_rational_bignum_canonical(void) {
    chk("Rational[2^70, 2^70]", "1");
    chk("FullForm[Rational[6, -2^70]]", "Rational[-3, 590295810358705651712]");
    chk_true("Rational[2^70 + 2, 4] === (2^69 + 1)/2");
    chk("Rational[0, 2^70]", "0");
}

static void test_b5_bigint_hash_matches_integer(void) {
    /* A 21-digit literal with leading zeros is the Integer 1 (parser
     * normalises); and any int64-sized BigInt hashes as its Integer, so the
     * hash agrees with SameQ. */
    chk("Head[000000000000000000001]", "Integer");
    chk("Length[Union[{1, 000000000000000000001}]]", "1");
    chk("<|1 -> a|>[000000000000000000001]", "a");
    chk_true("Hash[1] == Hash[000000000000000000001]");
    chk_true("Hash[2^62] == Hash[2^61 * 2]");
}

static void test_b6_bigint_exponent(void) {
    chk("(-1)^(2^70)", "1");
    chk("(-1)^(2^70 + 1)", "-1");
    chk("(-1)^(-2^70 - 1)", "-1");
    chk("1^(2^70)", "1");
    chk("1^(-2^70)", "1");
    chk("0^(2^70)", "0");
}

/* ---- C: machine <-> software boundary ------------------------------------ */

/* A machine number has a 53-bit mantissa and a normal exponent; a result
 * outside [DBL_MIN, DBL_MAX] is the same value at machine precision (53-bit
 * MPFR, MPFR's exponent range) -- never inf, a flushed 0, or a subnormal. */
static void test_c1_parser_degenerate_literals(void) {
    chk_true("1.*^400 == 10^400");
    chk_true("1.*^-400 == 10^-400");
    chk_true("1.*^-400 > 0");
    chk_true("-1.*^400 == -10^400");
    chk_true("Precision[1.*^400] == Precision[Exp[1000.]]");
    chk_true("1.*^-320 == 10^-320");                  /* not a 4-digit subnormal */
    chk("FreeQ[Hold[1.*^400], inf]", "True");
    chk("0.*^400", "0.0");
}

static void test_c2_plus_times_overflow_underflow(void) {
    chk_true("1.*^308 + 1.*^308 == 2 10^308");
    chk_true("10.^308 * 10 == 10^309");
    chk_true("1.*^-200 * 1.*^-200 == 10^-400");        /* was 0. */
    chk_true("2.^-1074 / 2 == 2^-1075");
    chk_true("1.5 + 10^400 == 10^400");
    chk_true("Precision[1.*^308 + 1.*^308] == Precision[Exp[1000.]]");
    chk("1. - 1.", "0.0");                              /* a true machine zero */
    chk("MachineNumberQ[1.5 * 2.]", "True");
    /* packed / listable path declines to the builtin, element for element */
    chk_true("{1.*^308, 1.} + {1.*^308, 1} == {2 10^308, 2}");
}

static void test_c4_negative_base_machine_power(void) {
    /* sin/cos of pi y reduced exactly: no 6.1*^-17 residue on any surface */
    chk("Sqrt[-1.]", "0.0 + 1.0*I");
    chk("Sqrt[-4.]", "0.0 + 2.0*I");
    chk("(-1.)^-0.5", "0.0 - 1.0*I");
    chk("Sqrt[Developer`ToPackedArray[{-1., 4.}]]", "{0.0 + 1.0*I, 2.0 + 0.0*I}");
    chk("Sqrt[NDArray[{-1., 4.}]]", "NDArray[{0.0 + 1.0*I, 2.0 + 0.0*I}]");
    chk("Compile[{{x, _Real}}, Sqrt[x]][-1.]", "0.0 + 1.0*I");
    chk_true("Sqrt[-2.] === Sqrt[Developer`ToPackedArray[{-2.}]][[1]]");
}

static void test_c6_machine_contagion(void) {
    chk("Precision[N[Pi, 30]^1.5]", "MachinePrecision");
    chk("Precision[1.5^N[Pi, 30]]", "MachinePrecision");
    chk("Precision[N[Pi, 30]^2]", "30.0");             /* exact exponent: no contagion */
    chk("Precision[Abs[Complex[N[Pi, 30], 1.5]]]", "MachinePrecision");
    chk_true("Complex[1.0, 2] === Complex[1.0, 2.0]");
    chk("Precision[Complex[1.5`30, 2.]]", "MachinePrecision");
    chk("Head[Complex[2^70, 0.]]", "Real");
    chk("FullForm[Complex[0, 1.5]]", "Complex[0, 1.5]");   /* exact zero part kept */
}

static void test_c7_compare_out_of_double_range(void) {
    chk("1.5 < 10^400", "True");
    chk("1.*^308 < 10^309", "True");
    chk("1.*^-300 > 1/10^400", "True");
    chk("-1.*^308 < -10^400", "False");
    chk("0. == 1/10^400", "False");
    chk_true("Sort[{N[10^500], N[10^400], 1.5}] == {1.5, 10^400, 10^500}");
    chk_true("Sort[{2^53 + 1, 2.^53, 2^53}] === {2^53, 2.^53, 2^53 + 1}");
    chk_true("TakeLargest[{N[10^400], N[10^500], 3.}, 1] == {10^500}");
}

static void test_c8_no_subnormal_machine_numbers(void) {
    chk_true("N[10^-320] == 10^-320");
    chk_true("N[3/2 2^-1074] == 3/2 2^-1074");
    chk_true("SetPrecision[2.^-1074, Infinity] === 2^-1074");
    chk("MachineNumberQ[2.^-1074]", "False");
    chk("MachineNumberQ[2.^-1022]", "True");            /* DBL_MIN is a machine number */
}

/* ---- Phase 4 (v0.6595): the remaining wrong answers ---------------------- */

static void test_c3_e5_buffer_paths_promote(void) {
    /* Packed / NDArray Plus, Times and the reductions decline a lane or a
     * scalar outside the machine range to the element path, which promotes. */
    chk_true("Developer`ToPackedArray[{1.*^308, 2.}] * 10 == {10^309, 20}");
    chk_true("Developer`ToPackedArray[{1.*^-200, 2.}] * 1.*^-200 == {10^-400, 2 10^-200}");
    chk_true("NDArray[{1.*^308, 2.}] + NDArray[{1.*^308, 2.}] == {2 10^308, 4}");
    chk_true("Developer`ToPackedArray[{1., 2.}] * 10^400 == {10^400, 2 10^400}");
    chk_true("Total[Developer`ToPackedArray[{1.*^308, 1.*^308}]] == 2 10^308");
    chk_true("Accumulate[Developer`ToPackedArray[{1.*^308, 1.*^308}]] == {10^308, 2 10^308}");
    chk_true("Mean[Developer`ToPackedArray[{1.*^308, 1.*^308}]] == 10^308");
    chk_true("Mean[{1.*^308, 1.*^308}] == 10^308");
    chk_true("RootMeanSquare[{1.*^308, 1.*^308}] == 10^308");
    chk_true("RootMeanSquare[{1.*^-200, 1.*^-200}] == 10^-200");
    chk_true("{1.*^308, 1.*^308} . {1.*^308, 1.*^308} == 2 10^616");
    chk_true("{1.*^-200, 1.*^-200} . {1.*^-200, 1.*^-200} == 2 10^-400");
    chk("{1., -1.} . {1., 1.}", "0.0");                      /* cancellation stays */
    chk("Developer`PackedArrayQ[Developer`ToPackedArray[{1., 2.}] * 10.]", "True");
}

static void test_e9_compile_underflow(void) {
    chk_true("Compile[{{x, _Real}}, x x][1.*^-200] == 10^-400");
    chk_true("Compile[{{v, _Real, 1}}, v v][{1.*^-200, 2.}] == {10^-400, 4}");
    chk_true("Compile[{{x, _Complex}}, x x][1.*^-200 I] == -10^-400");
    chk_true("Table[x x, {x, {1.*^-200, 2.}}] == {10^-400, 4}");
    chk("Compile[{{x, _Real}}, x x][3.]", "9.0");
}

/* E12: a _Real argument with no machine representation declines the compiled
 * run instead of entering the VM as inf. / 0. (found by check-machine-boundary). */
static void test_e12_compile_unrepresentable_arg(void) {
    chk("Quiet[Compile[{{x, _Real}, {y, _Real}}, x^y][0.5, 10^400]]", "Underflow[]");
    chk("Quiet[Compile[{{x, _Real}, {y, _Real}}, x^y][0.5, 1.*^400]]", "Underflow[]");
    chk("Quiet[Compile[{{x, _Real, 1}, {y, _Real}}, x^y][{1., 0.5}, 10^400]]",
        "{1.0, Underflow[]}");
    chk_true("Compile[{{x, _Real}}, 2. x][1/10^400] == 2 10^-400");
    chk_true("Compile[{{v, _Real, 1}}, 2. v][{10^400, 1.}] == {2 10^400, 2}");
    chk("Compile[{{x, _Real}}, 2. x][0]", "0.0");            /* exact zero compiles */
    chk_true("Compile[{{x, _Real}}, 2. x][1/3] == 2./3");
    chk_true("Compile[{{x, _Real}}, x^2, WorkingPrecision -> 30][10^400] == 10^800");
}

/* G1-G3: wrong answers on ordinary data, found by check-machine-boundary. */
static void test_g1_compile_integer_subtract(void) {
    chk("Compile[{{x, _Integer, 1}, {y, _Integer}}, Subtract[x, y]][{10, 3}, 2]", "{8, 1}");
    chk("Compile[{{v, _Integer, 1}, {w, _Integer, 1}}, Subtract[v, w]][{2^53 + 1, 3}, {0, 1}]",
        "{9007199254740993, 2}");                            /* stays an integer array */
    chk("Compile[{{y, _Integer}, {v, _Integer, 1}}, Subtract[y, v]][0, {2^53 + 1, 3}]",
        "{-9007199254740993, -3}");
    chk("Compile[{{x, _Integer, 1}, {y, _Integer}}, Subtract[x, y]][{3}, -2^63]",
        "{9223372036854775811}");                            /* overflow -> exact */
    chk("Compile[{{v, _Real, 1}, {y, _Real}}, Subtract[v, y]][{1., 3.}, 2.]", "{-1.0, 1.0}");
}

static void test_g2_inverse_jacobi_dn_packed(void) {
    chk_true("With[{v = {0.3, 0.5, 0.99, -0.2}}, InverseJacobiDN[Developer`ToPackedArray[v], 0.5]"
             " === (InverseJacobiDN[#, 0.5] & /@ v)]");
    chk_true("With[{v = {0.3, 0.5}}, Normal[InverseJacobiDN[NDArray[v], 0.5]]"
             " === (InverseJacobiDN[#, 0.5] & /@ v)]");
    chk_true("Im[InverseJacobiDN[Developer`ToPackedArray[{0.3}], 1.*^-10][[1]]] > 12");
    chk("Developer`PackedArrayQ[InverseJacobiDN[Developer`ToPackedArray[{0.99, 0.95}], 0.5]]",
        "True");                                             /* all-real stays packed */
}

static void test_g3_machine_mod(void) {
    /* the exact remainder, on every surface: m - n Floor[m/n] collapsed to 0.
     * when m/n rounded up to an integer */
    chk_true("Mod[Developer`ToPackedArray[{0.3, 0.5}], 1.*^-10] === {Mod[0.3, 1.*^-10], Mod[0.5, 1.*^-10]}");
    chk_true("Mod[0.5, 1.*^-10] == SetPrecision[Mod[SetPrecision[0.5, Infinity], SetPrecision[1.*^-10, Infinity]], MachinePrecision]");
    chk_true("Compile[{{a, _Real, 1}, {b, _Real}}, Mod[a, b]][{0.5}, 1.*^-10] === {Mod[0.5, 1.*^-10]}");
    chk("Mod[-7., 3.]", "2.0");
    chk("Mod[7., -3.]", "-2.0");
    chk("Mod[Developer`ToPackedArray[{-7., 7.}], -3.]", "{-1.0, -2.0}");
    chk("Mod[-2., 1.]", "0.0");
}

/* G5 + G6: every buffer / compiled kernel lane is checked by default (lane
 * modes, ndarray.h); a lane the scalar net would recompute declines. */
static void test_g56_lane_checks(void) {
    const char* pk = "Developer`ToPackedArray";
    char buf[512];
    /* G6: underflow flushed to 0. on the buffer surfaces */
    snprintf(buf, sizeof buf, "FresnelS[%s[{-1.*^-200, 0.5}]] === {FresnelS[-1.*^-200], FresnelS[0.5]}", pk);
    chk_true(buf);
    chk_true("Normal[FresnelS[NDArray[{-1.*^-200, 0.5}]]] === {FresnelS[-1.*^-200], FresnelS[0.5]}");
    snprintf(buf, sizeof buf, "Quiet[Zeta[%s[{1.*^308, 2.}], 1.*^308]] === Quiet[{Zeta[1.*^308, 1.*^308], Zeta[2., 1.*^308]}]", pk);
    chk_true(buf);
    chk_true("Quiet[SphericalBesselJ[NDArray[{1.*^308, 1.}], 10.]] === Quiet[{SphericalBesselJ[1.*^308, 10.], SphericalBesselJ[1., 10.]}]");
    /* G6: Compile -- a kernel lane (Csch), an overflowed intermediate, erfc */
    chk_true("Quiet[Compile[{{x, _Real}}, Csch[x]][-1.*^308] === Csch[-1.*^308]]");
    chk_true("Quiet[Compile[{{x, _Real}}, Sech[x]][1.*^308] === Sech[1.*^308]]");
    chk_true("Quiet[Compile[{{x, _Real}}, 1/Sinh[x]][800.] === 1/Sinh[800.]]");
    chk_true("Quiet[Compile[{{x, _Real}}, Erfc[x]][1.*^308] === Erfc[1.*^308]]");
    chk_true("Quiet[Compile[{{v, _Real, 1}}, Erfc[v]][{1.*^308, 1.}] === {Erfc[1.*^308], Erfc[1.]}]");
    /* G5: no subnormal machine number out of a buffer */
    snprintf(buf, sizeof buf, "ArcTanh[%s[{1.*^308, 0.5}]] === {ArcTanh[1.*^308], ArcTanh[0.5]}", pk);
    chk_true(buf);
    snprintf(buf, sizeof buf, "ArcCot[%s[{1.*^308, 0.5}]][[1]] === ArcCot[1.*^308]", pk);
    chk_true(buf);
    chk_true("MachineNumberQ[ArcCot[NDArray[{1.*^308, 0.5}]][[1]]] == False");
    /* exact zeros stay values and buffers stay packed (lane modes RANGE / NONE) */
    snprintf(buf, sizeof buf, "Developer`PackedArrayQ[Mod[%s[{2., 3., 4.}], 2.]]", pk);
    chk(buf, "True");
    snprintf(buf, sizeof buf, "Developer`PackedArrayQ[Log[2., %s[{1., 2., 4.}]]]", pk);
    chk(buf, "True");
    snprintf(buf, sizeof buf, "Developer`PackedArrayQ[BesselJ[1, %s[{0., 1., 2.}]]]", pk);
    chk(buf, "True");
    snprintf(buf, sizeof buf, "Developer`PackedArrayQ[LogGamma[%s[{1., 2., 3.}]]]", pk);
    chk(buf, "True");
    chk("Compile[{{x, _Real}}, Mod[x, 2.]][4.]", "0.0");
}

/* G4: one machine helper per elementary head -- a real argument takes the real
 * libm function on every surface (mth_machine_elem.h, mth_machine_log.h). */
static void test_g4_one_helper(void) {
    const char* tpl =
        "With[{v = {0.3, -0.7, 2.5, -3.1, 0.5, 17.25, -0.5}}, "
        "With[{s = %s /@ v}, {%s[Developer`ToPackedArray[v]] === s, Normal[%s[NDArray[v]]] === s, "
        "Compile[{{x, _Real, 1}}, %s[x]][v] === s, (Compile[{{x, _Real}}, %s[x]] /@ v) === s}]]";
    const char* heads[] = { "Sin", "Cos", "Tan", "Cot", "Sec", "Csc", "Sinh", "Cosh", "Tanh",
                            "Coth", "Sech", "Csch", "Exp", "ArcTan", "ArcCot", "ArcSinh",
                            "ArcCsch", "Log", "Log2", "Log10" };
    char buf[1024];
    for (size_t i = 0; i < sizeof heads / sizeof heads[0]; i++) {
        const char* h = heads[i];
        snprintf(buf, sizeof buf, tpl, h, h, h, h, h);
        chk(buf, "{True, True, True, True}");
    }
    chk("With[{w = {0.3, -0.7, 0.5, -0.5, 0.99, 1., -1.}}, With[{s = ArcSin /@ w}, "
        "{ArcSin[Developer`ToPackedArray[w]] === s, Compile[{{x, _Real, 1}}, ArcSin[x]][w] === s}]]",
        "{True, True}");
    chk_true("ArcTan[0.5] === Compile[{{x, _Real}}, ArcTan[x]][0.5]");
    chk_true("Binomial[5.5, 2.] == 12.375");
    chk_true("With[{v = {5.5, 7.25, 10.}}, Binomial[Developer`ToPackedArray[v], 2.] === (Binomial[#, 2.] & /@ v)]");
    chk_true("ExpIntegralE[7., Developer`ToPackedArray[{-1.5, 1.}]] === {ExpIntegralE[7., -1.5], ExpIntegralE[7., 1.]}");
    chk_true("Im[Compile[{{x, _Real}}, ExpIntegralE[7., x]][-1.5]] < 0");
}

/* G7 (IEEE signed zero on every surface), G8, G9, Quotient (exact floor). */
static void test_g789_decisions(void) {
    chk("ArcCot[-0.]", "-1.5708");
    chk("ArcCot[0.]", "1.5708");
    chk_true("ArcCot[Developer`ToPackedArray[{-0., 1.}]][[1]] === ArcCot[-0.]");
    chk_true("Compile[{{x, _Real}}, ArcCot[x]][-0.] === ArcCot[-0.]");
    chk("Arg[-0.]", "3.14159");
    chk("Arg[0.]", "0.0");
    chk("Arg[0]", "0");
    chk("Binomial[2., -1.*^308]", "0.0");                       /* G8 */
    chk("Log[1.*^308, 1.*^308]", "1.0");                         /* G9 */
    chk("Log[2, 2]", "1");
    chk("Log[x, x]", "1");
    chk("Quotient[0.5, 1.*^-10]", "4999999999");                 /* exact floor */
    chk_true("Quotient[Developer`ToPackedArray[{0.5, 7., -7.}], 1.*^-10] === Quotient[{0.5, 7., -7.}, 1.*^-10]");
    chk_true("Quotient[{0.5, 7., -7.}, 1.*^-10] === {4999999999, 69999999999, -70000000000}");
    chk_true("Compile[{{x, _Real}}, Quotient[x, 1.*^-10]][0.5] === 4999999999");
    chk_true("With[{m = 0.5, n = 1.*^-10}, SetPrecision[m, Infinity] == "
             "SetPrecision[n, Infinity] Quotient[m, n] + SetPrecision[Mod[m, n], Infinity]]");
    chk("{Quotient[7., 2.], Quotient[-7., 2.], Quotient[7., -2.]}", "{3, -4, -4}");
}

/* G10: undefined double -> int conversions and unreachable series caps. */
static void test_g10_huge_order(void) {
    chk_true("First[AbsoluteTiming[Quiet[StruveH[-1.*^308, 10.]]]] < 1");
    chk_true("First[AbsoluteTiming[Quiet[Compile[{{x, _Real}, {y, _Real}}, StruveH[x, y]][-1.*^308, 10.]]]] < 1");
    chk_true("First[AbsoluteTiming[Quiet[BesselK[1.*^308, 1.*^308]]]] < 5");
    chk_true("First[AbsoluteTiming[Quiet[HermiteH[-1.*^308, 10.]]]] < 1");
    chk_true("Quiet[HermiteH[-1.*^10, 10.]] =!= 0");             /* was an exact 0 */
    /* HermiteH for x > 0 by 2^n U(-n/2, 1/2, x^2): the 1F1 difference
     * cancelled (references: the integral representation for n < 0) */
    chk_true("Abs[HermiteH[-100., 10.]/1.312579565568286*^-138 - 1] < 1*^-13");
    chk_true("Abs[HermiteH[-1000., 10.]/2.572488113145133*^-1608 - 1] < 1*^-13");
    chk_true("Abs[HermiteH[-7/2 // N, 10.]/0.00002689700578763475 - 1] < 1*^-14");
    chk_true("Abs[HermiteH[2.5, 3.] - (2 3. HermiteH[1.5, 3.] - 2 1.5 HermiteH[0.5, 3.])] < 1*^-12");
}

static void test_d9_negative_literal_power(void) {
    /* -a^b is -(a^b) for every literal spelling of a. */
    chk("-2.`30^2 < 0", "True");
    chk("-2`10^2 < 0", "True");
    chk_true("-2.5*^3^2 == -6.25 10^6");
    chk("-2.^2", "-4.0");
    chk("FullForm[Hold[-2.`30^2]]", "Hold[Times[-1, Power[2.`30., 2]]]");
}

static void test_e4_quotient_mpfr(void) {
    chk("Quotient[N[10^30, 40] + 1/2, 1] - 10^30", "0");
    chk("Quotient[-N[7, 30], 2]", "-4");
    chk("Quotient[N[10, 30], 3, 1/2]", "3");
}

static void test_e2_machine_zero_times_mpfr(void) {
    chk("Exp[1000.] * 0.", "0.0");
    chk("MachineNumberQ[Exp[1000.] * 0.]", "True");
}

static void test_e10_tiny_iterator_step(void) {
    chk_true("Range[1.*^-200, 2.*^-200, 1.*^-200] == {10^-200, 2 10^-200}");
    chk("Length[Table[x, {x, 1.*^-200, 5.*^-200, 1.*^-200}]]", "5");
    chk("Length[Range[0., 1., 0.1]]", "11");                  /* unchanged */
    chk("Length[Table[x, {x, 0., 10., 0.01}]]", "1001");
}

/* ---- Phase 5 (v0.6596): fixed-precision correctness ---------------------- */

static void test_d3_accuracy_any_magnitude(void) {
    /* accuracy = precision - log10|x|, with log10 from the exponent (any range) */
    chk("Accuracy[N[10^-400, 20]]", "420.0");
    chk("Accuracy[N[10^500, 20]]", "-480.0");
    chk_true("Abs[Accuracy[1.*^400] - ($MachinePrecision - 400)] < 10^-6");
    chk("Accuracy[1.5]", "15.7785");                         /* unchanged */
}

static void test_d4_setaccuracy_magnitude(void) {
    chk_true("Precision[SetAccuracy[Pi 10^50, 10]] > 60");
    chk_true("Accuracy[SetAccuracy[Pi, 10]] >= 10");          /* never less than asked */
    chk_true("Accuracy[SetAccuracy[Pi 10^50, 10]] >= 10");
    chk("SetAccuracy[Pi/10^50, 10]", "0.*^-10");               /* no digit left: 0``10 */
    chk("Accuracy[SetAccuracy[Pi/10^50, 10]]", "10.0");
}

static void test_d8_precision_mark_is_decimal(void) {
    chk("FullForm[Hold[1.5`+nan]]", "Hold[Plus[1.5, nan]]");
    chk("Precision[1.5`30]", "30.0");
    chk("FullForm[Hold[2.5`10..]]", "Hold[Repeated[2.5`10.]]");
}

static void test_d7_nonpositive_precision(void) {
    chk("Quiet[N[1/3, -1]]", "N[1/3, -1]");
    chk("Check[N[1/3, 0], \"msg\"]", "\"msg\"");               /* the message fires */
    chk("Check[SetPrecision[1.5, 0], \"msg\"]", "\"msg\"");
    chk("N[1/3, 5]", "0.33333");
    chk("Block[{$MaxExtraPrecision = Infinity}, Precision[N[Sin[10^60], 20]]]", "20.0");
}

/* ---- D7b (v0.6597): $MinPrecision / $MaxPrecision ----------------------- */

static void test_d7b_min_max_precision(void) {
    chk("{$MinPrecision, $MaxPrecision}", "{0, Infinity}");
    /* requests are bounded (with a message) */
    chk("Quiet[Block[{$MaxPrecision = 20}, Precision[N[Pi, 50]]]]", "20.0");
    chk("Quiet[Block[{$MinPrecision = 30}, Precision[N[Pi, 10]]]]", "30.0");
    chk("Check[Block[{$MaxPrecision = 20}, N[Pi, 50]], \"msg\"]", "\"msg\"");
    chk("Quiet[Block[{$MaxPrecision = 20}, Precision[SetPrecision[1.5, 40]]]]", "20.0");
    /* every computed arbitrary-precision number obeys them */
    chk("Block[{$MaxPrecision = 20}, Precision[Sqrt[2.`30]]]", "20.0");
    chk("Block[{$MaxPrecision = 20}, Precision[1.5`40 + 0]]", "20.0");
    chk("Quiet[Block[{$MinPrecision = 25, $MaxPrecision = 25}, Precision[1.5`10 * 2.5`40]]]", "25.0");
    /* machine numbers are not arbitrary-precision numbers */
    chk("Block[{$MinPrecision = 30}, Precision[1.5 2.]]", "MachinePrecision");
    /* Block restores; assignment is validated */
    chk("Quiet[Block[{$MaxPrecision = 20}, N[Pi, 50]]]; Precision[N[Pi, 50]]", "50.0");
    chk("Quiet[$MaxPrecision = -3]; $MaxPrecision", "Infinity");
    chk("Quiet[$MaxPrecision = 30; $MinPrecision = 40]; {$MinPrecision, $MaxPrecision}", "{0, 30}");   /* min > max refused */
    chk("Quiet[$MinPrecision = 20; $MaxPrecision = 10]; {$MinPrecision, $MaxPrecision}", "{20, 30}");  /* max < min refused */
    chk("$MaxPrecision = Infinity; $MinPrecision = 0; {$MinPrecision, $MaxPrecision}", "{0, Infinity}");
    /* N's adaptive precision stops at the bound instead of spinning */
    chk("Block[{$MaxPrecision = 30}, First[AbsoluteTiming[N[Sin[10^60], 20]]] < 5]", "True");
}

/* ---- Phase 6: the policy items decided 2026-11-17 ------------------------ */

/* D5: below MPFR's exponent range the answer is Underflow[] with General::unfl,
 * never a zero that looks exact. A true zero stays a zero. */
static void test_d5_underflow(void) {
    chk("Quiet[Exp[-10.`30^10]]", "Underflow[]");
    chk("Quiet[Exp[-10^10 1.]]", "Underflow[]");                 /* machine path */
    chk("Quiet[1.`20*^-300000000 1.`20*^-300000000]", "Underflow[]");
    chk("Quiet[-1.`20*^-300000000 1.`20*^-300000000]", "-Underflow[]");
    chk("Check[Exp[-10.`30^10], \"msg\"]", "\"msg\"");               /* the message fires */
    chk_true("Quiet[Check[0. Exp[-10.], \"msg\"]] === 0.");      /* a true zero */
    chk_true("Quiet[Check[1.`30 - 1.`30, \"msg\"]] == 0");
    chk_true("Quiet[Check[Exp[-1000.], \"msg\"]] != 0");         /* in MPFR's range */
}

/* E1: a machine-precision value back inside the double range is a double. */
static void test_e1_demote(void) {
    chk_true("MachineNumberQ[(1.*^308 10.)/1.*^10]");
    chk_true("MachineNumberQ[(1.*^-300 1.*^-300) 1.*^300]");
    chk_true("MachineNumberQ[(1.*^308 10.) - 9.*^308]");
    chk_true("MachineNumberQ[Sqrt[1.*^400 1.]]");
    chk_true("MachineNumberQ[Exp[1000.] Exp[-999.]]");
    chk_true("MachineNumberQ[Re[Complex[Exp[1000.], Exp[1000.]] Exp[-999.]]] && Head[Complex[Exp[1000.], Exp[1000.]] Exp[-999.]] === Complex");
    chk_true("!MachineNumberQ[Im[Complex[Exp[1000.], 1.] Exp[-999.]]]");   /* 1.4*^-434 stays MPFR */
    chk_true("!MachineNumberQ[1.*^308 10.] && !MachineNumberQ[N[2, 30] 10^300]");
}

/* C9a / C9b: a machine real always shows a decimal point and *^ (it re-reads
 * as a Real); a 53-bit MPFR is machine precision and prints like one. */
static void test_c9_printing(void) {
    chk("1.*^15", "1.0*^15");
    chk("2.^-52", "2.22045*^-16");
    chk("123456789.", "1.23457*^8");
    chk("100000.", "100000.0");
    chk("-3.*^20", "-3.0*^20");
    chk("Exp[1000.]", "1.97007*^434");
    chk("InputForm[Exp[1000.]]", "1.970071114017047*^434");
    chk_true("ToExpression[ToString[InputForm[Exp[1000.]]]] === Exp[1000.]");
    chk_true("Head[ToExpression[ToString[1.*^15]]] === Real");
    chk("Precision[Exp[1000.]]", "MachinePrecision");
    chk("Precision[1.*^400]", "MachinePrecision");
    chk("Precision[Complex[Exp[1000.], 1.]]", "MachinePrecision");
    chk("Precision[1.`30*^400]", "30.0");
    chk("N[10^50, 20]", "1.0000000000000000000*^50");
    chk("TeXForm[1.*^15]", "1.\\times 10^{15}");
}

/* E3: past the exponent cap a Rational power stays symbolic; a negative
 * exponent moves to the reciprocal base, and a negative BigInt exponent prints
 * like a machine one (InputForm parentheses, OutputForm/TeX denominator). */
static void test_e3_rational_power_past_cap(void) {
    chk("InputForm[(2/3)^(-2^63)]", "(3/2)^9223372036854775808");
    chk("InputForm[(-2/3)^(-2^64 - 1)]", "(-3/2)^18446744073709551617");
    chk("InputForm[(2/3)^(2^64)]", "(2/3)^18446744073709551616");
    chk_true("(2/3)^(-2^63) === (3/2)^(2^63)");
    chk("InputForm[x^(-2^64)]", "x^(-18446744073709551616)");
    chk("ToString[x^(-2^64)]", "\"1/x^18446744073709551616\"");
    chk("ToString[TeXForm[x^(-2^64)]]", "\"\\\\frac{1}{x^{18446744073709551616}}\"");
    chk("InputForm[a x^(-2^64)]", "a*x^(-18446744073709551616)");
    chk_true("ToExpression[ToString[InputForm[x^(-2^64)]]] === x^(-2^64)");
}

/* F6 (found in phase 7): real Binomial through the Gamma quotient. A Gamma
 * factor past MPFR's range (Gamma[10^8 + 1]) gave Indeterminate for a small
 * binomial; integer-valued poles gave Indeterminate where the exact path has a
 * value. Expected values: the exact Binomial at the same integers. */
static void test_f6_binomial_real(void) {
    chk("InputForm[Binomial[1.*^8, 2.]]", "4.99999995*^15");
    chk("InputForm[Binomial[1.*^15, 2.]]", "4.999999999999995*^29");
    chk_true("Binomial[10^20, 2.] == N[Binomial[10^20, 2]]");
    chk("InputForm[Binomial[-5., 3.]]", "-35.0");
    chk("InputForm[Binomial[-1.*^3, 3.]]", "-167167000.0");
    chk("InputForm[Binomial[5, 7.]]", "0.0");
    chk("InputForm[Binomial[-5., -7.]]", "15.0");
    chk("InputForm[Binomial[-2., -5.]]", "-4.0");
    chk("InputForm[Binomial[-5., -2.]]", "0.0");
    chk("InputForm[Binomial[-5.5, 3.]]", "-44.6875");
    chk_true("Binomial[N[-5, 30], N[3, 30]] == -35");
    chk_true("Binomial[N[-5, 30], N[-7, 30]] == 15");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_a1_int64_min_rational_power);
    TEST(test_a2_quotient_huge_real);
    TEST(test_a5_power_past_gmp_range);
    TEST(test_a6_print_int64_min_exponent);
    TEST(test_b1_divide_cross_product_overflow);
    TEST(test_b2_gcd_lcm_int64_min);
    TEST(test_b3_negative_denominator);
    TEST(test_b4_rational_bignum_canonical);
    TEST(test_b5_bigint_hash_matches_integer);
    TEST(test_b6_bigint_exponent);
    TEST(test_c1_parser_degenerate_literals);
    TEST(test_c2_plus_times_overflow_underflow);
    TEST(test_c4_negative_base_machine_power);
    TEST(test_c6_machine_contagion);
    TEST(test_c7_compare_out_of_double_range);
    TEST(test_c8_no_subnormal_machine_numbers);
    TEST(test_c3_e5_buffer_paths_promote);
    TEST(test_e9_compile_underflow);
    TEST(test_e12_compile_unrepresentable_arg);
    TEST(test_g1_compile_integer_subtract);
    TEST(test_g2_inverse_jacobi_dn_packed);
    TEST(test_g3_machine_mod);
    TEST(test_g56_lane_checks);
    TEST(test_g4_one_helper);
    TEST(test_g789_decisions);
    TEST(test_g10_huge_order);
    TEST(test_d9_negative_literal_power);
    TEST(test_e4_quotient_mpfr);
    TEST(test_e2_machine_zero_times_mpfr);
    TEST(test_e10_tiny_iterator_step);
    TEST(test_d3_accuracy_any_magnitude);
    TEST(test_d4_setaccuracy_magnitude);
    TEST(test_d8_precision_mark_is_decimal);
    TEST(test_d7_nonpositive_precision);
    TEST(test_d7b_min_max_precision);
    TEST(test_d5_underflow);
    TEST(test_e1_demote);
    TEST(test_c9_printing);
    TEST(test_e3_rational_power_past_cap);
    TEST(test_f6_binomial_real);

    if (g_fail) {
        printf("%d numeric-tower review check(s) FAILED\n", g_fail);
        return 1;
    }
    printf("All numeric-tower review tests passed!\n");
    return 0;
}
