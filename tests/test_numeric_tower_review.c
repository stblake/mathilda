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
    TEST(test_d9_negative_literal_power);
    TEST(test_e4_quotient_mpfr);
    TEST(test_e2_machine_zero_times_mpfr);
    TEST(test_e10_tiny_iterator_step);

    if (g_fail) {
        printf("%d numeric-tower review check(s) FAILED\n", g_fail);
        return 1;
    }
    printf("All numeric-tower review tests passed!\n");
    return 0;
}
