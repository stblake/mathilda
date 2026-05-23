/* Unit tests for precision-aware numeric literal parsing.
 *
 * The parser routes real literals to one of:
 *   - EXPR_REAL  (machine double) when the implied precision of the
 *     typed mantissa is ≤ MachinePrecision (≈ 15.95 decimal digits);
 *   - EXPR_MPFR  with ceil(implied_precision × log2(10)) bits when
 *     the implied precision exceeds MachinePrecision.
 *
 * Implied precision is computed from the mantissa only (digits + '.',
 * excluding any 'e' / 'E' / '*^' exponent suffix):
 *
 *     implied = frac_digits + log10(|mantissa_value|)
 *
 * This matches Mathematica's behavior, where '1.0e22' (mantissa "1.0",
 * one fractional digit, magnitude 1) is MachinePrecision but
 * '1.0000000000000000000000*^22' (23 fractional digits) is high-precision.
 *
 * Coverage map:
 *   A  basic machine-precision literals
 *   B  high-precision literals (auto-promoted to MPFR)
 *   C  threshold boundary (just below / just above MachinePrecision)
 *   D  user-provided Mathematica examples (13.463 / 54.463 / 71.463)
 *   E  scientific notation (e/E)
 *   F  scaled scientific notation (*^)
 *   G  zero literals (always machine-precision)
 *   H  no-decimal reals (`1e10` etc.) stay machine
 *   I  pure integers untouched (Integer / BigInt)
 *   J  explicit backtick suffix still overrides (backwards-compat)
 *   K  Precision[] / Accuracy[] integration on user examples
 *   L  negative literals on both paths
 *   M  mantissa value preserved (round-trip through MPFR)
 */

#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* ------------------------------------------------------------------------
 *  Type-level assertions: parse, check the resulting Expr->type tag,
 *  then free. Frees on assertion failure are skipped since ASSERT exits.
 * ---------------------------------------------------------------------- */

static void assert_parses_as_real(const char* src) {
    Expr* p = parse_expression(src);
    ASSERT_MSG(p != NULL, "parse failed for: %s", src);
    ASSERT_MSG(p->type == EXPR_REAL,
        "expected EXPR_REAL (machine double), got type %d for: %s",
        (int)p->type, src);
    expr_free(p);
}

static void assert_parses_as_integer(const char* src) {
    Expr* p = parse_expression(src);
    ASSERT_MSG(p != NULL, "parse failed for: %s", src);
    ASSERT_MSG(p->type == EXPR_INTEGER,
        "expected EXPR_INTEGER, got type %d for: %s", (int)p->type, src);
    expr_free(p);
}

static void assert_parses_as_bigint(const char* src) {
    Expr* p = parse_expression(src);
    ASSERT_MSG(p != NULL, "parse failed for: %s", src);
    ASSERT_MSG(p->type == EXPR_BIGINT,
        "expected EXPR_BIGINT, got type %d for: %s", (int)p->type, src);
    expr_free(p);
}

#ifdef USE_MPFR
static void assert_parses_as_mpfr(const char* src) {
    Expr* p = parse_expression(src);
    ASSERT_MSG(p != NULL, "parse failed for: %s", src);
    ASSERT_MSG(p->type == EXPR_MPFR,
        "expected EXPR_MPFR (arbitrary-precision), got type %d for: %s",
        (int)p->type, src);
    expr_free(p);
}

/* Verify an MPFR literal has approximately `expected_bits` bits of
 * precision.  Tolerance handles small bit-quantization slop near the
 * boundary, but we keep it tight (≤ 2 bits) so accidental mis-binning
 * is still detected. */
static void assert_parses_as_mpfr_bits(const char* src,
                                       long expected_bits, long tol) {
    Expr* p = parse_expression(src);
    ASSERT_MSG(p != NULL, "parse failed for: %s", src);
    ASSERT_MSG(p->type == EXPR_MPFR,
        "expected EXPR_MPFR, got type %d for: %s", (int)p->type, src);
    long got = (long)mpfr_get_prec(p->data.mpfr);
    long delta = got - expected_bits;
    if (delta < 0) delta = -delta;
    ASSERT_MSG(delta <= tol,
        "for %s: expected ~%ld bits (tol %ld), got %ld",
        src, expected_bits, tol, got);
    expr_free(p);
}

/* Round-trip: MPFR literal must hold the same value as mpfr_set_str
 * at the parser's chosen precision — verifying the parser feeds the
 * literal text to MPFR directly (no lossy double round-trip). */
static void assert_mpfr_value_matches_str(const char* src) {
    Expr* p = parse_expression(src);
    ASSERT_MSG(p != NULL && p->type == EXPR_MPFR,
        "parse failed or wrong type for: %s", src);
    mpfr_prec_t bits = mpfr_get_prec(p->data.mpfr);
    mpfr_t ref;
    mpfr_init2(ref, bits);
    mpfr_set_str(ref, src, 10, MPFR_RNDN);
    int cmp = mpfr_cmp(p->data.mpfr, ref);
    ASSERT_MSG(cmp == 0,
        "MPFR value mismatch for %s (parsed ≠ direct mpfr_set_str @ %ld bits)",
        src, (long)bits);
    mpfr_clear(ref);
    expr_free(p);
}
#endif

/* ------------------------------------------------------------------------
 *  A.  Basic machine-precision literals
 * ---------------------------------------------------------------------- */

static void test_basic_machine_precision(void) {
    assert_parses_as_real("1.5");
    assert_parses_as_real("3.14");
    assert_parses_as_real("0.5");
    assert_parses_as_real("0.1");
    assert_parses_as_real("123.456");
    assert_parses_as_real("3.14159265");          /*  8 sig digits */
    assert_parses_as_real("3.141592653589");      /* 12 sig digits */
    assert_parses_as_real("0.123456789012345");   /* 15 frac digits, impl≈14.09 */
    assert_parses_as_real("123.");                /* trailing dot, no frac */
    assert_parses_as_real(".5");                  /* leading dot, no int */
}

static void test_basic_machine_precision_negative(void) {
    assert_parses_as_real("-1.5");
    assert_parses_as_real("-3.14");
    assert_parses_as_real("-0.0");
    assert_parses_as_real("-123.456");
    assert_parses_as_real("-29037945.290347");    /* user example 1, impl≈13.46 */
}

/* ------------------------------------------------------------------------
 *  B.  High-precision literals (auto-promoted to MPFR)
 * ---------------------------------------------------------------------- */

#ifdef USE_MPFR
static void test_promotes_to_mpfr_basic(void) {
    /* 16-digit fraction of 1.x — implied prec ≈ 16 → MPFR. */
    assert_parses_as_mpfr("1.2345678901234567");
    /* 18 fractional digits — implied ≈ 18.09 → MPFR. */
    assert_parses_as_mpfr("3.141592653589793238");
    /* 25-digit literal — implied ≈ 24.09 → MPFR. */
    assert_parses_as_mpfr("1.234567890123456789012345");
    /* 50-digit fraction → MPFR. */
    assert_parses_as_mpfr(
        "1.12345678901234567890123456789012345678901234567890");
}

static void test_promotes_to_mpfr_negative(void) {
    assert_parses_as_mpfr("-1.2345678901234567");
    assert_parses_as_mpfr("-3.141592653589793238462643383279");
}

static void test_promotes_to_mpfr_large_magnitude(void) {
    /* Large integer part + many decimals: the user's 54-digit example. */
    assert_parses_as_mpfr(
        "29037852093587905730945.29034875093457832094573984537498");
    assert_parses_as_mpfr(
        "-29037852093587905730945.29034875093457832094573984537498");
}

static void test_promotes_to_mpfr_small_magnitude(void) {
    /* Many leading zeros after the decimal — log10 is very negative.
     * The leading zeros do NOT count as precision: 3 sig digits at
     * implied prec ≈ 3 stays machine. We need lots of trailing
     * non-zero digits to cross the threshold. */
    assert_parses_as_real("0.0000000000000000000000000000000123");        /* 3 sig */
    assert_parses_as_mpfr(
        "0.000000000000000000000000000000012345678901234567890123456");   /* 26 sig */
}
#endif

/* ------------------------------------------------------------------------
 *  C.  Threshold boundary (MachinePrecision ≈ 15.95 digits)
 * ---------------------------------------------------------------------- */

#ifdef USE_MPFR
static void test_threshold_boundary(void) {
    /* implied ≈ 14.30 (15 frac digits, leading "1") → REAL. */
    assert_parses_as_real("1.23456789012345");
    /* implied ≈ 15.30 → REAL (just under). */
    assert_parses_as_real("1.234567890123456");
    /* implied ≈ 16.30 → MPFR (just over). */
    assert_parses_as_mpfr("1.2345678901234567");
    /* Borderline near magnitude ~10: implied ≈ 13.95 → REAL. */
    assert_parses_as_real("9.99999999999998");        /* 14 frac, log10≈1 */
    /* Mantissa 9.99... with 15 frac digits → implied ≈ 15.95-ish → REAL. */
    assert_parses_as_real("9.99999999999999");
    /* 16 frac digits with leading 9 → implied ≈ 16.95 → MPFR. */
    assert_parses_as_mpfr("9.999999999999999");
}
#endif

/* ------------------------------------------------------------------------
 *  D.  User-provided Mathematica examples
 *      Expected implied precisions: 13.463, 54.463, 71.463.
 *      Corresponding MPFR bits: 0 (machine), 181, 238.
 * ---------------------------------------------------------------------- */

static void test_user_example_1_machine(void) {
    /* implied = 6 + log10(2.9037945e7) ≈ 13.463 → MachinePrecision. */
    assert_parses_as_real("-29037945.290347");
    assert_eval_eq("Precision[-29037945.290347]", "MachinePrecision", 0);
}

#ifdef USE_MPFR
static void test_user_example_2_mpfr_54(void) {
    /* implied = 32 + log10(2.9037e22) ≈ 54.463 → ceil(54.463*log2(10))
     * = ceil(180.93) = 181 bits. */
    const char* src =
        "-29037852093587905730945.29034875093457832094573984537498";
    assert_parses_as_mpfr(src);
    assert_parses_as_mpfr_bits(src, 181, 1);
    /* Precision = 181 / log2(10) ≈ 54.487 → print starts with "54.". */
    assert_eval_startswith(
        "Precision[-29037852093587905730945.29034875093457832094573984537498]",
        "54.");
    /* Accuracy = Precision − log10|x| ≈ 54.487 − 22.463 = 32.024 → "32.". */
    assert_eval_startswith(
        "Accuracy[-29037852093587905730945.29034875093457832094573984537498]",
        "32.");
}

static void test_user_example_3_mpfr_71(void) {
    /* implied = 49 + log10(2.9037e22) ≈ 71.463 → ceil(71.463*log2(10))
     * = ceil(237.39) = 238 bits. */
    const char* src =
        "-29037852093587905730945.2903487509345783209457398453749803489530945034950";
    assert_parses_as_mpfr(src);
    assert_parses_as_mpfr_bits(src, 238, 1);
    assert_eval_startswith(
        "Precision[-29037852093587905730945."
        "2903487509345783209457398453749803489530945034950]",
        "71.");
    assert_eval_startswith(
        "Accuracy[-29037852093587905730945."
        "2903487509345783209457398453749803489530945034950]",
        "49.");
}
#endif

/* ------------------------------------------------------------------------
 *  E.  Scientific notation (e/E)
 *      Precision is determined by the MANTISSA only, NOT the final value.
 * ---------------------------------------------------------------------- */

static void test_scientific_low_precision_mantissa(void) {
    /* '1.0e22' has mantissa "1.0", implied precision 1 → REAL. */
    assert_parses_as_real("1.0e22");
    assert_parses_as_real("1.0E22");
    assert_parses_as_real("2.5e-3");
    assert_parses_as_real("-1.5e10");
    assert_parses_as_real("3.14e100");           /* 3 sig digits */
    /* The huge final magnitude must NOT trigger promotion. */
    assert_eval_eq("Precision[1.0e22]", "MachinePrecision", 0);
}

#ifdef USE_MPFR
static void test_scientific_high_precision_mantissa(void) {
    /* Many mantissa digits → MPFR despite exponent. */
    assert_parses_as_mpfr("1.123456789012345678e22");
    assert_parses_as_mpfr("1.12345678901234567890123456789e-5");
    /* Negative + high mantissa precision. */
    assert_parses_as_mpfr("-1.123456789012345678e22");
}
#endif

/* ------------------------------------------------------------------------
 *  F.  Scaled scientific notation (*^)
 *      The mantissa-only rule also applies to `*^` literals.
 * ---------------------------------------------------------------------- */

static void test_scaled_low_precision(void) {
    assert_eval_eq("Precision[1.0*^22]", "MachinePrecision", 0);
    /* 1.0*^22 evaluates to 1.×10^22 — must be machine precision. */
    Expr* p = parse_expression("1.0*^22");
    ASSERT_MSG(p != NULL, "parse failed");
    ASSERT_MSG(p->type == EXPR_REAL,
        "expected EXPR_REAL for 1.0*^22, got type %d", (int)p->type);
    expr_free(p);
}

#ifdef USE_MPFR
static void test_scaled_high_precision(void) {
    Expr* p = parse_expression("1.123456789012345678*^22");
    ASSERT_MSG(p != NULL, "parse failed");
    ASSERT_MSG(p->type == EXPR_MPFR,
        "expected EXPR_MPFR for 1.123456789012345678*^22, got type %d",
        (int)p->type);
    expr_free(p);
    /* And via Precision[]. */
    assert_eval_startswith("Precision[1.123456789012345678*^22]", "18.");
}
#endif

/* ------------------------------------------------------------------------
 *  G.  Zero literals  →  always machine precision.
 * ---------------------------------------------------------------------- */

static void test_zero_literals(void) {
    assert_parses_as_real("0.0");
    assert_parses_as_real("-0.0");
    assert_parses_as_real("0.000");
    assert_parses_as_real("0.0000000000000000000000000");
    assert_parses_as_real("0.0e10");
    /* Zero with many decimals must NOT become MPFR (no nonzero mantissa). */
    assert_eval_eq("Precision[0.000]", "MachinePrecision", 0);
    assert_eval_eq("Precision[0.0000000000000000000000]", "MachinePrecision", 0);
}

/* ------------------------------------------------------------------------
 *  H.  No-decimal reals (`1e10`, `1E5` etc.) stay machine.
 * ---------------------------------------------------------------------- */

static void test_no_decimal_reals(void) {
    /* Without a '.', precision computation is skipped — always REAL.
     * (frac_digits = 0, has_dot = false → never promote.) */
    assert_parses_as_real("1e10");
    assert_parses_as_real("1E10");
    assert_parses_as_real("5e-3");
    assert_parses_as_real("123e7");
}

/* ------------------------------------------------------------------------
 *  I.  Pure integers untouched: parser still produces Integer / BigInt.
 * ---------------------------------------------------------------------- */

static void test_integer_literals_untouched(void) {
    assert_parses_as_integer("0");
    assert_parses_as_integer("1");
    assert_parses_as_integer("12345");
    assert_parses_as_integer("-7");
    assert_parses_as_integer("9223372036854775806");   /* fits int64 */
    /* 20+ digits → BigInt promotion (existing behavior). */
    assert_parses_as_bigint("99999999999999999999");          /* 20 digits */
    assert_parses_as_bigint("100000000000000000000");
    assert_parses_as_bigint(
        "12345678901234567890123456789012345678901234567890");
}

/* ------------------------------------------------------------------------
 *  J.  Explicit backtick suffix still overrides implicit precision.
 * ---------------------------------------------------------------------- */

#ifdef USE_MPFR
static void test_backtick_suffix_still_works(void) {
    /* `3.14`50` → 50-digit MPFR (unchanged). */
    assert_eval_startswith("Precision[3.14`50]", "50.");
    /* `3.14``49` (double backtick) → 49-digit accuracy (unchanged). */
    assert_eval_startswith("Accuracy[3.14``49]", "49.");
    /* Bare backtick on integer → machine-precision Real. */
    Expr* p = parse_expression("3`");
    ASSERT_MSG(p != NULL, "parse failed for 3`");
    ASSERT_MSG(p->type == EXPR_REAL,
        "expected EXPR_REAL for 3`, got type %d", (int)p->type);
    expr_free(p);
    /* High implicit precision + explicit suffix — suffix wins. */
    assert_eval_startswith(
        "Precision[1.123456789012345678901234567890`20]", "20.");
}
#endif

/* ------------------------------------------------------------------------
 *  K.  Precision[] / Accuracy[] integration on user examples.
 *      Verifies the full end-to-end pipeline (parser → printer).
 * ---------------------------------------------------------------------- */

static void test_precision_machine_examples(void) {
    /* All machine-precision literals report MachinePrecision symbol. */
    assert_eval_eq("Precision[1.5]",     "MachinePrecision", 0);
    assert_eval_eq("Precision[3.14]",    "MachinePrecision", 0);
    assert_eval_eq("Precision[1.0e22]",  "MachinePrecision", 0);
    assert_eval_eq("Precision[-29037945.290347]", "MachinePrecision", 0);
    /* Integer / Rational accuracy = Infinity. */
    assert_eval_eq("Precision[123]",     "Infinity", 0);
    assert_eval_eq("Precision[1/3]",     "Infinity", 0);
}

#ifdef USE_MPFR
static void test_precision_mpfr_examples(void) {
    /* 24-digit mantissa → precision ≈ 24.09. */
    assert_eval_startswith(
        "Precision[1.234567890123456789012345]", "24.");
    /* Negative version same precision. */
    assert_eval_startswith(
        "Precision[-1.234567890123456789012345]", "24.");
    /* Accuracy = Precision − log10(|x|) ≈ 24.09 − 0.09 = 24.00. */
    assert_eval_startswith(
        "Accuracy[1.234567890123456789012345]", "24.");
}
#endif

/* ------------------------------------------------------------------------
 *  L.  Mantissa value preserved through MPFR round-trip.
 * ---------------------------------------------------------------------- */

#ifdef USE_MPFR
static void test_mpfr_value_roundtrip(void) {
    /* The parser must feed the original literal string to mpfr_set_str
     * at the chosen precision — not double-strtod first then convert. */
    assert_mpfr_value_matches_str("3.141592653589793238462643");
    assert_mpfr_value_matches_str(
        "29037852093587905730945.29034875093457832094573984537498");
    assert_mpfr_value_matches_str(
        "-29037852093587905730945."
        "2903487509345783209457398453749803489530945034950");
}
#endif

/* ------------------------------------------------------------------------
 *  M.  Frac-only literals (.5, .05) and very long fractional literals.
 * ---------------------------------------------------------------------- */

#ifdef USE_MPFR
static void test_frac_only_high_precision(void) {
    /* No integer part, many decimals. .12345...30digits has implied
     * precision ≈ 30 + log10(0.1234) = 30 + (-0.91) = 29.09 → MPFR. */
    assert_parses_as_mpfr(
        ".123456789012345678901234567890");
    /* The value should match mpfr_set_str at ~97 bits. */
    Expr* p = parse_expression(".123456789012345678901234567890");
    ASSERT_MSG(p != NULL && p->type == EXPR_MPFR, "parse failed");
    long bits = (long)mpfr_get_prec(p->data.mpfr);
    ASSERT_MSG(bits >= 95 && bits <= 100,
        "expected ~97 bits for 30-digit fraction, got %ld", bits);
    expr_free(p);
}

static void test_very_long_literal(void) {
    /* 100-digit fraction, integer part 1 → implied ≈ 100. */
    const char* src =
        "1."
        "1234567890123456789012345678901234567890"
        "1234567890123456789012345678901234567890"
        "12345678901234567890";
    assert_parses_as_mpfr(src);
    /* Implied precision ≈ 100.0 → bits = ceil(100 × log2(10)) = 333. */
    assert_parses_as_mpfr_bits(src, 333, 2);
}
#endif

/* ------------------------------------------------------------------------
 *  main()
 * ---------------------------------------------------------------------- */

int main(void) {
    symtab_init();
    core_init();

    /* A — basic machine precision */
    TEST(test_basic_machine_precision);
    TEST(test_basic_machine_precision_negative);

    /* E — scientific notation (no-decimal-promotion semantics) */
    TEST(test_scientific_low_precision_mantissa);

    /* F — scaled scientific (*^) */
    TEST(test_scaled_low_precision);

    /* G — zero literals */
    TEST(test_zero_literals);

    /* H — no-decimal reals */
    TEST(test_no_decimal_reals);

    /* I — pure integers */
    TEST(test_integer_literals_untouched);

    /* K — Precision[] on machine examples (works with USE_MPFR=0 too) */
    TEST(test_precision_machine_examples);

    /* D part 1 — user example 1 (machine path) */
    TEST(test_user_example_1_machine);

#ifdef USE_MPFR
    /* B — auto-promotion */
    TEST(test_promotes_to_mpfr_basic);
    TEST(test_promotes_to_mpfr_negative);
    TEST(test_promotes_to_mpfr_large_magnitude);
    TEST(test_promotes_to_mpfr_small_magnitude);

    /* C — threshold boundary */
    TEST(test_threshold_boundary);

    /* D part 2/3 — user examples 2 & 3 */
    TEST(test_user_example_2_mpfr_54);
    TEST(test_user_example_3_mpfr_71);

    /* E — high-precision mantissas under scientific notation */
    TEST(test_scientific_high_precision_mantissa);

    /* F — scaled scientific with high mantissa precision */
    TEST(test_scaled_high_precision);

    /* J — backtick suffix unchanged */
    TEST(test_backtick_suffix_still_works);

    /* K — MPFR Precision[] integration */
    TEST(test_precision_mpfr_examples);

    /* L — round-trip MPFR value preservation */
    TEST(test_mpfr_value_roundtrip);

    /* M — frac-only & 100-digit literals */
    TEST(test_frac_only_high_precision);
    TEST(test_very_long_literal);
#endif

    printf("All parser_precision_tests passed.\n");
    return 0;
}
