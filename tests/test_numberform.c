/*
 * Tests for NumberForm (and the minimal Row that its NumberFormat option uses).
 *
 * NumberForm is an inert PRINT WRAPPER: it stays in the expression tree and
 * only changes how the wrapped expression is rendered. Every test therefore
 * evaluates the NumberForm[...] call (which leaves the head in place) and then
 * checks the printed form via expr_to_string — exactly the REPL path.
 *
 * Assertions use the hard ASSERT_STR_EQ (exit(1) on failure), which survives a
 * Release/-DNDEBUG build, rather than the soft assert_eval_eq helper.
 *
 * Note on list spacing: Mathilda's list printer always separates elements with
 * ", " (comma + space). The Mathematica reference prints padded lists more
 * tightly; the per-number padding here is identical, and the visible extra
 * space in a padded list is purely that standard separator.
 */

#include "core.h"
#include "expr.h"
#include "symtab.h"
#include "eval.h"
#include "test_utils.h"
#include "print.h"
#include "parse.h"
#include <stdio.h>
#include <string.h>

/* Evaluate `input` and assert its printed (standard-form) output. */
static void check(const char* input, const char* expected) {
    Expr* p = parse_expression(input);
    ASSERT(p != NULL);
    Expr* e = evaluate(p);
    expr_free(p);
    char* s = expr_to_string(e);
    ASSERT(s != NULL);
    ASSERT_STR_EQ(s, expected);
    free(s);
    expr_free(e);
}

/* Evaluate `input` and assert its FullForm (raw-tree) output. */
static void check_fullform(const char* input, const char* expected) {
    Expr* p = parse_expression(input);
    ASSERT(p != NULL);
    Expr* e = evaluate(p);
    expr_free(p);
    char* s = expr_to_string_fullform(e);
    ASSERT(s != NULL);
    ASSERT_STR_EQ(s, expected);
    free(s);
    expr_free(e);
}

/* ---- precision spec: n, {n,f}, default ---- */
void test_precision(void) {
    check("NumberForm[N[Pi],10]", "3.141592654");
    check("NumberForm[1.23456,{3,4}]", "1.2300");
    check("NumberForm[N[E],10]", "2.718281828");
    check("NumberForm[N[E],3]", "2.72");
    check("NumberForm[1.23,2]", "1.2");
    /* default machine display is 6 significant digits */
    check("NumberForm[12345.6789]", "12345.7");
    check("NumberForm[N[E]]", "2.71828");
}

/* ---- the {n,f} pi table (double rounding: n sig figs, then f decimals) ---- */
void test_nf_table(void) {
    check("NumberForm[N[Pi],{1,1}]", "3.0");
    check("NumberForm[N[Pi],{1,5}]", "3.00000");
    check("NumberForm[N[Pi],{3,2}]", "3.14");
    check("NumberForm[N[Pi],{3,5}]", "3.14000");
    check("NumberForm[N[Pi],{4,2}]", "3.14");
    check("NumberForm[N[Pi],{4,5}]", "3.14200");
    check("NumberForm[N[Pi],{5,1}]", "3.1");
    check("NumberForm[N[Pi],{5,3}]", "3.142");
    check("NumberForm[N[Pi],{5,4}]", "3.1416");
    check("NumberForm[N[Pi],{5,5}]", "3.14160");
}

/* ---- scientific vs decimal + ScientificNotationThreshold ---- */
void test_scientific(void) {
    check("NumberForm[1234567890.,10]", "1.23456789*10^(9)");
    check("NumberForm[1.25*10^8,10]", "1.25*10^(8)");
    check("NumberForm[1.25*10^8,10,NumberMultiplier->\"*\"]", "1.25*10^(8)");
    check("NumberForm[{-0.00001234,2.468032,5234452.134},7]",
          "{-0.00001234, 2.468032, 5.234452*10^(6)}");
    check("NumberForm[{-0.00001234,2.468032,5234452.134},7,ScientificNotationThreshold->{-4,8}]",
          "{-1.234*10^(-5), 2.468032, 5234452.}");
}

/* ---- DigitBlock / NumberSeparator, incl. big integers ---- */
void test_digitblock(void) {
    check("NumberForm[10^9,10]", "1000000000");
    check("NumberForm[10^9,DigitBlock->3]", "1,000,000,000");
    check("NumberForm[30!,DigitBlock->3]",
          "265,252,859,812,191,058,636,308,480,000,000");
    check("NumberForm[30!,DigitBlock->3,NumberSeparator->\" \"]",
          "265 252 859 812 191 058 636 308 480 000 000");
    check("NumberForm[30!,DigitBlock->5,NumberSeparator->\" \"]",
          "265 25285 98121 91058 63630 84800 00000");
}

/* ---- NumberPoint ---- */
void test_numberpoint(void) {
    check("NumberForm[1.2345,3]", "1.23");
    check("NumberForm[1.2345,3,NumberPoint->\",\"]", "1,23");
}

/* ---- NumberSigns ---- */
void test_numbersigns(void) {
    check("NumberForm[{-1/3.,2/3.},5]", "{-0.33333, 0.66667}");
    check("NumberForm[{-1/3.,2/3.},5,NumberSigns->{\"-\",\"+\"}]",
          "{-0.33333, +0.66667}");
    check("NumberForm[{-1/3.,2/3.},5,NumberSigns->{\"minus \",\"plus \"}]",
          "{minus 0.33333, plus 0.66667}");
}

/* ---- ExponentStep / ExponentFunction (incl. a Null return) ---- */
void test_exponent_options(void) {
    check("NumberForm[1234567890.,10,ExponentStep->6]", "1234.56789*10^(6)");
    check("NumberForm[N[E^Range[10,50,10]],ExponentFunction->(3Quotient[#,3]&)]",
          "{22.0265*10^(3), 485.165*10^(6), 10.6865*10^(12), 235.385*10^(15), 5.18471*10^(21)}");
    check("NumberForm[N[E^Range[10,50,10]],ExponentFunction->(If[-10<#<10,Null,#]&)]",
          "{22026.5, 485165195., 1.06865*10^(13), 2.35385*10^(17), 5.18471*10^(21)}");
}

/* ---- NumberFormat (a string result, a mantissa-only result, a Row) ---- */
void test_numberformat(void) {
    check("NumberForm[{8.^5,11.^7,13.^9},NumberFormat->(Row[{#1,\"e\",#3}]&)]",
          "{32768.e, 1.94872e7, 1.06045e10}");
    check("NumberForm[{8.^5,11.^7,13.^9},NumberFormat->(#1&)]",
          "{32768., 1.94872, 1.06045}");
    check("NumberForm[{8.^5,11.^7,13.^9},NumberFormat->(#3&),ExponentFunction->(#&)]",
          "{4, 7, 10}");
}

/* ---- Row on its own ---- */
void test_row(void) {
    check("Row[{1,2,3}]", "123");
    check("Row[{1,2,3},\"-\"]", "1-2-3");
    check("Row[{\"a\",\"b\"}]", "ab");
}

/* ---- NumberMultiplier ---- */
void test_multiplier(void) {
    check("NumberForm[1.25*10^8,10,NumberMultiplier->\" x \"]", "1.25 x 10^(8)");
}

/* ---- DefaultPrintPrecision ---- */
void test_default_precision(void) {
    check("NumberForm[12345.6789]", "12345.7");
    check("NumberForm[12345.6789,DefaultPrintPrecision->8]", "12345.679");
}

/* ---- NumberPadding + SignPadding (Mathilda list spacing, see header) ---- */
void test_padding(void) {
    /* default: no padding */
    check("NumberForm[{-6.7,6.888,6.99999},4]", "{-6.7, 6.888, 7.}");
    check("NumberForm[{-6.7,6.888,6.99999},4,NumberPadding->{\" \",\"\"}]",
          "{  -6.7,  6.888,     7.}");
    check("NumberForm[{-6.7,6.888,6.99999},{4,3},NumberPadding->{\"\",\"0\"}]",
          "{-6.700, 6.888, 7.000}");
    check("NumberForm[{-1.2345,2.4680},{5,2},NumberPadding->{\" \",\" \"}]",
          "{  -1.23,    2.47}");
    check("NumberForm[{-1.2345,2.4680},{5,2},SignPadding->True,NumberPadding->{\" \",\" \"}]",
          "{-  1.23,    2.47}");
}

/* ---- reqsigz: fewer requested digits than integer digits ---- */
void test_reqsigz(void) {
    check("NumberForm[12345.6,3]", "12300.");
}

/* ---- high-precision (MPFR) re-rounding ---- */
void test_high_precision(void) {
    check("NumberForm[N[Pi,50],20]", "3.1415926535897932385");
    check("NumberForm[N[Pi,50],10]", "3.141592654");
    check("NumberForm[N[Pi,50],5]", "3.1416");
}

/* ---- works over lists, matrices, and mixed symbolic expressions ---- */
void test_structures(void) {
    check("NumberForm[{5.73141,6.93729,3.11996,1.98576,9.71174},4]",
          "{5.731, 6.937, 3.12, 1.986, 9.712}");
    check("NumberForm[{{3.83975,7.52103,0.67226},{9.16179,5.23369,4.95147},"
          "{3.53205,1.58439,1.16991}},2]",
          "{{3.8, 7.5, 0.67}, {9.2, 5.2, 5.}, {3.5, 1.6, 1.2}}");
    check("NumberForm[1/3.*Sin[x/7.],3]", "0.333 Sin[0.143 x]");
}

/* ---- zero / negative zero ---- */
void test_zero(void) {
    check("NumberForm[0.,3]", "0.");
    check("NumberForm[-0.,3]", "0.");
    check("NumberForm[0.,{3,2}]", "0.00");
}

/* ---- print-wrapper semantics: the head survives; FullForm shows it raw ---- */
void test_wrapper_semantics(void) {
    check_fullform("NumberForm[1.23,2]", "NumberForm[1.23, 2]");
    /* an intervening NumberForm blocks arithmetic (stays inert) */
    check("10*NumberForm[1.23,2]", "10 1.2");
    /* the first argument evaluates (NHoldRest, not HoldAll) */
    check("NumberForm[1 + 1.0,3]", "2.");
}

/* ---- SetOptions changes the defaults ---- */
void test_set_options(void) {
    check("NumberForm[12345.6789]", "12345.7");
    Expr* p = parse_expression("SetOptions[NumberForm, DefaultPrintPrecision->8]");
    Expr* e = evaluate(p); expr_free(p); expr_free(e);
    check("NumberForm[12345.6789]", "12345.679");
    /* restore so test ordering cannot leak */
    p = parse_expression("SetOptions[NumberForm, DefaultPrintPrecision->Automatic]");
    e = evaluate(p); expr_free(p); expr_free(e);
    check("NumberForm[12345.6789]", "12345.7");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_precision);
    TEST(test_nf_table);
    TEST(test_scientific);
    TEST(test_digitblock);
    TEST(test_numberpoint);
    TEST(test_numbersigns);
    TEST(test_exponent_options);
    TEST(test_numberformat);
    TEST(test_row);
    TEST(test_multiplier);
    TEST(test_default_precision);
    TEST(test_padding);
    TEST(test_reqsigz);
    TEST(test_high_precision);
    TEST(test_structures);
    TEST(test_zero);
    TEST(test_wrapper_semantics);
    TEST(test_set_options);

    printf("All NumberForm tests passed!\n");
    return 0;
}
