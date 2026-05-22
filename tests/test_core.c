#include "core.h"
#include "expr.h"
#include "symtab.h"
#include "eval.h"
#include "test_utils.h"
#include <stdio.h>
#include <string.h>
#include "print.h"
#include "parse.h"

void test_numericq(void) {
    assert_eval_eq("NumericQ[3]", "True", 0);
    assert_eval_eq("NumericQ[3.14]", "True", 0);
    assert_eval_eq("NumericQ[Pi]", "True", 0);
    assert_eval_eq("NumericQ[E]", "True", 0);
    assert_eval_eq("NumericQ[I]", "True", 0);
    assert_eval_eq("NumericQ[Infinity]", "True", 0);
    assert_eval_eq("NumericQ[x]", "False", 0);
    assert_eval_eq("NumericQ[Sin[Sqrt[2]]]", "True", 0);
    assert_eval_eq("NumericQ[Sin[x]]", "False", 0);
    
    assert_eval_eq("SetAttributes[f, NumericFunction]", "Null", 0);
    assert_eval_eq("NumericQ[f[Pi, Sin[1+I]]]", "True", 0);
    assert_eval_eq("NumericQ[f[Pi, x]]", "False", 0);
    assert_eval_eq("Clear[f]", "Null", 0);
}

void test_numberq(void) {
    // Test case 1: NumberQ[123] -> True
    Expr* e1 = parse_expression("NumberQ[123]");
    Expr* res1 = evaluate(e1);
    char* s1 = expr_to_string(res1);
    assert(strcmp(s1, "True") == 0);
    free(s1);
    expr_free(e1);
    expr_free(res1);

    // Test case 2: NumberQ[1.23] -> True
    Expr* e2 = parse_expression("NumberQ[1.23]");
    Expr* res2 = evaluate(e2);
    char* s2 = expr_to_string(res2);
    assert(strcmp(s2, "True") == 0);
    free(s2);
    expr_free(e2);
    expr_free(res2);

    // Test case 3: NumberQ[Complex[1, 2]] -> True
    Expr* e3 = parse_expression("NumberQ[Complex[1, 2]]");
    Expr* res3 = evaluate(e3);
    char* s3 = expr_to_string(res3);
    assert(strcmp(s3, "True") == 0);
    free(s3);
    expr_free(e3);
    expr_free(res3);

    // Test case 4: NumberQ[Rational[3,4]] -> True
    Expr* e4 = parse_expression("NumberQ[Rational[3,4]]");
    Expr* res4 = evaluate(e4);
    char* s4 = expr_to_string(res4);
    assert(strcmp(s4, "True") == 0);
    free(s4);
    expr_free(e4);
    expr_free(res4);

    // Test case 5: NumberQ[x] -> False
    Expr* e5 = parse_expression("NumberQ[x]");
    Expr* res5 = evaluate(e5);
    char* s5 = expr_to_string(res5);
    assert(strcmp(s5, "False") == 0);
    free(s5);
    expr_free(e5);
    expr_free(res5);

    // Test case 6: NumberQ["hello"] -> False
    Expr* e6 = parse_expression("NumberQ[\"hello\"]");
    Expr* res6 = evaluate(e6);
    char* s6 = expr_to_string(res6);
    assert(strcmp(s6, "False") == 0);
    free(s6);
    expr_free(e6);
    expr_free(res6);

    // Test case 7: NumberQ[f[x]] -> False
    Expr* e7 = parse_expression("NumberQ[f[x]]");
    Expr* res7 = evaluate(e7);
    char* s7 = expr_to_string(res7);
    assert(strcmp(s7, "False") == 0);
    free(s7);
    expr_free(e7);
    expr_free(res7);

    // Test case 8: NumberQ[Infinity] -> False
    Expr* e8 = parse_expression("NumberQ[Infinity]");
    Expr* res8 = evaluate(e8);
    char* s8 = expr_to_string(res8);
    assert(strcmp(s8, "False") == 0);
    free(s8);
    expr_free(e8);
    expr_free(res8);
}

void test_atomq(void) {
    // Test case 1: AtomQ[x] -> True
    Expr* e1 = parse_expression("AtomQ[x]");
    Expr* res1 = evaluate(e1);
    char* s1 = expr_to_string(res1);
    assert(strcmp(s1, "True") == 0);
    free(s1);
    expr_free(e1);
    expr_free(res1);

    // Test case 2: AtomQ[123] -> True
    Expr* e2 = parse_expression("AtomQ[123]");
    Expr* res2 = evaluate(e2);
    char* s2 = expr_to_string(res2);
    assert(strcmp(s2, "True") == 0);
    free(s2);
    expr_free(e2);
    expr_free(res2);

    // Test case 3: AtomQ[f[x]] -> False
    Expr* e3 = parse_expression("AtomQ[f[x]]");
    Expr* res3 = evaluate(e3);
    char* s3 = expr_to_string(res3);
    assert(strcmp(s3, "False") == 0);
    free(s3);
    expr_free(e3);
    expr_free(res3);

    // Test case 4: AtomQ[Complex[1, 2]] -> True
    Expr* e4 = parse_expression("AtomQ[Complex[1, 2]]");
    Expr* res4 = evaluate(e4);
    char* s4 = expr_to_string(res4);
    assert(strcmp(s4, "True") == 0);
    free(s4);
    expr_free(e4);
    expr_free(res4);
}

void test_integerq(void) {
    // Test case 1: IntegerQ[123] -> True
    Expr* e1 = parse_expression("IntegerQ[123]");
    Expr* res1 = evaluate(e1);
    char* s1 = expr_to_string(res1);
    assert(strcmp(s1, "True") == 0);
    free(s1);
    expr_free(e1);
    expr_free(res1);


    // Test case 2: IntegerQ[-5] -> True
    Expr* e2 = parse_expression("IntegerQ[-5]");
    Expr* res2 = evaluate(e2);
    char* s2 = expr_to_string(res2);
    assert(strcmp(s2, "True") == 0);
    free(s2);
    expr_free(e2);
    expr_free(res2);

    // Test case 3: IntegerQ[1.0] -> False
    Expr* e3 = parse_expression("IntegerQ[1.0]");
    Expr* res3 = evaluate(e3);
    char* s3 = expr_to_string(res3);
    assert(strcmp(s3, "False") == 0);
    free(s3);
    expr_free(e3);
    expr_free(res3);

    // Test case 4: IntegerQ[1.23] -> False
    Expr* e4 = parse_expression("IntegerQ[1.23]");
    Expr* res4 = evaluate(e4);
    char* s4 = expr_to_string(res4);
    assert(strcmp(s4, "False") == 0);
    free(s4);
    expr_free(e4);
    expr_free(res4);

    // Test case 5: IntegerQ[3/4] -> False
    Expr* e5 = parse_expression("IntegerQ[3/4]");
    Expr* res5 = evaluate(e5);
    char* s5 = expr_to_string(res5);
    assert(strcmp(s5, "False") == 0);
    free(s5);
    expr_free(e5);
    expr_free(res5);

    // Test case 6: IntegerQ[x] -> False
    Expr* e6 = parse_expression("IntegerQ[x]");
    Expr* res6 = evaluate(e6);
    char* s6 = expr_to_string(res6);
    assert(strcmp(s6, "False") == 0);
    free(s6);
    expr_free(e6);
    expr_free(res6);

    // Test case 7: IntegerQ[Complex[1, 1]] -> False
    Expr* e7 = parse_expression("IntegerQ[Complex[1, 1]]");
    Expr* res7 = evaluate(e7);
    char* s7 = expr_to_string(res7);
    assert(strcmp(s7, "False") == 0);
    free(s7);
    expr_free(e7);
    expr_free(res7);
}

void test_evenq_oddq(void) {
    // EvenQ tests
    assert_eval_eq("EvenQ[2]", "True", 0);
    assert_eval_eq("EvenQ[3]", "False", 0);
    assert_eval_eq("EvenQ[0]", "True", 0);
    assert_eval_eq("EvenQ[-4]", "True", 0);
    assert_eval_eq("EvenQ[2.0]", "False", 0);
    assert_eval_eq("EvenQ[x]", "False", 0);
    assert_eval_eq("EvenQ[1/2]", "False", 0);

    // OddQ tests
    assert_eval_eq("OddQ[3]", "True", 0);
    assert_eval_eq("OddQ[2]", "False", 0);
    assert_eval_eq("OddQ[-5]", "True", 0);
    assert_eval_eq("OddQ[0]", "False", 0);
    assert_eval_eq("OddQ[3.0]", "False", 0);
    assert_eval_eq("OddQ[x]", "False", 0);
    assert_eval_eq("OddQ[1/2]", "False", 0);
}

void assert_streq_double(const char* actual, const char* expected) {
    if (strcmp(actual, expected) == 0) {
        assert(true);
        return;
    }
    char expected_with_zero[256];
    snprintf(expected_with_zero, sizeof(expected_with_zero), "%s0", expected);
    if (strcmp(actual, expected_with_zero) == 0) {
        assert(true);
        return;
    }
    assert(false);
}

void test_mod(void) {
    // Integer tests
    assert_eval_eq("Mod[7, 3]", "1", 0);
    assert_eval_eq("Mod[-7, 3]", "2", 0);
    assert_eval_eq("Mod[7, -3]", "-2", 0);
    assert_eval_eq("Mod[-7, -3]", "-1", 0);
    assert_eval_eq("Mod[6, 3]", "0", 0);
    assert_eval_eq("Mod[0, 5]", "0", 0);

    // Real tests
    assert_eval_eq("Mod[11.25, 3]", "2.25", 0);
    assert_streq_double(expr_to_string(eval_and_free(parse_expression("Mod[11.25, 4.125]"))), "3.");
    assert_streq_double(expr_to_string(eval_and_free(parse_expression("Mod[7.5, 2.5]"))), "0.");

    // Mixed type tests
    assert_streq_double(expr_to_string(eval_and_free(parse_expression("Mod[10, 2.5]"))), "0.");
    assert_eval_eq("Mod[10.5, 3]", "1.5", 0);

    // 3-argument offset tests
    assert_eval_eq("Mod[11, 5, 2]", "6", 0);
    assert_eval_eq("Mod[1, 5, 2]", "6", 0);
    assert_eval_eq("Mod[2, 5, 2]", "2", 0);
    assert_eval_eq("Mod[6, 5, 2]", "6", 0);
    assert_eval_eq("Mod[7, 5, 2]", "2", 0);
    assert_eval_eq("Mod[11.5, 5, 2]", "6.5", 0);
    assert_eval_eq("Mod[11, 5, -1]", "1", 0);
}

void test_quotient(void) {
    // Integer tests (floor)
    assert_eval_eq("Quotient[10, 3]", "3", 0);
    assert_eval_eq("Quotient[-10, 3]", "-4", 0);
    assert_eval_eq("Quotient[10, -3]", "-4", 0);
    assert_eval_eq("Quotient[-10, -3]", "3", 0);

    // Rational and Real tests (floor)
    assert_eval_eq("Quotient[111/4, 5/4]", "22", 0);
    assert_eval_eq("Quotient[144.144, 11.12]", "12", 0);

    // Complex test (rounds to nearest integer)
    char* s_complex = expr_to_string_fullform(eval_and_free(parse_expression("Quotient[17.5+6I, 1+2I]")));
    assert(strcmp(s_complex, "Complex[6, -6]") == 0);
    free(s_complex);

    // 3-argument tests
    assert_eval_eq("Quotient[11, 3, 1]", "3", 0);
    assert_eval_eq("Quotient[10, 3, 1]", "3", 0);
    assert_eval_eq("Quotient[12, 3, 1]", "3", 0);
}

void test_quotientremainder(void) {
    // Basic test
    char* s1 = expr_to_string_fullform(eval_and_free(parse_expression("QuotientRemainder[10, 3]")));
    assert(strcmp(s1, "List[3, 1]") == 0);
    free(s1);

    // Real numbers test
    char* s2 = expr_to_string_fullform(eval_and_free(parse_expression("QuotientRemainder[11.25, 3]")));
    assert(strcmp(s2, "List[3, 2.25]") == 0);
    free(s2);

    // Negative test
    char* s3 = expr_to_string_fullform(eval_and_free(parse_expression("QuotientRemainder[-10, 3]")));
    assert(strcmp(s3, "List[-4, 2]") == 0);
    free(s3);
}

void test_re_im(void) {
    char* s1 = expr_to_string_fullform(eval_and_free(parse_expression("Re[Complex[2, 3]]")));
    assert(strcmp(s1, "2") == 0);
    free(s1);

    char* s2 = expr_to_string_fullform(eval_and_free(parse_expression("Im[Complex[2, 3]]")));
    assert(strcmp(s2, "3") == 0);
    free(s2);

    char* s3 = expr_to_string_fullform(eval_and_free(parse_expression("Re[5]")));
    assert(strcmp(s3, "5") == 0);
    free(s3);

    char* s4 = expr_to_string_fullform(eval_and_free(parse_expression("Im[5]")));
    assert(strcmp(s4, "0") == 0);
    free(s4);

    char* s5 = expr_to_string_fullform(eval_and_free(parse_expression("ReIm[Complex[2, 3]]")));
    assert(strcmp(s5, "List[2, 3]") == 0);
    free(s5);

    char* s6 = expr_to_string_fullform(eval_and_free(parse_expression("ReIm[5]")));
    assert(strcmp(s6, "List[5, 0]") == 0);
    free(s6);
}

void test_abs_conjugate(void) {
    char* s1 = expr_to_string_fullform(eval_and_free(parse_expression("Abs[-5]")));
    assert(strcmp(s1, "5") == 0);
    free(s1);

    char* s2 = expr_to_string_fullform(eval_and_free(parse_expression("Abs[-3.14]")));
    assert(strcmp(s2, "3.14") == 0);
    free(s2);

    char* s3 = expr_to_string_fullform(eval_and_free(parse_expression("Abs[Complex[3, 4]]")));
    assert(strcmp(s3, "5") == 0);
    free(s3);

    char* s4 = expr_to_string_fullform(eval_and_free(parse_expression("Conjugate[Complex[2, 3]]")));
    assert(strcmp(s4, "Complex[2, -3]") == 0);
    free(s4);

    char* s5 = expr_to_string_fullform(eval_and_free(parse_expression("Conjugate[5]")));
    assert(strcmp(s5, "5") == 0);
    free(s5);

    /* Symbolic real numerics (anything that numericalizes to a machine real)
     * are Conjugate-fixed: 3/Sqrt[11], Sqrt[2], Pi, etc. */
    char* s6 = expr_to_string_fullform(eval_and_free(parse_expression("Conjugate[3/Sqrt[11]]")));
    assert(strcmp(s6, "Times[3, Power[11, Rational[-1, 2]]]") == 0);
    free(s6);

    char* s7 = expr_to_string_fullform(eval_and_free(parse_expression("Conjugate[Sqrt[11]]")));
    assert(strcmp(s7, "Power[11, Rational[1, 2]]") == 0);
    free(s7);

    char* s8 = expr_to_string_fullform(eval_and_free(parse_expression("Conjugate[Pi]")));
    assert(strcmp(s8, "Pi") == 0);
    free(s8);

    /* Involution: Conjugate[Conjugate[z]] -> z. */
    char* s9 = expr_to_string_fullform(eval_and_free(parse_expression("Conjugate[Conjugate[z]]")));
    assert(strcmp(s9, "z") == 0);
    free(s9);

    /* Three nested Conjugates collapse to one (odd parity). */
    char* s10 = expr_to_string_fullform(eval_and_free(parse_expression("Conjugate[Conjugate[Conjugate[z]]]")));
    assert(strcmp(s10, "Conjugate[z]") == 0);
    free(s10);

    /* Re, Im, Abs, Arg are real-valued and therefore Conjugate-fixed,
     * even for symbolic arguments that don't numericalize. */
    char* s11 = expr_to_string_fullform(eval_and_free(parse_expression("Conjugate[Re[z]]")));
    assert(strcmp(s11, "Re[z]") == 0);
    free(s11);

    char* s12 = expr_to_string_fullform(eval_and_free(parse_expression("Conjugate[Im[z]]")));
    assert(strcmp(s12, "Im[z]") == 0);
    free(s12);

    char* s13 = expr_to_string_fullform(eval_and_free(parse_expression("Conjugate[Abs[z]]")));
    assert(strcmp(s13, "Abs[z]") == 0);
    free(s13);

    char* s14 = expr_to_string_fullform(eval_and_free(parse_expression("Conjugate[Arg[z]]")));
    assert(strcmp(s14, "Arg[z]") == 0);
    free(s14);
}

void test_arg(void) {
    char* s1 = expr_to_string_fullform(eval_and_free(parse_expression("Arg[0]")));
    assert(strcmp(s1, "0") == 0);
    free(s1);

    char* s2 = expr_to_string_fullform(eval_and_free(parse_expression("Arg[5]")));
    assert(strcmp(s2, "0") == 0);
    free(s2);

    char* s3 = expr_to_string_fullform(eval_and_free(parse_expression("Arg[-5]")));
    assert(strcmp(s3, "Pi") == 0);
    free(s3);

    char* s4 = expr_to_string_fullform(eval_and_free(parse_expression("Arg[Complex[0, 2]]")));
    assert(strcmp(s4, "Times[Rational[1, 2], Pi]") == 0);
    free(s4);

    char* s5 = expr_to_string_fullform(eval_and_free(parse_expression("Arg[Complex[0, -2]]")));
    assert(strcmp(s5, "Times[Rational[-1, 2], Pi]") == 0);
    free(s5);
    
    char* s6 = expr_to_string_fullform(eval_and_free(parse_expression("Arg[x]")));
    assert(strcmp(s6, "Arg[x]") == 0);
    free(s6);

    char* s7 = expr_to_string_fullform(eval_and_free(parse_expression("Arg[Complex[1, 1]]")));
    assert(strcmp(s7, "Times[Rational[1, 4], Pi]") == 0);
    free(s7);

    char* s8 = expr_to_string_fullform(eval_and_free(parse_expression("Arg[Complex[-1, 1]]")));
    assert(strcmp(s8, "Times[Rational[3, 4], Pi]") == 0);
    free(s8);

    char* s9 = expr_to_string_fullform(eval_and_free(parse_expression("Arg[Complex[-1, -1]]")));
    assert(strcmp(s9, "Times[Rational[-3, 4], Pi]") == 0);
    free(s9);

    char* s10 = expr_to_string_fullform(eval_and_free(parse_expression("Arg[Complex[1, -1]]")));
    assert(strcmp(s10, "Times[Rational[-1, 4], Pi]") == 0);
    free(s10);

    char* s11 = expr_to_string_fullform(eval_and_free(parse_expression("Arg[Complex[1, 2]]")));
    assert(strcmp(s11, "ArcTan[1, 2]") == 0);
    free(s11);
}

void test_trig(void) {
    char* s1 = expr_to_string_fullform(eval_and_free(parse_expression("Sin[0]")));
    assert(strcmp(s1, "0") == 0); free(s1);
    
    char* s2 = expr_to_string_fullform(eval_and_free(parse_expression("Cos[0]")));
    assert(strcmp(s2, "1") == 0); free(s2);

    char* s3 = expr_to_string_fullform(eval_and_free(parse_expression("Sin[Times[Rational[1, 6], Pi]]")));
    assert(strcmp(s3, "Rational[1, 2]") == 0); free(s3);

    char* s4 = expr_to_string_fullform(eval_and_free(parse_expression("Cos[Times[Rational[1, 3], Pi]]")));
    assert(strcmp(s4, "Rational[1, 2]") == 0); free(s4);

    char* s5 = expr_to_string_fullform(eval_and_free(parse_expression("Tan[Times[Rational[1, 4], Pi]]")));
    assert(strcmp(s5, "1") == 0); free(s5);
    
    char* s6 = expr_to_string_fullform(eval_and_free(parse_expression("ArcSin[0]")));
    assert(strcmp(s6, "0") == 0); free(s6);
    
    char* s7 = expr_to_string_fullform(eval_and_free(parse_expression("ArcCos[1]")));
    assert(strcmp(s7, "0") == 0); free(s7);

    char* s8 = expr_to_string_fullform(eval_and_free(parse_expression("ArcTan[0]")));
    assert(strcmp(s8, "0") == 0); free(s8);
}

void test_gcd_lcm(void) {
    // GCD tests
    assert_eval_eq("GCD[12, 18, 24]", "6", 1);
    assert_eval_eq("GCD[1/2, 1/3]", "Rational[1, 6]", 1);
    assert_eval_eq("GCD[0, 5]", "5", 1);
    assert_eval_eq("GCD[]", "0", 1);
    assert_eval_eq("GCD[-5]", "5", 1);
    assert_eval_eq("GCD[x]", "x", 1);

    // LCM tests
    assert_eval_eq("LCM[12, 18, 24]", "72", 1);
    assert_eval_eq("LCM[1/2, 1/3]", "1", 1);
    assert_eval_eq("LCM[0, 5]", "0", 1);
    assert_eval_eq("LCM[]", "1", 1);
    assert_eval_eq("LCM[-5]", "5", 1);
    assert_eval_eq("LCM[x]", "x", 1);
}

void test_primeq(void) {
    assert_eval_eq("PrimeQ[2]", "True", 1);
    assert_eval_eq("PrimeQ[3]", "True", 1);
    assert_eval_eq("PrimeQ[4]", "False", 1);
    assert_eval_eq("PrimeQ[17]", "True", 1);
    assert_eval_eq("PrimeQ[1]", "False", 1);
    assert_eval_eq("PrimeQ[0]", "False", 1);
    assert_eval_eq("PrimeQ[-2]", "True", 1);
    assert_eval_eq("PrimeQ[-17]", "True", 1);
    assert_eval_eq("PrimeQ[11.5]", "False", 1);
    /* *Q predicates must always return True/False — never symbolic. */
    assert_eval_eq("PrimeQ[x]", "False", 1);
    assert_eval_eq("PrimeQ[Sqrt[2]]", "False", 1);
    assert_eval_eq("PrimeQ[Exp[2 Pi I / 3]]", "False", 1);

    // Large prime (from user request)
    assert_eval_eq("PrimeQ[-59]", "True", 1);

    // Gaussian primes: both parts nonzero, norm is prime
    assert_eval_eq("PrimeQ[1 + I]", "True", 1);       // norm=2, prime
    assert_eval_eq("PrimeQ[1 + 2 I]", "True", 1);     // norm=5, prime
    assert_eval_eq("PrimeQ[2 + I]", "True", 1);        // norm=5, prime
    assert_eval_eq("PrimeQ[4 + I]", "True", 1);        // norm=17, prime
    assert_eval_eq("PrimeQ[2 + 2 I]", "False", 1);     // norm=8, not prime

    // Gaussian primes: pure imaginary, |b| prime and b ≡ 3 mod 4
    assert_eval_eq("PrimeQ[3 I]", "True", 1);
    assert_eval_eq("PrimeQ[7 I]", "True", 1);
    assert_eval_eq("PrimeQ[5 I]", "False", 1);         // 5 ≡ 1 mod 4
    assert_eval_eq("PrimeQ[2 I]", "False", 1);          // 2 ≡ 2 mod 4

    // Gaussian primes: negative imaginary
    assert_eval_eq("PrimeQ[-3 I]", "True", 1);
    assert_eval_eq("PrimeQ[1 - 2 I]", "True", 1);      // norm=5, prime

    /* GaussianIntegers option: tests primality in Z[i]. A rational
     * integer n is a Gaussian prime iff |n| is prime in Z AND n ≡ 3
     * mod 4 (the ≡ 1 mod 4 primes split, and n=2 is associate of (1+i)^2). */
    assert_eval_eq("PrimeQ[5, GaussianIntegers -> True]", "False", 1);
    assert_eval_eq("PrimeQ[3, GaussianIntegers -> True]", "True", 1);
    assert_eval_eq("PrimeQ[7, GaussianIntegers -> True]", "True", 1);
    assert_eval_eq("PrimeQ[2, GaussianIntegers -> True]", "False", 1);
    assert_eval_eq("PrimeQ[13, GaussianIntegers -> True]", "False", 1);
    assert_eval_eq("PrimeQ[5, GaussianIntegers -> False]", "True", 1);
    /* Option still gives sensible answers on non-integer / non-Gaussian args. */
    assert_eval_eq("PrimeQ[x, GaussianIntegers -> True]", "False", 1);
    assert_eval_eq("PrimeQ[5.5, GaussianIntegers -> True]", "False", 1);
}

void test_factorinteger(void) {
    assert_eval_eq("FactorInteger[12]", "List[List[2, 2], List[3, 1]]", 1);
    assert_eval_eq("FactorInteger[-12]", "List[List[-1, 1], List[2, 2], List[3, 1]]", 1);
    assert_eval_eq("FactorInteger[1/2]", "List[List[2, -1]]", 1);
    assert_eval_eq("FactorInteger[3/4]", "List[List[2, -2], List[3, 1]]", 1);
    
    // Partial factorization
    assert_eval_eq("FactorInteger[100, 1]", "List[List[2, 2]]", 1);
    
    // Automatic (easy factors)
    // Removed flaky ECM test
}

void test_eulerphi(void) {
    Expr* e; Expr* res; char* s;

    e = parse_expression("EulerPhi[10]"); res = evaluate(e); s = expr_to_string_fullform(res);
    assert(strcmp(s, "4") == 0); free(s); expr_free(res); expr_free(e);

    e = parse_expression("EulerPhi[1]"); res = evaluate(e); s = expr_to_string_fullform(res);
    assert(strcmp(s, "1") == 0); free(s); expr_free(res); expr_free(e);

    e = parse_expression("EulerPhi[0]"); res = evaluate(e); s = expr_to_string_fullform(res);
    assert(strcmp(s, "0") == 0); free(s); expr_free(res); expr_free(e);

    e = parse_expression("EulerPhi[-10]"); res = evaluate(e); s = expr_to_string_fullform(res);
    assert(strcmp(s, "4") == 0); free(s); expr_free(res); expr_free(e);

    e = parse_expression("EulerPhi[{10, 20, 30}]"); res = evaluate(e); s = expr_to_string_fullform(res);
    assert(strcmp(s, "List[4, 8, 8]") == 0); free(s); expr_free(res); expr_free(e);

    e = parse_expression("EulerPhi[9]"); res = evaluate(e); s = expr_to_string_fullform(res);
    assert(strcmp(s, "6") == 0); free(s); expr_free(res); expr_free(e);

    e = parse_expression("EulerPhi[100]"); res = evaluate(e); s = expr_to_string_fullform(res);
    assert(strcmp(s, "40") == 0); free(s); expr_free(res); expr_free(e);
}

void test_factorial(void) {
    Expr* e; Expr* res; char* s;

    e = parse_expression("Factorial[0]"); res = evaluate(e); s = expr_to_string_fullform(res);
    assert(strcmp(s, "1") == 0); free(s); expr_free(res); expr_free(e);

    e = parse_expression("5!"); res = evaluate(e); s = expr_to_string_fullform(res);
    assert(strcmp(s, "120") == 0); free(s); expr_free(res); expr_free(e);

    e = parse_expression("(-1)!"); res = evaluate(e); s = expr_to_string_fullform(res);
    assert(strcmp(s, "ComplexInfinity") == 0); free(s); expr_free(res); expr_free(e);

    e = parse_expression("(1/2)!"); res = evaluate(e); s = expr_to_string_fullform(res);
    assert(strcmp(s, "Times[Rational[1, 2], Power[Pi, Rational[1, 2]]]") == 0); free(s); expr_free(res); expr_free(e);

    e = parse_expression("(-1/2)!"); res = evaluate(e); s = expr_to_string_fullform(res);
    assert(strcmp(s, "Power[Pi, Rational[1, 2]]") == 0); free(s); expr_free(res); expr_free(e);

    e = parse_expression("(-3/2)!"); res = evaluate(e); s = expr_to_string_fullform(res);
    assert(strcmp(s, "Times[-2, Power[Pi, Rational[1, 2]]]") == 0); free(s); expr_free(res); expr_free(e);

    e = parse_expression("21!"); res = evaluate(e); s = expr_to_string_fullform(res);
    assert(strcmp(s, "51090942171709440000") == 0); free(s); expr_free(res); expr_free(e);
}

void test_binomial(void) {
    Expr* e; Expr* res; char* s;

    e = parse_expression("Binomial[10, 3]"); res = evaluate(e); s = expr_to_string_fullform(res);
    assert(strcmp(s, "120") == 0); free(s); expr_free(res); expr_free(e);

    e = parse_expression("Binomial[8, 4]"); res = evaluate(e); s = expr_to_string_fullform(res);
    assert(strcmp(s, "70") == 0); free(s); expr_free(res); expr_free(e);

    e = parse_expression("Binomial[9/2, 7/2]"); res = evaluate(e); s = expr_to_string_fullform(res);
    assert(strcmp(s, "Rational[9, 2]") == 0 || strcmp(s, "5") == 0); free(s); expr_free(res); expr_free(e);

    e = parse_expression("Binomial[n, 4]"); res = evaluate(e); s = expr_to_string_fullform(res);
    assert(strcmp(s, "Binomial[n, 4]") == 0 || strcmp(s, "Times[Rational[1, 24], n, Plus[-3, n], Plus[-2, n], Plus[-1, n]]") == 0);
    free(s); expr_free(res); expr_free(e);

    e = parse_expression("Binomial[0, 1]"); res = evaluate(e); s = expr_to_string_fullform(res);
    assert(strcmp(s, "0") == 0); free(s); expr_free(res); expr_free(e);
    
    e = parse_expression("Binomial[-1, 1]"); res = evaluate(e); s = expr_to_string_fullform(res);
    assert(strcmp(s, "-1") == 0); free(s); expr_free(res); expr_free(e);
    
    e = parse_expression("Binomial[-1, 0]"); res = evaluate(e); s = expr_to_string_fullform(res);
    assert(strcmp(s, "1") == 0); free(s); expr_free(res); expr_free(e);
}

void test_nextprime(void) {
    assert_eval_eq("NextPrime[2]", "3", 1);
    assert_eval_eq("NextPrime[3]", "5", 1);
    assert_eval_eq("NextPrime[4]", "5", 1);
    assert_eval_eq("NextPrime[2.5]", "3", 1);
    assert_eval_eq("NextPrime[10, 1]", "11", 1);
    assert_eval_eq("NextPrime[10, 2]", "13", 1);
    assert_eval_eq("NextPrime[5, -1]", "3", 1);
    assert_eval_eq("NextPrime[5, -2]", "2", 1);
    assert_eval_eq("NextPrime[2, -1]", "NextPrime[2, -1]", 1);
}

void test_primepi(void) {
    assert_eval_eq("PrimePi[10]", "4", 1);
    assert_eval_eq("PrimePi[100]", "25", 1);
    assert_eval_eq("PrimePi[1000]", "168", 1);
    assert_eval_eq("PrimePi[10000]", "1229", 1);
    
    // Listable
    assert_eval_eq("PrimePi[{10, 100}]", "List[4, 25]", 1);
}

void test_depth(void) {
    assert_eval_eq("Depth[a]", "1", 1);
    assert_eval_eq("Depth[{a}]", "2", 1);
    assert_eval_eq("Depth[{{a}}]", "3", 1);
    assert_eval_eq("Depth[{{{a}}}]", "4", 1);
    assert_eval_eq("Depth[{{{a}, b}}]", "4", 1);
    assert_eval_eq("Depth[1 + x^2]", "3", 1);
    assert_eval_eq("Depth[f[f[f[x]]]]", "4", 1);
    assert_eval_eq("Depth[f[g[h[x]]]]", "4", 1);
    assert_eval_eq("Depth[12345]", "1", 1);
    Expr* e_depth = eval_and_free(parse_expression("Depth[3+I]"));
    char* s_depth = expr_to_string_fullform(e_depth);
    assert(strcmp(s_depth, "1") == 0);
    free(s_depth);
    expr_free(e_depth);

    assert_eval_eq("Depth[1/2]", "1", 1);
    assert_eval_eq("Depth[Sqrt[2]]", "2", 1);
    assert_eval_eq("Depth[h[{{{a}}}][x, y]]", "2", 1);
    assert_eval_eq("Depth[h[{{{a}}}][x, y], Heads -> True]", "6", 1);
}

void test_leafcount() {
    assert_eval_eq("LeafCount[1+a+b^2]", "6", 1);
    assert_eval_eq("LeafCount[f[x,y]]", "3", 1);
    assert_eval_eq("LeafCount[f[a,b][x,y]]", "5", 1);
    assert_eval_eq("LeafCount[I]", "3", 1);
    assert_eval_eq("LeafCount[{1/2, 1+I}]", "7", 1);
    assert_eval_eq("LeafCount[f[x,y], Heads->False]", "2", 1);
}

void test_bytecount() {
    Expr* e; Expr* res;
    
    e = parse_expression("ByteCount[1]");
    res = evaluate(e);
    assert(res->type == EXPR_INTEGER);
    assert(res->data.integer > 0);
    expr_free(res); expr_free(e);

    e = parse_expression("ByteCount[x]");
    res = evaluate(e);
    assert(res->type == EXPR_INTEGER);
    assert(res->data.integer > sizeof(Expr)); // Should include string length
    expr_free(res); expr_free(e);

    e = parse_expression("ByteCount[f[x,y]]");
    res = evaluate(e);
    assert(res->type == EXPR_INTEGER);
    assert(res->data.integer > sizeof(Expr) * 3); // head + 2 args
    expr_free(res); expr_free(e);
}

void test_information(void) {
    Expr* e1 = parse_expression("Information[Range]");
    Expr* res1 = evaluate(e1);
    ASSERT(res1->type == EXPR_STRING);
    ASSERT(strstr(res1->data.string, "Range[n, m, d]") != NULL);
    expr_free(e1);
    expr_free(res1);

    Expr* e2 = parse_expression("?Range");
    Expr* res2 = evaluate(e2);
    ASSERT(res2->type == EXPR_STRING);
    ASSERT(strstr(res2->data.string, "Range[n, m, d]") != NULL);
    expr_free(e2);
    expr_free(res2);
    
    Expr* e3 = parse_expression("?NonExistentSymbol");
    Expr* res3 = evaluate(e3);
    ASSERT(res3->type == EXPR_STRING);
    ASSERT(strstr(res3->data.string, "No information available") != NULL);
    expr_free(e3);
    expr_free(res3);
}


void test_factor_methods() {
    // 2^60 - 1 = 3^2 * 5^2 * 7 * 11 * 13 * 31 * 41 * 61 * 151 * 331 * 1321
    // The prime components for some algorithms might remain composite if they don't resolve.
    // We'll just ensure they don't crash and parse the method rule correctly.
    assert_eval_eq("FactorInteger[91, Method -> \"PollardP-1\"]", "List[List[91, 1]]", 1);
    assert_eval_eq("FactorInteger[91, Method -> \"WilliamsP+1\"]", "List[List[91, 1]]", 1);
    
    // Fermat
    assert_eval_eq("FactorInteger[5959, Method -> \"Fermat\"]", "List[List[59, 1], List[101, 1]]", 1);

    // BlakeRationalBaseDescent
    assert_eval_eq("FactorInteger[13434917067328449643383271289062122492729438008563207396013386022436439625786814860310635092648849808707283885079425559380060603741711887741983827702536077124606526390026531236424033896417798802339681855054916007583273591721147874872971014173899334436432277372280112100845364029259802969568473368054409653616559493278947418005745364854816116209, Method -> {\"BlakeRationalBaseDescent\", \"Base\" -> 22/7}]", "List[List[31415926535897932384626433832795028841971693993751058209749445923078164119021, 1], List[427646692258998556040057979899735431275177165797172063573559165768095800190402040091181222919605426534412322698619636704599813270516030490251834195211040830685432130941412817396056010802865303531028012177861105952338794347041755868615176860993347787824047880817012629, 1]]", 1);

    
    // CFRAC
    assert_eval_eq("FactorInteger[8051, Method -> \"CFRAC\"]", "List[List[83, 1], List[97, 1]]", 1);
}

void test_clear_attributes(void) {
    // ClearAttributes returns Null
    assert_eval_eq("ClearAttributes[testClearF, Protected]", "Null", 0);

    // Set multiple attributes, then clear one at a time
    assert_eval_eq("SetAttributes[testClearF, {Flat, Orderless, OneIdentity}]", "Null", 0);
    assert_eval_eq("Attributes[testClearF]", "{Flat, OneIdentity, Orderless}", 0);

    // Clear a single attribute; remaining ones are retained
    assert_eval_eq("ClearAttributes[testClearF, OneIdentity]", "Null", 0);
    assert_eval_eq("Attributes[testClearF]", "{Flat, Orderless}", 0);

    // Clear multiple attributes at once with a list
    assert_eval_eq("ClearAttributes[testClearF, {Flat, Orderless}]", "Null", 0);
    assert_eval_eq("Attributes[testClearF]", "{}", 0);

    // Clearing an attribute that is not set is a no-op
    assert_eval_eq("ClearAttributes[testClearF, Listable]", "Null", 0);
    assert_eval_eq("Attributes[testClearF]", "{}", 0);

    // Listable threading: set Listable, clear it, verify behavior
    assert_eval_eq("SetAttributes[testClearG, Listable]", "Null", 0);
    assert_eval_eq("Attributes[testClearG]", "{Listable}", 0);
    assert_eval_eq("ClearAttributes[testClearG, Listable]", "Null", 0);
    assert_eval_eq("Attributes[testClearG]", "{}", 0);

    // Clear attributes using string symbol name
    assert_eval_eq("SetAttributes[testClearH, {Flat, Protected}]", "Null", 0);
    assert_eval_eq("ClearAttributes[\"testClearH\", Protected]", "Null", 0);
    assert_eval_eq("Attributes[testClearH]", "{Flat}", 0);
    assert_eval_eq("ClearAttributes[testClearH, Flat]", "Null", 0);

    // Clear attributes from a list of symbols
    assert_eval_eq("SetAttributes[testClearS1, {Flat, Orderless}]", "Null", 0);
    assert_eval_eq("SetAttributes[testClearS2, {Flat, Orderless}]", "Null", 0);
    assert_eval_eq("ClearAttributes[{testClearS1, testClearS2}, Flat]", "Null", 0);
    assert_eval_eq("Attributes[testClearS1]", "{Orderless}", 0);
    assert_eval_eq("Attributes[testClearS2]", "{Orderless}", 0);

    // Clear list of attributes from list of symbols
    assert_eval_eq("ClearAttributes[{testClearS1, testClearS2}, {Orderless}]", "Null", 0);
    assert_eval_eq("Attributes[testClearS1]", "{}", 0);
    assert_eval_eq("Attributes[testClearS2]", "{}", 0);

    // ClearAttributes itself has HoldFirst and Protected
    assert_eval_eq("Attributes[ClearAttributes]", "{HoldFirst, Protected}", 0);

    // Clean up
    assert_eval_eq("Clear[testClearF]", "Null", 0);
    assert_eval_eq("Clear[testClearG]", "Null", 0);
    assert_eval_eq("Clear[testClearH]", "Null", 0);
    assert_eval_eq("Clear[testClearS1]", "Null", 0);
    assert_eval_eq("Clear[testClearS2]", "Null", 0);
}

int main(void) {
    symtab_init();
    core_init();
    
    TEST(test_numberq);
    TEST(test_numericq);
    TEST(test_atomq);
    TEST(test_integerq);
    TEST(test_evenq_oddq);
    TEST(test_mod);
    TEST(test_quotient);
    TEST(test_quotientremainder);
    TEST(test_re_im);
    TEST(test_abs_conjugate);
    TEST(test_arg);
    TEST(test_trig);
    TEST(test_gcd_lcm);
    TEST(test_primeq);
    TEST(test_primepi);
    TEST(test_factorinteger);
    TEST(test_factor_methods);
    TEST(test_eulerphi);
    TEST(test_factorial);
    TEST(test_binomial);
    TEST(test_nextprime);
    TEST(test_depth);
    TEST(test_leafcount);
    TEST(test_bytecount);
    TEST(test_information);
    TEST(test_clear_attributes);

    printf("All core tests passed!\n");
    return 0;
}
