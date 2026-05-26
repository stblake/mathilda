#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include <string.h>
#include <stdlib.h>
#include "print.h"

void run_test(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* res_str = expr_to_string_fullform(res);
    if (strcmp(res_str, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, res_str);
        ASSERT_STR_EQ(res_str, expected);
    }
    free(res_str);
    expr_free(e);
    expr_free(res);
}

void test_block() {
    Expr* te = parse_expression("x = 1");
    Expr* tr = evaluate(te);
    expr_free(te); expr_free(tr);
    
    Expr* e1 = parse_expression("Block[{x = 2}, x + 1]");
    Expr* res1 = evaluate(e1);
    ASSERT(res1->type == EXPR_INTEGER);
    ASSERT(res1->data.integer == 3);
    
    Expr* e2 = parse_expression("x");
    Expr* res2 = evaluate(e2);
    ASSERT(res2->data.integer == 1);
    
    expr_free(e1); expr_free(res1); 
    expr_free(e2); expr_free(res2);
}

void test_module() {
    Expr* te = parse_expression("x = 1");
    Expr* tr = evaluate(te);
    expr_free(te); expr_free(tr);

    // Module should rename local x and not affect global x
    Expr* e1 = parse_expression("Module[{x = 10}, x + 5]");
    Expr* res1 = evaluate(e1);
    ASSERT(res1->data.integer == 15);
    
    Expr* e2 = parse_expression("x");
    Expr* res2 = evaluate(e2);
    ASSERT(res2->data.integer == 1);
    
    expr_free(e1); expr_free(res1);
    expr_free(e2); expr_free(res2);
}

void test_with() {
    Expr* te = parse_expression("x = 1");
    Expr* tr = evaluate(te);
    expr_free(te); expr_free(tr);

    Expr* e1 = parse_expression("With[{x = 100}, x / 2]");
    Expr* res1 = evaluate(e1);
    ASSERT(res1->data.integer == 50);
    
    Expr* e2 = parse_expression("x");
    Expr* res2 = evaluate(e2);
    ASSERT(res2->data.integer == 1);
    
    expr_free(e1); expr_free(res1);
    expr_free(e2); expr_free(res2);
}

void test_powermod() {
    run_test("PowerMod[2, 10, 3]", "1");
    run_test("PowerMod[3, 2, 7]", "2");
    run_test("PowerMod[3, -2, 7]", "4");
    run_test("PowerMod[2, {10, 11, 12, 13, 14}, 5]", "List[4, 3, 1, 2, 4]");
    run_test("PowerMod[3, -1, 7]", "5");
    run_test("PowerMod[3, -1, 6]", "PowerMod[3, -1, 6]"); // no inverse
}

/* Integer-exponent identities and edge cases that should hold for any modulus
 * and don't depend on a hard-coded reference value. */
void test_powermod_identities() {
    run_test("PowerMod[5, 0, 7]", "1");                 // a^0 = 1 (mod m)
    run_test("PowerMod[5, 0, 1]", "0");                 // mod 1 everything is 0
    run_test("PowerMod[7, 1, 11]", "7");                // a^1 = a (mod m)
    run_test("PowerMod[0, 5, 7]", "0");                 // 0^positive = 0
    run_test("PowerMod[123, 1, 100]", "23");            // reduction without power
    run_test("PowerMod[-5, 2, 7]", "4");                // negative base reduces correctly
    run_test("PowerMod[5, 2, -7]", "4");                // negative modulus treated as |m|
    /* The modulus 0 case has no sensible value; PowerMod stays unevaluated. */
    run_test("PowerMod[3, 2, 0]", "PowerMod[3, 2, 0]");
    /* Modular inverse path. */
    run_test("PowerMod[2, -1, 7]", "4");                // 2*4 = 8 ≡ 1 (mod 7)
    run_test("PowerMod[1, -1, 10]", "1");
    run_test("PowerMod[3, -1, 6]", "PowerMod[3, -1, 6]"); // gcd(3,6)=3, no inverse
}

/* Bignum integer-exponent inputs: a, b, m may each be far beyond int64.
 * Reference values were precomputed from the same algorithm in Mathematica /
 * Python's pow().  These exercise GMP's mpz_powm fast path. */
void test_powermod_bignum() {
    /* Bignum modulus, bignum base, bignum exponent — Fermat: a^(p-1) ≡ 1 */
    run_test("PowerMod[2, 10^18 + 8, 10^18 + 9]", "1");
    run_test("PowerMod[3, 10^18 + 8, 10^18 + 9]", "1");
    /* Bignum base mod bignum prime. */
    run_test("PowerMod[10^100, 2, 10^30 + 7]",
             "11764900000000000000000000");
    /* Very large modulus, moderate exponent. */
    run_test("PowerMod[2, 1000, 10^50 + 1]",
             "20435445070522076575670377779817606475351265541240");
    /* Bignum exponent with int64-fitting modulus — bignum exponent path. */
    run_test("PowerMod[2, 10^9, 10^18 + 9]", "904186674203936980");
    /* Mersenne-prime exponent of 2 — classical result mod the prime. */
    run_test("PowerMod[2, 2^31 - 1, 2^31 - 1]", "2");
    /* Bignum modular inverse: 12345678901234567890 mod (10^30 + 7). */
    run_test("PowerMod[12345678901234567890, -1, 10^30 + 7]",
             "973526894364094738699762122176");
    /* Negative bignum exponent (computes inverse, then powers it). */
    run_test("PowerMod[2, -10^9, 10^18 + 9]", "728548850360544027");
    /* Round-trip identity: PowerMod[PowerMod[a,b,m], 1, m] == a^b mod m. */
    run_test("PowerMod[PowerMod[7, 1234, 10^20 + 39], 1, 10^20 + 39]",
             "10656971757803551315");
}

/* Modular-root case (rational exponent p/q): PowerMod[a, p/q, m] returns an
 * x with x^q ≡ a^p (mod m) when one exists; otherwise stays unevaluated. */
void test_powermod_modular_root() {
    /* The bug from the original report. */
    run_test("PowerMod[3, 1/2, 2]", "1");

    /* Square roots mod small primes. */
    run_test("PowerMod[4, 1/2, 7]", "2");                   // 2^2 = 4
    run_test("PowerMod[2, 1/2, 7]", "3");                   // 3^2 ≡ 2 (mod 7)
    run_test("PowerMod[1, 1/2, 7]", "1");
    run_test("PowerMod[0, 1/2, 7]", "0");
    /* No square root: 3 is a QNR mod 7. */
    run_test("PowerMod[3, 1/2, 7]", "PowerMod[3, Rational[1, 2], 7]");

    /* Cube roots — coprime-exponent closed form (gcd(3, p-1) = 1). */
    run_test("PowerMod[2, 1/3, 5]", "3");                   // 3^3 = 27 ≡ 2 (mod 5)
    run_test("PowerMod[2, 1/3, 31]", "4");                  // 4^3 = 64 ≡ 2 (mod 31)
    /* Higher roots. */
    run_test("PowerMod[16, 1/4, 17]", "2");                 // 2^4 = 16 (mod 17)

    /* Composite modulus, square root via CRT. */
    run_test("PowerMod[100, 1/2, 17 * 19 * 23]", "10");     // 10^2 = 100 (mod 7429)
    run_test("PowerMod[4, 1/2, 6]", "2");                   // 2^2 = 4 (mod 6)
    /* 2 has no square root mod 19, so no root mod 17*19 either. */
    run_test("PowerMod[2, 1/2, 17 * 19]", "PowerMod[2, Rational[1, 2], 323]");

    /* Cube root mod a composite. */
    run_test("PowerMod[7, 1/3, 2 * 3 * 5 * 7 * 11]", "1183");

    /* Prime-power modulus (Hensel lift). */
    run_test("PowerMod[16, 1/4, 81]", "2");                 // 2^4 = 16 (mod 3^4)

    /* Tonelli-Shanks branch with p ≡ 1 (mod 4): 1000000007 ≡ 3 (mod 4) so
     * the p ≡ 3 (mod 4) shortcut is used; 4 has root 2 (and p-2). */
    run_test("PowerMod[4, 1/2, 1000000007]", "2");
    run_test("PowerMod[9, 1/2, 1000000007]", "3");

    /* Bignum modulus: square root mod a large prime via Tonelli-Shanks. */
    run_test("PowerMod[123, 1/2, 2^31 - 1]", "681925776");
    /* sqrt(2) mod (10^18 + 9): result was verified by squaring back. */
    run_test("PowerMod[2, 1/2, 10^18 + 9]", "742174169206529574");

    /* Rational exponent p/q with p != 1. */
    run_test("PowerMod[2, 3/2, 7]", "1");                   // sqrt(2^3=8 ≡ 1) = 1
    run_test("PowerMod[3, 2/3, 11]", "4");                  // cbrt(9) ≡ 4 (mod 11)

    /* Negative numerator: -1/2 means modular inverse, then square root. */
    run_test("PowerMod[4, -1/2, 11]", "5");                 // sqrt(4^-1 = 3) = 5

    /* Symbolic base or modulus: unevaluated. */
    run_test("PowerMod[a, 1/2, 5]", "PowerMod[a, Rational[1, 2], 5]");
    run_test("PowerMod[2, 1/2, m]", "PowerMod[2, Rational[1, 2], m]");

    /* m = 0 with rational exponent: unevaluated. */
    run_test("PowerMod[2, 1/2, 0]", "PowerMod[2, Rational[1, 2], 0]");

    /* Verification: PowerMod[PowerMod[a, 1/2, p], 2, p] == a (mod p). */
    run_test("PowerMod[PowerMod[2, 1/2, 7] ^ 2, 1, 7]", "2");
}

/* Listable threading: PowerMod has the Listable attribute, so list args
 * thread element-wise. */
void test_powermod_listable() {
    run_test("PowerMod[2, {1, 2, 3, 4}, 7]",
             "List[2, 4, 1, 2]");
    run_test("PowerMod[{2, 3, 5}, 2, 7]",
             "List[4, 2, 4]");
    run_test("PowerMod[2, 1, {3, 5, 7}]",
             "List[2, 2, 2]");
    /* Mixed integer / rational exponent in a list. */
    run_test("PowerMod[{4, 9, 16}, 1/2, 17]",
             "List[2, 3, 4]");
}

int main() {
    symtab_init();
    core_init();
    
    TEST(test_block);
    TEST(test_module);
    TEST(test_with);
    TEST(test_powermod);
    TEST(test_powermod_identities);
    TEST(test_powermod_bignum);
    TEST(test_powermod_modular_root);
    TEST(test_powermod_listable);

    printf("All modular tests passed!\n");
    return 0;
}
