
#include "parse.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include "symtab.h"
#include "core.h"
#include "eval.h"


void test_parse_atomics() {
    struct {
        const char* input;
        ExprType expected_type;
    } tests[] = {
        {"123", EXPR_INTEGER},
        {"-45", EXPR_INTEGER},
        {"3.14", EXPR_REAL},
	{"3.2E4", EXPR_REAL},
        {"x", EXPR_SYMBOL},
        {"\"text\"", EXPR_STRING}
    };

    for (size_t i = 0; i < sizeof(tests)/sizeof(tests[0]); i++) {
        Expr* e = parse_expression(tests[i].input);
        ASSERT(e != NULL);
        ASSERT(e->type == tests[i].expected_type);
        printf("PASS: %-10s → ", tests[i].input);
        expr_print(e);
        ASSERT_MSG(e->type == tests[i].expected_type, 
                 "For input '%s': expected type %d, got %d", 
                 tests[i].input, tests[i].expected_type, e->type);
        printf("\n");
        expr_free(e);
    }
}

void test_parse_functions() {
    const char* tests[] = {
        "f[]",
        "f[x]",
        "f[x,y]",
        "f[g[x]]",
        "f[1,\"two\"]"
    };

    for (size_t i = 0; i < sizeof(tests)/sizeof(tests[0]); i++) {
        Expr* e = parse_expression(tests[i]);
        ASSERT(e != NULL);
        ASSERT(e->type == EXPR_FUNCTION);
        printf("PASS: %-15s → ", tests[i]);
        expr_print(e);
        printf("\n");
        expr_free(e);
    }
}

void test_parse_lists() {
    const char* tests[] = {
        "{}",
        "{1}",
        "{1,2,3}",
        "{x,f[y],{}}"
    };

    for (size_t i = 0; i < sizeof(tests)/sizeof(tests[0]); i++) {
        Expr* e = parse_expression(tests[i]);
        ASSERT(e != NULL);
        ASSERT(e->type == EXPR_FUNCTION);
        ASSERT(strcmp(e->data.function.head->data.symbol.name, "List") == 0);
        printf("PASS: %-15s → ", tests[i]);
        expr_print(e);
        printf("\n");
        expr_free(e);
    }
}

void test_infix_operators() {
    const char* tests[] = {
        "a + b",
        "a * b + c",
        "a + b * c",
        "a ^ b * c",
        "a + b + c"
    };
    
    for (size_t i = 0; i < sizeof(tests)/sizeof(tests[0]); i++) {
        Expr* e = parse_expression(tests[i]);
        ASSERT(e != NULL);
        printf("Parsed: %s", tests[i]);
        expr_print(e);
        expr_free(e);
    }
}

void test_postfix_operators() {
    const char* tests[] = {
        "x // f",
        "4.13 // AtomQ",
        "a + b // g",
        "a // b // c"
    };

    for (size_t i = 0; i < sizeof(tests)/sizeof(tests[0]); i++) {
        Expr* e = parse_expression(tests[i]);
        ASSERT(e != NULL);
        printf("Parsed: %s -> ", tests[i]);
        expr_print(e);
        printf("\n");
        expr_free(e);
    }
}

void test_prefix_operators() {
    const char* tests[] = {
        "f @ x",
        "AtomQ @ 4.13",
        "g @ a + b",
        "c @ b @ a",
        "f @ x // g"
    };

    for (size_t i = 0; i < sizeof(tests)/sizeof(tests[0]); i++) {
        Expr* e = parse_expression(tests[i]);
        ASSERT(e != NULL);
        printf("Parsed: %s -> ", tests[i]);
        expr_print(e);
        printf("\n");
        expr_free(e);
    }
}

void test_fullform_parsing() {
    struct {
        const char* input;
        const char* expected;
    } tests[] = {
        {"f[]",          "f[]"},
        {"f[x]",         "f[x]"},
        {"f[x,y]",       "f[x, y]"},
        {"f[g[x], f[x]]",      "f[g[x],f[x]]"},
        {"{1,2,3}",      "List[1, 2, 3]"}
    };

    for (size_t i = 0; i < sizeof(tests)/sizeof(tests[0]); i++) {
        Expr* e = parse_expression(tests[i].input);
	printf("Parsed: %s", tests[i].input);
	expr_print(e);
        expr_free(e);
    }
}

void test_negative_numbers() {
    struct {
        const char* input;
        const char* expected;
    } tests[] = {
        {"-123", "-123"},
        {"-3.14", "-3.14"}
    };
    
    for (size_t i = 0; i < sizeof(tests)/sizeof(tests[0]); i++) {
        Expr* e = parse_expression(tests[i].input);
        ASSERT(e != NULL);
        printf("PASS: %-10s → ", tests[i].input);
        expr_print(e);
        printf("\n");
        expr_free(e);
    }
}

void test_parentheses() {
    struct {
        const char* input;
        const char* expected;
    } tests[] = {
        {"(1)", "1"},
        {"(x + y)", "Plus[x, y]"},
        {"(a * (b + c))", "Times[a, Plus[b, c]]"},
        {"- (x + y)", "Times[-1, Plus[x, y]]"},
        {"-(x + y)", "Times[-1, Plus[x, y]]"}	
    };
    
    for (size_t i = 0; i < sizeof(tests)/sizeof(tests[0]); i++) {
        Expr* e = parse_expression(tests[i].input);
        ASSERT(e != NULL);
        
        printf("PASS: %-15s → \n", tests[i].input);
	expr_print(e);
        expr_free(e);
    }
}

void test_parse_patterns() {
    const char* tests[] = {"_", "__", "___", "x_", "x__", "x_Integer", "_Integer"};
    for (size_t i = 0; i < sizeof(tests)/sizeof(tests[0]); i++) {
        Expr* e = parse_expression(tests[i]);
        ASSERT(e != NULL);
        printf("Parsed: %-15s -> ", tests[i]);
        expr_print(e);
        printf("\n");
        expr_free(e);
    }
}

void test_parse_assignments_and_equality() {
    const char* tests[] = {
        "x = 5",
        "f[x_] := x ^ 2",
        "a == b",
        "a === b",
        "a =!= b",
        "a != b",
        "a === b",
        "x = y = 2",
        "a != b",
        "a != b != c",
        "a && b",
        "a || b",
        "!a",
        "a || b && c",
        "!a || b",
        "a < b",
        "a > b",
        "a <= b",
        "a >= b",
        "a < b < c",
        "a <= b < c",
        "a; b",
        "a; b;"
    };

    for (size_t i = 0; i < sizeof(tests)/sizeof(tests[0]); i++) {
        Expr* e = parse_expression(tests[i]);
        ASSERT(e != NULL);
        printf("Parsed: %-25s -> ", tests[i]);
        expr_print(e);
        printf("\n");
        expr_free(e);
    }
}

void test_parse_part() {
    const char* tests[] = {
        "e[[1]]",
        "e[[1, 2]]",
        "e[[i, j, k, l]]",
        "f[x][[1]]",
        "(a+b)[[2]]",
        "e[[1]][[2]]"
    };

    for (size_t i = 0; i < sizeof(tests)/sizeof(tests[0]); i++) {
        Expr* e = parse_expression(tests[i]);
        ASSERT(e != NULL);
        printf("Parsed: %-25s -> ", tests[i]);
        expr_print(e);
        printf("\n");
        expr_free(e);
    }
}

void test_implicit_multiplication() {
    const char* tests[] = {
        "3x",
        "3 x",
        "x y",
        "2(x+y)",
        "(a)(b)"
    };
    
    for (size_t i = 0; i < sizeof(tests)/sizeof(tests[0]); i++) {
        Expr* e = parse_expression(tests[i]);
        ASSERT(e != NULL);
        printf("Parsed: %-15s -> ", tests[i]);
        expr_print(e);
        expr_free(e);
    }
}

void test_parse_precedence() {
    const char* input = "Times @@ Power @@@ %";
    Expr* e = parse_expression(input);
    ASSERT(e != NULL);
    char* s = expr_to_string_fullform(e);
    // Expected: Apply[Times, Apply[Power, Out[-1], List[1]]]
    // Since % parses as Out[-1]
    ASSERT_STR_EQ(s, "Apply[Times, Apply[Power, Out[-1], List[1]]]");
    free(s);
    expr_free(e);
}

void test_parse_comments() {
    const char* tests[] = {
        "123 (* this is a comment *)",
        "(* comment *) x",
        "f[ (* arg1 *) x, y (* arg2 *) ]",
        "(* multi \n line \n comment *) 42",
        "(* nested (* comment *) ! *) \"str\""
    };
    ExprType expected_types[] = {
        EXPR_INTEGER,
        EXPR_SYMBOL,
        EXPR_FUNCTION,
        EXPR_INTEGER,
        EXPR_STRING
    };

    for (size_t i = 0; i < sizeof(tests)/sizeof(tests[0]); i++) {
        Expr* e = parse_expression(tests[i]);
        ASSERT(e != NULL);
        ASSERT_MSG(e->type == expected_types[i],
                   "For input '%s': expected type %d, got %d",
                   tests[i], expected_types[i], e->type);
        expr_free(e);
    }
}



static void assert_parse_eq(const char* input, const char* expected) {
    Expr* p = parse_expression(input);
    if (!p) {
        printf("FAIL: parsed %s as NULL\n", input);
        return;
    }
    char* s = expr_to_string_fullform(p);
    if (strcmp(s, expected) != 0) {
        printf("FAIL: parsed %s\nExpected: %s\nActual:   %s\n", input, expected, s);
        assert(strcmp(s, expected) == 0);
    }
    free(s);
    expr_free(p);
}

/* Reported by Nasser, May 2026: `Timing[LUDecomposition[mat];]` printed a
 * spurious "Unexpected character: ']'" on stderr.  The parser already
 * substituted Null for the missing RHS of `;` (so the tree was correct),
 * but parse_primary was being invoked on the trailing delimiter first
 * and printing an error before unwinding.  These tests assert (a) the
 * parse tree shape and (b) that stderr stays silent. */
static void assert_parse_silent(const char* input, const char* expected) {
    fflush(stderr);
    int saved = dup(fileno(stderr));
    FILE* tmp = tmpfile();
    ASSERT(saved >= 0 && tmp != NULL);
    dup2(fileno(tmp), fileno(stderr));

    Expr* p = parse_expression(input);

    fflush(stderr);
    dup2(saved, fileno(stderr));
    close(saved);

    fseek(tmp, 0, SEEK_END);
    long sz = ftell(tmp);
    char buf[256] = {0};
    if (sz > 0) {
        fseek(tmp, 0, SEEK_SET);
        size_t r = fread(buf, 1, sizeof(buf) - 1, tmp);
        (void)r;
    }
    fclose(tmp);

    if (sz > 0) {
        fprintf(stderr, "FAIL: parse(%s) wrote to stderr: %s\n", input, buf);
        ASSERT(sz == 0);
    }
    ASSERT(p != NULL);
    char* s = expr_to_string_fullform(p);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: parse(%s)\n  Expected: %s\n  Actual:   %s\n",
                input, expected, s);
        ASSERT(strcmp(s, expected) == 0);
    }
    free(s);
    expr_free(p);
}

void test_parse_trailing_semicolon() {
    /* `expr;` followed immediately by a closing delimiter should produce
     * CompoundExpression[expr, Null] without writing to stderr.  Each
     * closer (`)`, `]`, `}`) goes through a different call site of
     * parse_expression_prec, so we cover all three. */
    assert_parse_silent("(1+1;)",        "CompoundExpression[Plus[1, 1], Null]");
    assert_parse_silent("f[x;]",         "f[CompoundExpression[x, Null]]");
    assert_parse_silent("{1+1;}",        "List[CompoundExpression[Plus[1, 1], Null]]");
    assert_parse_silent("g[x, y;]",      "g[x, CompoundExpression[y, Null]]");
    assert_parse_silent("{a, b;}",       "List[a, CompoundExpression[b, Null]]");
    /* Original repro from the bug report: AbsoluteTiming isn't a builtin
     * here so it stays symbolic, but the parse must still succeed
     * silently. */
    assert_parse_silent("Timing[LUDecomposition[mat];]",
        "Timing[CompoundExpression[LUDecomposition[mat], Null]]");
}

void test_parse_scaled_scientific() {
    /* Mathematica's `*^` scaled-scientific-notation suffix on a number
     * literal: mantissa *^ signed-integer-exponent. */
    assert_parse_eq("1.23*^4", "12300.0");
    assert_parse_eq("1.23E4", "12300.0");
    assert_parse_eq("1.23e4", "12300.0");
    assert_parse_eq("1.23*^-4", "0.000123");
    assert_parse_eq("123*^4", "1230000");
    assert_parse_eq("123*^-4", "Rational[123, 10000]");
    assert_parse_eq("100*^-2", "1");           /* reduces to integer */
    assert_parse_eq("5*^0", "5");              /* zero exponent */
    assert_parse_eq("2*^3", "2000");           /* simple integer scale */
    assert_parse_eq("-1.5*^2", "-150.0");      /* signed mantissa */
}

void test_parse_dots() {
    assert_parse_eq(".1", "0.1");
    assert_parse_eq("-.1", "-0.1");
    assert_parse_eq("+.1", "0.1");
    assert_parse_eq("x /. .1 -> 2", "ReplaceAll[x, Rule[0.1, 2]]");
    assert_parse_eq("x /. 1", "ReplaceAll[x, 1]");
    assert_parse_eq("x /.1", "Times[x, Power[0.1, -1]]");
    assert_parse_eq("1/.1", "Power[0.1, -1]");
    assert_parse_eq(".1 ..", "Repeated[0.1]");
    assert_parse_eq(".1 ...", "RepeatedNull[0.1]");
    assert_parse_eq("x / .1", "Times[x, Power[0.1, -1]]");
    assert_parse_eq("x * .1", "Times[x, 0.1]");
    assert_parse_eq("x + .1", "Plus[x, 0.1]");
    assert_parse_eq("x - .1", "Plus[x, Times[-1, 0.1]]");
    assert_parse_eq("x ^ .1", "Power[x, 0.1]");
    assert_parse_eq("x .1", "Times[x, 0.1]");
    // Let me check what `x .1` parses as in my code!
}

int main() {
    /* The parser builds expressions with the cached SYM_* symbol pointers
     * (e.g. expr_new_symbol(SYM_List)), which are only populated by
     * sym_names_init() — invoked from core_init(). Initialise the symbol
     * system before exercising the parser. */
    symtab_init(); core_init();
    TEST(test_parse_atomics);
    TEST(test_negative_numbers);
    TEST(test_implicit_multiplication);
    TEST(test_parentheses);
    TEST(test_parse_patterns);
    TEST(test_parse_functions);
    TEST(test_parse_lists);
    TEST(test_infix_operators);
    TEST(test_postfix_operators);
    TEST(test_prefix_operators);
    TEST(test_parse_assignments_and_equality);
    TEST(test_parse_part);
    TEST(test_parse_precedence);
    TEST(test_parse_comments);
    TEST(test_parse_trailing_semicolon);
    TEST(test_parse_dots);
    TEST(test_parse_scaled_scientific);

    printf("\nAll parser tests passed!\n");
    return 0;
}
