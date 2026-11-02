
#include "parse.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include <string.h>
#include <unistd.h>
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

/* Operator precedence & associativity ladder (get_operator in src/parse.c).
 *
 * get_operator is static, so it is exercised through the public parser: each
 * case pins the parse-tree shape that a specific precedence/associativity
 * relationship must produce. These lock the 1–10000 precedence ladder against
 * accidental drift (a mis-typed def.prec would regroup one of these). Every
 * expected form was cross-checked against Wolfram-Language semantics.
 *
 * Note: And and Or share a precedence level here (a historic Mathilda choice,
 * not WL-faithful), so mixed `&&`/`||` grouping is deliberately not asserted;
 * each is tested against unambiguous neighbours (comparisons, Not) instead. */
void test_operator_precedence_ladder() {
    /* --- relative precedence: the tighter operator groups on the inside --- */
    assert_parse_eq("a + b*c",      "Plus[a, Times[b, c]]");        /* * > + */
    assert_parse_eq("a*b + c",      "Plus[Times[a, b], c]");
    assert_parse_eq("a*b^c",        "Times[a, Power[b, c]]");       /* ^ > * */
    assert_parse_eq("a + b == c",   "Equal[Plus[a, b], c]");        /* + > == */
    assert_parse_eq("a && b == c",  "And[a, Equal[b, c]]");         /* == > && */
    assert_parse_eq("a == b && c",  "And[Equal[a, b], c]");
    assert_parse_eq("!a && b",      "And[Not[a], b]");              /* ! > && */
    assert_parse_eq("!a == b",      "Not[Equal[a, b]]");            /* == > ! */
    assert_parse_eq("a . b + c",    "Plus[Dot[a, b], c]");          /* . > + */
    assert_parse_eq("a*b . c",      "Times[a, Dot[b, c]]");         /* . > * */
    assert_parse_eq("2^3!",         "Power[2, Factorial[3]]");      /* ! > ^ */
    assert_parse_eq("f @ a + b",    "Plus[f[a], b]");               /* @ > + */
    assert_parse_eq("a + b // f",   "f[Plus[a, b]]");               /* + > // */
    assert_parse_eq("a -> b /; c",  "Rule[a, Condition[b, c]]");    /* /; > -> */
    assert_parse_eq("a /. b -> c",  "ReplaceAll[a, Rule[b, c]]");   /* -> > /. */
    assert_parse_eq("a ~~ b | c",   "StringExpression[a, Alternatives[b, c]]"); /* | > ~~ */
    assert_parse_eq("x = a -> b",   "Set[x, Rule[a, b]]");          /* -> > = */
    assert_parse_eq("a[[1]] + b",   "Plus[Part[a, 1], b]");         /* [[ ]] > + */
    assert_parse_eq("f[x]^2",       "Power[f[x], 2]");              /* f[x] > ^ */
    assert_parse_eq("a @* b @ c",   "Composition[a, b][c]");        /* @* > @ */
    assert_parse_eq("a ? b + c",    "Plus[PatternTest[a, b], c]");  /* ? > + */

    /* --- prefix unary minus: binds looser than Power/Times, tighter than +/- --- */
    assert_parse_eq("-a^2",         "Times[-1, Power[a, 2]]");
    assert_parse_eq("-f[x]",        "Times[-1, f[x]]");

    /* --- prefix unary plus: a one-argument Plus, as in WL (issue #101);
     *     binds tighter than binary +, looser than Times/Power --- */
    assert_parse_eq("+x",           "Plus[x]");
    assert_parse_eq("+ x",          "Plus[x]");
    assert_parse_eq("+a b",         "Plus[Times[a, b]]");
    assert_parse_eq("+a^2",         "Plus[Power[a, 2]]");
    assert_parse_eq("+x == y",      "Equal[Plus[x], y]");
    assert_parse_eq("+-x",          "Plus[Times[-1, x]]");
    assert_parse_eq("f[+x]",        "f[Plus[x]]");
    assert_parse_eq("{+x, -y}",     "List[Plus[x], Times[-1, y]]");
    assert_parse_eq("+5",           "5");                 /* signed literal */
    assert_parse_eq("++x",          "PreIncrement[x]");

    /* --- associativity --- */
    assert_parse_eq("a - b - c",    "Plus[a, Times[-1, b], Times[-1, c]]"); /* + left */
    assert_parse_eq("a^b^c",        "Power[a, Power[b, c]]");               /* ^ right */
    assert_parse_eq("a/b/c",        "Times[Times[a, Power[b, -1]], Power[c, -1]]"); /* / left */
    assert_parse_eq("a = b = c",    "Set[a, Set[b, c]]");                   /* = right */
    assert_parse_eq("a -> b -> c",  "Rule[a, Rule[b, c]]");                 /* -> right */
    assert_parse_eq("a // f // g",  "g[f[a]]");                             /* // left */
    assert_parse_eq("f @ g @ x",    "f[g[x]]");                             /* @ right */
    assert_parse_eq("a @@ b @@ c",  "Apply[a, Apply[b, c]]");              /* @@ right */
    assert_parse_eq("a && b && c",  "And[And[a, b], c]");                   /* && left (parse) */
    assert_parse_eq("a <> b <> c",  "StringJoin[a, b, c]");                 /* <> n-ary */
    assert_parse_eq("a . b . c",    "Dot[Dot[a, b], c]");                   /* . left */
    assert_parse_eq("a | b | c",    "Alternatives[a, b, c]");               /* | n-ary */
    assert_parse_eq("a ~~ b ~~ c",  "StringExpression[a, b, c]");           /* ~~ n-ary */

    /* `:` is Pattern (sym:p, binding tighter than ~~) or Optional (p:v, just
     * below Pattern), as in Mathematica (issue #98). */
    assert_parse_eq("x__ ~~ y:DigitCharacter..",
                    "StringExpression[Pattern[x, BlankSequence[]], Pattern[y, Repeated[DigitCharacter]]]");
    assert_parse_eq("x:_ ~~ \"b\"",  "StringExpression[Pattern[x, Blank[]], \"b\"]");
    assert_parse_eq("a ~~ x:_",      "StringExpression[a, Pattern[x, Blank[]]]");
    assert_parse_eq("x:a|b",         "Pattern[x, Alternatives[a, b]]");
    assert_parse_eq("x:p..",         "Pattern[x, Repeated[p]]");
    assert_parse_eq("x:_ /; c",      "Condition[Pattern[x, Blank[]], c]");
    assert_parse_eq("x_:v",          "Optional[Pattern[x, Blank[]], v]");
    assert_parse_eq("x:_:1",         "Optional[Pattern[x, Blank[]], 1]");
    assert_parse_eq("a:b:c",         "Optional[Pattern[a, b], c]");
    assert_parse_eq("a:b:c:d",       "Optional[Pattern[a, b], Pattern[c, d]]");
    assert_parse_eq("x:_ + y:_",     "Optional[Pattern[x, Plus[Blank[], y]], Blank[]]");
    assert_parse_eq("p:_ -> rhs",    "Rule[Pattern[p, Blank[]], rhs]");
    assert_parse_eq("x:_Integer:1 ~~ y",
                    "StringExpression[Optional[Pattern[x, Blank[Integer]], 1], y]");

    /* --- CompoundExpression is the loosest binary operator --- */
    assert_parse_eq("a ; b ; c",    "CompoundExpression[a, b, c]");

    /* --- one expression spanning many levels at once --- */
    assert_parse_eq("a + b c^d/e - f!",
        "Plus[a, Times[b, Times[Power[c, d], Power[e, -1]]], Times[-1, Factorial[f]]]");
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

/* Issue #102: `:` takes its LEFT operand at atom level -- the bare symbol
 * (Pattern) or blank token (Optional) just before it -- not the expression
 * built so far, so n_.*a:_g is Times[n_., a:_g], not Optional[n_.*a, _g].
 * Every expected FullForm below was produced by the reference kernel
 * (ToString[FullForm[ToExpression[src, InputForm, Hold]]]). */
void test_parse_colon_atom_level() {
    assert_parse_eq("n_.*a:_g", "Times[Optional[Pattern[n, Blank[]]], Pattern[a, Blank[g]]]");
    assert_parse_eq("n_*a:_g", "Times[Pattern[n, Blank[]], Pattern[a, Blank[g]]]");
    assert_parse_eq("x+a:_g", "Plus[x, Pattern[a, Blank[g]]]");
    assert_parse_eq("x + a:_g == b", "Plus[x, Pattern[a, Equal[Blank[g], b]]]");
    assert_parse_eq("x + a:_g + c", "Plus[x, Pattern[a, Plus[Blank[g], c]]]");
    assert_parse_eq("n a:_g", "Times[n, Pattern[a, Blank[g]]]");
    assert_parse_eq("n a:_g b", "Times[n, Pattern[a, Times[Blank[g], b]]]");
    assert_parse_eq("a:_g b", "Pattern[a, Times[Blank[g], b]]");
    assert_parse_eq("a.b:_g", "Dot[a, Pattern[b, Blank[g]]]");
    assert_parse_eq("x_:0", "Optional[Pattern[x, Blank[]], 0]");
    assert_parse_eq("a:b:c", "Optional[Pattern[a, b], c]");
    assert_parse_eq("a:b:c:d", "Optional[Pattern[a, b], Pattern[c, d]]");
    assert_parse_eq("x_:a:b", "Optional[Pattern[x, Blank[]], Pattern[a, b]]");
    assert_parse_eq("a:x_:0", "Optional[Pattern[a, Pattern[x, Blank[]]], 0]");
    assert_parse_eq("-a:_g", "Times[-1, Pattern[a, Blank[g]]]");
    assert_parse_eq("!a:_", "Not[Pattern[a, Blank[]]]");
    assert_parse_eq("2 a:_g", "Times[2, Pattern[a, Blank[g]]]");
    assert_parse_eq("b/a:_g", "Times[b, Power[Pattern[a, Blank[g]], -1]]");
    assert_parse_eq("x^a:_g", "Power[x, Pattern[a, Blank[g]]]");
    assert_parse_eq("f@a:_g", "f[Pattern[a, Blank[g]]]");
    assert_parse_eq("a:_g->1", "Rule[Pattern[a, Blank[g]], 1]");
    assert_parse_eq("a:_/;c", "Condition[Pattern[a, Blank[]], c]");
    assert_parse_eq("a:_g..", "Pattern[a, Repeated[Blank[g]]]");
    assert_parse_eq("a:b|c", "Pattern[a, Alternatives[b, c]]");
    assert_parse_eq("x|a:b", "Alternatives[x, Pattern[a, b]]");
    assert_parse_eq("x|a:b|c", "Alternatives[x, Pattern[a, Alternatives[b, c]]]");
    assert_parse_eq("a ~~ x:b", "StringExpression[a, Pattern[x, b]]");
    assert_parse_eq("x:a ~~ b", "StringExpression[Pattern[x, a], b]");
    assert_parse_eq("a:_+b", "Pattern[a, Plus[Blank[], b]]");
    assert_parse_eq("x->a:_g", "Rule[x, Pattern[a, Blank[g]]]");
    assert_parse_eq("{a:_g, b:_h}", "List[Pattern[a, Blank[g]], Pattern[b, Blank[h]]]");
    assert_parse_eq("f[a:_g, n_.*b:_h]", "f[Pattern[a, Blank[g]], Times[Optional[Pattern[n, Blank[]]], Pattern[b, Blank[h]]]]");
    assert_parse_eq("a:_g:1", "Optional[Pattern[a, Blank[g]], 1]");
    assert_parse_eq("n_.*a:_g:1", "Times[Optional[Pattern[n, Blank[]]], Optional[Pattern[a, Blank[g]], 1]]");
    assert_parse_eq("x__*a:_g", "Times[Pattern[x, BlankSequence[]], Pattern[a, Blank[g]]]");
    assert_parse_eq("a:_g?f", "Pattern[a, PatternTest[Blank[g], f]]");
    assert_parse_eq("x?f:_", "PatternTest[x, Pattern[f, Blank[]]]");
    assert_parse_eq("a:_g&&b", "Pattern[a, And[Blank[g], b]]");
    assert_parse_eq("x&&a:_g", "And[x, Pattern[a, Blank[g]]]");
    assert_parse_eq("x==a:_g", "Equal[x, Pattern[a, Blank[g]]]");
    assert_parse_eq("x>a:_g", "Greater[x, Pattern[a, Blank[g]]]");
    assert_parse_eq("a:_g;b", "CompoundExpression[Pattern[a, Blank[g]], b]");
    assert_parse_eq("a :> b", "RuleDelayed[a, b]");
    assert_parse_eq("a := b", "SetDelayed[a, b]");
    assert_parse_eq("a::b", "MessageName[a, \"b\"]");
    assert_parse_eq("x: _", "Pattern[x, Blank[]]");
    assert_parse_eq("x : _g * y", "Pattern[x, Times[Blank[g], y]]");
    assert_parse_eq("x_^n_.*a:_", "Times[Power[Pattern[x, Blank[]], Optional[Pattern[n, Blank[]]]], Pattern[a, Blank[]]]");
    assert_parse_eq("Times[a:_g]", "Times[Pattern[a, Blank[g]]]");
    assert_parse_eq("p:{___}", "Pattern[p, List[BlankNullSequence[]]]");
    assert_parse_eq("f[x:_:0]", "f[Optional[Pattern[x, Blank[]], 0]]");
    assert_parse_eq("n x_:0", "Times[n, Optional[Pattern[x, Blank[]], 0]]");
    assert_parse_eq("n_.*a:_g:1 + c", "Times[Optional[Pattern[n, Blank[]]], Optional[Pattern[a, Blank[g]], Plus[1, c]]]");
    assert_parse_eq("a:_g:1+c", "Optional[Pattern[a, Blank[g]], Plus[1, c]]");
    assert_parse_eq("x + a:b:c", "Plus[x, Optional[Pattern[a, b], c]]");
    assert_parse_eq("a:_:1 b", "Optional[Pattern[a, Blank[]], Times[1, b]]");
    assert_parse_eq("x + a:_g:1 == b", "Plus[x, Optional[Pattern[a, Blank[g]], Equal[1, b]]]");
    assert_parse_eq("x + a:_g:1 -> b", "Rule[Plus[x, Optional[Pattern[a, Blank[g]], 1]], b]");
    assert_parse_eq("x /; a:_g", "Condition[x, Pattern[a, Blank[g]]]");
    assert_parse_eq("x:a:b:c:d:e", "Optional[Pattern[x, a], Optional[Pattern[b, c], Pattern[d, e]]]");
    assert_parse_eq("a:_g:b:_h", "Optional[Pattern[a, Blank[g]], Pattern[b, Blank[h]]]");
    assert_parse_eq("2 a:_g:1", "Times[2, Optional[Pattern[a, Blank[g]], 1]]");
    assert_parse_eq("f[x] := n_*a:_g", "SetDelayed[f[x], Times[Pattern[n, Blank[]], Pattern[a, Blank[g]]]]");
    assert_parse_eq("a:_g:1?f", "Optional[Pattern[a, Blank[g]], PatternTest[1, f]]");
    assert_parse_eq("x^a:_:1", "Power[x, Optional[Pattern[a, Blank[]], 1]]");
    assert_parse_eq("{n_*a:_g, x}", "List[Times[Pattern[n, Blank[]], Pattern[a, Blank[g]]], x]");
    assert_parse_eq("a:_g + x_:0", "Optional[Pattern[a, Plus[Blank[g], Pattern[x, Blank[]]]], 0]");
    assert_parse_eq("a:f[b:c:d]", "Pattern[a, f[Optional[Pattern[b, c], d]]]");
    assert_parse_eq("a:{b:c, d:e:f}", "Pattern[a, List[Pattern[b, c], Optional[Pattern[d, e], f]]]");
    assert_parse_eq("a:(b:c)", "Pattern[a, Pattern[b, c]]");
    assert_parse_eq("a:b|c:d", "Optional[Pattern[a, Alternatives[b, c]], d]");
    assert_parse_eq("a:_?(x:_)", "Pattern[a, PatternTest[Blank[], Pattern[x, Blank[]]]]");
    assert_parse_eq("x_:a:b:c", "Optional[Pattern[x, Blank[]], Optional[Pattern[a, b], c]]");
    assert_parse_eq("x_:a:b:c:d", "Optional[Pattern[x, Blank[]], Optional[Pattern[a, b], Pattern[c, d]]]");
    assert_parse_eq("a:b..:c", "Optional[Pattern[a, Repeated[b]], c]");
    assert_parse_eq("x:_ ~~ y:_", "StringExpression[Pattern[x, Blank[]], Pattern[y, Blank[]]]");
    assert_parse_eq("x_:0 ~~ y", "StringExpression[Optional[Pattern[x, Blank[]], 0], y]");
    assert_parse_eq("a:x_+y_:1", "Optional[Pattern[a, Plus[Pattern[x, Blank[]], Pattern[y, Blank[]]]], 1]");
    assert_parse_eq("f[x_, y_:0 z]", "f[Pattern[x, Blank[]], Optional[Pattern[y, Blank[]], Times[0, z]]]");
    assert_parse_eq("f[x_, y_:0, z_:1]", "f[Pattern[x, Blank[]], Optional[Pattern[y, Blank[]], 0], Optional[Pattern[z, Blank[]], 1]]");
    assert_parse_eq("a:b->c:d", "Rule[Pattern[a, b], Pattern[c, d]]");
    assert_parse_eq("a:b:>c:d", "RuleDelayed[Pattern[a, b], Pattern[c, d]]");
    assert_parse_eq("x_:a+b:c", "Optional[Pattern[x, Blank[]], Plus[a, Pattern[b, c]]]");
    assert_parse_eq("x_h:0", "Optional[Pattern[x, Blank[h]], 0]");
    assert_parse_eq("x__:0", "Optional[Pattern[x, BlankSequence[]], 0]");
    assert_parse_eq("x___h:0", "Optional[Pattern[x, BlankNullSequence[h]], 0]");
    assert_parse_eq("_:0", "Optional[Blank[], 0]");
    assert_parse_eq("__:0", "Optional[BlankSequence[], 0]");
    assert_parse_eq("x_:0?f", "Optional[Pattern[x, Blank[]], PatternTest[0, f]]");
    assert_parse_eq("x:_:0 /; y", "Condition[Optional[Pattern[x, Blank[]], 0], y]");
    assert_parse_eq("h[a:_:0, b:_.]", "h[Optional[Pattern[a, Blank[]], 0], Pattern[b, Optional[Blank[]]]]");
    assert_parse_eq("x_:a:b|c", "Optional[Pattern[x, Blank[]], Pattern[a, Alternatives[b, c]]]");
    assert_parse_eq("a:-b:c", "Optional[Pattern[a, Times[-1, b]], c]");
    assert_parse_eq("a:!b:c", "Optional[Pattern[a, Not[b]], c]");
    assert_parse_eq("f[_, n_.*a:_g:2]", "f[Blank[], Times[Optional[Pattern[n, Blank[]]], Optional[Pattern[a, Blank[g]], 2]]]");
    assert_parse_eq("p:x_*y_", "Pattern[p, Times[Pattern[x, Blank[]], Pattern[y, Blank[]]]]");
    assert_parse_eq("p:x_*y_:1", "Optional[Pattern[p, Times[Pattern[x, Blank[]], Pattern[y, Blank[]]]], 1]");
    assert_parse_eq("x:a.b:c", "Optional[Pattern[x, Dot[a, b]], c]");
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
    assert_parse_eq("x - .1", "Plus[x, -0.1]");   /* a literal is negated in place, as the reference */
    assert_parse_eq("x ^ .1", "Power[x, 0.1]");
    assert_parse_eq("x .1", "Times[x, 0.1]");
    // Let me check what `x .1` parses as in my code!
}

/* Helper: assert that parse_next_expression yields exactly `expected_count`
 * statements from `buffer`, and that the last one's FullForm is `last`.
 * Mirrors how mathilda_run_file (Get) iterates a multi-statement file. */
static void assert_statements(const char* buffer, int expected_count,
                              const char* last) {
    const char* ptr = buffer;
    int count = 0;
    char* last_form = NULL;
    for (;;) {
        Expr* e = parse_next_expression(&ptr);
        if (!e) break;
        count++;
        free(last_form);
        last_form = expr_to_string_fullform(e);
        expr_free(e);
    }
    if (count != expected_count) {
        printf("FAIL: %s -> %d statements (expected %d)\n",
               buffer, count, expected_count);
        assert(count == expected_count);
    }
    if (last && (!last_form || strcmp(last_form, last) != 0)) {
        printf("FAIL: %s -> last %s (expected %s)\n",
               buffer, last_form ? last_form : "(none)", last);
        assert(last_form && strcmp(last_form, last) == 0);
    }
    free(last_form);
}

/* GitHub issue #20: two top-level expressions on separate lines with no
 * trailing `;` were merged into one via implicit multiplication, so
 * Get["file.m"] returned the wrong value. A bracket-depth-0 line break must
 * terminate the statement (Mathematica semantics), while line breaks inside
 * brackets stay insignificant. */
void test_parse_newline_separator() {
    /* Bare juxtaposition across a newline: two statements, not a product. */
    assert_statements("a\nb\n", 2, "b");
    /* Each statement is itself an implicit product on one line. */
    assert_statements("a b\nc d\n", 2, "Times[c, d]");
    /* Real repro: Series then Integrate — must stay two statements. */
    assert_statements("Series[Sin[x],{x,0,5}]\nIntegrate[Sin[x],x]\n", 2,
                      "Integrate[Sin[x], x]");
    /* Explicit `;` and bare newlines mix freely. */
    assert_statements("x = 3;\nx^2\n", 2, "Power[x, 2]");

    /* Line breaks inside brackets do NOT separate: one statement each. */
    assert_statements("f[a,\n b,\n c]\n", 1, "f[a, b, c]");
    assert_statements("{1,\n 2,\n 3}\n", 1, "List[1, 2, 3]");
    assert_statements("(a\n b)\n", 1, "Times[a, b]");
    /* Juxtaposition inside braces across a newline stays implicit Times. */
    assert_statements("{a\n b}\n", 1, "List[Times[a, b]]");

    /* A line ending in an operator continues onto the next line. */
    assert_statements("1 +\n2\n", 1, "Plus[1, 2]");
    assert_statements("x =\n 5\n", 1, "Set[x, 5]");
}

static void test_parentheses_break_comparison_chains(void) {
    /* A chain of comparisons folds into one variadic Inequality, and the fold used to decide by
     * INSPECTING the built subtree -- which cannot tell `a >= b == c` (a chain) from `(a >= b) == c`
     * (an Equal whose left operand happens to be a comparison). So parentheses were ignored on the
     * LEFT operand, and `(2.0 >= 2.) == (1.0 > 0.)` evaluated to `2.0 == True` instead of True.
     *
     * The right operand was never affected, being parsed as its own subexpression, which is why only
     * one side looked wrong and why the bug survived: every symmetric-looking test passed. */
    assert_eval_eq("FullForm[Hold[(a >= b) == c]]", "Hold[Equal[GreaterEqual[a, b], c]]", 0);
    assert_eval_eq("FullForm[Hold[(a > b) == (c > d)]]",
                   "Hold[Equal[Greater[a, b], Greater[c, d]]]", 0);
    /* The value form, which is what surfaced it: both must be True, not `2.0 == True`. */
    assert_eval_eq("(2.0 >= 2.) == (1.0 > 0.)", "True", 0);
    assert_eval_eq("(1.0 >= 2.) == (0.0 > 0.)", "True", 0);
    /* Nesting a whole chain in parentheses must not extend it either -- the same provenance
     * question, reached through extend_inequality rather than the fold. */
    assert_eval_eq("FullForm[Hold[(a < b < c) == d]]",
                   "Hold[Equal[Inequality[a, Less, b, Less, c], d]]", 0);
}

static void test_unparenthesised_chains_still_chain(void) {
    /* The other half: the fix must not stop real chains from folding. Without these the fix could
     * "work" by simply never chaining, which would be a different bug with the same test passing. */
    assert_eval_eq("FullForm[Hold[a >= b == c]]", "Hold[Inequality[a, GreaterEqual, b, Equal, c]]",
                   0);
    assert_eval_eq("FullForm[Hold[1 < 2 < 3]]", "Hold[Inequality[1, Less, 2, Less, 3]]", 0);
    assert_eval_eq("FullForm[Hold[a < b <= c == d > e]]",
                   "Hold[Inequality[a, Less, b, LessEqual, c, Equal, d, Greater, e]]", 0);
    /* And chain SEMANTICS are unchanged: a chain is True only when every adjacent pair holds. */
    assert_eval_eq("{1 < 2 < 3, 1 < 3 < 2, 2 == 2 == 2, 3 > 2 > 1}",
                   "{True, False, True, True}", 0);
}

static void test_variadic_comparison_heads(void) {
    /* SameQ, UnsameQ and Unequal are variadic predicates: a run of one of them
     * collects into a single call, as in Mathematica. Parenthesising still
     * nests, and a different operator does not join. */
    assert_eval_eq("FullForm[Hold[a === b === c]]", "Hold[SameQ[a, b, c]]", 0);
    assert_eval_eq("FullForm[Hold[a === b === c === d]]", "Hold[SameQ[a, b, c, d]]", 0);
    assert_eval_eq("FullForm[Hold[a =!= b =!= c]]", "Hold[UnsameQ[a, b, c]]", 0);
    assert_eval_eq("FullForm[Hold[a != b != c]]", "Hold[Unequal[a, b, c]]", 0);
    assert_eval_eq("FullForm[Hold[(a === b) === c]]", "Hold[SameQ[SameQ[a, b], c]]", 0);
    assert_eval_eq("FullForm[Hold[a === (b === c)]]", "Hold[SameQ[a, SameQ[b, c]]]", 0);
    assert_eval_eq("FullForm[Hold[a === b =!= c]]", "Hold[UnsameQ[SameQ[a, b], c]]", 0);
    /* Unequal mixed with an ordering head joins the Inequality. */
    assert_eval_eq("FullForm[Hold[a < b != c]]", "Hold[Inequality[a, Less, b, Unequal, c]]", 0);
    assert_eval_eq("FullForm[Hold[a != b < c]]", "Hold[Inequality[a, Unequal, b, Less, c]]", 0);
    assert_eval_eq("FullForm[Hold[a != b != c < d]]",
                   "Hold[Inequality[a, Unequal, b, Unequal, c, Less, d]]", 0);
    assert_eval_eq("FullForm[Hold[(a != b) < c]]", "Hold[Less[Unequal[a, b], c]]", 0);
}

/* A string literal, and a symbol name, must survive the lexer at ANY length.
 * Both buffers were a fixed char[256] whose overflow was silently dropped, so
 * every literal past 255 characters lost its tail with no diagnostic: a wrong
 * answer outright, and, for the common case of a literal holding an expression
 * (ToExpression, a generated .m file), a bogus syntax error at the cut point.
 * A string built at RUN time (StringJoin) was never affected, which is how it
 * survived -- so the lengths below must come from source literals. */
void test_long_string_literal_not_truncated() {
    char src[2600];
    const size_t lens[] = {200, 255, 256, 257, 300, 1000, 2048};
    for (size_t k = 0; k < sizeof(lens) / sizeof(lens[0]); k++) {
        size_t n = lens[k];
        src[0] = '"';
        for (size_t i = 0; i < n; i++) src[1 + i] = (char)('a' + (i % 26));
        src[1 + n] = '"';
        src[2 + n] = '\0';
        Expr* e = parse_expression(src);
        ASSERT(e != NULL);
        ASSERT(e->type == EXPR_STRING);
        ASSERT_MSG(strlen(e->data.string) == n,
                   "string literal of %zu chars parsed to %zu",
                   n, strlen(e->data.string));
        /* the TAIL must be intact, not merely the length */
        ASSERT(e->data.string[n - 1] == (char)('a' + ((n - 1) % 26)));
        expr_free(e);
    }
    /* escapes past the old boundary still decode */
    {
        char esc[700];
        size_t p = 0;
        esc[p++] = '"';
        for (size_t i = 0; i < 300; i++) esc[p++] = 'z';
        esc[p++] = '\\'; esc[p++] = 'n';
        for (size_t i = 0; i < 300; i++) esc[p++] = 'z';
        esc[p++] = '"'; esc[p] = '\0';
        Expr* e = parse_expression(esc);
        ASSERT(e != NULL && e->type == EXPR_STRING);
        ASSERT(strlen(e->data.string) == 601);
        ASSERT(e->data.string[300] == '\n');
        expr_free(e);
    }
    /* a long context-qualified symbol name must not alias onto a shorter one */
    {
        char a[600], b[600];
        size_t p = 0;
        for (size_t i = 0; i < 300; i++) { a[p] = 'q'; b[p] = 'q'; p++; }
        a[p] = 'A'; b[p] = 'B'; p++;
        a[p] = '\0'; b[p] = '\0';
        Expr* ea = parse_expression(a);
        Expr* eb = parse_expression(b);
        ASSERT(ea != NULL && eb != NULL);
        ASSERT(ea->type == EXPR_SYMBOL && eb->type == EXPR_SYMBOL);
        ASSERT_MSG(!expr_eq(ea, eb),
                   "two 301-char symbol names differing in the last character "
                   "were parsed as the same symbol");
        expr_free(ea); expr_free(eb);
    }
    printf("PASS: long string literals and symbol names are not truncated\n");
}

/* A '[' with no head in front of it is a syntax error. It used to parse to
 * an EXPR_FUNCTION with a NULL head, and `f[x,[x]]` then segfaulted in
 * evaluate_step. Part (`g[[x]]`) and curried calls (`f[x][y]`) are
 * postfix on a real head and must be unaffected. */
void test_parse_headless_bracket() {
    const char* bad[] = { "f[x,[x]]", "[x]", "{[x]}", "x + [y]", "a; [y]" };
    for (size_t i = 0; i < sizeof(bad)/sizeof(bad[0]); i++) {
        Expr* e = parse_expression(bad[i]);
        ASSERT_MSG(e == NULL, "expected syntax error for '%s'", bad[i]);
        if (e) expr_free(e);
    }
    const char* good[] = { "g[[x]]", "f[x][y]", "f[x, Sqrt[x]]" };
    for (size_t i = 0; i < sizeof(good)/sizeof(good[0]); i++) {
        Expr* e = parse_expression(good[i]);
        ASSERT_MSG(e != NULL && e->type == EXPR_FUNCTION
                   && e->data.function.head != NULL,
                   "expected a headed function for '%s'", good[i]);
        expr_free(e);
    }
}

/* ToExpression reads a string as a sequence of line-oriented statements --
 * issue #108, where "a = 2\n b = 3" was a syntax error. Every expected form
 * below was taken from the reference implementation (ToString[FullForm[
 * ToExpression[s, InputForm, Hold]]]): a line break ends a complete statement
 * whatever the next line starts with, an incomplete one continues, a `;` at
 * the end of a line takes an implicit Null, and each blank or comment-only
 * line between statements is a Null statement. */
/* A backslash at the end of a line inside a string literal is a line
 * continuation: both characters are dropped.  An escaped backslash before a
 * newline is a backslash and a real newline; before CRLF only the CR goes. */
void test_string_line_continuation() {
    /* C source "ab\\\ncd" is the Mathilda source "ab\<LF>cd". */
    assert_eval_eq("\"ab\\\ncd\"", "\"abcd\"", 0);
    assert_eval_eq("StringLength[\"ab\\\ncd\"]", "4", 0);
    assert_eval_eq("\"a\\\n\\\nb\"", "\"ab\"", 0);
    assert_eval_eq("\"tab\\\n   cd\"", "\"tab   cd\"", 0);
    assert_eval_eq("\"x\\\\\ny\"", "\"x\\\\\\ny\"", 0);
    assert_eval_eq("\"ab\\\r\ncd\"", "\"ab\\ncd\"", 0);
    /* Through ToExpression too. */
    assert_eval_eq("ToExpression[\"\\\"p\\\\\nq\\\"\"]", "\"pq\"", 0);
}

/* The statement reader (parse_next_expression, used by Get and -file): the
 * FullForm of each statement it yields from `src`, joined by " | ". */
static void statements_of(const char* src, char* out, size_t cap) {
    out[0] = '\0';
    const char* p = src;
    for (int k = 0; k < 16; k++) {
        Expr* e = parse_next_expression(&p);
        if (!e) break;
        char* t = expr_to_string_fullform(e);
        if (k) strncat(out, " | ", cap - strlen(out) - 1);
        strncat(out, t, cap - strlen(out) - 1);
        free(t);
        expr_free(e);
    }
}

static void assert_statement_list(const char* src, const char* expected) {
    char got[2048];
    statements_of(src, got, sizeof got);
    if (strcmp(got, expected) != 0) {
        printf("FAIL statements: %s\n  Expected: %s\n  Actual:   %s\n", src, expected, got);
        ASSERT(0);
    }
}

/* A top-level `;` ends a statement unless what follows can only continue an
 * expression: then the `;` is CompoundExpression with an implicit Null and that
 * compound is the operator's left operand.  Every row's expected shape is the
 * reference reader's (a reference CompoundExpression[a, s2] of two statements
 * is the two statements here, which evaluate identically). */
void test_statement_continuation_after_semicolon() {
    static const char* const cases[][2] = {
        { "a; // f", "f[CompoundExpression[a, Null]]" },
        { "a; /. x -> 1", "ReplaceAll[CompoundExpression[a, Null], Rule[x, 1]]" },
        { "a; //. x -> 1", "ReplaceRepeated[CompoundExpression[a, Null], Rule[x, 1]]" },
        { "a; -> b", "Rule[CompoundExpression[a, Null], b]" },
        { "a; :> b", "RuleDelayed[CompoundExpression[a, Null], b]" },
        { "a; == b", "Equal[CompoundExpression[a, Null], b]" },
        { "a; === b", "SameQ[CompoundExpression[a, Null], b]" },
        { "a; && b", "And[CompoundExpression[a, Null], b]" },
        { "a; || b", "Or[CompoundExpression[a, Null], b]" },
        { "a; &", "Function[CompoundExpression[a, Null]]" },
        { "a; @@ b", "Apply[CompoundExpression[a, Null], b]" },
        { "a; /@ b", "Map[CompoundExpression[a, Null], b]" },
        { "a; @ b", "CompoundExpression[a, Null][b]" },
        { "a; ^ b", "Power[CompoundExpression[a, Null], b]" },
        { "a; * b", "Times[CompoundExpression[a, Null], b]" },
        { "a; > b", "Greater[CompoundExpression[a, Null], b]" },
        { "a; >= b", "GreaterEqual[CompoundExpression[a, Null], b]" },
        { "a; < b", "Less[CompoundExpression[a, Null], b]" },
        { "a; <= b", "LessEqual[CompoundExpression[a, Null], b]" },
        { "a; <> b", "StringJoin[CompoundExpression[a, Null], b]" },
        { "a; ~~ b", "StringExpression[CompoundExpression[a, Null], b]" },
        { "a; | b", "Alternatives[CompoundExpression[a, Null], b]" },
        { "a; := b", "SetDelayed[CompoundExpression[a, Null], b]" },
        { "a; = b", "Set[CompoundExpression[a, Null], b]" },
        { "a; /; b /; c", "Condition[Condition[CompoundExpression[a, Null], b], c]" },
        { "a; /; b; /; c", "Condition[CompoundExpression[Condition[CompoundExpression[a, Null], b], Null], c]" },
        { "a; -b", "a | Times[-1, b]" },
        { "a; !b", "a | Not[b]" },
        { "a;;b", "Span[a, b]" },
        { "a; <|x->1|>", "a | Association[Rule[x, 1]]" },
        { "a; ::x", "MessageName[CompoundExpression[a, Null], \"x\"]" },
        { "a; .5", "a | 0.5" },
        { "a; +b", "a | Plus[b]" },
        { "a; b", "a | b" },
    };
    for (size_t i = 0; i < sizeof cases / sizeof cases[0]; i++)
        assert_statement_list(cases[i][0], cases[i][1]);
    /* The FeynCalc shape (SharedTools.m): a definition, `;/;` guard, and the
     * statement after it still read as its own statement. */
    assert_statement_list("FCFactorOut[expr_, pref_] := pref f[expr/pref];/; !MemberQ[{List}, Head[expr]];\nnext = 1",
                      "Condition[CompoundExpression[SetDelayed[FCFactorOut[Pattern[expr, Blank[]], "
                      "Pattern[pref, Blank[]]], Times[pref, f[Times[expr, Power[pref, -1]]]]], Null], "
                      "Not[MemberQ[List[List], Head[expr]]]] | Set[next, 1]");
    assert_statement_list("a; -= 2", "SubtractFrom[CompoundExpression[a, Null], 2]");
    assert_statement_list("a; += 2", "AddTo[CompoundExpression[a, Null], 2]");
    /* Statements still split at `;` -- the context-prologue guarantee. */
    assert_statement_list("BeginPackage[\"P`\"]; foo; EndPackage[]",
                      "BeginPackage[\"P`\"] | foo | EndPackage[]");
    assert_statement_list("a := b; c := d", "SetDelayed[a, b] | SetDelayed[c, d]");
    assert_statement_list("a;\n\nb;", "a | b");
    /* A guard on the next line is a syntax error (the reference rejects it too):
     * the reader yields the first statement, then stops. */
    assert_statement_list("f[x_] := g[x];\n/; x > 0", "SetDelayed[f[Pattern[x, Blank[]]], g[x]]");
}

/* The same continuation inside one expression (ToExpression, a bracketed
 * body): the infix `;` must decide "implicit Null" BEFORE parsing an operand.
 * Parsing first let a prefix minus swallow the `-` of `->` / `-=`, so
 * `a; -> c` came back as Greater[..., c] and `a; -= c` as Set[..., c]. */
void test_infix_semicolon_before_operator() {
    assert_parse_eq("a; -> c", "Rule[CompoundExpression[a, Null], c]");
    assert_parse_eq("a; -= c", "SubtractFrom[CompoundExpression[a, Null], c]");
    assert_parse_eq("a; += c", "AddTo[CompoundExpression[a, Null], c]");
    assert_parse_eq("a; != c", "Unequal[CompoundExpression[a, Null], c]");
    assert_parse_eq("a; >= c", "GreaterEqual[CompoundExpression[a, Null], c]");
    assert_parse_eq("a; || c", "Or[CompoundExpression[a, Null], c]");
    assert_parse_eq("a; /; c", "Condition[CompoundExpression[a, Null], c]");
    assert_parse_eq("f[a; -> c]", "f[Rule[CompoundExpression[a, Null], c]]");
    assert_parse_eq("(a; b; -> c)", "Rule[CompoundExpression[a, b, Null], c]");
    assert_parse_eq("a; b; c", "CompoundExpression[a, b, c]");
    assert_parse_eq("(a; b); c", "CompoundExpression[CompoundExpression[a, b], c]");
    assert_parse_eq("a; (b; c)", "CompoundExpression[a, CompoundExpression[b, c]]");
    assert_parse_eq("a; b; c;", "CompoundExpression[a, b, c, Null]");
    /* Operands that DO begin an expression are still the next statement. */
    assert_parse_eq("a; -c", "CompoundExpression[a, Times[-1, c]]");
    assert_parse_eq("a; !c", "CompoundExpression[a, Not[c]]");
    assert_parse_eq("a; -1", "CompoundExpression[a, -1]");
    assert_parse_eq("a; <|x -> 1|>", "CompoundExpression[a, Association[Rule[x, 1]]]");
    assert_statement_list("a; != c", "Unequal[CompoundExpression[a, Null], c]");
}

/* `<<file` is Get["file"]; the name after `<<`, `>>` and `>>>` is a file
 * token, not an expression: letters, digits, `` !$*-./:?\_`~ `` and balanced
 * [...], ended by whitespace or ; , ) ] } -- any other character is a syntax
 * error.  Get binds tighter than everything: a<<f is Times[a, Get["f"]]. */
static void assert_parse_fails(const char* input);

void test_get_put_file_tokens(void) {
    assert_parse_eq("<<f", "Get[\"f\"]");
    assert_parse_eq("<<a/b.m", "Get[\"a/b.m\"]");
    assert_parse_eq("<<Foo`", "Get[\"Foo`\"]");
    assert_parse_eq("<<  f.m  ", "Get[\"f.m\"]");
    assert_parse_eq("<< \"x y.m\"", "Get[\"x y.m\"]");
    assert_parse_eq("<<~/x-1.m", "Get[\"~/x-1.m\"]");
    assert_parse_eq("<<f[x]", "Get[\"f[x]\"]");
    assert_parse_eq("<<f//g", "Get[\"f//g\"]");
    assert_parse_eq("<<f!g", "Get[\"f!g\"]");
    assert_parse_eq("<<f;g", "CompoundExpression[Get[\"f\"], g]");
    assert_parse_eq("{<<f, 2}", "List[Get[\"f\"], 2]");
    assert_parse_eq("(<<f)", "Get[\"f\"]");
    assert_parse_eq("f[<<g]", "f[Get[\"g\"]]");
    assert_parse_eq("a<<f", "Times[a, Get[\"f\"]]");
    assert_parse_eq("x=<<f", "Set[x, Get[\"f\"]]");
    assert_parse_eq("<<\"f\"[x]", "Get[\"f\"][x]");
    assert_parse_eq("<<\"a\" <> \"b\"", "StringJoin[Get[\"a\"], \"b\"]");
    assert_parse_fails("<<f+1");
    assert_parse_fails("<<f[g");
    assert_parse_fails("<<f)g");
    assert_parse_fails("<<");
    /* >> and >>> read their file name with the same token rules. */
    assert_parse_eq("x>>f[1]", "Put[x, \"f[1]\"]");
    assert_parse_eq("x>>f*g", "Put[x, \"f*g\"]");
    assert_parse_eq("x>>f`g", "Put[x, \"f`g\"]");
    assert_parse_eq("x>>>f?g", "PutAppend[x, \"f?g\"]");
    assert_parse_eq("x>>f;g", "CompoundExpression[Put[x, \"f\"], g]");
    assert_parse_eq("{x>>f, 1}", "List[Put[x, \"f\"], 1]");
    assert_parse_fails("x>>f+g");
    /* The statement reader: <<file at the start of a line is one statement. */
    assert_statement_list("<<f\ng", "Get[\"f\"] | g");
    assert_statement_list("a; <<f", "a | Get[\"f\"]");
}

/* Get stops at a syntax error -- reporting Syntax::sntx through the message
 * funnel (Quiet / Check see it) with the line number -- and returns Null; the
 * parser's own diagnostics are not printed on top. */
void test_get_syntax_error_reported() {
    char path[128], src[512];
    snprintf(path, sizeof path, "/tmp/Mathilda_parse_get_%d.m", (int)getpid());
    FILE* fp = fopen(path, "w");
    ASSERT(fp != NULL);
    fputs("gsA = 1;\ngsB = (2 + ;\ngsC = 3;\n", fp);
    fclose(fp);
    snprintf(src, sizeof src, "{Quiet[Get[\"%s\"]], gsA, gsB, gsC}", path);
    assert_eval_eq(src, "{Null, 1, gsB, gsC}", 0);
    snprintf(src, sizeof src, "Check[Get[\"%s\"], \"caught\"]", path);
    assert_eval_eq(src, "\"caught\"", 0);
    /* An unfinished expression at end of file is an error too, not "complete". */
    fp = fopen(path, "w");
    ASSERT(fp != NULL);
    fputs("gsD = 1;\ngsE = (1 + \n", fp);
    fclose(fp);
    snprintf(src, sizeof src, "{Quiet[Get[\"%s\"]], gsD, gsE}", path);
    assert_eval_eq(src, "{Null, 1, gsE}", 0);
    /* The `;/;` form loads, and the statements after it run. */
    fp = fopen(path, "w");
    ASSERT(fp != NULL);
    fputs("gsF[x_] := gsG[x];/; x > 0;\ngsH = 5;\n", fp);
    fclose(fp);
    snprintf(src, sizeof src, "Get[\"%s\"]; {gsF[2], gsH}", path);
    assert_eval_eq(src, "{gsG[2], 5}", 0);
    /* With the guard on the next line: the first definition stands, then a
     * syntax error stops the read (gsJ stays unset). */
    fp = fopen(path, "w");
    ASSERT(fp != NULL);
    fputs("gsI[x_] := gsG[x];\n/; x > 0;\ngsJ = 5;\n", fp);
    fclose(fp);
    snprintf(src, sizeof src, "{Quiet[Get[\"%s\"]], gsI[2], gsJ}", path);
    assert_eval_eq(src, "{Null, gsG[2], gsJ}", 0);
    unlink(path);
    /* ToExpression / SyntaxQ parse quietly: Quiet leaks nothing, Check sees the error. */
    assert_eval_eq("Quiet[ToExpression[\"(1+\"]]", "$Failed", 0);
    assert_eval_eq("Check[ToExpression[\"(1+\"], \"caught\"]", "\"caught\"", 0);
    assert_eval_eq("SyntaxQ[\"(1+\"]", "False", 0);
}

/* ToExpression evaluates statement by statement: Begin["ctx`"] on one line
 * governs the symbols of the lines after it; a syntax error further on gives
 * $Failed after the earlier statements ran; the Hold form only parses. */
void test_toexpression_sequencing() {
    assert_eval_eq("ToExpression[\"Begin[\\\"TQ`Private`\\\"]\\ntqz = 1\\nEnd[]\"]; "
                   "{Names[\"TQ`Private`*\"], Names[\"Global`tqz\"]}", "{{\"TQ`Private`tqz\"}, {}}", 0);
    assert_eval_eq("tqr = Quiet[ToExpression[\"tqa = 1\\ntqc = (2 + \\ntqd = 3\"]]; {tqr, tqa, tqc, tqd}",
                   "{$Failed, 1, tqc, tqd}", 0);
    assert_eval_eq("ToExpression[\"tqf = 1\\ntqg = 2\", InputForm, Hold]", "Hold[tqf = 1, tqg = 2]", 0);
    assert_eval_eq("{ToExpression[\"1+1\"], ToExpression[\"\"], ToExpression[\"q1;\\nq2\"]}", "{2, Null, q2}", 0);
    /* An abort stops the read. */
    assert_eval_eq("CheckAbort[ToExpression[\"tqh = 1\\nAbort[]\\ntqi = 2\"], \"ab\"]; {tqh, tqi}", "{1, tqi}", 0);
}

/* An unclosed "(" is a syntax error, not an expression ending at end of input. */
void test_unclosed_paren_rejected() {
    assert_eval_eq("Quiet[ToExpression[\"(2 + x\"]]", "$Failed", 0);
    assert_eval_eq("SyntaxQ[\"(2 + x\"]", "False", 0);
    assert_eval_eq("SyntaxQ[\"(2 + x)\"]", "True", 0);
    assert_eval_eq("Quiet[ToExpression[\"f[(1)\"]]", "$Failed", 0);
}

/* Reading a statement creates its symbols (outside kernel module loads): the
 * loader idiom v = Names["Global`*"] lists v itself. */
void test_reading_creates_symbols() {
    assert_eval_eq("ToExpression[\"Global`rcsv = Names[\\\"Global`*\\\"]\"]; MemberQ[rcsv, \"rcsv\"]", "True", 0);
    assert_eval_eq("ToExpression[\"Function[rcsArg, rcsArg + 1][1]\"]; Names[\"Global`rcsArg\"]",
                   "{\"rcsArg\"}", 0);
}

void test_toexpression_line_statements() {
    assert_eval_eq("Quiet[ToExpression[\"a = 2\\n b = 3\", InputForm, Hold]]", "Hold[Set[a, 2], Set[b, 3]]", 1);
    assert_eval_eq("Quiet[ToExpression[\"1+\\n2\", InputForm, Hold]]", "Hold[Plus[1, 2]]", 1);
    assert_eval_eq("Quiet[ToExpression[\"f[1,\\n2]\\n\\n g[3]\\n\", InputForm, Hold]]", "Hold[f[1, 2], Null, g[3]]", 1);
    assert_eval_eq("Quiet[ToExpression[\"a;\\n b\", InputForm, Hold]]", "Hold[CompoundExpression[a, Null], b]", 1);
    assert_eval_eq("Quiet[ToExpression[\"a b\", InputForm, Hold]]", "Hold[Times[a, b]]", 1);
    assert_eval_eq("Quiet[ToExpression[\"(a\\n b)\", InputForm, Hold]]", "Hold[Times[a, b]]", 1);
    assert_eval_eq("Quiet[ToExpression[\"{a\\n b}\", InputForm, Hold]]", "Hold[List[Times[a, b]]]", 1);
    assert_eval_eq("Quiet[ToExpression[\"\\n\", InputForm, Hold]]", "Hold[Null]", 1);
    assert_eval_eq("Quiet[ToExpression[\"\", InputForm, Hold]]", "Hold[Null]", 1);
    assert_eval_eq("Quiet[ToExpression[\"a\\n\", InputForm, Hold]]", "Hold[a]", 1);
    assert_eval_eq("Quiet[ToExpression[\"a\\n b +\", InputForm, Hold]]", "$Failed", 1);
    assert_eval_eq("Quiet[ToExpression[\"a // f\\n b\", InputForm, Hold]]", "Hold[f[a], b]", 1);
    assert_eval_eq("Quiet[ToExpression[\"{\\\"a\\n b\\\"}\", InputForm, Hold]]", "Hold[List[\"a\\n b\"]]", 1);
    assert_eval_eq("Quiet[ToExpression[\"a (* c\\n *) b\", InputForm, Hold]]", "Hold[Times[a, b]]", 1);
    assert_eval_eq("Quiet[ToExpression[\"a\\n (* c *) b\", InputForm, Hold]]", "Hold[a, b]", 1);
    assert_eval_eq("Quiet[ToExpression[\"a\\n\\nb\", InputForm, Hold]]", "Hold[a, Null, b]", 1);
    assert_eval_eq("Quiet[ToExpression[\"a\\n\\n\\nb\", InputForm, Hold]]", "Hold[a, Null, Null, b]", 1);
    assert_eval_eq("Quiet[ToExpression[\"a\\n \\nb\", InputForm, Hold]]", "Hold[a, Null, b]", 1);
    assert_eval_eq("Quiet[ToExpression[\"\\na\\nb\", InputForm, Hold]]", "Hold[Null, a, b]", 1);
    assert_eval_eq("Quiet[ToExpression[\"\\n\\na\", InputForm, Hold]]", "Hold[Null, Null, a]", 1);
    assert_eval_eq("Quiet[ToExpression[\"a\\n\\n\", InputForm, Hold]]", "Hold[a, Null]", 1);
    assert_eval_eq("Quiet[ToExpression[\"a\\r\\nb\", InputForm, Hold]]", "Hold[a, Null, b]", 1);
    assert_eval_eq("Quiet[ToExpression[\"a ;b\\nc\", InputForm, Hold]]", "Hold[CompoundExpression[a, b], c]", 1);
    assert_eval_eq("Quiet[ToExpression[\"a\\n)\", InputForm, Hold]]", "$Failed", 1);
    assert_eval_eq("Quiet[ToExpression[\"a\\n\\n(* x *)\\n\\nb\", InputForm, Hold]]", "Hold[a, Null, Null, Null, b]", 1);
    assert_eval_eq("Quiet[ToExpression[\"a\\n(* x *)\\nb\", InputForm, Hold]]", "Hold[a, Null, b]", 1);
    assert_eval_eq("Quiet[ToExpression[\"a\\n  \\n  \\nb\", InputForm, Hold]]", "Hold[a, Null, Null, b]", 1);
    assert_eval_eq("Quiet[ToExpression[\"a\\n+b\", InputForm, Hold]]", "Hold[a, Plus[b]]", 1);
    assert_eval_eq("Quiet[ToExpression[\"a\\n-b\", InputForm, Hold]]", "Hold[a, Times[-1, b]]", 1);
    assert_eval_eq("Quiet[ToExpression[\"a +\\n\\nb\", InputForm, Hold]]", "Hold[Plus[a, b]]", 1);
    assert_eval_eq("Quiet[ToExpression[\"a;\\n\\nb\", InputForm, Hold]]", "Hold[CompoundExpression[a, Null], Null, b]", 1);
    assert_eval_eq("Quiet[ToExpression[\"f\\n[x]\", InputForm, Hold]]", "$Failed", 1);
    assert_eval_eq("Quiet[ToExpression[\"a\\n// f\", InputForm, Hold]]", "$Failed", 1);
    assert_eval_eq("Quiet[ToExpression[\"a\\n;b\", InputForm, Hold]]", "$Failed", 1);
    assert_eval_eq("Quiet[ToExpression[\"a;\", InputForm, Hold]]", "Hold[CompoundExpression[a, Null]]", 1);
    assert_eval_eq("Quiet[ToExpression[\"a;\\n\", InputForm, Hold]]", "Hold[CompoundExpression[a, Null]]", 1);
    assert_eval_eq("Quiet[ToExpression[\";\", InputForm, Hold]]", "$Failed", 1);
    assert_eval_eq("Quiet[ToExpression[\"a = \\n 2\", InputForm, Hold]]", "Hold[Set[a, 2]]", 1);
    assert_eval_eq("Quiet[ToExpression[\"  a  \\n  b  \", InputForm, Hold]]", "Hold[a, b]", 1);
    assert_eval_eq("Quiet[ToExpression[\"a\\tb\", InputForm, Hold]]", "Hold[Times[a, b]]", 1);
    assert_eval_eq("Quiet[ToExpression[\"x=1;y=2\\nz\", InputForm, Hold]]", "Hold[CompoundExpression[Set[x, 1], Set[y, 2]], z]", 1);
    assert_eval_eq("Quiet[ToExpression[\"a (* \\n *)\\n b\", InputForm, Hold]]", "Hold[a, b]", 1);
    assert_eval_eq("Quiet[ToExpression[\"a; (* c *)\\n b\", InputForm, Hold]]", "Hold[CompoundExpression[a, Null], b]", 1);
    assert_eval_eq("Quiet[ToExpression[\"a;;\\nb\", InputForm, Hold]]", "Hold[Span[a, All], b]", 1);
    assert_eval_eq("Quiet[ToExpression[\"f[x\\n] b\", InputForm, Hold]]", "Hold[Times[f[x], b]]", 1);
    assert_eval_eq("Quiet[ToExpression[\"{a}\\n b\", InputForm, Hold]]", "Hold[List[a], b]", 1);
    assert_eval_eq("Quiet[ToExpression[\"f[x_] :=\\n x^2\\ng[2]\", InputForm, Hold]]", "Hold[SetDelayed[f[Pattern[x, Blank[]]], Power[x, 2]], g[2]]", 1);
    assert_eval_eq("Quiet[ToExpression[\"If[a,\\n b,\\n c]\\nd\", InputForm, Hold]]", "Hold[If[a, b, c], d]", 1);
    assert_eval_eq("Quiet[ToExpression[\"(* only *)\", InputForm, Hold]]", "Hold[Null]", 1);
    assert_eval_eq("Quiet[ToExpression[\"a\\rb\", InputForm, Hold]]", "Hold[a, b]", 1);
    assert_eval_eq("Quiet[ToExpression[\"x = 1 (* t *)\\n\", InputForm, Hold]]", "Hold[Set[x, 1]]", 1);
    assert_eval_eq("Quiet[ToExpression[\"<|a -> 1,\\n b -> 2|>\\nc\", InputForm, Hold]]", "Hold[Association[Rule[a, 1], Rule[b, 2]], c]", 1);
    assert_eval_eq("Quiet[ToExpression[\"a[[1\\n]]\\nb\", InputForm, Hold]]", "Hold[Part[a, 1], b]", 1);
    assert_eval_eq("Quiet[ToExpression[\"-x\\n!y\", InputForm, Hold]]", "Hold[Times[-1, x], Not[y]]", 1);
    /* StandardForm: line breaks are whitespace. */
    assert_eval_eq("ToExpression[\"a\\nb\", StandardForm, Hold]", "Hold[Times[a, b]]", 1);
    assert_eval_eq("ToExpression[\"a\\n\\nb\\n\", StandardForm, Hold]", "Hold[Times[a, b]]", 1);
    /* Without a head each statement is evaluated in turn; the last result
     * is returned. */
    assert_eval_eq("ToExpression[\"tx0 = 2\\n ty0 = tx0 + 1\"]", "3", 0);
    assert_eval_eq("{tx0, ty0}", "{2, 3}", 0);
    assert_eval_eq("ToExpression[\"a\\n\\n\"]", "Null", 0);
    assert_eval_eq("ToExpression[\"a\\nb\", InputForm, f]", "f[a, b]", 1);
    assert_eval_eq("ToExpression[{\"a\\nb\", \"c\"}, InputForm, Hold]",
                   "List[Hold[a, b], Hold[c]]", 1);
    /* A syntax error is a message Check can see. */
    assert_eval_eq("Check[ToExpression[\"a\\n// f\"], err, ToExpression::sntx]", "err", 0);
}

/* Empty argument / element slots are Null, as in Mathematica: `f[,x]` is
 * f[Null, x] (the leading comma used to be swallowed, giving f[x]), and
 * `f[x,]` / `f[x,,y]` used to be syntax errors. `f[]` / `{}` stay empty. */
void test_parse_empty_slots_are_null() {
    struct { const char* input; const char* fullform; } tests[] = {
        {"f[,x]",   "f[Null, x]"},
        {"f[x,]",   "f[x, Null]"},
        {"f[x,,y]", "f[x, Null, y]"},
        {"f[,]",    "f[Null, Null]"},
        {"f[]",     "f[]"},
        {"f[ ]",    "f[]"},
        {"{,1}",    "List[Null, 1]"},
        {"{1,}",    "List[1, Null]"},
        {"{a,,b}",  "List[a, Null, b]"},
        {"{}",      "List[]"},
    };
    for (size_t i = 0; i < sizeof(tests)/sizeof(tests[0]); i++) {
        Expr* e = parse_expression(tests[i].input);
        ASSERT_MSG(e != NULL, "parse failed for '%s'", tests[i].input);
        char* got = expr_to_string_fullform(e);
        ASSERT_MSG(strcmp(got, tests[i].fullform) == 0,
                   "'%s': expected %s, got %s", tests[i].input,
                   tests[i].fullform, got);
        free(got);
        expr_free(e);
    }
    /* A missing closer is still a syntax error. */
    const char* bad[] = { "f[x y", "{1, 2", "f[x,", "{," };
    for (size_t i = 0; i < sizeof(bad)/sizeof(bad[0]); i++) {
        Expr* e = parse_expression(bad[i]);
        ASSERT_MSG(e == NULL, "expected syntax error for '%s'", bad[i]);
        if (e) expr_free(e);
    }
}

/* The input must be rejected as a syntax error (parse_expression -> NULL). */
static void assert_parse_fails(const char* input) {
    Expr* p = parse_expression(input);
    if (p) {
        char* s = expr_to_string_fullform(p);
        printf("FAIL: expected a syntax error for %s, parsed as %s\n", input, s);
        free(s);
        expr_free(p);
        assert(0);
    }
}

/* `..` / `...` (Repeated / RepeatedNull) tokenization -- issue #100.
 *
 * `..` is a single token that wins over a decimal point by longest munch, so
 * `0..` is Repeated[0], not the real `0.` followed by a stray `.`.  Every
 * expected form here was taken from Mathematica 14 (ToExpression[s,
 * InputForm, Hold] // FullForm), including the syntax errors: WL forms the
 * Optional token `_.` only from a BARE single blank, so `_...` is (_.).. and
 * `_..` does not parse, while a headed or sequence blank takes the dots as
 * the repeat operator. */
void test_parse_repeated_tokens() {
    /* Integers: the issue itself. */
    assert_parse_eq("0..",          "Repeated[0]");
    assert_parse_eq("0...",         "RepeatedNull[0]");
    assert_parse_eq("-0..",         "Repeated[0]");
    assert_parse_eq("-1..",         "Repeated[-1]");
    assert_parse_eq("12..",         "Repeated[12]");
    assert_parse_eq("123456789012345678901234567890..",
                    "Repeated[123456789012345678901234567890]");
    assert_parse_eq("0.....",       "Repeated[RepeatedNull[0]]");
    assert_parse_eq("1*^2..",       "Repeated[100]");
    /* Reals keep their decimal point; only a dot followed by a dot is the
     * operator. */
    assert_parse_eq("1.5..",        "Repeated[1.5]");
    assert_parse_eq("1.5...",       "RepeatedNull[1.5]");
    assert_parse_eq("1. ..",        "Repeated[1.0]");
    assert_parse_eq(".1..",         "Repeated[0.1]");
    assert_parse_eq("1.2",          "1.2");
    assert_parse_eq("1. 2",         "Times[1.0, 2]");
    /* Precision suffix: strtod must not swallow the first dot of `..`. */
    {
        Expr* p = parse_expression("2.5`10..");
        assert(p && p->type == EXPR_FUNCTION);
        char* h = expr_to_string(p->data.function.head);
        assert(strcmp(h, "Repeated") == 0);
        free(h);
        expr_free(p);
    }
    /* A repeated number inside larger expressions. */
    assert_parse_eq("{0..}",        "List[Repeated[0]]");
    assert_parse_eq("f[0.., 1...]", "f[Repeated[0], RepeatedNull[1]]");
    assert_parse_eq("0..|1",        "Alternatives[Repeated[0], 1]");
    assert_parse_eq("x:0..",        "Pattern[x, Repeated[0]]");
    assert_parse_eq("0..->a",       "Rule[Repeated[0], a]");
    assert_parse_eq("0.. :> 1",     "RuleDelayed[Repeated[0], 1]");
    assert_parse_eq("x_:0..",       "Optional[Pattern[x, Blank[]], Repeated[0]]");
    assert_parse_eq("1..2",         "Times[Repeated[1], 2]");
    assert_parse_eq("1...2",        "Times[RepeatedNull[1], 2]");
    assert_parse_eq("1..^2",        "Power[Repeated[1], 2]");
    assert_parse_eq("-3^2",         "Times[-1, Power[3, 2]]");
    assert_parse_eq("{0,0} /. {0..} -> z",
                    "ReplaceAll[List[0, 0], Rule[List[Repeated[0]], z]]");
    /* Bare single blanks: `_.` is Optional and binds first. */
    assert_parse_eq("_.",           "Optional[Blank[]]");
    assert_parse_eq("x_.",          "Optional[Pattern[x, Blank[]]]");
    assert_parse_eq("_...",         "Repeated[Optional[Blank[]]]");
    assert_parse_eq("x_...",        "Repeated[Optional[Pattern[x, Blank[]]]]");
    assert_parse_eq("_....",        "RepeatedNull[Optional[Blank[]]]");
    assert_parse_eq("x_. ..",       "Repeated[Optional[Pattern[x, Blank[]]]]");
    assert_parse_eq("_. ..",        "Repeated[Optional[Blank[]]]");
    /* Headed and sequence blanks: the dots are the repeat operator. */
    assert_parse_eq("_Integer..",   "Repeated[Blank[Integer]]");
    assert_parse_eq("_Integer...",  "RepeatedNull[Blank[Integer]]");
    assert_parse_eq("x_Integer..",  "Repeated[Pattern[x, Blank[Integer]]]");
    assert_parse_eq("x_Integer...", "RepeatedNull[Pattern[x, Blank[Integer]]]");
    assert_parse_eq("_h..",         "Repeated[Blank[h]]");
    assert_parse_eq("__..",         "Repeated[BlankSequence[]]");
    assert_parse_eq("___...",       "RepeatedNull[BlankNullSequence[]]");
    assert_parse_eq("x__..",        "Repeated[Pattern[x, BlankSequence[]]]");
    assert_parse_eq("x___..",       "Repeated[Pattern[x, BlankNullSequence[]]]");
    assert_parse_eq("x__h..",       "Repeated[Pattern[x, BlankSequence[h]]]");
    assert_parse_eq("__h...",       "RepeatedNull[BlankSequence[h]]");
    assert_parse_eq("{_Integer..}", "List[Repeated[Blank[Integer]]]");
    assert_parse_eq("f[x_Real.., y___]",
                    "f[Repeated[Pattern[x, Blank[Real]]], Pattern[y, BlankNullSequence[]]]");
    /* Other operands, and Dot left intact. */
    assert_parse_eq("a..",          "Repeated[a]");
    assert_parse_eq("a...",         "RepeatedNull[a]");
    assert_parse_eq("(a|b)..",      "Repeated[Alternatives[a, b]]");
    assert_parse_eq("{a,b}..",      "Repeated[List[a, b]]");
    assert_parse_eq("\"a\"..",      "Repeated[\"a\"]");
    assert_parse_eq("#..",          "Repeated[Slot[1]]");
    assert_parse_eq("a . b",        "Dot[a, b]");
    assert_parse_eq("a.b.c",        "Dot[Dot[a, b], c]");
    /* Syntax errors in WL too. */
    assert_parse_fails("_..");
    assert_parse_fails("x_..");
    assert_parse_fails("0....");
    assert_parse_fails("a....");
}

/* Printed repeat patterns must read back as the same expression. */
void test_repeated_print_roundtrip() {
    const char* cases[] = {
        "Repeated[0]", "RepeatedNull[0]", "Repeated[1.]", "RepeatedNull[2.5]",
        "Repeated[-1]", "Repeated[_]", "RepeatedNull[_Integer]",
        "Repeated[Optional[x_]]", "Repeated[x__]", "Repeated[a | b]",
        "Repeated[PatternSequence[a, b]]", "Repeated[RepeatedNull[0]]", NULL
    };
    for (int i = 0; cases[i]; i++) {
        Expr* e = parse_expression(cases[i]);
        assert(e);
        char* printed = expr_to_string(e);
        Expr* back = parse_expression(printed);
        if (!back || !expr_eq(e, back)) {
            printf("FAIL: %s printed as %s, which does not read back\n", cases[i], printed);
            assert(0);
        }
        free(printed);
        expr_free(back);
        expr_free(e);
    }
}


/* Named characters in source code (src/parse.c src_nonascii / u8_char_class).
 * Every expectation was checked against the reference parser
 * (ToExpression[s, InputForm, Hold]); the reference prints a non-ASCII name as
 * \[Mu] in FullForm where Mathilda prints the character itself.  A non-ASCII
 * letter is stored as UTF-8, so the escape, the hex escape and the raw
 * character are the same symbol. */
#define MU "\xce\xbc"
void test_named_characters_in_source() {
    /* Letters form identifiers, whatever the spelling. */
    assert_parse_eq("\\[Mu]+1",            "Plus[" MU ", 1]");
    assert_parse_eq(MU "+1",               "Plus[" MU ", 1]");
    assert_parse_eq("\\:03bc+1",           "Plus[" MU ", 1]");
    assert_parse_eq("a\\[Mu]1",            "a" MU "1");
    assert_parse_eq("\\[Mu]$",             MU "$");
    assert_parse_eq("\\[Alpha] \\[Beta]",  "Times[\xce\xb1, \xce\xb2]");
    assert_parse_eq("1\\[Mu]",             "Times[1, " MU "]");
    assert_parse_eq("f[\\[Mu]_, x_\\[Mu]]",
                    "f[Pattern[" MU ", Blank[]], Pattern[x, Blank[" MU "]]]");
    assert_parse_eq("\\[CapitalSigma][x]", "\xce\xa3[x]");
    /* A constant's letter alone is the constant; inside a name it is a letter. */
    assert_parse_eq("\\[Pi]",              "Pi");
    assert_parse_eq("\xcf\x80",            "Pi");
    assert_parse_eq("2\\[Pi]",             "Times[2, Pi]");
    assert_parse_eq("a\\[Pi]",             "a\xcf\x80");
    assert_parse_eq("\\[ImaginaryI]",      "I");
    assert_parse_eq("\\[Infinity]",        "Infinity");
    assert_parse_eq("\\[ExponentialE]",    "E");
    assert_parse_eq("\\[Degree]",          "Degree");
    /* Operators with an ASCII spelling reuse it. */
    assert_parse_eq("a\\[Rule]b",          "Rule[a, b]");
    assert_parse_eq("a\\[Rule]b\\[Rule]c", "Rule[a, Rule[b, c]]");
    assert_parse_eq("a\\[RuleDelayed]b",   "RuleDelayed[a, b]");
    assert_parse_eq("\\[Theta]\\[Rule]\\[Phi]", "Rule[\xce\xb8, \xcf\x95]");
    assert_parse_eq("a\\[Equal]b",         "Equal[a, b]");
    assert_parse_eq("a\\[NotEqual]b",      "Unequal[a, b]");
    assert_parse_eq("a < b \\[LessEqual] c", "Inequality[a, Less, b, LessEqual, c]");
    assert_parse_eq("a\\[Times]b",         "Times[a, b]");
    assert_parse_eq("a\\[And]b",           "And[a, b]");
    assert_parse_eq("\\[Not]a",            "Not[a]");
    /* Operators without one: precedence and grouping from the reference. */
    assert_parse_eq("a\\[CenterDot]b",     "CenterDot[a, b]");
    assert_parse_eq("a\\[CenterDot]b\\[CenterDot]c", "CenterDot[a, b, c]");
    assert_parse_eq("(a\\[CenterDot]b)\\[CenterDot]c", "CenterDot[CenterDot[a, b], c]");
    assert_parse_eq("a + b \\[CenterDot] c", "Plus[a, CenterDot[b, c]]");
    assert_parse_eq("a\\[CenterDot]b^2",   "CenterDot[a, Power[b, 2]]");
    assert_parse_eq("a \\[Cross] b \\[CenterDot] c", "CenterDot[Cross[a, b], c]");
    assert_parse_eq("a\\[Perpendicular]b\\[Perpendicular]c",
                    "Perpendicular[Perpendicular[a, b], c]");
    assert_parse_eq("x\\[Element]y",       "Element[x, y]");
    assert_parse_eq("a \\[Element] b \\[Element] c", "Element[a, b, c]");
    assert_parse_eq("a\\[LeftRightArrow]b", "LeftRightArrow[a, b]");
    assert_parse_eq("a\\[Function]b",      "Function[a, b]");
    assert_parse_eq("\\[Del]f",            "Del[f]");
    assert_parse_eq("\\[Sqrt]x",           "Sqrt[x]");
    assert_parse_eq("a\\[Transpose]",      "Transpose[a]");
    assert_parse_eq("x^2\\[Transpose]",    "Power[x, Transpose[2]]");
    /* Whitespace characters. */
    assert_parse_eq("a\\[InvisibleSpace]b", "Times[a, b]");
    assert_parse_eq("a\\[VeryThinSpace]b", "Times[a, b]");
    assert_parse_eq("{a \\[IndentingNewLine] b}", "List[Times[a, b]]");
    assert_eval_eq("Quiet[ToExpression[\"a \\[IndentingNewLine] b\", InputForm, Hold]]",
                   "Hold[a, b]", 1);
    /* Through ToExpression (the FeynCalc test files' path), and contexts. */
    assert_eval_eq("ToExpression[\"{\\[Mu], \\[Nu]} /. \\[Mu] \\[Rule] 1\"]",
                   "{1, \xce\xbd}", 0);
    assert_eval_eq("Global`\\[Mu] === " MU, "True", 0);
    assert_eval_eq("Context[\\[Mu]]", "\"Global`\"", 0);
    assert_eval_eq("ToString[a \\[Rule] \\[Mu], InputForm]", "\"a -> " MU "\"", 0);
    /* Not source characters at all. */
    ASSERT(parse_expression("\\[PartialD]x") == NULL);
}
#undef MU

/* Flat operator chains flatten only what the chain itself built (src/parse.c
 * left_flat_head / chain_head); a parenthesised or explicitly written Plus /
 * Times, or a prefix minus, stays one operand.  Checked against the reference. */
void test_flat_chain_provenance() {
    assert_parse_eq("a b c",          "Times[a, b, c]");
    assert_parse_eq("(a b) c",        "Times[Times[a, b], c]");
    assert_parse_eq("Times[a, b] c",  "Times[Times[a, b], c]");
    assert_parse_eq("Times[a, b]*c",  "Times[Times[a, b], c]");
    assert_parse_eq("a (b c)",        "Times[a, Times[b, c]]");
    assert_parse_eq("(a + b) + c",    "Plus[Plus[a, b], c]");
    assert_parse_eq("Plus[a, b] + c", "Plus[Plus[a, b], c]");
    assert_parse_eq("a + b - c",      "Plus[a, b, Times[-1, c]]");
    assert_parse_eq("a - b + c",      "Plus[a, Times[-1, b], c]");
    assert_parse_eq("(a + b) - c",    "Plus[Plus[a, b], Times[-1, c]]");
    assert_parse_eq("(a<>b)<>c",      "StringJoin[StringJoin[a, b], c]");
    assert_parse_eq("(a**b)**c",      "NonCommutativeMultiply[NonCommutativeMultiply[a, b], c]");
    /* prefix minus binds tighter than Times and is not flattened into it */
    assert_parse_eq("-a b",           "Times[Times[-1, a], b]");
    assert_parse_eq("-x y z",         "Times[Times[-1, x], y, z]");
    assert_parse_eq("-(x + y) z",     "Times[Times[-1, Plus[x, y]], z]");
    assert_parse_eq("-(x + y)^2",     "Times[-1, Power[Plus[x, y], 2]]");
    /* binary minus: the -1 joins a bare product operand only */
    assert_parse_eq("a - b c",        "Plus[a, Times[-1, b, c]]");
    assert_parse_eq("a - 2 b",        "Plus[a, Times[-1, 2, b]]");
    assert_parse_eq("a - (b c)",      "Plus[a, Times[-1, Times[b, c]]]");
    assert_parse_eq("a - (x + y) z",  "Plus[a, Times[-1, Plus[x, y], z]]");
    assert_parse_eq("a - b/c",        "Plus[a, Times[-1, Times[b, Power[c, -1]]]]");
    /* An e/E exponent needs digits after it (a Mathilda extension; the
     * reference has only *^): 2ep is 2 ep and 3E is 3 E, as it reads them. */
    assert_parse_eq("2ep",            "Times[2, ep]");
    assert_parse_eq("3E",             "Times[3, E]");
    assert_parse_eq("2e-x",           "Plus[Times[2, e], Times[-1, x]]");
    assert_parse_eq("2e5",            "200000.0");
    assert_parse_eq("1.5e-3",         "0.0015");
}

/* Linear box syntax: \!\( ... \) groups with \*Box[...] items read as the
 * expressions they display, both in source and inside a string handed to
 * ToExpression (FeynCalc's tests write powers and fractions this way). The
 * control characters in a string are the reference's private-use code points.
 * Every expected value is the reference's. */
static void test_linear_box_syntax(void) {
    assert_eval_eq("ToCharacterCode[\"\\!\\(\\)\\*\"]", "List[63425, 63433, 63424, 63432]", 1);
    assert_eval_eq("StringLength[\"\\(\\)\"]", "2", 1);
    assert_eval_eq("ToExpression[\"\\!\\(\\*SuperscriptBox[\\(mb\\), \\(2\\)]\\)-x\"]",
                   "Plus[Power[mb, 2], Times[-1, x]]", 1);
    assert_eval_eq("ToExpression[\"a + \\!\\(\\*FractionBox[\\(x + 1\\), \\(y\\)]\\) + "
                   "\\!\\(\\*SqrtBox[\\(z \\(-1\\)\\)]\\) + \\!\\(\\*SubscriptBox[\\(p\\), \\(2\\)]\\)\"]",
                   "Plus[a, Times[Plus[1, x], Power[y, -1]], Power[Times[-1, z], Rational[1, 2]], Subscript[p, 2]]", 1);
    assert_eval_eq("ToExpression[\"\\!\\(\\*SuperscriptBox[\\(mb\\), \\(2\\)] - \\*SqrtBox[\\(4\\)]\\)\"]",
                   "Plus[-2, Power[mb, 2]]", 1);
    /* juxtaposed groups multiply */
    assert_eval_eq("ToExpression[\"\\!\\(\\*SuperscriptBox[\\(mb\\), \\(2\\)]\\) \\!\\(\\*SuperscriptBox[\\(la\\), \\(2\\)]\\)\"]",
                   "Times[Power[la, 2], Power[mb, 2]]", 1);
    /* nested boxes; RowBox of strings; a bare group */
    assert_eval_eq("ToExpression[\"\\!\\(\\*SuperscriptBox[\\(a\\), \\(\\*FractionBox[\\(1\\), \\(2\\)]\\)]\\)\"]",
                   "Power[a, Rational[1, 2]]", 1);
    assert_eval_eq("ToExpression[\"\\!\\(\\*RowBox[{\\\"a\\\", \\\"+\\\", \\\"b\\\"}]\\)\"]", "Plus[a, b]", 1);
    assert_eval_eq("{ToExpression[\"\\!\\(x + y\\)\"], ToExpression[\"\\(x\\)\"]}", "List[Plus[x, y], x]", 1);
    /* the same syntax directly in source */
    assert_eval_eq("\\!\\(\\*SuperscriptBox[\\(q\\), \\(3\\)]\\) + 1", "Plus[1, Power[q, 3]]", 1);
    /* "\\ " is a space token in box text (FeynCalc's PaVeDenner-ID7) */
    assert_eval_eq("ToExpression[\"\\!\\(\\*FractionBox[\\(1\\), \\(3  p^\\\\ 2\\)]\\)\"]",
                   "Times[Rational[1, 3], Power[p, -2]]", 1);
    assert_eval_eq("ToExpression[\"\\!\\(a\\\\ b\\)\"]", "Times[a, b]", 1);
}

/* ToExpression of a non-string, as the reference treats it. */
static void test_toexpression_nonstring(void) {
    assert_eval_eq("{ToExpression[4], ToExpression[x], ToExpression[1/2], ToExpression[{1, \"2+2\"}]}",
                   "List[4, x, Rational[1, 2], List[1, 4]]", 1);
    assert_eval_eq("Quiet[ToExpression[f[x]]]", "$Failed", 1);
    assert_eval_eq("Quiet[ToExpression[3, InputForm]]", "ToExpression[3, InputForm]", 1);
}

int main() {
    /* The parser builds expressions with the cached SYM_* symbol pointers
     * (e.g. expr_new_symbol(SYM_List)), which are only populated by
     * sym_names_init() — invoked from core_init(). Initialise the symbol
     * system before exercising the parser. */
    symtab_init(); core_init();
    TEST(test_parse_atomics);
    TEST(test_variadic_comparison_heads);
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
    TEST(test_operator_precedence_ladder);
    TEST(test_parse_comments);
    TEST(test_parse_trailing_semicolon);
    TEST(test_parse_newline_separator);
    TEST(test_parse_dots);
    TEST(test_parse_repeated_tokens);
    TEST(test_repeated_print_roundtrip);
    TEST(test_parse_scaled_scientific);
    TEST(test_parentheses_break_comparison_chains);
    TEST(test_unparenthesised_chains_still_chain);
    TEST(test_long_string_literal_not_truncated);
    TEST(test_parse_headless_bracket);
    TEST(test_parse_empty_slots_are_null);
    TEST(test_parse_colon_atom_level);
    TEST(test_toexpression_line_statements);
    TEST(test_string_line_continuation);
    TEST(test_statement_continuation_after_semicolon);
    TEST(test_infix_semicolon_before_operator);
    TEST(test_get_put_file_tokens);
    TEST(test_get_syntax_error_reported);
    TEST(test_toexpression_sequencing);
    TEST(test_unclosed_paren_rejected);
    TEST(test_reading_creates_symbols);
    TEST(test_named_characters_in_source);
    TEST(test_flat_chain_provenance);
    TEST(test_linear_box_syntax);
    TEST(test_toexpression_nonstring);

    printf("\nAll parser tests passed!\n");
    return 0;
}
