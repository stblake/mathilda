
#include "test_utils.h"
#include "match.h"
#include "expr.h"
#include "parse.h"
#include "symtab.h"
#include "core.h"
#include <stdio.h>

static Expr* parse_or_die(const char* s) {
    Expr* e = parse_expression(s);
    if (!e) {
        fprintf(stderr, "Failed to parse: %s\n", s);
        exit(1);
    }
    return e;
}

static bool test_match_q(const char* expr_str, const char* pat_str) {
    MatchEnv* env = env_new();
    Expr* e = parse_or_die(expr_str);
    Expr* p = parse_or_die(pat_str);
    bool res = match(e, p, env);
    expr_free(e);
    expr_free(p);
    env_free(env);
    return res;
}

void test_match_extensive_literals() {
    ASSERT(test_match_q("1", "1") == true);
    ASSERT(test_match_q("1", "2") == false);
    ASSERT(test_match_q("3.14", "3.14") == true);
    ASSERT(test_match_q("3.14", "3.15") == false);
    ASSERT(test_match_q("x", "x") == true);
    ASSERT(test_match_q("x", "y") == false);
    ASSERT(test_match_q("\"str\"", "\"str\"") == true);
    ASSERT(test_match_q("\"str\"", "\"other\"") == false);
}

void test_match_extensive_blanks() {
    ASSERT(test_match_q("1", "_") == true);
    ASSERT(test_match_q("f[x]", "_") == true);
    ASSERT(test_match_q("f[x]", "f[_]") == true);
    ASSERT(test_match_q("f[x, y]", "f[_]") == false);
    ASSERT(test_match_q("f[x, y]", "f[_, _]") == true);
}

void test_match_extensive_typed_blanks() {
    ASSERT(test_match_q("1", "_Integer") == true);
    ASSERT(test_match_q("1.1", "_Integer") == false);
    ASSERT(test_match_q("1.1", "_Real") == true);
    ASSERT(test_match_q("x", "_Symbol") == true);
    ASSERT(test_match_q("f[x]", "_f") == true);
    ASSERT(test_match_q("g[x]", "_f") == false);
}

void test_match_extensive_named_patterns() {
    MatchEnv* env = env_new();
    Expr* e = parse_or_die("f[1, 2]");
    Expr* p = parse_or_die("f[x_, y_]");
    ASSERT(match(e, p, env) == true);
    ASSERT(expr_eq(env_get(env, "x"), parse_or_die("1")));
    ASSERT(expr_eq(env_get(env, "y"), parse_or_die("2")));
    expr_free(e); expr_free(p); env_free(env);

    env = env_new();
    e = parse_or_die("f[1, 1]");
    p = parse_or_die("f[x_, x_]");
    ASSERT(match(e, p, env) == true);
    expr_free(e); expr_free(p); env_free(env);

    env = env_new();
    e = parse_or_die("f[1, 2]");
    p = parse_or_die("f[x_, x_]");
    ASSERT(match(e, p, env) == false);
    expr_free(e); expr_free(p); env_free(env);
}

void test_match_extensive_alternatives() {
    ASSERT(test_match_q("1", "1 | 2") == true);
    ASSERT(test_match_q("2", "1 | 2") == true);
    ASSERT(test_match_q("3", "1 | 2") == false);
    ASSERT(test_match_q("f[1]", "f[1 | 2]") == true);
    ASSERT(test_match_q("f[3]", "f[1 | 2]") == false);
}

void test_match_extensive_conditions() {
    // Condition[x_, x === 1]
    ASSERT(test_match_q("1", "x_ /; x === 1") == true);
    ASSERT(test_match_q("2", "x_ /; x === 1") == false);
}

void test_match_extensive_sequences() {
    ASSERT(test_match_q("f[1, 2, 3]", "f[__]") == true);
    ASSERT(test_match_q("f[1, 2, 3]", "f[___]") == true);
    ASSERT(test_match_q("f[]", "f[___]") == true);
    ASSERT(test_match_q("f[]", "f[__]") == false);
    ASSERT(test_match_q("f[1, 2, 3]", "f[_, __]") == true);
    ASSERT(test_match_q("f[1, 2, 3]", "f[__, _]") == true);
    
    // Typed sequences
    ASSERT(test_match_q("f[1, 2, 3]", "f[__Integer]") == true);
    ASSERT(test_match_q("f[1, 2.5, 3]", "f[__Integer]") == false);

    // Sequences in Condition
    ASSERT(test_match_q("f[1, 2]", "f[x__] /; Length[{x}] === 2") == true);
    
    // Longest and Shortest constraints (tested via ReplaceAll)
    Expr* e1 = evaluate(parse_or_die("ReplaceAll[{a,b,c,d,e,f}, {x__, y__} :> {x}]"));
    Expr* r1 = parse_or_die("{a}");
    ASSERT(expr_eq(e1, r1) == true);
    expr_free(e1); expr_free(r1);

    Expr* e2 = evaluate(parse_or_die("ReplaceAll[{a,b,c,d,e,f}, {x__, Longest[y__]} :> {x}]"));
    Expr* r2 = parse_or_die("{a}");
    ASSERT(expr_eq(e2, r2) == true);
    expr_free(e2); expr_free(r2);

    Expr* e3 = evaluate(parse_or_die("ReplaceAll[{a,b,c,d,e,f}, {Longest[x__], y__} :> {y}]"));
    Expr* r3 = parse_or_die("{f}");
    ASSERT(expr_eq(e3, r3) == true);
    expr_free(e3); expr_free(r3);

    Expr* e4 = evaluate(parse_or_die("ReplaceAll[{a,b,c,d,e,f}, {Shortest[x__], y__} :> {x}]"));
    Expr* r4 = parse_or_die("{a}");
    ASSERT(expr_eq(e4, r4) == true);
    expr_free(e4); expr_free(r4);
}

void test_match_extensive_named_sequences() {
    printf("  Starting test_match_extensive_named_sequences...\n");
    MatchEnv* env = env_new();
    Expr* e = parse_or_die("f[1, 2, 3]");
    Expr* p = parse_or_die("f[x_, y__]");
    printf("  Parsed e and p\n");
    bool m = match(e, p, env);
    printf("  Match result: %d\n", m);
    ASSERT(m == true);
    
    printf("  Checking x...\n");
    Expr* expected_x = parse_or_die("1");
    ASSERT(expr_eq(env_get(env, "x"), expected_x));
    expr_free(expected_x);
    
    printf("  Checking y...\n");
    // y should be Sequence[2, 3]
    Expr* val_y = env_get(env, "y");
    ASSERT(val_y != NULL);
    printf("  y found in env, type=%d\n", val_y->type);
    ASSERT(val_y->type == EXPR_FUNCTION && strcmp(val_y->data.function.head->data.symbol.name, "Sequence") == 0);
    ASSERT(val_y->data.function.arg_count == 2);
    
    Expr* expected_y = parse_or_die("Sequence[2, 3]");
    ASSERT(expr_eq(val_y, expected_y));
    expr_free(expected_y);
    
    printf("  Cleaning up...\n");
    expr_free(e); expr_free(p); env_free(env);
    printf("  Done.\n");
}

void test_match_extensive_repeated() {
    ASSERT(test_match_q("{f[a,b], f[a,c], f[a,d]}", "{f[x_, _]..}") == true);
    ASSERT(test_match_q("{f[a,b], f[a,c], f[a,d]}", "{f[_, x_]..}") == false);
    ASSERT(test_match_q("{4,5,6}", "{Repeated[x_Integer]}") == false);
    ASSERT(test_match_q("{4,4,4}", "{Repeated[x_Integer]}") == true);
    ASSERT(test_match_q("{4,5,6}", "{Repeated[_Integer]}") == true);
}

void test_match_extensive_repeatednull() {
    ASSERT(test_match_q("{f[a,b], f[a,c], f[a,d]}", "{f[x_, _]...}") == true);
    ASSERT(test_match_q("{f[a,b], f[a,c], f[a,d]}", "{f[_, x_]...}") == false);
    
    ASSERT(test_match_q("{4,5,6}", "{RepeatedNull[x_Integer]}") == false);
    ASSERT(test_match_q("{4,4,4}", "{RepeatedNull[x_Integer]}") == true);
    ASSERT(test_match_q("{4,5,6}", "{RepeatedNull[_Integer]}") == true);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_match_extensive_literals);
    TEST(test_match_extensive_blanks);
    TEST(test_match_extensive_typed_blanks);
    TEST(test_match_extensive_named_patterns);
    TEST(test_match_extensive_sequences);
    TEST(test_match_extensive_named_sequences);
    TEST(test_match_extensive_alternatives);
    TEST(test_match_extensive_conditions);
    TEST(test_match_extensive_repeated);
    TEST(test_match_extensive_repeatednull);
    printf("All extensive match tests passed!\n");
    return 0;
}
