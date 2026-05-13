#include "test_utils.h"
#include "match.h"
#include "expr.h"
#include "parse.h"
#include "symtab.h"
#include "core.h"

static Expr* parse_or_die(const char* s) {
    Expr* e = parse_expression(s);
    if (!e) {
        fprintf(stderr, "Failed to parse: %s\n", s);
        exit(1);
    }
    return e;
}

void test_match_basic() {
    MatchEnv* env = env_new();
    Expr* e1 = parse_or_die("123");
    Expr* p1 = parse_or_die("123");
    ASSERT(match(e1, p1, env) == true);
    expr_free(e1); expr_free(p1);

    Expr* e2 = parse_or_die("x");
    Expr* p2 = parse_or_die("y");
    ASSERT(match(e2, p2, env) == false);
    expr_free(e2); expr_free(p2);

    Expr* e3 = parse_or_die("f[1]");
    Expr* p3 = parse_or_die("f[1]");
    ASSERT(match(e3, p3, env) == true);
    expr_free(e3); expr_free(p3);
    
    env_free(env);
}

void test_match_blank() {
    MatchEnv* env = env_new();
    Expr* e1 = parse_or_die("f[1, g[2]]");
    Expr* p1 = parse_or_die("Blank[]");
    ASSERT(match(e1, p1, env) == true);
    expr_free(e1); expr_free(p1);

    Expr* e2 = parse_or_die("f[1, 2]");
    Expr* p2 = parse_or_die("f[Blank[], Blank[]]");
    ASSERT(match(e2, p2, env) == true);
    expr_free(e2); expr_free(p2);

    Expr* e3 = parse_or_die("f[1, 2]");
    Expr* p3 = parse_or_die("f[Blank[], Blank[], Blank[]]");
    ASSERT(match(e3, p3, env) == false);
    expr_free(e3); expr_free(p3);
    
    env_free(env);
}

void test_match_pattern() {
    MatchEnv* env = env_new();
    Expr* e1 = parse_or_die("f[42]");
    Expr* p1 = parse_or_die("f[Pattern[x, Blank[]]]"); // f[x_]
    ASSERT(match(e1, p1, env) == true);
    
    Expr* val = env_get(env, "x");
    ASSERT(val != NULL);
    ASSERT(val->type == EXPR_INTEGER && val->data.integer == 42);
    
    expr_free(e1); expr_free(p1);
    env_free(env);
}

void test_match_pattern_same_symbol() {
    MatchEnv* env = env_new();
    Expr* e1 = parse_or_die("f[42, 42]");
    Expr* p1 = parse_or_die("f[Pattern[x, Blank[]], Pattern[x, Blank[]]]"); // f[x_, x_]
    ASSERT(match(e1, p1, env) == true);
    expr_free(e1); expr_free(p1);

    env_free(env);
    
    env = env_new();
    Expr* e2 = parse_or_die("f[42, 43]");
    Expr* p2 = parse_or_die("f[Pattern[x, Blank[]], Pattern[x, Blank[]]]");
    ASSERT(match(e2, p2, env) == false); // Should fail because x is bound to 42, then compared to 43
    expr_free(e2); expr_free(p2);
    
    env_free(env);
}

void test_replace_bindings() {
    MatchEnv* env = env_new();
    Expr* e1 = parse_or_die("f[1]");
    Expr* p1 = parse_or_die("f[Pattern[x, Blank[]]]");
    ASSERT(match(e1, p1, env) == true);
    
    Expr* rhs = parse_or_die("g[x, x]");
    Expr* res = replace_bindings(rhs, env);
    
    Expr* expected = parse_or_die("g[1, 1]");
    ASSERT(expr_eq(res, expected));
    
    expr_free(e1); expr_free(p1); expr_free(rhs); expr_free(res); expr_free(expected);
    env_free(env);
}

void test_match_deeply_nested() {
    MatchEnv* env = env_new();
    Expr* e = parse_or_die("f[g[h[i[j[k[42]]]]]]");
    Expr* p = parse_or_die("f[g[h[i[j[Pattern[x, Blank[]]]]]]]");
    ASSERT(match(e, p, env) == true);
    
    Expr* val = env_get(env, "x");
    ASSERT(val != NULL);
    
    Expr* expected = parse_or_die("k[42]");
    ASSERT(expr_eq(val, expected));
    
    expr_free(e); expr_free(p); expr_free(expected);
    env_free(env);
}

void test_match_multiple_bindings() {
    MatchEnv* env = env_new();
    Expr* e = parse_or_die("f[1, g[2, 3], h[4, 5, 6]]");
    Expr* p = parse_or_die("f[Pattern[a, Blank[]], g[Pattern[b, Blank[]], Pattern[c, Blank[]]], Pattern[d, Blank[]]]");
    ASSERT(match(e, p, env) == true);
    
    Expr* val_a = env_get(env, "a");
    ASSERT(val_a != NULL && val_a->type == EXPR_INTEGER && val_a->data.integer == 1);
    
    Expr* val_b = env_get(env, "b");
    ASSERT(val_b != NULL && val_b->type == EXPR_INTEGER && val_b->data.integer == 2);
    
    Expr* val_c = env_get(env, "c");
    ASSERT(val_c != NULL && val_c->type == EXPR_INTEGER && val_c->data.integer == 3);
    
    Expr* val_d = env_get(env, "d");
    ASSERT(val_d != NULL);
    Expr* expected_d = parse_or_die("h[4, 5, 6]");
    ASSERT(expr_eq(val_d, expected_d));
    
    expr_free(e); expr_free(p); expr_free(expected_d);
    env_free(env);
}

void test_match_patterntest() {
    symtab_init();
    core_init();
    MatchEnv* env = env_new();
    
    // _?NumberQ matching 123
    Expr* e1 = parse_or_die("123");
    Expr* p1 = parse_or_die("_?NumberQ");
    ASSERT(match(e1, p1, env) == true);
    expr_free(e1); expr_free(p1);

    // _?NumberQ matching x
    Expr* e2 = parse_or_die("x");
    Expr* p2 = parse_or_die("_?NumberQ");
    ASSERT(match(e2, p2, env) == false);
    expr_free(e2); expr_free(p2);
    
    // _?(Mod[#, 2] == 0 &) matching 4
    Expr* e3 = parse_or_die("4");
    Expr* p3 = parse_or_die("_?(Mod[#, 2] == 0 &)");
    ASSERT(match(e3, p3, env) == true);
    expr_free(e3); expr_free(p3);

    // _?(Mod[#, 2] == 0 &) matching 3
    Expr* e4 = parse_or_die("3");
    Expr* p4 = parse_or_die("_?(Mod[#, 2] == 0 &)");
    ASSERT(match(e4, p4, env) == false);
    expr_free(e4); expr_free(p4);

    // MatchQ[Range[10], {___, _?PrimeQ, ___}]
    Expr* e5 = parse_or_die("MatchQ[Range[10], {___, _?PrimeQ, ___}]");
    Expr* res5 = evaluate(e5);
    ASSERT(res5->type == EXPR_SYMBOL && strcmp(res5->data.symbol, "True") == 0);
    expr_free(e5); expr_free(res5);

    // MatchQ[Range[10], {___, _?(Mod[#, 2] == 0 &), ___}]
    Expr* e6 = parse_or_die("MatchQ[Range[10], {___, _?(Mod[#, 2] == 0 &), ___}]");
    Expr* res6 = evaluate(e6);
    ASSERT(res6->type == EXPR_SYMBOL && strcmp(res6->data.symbol, "True") == 0);
    expr_free(e6); expr_free(res6);

    env_free(env);
}

int main() {
    symtab_init();
    TEST(test_match_basic);
    TEST(test_match_blank);
    TEST(test_match_pattern);
    TEST(test_match_pattern_same_symbol);
    TEST(test_replace_bindings);
    TEST(test_match_deeply_nested);
    TEST(test_match_multiple_bindings);
    TEST(test_match_patterntest);
    printf("All match tests passed!\n");
    return 0;
}
