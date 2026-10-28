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
    ASSERT(res5->type == EXPR_SYMBOL && strcmp(res5->data.symbol.name, "True") == 0);
    expr_free(e5); expr_free(res5);

    // MatchQ[Range[10], {___, _?(Mod[#, 2] == 0 &), ___}]
    Expr* e6 = parse_or_die("MatchQ[Range[10], {___, _?(Mod[#, 2] == 0 &), ___}]");
    Expr* res6 = evaluate(e6);
    ASSERT(res6->type == EXPR_SYMBOL && strcmp(res6->data.symbol.name, "True") == 0);
    expr_free(e6); expr_free(res6);

    env_free(env);
}

/* A failing /; condition (or later argument) must backtrack into the other
 * assignments of a nested Orderless match, not commit to its first split
 * (v0.451).  The first assignment of a_ f_ against z x^k binds a -> x^k
 * whenever x^k sorts first; FreeQ[a, x] then fails and a -> z must be tried. */
void test_match_condition_backtracks_into_nested() {
    assert_eval_eq("MatchQ[g[z x^k], g[a_ f_] /; FreeQ[a, x]]", "True", 0);
    assert_eval_eq("MatchQ[g[z x^k], g[(a_ f_) /; FreeQ[a, x]]]", "True", 0);
    assert_eval_eq("ReplaceAll[g[z x^k], g[a_ f_] /; FreeQ[a, x] :> {a, f}]", "{z, x^k}", 0);
    assert_eval_eq("Module[{tt}, tt[a_ f_, v_] /; FreeQ[a, v] := {a, f}; tt[z x^k, x]]", "{z, x^k}", 0);
    /* a later ARGUMENT failing also re-splits the nested match */
    assert_eval_eq("MatchQ[h[z x^k, z], h[a_ f_, a_]]", "True", 0);
    assert_eval_eq("MatchQ[g[z x^k], g[a_ f_] /; FreeQ[{a, f}, z]]", "False", 0);
    /* the integrand that used to spin: the constant factor pulls out.  The
     * power-rule denominator 1 + 1/(-1 + a) is combined by Cancel, which no
     * longer declines on a generator whose EXPONENT is rational in a variable
     * (y^(1/(a - 1)) is one opaque generator); the reference's form is
     * ((-1 + a) x^(a/(-1 + a)))/a, the same antiderivative. */
    assert_eval_eq("Integrate[x^(1/(a-1)) y^(1/(a-1)), x]",
                   "(-x^(1 + 1/(-1 + a)) y^(1/(-1 + a)) + a x^(1 + 1/(-1 + a)) y^(1/(-1 + a)))/a", 0);
}

/* A Condition over a OneIdentity head with an Optional element: the subject
 * need not carry the head, because the Optional takes its default (h[g]
 * against Times[a_., h[b_]] as a -> 1), exactly as the unguarded pattern
 * matches. The guarded form used to fail the head test outright, so a
 * package's per-term structure check (a_. h[b_] /; ...) rejected every
 * coefficient-free term. */
void test_match_condition_oneidentity_optional_default() {
    assert_eval_eq("MatchQ[h[g], a_. h[b_] /; FreeQ[{a, b}, h]]", "True", 0);
    assert_eval_eq("Replace[h[g], a_. h[b_] /; b === g :> {a, b}]", "{1, g}", 0);
    assert_eval_eq("Replace[h[k], a_. h[b_] /; b === g :> {a, b}]", "h[k]", 0);
    assert_eval_eq("MatchQ[3 h[g], a_. h[b_] /; a == 3]", "True", 0);
    assert_eval_eq("Replace[x, a_. + y_ /; y === x :> {a, y}]", "{0, x}", 0);
    assert_eval_eq("Cases[{h[1], 2 h[3], h[4] h[5], k[1]}, a_. h[b_] /; OddQ[b] :> {a, b}]",
                   "{{1, 1}, {2, 3}, {h[4], 5}}", 0);
    /* a head that is not OneIdentity still fails on a mismatch */
    assert_eval_eq("MatchQ[g[x], f[x_] /; True]", "False", 0);
    /* and as a DownValue */
    assert_eval_eq("qq[a_. h[b_] /; b > 1] := {a, b}; {qq[h[2]], qq[5 h[3]], qq[h[0]]}",
                   "{{1, 2}, {5, 3}, qq[h[0]]}", 0);
}

/* A TimeConstrained timeout restores registered nesting counters (v0.451): an
 * aborted nested DSolve must not leave DSolve's depth raised, which silently
 * disabled the top-level-only methods for the rest of the session.  After the
 * abort, the top-level-only DSolve`InverseOrientation still answers. */
void test_timeconstrained_restores_dsolve_depth() {
    /* an equation the cascade grinds on: the 2 s budget aborts it mid-recursion */
    assert_eval_eq("TimeConstrained[DSolve[(a + x y[x]^2) Sin[x] y'[x] + b Cos[y[x]] y[x] == x^3 y[x]^7, y, x], 2, ok]", "ok", 0);
    assert_eval_eq("Head[DSolve[y'[x] == 1/(x y[x] + x^3 y[x]^3), y, x]]", "List", 0);
}

/* A call with a pattern head (g_[y__]) matches ONE operand of a Flat call,
 * never a group h[...] of several (the reference's semantics). */
static void test_match_flat_pattern_head_single_operand(void) {
    assert_eval_eq("ClearAll[ff, q, r]; SetAttributes[ff, {Flat, Orderless}]; "
                   "{ReplaceList[ff[a[1], b[2], c[3]], ff[g_[y__], c[3]] :> {g, {y}}], "
                   "ff[a[1], b[2], c[3]] /. ff[g_[y__], c[3]] :> {g, {y}}, "
                   "ReplaceList[ff[a[1], b[2], c[3]], ff[x_, c[3]] :> x]}",
                   "{{}, ff[b[2], {a, {1}}], {ff[a[1], b[2]]}}", 0);
    assert_eval_eq("SetAttributes[q, Orderless]; "
                   "q /: q[z_, x_] f_[y__] := r[z, f]; FreeQ[q[a, b] e[a] g[b], r[_, Times]]",
                   "True", 0);
}

/* A failing c1_. Cosh[a_] + c2_. Cosh[b_] + r___ /; cond over a wide sum is
 * polynomial: neither the Optional-coefficient element (whose required
 * Cosh[a_] cannot take a Plus group) nor the leftover retry (r___ already
 * catches the rest) enumerates subsets.  Both were 2^n -- 20 s for TrigReduce
 * on a 16-term Cosh/Sinh sum; at 30 terms either alone would never finish. */
static void test_match_flat_orderless_no_subset_blowup(void) {
    assert_eval_eq("ws = Sum[c[k] Cosh[k u], {k, 30}]; "
                   "TimeConstrained[{MatchQ[ws, c1_. Cosh[a_] + c2_. Cosh[b_] + r___ /; SameQ[a + b, 0]], "
                   "ws /. (c1_. Cosh[a_] + c2_. Cosh[b_] + r___ /; SameQ[a + b, 0]) :> 0}, 10, $Aborted] === {False, ws}",
                   "True", 0);
    assert_eval_eq("3 Cosh[x] + 2 Cosh[-x] + y /. c1_. Cosh[a_] + c2_. Cosh[b_] + r___ "
                   "/; SameQ[Expand[a + b], 0] :> (c1 + c2) Cosh[a] + r",
                   "y + 5 Cosh[x]", 0);
    /* A condition that reads r still gets the leftover retry. */
    assert_eval_eq("{a + b + c + d /. (a + r___ /; Length[{r}] == 2) :> f[r], "
                   "a + b + c + d /. (a + r__ /; Length[{r}] == 1) :> f[r], "
                   "2 x y /. c_. x^n_. :> {c, n}}",
                   "{d + f[b, c], c + d + f[b], {2 y, 1}}", 0);
}

/* Issue #116: every pattern element used to be matched with "the rest of this
 * level" as a continuation, so match depth grew per ARGUMENT and hit the
 * $RecursionLimit guard at 512 arguments -- MatchQ[e, e] was False for a
 * 512-term sum.  A leading run of elements with a single way to match
 * (literals; _, x_, x_h in an ordered non-Flat head) is now consumed
 * iteratively, at any width. */
static void test_match_wide_call_no_depth_limit(void) {
    assert_eval_eq("good = Sum[f[i], {i, 511}]; bad = good + f[512]; "
                   "{MatchQ[good, good], MatchQ[bad, bad], MatchQ[bad, bad + f[513]]}",
                   "{True, True, False}", 0);
    assert_eval_eq("w = Array[x, 2000]; {MatchQ[Times @@ w, Times @@ w], MatchQ[g @@ w, g @@ w], "
                   "MatchQ[g @@ Range[2000], g @@ Table[_Integer, 2000]], "
                   "MatchQ[g @@ Range[2000], g @@ Append[Table[_, 1999], _String]], "
                   "Replace[g @@ Range[2000], g[a_, b_, Sequence @@ Table[_, 1996], c_, d_] :> {a, b, c, d}]}",
                   "{True, True, True, False, {1, 2, 1999, 2000}}", 0);
    /* repeated names still constrain; leading literals keep enumeration order */
    assert_eval_eq("{MatchQ[g[1, 1, 2], g[x_, x_, _]], MatchQ[g[1, 2, 2], g[x_, x_, _]], "
                   "ReplaceList[o[a, b, c] /. o -> Plus, a + x_ + y_ :> {x, y}]}",
                   "{True, False, {{b, c}, {c, b}}}", 0);
}

/* Issue #115: FreeQ also finds a Flat call form as a SUBSET of a longer call,
 * as ReplaceAll does; MatchQ, MemberQ, Count and Cases stay whole-expression. */
static void test_freeq_flat_subset(void) {
    assert_eval_eq("e = a f[1] f[2]; {FreeQ[e, x_f y_f], FreeQ[e, x_f y_f, {0}], "
                   "FreeQ[e, x_f y_f, {0, Infinity}], FreeQ[e, x_f y_f, {1}], "
                   "MatchQ[e, x_f y_f], MemberQ[{e}, x_f y_f], Count[e, x_f y_f, {0}]}",
                   "{False, False, False, True, False, False, 0}", 0);
    assert_eval_eq("{FreeQ[a + b + c, a + b], FreeQ[a + b + c, a + d], FreeQ[d**a**b**c, a**b], "
                   "FreeQ[d**a**b**c, b**a], FreeQ[a && b && c, a && c]}",
                   "{False, True, False, True, True}", 0);
    assert_eval_eq("{FreeQ[f[g[a]], a, Heads -> False], FreeQ[f[g[a]], g, Heads -> False], "
                   "FreeQ[f[g[a]], g[a], {2, Infinity}], FreeQ[f[g[a]], a, {2, Infinity}]}",
                   "{False, True, True, False}", 0);
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
    TEST(test_match_condition_backtracks_into_nested);
    TEST(test_match_condition_oneidentity_optional_default);
    TEST(test_match_flat_pattern_head_single_operand);
    TEST(test_match_flat_orderless_no_subset_blowup);
    TEST(test_timeconstrained_restores_dsolve_depth);
    TEST(test_match_wide_call_no_depth_limit);
    TEST(test_freeq_flat_subset);
    printf("All match tests passed!\n");
    return 0;
}
