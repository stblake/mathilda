#include <stdio.h>
#include <assert.h>
#include <string.h>
#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "replace.h"
#include "core.h"
#include "print.h"

void run_test(const char* input, const char* expected) {
    Expr* expr = parse_expression(input);
    assert(expr != NULL);
    Expr* evaled = evaluate(expr);
    char* s = expr_to_string(evaled);
    if (strcmp(s, expected) != 0) {
        printf("FAIL: %s\nExpected: %s\nGot:      %s\n", input, expected, s);
        assert(false);
    } else {
        printf("PASS: %s -> %s\n", input, s);
    }
    free(s);
    expr_free(expr);
    expr_free(evaled);
}

/* A Rational or Complex number is an atom: replacement, Cases, Count,
 * Position and DeleteCases never reach its parts or its head, while FreeQ and
 * MemberQ (with Heads -> True) see its head at its own level. Every expected
 * value is the reference's. {1/2} /. 2 -> 3 was {1/3}; FeynCalc's
 * MomentumCombine gates on FreeQ[expr, Complex]. */
static void test_number_atoms(void) {
    run_test("{{1/2} /. 2 -> 3, {1/2} /. 1 -> 5, {1/2} /. Rational -> f, {I} /. Complex -> h, {1/2} //. 2 -> 3}",
             "{{1/2}, {1/2}, {1/2}, {I}, {1/2}}");
    run_test("{1/2} /. Rational[a_, b_] -> g[a, b]", "{g[1, 2]}");
    run_test("Replace[{1/2}, Rational -> f, Infinity, Heads -> True]", "{1/2}");
    run_test("{Count[{1/2}, 2, Infinity], Cases[{1/2, I}, 2 | 0 | 1, Infinity], Position[{1/2}, 2], Position[{I}, Complex]}",
             "{0, {}, {}, {}}");
    run_test("{DeleteCases[{1/2, x}, 2, Infinity], DeleteCases[{I}, 0, Infinity]}", "{{1/2, x}, {I}}");
    run_test("{FreeQ[I x, Complex], FreeQ[f[1/2], Rational], FreeQ[1/2, Rational, Heads -> False], FreeQ[f[I], Complex, 1], "
             "FreeQ[f[I], Complex, {2}], FreeQ[f[1/2], 2]}",
             "{False, False, True, False, True, True}");
    run_test("{MemberQ[{1/2}, 2, Infinity], MemberQ[f[1/2], Rational, {1}, Heads -> True], MemberQ[{1/2}, Rational, Infinity]}",
             "{False, True, False}");
}

int main() {
    symtab_init();
    core_init();
    test_number_atoms();
    
    // ReplacePart[expr, i -> new]
    run_test("ReplacePart[{a, b, c, d}, 2 -> x]", "{a, x, c, d}");
    
    // Negative indices
    run_test("ReplacePart[{5, 4, 3, 2, 1}, _?PrimeQ -> \"p\"]", "{5, \"p\", \"p\", 2, \"p\"}");
    
    // ReplacePart[expr, {i, j} -> new]
    run_test("ReplacePart[{{a, b}, {c, d}}, {2, 1} -> x]", "{{a, b}, {x, d}}");
    
    // ReplacePart[expr, {{i1, j1}, {i2, j2}} -> new]
    run_test("ReplacePart[{{a, b}, {c, d}}, {{1, 1}, {2, 2}} -> x]", "{{x, b}, {c, x}}");
    
    // ReplacePart[expr, {i1 -> new1, i2 -> new2}]
    run_test("ReplacePart[{a, b, c, d}, {1 -> x, 3 -> y}]", "{x, b, y, d}");
    
    // Check replacing head
    run_test("ReplacePart[f[x], 0 -> g]", "g[x]");
    
    // Replace
    run_test("Replace[x^2, x^2 -> a + b]", "a + b");
    run_test("Replace[1 + x^2, x^2 -> a + b]", "1 + x^2");
    run_test("Replace[x, {{x -> a}, {x -> b}}]", "{a, b}");
    run_test("Replace[f[2], f[x_] :> x]", "2");
    run_test("Replace[f[2], {f[x_] :> x, f[2] :> \"two\"}]", "2");
    run_test("Replace[f[2], {{f[x_] :> x}, {f[2] :> \"two\"}}]", "{2, \"two\"}");
    run_test("Replace[x, {y -> 2, z -> 3}]", "x");
    run_test("Replace[1 + x^2, x^2 -> a + b, {1}]", "1 + a + b");
    run_test("Replace[{1, 3, 2, x, 6, Pi}, _?PrimeQ -> \"prime\", {1}]", "{1, \"prime\", \"prime\", x, 6, Pi}");
    run_test("Replace[{1, 3, 2, x, 6, Pi}, t_ /; Mod[t, 3] == 0 -> \"3k\", {1}]", "{1, \"3k\", 2, x, \"3k\", Pi}");
    run_test("Replace[{1, Pi, {2 z, Cos[x]}}, x_ -> f[x], {-1}]", "{f[1], f[Pi], {f[2] f[z], Cos[f[x]]}}");
    run_test("Replace[f[f[f[f[x]]]], f[x_] :> g[x], 3]", "f[g[g[g[x]]]]");
    run_test("Replace[f[f[f[f[x]]]], f[x_] :> g[x], {0, 2}]", "g[g[g[f[x]]]]");
    run_test("Replace[f[f[f[f[x]]]], f[x_] :> g[x], All]", "g[g[g[g[x]]]]");
    run_test("Replace[f[1, 2, f[3]], f -> g, All]", "f[1, 2, f[3]]");
    run_test("Replace[f[1, 2, f[3]], f -> g, All, Heads -> True]", "g[1, 2, g[3]]");
    run_test("Replace[f[1, 2, f[3]], f -> g, {1}, Heads -> True]", "g[1, 2, f[3]]");
    run_test("Replace[f[1, 2, f[3]], f -> g, {2}, Heads -> True]", "f[1, 2, g[3]]");
    run_test("Replace[f[1, 2, f[3]], Heads -> True]", "f[1, 2, f[3]]");
    
    // ReplaceAll
    run_test("{x, x^2, y, z} /. x -> 1", "{1, 1, y, z}");
    run_test("{x, x^2, y, z} /. x -> {a, b}", "{{a, b}, {a^2, b^2}, y, z}");
    run_test("Sin[x] /. Sin -> Cos", "Cos[x]");
    run_test("{a, b, c} /. List -> f", "f[a, b, c]");
    run_test("1 + x^2 + x^4 /. x^p_ :> f[p]", "1 + f[2] + f[4]");
    run_test("x /. {x -> 1, x -> 3, x -> 7}", "1");
    run_test("x /. {{x -> 1}, {x -> 3}, {x -> 7}}", "{1, 3, 7}");
    run_test("{f[2], f[x, y], h[], f[]} /. f[x_] -> \"OK\"", "{\"OK\", f[x, y], h[], f[]}");
    run_test("{f[2], f[x, y], h[], f[]} /. f[x__] -> \"OK\"", "{\"OK\", \"OK\", h[], f[]}");
    run_test("{f[2], f[x, y], h[], f[]} /. f[x___] -> \"OK\"", "{\"OK\", \"OK\", h[], \"OK\"}");
    run_test("{1, 3, 2, x, 6, Pi} /. _Integer :> \"int\"", "{\"int\", \"int\", \"int\", x, \"int\", Pi}");
    run_test("{1, 3, 2, x, 6, Pi} /. _?PrimeQ -> \"prime\"", "{1, \"prime\", \"prime\", x, 6, Pi}");
    run_test("{1, 3, 2, x, 6, Pi} /. t_ /; Mod[t, 3] == 0 -> \"3k\"", "{1, \"3k\", 2, x, \"3k\", Pi}");
    run_test("{{a, 1}, {b, 2}} /. {x_, y_} :> Subscript[x, y]", "Subscript[{a, 1}, {b, 2}]");
    run_test("{{a, 1}, {b, 2}} /. {x_Symbol, y_} :> Subscript[x, y]", "{Subscript[a, 1], Subscript[b, 2]}");
    run_test("Hold[x + x] /. x -> 7", "Hold[7 + 7]");
    run_test("Hold[x + x] /. x :> 2^2", "Hold[2^2 + 2^2]");
    run_test("Hold[x + x] /. x -> 2^2", "Hold[4 + 4]"); // Note: Power evaluated outside Hold.
    run_test("x /. {}", "x");
    run_test("x /. {{x -> 1}, {y -> 2}}", "{1, x}");
    run_test("{a, b, c} /. {a -> b, b -> d}", "{b, d, c}");
    run_test("{a, b, c} /. a -> b /. b -> d", "{d, d, c}");

    // ReplaceRepeated
    run_test("{a, b, c} //. {a -> b, b -> d}", "{d, d, c}");
    run_test("{1, 3, 4, 2, 3, 4, 5, 4, 2, 4, 6, 2, 1} //. {a___, x_, b___, x_, c___} -> {a, x, b, c}", "{1, 3, 4, 2, 5, 6}"); // Delete duplicates test case

    // Repeated
    run_test("{{}, {a,a}, {a,b}, {a,a,a}, {a}} /. {a..} -> x", "{{}, x, {a, b}, x, x}");
    run_test("{{}, {f[a],f[b]}, {f[a]}, {f[a,b]}, {f[a],g[b]}} /. {f[_]..} -> x", "{{}, x, x, {f[a, b]}, {f[a], g[b]}}");
    run_test("{f[], f[a,a], f[a,b], f[a,a,a]} /. f[a..] -> x", "{f[], x, f[a, b], x}");
    run_test("{{}, {a}, {a,a}, {a,a,a}, {a,a,a,a}} /. {Repeated[a, 3]} -> x", "{{}, x, x, x, {a, a, a, a}}");
    run_test("{{}, {a}, {a,a}, {a,a,a}, {a,a,a,a}} /. {Repeated[a, {2, 3}]} -> x", "{{}, {a}, x, x, {a, a, a, a}}");
    run_test("{{}, {a}, {a,b}, {a,b,c}, {a,b,c,d}} /. {Repeated[_, {3}]} -> x", "{{}, {a}, {a, b}, x, {a, b, c, d}}");
    run_test("{{}, {a}, {a,a}, {a,a,a}, {a,b}, {b,b}} /. {Repeated[_, {0, 2}]} -> x", "{x, x, x, {a, a, a}, x, x}");
    run_test("{{}, {a}, {a,a}, {a,a,a}, {a,b}, {b,b}} /. {Repeated[z_, {0, 3}]} -> x", "{x, x, x, x, {a, b}, x}");
    
    // RepeatedNull
    run_test("{{}, {a,a}, {a,b}, {a,a,a}, {a}} /. {a...} -> x", "{x, x, {a, b}, x, x}");
    run_test("{{}, {f[a],f[b]}, {f[a]}, {f[a,b]}, {f[a],g[b]}} /. {f[_]...} -> x", "{x, x, x, {f[a, b]}, {f[a], g[b]}}");
    run_test("{f[], f[a,a], f[a,b], f[a,a,a]} /. f[a...] -> x", "{x, x, f[a, b], x}");
    run_test("{{}, {a}, {a,a}, {a,a,a}, {a,a,a,a}} /. {RepeatedNull[a, 3]} -> x", "{x, x, x, x, {a, a, a, a}}");
    run_test("{{}, {a}, {a,a}, {a,a,a}, {a,a,a,a}} /. {RepeatedNull[a, {2, 3}]} -> x", "{{}, {a}, x, x, {a, a, a, a}}");
    run_test("{{}, {a}, {a,a}, {a,a,a}, {a,b}, {b,b}} /. {RepeatedNull[_, 2]} -> x", "{x, x, x, {a, a, a}, x, x}");
    run_test("{{}, {a}, {a,a}, {a,a,a}, {a,b}, {b,b}} /. {RepeatedNull[z_, 3]} -> x", "{x, x, x, x, {a, b}, x}");

    // Longest and Shortest
    run_test("{a,b,c,d,e,f} /. {x__, y__} :> {y}", "{b, c, d, e, f}"); // Default Shortest behavior
    run_test("{a,b,c,d,e,f} /. {Longest[x__], y__} :> {y}", "{f}");
    run_test("{a,b,c,d,e,f} /. {Shortest[x__], y__} :> {x}", "{a}");
    run_test("{a,b,c,d,e,f} /. {x__, Longest[y__]} :> {x}", "{a}");


    // ReplaceList
    run_test("ReplaceList[{a,b,c,d,e,f},{x__,y__}->{{x},{y}}]", "{{{a}, {b, c, d, e, f}}, {{a, b}, {c, d, e, f}}, {{a, b, c}, {d, e, f}}, {{a, b, c, d}, {e, f}}, {{a, b, c, d, e}, {f}}}");
    run_test("ReplaceList[x,{x->a,x->b,x->c}]", "{a, b, c}");
    run_test("ReplaceList[a+b+c,x_+y_->g[x,y]]", "{g[a, b + c], g[b, a + c], g[c, a + b], g[a + b, c], g[a + c, b], g[b + c, a]}");
    run_test("ReplaceList[{a,b,c,d},{___,x__,___}->{x}]", "{{a}, {a, b}, {a, b, c}, {a, b, c, d}, {b}, {b, c}, {b, c, d}, {c}, {c, d}, {d}}");
    run_test("ReplaceList[{a,b,b,b,c,c,a},{___,x_,x_,___}->x]", "{b, b, c}");
    run_test("ReplaceList[{a,b,c,a,d,b,d},{___,x_,y__,x_,___}->{x,{y}}]", "{{a, {b, c}}, {b, {c, a, d}}, {d, {b}}}");
    run_test("ReplaceList[a+b+c,x_+y_:>{x,y}]", "{{a, b + c}, {b, a + c}, {c, a + b}, {a + b, c}, {a + c, b}, {b + c, a}}");

    /* Flat-head leftover matching: a Flat (Plus/Times/...) rule rewrites a
     * SUBSET of a longer operand sequence and keeps the rest (PR #87 review). */
    run_test("1 + x + y /. x + y -> z", "1 + z");
    run_test("a + b + c /. x_ + y_ -> f[x, y]", "f[a, b + c]");   /* group absorption, unchanged */
    run_test("modr[1, 5] + modr[2, 5] /. modr[a_, p_] + modr[b_, p_] :> modr[Mod[a + b, p], p]", "modr[3, 5]");
    run_test("modr[1, 5] + modr[2, 5] + xr /. modr[a_, p_] + modr[b_, p_] :> modr[Mod[a + b, p], p]", "xr + modr[3, 5]");

    /* Ordered (non-Orderless) Flat heads: the matched run may sit anywhere --
     * $FlatPre$___ + $FlatRest$___ -- leftmost run first, one rewrite per
     * application (ReplaceAll does not re-traverse the result). */
    run_test("x ** y ** z /. y ** z -> w", "x ** w");
    run_test("x ** y ** z /. x ** y -> w", "w ** z");
    run_test("a . b . c . d /. b . c -> w", "a . w . d");
    run_test("SetAttributes[ffr, Flat]; ffr[x, y, z] /. ffr[y, z] -> w", "ffr[x, w]");
    run_test("ffr[x, y, z, q] /. ffr[y, z] -> w", "ffr[x, w, q]");
    run_test("ffr[a, b, a, b] /. ffr[a, b] -> w", "ffr[w, a, b]");
    run_test("ffr[x, a, b, a, b] /. ffr[a, b] -> w", "ffr[x, w, a, b]");
    run_test("ffr[a, b, a, b] //. ffr[a, b] -> w", "ffr[w, w]");
    run_test("ffr[x, b, a] /. ffr[a, b] -> w", "ffr[x, b, a]");          /* order kept */
    run_test("ffr[x, y, z] /. ffr[a_, b_] -> g[a, b]", "g[x, ffr[y, z]]");  /* full match first */
    run_test("ffr[x, gr[1], gr[2], y] /. ffr[gr[i_], gr[j_]] :> gr[i + j]", "ffr[x, gr[3], y]");
    run_test("ffr[x, gr[1], gr[2], y] /. HoldPattern[ffr[gr[i_], gr[j_]] /; i < j] :> gr[i + j]", "ffr[x, gr[3], y]");
    run_test("ffr[x, gr[2], gr[1], y] /. HoldPattern[ffr[gr[i_], gr[j_]] /; i < j] :> gr[i + j]", "ffr[x, gr[2], gr[1], y]");
    run_test("MatchQ[ffr[x, y, z], ffr[y, z]]", "False");                  /* MatchQ stays whole */
    /* a leading blank (no single-operand pre-screen) still finds the run */
    run_test("ffr[1, -2, 3, 4] /. ffr[a_?Negative, b_] :> gr[a, b]", "ffr[1, gr[-2, 3], 4]");
    run_test("ffr[1, 2, 2, 3] /. ffr[x_, x_] :> hr[x]", "ffr[1, hr[2], 3]");
    run_test("(p && q && r && s) /. (r && s) -> t", "p && q && t");     /* And: Flat, ordered */

    /* Offset / span screens (anchor = first single-operand element; bound
     * names fix the group span).  These must leave every result unchanged. */
    run_test("ffr[q, a, b, a, b] /. ffr[x_, x_] :> sq[x]", "ffr[q, sq[ffr[a, b]]]");
    run_test("(p && q && !q && r) /. And[x_, Not[x_]] :> False", "False");
    run_test("(p && q && r && !s) /. And[x_, Not[x_]] :> False", "p && q && r && !s");
    run_test("ffr[a, b, c, b, c] /. ffr[x__, x__] :> sq[x]", "ffr[a, sq[b, c]]");
    run_test("ffr[k, a, b, c, op[3], m] /. ffr[x_, y_, op[j_]] :> h[x, y, j]",
             "ffr[h[k, ffr[a, b, c], 3], m]");
    run_test("ffr[k, a, b, c, d] /. ffr[b, x_, d] :> h[x]", "ffr[k, a, h[c]]");
    /* an anchor with an Optional default may match a non-call operand */
    run_test("ffr[u, x, v] /. ffr[y_, x^n_.] :> g[y, n]", "ffr[g[u, 1], v]");
    /* a Condition inside the anchor is not used to screen */
    run_test("ffr[z, a, op[1], b, op[2], c] /. ffr[x_, op[j_ /; j > 1]] :> h[x, j]",
             "ffr[h[ffr[z, a, op[1], b], 2], c]");
    run_test("MatchQ[ffr[q, a, a], ffr[x_, x_]]", "False");
    /* PatternSequence is spliced, not a single operand: it must not pre-screen */
    run_test("ffr[x, a, b, y] /. ffr[PatternSequence[a, b]] -> w", "ffr[x, w, y]");
    /* scale: a failing And[x_, Not[x_]] over a long chain is O(n), not O(n^3) */
    run_test("Length[(And @@ Array[av, 1500]) /. And[x_, Not[x_]] :> False]", "1500");

    /* ReplaceRepeated MaxIterations: the default, a limit hit (rrlim, the
     * expression reached is returned), a limit not hit, Infinity, a bad value
     * (iopnf, unevaluated), an unknown option (optx), and SetOptions. */
    run_test("Options[ReplaceRepeated]", "{MaxIterations -> 65536}");
    run_test("ReplaceRepeated[f[a], {f[a] -> f[b]}, MaxIterations -> 20]", "f[b]");
    run_test("Quiet[ReplaceRepeated[0, x_Integer :> x + 1, MaxIterations -> 5]]", "5");
    run_test("Check[ReplaceRepeated[0, x_Integer :> x + 1, MaxIterations -> 5], \"msg\"]", "\"msg\"");
    run_test("Quiet[ReplaceRepeated[0, x_Integer :> x + 1, MaxIterations -> 0]]", "0");
    /* g[0] is reached on the 3rd scan, but confirming it needs a 4th. */
    run_test("Quiet[ReplaceRepeated[g[3], g[n_ /; n > 0] :> g[n - 1], MaxIterations -> 3]]", "g[0]");
    run_test("ReplaceRepeated[g[3], g[n_ /; n > 0] :> g[n - 1], MaxIterations -> 4]", "g[0]");
    run_test("ReplaceRepeated[f[1], f[n_ /; n < 10] :> f[n + 1], MaxIterations -> Infinity]", "f[10]");
    run_test("ReplaceRepeated[{a}, {{a -> b}, {a -> c}}, MaxIterations -> 2]", "{{b}, {c}}");
    run_test("Quiet[ReplaceRepeated[x, x -> y, MaxIterations -> -1]]",
             "ReplaceRepeated[x, x -> y, MaxIterations -> -1]");
    run_test("Check[ReplaceRepeated[x, x -> y, MaxIterations -> \"a\"], \"msg\"]", "\"msg\"");
    run_test("Check[ReplaceRepeated[x, {x -> y}, Method -> 1], \"msg\"]", "\"msg\"");
    run_test("SetOptions[ReplaceRepeated, MaxIterations -> 4]; r = Quiet[0 //. x_Integer :> x + 1]; "
             "SetOptions[ReplaceRepeated, MaxIterations -> 65536]; r", "4");

    /* Dispatch: normalisation, Normal, and every Replace* form accepting it. */
    run_test("Dispatch[a -> 1]", "Dispatch[{a -> 1}]");
    run_test("Dispatch[Dispatch[{a -> 1}]]", "Dispatch[{a -> 1}]");
    run_test("Dispatch[3]", "Dispatch[3]");
    run_test("Normal[Dispatch[{a -> 1, b -> 2}]]", "{a -> 1, b -> 2}");
    run_test("dsp = Dispatch[{a -> 1, b -> 2, f[x_] :> x^2}]; {a, b, f[3], c} /. dsp", "{1, 2, 9, c}");
    run_test("ReplaceAll[g[a], dsp]", "g[1]");
    run_test("{Replace[a, dsp], Replace[{a}, dsp, {1}], ReplaceList[a, dsp]}", "{1, {1}, {1}}");
    run_test("f[f[2]] //. dsp", "16");
    run_test("{a} /. Dispatch[{{a -> 1}, {a -> 2}}]", "{{1}, {2}}");
    /* The index never changes which rule wins: keyed rules, wild rules, guards,
     * Optional arguments (x^n_. matches a bare x), Verbatim, heads. */
    run_test("{x, x^2, x^3} /. Dispatch[{x^n_. :> p[n]}]", "{p[1], p[2], p[3]}");
    run_test("{f[1], f[2], h[0], h[5]} /. Dispatch[{f[1] -> one, z_Integer :> int[z], h[q_] /; q > 1 -> big}]",
             "{one, f[int[2]], h[int[0]], big}");
    run_test("{1, 2, a} /. Dispatch[{_Integer -> i, 1 -> one}]", "{i, i, a}");
    run_test("{x_, y} /. Dispatch[{Verbatim[x_] -> v}]", "{v, y}");
    run_test("a + b + z /. Dispatch[{a + b -> s}]", "s + z");
    /* Randomised: Dispatch and the plain rule list agree. */
    run_test("SeedRandom[3]; pool = {a, b, 1, f[x_], f[1], g[a, y_], h[q_] /; q > 1, x_Integer, _g, a | b, "
             "x^n_., p_?EvenQ, a + z_, 1/2, u_[1], h[5], {1, 2}, x^2, f[f[a]]}; "
             "ex = {a, b, 1, 2, f[1], f[2], f[a], g[a, 3], h[0], h[5], x, x^2, a + b + z, 1/2, k[1], {f[f[a]]}, {1, 2}, h[5]}; "
             "And @@ Table[rs = Table[RandomChoice[pool] -> RandomInteger[{100, 999}], {RandomInteger[{1, 10}]}]; "
             "(ex /. rs) === (ex /. Dispatch[rs]), {300}]", "True");
    run_test("Attributes[Dispatch]", "{Protected}");

    /* Pattern-free compound left-hand sides are keyed by the whole expression,
     * not only by their head, and first-match order still holds between a
     * literal rule and a pattern rule that both match. */
    run_test("{f[1], f[2], g[f[1]]} /. Dispatch[{f[x_] :> p[x], f[1] -> lit}]", "{p[1], p[2], g[p[1]]}");
    run_test("{f[1], f[2], g[f[1]]} /. Dispatch[{f[1] -> lit, f[x_] :> p[x]}]", "{lit, p[2], g[lit]}");
    run_test("{f[1], f[2]} /. Dispatch[{_f :> w, f[1] -> lit}]", "{w, w}");
    run_test("{f[1], f[2]} /. Dispatch[{f[1] /; False -> lit, f[1] -> lit2}]", "{lit2, f[2]}");
    run_test("{a + b, a + b + c, h[a + b + c]} /. Dispatch[{a + b -> s}]", "{s, c + s, h[c + s]}");
    run_test("{{1, 2}, {1, 2, 3}, Range[3]} /. Dispatch[{{1, 2, 3} -> lit, {1, 2} -> two}]", "{two, lit, lit}");
    run_test("{Range[100], {1}} /. Dispatch[{Range[100] -> big}]", "{big, {1}}");
    run_test("{Range[5.]} /. Dispatch[{{x_, y__} :> first[x]}]", "{first[1.0]}");
    run_test("{x^2, x^3, Sin[x^2]} /. Dispatch[{x^2 -> sq, Sin[x^2] -> s2}]", "{sq, x^3, s2}");
    run_test("{h[\"a\"], h[1.5], h[1], f[1][2], f[1]} /. Dispatch[{h[\"a\"] -> ha, h[1.5] -> hr, h[1] -> hi, f[1][2] -> c12}]",
             "{ha, hr, hi, c12, f[1]}");
    run_test("f[1] //. Dispatch[{f[1] -> f[2], f[2] -> done}]", "done");
    /* ... and that is what makes a table of them fast: 4000 f[i] rules against
     * 2000 subexpressions. Loose bound (a linear scan is ~100x slower), so the
     * test is not timing-flaky. */
    run_test("rls = Table[f[i] -> i^2, {i, 4000}]; dsp4 = Dispatch[rls]; SeedRandom[1]; "
             "ex4 = Table[f[RandomInteger[{1, 4000}]], {2000}]; "
             "tl = First[AbsoluteTiming[r1 = ex4 /. rls]]; td = First[AbsoluteTiming[r2 = ex4 /. dsp4]]; "
             "{r1 === r2, tl > 10 td}", "{True, True}");

    printf("All Replace tests passed!\n");
    symtab_clear();
    return 0;
}
