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

int main() {
    symtab_init();
    core_init();
    
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
    run_test("modr[1, 5] + modr[2, 5] + xr /. modr[a_, p_] + modr[b_, p_] :> modr[Mod[a + b, p], p]", "modr[3, 5] + xr");

    printf("All Replace tests passed!\n");
    symtab_clear();
    return 0;
}
