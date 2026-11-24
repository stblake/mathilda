/* Builtin-audit regressions: functional programming (review/issues/; tests/test_audit.h). */
#include "test_audit.h"

/* FixedPointList#1 (high): "FixedPointList segfaults on long symbolic
 * iterations (e.g. FixedPointList[f, x, 100000]); run time grows
 * quadratically".  SameQ now sits behind the memoised-hash pre-check, so each
 * step is O(1) for a growing iterate. */
static void test_fixedpointlist_long_symbolic(void) {
    audit_chk("FixedPointList#1", "Length[FixedPointList[f, x, 100000]]", "100001");
    audit_true("FixedPointList#1",
               "Last[FixedPointList[f, x, 100000]] === Nest[f, x, 100000]");
    audit_true("FixedPointList#1",
               "FixedPoint[g, x, 100000] === Nest[g, x, 100000]");
    /* Convergence is still detected by SameQ. */
    audit_chk("FixedPointList#1", "FixedPointList[Floor[#/2] &, 20]", "{20, 10, 5, 2, 1, 0, 0}");
    audit_chk("FixedPointList#1", "FixedPoint[Floor[#/2] &, 2^100]", "0");
    audit_chk("FixedPointList#1", "FixedPointList[# /. a -> b &, {a, a}]", "{{a, a}, {b, b}, {b, b}}");
    /* Freeing the chain: expr_free loops on the last argument, so releasing a
     * 2*10^6-deep chain needs constant C stack (recursion overflowed the 8 MB
     * default at ~10^5 levels once the list was dropped). */
    audit_chk("FixedPointList#1",
              "Module[{c = Nest[g, y, 2000000]}, c = 0; Length[FixedPointList[f, x, 300000]]]",
              "300001");
    /* 0. and -0. are SameQ, so they hash alike and still converge. */
    audit_chk("FixedPointList#1", "Length[FixedPointList[-# &, 0., 10]]", "2");
    audit_true("FixedPointList#1", "Hash[0.] === Hash[-0.]");
    audit_chk("FixedPointList#1", "FixedPointList[f, x, SameTest -> (True &)]", "{x, f[x]}");
}

/* TakeLargest#1 (high): "TakeLargest ranks symbolic numeric constants (Pi, E,
 * -Pi, Sqrt[2]) above every explicit number".  Keys rank in certified numeric
 * order (src/numeric_order.c), not canonical order. */
static void test_takelargest_numeric_constants(void) {
    audit_true("TakeLargest#1", "TakeLargest[{1, Pi, 4}, 1] === {4}");
    audit_true("TakeLargest#1", "TakeLargest[{-Pi, 1}, 1] === {1}");
    audit_true("TakeLargest#1", "TakeLargest[{E, 3, Pi}, 3] === {Pi, 3, E}");
    audit_true("TakeLargest#1", "TakeSmallest[{Pi, 3, Sqrt[2], 1/2}, 2] === {1/2, Sqrt[2]}");
    /* close values are certified, not read through doubles */
    audit_true("TakeLargest#1", "TakeLargest[{314159265358979/10^14, Pi}, 1] === {Pi}");
    audit_true("TakeLargest#1", "RankedMin[{Pi, 3, E}, 2] === 3");
    audit_true("TakeLargest#1", "RankedMax[{Pi, 314159265358979/10^14, E}, 1] === Pi");
    /* ties settle by position; machine data unchanged */
    audit_true("TakeLargest#1", "TakeLargest[{3, 3., 1}, 2] === {3, 3.}");
    audit_true("TakeLargest#1", "TakeLargest[{5, 1, 4, 2}, 2] === {5, 4}");
    audit_true("TakeLargest#1",
               "Module[{v = RandomReal[1, 500]}, TakeLargest[v, 5] === TakeLargest[Developer`FromPackedArray[v], 5]]");
}

/* TakeLargest#2 (high): "TakeLargest treats -Infinity as larger than finite
 * numbers". */
static void test_takelargest_minus_infinity(void) {
    audit_true("TakeLargest#2", "TakeLargest[{-Infinity, 3}, 1] === {3}");
    audit_true("TakeLargest#2", "TakeLargest[{-Infinity, Infinity, 3}, 2] === {Infinity, 3}");
    audit_true("TakeLargest#2", "TakeSmallest[{-Infinity, Infinity, 3}, 1] === {-Infinity}");
}

/* TakeLargestBy#1 (high): "Exact numeric keys (Pi, E, Sqrt[2]) are ranked by
 * canonical order, not numeric value". */
static void test_takelargestby_numeric_keys(void) {
    audit_true("TakeLargestBy#1", "TakeLargestBy[{1, 2, 3}, {Pi, 3, E}[[#]] &, 2] === {1, 2}");
    audit_true("TakeLargestBy#1", "TakeLargestBy[{1, 2, 3}, {Sqrt[2], 1.5, 1}[[#]] &, 1] === {2}");
    audit_true("TakeLargestBy#1", "TakeLargestBy[{Pi, 3, E}, N, 2] === {Pi, 3}");
}

/* TakeSmallestBy#1 (high): "TakeSmallestBy ranks exact numeric keys like
 * Sqrt[2] and E by canonical order, not numeric value". */
static void test_takesmallestby_numeric_keys(void) {
    audit_true("TakeSmallestBy#1", "TakeSmallestBy[{Sqrt[2], 1.5, 1}, Identity, 2] === {1, Sqrt[2]}");
    audit_true("TakeSmallestBy#1", "TakeSmallestBy[{E, 2.8, 3}, Identity, 1] === {E}");
    audit_true("TakeSmallestBy#1", "TakeSmallestBy[{Pi, 3, 4}, Identity, 1] === {3}");
}

/* Function#1 / Slot#1 (high): "#0 (the pure function itself) is never
 * substituted, so recursive pure functions fail". */
static void test_function_slot_zero(void) {
    audit_chk("Function#1", "If[# == 0, 1, # #0[# - 1]] &[5]", "120");
    audit_chk("Slot#1", "If[#1 <= 1, 1, #1 #0[#1 - 1]] &[6]", "720");
    audit_true("Slot#1", "(f[#0] &)[a] === f[f[#0] &]");
    audit_true("Slot#1", "(#0 &)[] === (#0 &)");
    /* sanity: ordinary slots */
    audit_chk("Slot#1", "(#1 + #2 &)[3, 4]", "7");
    audit_chk("Slot#1", "Map[#^2 &, {1, 2, 3}]", "{1, 4, 9}");
}

/* Function#2 (high): "Nested Function captures a substituted argument with
 * the same name as its parameter (no lexical renaming)".  The inner
 * parameter is alpha-renamed, so the free y stays free. */
static void test_function_nested_capture(void) {
    audit_run("ClearAll[y, x, z]");
    audit_true("Function#2", "Function[x, Function[y, x + y]][y][1] === 1 + y");
    audit_true("Function#2", "Function[x, Function[{y}, x + y]][y][3] === 3 + y");
    audit_true("Function#2", "Function[x, Function[y, x y]][2 y][5] === 10 y");
    audit_true("Function#2", "Function[x, Function[y, x + y]][y] =!= Function[y, y + y]");
    /* sanity: no clash, no rename; a rebinding still shadows */
    audit_true("Function#2", "Function[x, Function[y, x + y]][a] === Function[y, a + y]");
    audit_true("Function#2", "Function[x, Function[x, x]][1] === Function[x, x]");
    audit_true("Function#2", "Function[x, (x + #) &][a] === (a + #1 &)");
}

/* Function#3 (high): "Function substituting into Module/With captures the
 * argument in the inner local variable". */
static void test_function_scope_capture(void) {
    audit_run("ClearAll[y, x]");
    audit_true("Function#3", "Function[x, With[{y = 2}, x + y]][y] === 2 + y");
    audit_true("Function#3", "Module[{r = Function[x, Module[{y}, x + y]][y]}, "
                             "Length[r] == 2 && MemberQ[r, y] && Count[r, y] == 1]");
    audit_true("Function#3", "(With[{y = 3}, # + y] &)[y] === 3 + y");
    /* a Module / With that binds the parameter shadows it in the body only */
    audit_chk("Function#3", "Function[x, Module[{x = x + 1}, x]][5]", "6");
    audit_chk("Function#3", "Function[{x, z}, With[{x = z}, x + z]][1, 2]", "4");
    audit_chk("Function#3", "Function[x, Module[{y = x}, y + 1]][5]", "6");
}

/* Apply#1 / Map#1 (high): "level spec {m, Infinity} silently falls back"
 * (Apply to {0}, Map to level 1).  One shared parser now serves every
 * level-spec head. */
static void test_level_spec_infinity(void) {
    audit_run("ClearAll[f, a, b, c, d]");
    audit_chk("Apply#1", "Apply[f, {{a, {b}}, {c, d}}, {0, Infinity}]", "f[f[a, f[b]], f[c, d]]");
    audit_chk("Apply#1", "Apply[f, {{a, {b}}, {c, d}}, {1, Infinity}]", "{f[a, f[b]], f[c, d]}");
    audit_chk("Apply#1", "Apply[f, {{a}}, All]", "f[f[a]]");
    audit_chk("Map#1", "Map[f, {{a, b}, {c}}, {0, Infinity}]", "f[{f[{f[a], f[b]}], f[{f[c]}]}]");
    audit_chk("Map#1", "Map[f, {{a, b}, {c}}, {2, Infinity}]", "{{f[a], f[b]}, {f[c]}}");
    audit_chk("Map#1", "Map[f, {a, b}, {0, Infinity}]", "f[{f[a], f[b]}]");
    audit_chk("Map#1", "Map[f, {a, {b}}, -1]", "{f[a], f[{f[b]}]}");
    /* the other level-spec heads agree */
    audit_chk("Map#1", "Level[{a, {b, c}}, -1]", "{a, b, c, {b, c}}");
    audit_chk("Map#1", "Level[{a, {b, c}}, {-1}]", "{a, b, c}");
    audit_chk("Map#1", "Cases[{1, {2}}, _Integer, {-1}]", "{1, 2}");
    audit_chk("Map#1", "Cases[{1, {2}}, _List, All]", "{{2}, {1, {2}}}");
    audit_chk("Map#1", "Count[{1, {2, 3}}, _Integer, {1, Infinity}]", "3");
    audit_chk("Map#1", "Replace[{a, {a}}, a -> b, {0, Infinity}]", "{b, {b}}");
    audit_chk("Map#1", "MapIndexed[f, {a, {b}}, {1, Infinity}]", "{f[a, {1}], f[{f[b, {2, 1}]}, {2}]}");
    audit_chk("Map#1", "FreeQ[{a, {b}}, b, {2, Infinity}]", "False");
    /* a malformed spec says Head::level and stays unevaluated */
    audit_chk("Map#1", "Check[Map[f, {a}, x], fired, Map::level]", "fired");
    audit_chk("Map#1", "Check[Apply[f, {a}, x], fired, Apply::level]", "fired");
    audit_chk("Map#1", "Check[Level[{a}, x], fired, Level::level]", "fired");
    /* sanity: defaults and an option in the level slot */
    audit_chk("Map#1", "Map[f, {a, b}]", "{f[a], f[b]}");
    audit_chk("Map#1", "Apply[f, {a, b}]", "f[a, b]");
    audit_chk("Map#1", "Cases[{1, a, 2}, _Integer]", "{1, 2}");
    audit_chk("Map#1", "Map[f, g[a], Heads -> True]", "f[g][f[a]]");
}

/* S14 (Wave 3): the shared level-spec parser takes BigInt and -Infinity
 * bounds; the level forms of AllTrue/AnyTrue/NoneTrue/FirstCase; Apply below
 * a held argument; the ::level / ::intnm messages. */
static void test_s14_level_specs(void) {
    audit_run("ClearAll[f, g, a, b, c, x, q]");
    /* Map#2: a level past any depth is Infinity */
    audit_true("Map#2", "Map[f, {{a}}, 10^20] === {f[{f[a]}]}");
    audit_true("Map#2", "Map[f, {{a}}, {1, 2^70}] === {f[{f[a]}]}");
    audit_true("Map#2", "Map[f, {{a}}, -10^20] === {{a}}");
    /* Apply#3: the new f[...] below Hold is not evaluated */
    audit_chk("Apply#3", "Apply[f, Hold[g[1 + 1]], {1}]", "Hold[f[1 + 1]]");
    audit_chk("Apply#3", "Apply[f, Hold[{1 + 1}], {1}]", "Hold[f[1 + 1]]");
    audit_true("Apply#3", "Apply[Plus, {{1, 2}, {3, 4}}, {1}] === {3, 7}");
    audit_true("Apply#3", "(Plus @@@ {{1, 2}, {3, x}}) === {3, 3 + x}");
    /* AllTrue#1 / AnyTrue#1 / NoneTrue#1: the level form (levels 1..n for n) */
    audit_chk("AllTrue#1", "AllTrue[{{2}, {3}}, EvenQ, {2}]", "False");
    audit_chk("AllTrue#1", "AllTrue[{{2}, {4}}, EvenQ, {2}]", "True");
    audit_chk("AllTrue#1", "AllTrue[{2, 4}, EvenQ, 1]", "True");
    audit_chk("AllTrue#1", "AllTrue[{{2, 4}, {6}}, IntegerQ[#] || ListQ[#] &, 2]", "True");
    audit_chk("AnyTrue#1", "AnyTrue[{{1, 3}, {2}}, EvenQ, 2]", "True");
    audit_chk("AnyTrue#1", "AnyTrue[{{1, 3}, {2}}, EvenQ, {2}]", "True");
    audit_chk("AnyTrue#1", "AnyTrue[{1, {2}}, EvenQ, Infinity]", "True");
    audit_chk("AnyTrue#1", "AnyTrue[{1, {3}}, EvenQ, Infinity]", "False");
    audit_chk("AnyTrue#1", "AnyTrue[{1, 2}, EvenQ, 1]", "True");
    audit_chk("NoneTrue#1", "NoneTrue[{{1, 2}, {3, 5}}, EvenQ, 2]", "False");
    audit_chk("NoneTrue#1", "NoneTrue[{{1, 3}, {3, 5}}, EvenQ, 2]", "True");
    audit_chk("NoneTrue#1", "NoneTrue[{{1, 2}, {3, 5}}, EvenQ, {2}]", "False");
    audit_chk("AllTrue#1", "Check[AllTrue[{1}, EvenQ, x], fired, AllTrue::level]", "fired");
    /* FirstCase#2: the level spec (fourth argument) */
    audit_chk("FirstCase#2", "FirstCase[{{1, 2}, {3, a}}, _Symbol, None, Infinity]", "a");
    audit_chk("FirstCase#2", "FirstCase[{{1, 2}, {3, a}}, _Symbol, None]", "None");
    audit_chk("FirstCase#2", "FirstCase[{{1, 2}, {3, a}}, _Integer, None, {2}]", "1");
    audit_chk("FirstCase#2", "Check[FirstCase[{1}, _, None, x], fired, FirstCase::level]", "fired");
    /* MapThread#1: a bad level is ::intnm, unevaluated */
    audit_chk("MapThread#1", "Check[MapThread[f, {{1, 2}, {3, 4}}, -1], fired, MapThread::intnm]",
              "fired");
    audit_chk("MapThread#1", "Check[MapThread[f, {{1, 2}, {3, 4}}, 1.0], fired, MapThread::intnm]",
              "fired");
    audit_true("MapThread#1", "Quiet[Head[MapThread[f, {{1, 2}, {3, 4}}, -1]]] === MapThread");
    /* Scan#5: ::level */
    audit_chk("Scan#5", "Check[Scan[Print, {a}, x], fired, Scan::level]", "fired");
}

/* Wave 3 S12/S13: declared argument counts and option checks
 * (src/arity.h, src/arity_table.c).  The evaluator checks a call
 * against its head's declared arity before the C code runs and
 * prints argx/argr/argrx/argt(u)/argb(u)/argm(u), nonopt or optx;
 * the call stays unevaluated. */
static void test_wave3_arity_options(void) {
    audit_chk("Map#3", "Quiet[Check[Map[f, {a, b}, x], fired, Map::level]]", "fired");
    audit_chk("Map#3", "Quiet[Check[Map[f, {a, b}, {1, 1, 1}], fired, Map::level]]", "fired");
    audit_chk("Map#3", "Quiet[Check[Map[f, {a, b}, 1, 2], fired, Map::nonopt]]", "fired");
    audit_chk("Map#3", "Quiet[Check[Map[], fired, Map::argb]]", "fired");
    audit_chk("Map#3", "Map[f, a -> b]", "f[a] -> f[b]");
    audit_chk("NestList#2", "Quiet[Check[NestList[f, x], fired, NestList::argrx]]", "fired");
    audit_chk("NestList#2", "Quiet[NestList[f, x]]", "NestList[f, x]");
    audit_chk("Thread#2", "Quiet[Check[Thread[f[{a}], List, 1, 2], fired, Thread::argb]]", "fired");
    audit_chk("Thread#2", "Quiet[Thread[f[{a}], List, 1, 2]]", "Thread[f[{a}], List, 1, 2]");
    audit_chk("Distribute#1", "Quiet[Check[Distribute[], fired, Distribute::argb]]", "fired");
    audit_chk("Distribute#1", "Quiet[Distribute[]]", "Distribute[]");
    audit_chk("InterpolatingPolynomial#3", "Quiet[Check[InterpolatingPolynomial[{1, 4, 9}, x, Foo -> 1], fired, InterpolatingPolynomial::optx]]", "fired");
    audit_chk("InterpolatingPolynomial#3", "Quiet[InterpolatingPolynomial[{1, 4, 9}, x, Foo -> 1]]", "InterpolatingPolynomial[{1, 4, 9}, x, Foo -> 1]");
    audit_chk("LengthWhile#1", "Quiet[Check[LengthWhile[{1, 2}], fired, LengthWhile::argr]]", "fired");
    audit_chk("LengthWhile#1", "Quiet[LengthWhile[{1, 2}]]", "LengthWhile[{1, 2}]");
    audit_chk("NestWhile#1", "Quiet[Check[NestWhile[f, x], fired, NestWhile::argb]]", "fired");
    audit_chk("NestWhile#1", "Quiet[NestWhile[f, x]]", "NestWhile[f, x]");
    audit_chk("NestWhileList#2", "Quiet[Check[NestWhileList[f, x], fired, NestWhileList::argb]]", "fired");
    audit_chk("NestWhileList#2", "Quiet[NestWhileList[f, x]]", "NestWhileList[f, x]");
    audit_chk("TakeWhile#1", "Quiet[Check[TakeWhile[{1, 2}, EvenQ, 3], fired, TakeWhile::argrx]]", "fired");
    audit_chk("TakeWhile#1", "Quiet[TakeWhile[{1, 2}, EvenQ, 3]]", "TakeWhile[{1, 2}, EvenQ, 3]");
    audit_chk("MapAll#1", "Quiet[Check[MapAll[f, {a, b}, 3], fired, MapAll::nonopt]]", "fired");
    audit_chk("MapAll#1", "Quiet[MapAll[f, {a, b}, 3]]", "MapAll[f, {a, b}, 3]");
    audit_chk("MapAll#2", "Quiet[Check[MapAll[f, {a, b}, Foo -> 1], fired, MapAll::optx]]", "fired");
    audit_chk("MapAll#2", "Quiet[MapAll[f, {a, b}, Foo -> 1]]", "MapAll[f, {a, b}, Foo -> 1]");
    audit_chk("Apply#5", "Quiet[Check[Apply[f, {a, b}, {1}, 7], fired, Apply::nonopt]]", "fired");
    audit_chk("Apply#5", "Quiet[Apply[f, {a, b}, {1}, 7]]", "Apply[f, {a, b}, {1}, 7]");
    audit_chk("Scan#3", "Quiet[Check[Scan[s12p, {a}, 1, 2], fired, Scan::nonopt]]", "fired");
    audit_chk("Scan#3", "Quiet[Scan[s12p, {a}, 1, 2]]", "Scan[s12p, {a}, 1, 2]");
    audit_chk("Scan#4", "Quiet[Check[Scan[s12p, {a}, Foo -> True], fired, Scan::optx]]", "fired");
    audit_chk("Scan#4", "Quiet[Scan[s12p, {a}, Foo -> True]]", "Scan[s12p, {a}, Foo -> True]");
    audit_chk("MapAll#2", "MapAll[f, g[a], Heads -> True]", "f[f[g][f[a]]]");
    audit_chk("Scan#4", "Reap[Scan[Sow, f[a], Heads -> True]][[2]]", "{{f, a}}");
}

int main(void) {
    audit_init();
    TEST(test_fixedpointlist_long_symbolic);
    TEST(test_takelargest_numeric_constants);
    TEST(test_takelargest_minus_infinity);
    TEST(test_takelargestby_numeric_keys);
    TEST(test_takesmallestby_numeric_keys);
    TEST(test_function_slot_zero);
    TEST(test_function_nested_capture);
    TEST(test_function_scope_capture);
    TEST(test_level_spec_infinity);
    TEST(test_s14_level_specs);
    TEST(test_wave3_arity_options);
    return audit_finish("audit_functional_programming");
}
