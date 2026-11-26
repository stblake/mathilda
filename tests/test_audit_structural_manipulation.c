/* Builtin-audit regressions: structural manipulation (review/issues/; tests/test_audit.h). */
#include "test_audit.h"

/* ArrayFlatten#1 (critical): "ArrayFlatten segfaults when a block is nested
 * more than 64 levels deep (stack buffer overflow in af_tensor_dims)".  A
 * block needs only its first r dimensions; deeper levels are content. */
static void test_array_flatten_deep_block(void) {
    audit_true("ArrayFlatten#1",
               "ArrayFlatten[{{Nest[List, 1, 100]}}] === Nest[List, 1, 100]");
    audit_true("ArrayFlatten#1",
               "ArrayFlatten[{{Nest[List, 1, 5000]}}] === Nest[List, 1, 5000]");
    audit_true("ArrayFlatten#1",
               "ArrayFlatten[{{{{Nest[List, 1, 200]}}}}, 1] === {{{Nest[List, 1, 200]}}}");
    /* Content jagged below the r flattened levels is carried along. */
    audit_chk("ArrayFlatten#1", "ArrayFlatten[{{{{1, {2}}, {3, 4}}}}]", "{{1, {2}}, {3, 4}}");
    /* Ordinary block matrices and scalar padding are unchanged. */
    audit_chk("ArrayFlatten#1", "ArrayFlatten[{{{{1, 2}, {3, 4}}, 0}, {0, {{5}}}}]",
              "{{1, 2, 0}, {3, 4, 0}, {0, 0, 5}}");
    audit_chk("ArrayFlatten#1", "ArrayFlatten[{{{{1, 2}}, {{3}, {4}}}}]",
              "ArrayFlatten[{{{{1, 2}}, {{3}, {4}}}}]");
}

/* Depth#1 (high): "Depth segfaults on a malformed Rule option (Rule[] or
 * Rule[Heads])".  Now Depth::nonopt and the call stays unevaluated. */
static void test_depth_malformed_option(void) {
    audit_chk("Depth#1", "Check[Depth[x, Rule[]], fired, Depth::nonopt]", "fired");
    audit_chk("Depth#1", "Quiet[Depth[x, Rule[]]]", "Depth[x, Rule[]]");
    audit_chk("Depth#1", "Check[Depth[x, Rule[Heads]], fired, Depth::nonopt]", "fired");
    audit_chk("Depth#1", "Quiet[Depth[x, y]]", "Depth[x, y]");
    audit_chk("Depth#1", "Depth[h[{{{a}}}][x, y], Heads -> True]", "6");
    audit_chk("Depth#1", "Depth[h[{{{a}}}][x, y], Heads -> False]", "2");
    audit_chk("Depth#1", "Depth[a -> b]", "2");
}

/* Depth#2 (medium): "Depth crashes (stack overflow) on deeply nested
 * expressions".  The walk is iterative; LeafCount crashed the same way. */
static void test_depth_deep_nesting(void) {
    audit_run("auditDeep = Nest[List, 1, 200000]");
    audit_chk("Depth#2", "Depth[auditDeep]", "200001");
    audit_chk("Depth#2", "LeafCount[auditDeep]", "200001");
    audit_chk("Depth#2", "LeafCount[auditDeep, Heads -> False]", "1");
    audit_run("auditDeepH = Nest[g[#][x] &, y, 100000]");
    audit_chk("Depth#2", "Depth[auditDeepH]", "2");
    audit_chk("Depth#2", "Depth[auditDeepH, Heads -> True]", "200001");
    audit_run("Clear[auditDeep, auditDeepH]");
    audit_chk("Depth#2", "Depth[{{1, 2}, {3, {4}}}]", "4");
    audit_chk("Depth#2", "Depth[<|a -> <|b -> 1|>|>]", "3");
    audit_chk("Depth#2", "Depth[{}]", "2");
    audit_chk("Depth#2", "Depth[{1/2, 1 + I}]", "2");
    audit_chk("Depth#2", "LeafCount[1 + a + b^2]", "6");
}

/* An impossible size aborts with General::nomem: CheckAbort sees the abort and
 * Check sees the (routed) message. */
static void nomem_chk(const char* id, const char* expr) {
    char buf[512];
    snprintf(buf, sizeof buf, "Quiet[CheckAbort[%s, aborted]]", expr);
    audit_chk(id, buf, "aborted");
    snprintf(buf, sizeof buf, "Quiet[Check[CheckAbort[%s, 0], fired, General::nomem]]", expr);
    audit_chk(id, buf, "fired");
}

/* PadLeft#1 (high): a huge machine-integer length crashed the kernel
 * (unchecked malloc; -n overflowed at INT64_MIN).  PadRight shares the
 * driver, the packed rank-1 path and ArrayPad have their own buffers. */
static void test_pad_huge_length(void) {
    nomem_chk("PadLeft#1", "PadLeft[{1}, 10^15]");
    nomem_chk("PadLeft#1", "PadLeft[{1}, 2^62]");
    nomem_chk("PadLeft#1", "PadLeft[{1}, -2^63]");
    nomem_chk("PadLeft#1", "PadRight[{1}, 10^15]");
    nomem_chk("PadLeft#1", "PadRight[{{1}}, {10^10, 10^10}]");
    nomem_chk("PadLeft#1", "PadLeft[Range[100], 2^62]");          /* packed */
    nomem_chk("PadLeft#1", "PadRight[N[Range[100]], 2^62, 0.]");  /* packed */
    nomem_chk("PadLeft#1", "ArrayPad[{1}, 10^15]");
    nomem_chk("PadLeft#1", "ArrayPad[{1}, {0, 2^62}]");
    nomem_chk("PadLeft#1", "ArrayPad[Range[100], {0, 2^62}]");     /* packed */
    audit_chk("PadLeft#1", "PadLeft[{1, 2}, 4]", "{0, 0, 1, 2}");
    audit_chk("PadLeft#1", "PadRight[{1, 2}, -4]", "{0, 0, 1, 2}");
    audit_chk("PadLeft#1", "PadRight[{{1}, {2, 3}}, {2, 3}, x]", "{{1, x, x}, {2, 3, x}}");
    audit_chk("PadLeft#1", "ArrayPad[{1, 2}, 1]", "{0, 1, 2, 0}");
    audit_chk("PadLeft#1", "ArrayPad[Range[5], {-2, 1}]", "{3, 4, 5, 0}");
    audit_chk("PadLeft#1", "ArrayPad[Range[3], -2^63]", "{}");
}

/* Drop#1 (high): Drop[{a,b,c},{1,2^63-1}] segfaulted -- the index count came
 * from the raw spec before the range check (also Take, and the packed path). */
static void test_drop_huge_sequence_spec(void) {
    audit_chk("Drop#1", "Drop[{a, b, c}, {1, 9223372036854775807}]",
              "Drop[{a, b, c}, {1, 9223372036854775807}]");
    audit_chk("Drop#1", "Drop[{a, b, c}, {1, 2^62}]", "Drop[{a, b, c}, {1, 4611686018427387904}]");
    audit_chk("Drop#1", "Take[{a, b, c}, {1, 9223372036854775807}]",
              "Take[{a, b, c}, {1, 9223372036854775807}]");
    audit_chk("Drop#1", "Take[{a, b, c}, {3, -9223372036854775807 - 1, -1}]",
              "Take[{a, b, c}, {3, -9223372036854775808, -1}]");
    audit_chk("Drop#1", "Take[{a, b, c}, {1, -1, -9223372036854775807 - 1}]", "{}");
    audit_chk("Drop#1", "Head[Take[Range[10], {1, 2^63 - 1}]]", "Take");
    audit_chk("Drop#1", "Drop[{a, b, c, d, e}, {1, -1, 2}]", "{b, d}");
    audit_chk("Drop#1", "Take[{a, b, c, d, e}, {5, 1, -2}]", "{e, c, a}");
    audit_chk("Drop#1", "Take[{a, b, c, d, e}, {1, 4, 2}]", "{a, c}");
    audit_chk("Drop#1", "Take[Range[10], {2, 9, 3}]", "{2, 5, 8}");
}

/* Expand#1 (high): Expand[e, Trig -> False] read the option as the pattern
 * argument, so nothing was "free of" it and nothing expanded. The default
 * options are now options; a non-default (unimplemented) one leaves the call. */
static void test_expand_option_not_pattern(void) {
    audit_true("Expand#1", "Expand[(a+b)^2, Trig -> False] === a^2 + 2 a b + b^2");
    audit_true("Expand#1", "Expand[(a+b)^2, Modulus -> 0] === Expand[(a+b)^2]");
    audit_true("Expand#1",
               "Expand[(a+b)^2 (x+1)^2, x, Trig -> False] === Expand[(a+b)^2 (x+1)^2, x]");
    audit_true("Expand#1", "Expand[(a+b)^2 + (x+1)^2, x] === 1 + (a+b)^2 + 2 x + x^2");
}

/* Flatten#1 (high): Flatten[expr] flattens subexpressions with expr's own head
 * (the default head was always List). */
static void test_flatten_default_head(void) {
    audit_true("Flatten#1", "Flatten[f[a, f[b, c], f[f[d]]]] === f[a, b, c, d]");
    audit_true("Flatten#1", "Flatten[f[a, {b}, f[c]], 1] === f[a, {b}, c]");
    /* Lists are unchanged; a List does not flatten a nested f. */
    audit_true("Flatten#1", "Flatten[{a, {b, {c}}, f[{d}]}] === {a, b, c, f[{d}]}");
    audit_true("Flatten#1", "Flatten[{{a}, {{b}}}, 1] === {a, {b}}");
}

/* CoefficientList#1 (high): terms with a non-polynomial power of the variable
 * (1/x, x^n, Sqrt[x], x^2., 1/(1+x)) were silently dropped. Now
 * CoefficientList::poly and unevaluated. */
static void test_coefficientlist_nonpolynomial(void) {
    static const char* bad[] = {
        "1/x + x", "x^n + 1", "Sqrt[x] + 1", "x^2.0 + 1", "1/(1 + x)",
    };
    char buf[160];
    for (size_t i = 0; i < sizeof bad / sizeof bad[0]; i++) {
        snprintf(buf, sizeof buf, "Quiet[Check[CoefficientList[%s, x], fired, CoefficientList::poly]]", bad[i]);
        audit_chk("CoefficientList#1", buf, "fired");
    }
    audit_chk("CoefficientList#1",
              "Quiet[Check[CoefficientList[(1 + x)^2 + y/x, {x, y}], fired, CoefficientList::poly]]",
              "fired");
    audit_true("CoefficientList#1", "Head[Quiet[CoefficientList[1/x + x, x]]] === CoefficientList");
    /* Polynomials, and opaque x-dependent coefficients, are unchanged. */
    audit_true("CoefficientList#1", "CoefficientList[(1 + x)^3, x] === {1, 3, 3, 1}");
    audit_true("CoefficientList#1", "CoefficientList[Sin[x] + x, x] === {Sin[x], 1}");
    audit_true("CoefficientList#1",
               "CoefficientList[(1 + x)^3 + y x, {x, y}] === {{1, 0}, {3, 1}, {3, 0}, {1, 0}}");
}

/* Delete#1, Delete#2 (high): every position names a part of the ORIGINAL
 * expression -- negative indices are resolved first, and a repeated or aliased
 * position deletes once. */
static void test_delete_positions_original(void) {
    audit_true("Delete#1", "Delete[{a, b, c}, {{-3}, {2}}] === {c}");
    audit_true("Delete#1", "Delete[{a, b, c}, {{-3}, {2}}] === Delete[{a, b, c}, {{1}, {2}}]");
    audit_true("Delete#1", "Delete[{a, {b, c}, d}, {{2, -1}, {-1}}] === {a, {b}}");
    audit_true("Delete#2", "Delete[{a, b, c, d}, {{2}, {2}}] === {a, c, d}");
    audit_true("Delete#2", "Delete[{a, b, c}, {{-1}, {3}}] === {a, b}");
    audit_true("Delete#2", "Delete[{a, b, c, d}, {{-2}, {3}, {1}}] === {b, d}");
    /* Ordinary lists of positions are unchanged. */
    audit_true("Delete#1", "Delete[Range[10], {{1}, {3}, {5}}] === {2, 4, 6, 7, 8, 9, 10}");
    audit_true("Delete#1", "Delete[{a, {b, c}, d}, {{2, 1}, {2}}] === {a, d}");
}

/* Insert#1 (high): several positions refer to the original expression, so
 * negative slots are resolved before any insertion. */
static void test_insert_positions_original(void) {
    audit_true("Insert#1", "Insert[{a, b, c}, x, {{-1}, {-2}}] === {a, b, x, c, x}");
    audit_true("Insert#1", "Insert[{a, b, c, d}, x, {{-2}, {-4}}] === {a, x, b, c, x, d}");
    audit_true("Insert#1", "Insert[{a, b, c}, x, {{4}, {-4}}] === {x, a, b, c, x}");
    audit_true("Insert#1", "Insert[{a, b, c, d, e}, x, {{1}, {-1}}] === {x, a, b, c, d, e, x}");
    audit_true("Insert#1", "Insert[{a, {b, c}}, x, {{2, -1}, {-1}}] === {a, {b, c, x}, x}");
    /* A repeated slot inserts once per listing. */
    audit_true("Insert#1", "Insert[{a, b}, x, {{1}, {1}}] === {x, x, a, b}");
}

/* Partition#1 (high): Partition[m, {n1, n2}] is an N1 x N2 grid of n1 x n2
 * blocks (the row index of a block moves below the grid levels). */
static void test_partition_multilevel_blocks(void) {
    audit_true("Partition#1", "Partition[{{1,2,3},{4,5,6}}, {2,2}] === {{{{1, 2}, {4, 5}}}}");
    audit_true("Partition#1", "Dimensions[Partition[{{1,2,3},{4,5,6}}, {2,2}]] === {1, 1, 2, 2}");
    audit_true("Partition#1",
               "Partition[Partition[Range[16], 4], {2, 2}] === "
               "{{{{1, 2}, {5, 6}}, {{3, 4}, {7, 8}}}, {{{9, 10}, {13, 14}}, {{11, 12}, {15, 16}}}}");
    audit_true("Partition#1",
               "Partition[{{1,2,3},{4,5,6},{7,8,9}}, {2,2}, {1,1}][[1]] === "
               "{{{1, 2}, {4, 5}}, {{2, 3}, {5, 6}}}");
    /* Block (i, j, k) of a rank-3 partition is the matching sub-cube. */
    audit_true("Partition#1",
               "With[{p = Partition[Array[a, {4, 4, 4}], {2, 2, 2}]}, "
               "Dimensions[p] === {2, 2, 2, 2, 2, 2} && "
               "p[[2, 1, 2]] === Array[a, {2, 2, 2}, {3, 1, 3}]]");
    /* One-level partitions are unchanged. */
    audit_true("Partition#1", "Partition[Range[7], 3] === {{1, 2, 3}, {4, 5, 6}}");
    audit_true("Partition#1", "Partition[Range[5], 3, 1] === {{1, 2, 3}, {2, 3, 4}, {3, 4, 5}}");
}

/* CoefficientRules#1, MonomialList#1, MonomialList#2 (high): exponents were C
 * ints -- the FLINT read-off truncated x^(2^40) to x^0 and the generic walk
 * wrapped x^2147483647 x to a negative exponent. Now full int64, checked. */
static void test_monomial_exponent_width(void) {
    audit_true("CoefficientRules#1",
               "CoefficientRules[x^(2^40) + 1, x] === {{1099511627776} -> 1, {0} -> 1}");
    audit_true("CoefficientRules#1",
               "CoefficientRules[x^(2^32) + 2, x] === {{4294967296} -> 1, {0} -> 2}");
    audit_true("CoefficientRules#1",
               "FromCoefficientRules[CoefficientRules[x^(2^40) + 1, x], x] === x^(2^40) + 1");
    audit_true("CoefficientRules#1", "MonomialList[x^(2^40) + 1, x] === {x^1099511627776, 1}");
    audit_true("MonomialList#1", "MonomialList[x^2147483647 * x, x] === {x^2147483648}");
    audit_true("MonomialList#1", "CoefficientRules[x^2147483647 * x, x] === {{2147483648} -> 1}");
    audit_true("MonomialList#2",
               "MonomialList[x^2147483647 + x^2147483647 * x^2, x] === {x^2147483649, x^2147483647}");
    audit_true("MonomialList#2",
               "MonomialList[a x^2147483647 + b x^2147483649, x] === {b x^2147483649, a x^2147483647}");
    /* Ordinary monomial orders are unchanged. */
    audit_true("MonomialList#1",
               "CoefficientRules[(x + y)^3, {x, y}] === {{3, 0} -> 1, {2, 1} -> 3, {1, 2} -> 3, {0, 3} -> 1}");
    audit_true("MonomialList#1",
               "MonomialList[x^3 y + x y^5 + 2, {x, y}, \"DegreeReverseLexicographic\"] === {x y^5, x^3 y, 2}");
}

/* Eliminate#1 (high): Gaussian-rational coefficients were "not polynomial".
 * i is adjoined as a parameter g with g^2 + 1 = 0 (smallest in the order). */
static void test_eliminate_gaussian(void) {
    audit_true("Eliminate#1", "Eliminate[{x == I y, y == 2}, y] === (x == 2 I)");
    /* (x - I)^2 == -1  <=>  x^2 == 2 I x */
    audit_true("Eliminate#1", "Eliminate[{x == y + I, y^2 == -1}, y] === (x^2 == 2 I x)");
    audit_true("Eliminate#1", "Eliminate[{x == y, y == 2}, y] === (x == 2)");
}

/* FactorTermsList#1 (high): E^x was no generator at all (collect_variables
 * stops at a numeric base), so the whole sum was the "numerical factor". */
static void test_factortermslist_exp_kernel(void) {
    audit_true("FactorTermsList#1", "FactorTermsList[2 E^x + 4] === {2, 2 + E^x}");
    audit_true("FactorTermsList#1", "FactorTermsList[2 Exp[y] + 4] === {2, 2 + E^y}");
    audit_true("FactorTermsList#1", "FactorTermsList[2 E^(2 x) + 4 E^x] === {2, 2 E^x + E^(2 x)}");
    audit_true("FactorTermsList#1", "FactorTermsList[6 E^x y + 4 y, y] === {2, 2 + 3 E^x, y}");
    /* Unchanged kernels. */
    audit_true("FactorTermsList#1", "FactorTermsList[2 Sin[x] + 4] === {2, 2 + Sin[x]}");
    audit_true("FactorTermsList#1", "FactorTermsList[2 a^x + 4] === {2, 2 + a^x}");
    audit_true("FactorTermsList#1", "FactorTermsList[6 x^2 + 4 x] === {2, 2 x + 3 x^2}");
}

/* Level#1 (high): Level[expr, levelspec, f] ignored f. */
static void test_level_head_argument(void) {
    audit_true("Level#1", "Level[a + b, {1}, g] === g[a, b]");
    audit_true("Level#1", "Level[a + b, {1}, g, Heads -> True] === g[Plus, a, b]");
    audit_true("Level#1", "Level[a + b, {1}, Heads -> True] === {Plus, a, b}");
    audit_true("Level#1", "Level[{{a}, b}, -1] === {a, {a}, b}");
}

/* PolynomialMod#1 (high): a BigInt modulus was ignored (the integer branch
 * tested EXPR_INTEGER only). PolynomialExtendedGCD#1 (high) rides on it: its
 * big-prime fallback reduces through PolynomialMod. */
static void test_polynomialmod_bigint_modulus(void) {
    audit_true("PolynomialMod#1",
               "PolynomialMod[10^25 x + 3, 10^20 + 3] === 3 + Mod[10^25, 10^20 + 3] x");
    audit_true("PolynomialMod#1",
               "PolynomialMod[10^25 x + 3, 10^20 + 3] === 3 + 99999999999999700003 x");
    audit_true("PolynomialMod#1",
               "PolynomialMod[10^25 x + 3, -(10^20 + 3)] === 3 + 99999999999999700003 x");
    audit_true("PolynomialMod#1", "PolynomialMod[x/2 + 3, 2^64 + 13] === 3 + (2^63 + 7) x");
    audit_true("PolynomialMod#1", "PolynomialMod[10^25 x + 3, 7] === 3 + 3 x");
    audit_true("PolynomialMod#1", "PolynomialMod[x^2 + 3 x, x + 1] === -2");
    audit_true("PolynomialExtendedGCD#1",
               "With[{p = 18446744073709551629}, "
               "With[{r = PolynomialExtendedGCD[x^2 + 1, x + 1, x, Modulus -> p]}, "
               "r[[1]] === 1 && FreeQ[r, _Rational] && "
               "PolynomialMod[r[[2, 1]] (x^2 + 1) + r[[2, 2]] (x + 1) - 1, p] === 0]]");
    audit_true("PolynomialExtendedGCD#1",
               "PolynomialExtendedGCD[x^2 + 1, x + 1, x, Modulus -> 7] === {1, {4, 4 + 3 x}}");
}

/* PolynomialMod#2 (medium, run-dependent): the Gaussian branch read both parts'
 * data.integer unchecked; a Rational part (an EXPR_FUNCTION) contributed the
 * bits of a heap pointer, so the real part varied with the allocator.
 * 3/2 = 3 * 2^-1 = 3 * 3 = 4 (mod 5). */
static void test_polynomialmod_gaussian_rational(void) {
    for (int i = 0; i < 5; i++)
        audit_true("PolynomialMod#2", "PolynomialMod[3/2 + 3 I, 5] === 4 + 3 I");
    audit_true("PolynomialMod#2", "PolynomialMod[1/2 + 3 I, 5] === 3 + 3 I");
    audit_true("PolynomialMod#2", "PolynomialMod[(3/2 + 3 I) x + 7, 5] === 2 + (4 + 3 I) x");
    audit_true("PolynomialMod#2", "PolynomialMod[7 + 8 I, 5] === 2 + 3 I");
    audit_true("PolynomialMod#2", "PolynomialMod[x/2 + 1, 5] === 1 + 3 x");
}

/* RankedMax#1: already fixed by v0.6715 (order statistics rank by exact value);
 * pinned. */
static void test_rankedmax_exact_bigints(void) {
    audit_true("RankedMax#1", "RankedMax[{10^30 + 1, 10^30, Pi}, 1] === 10^30 + 1");
    audit_true("RankedMax#1", "RankedMax[{2^53 + 1, 2^53, -Infinity}, -2] === 2^53");
}

/* S14 (Wave 3): level specs and level arguments in the structural heads. */
static void test_s14_levels_and_shapes(void) {
    audit_run("ClearAll[a, b, c, d, e, f, h, x]");
    /* Level#4: -Infinity as a depth bound includes level 0 */
    audit_true("Level#4", "Level[{a, {b, {c}}}, {-Infinity, -1}] === Level[{a, {b, {c}}}, {-4, -1}]");
    audit_true("Level#4", "Level[{a, {b, {c}}}, {-Infinity, -1}] === "
                          "{a, b, c, {c}, {b, {c}}, {a, {b, {c}}}}");
    audit_true("Level#4", "Level[{a, {b}}, {1, -Infinity}] === {}");
    /* Flatten#2: the list-of-level-lists form */
    audit_true("Flatten#2", "Flatten[{{a, b}, {c, d}}, {{2}, {1}}] === {{a, c}, {b, d}}");
    audit_true("Flatten#2", "Flatten[{{a, b}, {c, d}}, {{1, 2}}] === {a, b, c, d}");
    audit_true("Flatten#2", "Flatten[{{a, b, c}, {d, e}, {f}}, {{2}, {1}}] === {{a, d, f}, {b, e}, {c}}");
    audit_true("Flatten#2", "Flatten[Array[h, {2, 3, 4}], {{3}, {1}, {2}}] === "
                            "Transpose[Array[h, {2, 3, 4}], {2, 3, 1}]");
    audit_true("Flatten#2", "Flatten[Array[h, {2, 3, 4}], {{1, 3}, {2}}] === "
                            "Flatten[Transpose[Array[h, {2, 3, 4}], {1, 3, 2}], 1]");
    /* Flatten#4: a symbol or Real level is Flatten::flev */
    audit_chk("Flatten#4", "Check[Flatten[{{a}}, x], fired, Flatten::flev]", "fired");
    audit_chk("Flatten#4", "Check[Flatten[{{a}}, 1.5], fired, Flatten::flev]", "fired");
    audit_true("Flatten#4", "Flatten[{{a}, {{b}}}, 2^70] === {a, b}");
    /* ConjugateTranspose#1: a ragged list is ::nmtx */
    audit_chk("ConjugateTranspose#1",
              "Check[ConjugateTranspose[{{1, I}, {3}}], fired, ConjugateTranspose::nmtx]", "fired");
    audit_chk("ConjugateTranspose#1",
              "Check[ConjugateTranspose[{{1, I}, 3}], fired, ConjugateTranspose::nmtx]", "fired");
    audit_true("ConjugateTranspose#1", "ConjugateTranspose[{1, I}] === {1, -I}");
    audit_true("ConjugateTranspose#1", "ConjugateTranspose[{{1, I}, {2, 3}}] === {{1, 2}, {-I, 3}}");
    /* Transpose#1: m <-> n swaps two levels */
    audit_true("Transpose#1", "Transpose[{{a, b}, {c, d}}, 1 <-> 2] === {{a, c}, {b, d}}");
    audit_true("Transpose#1", "Dimensions[Transpose[Array[h, {2, 3, 4}], 1 <-> 3]] === {4, 3, 2}");
    audit_true("Transpose#1", "Transpose[Array[h, {2, 3, 4}], 1 <-> 3] === "
                              "Transpose[Array[h, {2, 3, 4}], {3, 2, 1}]");
    /* Reverse#2: a bad level is Reverse::ilsmp */
    audit_chk("Reverse#2", "Check[Reverse[{1, 2}, 0], fired, Reverse::ilsmp]", "fired");
    audit_chk("Reverse#2", "Check[Reverse[{1, 2}, -1], fired, Reverse::ilsmp]", "fired");
    audit_chk("Reverse#2", "Check[Reverse[{1, 2}, x], fired, Reverse::ilsmp]", "fired");
    audit_chk("Reverse#2", "Check[Reverse[{1, 2}, 1.5], fired, Reverse::ilsmp]", "fired");
    audit_chk("Reverse#2", "Check[Reverse[{1, 2}, {1, x}], fired, Reverse::ilsmp]", "fired");
    audit_true("Reverse#2", "Reverse[{{1, 2}, {3}}, {1, 2}] === {{3}, {2, 1}}");
    /* Dimensions#2: no 64-level cap */
    audit_chk("Dimensions#2", "Length[Dimensions[Nest[List, 1, 70]]]", "70");
    audit_chk("Dimensions#2", "Length[Dimensions[Nest[List, 1, 70], 66]]", "66");
    audit_true("Dimensions#2", "Dimensions[{{1, 2}, {3}}] === {2}");
    audit_true("Dimensions#2", "Dimensions[{{{1}}}, Infinity] === {1, 1, 1}");
    /* RotateLeft#3 / RotateRight#4: the docstrings name every form */
    audit_chk("RotateLeft#3", "StringContainsQ[Information[\"RotateLeft\"], \"{n1, n2, ...}\"]",
              "True");
    audit_chk("RotateRight#4", "StringContainsQ[Information[\"RotateRight\"], \"RotateRight[expr]\"]",
              "True");
}

/* Wave 3 S12/S13: declared argument counts and option checks
 * (src/arity.h, src/arity_table.c).  The evaluator checks a call
 * against its head's declared arity before the C code runs and
 * prints argx/argr/argrx/argt(u)/argb(u)/argm(u), nonopt or optx;
 * the call stays unevaluated. */
static void test_wave3_arity_options(void) {
    audit_chk("Level#5", "Quiet[Check[Level[a + b], fired, Level::argtu]]", "fired");
    audit_chk("Head#2", "Quiet[Check[Head[x, f, g], fired, Head::argt]]", "fired");
    audit_chk("Head#2", "Quiet[Head[x, f, g]]", "Head[x, f, g]");
    audit_chk("Length#3", "Quiet[Check[Length[a, b], fired, Length::argx]]", "fired");
    audit_chk("Length#3", "Quiet[Length[a, b]]", "Length[a, b]");
    audit_chk("Sort#1", "Quiet[Check[Sort[{3, 1, 2}, Less, 1], fired, Sort::argt]]", "fired");
    audit_chk("Sort#1", "Quiet[Sort[{3, 1, 2}, Less, 1]]", "Sort[{3, 1, 2}, Less, 1]");
    audit_chk("Append#3", "Quiet[Check[Append[{a}, b, c], fired, Append::argt]]", "fired");
    audit_chk("Append#3", "Quiet[Append[{a}, b, c]]", "Append[{a}, b, c]");
    audit_chk("Collect#4", "Quiet[Check[Collect[x], fired, Collect::argtu]]", "fired");
    audit_chk("Collect#4", "Quiet[Collect[x]]", "Collect[x]");
    audit_chk("Delete#7", "Quiet[Check[Delete[{a}, 1, 2], fired, Delete::argt]]", "fired");
    audit_chk("Delete#7", "Quiet[Delete[{a}, 1, 2]]", "Delete[{a}, 1, 2]");
    audit_chk("DeleteDuplicates#3", "Quiet[Check[DeleteDuplicates[{1, 2, 3}, Equal, 3], fired, DeleteDuplicates::argt]]", "fired");
    audit_chk("DeleteDuplicates#3", "Quiet[DeleteDuplicates[{1, 2, 3}, Equal, 3]]", "DeleteDuplicates[{1, 2, 3}, Equal, 3]");
    audit_chk("First#3", "Quiet[Check[First[{1}, 2, 3], fired, First::argt]]", "fired");
    audit_chk("First#3", "Quiet[First[{1}, 2, 3]]", "First[{1}, 2, 3]");
    audit_chk("Insert#6", "Quiet[Check[Insert[{a, b}, x, 1, 2], fired, Insert::argt]]", "fired");
    audit_chk("Insert#6", "Quiet[Insert[{a, b}, x, 1, 2]]", "Insert[{a, b}, x, 1, 2]");
    audit_chk("Most#4", "Quiet[Check[Most[{1, 2}, {3}], fired, Most::argx]]", "fired");
    audit_chk("Most#4", "Quiet[Most[{1, 2}, {3}]]", "Most[{1, 2}, {3}]");
    audit_chk("Partition#5", "Quiet[Check[Partition[{1, 2, 3}], fired, Partition::argbu]]", "fired");
    audit_chk("Partition#5", "Quiet[Partition[{1, 2, 3}]]", "Partition[{1, 2, 3}]");
    audit_chk("Prepend#4", "Quiet[Check[Prepend[{a}, b, c], fired, Prepend::argt]]", "fired");
    audit_chk("Prepend#4", "Quiet[Prepend[{a}, b, c]]", "Prepend[{a}, b, c]");
    audit_chk("Rest#4", "Quiet[Check[Rest[{a, b}, {c}], fired, Rest::argx]]", "fired");
    audit_chk("Rest#4", "Quiet[Rest[{a, b}, {c}]]", "Rest[{a, b}, {c}]");
    audit_chk("Reverse#3", "Quiet[Check[Reverse[{1, 2}, 1, 2], fired, Reverse::argt]]", "fired");
    audit_chk("Reverse#3", "Quiet[Reverse[{1, 2}, 1, 2]]", "Reverse[{1, 2}, 1, 2]");
    audit_chk("Tally#3", "Quiet[Check[Tally[{a}, b, c], fired, Tally::argt]]", "fired");
    audit_chk("Tally#3", "Quiet[Tally[{a}, b, c]]", "Tally[{a}, b, c]");
    audit_chk("Variables#2", "Quiet[Check[Variables[x, y], fired, Variables::argx]]", "fired");
    audit_chk("Variables#2", "Quiet[Variables[x, y]]", "Variables[x, y]");
    audit_chk("AppendTo#2", "Quiet[Check[AppendTo[s12x, 1, 2], fired, AppendTo::argrx]]", "fired");
    audit_chk("AppendTo#2", "Quiet[AppendTo[s12x, 1, 2]]", "AppendTo[s12x, 1, 2]");
    audit_chk("ExpandDenominator#1", "Quiet[Check[ExpandDenominator[], fired, ExpandDenominator::argx]]", "fired");
    audit_chk("ExpandDenominator#1", "Quiet[ExpandDenominator[]]", "ExpandDenominator[]");
    audit_chk("FactorTermsList#3", "Quiet[Check[FactorTermsList[x, y, z], fired, FactorTermsList::argt]]", "fired");
    audit_chk("FactorTermsList#3", "Quiet[FactorTermsList[x, y, z]]", "FactorTermsList[x, y, z]");
    audit_chk("FlattenAt#4", "Quiet[Check[FlattenAt[{a, {b}}], fired, FlattenAt::argr]]", "fired");
    audit_chk("FlattenAt#4", "Quiet[FlattenAt[{a, {b}}]]", "FlattenAt[{a, {b}}]");
    audit_chk("FromCoefficientRules#2", "Quiet[Check[FromCoefficientRules[{{1} -> 1}], fired, FromCoefficientRules::argr]]", "fired");
    audit_chk("FromCoefficientRules#2", "Quiet[FromCoefficientRules[{{1} -> 1}]]", "FromCoefficientRules[{{1} -> 1}]");
    audit_chk("OrderedQ#1", "Quiet[Check[OrderedQ[{1, 2}, Less, 3], fired, OrderedQ::argt]]", "fired");
    audit_chk("OrderedQ#1", "Quiet[OrderedQ[{1, 2}, Less, 3]]", "OrderedQ[{1, 2}, Less, 3]");
    audit_chk("RotateLeft#2", "Quiet[Check[RotateLeft[{a, b, c}, 1, 2], fired, RotateLeft::argt]]", "fired");
    audit_chk("RotateLeft#2", "Quiet[RotateLeft[{a, b, c}, 1, 2]]", "RotateLeft[{a, b, c}, 1, 2]");
    audit_chk("Split#2", "Quiet[Check[Split[{1, 1, 2}, SameQ, 3], fired, Split::argt]]", "fired");
    audit_chk("Split#2", "Quiet[Split[{1, 1, 2}, SameQ, 3]]", "Split[{1, 1, 2}, SameQ, 3]");
    audit_chk("PolynomialGCD#1", "Quiet[Check[PolynomialGCD[x^2 - 1, x - 1, Foo -> 1], fired, PolynomialGCD::optx]]", "fired");
    audit_chk("PolynomialGCD#1", "Quiet[PolynomialGCD[x^2 - 1, x - 1, Foo -> 1]]", "PolynomialGCD[-1 + x^2, -1 + x, Foo -> 1]");
    audit_chk("PolynomialLCM#1", "Quiet[Check[PolynomialLCM[x + 1, x - 1, Foo -> 3], fired, PolynomialLCM::optx]]", "fired");
    audit_chk("PolynomialLCM#1", "Quiet[PolynomialLCM[x + 1, x - 1, Foo -> 3]]", "PolynomialLCM[1 + x, -1 + x, Foo -> 3]");
    audit_chk("GroebnerBasis#2", "Quiet[Check[GroebnerBasis[{x^2 - 1, y}, {x, y}, Foo -> 1], fired, GroebnerBasis::optx]]", "fired");
    audit_chk("IntersectingQ#2", "Quiet[Check[IntersectingQ[{a}, {a}, {a}], fired, IntersectingQ::nonopt]]", "fired");
    audit_chk("IntersectingQ#2", "Quiet[IntersectingQ[{a}, {a}, {a}]]", "IntersectingQ[{a}, {a}, {a}]");
    audit_chk("PolynomialReduce#2", "Quiet[Check[PolynomialReduce[x^2 + 1, {x - 1}, x, Foo -> 1], fired, PolynomialReduce::optx]]", "fired");
    audit_chk("Last#2", "Quiet[Check[Last[{}], fired, Last::nolast]]", "fired");
    audit_chk("Last#2", "Quiet[Last[{}]]", "Last[{}]");
    audit_chk("Last#2", "Quiet[Check[Last[s12a], fired, Last::normal]]", "fired");
    audit_chk("Last#2", "Quiet[Last[s12a]]", "Last[s12a]");
    audit_chk("Last#2", "Quiet[Check[Last[], fired, Last::argt]]", "fired");
    audit_chk("Last#2", "Quiet[Last[]]", "Last[]");
    audit_chk("Append#3", "Quiet[Check[Append[s12a, y], fired, Append::normal]]", "fired");
    audit_chk("Append#3", "Quiet[Append[s12a, y]]", "Append[s12a, y]");
    audit_chk("PolynomialGCD#1", "PolynomialGCD[x^2 - 1, x - 1, Modulus -> 2]", "1 + x");
    audit_chk("FactorList#2", "FactorList[x^2 - 1, {Modulus -> 2}]", "{{1, 1}, {1 + x, 2}}");
    audit_chk("Last#2", "Last[{}, d]", "d");
    audit_chk("Append#3", "Append[f[a], b]", "f[a, b]");
}

/* ================================================================== */
/* Wave 3, shared fix S11: Indeterminate / infinity / pole values     */
/* (src/infinity_rules.c and the per-head rules it does not cover).   */
/* ================================================================== */

static void test_s11_special_values(void) {
    /* ConjugateTranspose#2: conjugates of infinite entries */
    audit_true("ConjugateTranspose#2", "Check[ConjugateTranspose[{{Infinity, I Infinity}}], \"m\"] === {{Infinity}, {-I Infinity}}");
    /* MinMax#2: an Indeterminate element */
    audit_true("MinMax#2", "MinMax[{Indeterminate, 1}] === {Indeterminate, Indeterminate}");
    audit_true("MinMax#2", "Max[1, Indeterminate] === Indeterminate");
    audit_true("MinMax#2", "Min[{2, Indeterminate}] === Indeterminate");
    audit_true("MinMax#2", "MinMax[{3, 1, 2}] === {1, 3}");
}

/* ================================================================== */
/* Wave 4, G2-structural (batch 1): extremum of Intervals, Span and   */
/* Take ranges, argument validation, Listable tdlen settling.         */
/* ================================================================== */

static void test_wave4_structural_batch1(void) {
    /* Min/Max read an Interval argument's extreme endpoint */
    audit_chk("Max#1", "Max[Interval[{1,3}], 5]", "5");
    audit_chk("Max#1", "Max[Interval[{1,3}], Interval[{2,7}]]", "7");
    audit_chk("Max#1", "Min[Interval[{1,3}], 0]", "0");
    audit_chk("Max#1", "Max[Interval[{1,3}]]", "3");
    /* a Span running off the end is Part::take, unevaluated */
    audit_chk("Part#1", "Quiet[Check[{a, b, c, d, e}[[2;;10]], fired, Part::take]]", "fired");
    audit_chk("Part#1", "Quiet[{a, b, c, d, e}[[2;;10]]]", "{a, b, c, d, e}[[Span[2, 10]]]");
    audit_chk("Part#1", "Quiet[Check[<|\"a\" -> 1, \"b\" -> 2|>[[1;;5]], fired, Part::take]]", "fired");
    audit_chk("Part#1", "{a, b, c}[[3;;2]]", "{}");
    audit_chk("Part#1", "{a, b, c}[[2;;3]]", "{b, c}");
    /* Expand's Modulus and Trig options */
    audit_chk("Expand#2", "Expand[(1+x)^3, Modulus -> 2]", "1 + x + x^2 + x^3");
    audit_chk("Expand#3", "Expand[Cos[x + y], Trig -> True]", "Cos[x] Cos[y] - Sin[x] Sin[y]");
    /* Delete of an association entry by key */
    audit_chk("Delete#3", "Delete[<|a->1,b->2|>,Key[a]]", "<|b -> 2|>");
    audit_chk("Delete#3", "Delete[<|\"a\"->1,\"b\"->2|>,\"a\"]", "<|\"b\" -> 2|>");
    audit_chk("Delete#3", "Delete[<|a->1,b->2|>,{Key[a]}]", "<|b -> 2|>");
    audit_chk("DeleteDuplicates#2", "DeleteDuplicates[<|a->1,b->1,c->2|>, Equal]", "<|a -> 1, c -> 2|>");
    /* Listable threading over lists of unequal length: tdlen, unevaluated */
    audit_chk("Exponent#1", "Quiet[Exponent[{x^2, x^3}, {x, y, z}]]", "Exponent[{x^2, x^3}, {x, y, z}]");
    audit_chk("Exponent#1", "Quiet[Check[Exponent[{x^2, x^3}, {x, y, z}], fired, Thread::tdlen]]", "fired");
    audit_chk("Exponent#1", "Quiet[Exponent[{x^2, x^3}, {x, y, z}, List]]", "Exponent[{x^2, x^3}, {x, y, z}, List]");
    audit_chk("PolynomialLCM#2", "Quiet[PolynomialLCM[{x, x^2}, {x^3}]]", "PolynomialLCM[{x, x^2}, {x^3}]");
    /* Insert::psl for a position that is no index */
    audit_chk("Insert#3", "Quiet[Check[Insert[{a, b}, x, y], fired, Insert::psl]]", "fired");
    audit_chk("Insert#3", "Quiet[Insert[{a, b}, x, y]]", "Insert[{a, b}, x, y]");
    audit_chk("Insert#3", "Quiet[Insert[{a, b}, x, 1.5]]", "Insert[{a, b}, x, 1.5]");
    audit_chk("Insert#3", "Quiet[Insert[{a, b}, x, 10^30]]", "Insert[{a, b}, x, 1000000000000000000000000000000]");
    audit_chk("Insert#3", "Quiet[Insert[{a, b}, x, {{1}, 2}]]", "Insert[{a, b}, x, {{1}, 2}]");
    audit_chk("Insert#3", "Quiet[Insert[{a, b}, x, {1, 1}]]", "Insert[{a, b}, x, {1, 1}]");
    audit_chk("Insert#3", "Insert[{a, b}, x, {{1}, {2}}]", "{x, a, x, b}");
    audit_chk("Intersection#1", "Intersection[{1,2,3}, {2,3,4}, SameTest -> Automatic]", "{2, 3}");
    audit_chk("Intersection#1", "Union[{3,1}, {2,1}, SameTest -> Automatic]", "{1, 2, 3}");
    /* Partition::ilsmp for a block size or offset that is not a positive integer */
    audit_chk("Partition#2", "Quiet[Check[Partition[{1,2,3}, 0], fired, Partition::ilsmp]]", "fired");
    audit_chk("Partition#2", "Quiet[Partition[{1,2,3}, x]]", "Partition[{1, 2, 3}, x]");
    audit_chk("Partition#2", "Quiet[Partition[{1,2}, 10^30]]", "Partition[{1, 2}, 1000000000000000000000000000000]");
    audit_chk("Partition#2", "Quiet[Partition[Range[5], UpTo[0]]]", "Partition[{1, 2, 3, 4, 5}, UpTo[0]]");
    audit_chk("Partition#2", "Partition[Range[5], UpTo[2]]", "{{1, 2}, {3, 4}, {5}}");
    audit_chk("Partition#3", "Quiet[Partition[{1,2,3,4}, 2, x]]", "Partition[{1, 2, 3, 4}, 2, x]");
    audit_chk("Partition#3", "Quiet[Check[Partition[{1,2,3,4}, 2, 0], fired, Partition::ilsmp]]", "fired");
    /* PolynomialQuotient of non-polynomials: ::poly, unevaluated */
    audit_chk("PolynomialQuotient#1", "Quiet[PolynomialQuotient[x^2.5, x, x]]", "PolynomialQuotient[x^2.5, x, x]");
    audit_chk("PolynomialQuotient#1", "Quiet[Check[PolynomialQuotient[1/x + x^2, x+1, x], fired, PolynomialQuotient::poly]]", "fired");
    audit_chk("PolynomialQuotient#1", "Quiet[PolynomialQuotient[Sqrt[x] + x^2, x+1, x]]", "PolynomialQuotient[Sqrt[x] + x^2, 1 + x, x]");
    audit_chk("PolynomialQuotient#2", "Quiet[PolynomialQuotient[x^2, Sqrt[x], x]]", "PolynomialQuotient[x^2, Sqrt[x], x]");
    audit_chk("PolynomialQuotient#2", "Quiet[PolynomialQuotient[x^2+1, 1/x, x]]", "PolynomialQuotient[1 + x^2, 1/x, x]");
    audit_chk("PolynomialQuotient#3", "Quiet[PolynomialQuotient[{x^2,x^3}, x+1, x]]", "PolynomialQuotient[{x^2, x^3}, 1 + x, x]");
    audit_chk("PolynomialQuotient#3", "PolynomialQuotient[x^3 + a x, x + b, x]", "a + b^2 - b x + x^2");
    audit_chk("PolynomialQuotient#4", "Check[PolynomialQuotient[x^4+1, 2x+1, x, Modulus->4], fired]", "PolynomialQuotient[1 + x^4, 1 + 2 x, x, Modulus -> 4]");
    audit_chk("PolynomialQuotient#4", "PolynomialQuotient[x^4+1, x+1, x, Modulus->4]", "3 + x + 3 x^2 + x^3");
    /* Reverse::normal for an atom */
    audit_chk("Reverse#1", "Quiet[Check[Reverse[x], fired, Reverse::normal]]", "fired");
    audit_chk("Reverse#1", "Quiet[Reverse[x]]", "Reverse[x]");
    audit_chk("Reverse#1", "Quiet[Reverse[5]]", "Reverse[5]");
    /* Take::take / Drop::drop for positions off the end */
    audit_chk("Take#1", "Quiet[Check[Take[{a,b,c},5], fired, Take::take]]", "fired");
    audit_chk("Take#1", "Quiet[Check[Take[{a,b,c},{-4}], fired, Take::take]]", "fired");
    audit_chk("Take#1", "Quiet[Check[Take[Range[3],5], fired, Take::take]]", "fired");
    audit_chk("Take#1", "Quiet[Take[{a,b,c},5]]", "Take[{a, b, c}, 5]");
    audit_chk("Take#1", "Quiet[Check[Drop[{a,b,c},5], fired, Drop::drop]]", "fired");
    /* Tally of an atom */
    audit_chk("Tally#1", "Quiet[Check[Tally[5], fired, Tally::list]]", "fired");
    audit_chk("Tally#1", "Quiet[Tally[5]]", "Tally[5]");
    /* Union with SameTest compares with every kept element */
    audit_chk("Union#1", "Union[Range[6], SameTest->(Mod[#1-#2,3]==0&)]", "{1, 2, 3}");
    /* the Constant symbols are not variables */
    audit_chk("Variables#1", "Variables[EulerGamma x + Catalan + Degree y]", "{x, y}");
    audit_chk("Variables#1", "Variables[GoldenRatio + z]", "{z}");
    /* Coefficient of a power series */
    audit_chk("Coefficient#1", "Coefficient[Series[Exp[x], {x, 0, 5}], x, 3]", "1/6");
    /* CoefficientList of a list of polynomials, padded to one shape */
    audit_chk("CoefficientList#2", "CoefficientList[{1+x, 1+x^2}, x]", "{{1, 1, 0}, {1, 0, 1}}");
    /* HornerForm::poly through the funnel */
    audit_chk("HornerForm#1", "Quiet[Check[HornerForm[1 + x^a, x], failed]]", "failed");
    /* Resultant: the zero polynomial, and Modulus */
    audit_chk("Resultant#1", "Resultant[0,0,x]", "0");
    audit_chk("Resultant#1", "Resultant[0,x+1,x]", "0");
    audit_chk("Resultant#2", "Resultant[x^2+1, x^3-1, x, Modulus->2]", "0");
    audit_chk("Resultant#2", "Resultant[x^2+1, x^3-1, x, Modulus->3]", "2");
    audit_chk("Resultant#2", "Options[Resultant]", "{Method -> Automatic, Modulus -> 0}");
    /* PolynomialQ: bigint degrees, compound variables */
    audit_chk("PolynomialQ#1", "PolynomialQ[x^(10^20), x]", "True");
    audit_chk("PolynomialQ#1", "PolynomialQ[x^(2^63) + 1, x]", "True");
    audit_chk("PolynomialQ#2", "PolynomialQ[x[1]^2 + x[2], x[1]]", "True");
    audit_chk("PolynomialQ#2", "PolynomialQ[Cos[x] + Sin[x]^2, Sin[x]]", "True");
    audit_chk("PolynomialQ#2", "PolynomialQ[Sin[x]^(1/2), Sin[x]]", "False");
    /* a Throw from the test ends Split at once */
    audit_chk("Split#1", "k = 0; Catch[Split[{1,2,3,4}, (k++; Throw[#1]; True)&]]", "1");
    audit_chk("Split#1", "k", "1");
}

/* Wave 4, G2-structural (batch 2): Gaussian Factor in two variables,
 * Collect on symbolic powers, MinMax padding, PadRight margins, FlattenAt
 * splicing, indexed PrependTo targets, Catenate of associations. */
static void test_wave4_structural_batch2(void) {
    audit_chk("Factor#1", "Factor[2 + 2 y - I x - I x y + x^2 + x^2 y]", "(I + x) (-2*I + x) (1 + y)");
    audit_chk("Factor#1", "Factor[2 - I x + x^2]", "(I + x) (-2*I + x)");
    audit_chk("Factor#1", "Factor[x^2 + 2 I x y - y^2]", "(x + I y)^2");
    audit_chk("Factor#1", "Factor[x^2 + y^2]", "x^2 + y^2");
    audit_chk("Collect#1", "Collect[a x^n + b x^n + c, x]", "c + (a + b) x^n");
    audit_chk("Collect#1", "Collect[a x^n + b x^(2 n) + c x^n, x]", "(a + c) x^n + b x^(2 n)");
    audit_chk("Collect#2", "Collect[a x^(n+1) + b x^(n+1), x]", "(a + b) x^(1 + n)");
    audit_chk("Collect#2", "Collect[a x y + b x y + x, x]", "x (1 + a y + b y)");
    audit_chk("MinMax#1", "MinMax[{1, 2, 3}, 1]", "{0, 4}");
    audit_chk("MinMax#1", "MinMax[{1, 2, 3}, {1, 2}]", "{0, 5}");
    audit_chk("MinMax#1", "MinMax[{1, 2, 3}, Scaled[1/2]]", "{0, 4}");
    audit_chk("PadRight#1", "Quiet[PadRight[{a,b},4,x,m]]", "PadRight[{a, b}, 4, x, m]");
    audit_chk("PadRight#1", "Quiet[Check[PadRight[{a,b},4,x,2.0], fired, PadRight::ilsmp]]", "fired");
    audit_chk("PadRight#1", "PadRight[{a,b},4,x,1]", "{x, a, b, x}");
    audit_chk("FlattenAt#1", "FlattenAt[HoldComplete[{a,b}], 1]", "HoldComplete[a, b]");
    audit_chk("FlattenAt#1", "FlattenAt[{HoldComplete[{a,b}]}, {1,1}]", "{HoldComplete[a, b]}");
    audit_chk("FlattenAt#1", "FlattenAt[{a, {b, c}, {d}}, {{2}, {3}}]", "{a, b, c, d}");
    audit_chk("FlattenAt#2", "Quiet[r4f = {FlattenAt[{a,{b}}, {{}}]}]", "{FlattenAt[{a, {b}}, {{}}]}");
    audit_chk("FlattenAt#2", "Quiet[Check[FlattenAt[{a,{b}}, {{}}], fired, FlattenAt::psl]]", "fired");
    audit_chk("PrependTo#2", "Quiet[Check[PrependTo[h4p[1], z], fired, PrependTo::rvalue]]", "fired");
    audit_chk("PrependTo#2", "h4p[1]", "h4p[1]");
    audit_chk("PrependTo#2", "v4p = {1}; PrependTo[v4p, 0]", "{0, 1}");
    audit_chk("Catenate#1", "Catenate[{<|\"a\" -> 1, \"b\" -> 2|>, <|\"b\" -> 3, \"c\" -> 4|>}]", "{1, 2, 3, 4}");
    audit_chk("PolynomialMod#3", "PolynomialMod[(3+5I) x + 7, 2]", "1 + (1 + I) x");
}

/* Wave 4, G2-structural (batch 3): polynomial argument checks, list forms,
 * ExpandAll denominators and connectives, FactorTerms of radicals,
 * inexact PolynomialReduce, Eliminate of a solved variable. */
static void test_wave4_structural_batch3(void) {
    audit_chk("PolynomialExtendedGCD#3", "PolynomialExtendedGCD[{x^2-1, x^2+x}, x - 1, x]",
              "PolynomialExtendedGCD[{-1 + x^2, x + x^2}, -1 + x, x]");
    audit_chk("PolynomialExtendedGCD#3", "PolynomialExtendedGCD[x^2 - 1, x - 1, x]", "{-1 + x, {0, 1}}");
    audit_chk("PolynomialExtendedGCD#2", "PolynomialExtendedGCD[x^2 + 1, x + 1, x, Modulus -> 4]",
              "PolynomialExtendedGCD[1 + x^2, 1 + x, x, Modulus -> 4]");
    audit_chk("PolynomialExtendedGCD#2", "PolynomialExtendedGCD[x^2 + 1, x + 1, x, Modulus -> 3]", "{1, {2, 2 + x}}");
    audit_chk("CoefficientRules#2", "CoefficientRules[{x + y, x y}, {x, y}]",
              "{{{1, 0} -> 1, {0, 1} -> 1}, {{1, 1} -> 1}}");
    audit_chk("CoefficientRules#2", "CoefficientRules[{x + y, x y}]",
              "{{{1, 0} -> 1, {0, 1} -> 1}, {{1, 1} -> 1}}");
    audit_chk("CoefficientArrays#1", "Quiet[CoefficientArrays[Sin[x] + x, x]]", "CoefficientArrays[x + Sin[x], x]");
    audit_chk("CoefficientArrays#1", "Quiet[Check[CoefficientArrays[Sin[x] + x, x], fired, CoefficientArrays::poly]]", "fired");
    audit_chk("ExpandAll#1", "ExpandAll[1/((x+1)(x+2))]", "1/(2 + 3 x + x^2)");
    audit_chk("ExpandAll#1", "ExpandAll[a/((x+1)(x+2)) + b]", "b + a/(2 + 3 x + x^2)");
    audit_chk("ExpandAll#2", "ExpandAll[f[(x+y)^2] && g[(x+y)^2]]", "f[x^2 + 2 x y + y^2] && g[x^2 + 2 x y + y^2]");
    audit_chk("FactorTerms#1", "FactorTerms[2 Sqrt[x] + 4]", "2 (2 + Sqrt[x])");
    audit_chk("FactorTerms#1", "FactorTerms[6 x^(3/2) + 4 x]", "2 (2 x + 3 x^(3/2))");
    audit_chk("FactorTermsList#2", "FactorTermsList[2 Sqrt[x] + 4]", "{2, 2 + Sqrt[x]}");
    audit_chk("FactorTermsList#2", "FactorTermsList[6 x^(1/3) + 9]", "{3, 3 + 2 x^(1/3)}");
    audit_chk("FactorTermsList#2", "FactorTermsList[4 Sqrt[x] + 6 x^(3/2)]", "{2, 2 Sqrt[x] + 3 x^(3/2)}");
    audit_chk("PolynomialReduce#1", "PolynomialReduce[x^2 + 1.5, {x - 0.5}, x]", "{{0.5 + x}, 1.75}");
    audit_chk("PolynomialReduce#1", "PolynomialReduce[x^2 + 1, {x - 1}, x]", "{{1 + x}, 2}");
    audit_chk("Eliminate#3", "Quiet[Eliminate[{x == E^y, y == 0}, y]]", "x == 1");
    audit_chk("FactorSquareFree#1", "FactorSquareFree[E^(2x)+2E^x+1]", "(1 + E^x)^2");
    audit_chk("FactorSquareFree#1", "FactorSquareFree[x^(2n)+2x^n+1]", "(1 + x^n)^2");
    audit_chk("FactorSquareFree#1", "FactorSquareFree[x^(1/2) + 2 x^(1/4) + 1]", "(1 + x^(1/4))^2");
    audit_chk("FactorSquareFree#1", "FactorSquareFree[x^2 + 2 x + 1]", "(1 + x)^2");
    audit_chk("Eliminate#2", "Quiet[Eliminate[{x == 0.5 y, y == 2}, y]]", "x == 1.0");
}

/* Wave 5: `expr` raises Head::tag through the funnel (Check sees it). */
static void w5s_msg(const char* id, const char* expr, const char* tag) {
    char buf[1024];
    snprintf(buf, sizeof buf, "Quiet[Check[%s, fired, %s]]", expr, tag);
    audit_chk(id, buf, "fired");
}

/* Wave 5 G2: argument counts, ::normal for atoms, ::heads, Part/Extract
 * diagnostics, Most/Rest ::norest, Take/Drop ::seqs/::partd. */
static void test_wave5_structural_messages(void) {
    w5s_msg("Expand#5", "Expand[]", "Expand::argt");
    audit_chk("Expand#4", "Options[Expand]", "{Modulus -> 0, Trig -> False}");
    w5s_msg("Factor#2", "Factor[]", "Factor::argx");
    w5s_msg("Flatten#3", "Flatten[x]", "Flatten::normal");
    audit_chk("Flatten#3", "Quiet[Flatten[x]]", "Flatten[x]");
    audit_chk("Join#1", "Join[]", "{}");
    w5s_msg("Join#3", "Join[{a}, f[b]]", "Join::heads");
    audit_chk("Join#3", "Quiet[Join[{a}, f[b]]]", "Join[{a}, f[b]]");
    w5s_msg("Max#2", "Max[1, 2 I]", "Max::nord");
    w5s_msg("Max#2", "Max[1, 2 + 3 I]", "Max::nord");
    w5s_msg("Min#3", "Min[I, 2]", "Min::nord");
    audit_chk("Min#3", "Quiet[Min[I, 2]]", "Min[I, 2]");
    w5s_msg("Part#2", "{a, b, c}[[5]]", "Part::partw");
    w5s_msg("Part#2", "x[[1]]", "Part::partd");
    w5s_msg("Part#2", "{a, b, c}[[1.5]]", "Part::pkspec1");
    audit_chk("Part#2", "Quiet[{a, b, c}[[5]]]", "{a, b, c}[[5]]");
    w5s_msg("Complement#2", "Complement[{a, b}, f[a]]", "Complement::heads");
    audit_chk("Complement#3", "Complement[{a, b, c}, {b}, SameTest :> SameQ]", "{a, c}");
    w5s_msg("Delete#6", "Delete[{a, b}, x]", "Delete::psl");
    audit_chk("Delete#6", "Quiet[Delete[{a, b}, x]]", "Delete[{a, b}, x]");
    w5s_msg("Drop#2", "Drop[x, 1]", "Drop::normal");
    w5s_msg("Drop#2", "Drop[{a, b, c}, 10]", "Drop::drop");
    w5s_msg("Insert#4", "Insert[x, y, 1]", "Insert::normal");
    w5s_msg("Insert#4", "Insert[\"abc\", \"x\", 2]", "Insert::normal");
    w5s_msg("Intersection#2", "Intersection[{1}, f[1]]", "Intersection::heads");
    w5s_msg("Most#2", "Most[{}]", "Most::norest");
    w5s_msg("Most#3", "Most[x]", "Most::normal");
    w5s_msg("Most#3", "Most[5]", "Most::normal");
    audit_chk("Most#3", "Quiet[Most[x]]", "Most[x]");
    w5s_msg("Rest#2", "Rest[{}]", "Rest::norest");
    w5s_msg("Rest#3", "Rest[x]", "Rest::normal");
    w5s_msg("Partition#4", "Partition[x, 2]", "Partition::normal");
    audit_chk("Partition#4", "Quiet[Partition[\"abcd\", 2]]", "Partition[\"abcd\", 2]");
    w5s_msg("Take#2", "Take[x, 1]", "Take::normal");
    w5s_msg("Take#2", "Take[{a, b, c}, 1.5]", "Take::seqs");
    w5s_msg("Take#2", "Take[{a, b, c}, 1, 1]", "Take::partd");
    audit_chk("Take#2", "Take[{{1, 2}, {3, 4}}, 1, 1]", "{{1}}");
    w5s_msg("Union#3", "Union[{1, 2}, f[3]]", "Union::heads");
    audit_chk("Union#2", "{Union[{0., -0.}], Union[{-0., 0., 0.}]}", "{{0.0}, {0.0}}");
    w5s_msg("Commonest#2", "Commonest[x]", "Commonest::list");
    audit_chk("Commonest#2", "Quiet[Commonest[5]]", "Commonest[5]");
    w5s_msg("Extract#1", "Extract[{a, b}, {3}]", "Extract::partw");
    w5s_msg("Extract#1", "Extract[x, {1}]", "Extract::partd");
    w5s_msg("Extract#1", "Extract[{a, b}, {1, 1}]", "Extract::partd");
    w5s_msg("FlattenAt#3", "FlattenAt[{a, b}, 5]", "FlattenAt::partw");
    audit_chk("FlattenAt#3", "Quiet[FlattenAt[{a, b}, 5]]", "FlattenAt[{a, b}, 5]");
    w5s_msg("FactorTerms#4", "FactorTerms[]", "FactorTerms::argt");
    w5s_msg("FactorTerms#4", "FactorTerms[x, y, z]", "FactorTerms::nonopt");
    w5s_msg("PadLeft#2", "PadLeft[x, 3]", "PadLeft::normal");
    w5s_msg("PolynomialQ#3", "PolynomialQ[]", "PolynomialQ::argt");
    w5s_msg("PolynomialQ#3", "PolynomialQ[x, x, y]", "PolynomialQ::argt");
    w5s_msg("PrependTo#3", "PrependTo[1, 2]", "PrependTo::rvalue");
    w5s_msg("PrependTo#3", "PrependTo[{1, 2}, 3]", "PrependTo::rvalue");
    w5s_msg("RotateLeft#1", "RotateLeft[x]", "RotateLeft::normal");
    w5s_msg("RotateLeft#1", "RotateLeft[\"abc\"]", "RotateLeft::normal");
    w5s_msg("RotateRight#2", "RotateRight[x]", "RotateRight::normal");
    w5s_msg("RotateRight#3", "RotateRight[{a, b, c}, 1.5]", "RotateRight::rspec");
    w5s_msg("RotateRight#3", "RotateRight[{a, b}, 1, 2]", "RotateRight::argt");
    audit_chk("RotateRight#1", "RotateRight[{a, b, c}, -2^63]", "{c, a, b}");
    audit_chk("RotateRight#1", "RotateRight[{a, b, c}, -2^63 + 1]", "{b, c, a}");
    audit_chk("RotateRight#1", "RotateLeft[{a, b, c}, -2^63]", "{b, c, a}");
    w5s_msg("HornerForm#2", "HornerForm[x + 1/x]", "HornerForm::poly");
    audit_chk("HornerForm#2", "Quiet[HornerForm[x + 1/x]]", "HornerForm[1/x + x]");
}

/* Wave 5 G2: polynomial options and edge forms. */
static void test_wave5_structural_algebra(void) {
    audit_chk("ExpandAll#3", "ExpandAll[(a + b)^2, Modulus -> 2]", "a^2 + b^2");
    audit_chk("ExpandAll#3", "ExpandAll[f[(a + b)^2], Trig -> False]", "f[a^2 + 2 a b + b^2]");
    audit_chk("Coefficient#2", "Coefficient[(x + 1)^3, x, 10^20]", "0");
    audit_chk("Coefficient#2", "Coefficient[x^(10^20) + x, x, 10^20]", "1");
    audit_chk("CoefficientList#3", "CoefficientList[x^2 + y, {x, y}, {4, 3}]",
              "{{0, 1, 0}, {0, 0, 0}, {1, 0, 0}, {0, 0, 0}}");
    audit_chk("CoefficientList#3", "CoefficientList[1 + 2 x + 3 x^2, x, {2}]", "{1, 2}");
    audit_chk("Collect#3", "Collect[a x^2.5 + b x^2.5, x]", "(a + b) x^2.5");
    audit_chk("Exponent#2", "{Exponent[x^4, 2 x^2], Exponent[x^4, x^2/2], Exponent[x^4, I x^2]}",
              "{2, 2, 2}");
    audit_chk("PolynomialGCD#2", "PolynomialGCD[3/2, 9/4]", "3/4");
    audit_chk("PolynomialGCD#2", "PolynomialGCD[3/2 x, 9/4 x]", "3/4 x");
    audit_chk("Discriminant#1", "Discriminant[x^2 + 3 x + 1, x, Modulus -> 5]", "0");
    audit_chk("Discriminant#1", "Discriminant[x^2 + 3 x + 1, x, Modulus -> 7]", "5");
    audit_chk("ExpandNumerator#1", "ExpandNumerator[(a + b)^2 (1 + x)^(-100000000000000000000)]",
              "(a^2 + 2 a b + b^2)/(1 + x)^100000000000000000000");
    audit_chk("FactorSquareFreeList#1",
              "FactorSquareFreeList[x^2 - 2 Sqrt[2] x + 2, Extension -> Sqrt[2]]",
              "{{1, 1}, {-Sqrt[2] + x, 2}}");
    audit_chk("FactorSquareFreeList#1",
              "FactorSquareFreeList[x^4 - 4 x^2 + 4, Extension -> Sqrt[2]]",
              "{{1, 1}, {-2 + x^2, 2}}");
    audit_chk("FactorTerms#3", "FactorTerms[(2 x + 2)/(4 y + 4)]", "(1 + x)/(2 (1 + y))");
    audit_chk("FromCoefficientRules#1", "FromCoefficientRules[{{10^20} -> 1}, {x}]",
              "x^100000000000000000000");
    audit_chk("GroebnerBasis#3",
              "Check[GroebnerBasis[{x^2 - 1, y^2 - 1}, {x, y}, Modulus -> 0], fired]",
              "{-1 + y^2, -1 + x^2}");
    w5s_msg("PolynomialExtendedGCD#4", "PolynomialExtendedGCD[x^2 + 1, x + 1, x, Modulus -> a]",
            "PolynomialExtendedGCD::bmod");
    w5s_msg("PolynomialExtendedGCD#6", "PolynomialExtendedGCD[x^2 - 1, x - 1, 3]",
            "PolynomialExtendedGCD::ivar");
    audit_chk("PolynomialExtendedGCD#6", "Quiet[PolynomialExtendedGCD[x^2 - 1, x - 1, 3]]",
              "PolynomialExtendedGCD[-1 + x^2, -1 + x, 3]");
    audit_chk("PowerExpand#2", "PowerExpand[Sqrt[x y], {}]", "Sqrt[x y]");
    audit_chk("PowerExpand#2", "PowerExpand[Log[12], {}]", "Log[12]");
    audit_chk("PowerExpand#2", "PowerExpand[Sqrt[x y], {x}]", "Sqrt[x] Sqrt[y]");
}

/* Wave 5 G2: docstrings. */
static void test_wave5_structural_docs(void) {
    audit_true("Head#3", "StringFreeQ[Information[Head], \"BigInt\"] && Head[2^100] === Integer");
    audit_true("First#4", "StringContainsQ[Information[First], \"First[expr, def]\"] && First[{}, d] === d");
    audit_true("Depth#3", "StringContainsQ[Information[Depth], \"maximum number of indices\"]");
    audit_chk("Prepend#5", "Attributes[Prepend]", "{Protected}");
}

int main(void) {
    audit_init();
    TEST(test_expand_option_not_pattern);
    TEST(test_flatten_default_head);
    TEST(test_coefficientlist_nonpolynomial);
    TEST(test_delete_positions_original);
    TEST(test_insert_positions_original);
    TEST(test_partition_multilevel_blocks);
    TEST(test_monomial_exponent_width);
    TEST(test_eliminate_gaussian);
    TEST(test_factortermslist_exp_kernel);
    TEST(test_level_head_argument);
    TEST(test_polynomialmod_bigint_modulus);
    TEST(test_polynomialmod_gaussian_rational);
    TEST(test_rankedmax_exact_bigints);
    TEST(test_array_flatten_deep_block);
    TEST(test_depth_malformed_option);
    TEST(test_depth_deep_nesting);
    TEST(test_pad_huge_length);
    TEST(test_drop_huge_sequence_spec);
    TEST(test_s14_levels_and_shapes);
    TEST(test_wave3_arity_options);
    TEST(test_s11_special_values);
    TEST(test_wave4_structural_batch1);
    TEST(test_wave4_structural_batch2);
    TEST(test_wave4_structural_batch3);
    TEST(test_wave5_structural_messages);
    TEST(test_wave5_structural_algebra);
    TEST(test_wave5_structural_docs);
    return audit_finish("audit_structural_manipulation");
}
