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
    return audit_finish("audit_structural_manipulation");
}
