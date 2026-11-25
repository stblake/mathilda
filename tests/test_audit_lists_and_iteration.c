/* Builtin-audit regressions: lists and iteration (review/issues/; tests/test_audit.h). */
#include "test_audit.h"

/* UnitVector#1 (critical): UnitVector[opt] read args[1] past the end. */
static void test_unit_vector_option_only(void) {
    audit_chk("UnitVector#1",
              "Check[UnitVector[WorkingPrecision -> MachinePrecision], fired, UnitVector::argt]",
              "fired");
    audit_chk("UnitVector#1", "Quiet[UnitVector[WorkingPrecision -> MachinePrecision]]",
              "UnitVector[WorkingPrecision -> MachinePrecision]");
    audit_chk("UnitVector#1", "UnitVector[3, 2]", "{0, 1, 0}");
    audit_true("UnitVector#1", "UnitVector[2, WorkingPrecision -> MachinePrecision] === {0., 1.}");
    /* A length past memory aborts with General::nomem, never mallocs it. */
    audit_chk("UnitVector#1", "Quiet[CheckAbort[UnitVector[10^15, 1], aborted]]", "aborted");
    audit_chk("UnitVector#1",
              "Quiet[Check[CheckAbort[UnitVector[10^15, 1], 0], fired, General::nomem]]",
              "fired");
}

/* Tuples#1 (high): "Tuples[{a}, 10^6] crashes with SIGSEGV (unbounded
 * recursion in tuples_rec)".  The index tuples come from an odometer loop. */
static void test_tuples_long_tuple(void) {
    audit_chk("Tuples#1", "Length[Tuples[{a}, 10^6]]", "1");
    audit_chk("Tuples#1", "Length[First[Tuples[{a}, 10^6]]]", "1000000");
    audit_true("Tuples#1", "First[Tuples[{a}, 10^6]] === ConstantArray[a, 10^6]");
    /* Lexicographic order and the other forms are unchanged. */
    audit_chk("Tuples#1", "Tuples[{0, 1}, 3]",
              "{{0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1}, {1, 0, 0}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1}}");
    audit_chk("Tuples#1", "Tuples[{{a, b}, {1, 2, 3}}]",
              "{{a, 1}, {a, 2}, {a, 3}, {b, 1}, {b, 2}, {b, 3}}");
    audit_chk("Tuples#1", "Tuples[{{a, b}, {}}]", "{}");
    audit_chk("Tuples#1", "Tuples[{}]", "{{}}");
    audit_chk("Tuples#1", "Tuples[{a, b}, {2, 1}]",
              "{{{a}, {a}}, {{a}, {b}}, {{b}, {a}}, {{b}, {b}}}");
}

/* Array#1 (high): Array[f, 10^15] segfaulted (unchecked malloc of the result
 * slots, which also wrapped at 2^61).  Every level's slot count is checked. */
static void test_array_huge_dimension(void) {
    audit_chk("Array#1", "Quiet[CheckAbort[Array[f, 10^15], aborted]]", "aborted");
    audit_chk("Array#1", "Quiet[Check[CheckAbort[Array[f, 10^15], 0], fired, General::nomem]]",
              "fired");
    audit_chk("Array#1", "Quiet[CheckAbort[Array[f, 2^61], aborted]]", "aborted");
    audit_chk("Array#1", "Quiet[CheckAbort[Array[f, {10^15, 0}], aborted]]", "aborted");
    audit_chk("Array#1", "Quiet[CheckAbort[Array[f, {2^32, 2^32}], aborted]]", "aborted");
    audit_chk("Array#1", "Array[f, 3]", "{f[1], f[2], f[3]}");
    audit_chk("Array#1", "Array[f, {2, 2}, 0]", "{{f[0, 0], f[0, 1]}, {f[1, 0], f[1, 1]}}");
    audit_chk("Array#1", "Array[f, {0, 10^15}]", "{}");
    /* The same class in ArrayReshape (List and packed paths). */
    audit_chk("Array#1", "Quiet[CheckAbort[ArrayReshape[{1}, {10^15, 0}], aborted]]", "aborted");
    audit_chk("Array#1", "Quiet[CheckAbort[ArrayReshape[Range[100], 2^61], aborted]]", "aborted");
    audit_chk("Array#1", "ArrayReshape[Range[6], {2, 3}]", "{{1, 2, 3}, {4, 5, 6}}");
    audit_chk("Array#1", "ArrayReshape[{1, 2}, {2, 2}]", "{{1, 2}, {0, 0}}");
}

/* Range#1 (high): "Range[0., 100., 0.01] drops the final element (10000
 * instead of 10001 elements)".  The count is Floor[(imax - imin)/di] + 1 with
 * a relative slack, and element k is imin + k di (never past imax). */
static void test_range_real_count(void) {
    audit_chk("Range#1", "Length[Range[0., 100., 0.01]]", "10001");
    audit_true("Range#1", "Last[Range[0., 100., 0.01]] === 100.");
    audit_true("Range#1", "Range[0., 0.3, 0.1] === {0., 0.1, 0.2, 0.3}");
    audit_true("Range#1", "Last[Range[0., 1., 0.1]] === 1.");
    audit_true("Range#1", "Range[0.1, 1., 0.3] === {0.1, 0.4, 0.7, 0.1 + 3*0.3}");
    audit_true("Range#1", "Length[Range[0., 0.2999999999999, 0.1]] === 3");
    audit_true("Range#1", "Range[1., 0., -0.5] === {1., 0.5, 0.}");
    audit_true("Range#1", "Range[1, 2, 0.25] === {1., 1.25, 1.5, 1.75, 2.}");
    audit_true("Range#1", "Range[0, 1, 1/3] === {0, 1/3, 2/3, 1}");
    audit_true("Range#1", "Range[5., 1., 1.] === {}");
    audit_true("Range#1", "Length[Range[-1.*^6, 0., 0.1]] === 10000001");
}

/* ArrayPad#1 (high): "ArrayPad silently erases the array's data when the
 * per-level spec has more levels than the array". */
static void test_arraypad_spec_too_deep(void) {
    audit_chk("ArrayPad#1", "Check[ArrayPad[{1, 2}, {{1, 1}, {1, 1}}], fired, ArrayPad::depth]", "fired");
    audit_true("ArrayPad#1",
               "MatchQ[Quiet[ArrayPad[{1, 2}, {{1, 1}, {1, 1}}]], HoldPattern[ArrayPad[{1, 2}, {{1, 1}, {1, 1}}]]]");
    audit_chk("ArrayPad#1",
              "Check[ArrayPad[{{1, 2}, {3, 4}}, {{1, 1}, {1, 1}, {1, 1}}], fired, ArrayPad::depth]", "fired");
    /* a ragged level is not an array level */
    audit_chk("ArrayPad#1", "Check[ArrayPad[{1, {2, 3}}, {{1, 1}, {1, 1}}], fired, ArrayPad::depth]", "fired");
    audit_true("ArrayPad#1", "ArrayPad[{1, 2}, 1] === {0, 1, 2, 0}");
    audit_true("ArrayPad#1", "ArrayPad[{1, 2}, {{1, 2}}] === {0, 1, 2, 0, 0}");
    audit_true("ArrayPad#1",
               "ArrayPad[{{1, 2}, {3, 4}}, {{1, 0}, {0, 1}}] === {{0, 0, 0}, {1, 2, 0}, {3, 4, 0}}");
    audit_true("ArrayPad#1",
               "ArrayPad[{{1, 2}, {3, 4}}, 1] === {{0, 0, 0, 0}, {0, 1, 2, 0}, {0, 3, 4, 0}, {0, 0, 0, 0}}");
}

/* Permutations#1 (high): "Permutations[list, n] returns permutations of
 * exactly n elements instead of at most n". */
static void test_permutations_at_most_n(void) {
    audit_chk("Permutations#1", "Permutations[{a, b, c}, 2]",
              "{{}, {a}, {b}, {c}, {a, b}, {a, c}, {b, a}, {b, c}, {c, a}, {c, b}}");
    audit_chk("Permutations#1", "Permutations[{a, b}, Infinity]", "{{}, {a}, {b}, {a, b}, {b, a}}");
    audit_chk("Permutations#1", "Permutations[{a, a, b}, 2]",
              "{{}, {a}, {b}, {a, a}, {a, b}, {b, a}}");
    /* {n} is the exactly-n form */
    audit_chk("Permutations#1", "Permutations[{a, b, c}, {2}]",
              "{{a, b}, {a, c}, {b, a}, {b, c}, {c, a}, {c, b}}");
    audit_chk("Permutations#1", "Permutations[{a, b}, 0]", "{{}}");
    audit_chk("Permutations#1", "Length[Permutations[Range[4], 9]]", "65");
    audit_chk("Permutations#1", "Permutations[{a, b}]", "{{a, b}, {b, a}}");
}

/* ArrayPad#2 (medium): amounts for fewer levels than the rank pad those
 * levels with full-shape slices. */
static void test_s14_arraypad_levels(void) {
    audit_true("ArrayPad#2", "ArrayPad[{{1, 2}, {3, 4}}, {{1, 0}}] === {{0, 0}, {1, 2}, {3, 4}}");
    audit_true("ArrayPad#2", "ArrayPad[{{1, 2}, {3, 4}}, {{1, 1}}, 9] === "
                             "{{9, 9}, {1, 2}, {3, 4}, {9, 9}}");
    audit_true("ArrayPad#2", "ArrayPad[{{1, 2}, {3, 4}}, {{0, 0}, {1, 0}}] === {{0, 1, 2}, {0, 3, 4}}");
}

/* Wave 3 S12/S13: declared argument counts and option checks
 * (src/arity.h, src/arity_table.c).  The evaluator checks a call
 * against its head's declared arity before the C code runs and
 * prints argx/argr/argrx/argt(u)/argb(u)/argm(u), nonopt or optx;
 * the call stays unevaluated. */
static void test_wave3_arity_options(void) {
    audit_chk("GroupElements#2", "Quiet[Check[GroupElements[], fired, GroupElements::argt]]", "fired");
    audit_chk("GroupElements#2", "Quiet[GroupElements[]]", "GroupElements[]");
    audit_chk("PermutationCycles#2", "Quiet[Check[PermutationCycles[{2, 1}, h, 3], fired, PermutationCycles::argt]]", "fired");
    audit_chk("PermutationCycles#2", "Quiet[PermutationCycles[{2, 1}, h, 3]]", "PermutationCycles[{2, 1}, h, 3]");
    audit_chk("Signature#1", "Quiet[Check[Signature[{1, 2}, {3}], fired, Signature::argx]]", "fired");
    audit_chk("Signature#1", "Quiet[Signature[{1, 2}, {3}]]", "Signature[{1, 2}, {3}]");
    audit_chk("UnitVector#3", "Quiet[Check[UnitVector[3, 1, Foo -> 2], fired, UnitVector::optx]]", "fired");
    audit_chk("UnitVector#3", "Quiet[UnitVector[3, 1, Foo -> 2]]", "UnitVector[3, 1, Foo -> 2]");
    audit_chk("UnitVector#2", "UnitVector[3, 1, {WorkingPrecision -> MachinePrecision}]", "{1.0, 0.0, 0.0}");
}

/* ================================================================== */
/* Wave 3, shared fix S11: Indeterminate / infinity / pole values     */
/* (src/infinity_rules.c and the per-head rules it does not cover).   */
/* ================================================================== */

static void test_s11_special_values(void) {
    /* CosineDistance#1: an Infinity component gives no finite distance */
    audit_true("CosineDistance#1", "Quiet[CosineDistance[{1, Infinity}, {1, 2}]] === Indeterminate");
    /* EuclideanDistance#2: infinite components */
    audit_true("EuclideanDistance#2", "EuclideanDistance[{Infinity, 0}, {0, 0}] === Infinity");
    /* EuclideanDistance#3: no spurious message for -Infinity */
    audit_true("EuclideanDistance#3", "Check[EuclideanDistance[{-Infinity}, {0}], \"m\"] === Infinity");
    /* EuclideanDistance#4: an Indeterminate component */
    audit_true("EuclideanDistance#4", "EuclideanDistance[{Indeterminate}, {0}] === Indeterminate");
    /* ManhattanDistance#2: infinite components */
    audit_true("ManhattanDistance#2", "ManhattanDistance[{Infinity, 0}, {0, 0}] === Infinity");
    audit_true("ManhattanDistance#2", "Check[ManhattanDistance[{-Infinity}, {0}], \"m\"] === Infinity");
    /* SquaredEuclideanDistance#2: infinite components */
    audit_true("SquaredEuclideanDistance#2", "SquaredEuclideanDistance[{Infinity, 1}, {0, 1}] === Infinity");
    audit_true("SquaredEuclideanDistance#2", "SquaredEuclideanDistance[{ComplexInfinity}, {0}] === Infinity");
}

/* ================================================================== */
/* Wave 4, G2-structural (lists and iteration): symbolic lattices,    */
/* Table iterator validation, scalar distances, named groups, ...     */
/* ================================================================== */

static void test_wave4_lists_iteration(void) {
    /* symbolic bounds a numeric distance apart */
    audit_chk("Range#3", "Range[a, a+3]", "{a, 1 + a, 2 + a, 3 + a}");
    audit_chk("Range#3", "Range[0, 3 b, b]", "{0, b, 2 b, 3 b}");
    audit_chk("Range#3", "Range[x]", "Range[x]");
    audit_chk("Table#3", "Table[i, {i, a, a+3}]", "{a, 1 + a, 2 + a, 3 + a}");
    audit_chk("Table#3", "Table[i^2, {i, a, a + 1}]", "{a^2, (1 + a)^2}");
    /* a bad iterator anywhere leaves the whole multi-iterator Table */
    audit_chk("Table#1", "Quiet[Table[i, {i, 2}, {j, 1, 2, 0}]]", "Table[i, {i, 2}, {j, 1, 2, 0}]");
    audit_chk("Table#1", "Quiet[Check[Table[i, {i, 2}, {j, 1, 2, 0}], fired, Table::iterb]]", "fired");
    audit_chk("Table#1", "Quiet[Table[x, {i, 2}, a]]", "Table[x, {i, 2}, a]");
    audit_chk("Table#1", "Quiet[Check[Table[x, {i, 2}, a], fired, Table::itform]]", "fired");
    audit_chk("Table#1", "Quiet[Table[i, {i, 1, 3}, {}]]", "Table[i, {i, 1, 3}, {}]");
    audit_chk("Table#1", "Table[i j, {i, 3}, {j, i}]", "{{1}, {2, 4}, {3, 6, 9}}");
    audit_chk("Table#1", "Table[x, {i, 2}, 3]", "{{x, x, x}, {x, x, x}}");
    audit_chk("Table#2", "Quiet[Table[i, i, {i, 3}]]", "Table[i, i, {i, 3}]");
    audit_chk("Table#2", "Quiet[Table[Sequence[i, i], {i, 3}]]", "Table[i, i, {i, 3}]");
    /* Array with a head */
    audit_chk("Array#2", "Array[f, 3, 1, g]", "g[f[1], f[2], f[3]]");
    audit_chk("Array#2", "Array[f, {2, 2}, 1, g]", "g[g[f[1, 1], f[1, 2]], g[f[2, 1], f[2, 2]]]");
    audit_chk("Array#2", "Array[#^2 &, 5, 1, Plus]", "55");
    /* ArrayReshape reads only a List as an array */
    audit_chk("ArrayReshape#1", "Quiet[ArrayReshape[f[a,b,c,d],{2,2}]]", "ArrayReshape[f[a, b, c, d], {2, 2}]");
    audit_chk("ArrayReshape#1", "ArrayReshape[{f[a,b],c},{3}]", "{f[a, b], c, 0}");
    /* Cycles past the image-array degree cap */
    audit_chk("Cycles#1", "Cycles[{{10^9,5}}]", "Cycles[{{5, 1000000000}}]");
    audit_chk("Cycles#1", "Quiet[Check[Cycles[{{10^9,5},{5,7}}], fired, Cycles::reppoint]]", "fired");
    audit_chk("Cycles#1", "Cycles[{{10^9,5},{3,7,2}}]", "Cycles[{{2, 3, 7}, {5, 1000000000}}]");
    /* distances between scalars */
    audit_chk("EuclideanDistance#1", "EuclideanDistance[-3, 7/2]", "13/2");
    audit_chk("EuclideanDistance#1", "EuclideanDistance[I, 0]", "1");
    audit_chk("EuclideanDistance#1", "EuclideanDistance[1.5+2.I, 0]", "2.5");
    audit_chk("EuclideanDistance#1", "EuclideanDistance[{1,2},{4,6}]", "5");
    audit_chk("ManhattanDistance#1", "ManhattanDistance[1/2, 1/3]", "1/6");
    audit_chk("ManhattanDistance#1", "ManhattanDistance[I, 0]", "1");
    audit_chk("SquaredEuclideanDistance#1", "SquaredEuclideanDistance[1/2, 0]", "1/4");
    audit_chk("SquaredEuclideanDistance#1", "SquaredEuclideanDistance[1+I, 0]", "2");
    audit_chk("SquaredEuclideanDistance#1", "SquaredEuclideanDistance[a + b, c]", "Abs[a + b - c]^2");
    audit_chk("SquaredEuclideanDistance#1", "SquaredEuclideanDistance[Sqrt[2], 0]", "2");
    /* the named permutation groups */
    audit_chk("GroupElements#1", "GroupElements[SymmetricGroup[3]]",
              "{Cycles[{}], Cycles[{{2, 3}}], Cycles[{{1, 2}}], Cycles[{{1, 2, 3}}], Cycles[{{1, 3, 2}}], Cycles[{{1, 3}}]}");
    audit_chk("GroupElements#1", "GroupElements[CyclicGroup[4]]",
              "{Cycles[{}], Cycles[{{1, 2, 3, 4}}], Cycles[{{1, 3}, {2, 4}}], Cycles[{{1, 4, 3, 2}}]}");
    audit_chk("GroupElements#1", "GroupElements[AlternatingGroup[3]]",
              "{Cycles[{}], Cycles[{{1, 2, 3}}], Cycles[{{1, 3, 2}}]}");
    audit_chk("GroupElements#1", "Length[GroupElements[DihedralGroup[3]]]", "6");
    audit_chk("GroupElements#1", "Length[GroupElements[AlternatingGroup[4]]]", "12");
    audit_chk("GroupOrder#1", "GroupOrder[SymmetricGroup[4]]", "24");
    audit_chk("GroupOrder#1", "GroupOrder[AlternatingGroup[5]]", "60");
    audit_chk("GroupOrder#1", "GroupOrder[CyclicGroup[7]]", "7");
    audit_chk("GroupOrder#1", "GroupOrder[DihedralGroup[5]]", "10");
    audit_chk("GroupOrder#1", "GroupOrder[SymmetricGroup[25]]", "15511210043330985984000000");
    /* PermutationCycles with a head */
    audit_chk("PermutationCycles#1", "PermutationCycles[{2,3,1}, h]", "h[{{1, 2, 3}}]");
    /* Permutations::nninfseq */
    audit_chk("Permutations#2", "Quiet[Permutations[{a,b,c},x]]", "Permutations[{a, b, c}, x]");
    audit_chk("Permutations#2", "Quiet[Check[Permutations[{a,b,c},1.5], fired, Permutations::nninfseq]]", "fired");
    audit_chk("Permutations#2", "Quiet[Permutations[{a,b,c},{-1}]]", "Permutations[{a, b, c}, {-1}]");
    audit_chk("Permutations#2", "Quiet[Permutations[{a,b,c},{1,3,0}]]", "Permutations[{a, b, c}, {1, 3, 0}]");
    audit_chk("Permutations#2", "Length[Permutations[{a,b,c},2^70]]", "16");
    /* Riffle's three-argument form (already fixed) */
    audit_chk("Riffle#1", "Riffle[{a,b,c},x,2]", "{a, x, b, x, c}");
    audit_chk("Riffle#1", "Riffle[{a,b,c,d,e},x,{1,-1,2}]", "{x, a, x, b, x, c, x, d, x, e, x}");
    /* Tuples with a zero dimension */
    audit_chk("Tuples#2", "Tuples[{a,b},{0}]", "{{}}");
    audit_chk("Tuples#2", "Tuples[{a,b},{2,0}]", "{{{}, {}}}");
}

/* Wave 5: `expr` raises Head::tag through the funnel (Check sees it). */
static void w5l_msg(const char* id, const char* expr, const char* tag) {
    char buf[1024];
    snprintf(buf, sizeof buf, "Quiet[Check[%s, fired, %s]]", expr, tag);
    audit_chk(id, buf, "fired");
}

/* Wave 5 G2: iterator bounds, real ranges, sequence distances, specs. */
static void test_wave5_lists_batch(void) {
    w5l_msg("Range#4", "Range[x]", "Range::range");
    w5l_msg("Range#4", "Range[1, 2, 0]", "Range::range");
    audit_chk("Range#4", "Quiet[Range[1, 2, 0]]", "Range[1, 2, 0]");
    w5l_msg("Table#4", "Table[i, {i, a, b}]", "Table::iterb");
    w5l_msg("Table#4", "Table[i, {i, 1, 10, 0}]", "Table::iterb");
    audit_chk("Table#4", "Quiet[Table[i, {i, 1, 10, 0}]]", "Table[i, {i, 1, 10, 0}]");
    audit_chk("Table#5", "Last[Table[i, {i, 0., 1., 0.01}]]", "1.0");
    audit_true("Table#5", "Table[i, {i, 1, 2, 0.1}] === Range[1, 2, 0.1]");
    audit_true("Table#5", "Table[i, {i, 0., 1., 0.01}] === Range[0., 1., 0.01]");
    audit_true("Table#5", "Table[Sin[i], {i, 0., 1., 0.25}] === Sin[Range[0., 1., 0.25]]");
    audit_true("Table#6", "StringContainsQ[Information[Table], \"Table[expr, {i, imin, imax, di}]\"]");
    audit_chk("ConstantArray#1", "ConstantArray[c, {}]", "c");
    w5l_msg("ConstantArray#2", "ConstantArray[c, -1]", "ConstantArray::ilsmn");
    w5l_msg("ConstantArray#2", "ConstantArray[c, 2.0]", "ConstantArray::ilsmn");
    w5l_msg("ConstantArray#2", "ConstantArray[c, 10^20]", "ConstantArray::ilsmn");
    audit_chk("CosineDistance#2", "CosineDistance[\"abc\", \"abd\"]", "CosineDistance[\"abc\", \"abd\"]");
    w5l_msg("Cycles#2", "Cycles[{{1, -10^20}}]", "Cycles::pospoint");
    w5l_msg("Cycles#3", "Cycles[{{a, -1}}]", "Cycles::pospoint");
    w5l_msg("Cycles#3", "Cycles[{{-1, a}}]", "Cycles::pospoint");
    audit_chk("EditDistance#2", "EditDistance[\"abc\", \"ABC\", IgnoreCase -> True]", "0");
    audit_chk("EditDistance#2", "EditDistance[\"abc\", \"ABC\"]", "3");
    audit_chk("HammingDistance#2", "HammingDistance[\"abc\", \"ABC\", IgnoreCase -> True]", "0");
    w5l_msg("HammingDistance#3", "HammingDistance[{1, 2}, {1, 2, 3}]", "HammingDistance::idim");
    audit_chk("HammingDistance#3", "Quiet[HammingDistance[{1, 2}, {1, 2, 3}]]",
              "HammingDistance[{1, 2}, {1, 2, 3}]");
    audit_chk("FindClusters#2", "FindClusters[{1 -> a, 2 -> b, 10 -> c}]", "{{a, b}, {c}}");
    audit_chk("FindClusters#2", "FindClusters[{1, 2, 10} -> {a, b, c}]", "{{a, b}, {c}}");
    audit_chk("Nearest#1", "Nearest[{1/10^25, 1}, 0]", "{1/10000000000000000000000000}");
    w5l_msg("Permutations#3", "Permutations[x]", "Permutations::normal");
    w5l_msg("Permutations#3", "Permutations[{a, b}, 1, 2]", "Permutations::argt");
    audit_chk("Subdivide#1", "Subdivide[{1, 2}, 2]", "{{0, 0}, {1/2, 1}, {1, 2}}");
    audit_true("Subdivide#1", "Subdivide[{1, 2}, 2] === Subdivide[{0, 0}, {1, 2}, 2]");
    w5l_msg("Tuples#3", "Tuples[{a, b}, -1]", "Tuples::intnm");
    audit_chk("Tuples#3", "Quiet[Tuples[{a, b}, -1]]", "Tuples[{a, b}, -1]");
    audit_chk("UnitVector#5",
              "UnitVector[3, 1, WorkingPrecision -> MachinePrecision, WorkingPrecision -> Infinity]",
              "{1.0, 0.0, 0.0}");
}

int main(void) {
    audit_init();
    TEST(test_unit_vector_option_only);
    TEST(test_tuples_long_tuple);
    TEST(test_array_huge_dimension);
    TEST(test_range_real_count);
    TEST(test_arraypad_spec_too_deep);
    TEST(test_permutations_at_most_n);
    TEST(test_s14_arraypad_levels);
    TEST(test_wave3_arity_options);
    TEST(test_s11_special_values);
    TEST(test_wave4_lists_iteration);
    TEST(test_wave5_lists_batch);
    return audit_finish("audit_lists_and_iteration");
}
