/* test_sparse_dot.c -- native Dot / MatrixPower / Inner over SparseArray
 * (src/sparse/sparse_dot.c).
 *
 * Covers the answer KIND of every shape (dense, sparse, scalar), the
 * background of a sparse product, the machine / exact / boxed value paths
 * (including int64 overflow), shape errors, chains, and -- at sizes that take
 * the threaded kernels -- agreement with the dense product.
 *
 * Failures are counted rather than asserted: the CMake Release build passes
 * -DNDEBUG, which would no-op assert(). */

#include "test_utils.h"
#include "symtab.h"
#include "core.h"

static int failures = 0;

/* Evaluate `input`; fail unless it prints exactly `expected`. */
static void chk(const char* input, const char* expected) {
    struct Expr* p = parse_expression(input);
    if (!p) { fprintf(stderr, "FAIL(parse): %s\n", input); failures++; return; }
    struct Expr* e = evaluate(p);
    expr_free(p);
    char* s = expr_to_string(e);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  Expected: %s\n  Actual:   %s\n", input, expected, s);
        failures++;
    }
    free(s);
    expr_free(e);
}

/* The InputForm (canonical form) of a SparseArray-valued input. */
static void chk_if(const char* input, const char* expected) {
    char buf[4096];
    snprintf(buf, sizeof buf, "InputForm[%s]", input);
    chk(buf, expected);
}

static void chk_true(const char* input) { chk(input, "True"); }

/* ------------------------------------------------------------------------ */

static void setup(void) {
    chk("s = SparseArray[{{1,1}->2, {2,3}->3, {3,2}->-1}, {3,3}]; "
        "r = SparseArray[{{1,1}->1.5, {3,3}->2.}, {3,3}]; "
        "sv = SparseArray[{1->2, 3->5}, 3]; "
        "z = SparseArray[{}, {3,3}]; "
        "c = SparseArray[{{1,2}->1.+2.I}, {3,3}]; Null", "Null");
}

/* s . v, v . s, s . m, m . s: dense answers */
static void test_dense_answers(void) {
    chk("s . {1, 2, 3}", "{2, 9, -2}");
    chk("{1, 2, 3} . s", "{2, -3, 6}");
    chk("s . {{1, 2}, {3, 4}, {5, 6}}", "{{2, 4}, {15, 18}, {-3, -4}}");
    chk("{{1, 2, 3}} . s", "{{2, -3, 6}}");
    chk("{{1, 2, 3}, {0, 1, 0}} . s", "{{2, -3, 6}, {0, 0, 3}}");
    chk("s . {1., 2., 3.}", "{2.0, 9.0, -2.0}");
    /* a machine operand makes every element machine, empty rows included */
    chk("SparseArray[{{1,1}->2},{2,2}] . {1., 2.}", "{2.0, 0.0}");
    chk("{1., 2.} . SparseArray[{{1,1}->2},{2,2}]", "{2.0, 0.0}");
    chk("SparseArray[{{1,1}->2.},{2,2}] . {{1},{2}}", "{{2.0}, {0.0}}");
    chk("Head[s . {1, 2, 3}]", "List");
    /* the answer is a packed buffer; a visible NDArray stays visible */
    chk("NDArrayQ[s . {1., 2., 3.}]", "True");
    chk("Head[s . NDArray[{1., 2., 3.}]]", "NDArray");
    /* exact numbers that are not Integers meet a machine number: machine */
    chk("SparseArray[{{1,1}->1/2},{2,2}] . {1., 2.}", "{0.5, 0.0}");
    chk("SparseArray[{{1,1}->2.},{2,2}] . {1/2, 2}", "{1.0, 0.0}");
    /* rank-1 sparse against a matrix, and a matrix against a sparse vector */
    chk("sv . {{1, 2}, {3, 4}, {5, 6}}", "{27, 34}");
    chk("{{1, 2, 3}, {4, 5, 6}} . sv", "{17, 38}");
    /* complex */
    chk("c . {1, 2, 3} == {2. + 4. I, 0, 0}", "True");
}

/* s1 . s2 and the vector forms: sparse answers */
static void test_sparse_answers(void) {
    chk_if("s . s", "SparseArray[Automatic, {3, 3}, 0, {1, {{0, 1, 2, 3}, {{1}, {2}, {3}}}, {4, -3, -3}}]");
    chk_if("s . r", "SparseArray[Automatic, {3, 3}, 0.0, {1, {{0, 1, 2, 2}, {{1}, {3}}}, {3.0, 6.0}}]");
    chk_if("r . s", "SparseArray[Automatic, {3, 3}, 0.0, {1, {{0, 1, 1, 2}, {{1}, {2}}}, {3.0, -2.0}}]");
    chk_if("s . sv", "SparseArray[Automatic, {3}, 0, {1, {{0, 2}, {{1}, {2}}}, {4, 15}}]");
    chk_if("sv . s", "SparseArray[Automatic, {3}, 0, {1, {{0, 2}, {{1}, {2}}}, {4, -5}}]");
    chk("sv . sv", "29");
    chk("sv . {1, 2, 3}", "17");
    chk("{1, 2, 3} . sv", "17");
    chk("SparseArray[{1->1.}, 3] . SparseArray[{2->1}, 3]", "0.0");
    chk("SparseArray[{1->1}, 3] . SparseArray[{2->1}, 3]", "0");
    /* entries are structural: a product that cancels stays explicit */
    chk_if("SparseArray[{{1,1}->1,{1,2}->1},{2,2}] . SparseArray[{{1,1}->1,{2,1}->-1},{2,2}]",
           "SparseArray[Automatic, {2, 2}, 0, {1, {{0, 1, 1}, {{1}}}, {0}}]");
    /* no explicit product: exact zero background */
    chk_if("SparseArray[{{1,2}->3},{2,2}] . SparseArray[{{1,1}->2.},{2,2}]",
           "SparseArray[Automatic, {2, 2}, 0, {1, {{0, 0, 0}, {}}, {}}]");
    /* columns come out sorted within each row */
    chk_if("SparseArray[{{1,1}->1,{1,2}->1},{2,2}] . SparseArray[{{1,2}->1,{2,1}->2},{2,2}]",
           "SparseArray[Automatic, {2, 2}, 0, {1, {{0, 2, 2}, {{1}, {2}}}, {2, 1}}]");
    /* complex: the zero background is the Real 0. */
    chk_if("c . s", "SparseArray[Automatic, {3, 3}, 0.0, {1, {{0, 1, 1, 1}, {{3}}}, {3.0 + 6.0*I}}]");
    chk_if("c . c", "SparseArray[Automatic, {3, 3}, 0, {1, {{0, 0, 0, 0}, {}}, {}}]");
    /* rationals stay exact */
    chk_if("SparseArray[{{1,1}->1/2},{2,2}] . SparseArray[{{1,2}->2/3},{2,2}]",
           "SparseArray[Automatic, {2, 2}, 0, {1, {{0, 1, 1}, {{2}}}, {1/3}}]");
}

/* an operand without explicit entries: an empty SparseArray, whatever the
 * other operand; 0. when either operand carries machine numbers */
static void test_empty(void) {
    chk_if("z . z", "SparseArray[Automatic, {3, 3}, 0, {1, {{0, 0, 0, 0}, {}}, {}}]");
    chk_if("z . {1, 2, 3}", "SparseArray[Automatic, {3}, 0, {1, {{0, 0}, {}}, {}}]");
    chk_if("z . {1., 2., 3.}", "SparseArray[Automatic, {3}, 0.0, {1, {{0, 0}, {}}, {}}]");
    chk_if("{1., 2., 3.} . z", "SparseArray[Automatic, {3}, 0.0, {1, {{0, 0}, {}}, {}}]");
    chk_if("z . r", "SparseArray[Automatic, {3, 3}, 0.0, {1, {{0, 0, 0, 0}, {}}, {}}]");
    chk_if("SparseArray[{}, {3, 3}, 0.] . s", "SparseArray[Automatic, {3, 3}, 0.0, {1, {{0, 0, 0, 0}, {}}, {}}]");
    chk_if("SparseArray[{{1,1}->x},{3,3}] . z", "SparseArray[Automatic, {3, 3}, 0, {1, {{0, 0, 0, 0}, {}}, {}}]");
    chk("SparseArray[{}, {2}] . {1, 2}", "0");
    chk("SparseArray[{}, {2}, 0.] . SparseArray[{}, {2}]", "0.0");
}

/* exact and symbolic values, and int64 overflow */
static void test_boxed(void) {
    chk("SparseArray[{{1,1}->2},{2,2}] . {x, y}", "{2 x, 0}");
    chk("SparseArray[{{1,1}->x},{2,2}] . {1, 2}", "{x, 0}");
    chk_if("s . SparseArray[{{1,1}->x}, {3,3}]",
           "SparseArray[Automatic, {3, 3}, 0, {1, {{0, 1, 1, 1}, {{1}}}, {2*x}}]");
    chk_if("SparseArray[{{1,1}->x},{2,2}] . SparseArray[{{1,1}->y},{2,2}]",
           "SparseArray[Automatic, {2, 2}, 0, {1, {{0, 1, 1}, {{1}}}, {x*y}}]");
    chk("SparseArray[{1->1+I},2] . {1, 0}", "1 + I");
    chk("SparseArray[{{1,1}->2^62, {1,2}->2^62},{2,2}] . {1, 1}", "{9223372036854775808, 0}");
    chk_if("SparseArray[{{1,1}->2^62, {1,2}->2^62},{2,2}] . SparseArray[{{1,1}->2^62, {1,2}->2^62},{2,2}]",
           "SparseArray[Automatic, {2, 2}, 0, {1, {{0, 2, 2}, {{1}, {2}}}, "
           "{21267647932558653966460912964485513216, 21267647932558653966460912964485513216}}]");
    chk("SparseArray[{{1,1}->2^70},{2,2}] . {1., 2.} == {2.^70, 0}", "True");
}

/* mismatched shapes: Dot::dotsh with the arrays as written, call unevaluated */
static void test_shapes(void) {
    chk("Quiet[Check[s . {1, 2}, $Failed]]", "$Failed");
    chk("Head[Quiet[s . {1, 2}]]", "Dot");
    /* a large mismatched product is reported, never densified */
    chk("big = SparseArray[{{1,1}->1., {100000,100000}->2.}, {100000, 100000}]; Null", "Null");
    chk("Quiet[Check[big . {1, 2}, $Failed]]", "$Failed");
    chk("Quiet[Check[big . big . s, $Failed]]", "$Failed");
    /* a nonzero background: the answer is dense */
    chk("SparseArray[{{1,1}->1}, {3,3}, 1] . s", "{{2, -1, 3}, {2, -1, 3}, {2, -1, 3}}");
    chk("SparseArray[{{1,1}->1}, {3,3}, 1] . {1, 2, 3}", "{6, 6, 6}");
    chk("Dot[s]", "SparseArray[<3>, {3, 3}]");
}

/* chains */
static void test_chains(void) {
    chk("s . s . {1, 2, 3}", "{4, -6, -9}");
    chk("Dot[{1, 2, 3}, s, s]", "{4, -6, -9}");
    chk_if("s . s . s", "SparseArray[Automatic, {3, 3}, 0, {1, {{0, 1, 2, 3}, {{1}, {3}, {2}}}, {8, -9, 3}}]");
    chk("Dot[{1, 2, 3}, s, {1, 1, 1}]", "5");
    chk("Dot[s, {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, s] == Normal[s . s]", "True");
}

static void test_matrixpower(void) {
    chk_if("MatrixPower[s, 3]", "SparseArray[Automatic, {3, 3}, 0, {1, {{0, 1, 2, 3}, {{1}, {3}, {2}}}, {8, -9, 3}}]");
    chk_if("MatrixPower[s, 1]", "SparseArray[Automatic, {3, 3}, 0, {1, {{0, 1, 2, 3}, {{1}, {3}, {2}}}, {2, 3, -1}}]");
    chk_if("MatrixPower[r, 2]", "SparseArray[Automatic, {3, 3}, 0.0, {1, {{0, 1, 1, 2}, {{1}, {3}}}, {2.25, 4.0}}]");
    chk_if("MatrixPower[SparseArray[{{1,2}->1},{2,2}], 2]",
           "SparseArray[Automatic, {2, 2}, 0, {1, {{0, 0, 0}, {}}, {}}]");
    chk_true("Normal[MatrixPower[s, 7]] === MatrixPower[Normal[s], 7]");
    chk("MatrixPower[s, 2, {1, 2, 3}]", "{4, -6, -9}");
    chk("MatrixPower[s, 2, SparseArray[{1->1}, 3]]", "{4, 0, 0}");
    chk("MatrixPower[s, 0, {1, 2, 3}]", "{1, 2, 3}");
    /* nonzero background and negative powers: dense */
    chk("MatrixPower[SparseArray[{{1,2}->1},{2,2},1], 2]", "{{2, 2}, {2, 2}}");
    chk("MatrixPower[s, -1]", "{{1/2, 0, 0}, {0, 0, -1}, {0, 1/3, 0}}");
    chk("Quiet[Check[MatrixPower[SparseArray[{1->1},3], 2], $Failed]]", "$Failed");
    chk("Quiet[Check[MatrixPower[s, 2, {1, 2}], $Failed]]", "$Failed");
}

static void test_inner(void) {
    chk("Inner[Times, s, {1, 2, 3}, Plus]", "{2, 9, -2}");
    chk("Inner[Times, {1, 2, 3}, s, Plus]", "{2, -3, 6}");
    chk("Inner[Times, s, s, Plus]", "{{4, 0, 0}, {0, -3, 0}, {0, 0, -3}}");
    chk("Inner[Times, SparseArray[{1->2},3], {1, 2, 3}, Plus]", "2");
}

/* sizes that take the threaded kernels, against the dense product */
static void test_large(void) {
    /* 0. backgrounds, so the dense references are packed */
    chk("SeedRandom[7]; n = 2000; "
        "A = SparseArray[Thread[RandomInteger[{1, n}, {120000, 2}] -> RandomReal[{-1, 1}, 120000]], {n, n}, 0.]; "
        "B = SparseArray[Thread[RandomInteger[{1, n}, {120000, 2}] -> RandomReal[{-1, 1}, 120000]], {n, n}, 0.]; "
        "I1 = SparseArray[Thread[RandomInteger[{1, n}, {120000, 2}] -> RandomInteger[{-9, 9}, 120000]], {n, n}]; "
        "x = RandomReal[{-1, 1}, n]; M = RandomReal[{-1, 1}, {n, 4}]; W = RandomReal[{-1, 1}, {4, n}]; "
        "dA = Normal[A]; dB = Normal[B]; Null", "Null");
    chk_true("Max[Abs[A . x - dA . x]] < 10.^-10");
    chk_true("Max[Abs[x . A - x . dA]] < 10.^-10");
    chk_true("Max[Abs[A . M - dA . M]] < 10.^-10");
    chk_true("Max[Abs[W . A - W . dA]] < 10.^-10");
    chk_true("Head[A . B] === SparseArray && Max[Abs[Normal[A . B] - dA . dB]] < 10.^-10");
    chk_true("Max[Abs[Normal[A . Transpose[B]] - dA . Transpose[dB]]] < 10.^-10");
    /* exact integers through the threaded SpGEMM: the same numbers as the
     * (exactly representable) machine product */
    chk_true("Max[Abs[Normal[I1 . I1] - Normal[N[I1]] . Normal[N[I1]]]] == 0.");
    chk_true("Head[First[ArrayRules[I1 . I1]][[2]]] === Integer");
    chk_true("Max[Abs[MatrixPower[A, 3] . x - dA . (dA . (dA . x))]] < 10.^-10");
    chk_true("Max[Abs[MatrixPower[A, 3, x] - dA . (dA . (dA . x))]] < 10.^-10");
    chk_true("Max[Abs[A . A . x - dA . (dA . x)]] < 10.^-10");
}

int main(void) {
    symtab_init();
    core_init();
    setup();
    test_dense_answers();
    test_sparse_answers();
    test_empty();
    test_boxed();
    test_shapes();
    test_chains();
    test_matrixpower();
    test_inner();
    test_large();
    if (failures) {
        fprintf(stderr, "\n%d sparse Dot test(s) FAILED\n", failures);
        return 1;
    }
    printf("All sparse Dot tests passed.\n");
    return 0;
}
