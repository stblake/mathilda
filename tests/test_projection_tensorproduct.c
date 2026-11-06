/* test_projection_tensorproduct.c -- Projection, TensorProduct (new heads,
 * src/linalg/projection.c, tensorproduct.c), the SparseArray matrix
 * predicates in O(nnz) (src/sparse/sparse_la_pred.c) and Cross of sparse
 * vectors. Expected values were taken from the reference system
 * (wolframscript). Failures are counted (Release builds pass -DNDEBUG). */

#include "test_utils.h"
#include "symtab.h"
#include "core.h"

static int failures = 0;

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
static void chk_if(const char* input, const char* expected) {
    char buf[4096];
    snprintf(buf, sizeof buf, "InputForm[%s]", input);
    chk(buf, expected);
}
static void chk_true(const char* input) { chk(input, "True"); }
static void chk_msg(const char* input, const char* tag) {
    char buf[4096];
    snprintf(buf, sizeof buf, "Quiet[Check[%s; \"none\", \"msg\", %s]]", input, tag);
    chk(buf, "\"msg\"");
}

static void test_projection(void) {
    chk("Projection[{1,2,3},{1,0,1}]", "{2, 0, 2}");
    chk_if("Projection[{1,I,3},{1,0,I}]", "{1/2 - (3*I)/2, 0, 3/2 + I/2}");
    chk_if("Projection[{a,b},{c,d}]",
           "{(c*(a*Conjugate[c] + b*Conjugate[d]))/(c*Conjugate[c] + d*Conjugate[d]), "
           "(d*(a*Conjugate[c] + b*Conjugate[d]))/(c*Conjugate[c] + d*Conjugate[d])}");
    chk("Projection[{1.,2.,3.},{1,0,1}]", "{2.0, 0.0, 2.0}");
    chk("Projection[{1,2},{3,4},Dot]", "{33/25, 44/25}");
    chk("Projection[{1,2},{3,4},#1 . #2 + 1 &]", "{18/13, 24/13}");
    chk("Projection[{1,2},{0,0}]", "{0, 0}");
    chk("Projection[{1.,2.},{0.,0.}]", "{0.0, 0.0}");
    chk_if("Projection[x, y]", "(y*Conjugate[y] . x)/Conjugate[y] . y");
    chk_msg("Projection[{1,2},{1,2,3}]", "Projection::prnv");
    chk("Head[Quiet[Projection[{1,2},{1,2,3}]]]", "Projection");
    /* Precision of a vector is one number, the least over its entries */
    chk("Precision[Projection[N[{1,2},30],{1,1}]]", "30.103");
    chk_true("Chop[Projection[{1.+I, 2.}, {1., I}] - Projection[{1+I, 2}, {1, I}]] == {0, 0}");
    /* sparse: a sparse v gives a sparse answer, a sparse u a dense one */
    chk("Head[Projection[SparseArray[{1->3,3->4},4], {1.,2.,3.,4.}]]", "List");
    chk("Projection[SparseArray[{1->3,3->4},4], {1.,2.,3.,4.}]", "{0.5, 1.0, 1.5, 2.0}");
    chk("Head[Projection[{1.,2.,3.,4.}, SparseArray[{1->3,3->4},4]]]", "SparseArray");
    /* every machine surface takes the one kernel: List, packed, NDArray, Compile */
    chk("Projection[NDArray[{1.,2.}], NDArray[{1.,1.}]]", "NDArray[{1.5, 1.5}]");
    chk("Lookup[CompileDiagnostics[{{u, _Real, 1}, {v, _Real, 1}}, Projection[u, v]], \"Compiled\"]", "True");
    chk("Lookup[CompileDiagnostics[{{u, _Complex, 1}, {v, _Complex, 1}}, Projection[u, v]], \"Compiled\"]", "True");
    chk_true("With[{u = RandomReal[{-1, 1}, 50], v = RandomReal[{-1, 1}, 50]}, "
             "Compile[{{a, _Real, 1}, {b, _Real, 1}}, Projection[a, b]][u, v] === Projection[u, v] && "
             "Normal[Projection[NDArray[u], NDArray[v]]] === Projection[u, v] && "
             "Projection[ToNDArray[u], ToNDArray[v]] === Projection[u, v]]");
}

static void test_tensorproduct(void) {
    chk("TensorProduct[{1,2},{3,4}]", "{{3, 4}, {6, 8}}");
    chk_if("TensorProduct[{{1,2},{3,4}},{a,b}]", "{{{a, b}, {2*a, 2*b}}, {{3*a, 3*b}, {4*a, 4*b}}}");
    chk("TensorProduct[2, {1,2}]", "{2, 4}");
    chk("TensorProduct[{1,2}]", "{1, 2}");
    chk("TensorProduct[]", "1");
    chk("TensorProduct[x, {1,2}]", "TensorProduct[x, {1, 2}]");
    chk("TensorProduct[a, TensorProduct[b, c]]", "TensorProduct[a, b, c]");
    chk("TensorProduct[{1.,2.},{3,4}]", "{{3.0, 4.0}, {6.0, 8.0}}");
    chk("Dimensions[TensorProduct[{1,2},{{1,2,3},{4,5,6}},{7,8}]]", "{2, 2, 3, 2}");
    chk_if("TensorProduct[SparseArray[{1->3,3->4},4], SparseArray[{2->1},3]]",
           "SparseArray[Automatic, {4, 3}, 0, {1, {{0, 1, 1, 2, 2}, {{2}, {2}}}, {3, 4}}]");
    chk_if("TensorProduct[SparseArray[{1->3,3->4},4], {1,2}]",
           "SparseArray[Automatic, {4, 2}, 0, {1, {{0, 2, 2, 4, 4}, {{1}, {2}, {1}, {2}}}, {3, 6, 4, 8}}]");
    chk("Lookup[CompileDiagnostics[{{u, _Real, 1}, {v, _Real, 2}}, TensorProduct[u, v]], \"Compiled\"]", "True");
    chk_true("With[{u = RandomReal[1, 4], v = RandomReal[1, {3, 2}]}, "
             "Compile[{{a, _Real, 1}, {b, _Real, 2}}, TensorProduct[a, b]][u, v] == TensorProduct[u, v]]");
}

static void test_predicates(void) {
    chk("s = SparseArray[{{1,1}->2.,{1,2}->1.,{2,1}->1.,{3,3}->4.}, {3,3}]; Null", "Null");
    chk("{SymmetricMatrixQ[s], HermitianMatrixQ[s], DiagonalMatrixQ[s], UpperTriangularMatrixQ[s], "
        "LowerTriangularMatrixQ[s], SquareMatrixQ[s]}", "{True, True, False, False, False, True}");
    chk("h = SparseArray[{{1,1}->2.,{1,2}->1.+I,{2,1}->1.-I}, {2,2}]; {SymmetricMatrixQ[h], HermitianMatrixQ[h]}",
        "{False, True}");
    chk("u = SparseArray[{{1,1}->1,{1,3}->2,{2,2}->3}, {3,3}]; {UpperTriangularMatrixQ[u], LowerTriangularMatrixQ[u], "
        "DiagonalMatrixQ[u]}", "{True, False, False}");
    chk("SquareMatrixQ[SparseArray[{{1,1}->1}, {2,3}]]", "False");
    chk("SymmetricMatrixQ[SparseArray[{{1,2}->a,{2,1}->a}, {2,2}]]", "True");
    chk("SymmetricMatrixQ[SparseArray[{{1,2}->a,{2,1}->b}, {2,2}]]", "False");
    /* an explicit stored zero (sparse arithmetic keeps it) is a zero entry */
    chk("z = SparseArray[{{1,1}->1,{1,2}->1,{2,2}->1}, {2,2}] - SparseArray[{{1,2}->1}, {2,2}]; "
        "{DiagonalMatrixQ[z], SymmetricMatrixQ[z], LowerTriangularMatrixQ[z]}", "{True, True, True}");
    /* far past the densify cap: answered from the pattern */
    chk("big = SparseArray[{Band[{1,1}] -> 2., Band[{1,2}] -> -1., Band[{2,1}] -> -1.}, {100000, 100000}]; "
        "{SymmetricMatrixQ[big], HermitianMatrixQ[big], DiagonalMatrixQ[big], UpperTriangularMatrixQ[big], "
        "SquareMatrixQ[big]}", "{True, True, False, False, True}");
    chk("PositiveDefiniteMatrixQ[big]", "True");
    chk("PositiveDefiniteMatrixQ[-big]", "False");
    chk_true("With[{m = SparseArray[{Band[{1,1}] -> 1., Band[{1,2}] -> -3.}, {500, 500}]}, "
             "PositiveDefiniteMatrixQ[m] === PositiveDefiniteMatrixQ[Normal[m]] === False]");
    chk_true("With[{m = SparseArray[{Band[{1,1}] -> 3., Band[{1,2}] -> -1.}, {500, 500}]}, "
             "PositiveDefiniteMatrixQ[m] === PositiveDefiniteMatrixQ[Normal[m]] === True]");
    /* agreement with the dense predicates */
    chk_true("With[{m = SparseArray[{{1,1}->2,{2,3}->1,{3,2}->1}, {4,4}]}, "
             "Table[p[m] === p[Normal[m]], {p, {SymmetricMatrixQ, HermitianMatrixQ, DiagonalMatrixQ, "
             "UpperTriangularMatrixQ, LowerTriangularMatrixQ, SquareMatrixQ}}] === ConstantArray[True, 6]]");
}

static void test_cross(void) {
    chk("Cross[SparseArray[{1->1},3], SparseArray[{2->1},3]]", "{0, 0, 1}");
    chk_true("Cross[SparseArray[{1->1.+I},3], {0., 1., 0.}] === Cross[{1.+I, 0, 0}, {0., 1., 0.}]");
    /* LowerTriangularMatrixQ (new dense builtin) -- reference values */
    chk("m = {{1,0,0},{2,3,0},{4,5,6}}; {LowerTriangularMatrixQ[m], LowerTriangularMatrixQ[m, 1], "
        "LowerTriangularMatrixQ[m, -1], LowerTriangularMatrixQ[{{1,2},{0,3}}], LowerTriangularMatrixQ[{{1,0,0},{2,3,0}}], "
        "LowerTriangularMatrixQ[{{1,10^-12},{1,1.}}, Tolerance -> 10^-10], LowerTriangularMatrixQ[x], "
        "LowerTriangularMatrixQ[{1,2}]}", "{True, True, False, False, True, True, False, False}");
    chk_msg("LowerTriangularMatrixQ[{{1,0},{2,3}}, 1, 2]", "LowerTriangularMatrixQ::nonopt");
    chk_msg("LowerTriangularMatrixQ[{{1,0},{2,3}}, x]", "LowerTriangularMatrixQ::nmint");
    chk_msg("UpperTriangularMatrixQ[{{1,0},{2,3}}, x]", "UpperTriangularMatrixQ::nmint");
}

int main(void) {
    symtab_init();
    core_init();
    test_mute_messages();
    test_projection();
    test_tensorproduct();
    test_predicates();
    test_cross();
    if (failures) { fprintf(stderr, "%d failure(s)\n", failures); return 1; }
    printf("projection_tensorproduct_tests: all passed\n");
    return 0;
}
