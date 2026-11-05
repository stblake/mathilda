/* test_sparse_la_solve.c -- native LinearSolve / Det / Inverse over
 * SparseArray, the Krylov solvers, and LinearSolveFunction
 * (src/sparse/sparse_la_*.c, src/linalg/linsolvefn.c).
 *
 * Expected values were taken from the reference system (wolframscript) where
 * an answer is unique; a numeric answer that depends on the backend
 * (UMFPACK, CHOLMOD, LAPACK, a Krylov iteration) is checked through its
 * residual instead. Covers: machine real and complex, exact, MPFR and
 * symbolic entries; dense, packed, NDArray and SparseArray right-hand sides;
 * the reference's sparse-only messages (sing, sqmat, npdef, kryme, krypre,
 * krymit, kryinfp, sing1, sqmat1, lslc, notrhsm); every Krylov method and
 * preconditioner; sizes above the densify cap.
 *
 * Failures are counted rather than asserted (Release builds pass -DNDEBUG). */

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

static void chk_true(const char* input) { chk(input, "True"); }

/* `input` must leave a message with this tag (Check returns "msg") */
static void chk_msg(const char* input, const char* tag) {
    char buf[4096];
    snprintf(buf, sizeof buf, "Quiet[Check[%s; \"none\", \"msg\", %s]]", input, tag);
    chk(buf, "\"msg\"");
}

static void setup(void) {
    chk("s = SparseArray[{{1,1}->4.,{2,2}->3.,{1,2}->1.}]; "
        "c = SparseArray[{{1,1}->2.+I,{2,2}->3.,{1,2}->1.}]; "
        "sing = SparseArray[{{1,1}->1.,{2,2}->1.}, {3,3}]; "
        "rect = SparseArray[{{1,1}->1.,{2,3}->2.}]; "
        "lap[n_] := SparseArray[{Band[{1,1}] -> 2., Band[{1,2}] -> -1., Band[{2,1}] -> -1.}, {n, n}]; "
        "nonsym[n_] := SparseArray[{Band[{1,1}] -> 4., Band[{1,2}] -> -1., Band[{2,1}] -> -2., "
        "  Band[{1,3}] -> 0.5}, {n, n}]; "
        "res[a_, x_, b_] := Max[Abs[Flatten[{a . x - b}]]]/Max[1, Max[Abs[Flatten[{b}]]]]; "
        "Null", "Null");
}

/* ---- direct machine solves --------------------------------------------- */

static void test_direct(void) {
    /* reference: {0.16666666666666669, 0.3333333333333333} */
    chk("LinearSolve[s, {1., 1.}]", "{0.166667, 0.333333}");
    chk("Head[LinearSolve[s, {1., 1.}]]", "List");
    chk("LinearSolve[s, {{1., 0.}, {0., 1.}}]", "{{0.25, -0.0833333}, {0.0, 0.333333}}");
    /* a sparse right-hand side: a dense answer, as in the reference */
    chk("Head[LinearSolve[s, SparseArray[{1 -> 1.}, 2]]]", "List");
    chk("LinearSolve[s, SparseArray[{1 -> 1.}, 2]]", "{0.25, 0.0}");
    chk("Head[LinearSolve[s, SparseArray[{{1,1} -> 1.}, {2, 2}]]]", "List");
    /* exact right-hand side, machine matrix */
    chk("LinearSolve[s, {1, 1}]", "{0.166667, 0.333333}");
    /* complex: agrees with the dense matrix */
    chk("LinearSolve[c, {1., I}]", "{0.333333 - 0.333333*I, 0.0 + 0.333333*I}");
    chk_true("Chop[LinearSolve[c, {1., I}] - LinearSolve[Normal[c], {1., I}]] == {0, 0}");
    /* a complex right-hand side on a real matrix */
    chk_true("Chop[LinearSolve[s, {1., I}] - LinearSolve[Normal[s], {1., I}]] == {0, 0}");
    /* NDArray right-hand side: an NDArray answer */
    chk("Head[LinearSolve[s, NDArray[{1., 1.}]]]", "NDArray");
    /* large: above the dense threshold, nonsymmetric and SPD */
    chk_true("With[{a = nonsym[3000], b = N[Range[3000]]}, res[a, LinearSolve[a, b], b] < 10^-10]");
    chk_true("With[{a = lap[3000], b = ConstantArray[1., 3000]}, res[a, LinearSolve[a, b], b] < 10^-6]");
    chk_true("With[{a = lap[300], b = N[Range[300]]}, "
             "res[a, LinearSolve[a, Transpose[{b, 2 b}]], Transpose[{b, 2 b}]] < 10^-8]");
    /* far above the densify cap (2^27 elements): sparse all the way */
    chk_true("With[{a = lap[20000], b = ConstantArray[1., 20000]}, "
             "x = LinearSolve[a, b]; Length[x] == 20000 && Max[Abs[a . x - b]] < 10^-4]");
    chk_true("With[{a = nonsym[20000], b = ConstantArray[1., 20000]}, "
             "x = LinearSolve[a, b]; Max[Abs[a . x - b]] < 10^-9]");
}

/* ---- reference semantics that differ from the dense matrix ---------------- */

static void test_messages(void) {
    /* singular machine sparse: ::sing, unevaluated (dense: particular solution) */
    chk_msg("LinearSolve[sing, {1., 1., 0.}]", "LinearSolve::sing");
    chk("Head[Quiet[LinearSolve[sing, {1., 1., 0.}]]]", "LinearSolve");
    chk("Quiet[LinearSolve[Normal[sing], {1., 1., 0.}]]", "{1.0, 1.0, 0.0}");
    /* rectangular machine sparse: ::sqmat */
    chk_msg("LinearSolve[rect, {1., 1.}]", "LinearSolve::sqmat");
    chk("Head[Quiet[LinearSolve[rect, {1., 1.}]]]", "LinearSolve");
    /* exact sparse behaves like the dense matrix (reference) */
    chk("LinearSolve[SparseArray[{{1,1}->1,{2,2}->1}, {3,3}], {1,1,0}]", "{1, 1, 0}");
    chk("LinearSolve[SparseArray[{{1,1}->1,{2,2}->2,{3,1}->1}, {3,2}], {1,2,1}]", "{1, 1}");
    chk_msg("LinearSolve[SparseArray[{{1,1}->1,{2,2}->1}, {3,3}], {1,1,1}]", "LinearSolve::nosol");
    chk("LinearSolve[SparseArray[{{1,1}->a,{2,2}->b}], {1,1}]", "{1/a, 1/b}");
    /* Cholesky on a matrix that is not positive definite */
    chk_msg("LinearSolve[SparseArray[{{1,1}->-4.,{2,2}->3.}], {1., 1.}, Method -> \"Cholesky\"]",
            "LinearSolve::npdef");
    chk("Head[Quiet[LinearSolve[SparseArray[{{1,1}->-4.,{2,2}->3.}], {1., 1.}, Method -> \"Cholesky\"]]]",
        "LinearSolve");
}

/* ---- methods ---------------------------------------------------------------- */

static void test_methods(void) {
    chk_true("With[{a = lap[50], b = N[Range[50]]}, res[a, LinearSolve[a, b, Method -> \"Multifrontal\"], b] < 10^-10]");
    chk_true("With[{a = lap[50], b = N[Range[50]]}, res[a, LinearSolve[a, b, Method -> \"Cholesky\"], b] < 10^-10]");
    chk_true("With[{a = lap[500], b = N[Range[500]]}, res[a, LinearSolve[a, b, Method -> \"Cholesky\"], b] < 10^-8]");
    chk_true("With[{a = SparseArray[{{1,1}->4.,{2,2}->3.,{1,2}->1.+I,{2,1}->1.-I}], b = {1., 2.}}, "
             "res[a, LinearSolve[a, b, Method -> \"Cholesky\"], b] < 10^-14]");
    /* every Krylov method and preconditioner */
    const char* meths[] = { "ConjugateGradient", "BiCGSTAB", "GMRES" };
    const char* precs[] = { "\"ILU0\"", "\"ILUT\"", "\"ILUTP\"", "None", "Automatic" };
    char buf[1024];
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 5; j++) {
            snprintf(buf, sizeof buf,
                     "With[{a = lap[200], b = N[Range[200]]}, res[a, LinearSolve[a, b, Method -> "
                     "{\"Krylov\", Method -> \"%s\", \"Preconditioner\" -> %s}], b] < 10^-8]",
                     meths[i], precs[j]);
            chk_true(buf);
        }
    /* nonsymmetric: BiCGSTAB and GMRES */
    chk_true("With[{a = nonsym[2000], b = N[Range[2000]]}, res[a, LinearSolve[a, b, Method -> {\"Krylov\", Method -> \"GMRES\"}], b] < 10^-9]");
    chk_true("With[{a = nonsym[2000], b = N[Range[2000]]}, res[a, LinearSolve[a, b, Method -> \"Krylov\"], b] < 10^-9]");
    /* complex Krylov */
    chk_true("With[{a = SparseArray[{Band[{1,1}] -> 4. + I, Band[{1,2}] -> -1., Band[{2,1}] -> -1. + 0.5 I}, {300, 300}], "
             "b = N[Range[300]] I}, res[a, LinearSolve[a, b, Method -> {\"Krylov\", Method -> \"GMRES\"}], b] < 10^-9 && "
             "res[a, LinearSolve[a, b, Method -> {\"Krylov\", Method -> \"BiCGSTAB\"}], b] < 10^-9]");
    /* a user preconditioner (Jacobi as a function) */
    chk_true("With[{a = lap[100], b = N[Range[100]]}, res[a, LinearSolve[a, b, Method -> "
             "{\"Krylov\", \"Preconditioner\" -> (#/2. &)}], b] < 10^-8]");
    /* Tolerance and a starting vector */
    chk_true("With[{a = lap[100], b = N[Range[100]]}, res[a, LinearSolve[a, b, Method -> "
             "{\"Krylov\", Tolerance -> 10^-4, \"StartingVector\" -> ConstantArray[1., 100]}], b] < 10^-3]");
    /* messages */
    chk_msg("LinearSolve[lap[20], N[Range[20]], Method -> {\"Krylov\", Method -> \"Bogus\"}]", "LinearSolve::kryme");
    chk_msg("LinearSolve[lap[20], N[Range[20]], Method -> {\"Krylov\", \"Preconditioner\" -> \"Bogus\"}]", "LinearSolve::krypre");
    chk_msg("LinearSolve[lap[200], N[Range[200]], Method -> {\"Krylov\", \"MaxIterations\" -> 3, "
            "Method -> \"ConjugateGradient\", \"Preconditioner\" -> None}]", "LinearSolve::krymit");
    chk("Length[Quiet[LinearSolve[lap[200], N[Range[200]], Method -> {\"Krylov\", \"MaxIterations\" -> 3, "
        "Method -> \"ConjugateGradient\", \"Preconditioner\" -> None}]]]", "200");
    /* a dense matrix takes the sparse methods too; exact input says kryinfp
     * and still answers exactly (reference) */
    chk_true("With[{a = Normal[lap[60]], b = N[Range[60]]}, res[a, LinearSolve[a, b, Method -> \"Krylov\"], b] < 10^-9]");
    chk_msg("LinearSolve[{{1, 2}, {3, 4}}, {1, 1}, Method -> \"Krylov\"]", "LinearSolve::kryinfp");
    chk("Quiet[LinearSolve[{{1, 2}, {3, 4}}, {1, 1}, Method -> \"Krylov\"]]", "{-1, 1}");
    /* exact input: "Multifrontal" refuses (::exopt1), "Cholesky" numericalises */
    chk_msg("LinearSolve[{{1, 2}, {3, 4}}, {1, 1}, Method -> \"Multifrontal\"]", "LinearSolve::exopt1");
    chk("Head[Quiet[LinearSolve[SparseArray[{{1,1}->4,{2,2}->3}], {1,1}, Method->\"Multifrontal\"]]]", "LinearSolve");
    chk("LinearSolve[SparseArray[{{1,1}->4,{2,2}->3}], {1,1}, Method->\"Cholesky\"]", "{0.25, 0.333333}");
    /* method preconditions: Cholesky needs a Hermitian matrix (::herm, also
     * for machine input -- never a one-triangle factor of a nonsymmetric
     * matrix), Krylov explicit numbers (::krynum) */
    chk_msg("LinearSolve[SparseArray[{{1,1}->4,{2,2}->3,{1,2}->aa}], {1,1}, Method->\"Cholesky\"]", "LinearSolve::herm");
    chk_msg("LinearSolve[{{4., 1.}, {0., 3.}}, {1.,1.}, Method->\"Cholesky\"]", "LinearSolve::herm");
    chk_msg("LinearSolve[SparseArray[{Band[{1,1}] -> 4., Band[{1,2}] -> 1.}, {300, 300}], N[Range[300]], Method->\"Cholesky\"]",
            "LinearSolve::herm");
    chk_msg("LinearSolve[SparseArray[{{1,1}->4,{2,2}->3,{1,2}->aa}], {1,1}, Method->\"Krylov\"]", "LinearSolve::krynum");
    chk("Head[Quiet[LinearSolve[{{4, aa}, {0, 3}}, {1,1}, Method->\"Krylov\"]]]", "LinearSolve");
    /* machine contagion: integer or empty sparse with a machine b */
    chk_msg("LinearSolve[SparseArray[{}, {2,2}], {1.,1.}]", "LinearSolve::sing");
    chk_msg("LinearSolve[SparseArray[{{1,1}->1,{2,3}->2}], {1.,1.}]", "LinearSolve::sqmat");
    chk("LinearSolve[SparseArray[{{1,1}->4,{2,2}->3,{1,2}->1}], {1.,1.}]", "{0.166667, 0.333333}");
    /* a SparseArray mixing machine reals and complexes (kept boxed) */
    chk_true("With[{c = SparseArray[{Band[{1,1}] -> 4. + I, Band[{1,2}] -> -1., Band[{2,1}] -> -1. + 0.5 I}, {3000, 3000}], "
             "b = N[Range[3000]] I}, res[c, LinearSolve[c, b], b] < 10^-10]");
    /* dense mixed real/complex and exact-zero matrices take the machine path */
    chk_true("With[{d = Table[If[i == j, 4. + I, If[Abs[i - j] == 1, -1., 0]], {i, 200}, {j, 200}]}, "
             "TimeConstrained[res[d, LinearSolve[d, N[Range[200]]], N[Range[200]]] < 10^-12, 5, False]]");
}

/* ---- Det / Inverse --------------------------------------------------------- */

static void test_det_inverse(void) {
    chk("Det[s]", "12.0");
    chk("Inverse[s]", "{{0.25, -0.0833333}, {0.0, 0.333333}}");
    chk("Head[Inverse[s]]", "List");
    /* above the dense threshold: the sparse factor; compare with the dense one */
    chk_true("With[{a = nonsym[300]}, Abs[Det[a]/Det[Normal[a]] - 1] < 10^-10]");
    chk_true("With[{a = lap[200]}, Abs[Det[a] - 201] < 10^-8]");
    chk_true("With[{a = nonsym[200]}, Max[Abs[Inverse[a] . Normal[a] - IdentityMatrix[200]]] < 10^-12]");
    chk("Head[Inverse[nonsym[200]]]", "List");
    chk_true("With[{a = SparseArray[{Band[{1,1}] -> 4. + I, Band[{1,2}] -> -1.}, {200, 200}]}, "
             "Abs[Det[a]/Det[Normal[a]] - 1] < 10^-10]");
    /* singular */
    chk("Det[SparseArray[{{1,1}->1.}, {200, 200}]]", "0.0");
    chk_msg("Inverse[SparseArray[{{1,1}->1.}, {200, 200}]]", "Inverse::sing");
    /* exact sparse goes through the exact dense path */
    chk("Det[SparseArray[{{1,1}->4,{2,2}->3,{1,2}->1}]]", "12");
    chk("Inverse[SparseArray[{{1,1}->4,{2,2}->3,{1,2}->1}]]", "{{1/4, -1/12}, {0, 1/3}}");
}

/* ---- LinearSolveFunction ---------------------------------------------------- */

static void test_lsf(void) {
    chk("f = LinearSolve[{{4., 1.}, {2., 3.}}]; Head[f]", "LinearSolveFunction");
    chk("Length[f]", "2");
    chk("f[{1., 2.}]", "{0.1, 0.6}");
    chk("f[{{1., 0.}, {0., 1.}}]", "{{0.3, -0.1}, {-0.2, 0.4}}");
    chk("f[{1., 2.}, \"T\"]", "{-0.1, 0.7}");
    chk("f[{1., 2.}, \"N\"]", "{0.1, 0.6}");
    chk("f[SparseArray[{1 -> 1.}, 2]]", "{0.3, -0.2}");
    chk("f[\"MatrixDimensions\"]", "{2, 2}");
    chk("f[\"ConditionNumber\"]", "3.0");
    chk("f[\"Properties\"]", "{\"ConditionNumber\", \"getL\", \"getPermutations\", \"getU\", \"MatrixDimensions\", \"Properties\"}");
    chk_msg("f[{1., 2., 3.}]", "LinearSolveFunction::lslc");
    chk_msg("f[{1., 2.}, \"X\"]", "LinearSolve::notrhsm");
    /* complex: the four modes */
    chk("g = LinearSolve[{{1. + I, 2.}, {0., 3. - I}}]; Null", "Null");
    chk_true("Chop[g[{1., 2.}] - LinearSolve[{{1. + I, 2.}, {0., 3. - I}}, {1., 2.}]] == {0, 0}");
    chk_true("Chop[g[{1., 2.}, \"T\"] - LinearSolve[Transpose[{{1. + I, 2.}, {0., 3. - I}}], {1., 2.}]] == {0, 0}");
    chk_true("Chop[g[{1., 2.}, \"C\"] - LinearSolve[ConjugateTranspose[{{1. + I, 2.}, {0., 3. - I}}], {1., 2.}]] == {0, 0}");
    chk_true("Chop[g[{1., 2.}, \"J\"] - LinearSolve[Conjugate[{{1. + I, 2.}, {0., 3. - I}}], {1., 2.}]] == {0, 0}");
    /* exact, singular, rectangular (reference: still a LinearSolveFunction) */
    chk("e = LinearSolve[{{1, 2}, {3, 4}}]; e[{1, 1}]", "{-1, 1}");
    chk_msg("LinearSolve[{{1., 2.}, {2., 4.}}]", "LinearSolve::sing1");
    chk("Head[Quiet[LinearSolve[{{1., 2.}, {2., 4.}}]]]", "LinearSolveFunction");
    chk_msg("LinearSolve[{{1., 2., 3.}, {4., 5., 6.}}]", "LinearSolve::sqmat1");
    /* exact / symbolic singular matrices say sing1 too; badly conditioned luc */
    chk_msg("LinearSolve[SparseArray[{}, {3, 3}]]", "LinearSolve::sing1");
    chk_msg("LinearSolve[{{a, b}, {2 a, 2 b}}]", "LinearSolve::sing1");
    chk_msg("LinearSolve[N[HilbertMatrix[14]]]", "LinearSolve::luc");
    chk("Quiet[LinearSolve[{{1, 2, 3}, {4, 5, 6}}]][{1, 2}]", "{-1/3, 2/3, 0}");
    /* sparse: reuse the factor many times */
    chk("h = LinearSolve[lap[3000]]; Head[h]", "LinearSolveFunction");
    chk_true("And @@ Table[res[lap[3000], h[N[Range[3000]] k], N[Range[3000]] k] < 10^-6, {k, 1, 5}]");
    chk_true("res[lap[3000], h[N[Range[3000]], \"T\"], N[Range[3000]]] < 10^-6");
    chk_true("Head[LinearSolve[lap[20000]]] === LinearSolveFunction");
    /* Krylov LinearSolveFunction solves afresh */
    chk_true("With[{k = LinearSolve[lap[100], Method -> \"Krylov\"]}, res[lap[100], k[N[Range[100]]], N[Range[100]]] < 10^-8]");
    /* many objects: the cache evicts and refactors transparently */
    chk_true("fs = Table[LinearSolve[lap[150] + k IdentityMatrix[150, SparseArray]], {k, 12}]; "
             "And @@ Table[res[lap[150] + k IdentityMatrix[150, SparseArray], fs[[k]][N[Range[150]]], N[Range[150]]] < 10^-10, {k, 12}]");
}

/* ---- MatrixExp of a SparseArray (sparse_la_func.c) ------------------------- */

static void test_matrixexp(void) {
    /* a machine SparseArray answers a SparseArray (reference head); exact a List */
    chk("Head[MatrixExp[SparseArray[{{1,1}->1.,{2,2}->2.,{1,2}->0.5}, {2,2}]]]", "SparseArray");
    chk("Head[MatrixExp[SparseArray[{{1,1}->1,{2,2}->2}]]]", "List");
    /* an arbitrary-precision SparseArray: SparseArray (reference); the action
     * keeps the precision (dense path) */
    chk("Head[MatrixExp[N[SparseArray[{{1,1}->1,{2,2}->2}],30]]]", "SparseArray");
    chk("Precision[MatrixExp[N[SparseArray[{{1,1}->1,{2,2}->2}],30], {1,1}]]", "{30.103, 30.103}");
    chk_true("With[{s = SparseArray[{{1,1}->1.,{2,2}->2.,{1,2}->0.5}, {2,2}]}, "
             "Max[Abs[Normal[MatrixExp[s]] - MatrixExp[Normal[s]]]] < 10^-14]");
    /* the action on a vector, on SpMV: real, complex, against the dense action */
    chk("lap2[k_] := Module[{t = lap[k], id = IdentityMatrix[k, SparseArray]}, "
        "KroneckerProduct[t, id] + KroneckerProduct[id, t]]; Null", "Null");
    chk_true("With[{a = -0.1 lap2[30], v = N[Range[900]]/900}, "
             "Max[Abs[MatrixExp[a, v] - MatrixExp[Normal[a]] . v]] < 10^-13]");
    chk_true("With[{c = SparseArray[{Band[{1,1}] -> 1. + I, Band[{1,2}] -> -0.5}, {200, 200}], v = N[Range[200]]}, "
             "Max[Abs[MatrixExp[c, v] - MatrixExp[Normal[c]] . v]]/Max[Abs[MatrixExp[Normal[c]] . v]] < 10^-13]");
    chk("Head[MatrixExp[-0.1 lap2[30], N[Range[900]]]]", "List");
    /* far past the densify cap */
    chk_true("With[{a = -0.1 lap2[300], v = N[Range[90000]]/90000}, "
             "x = MatrixExp[a, v]; Length[x] == 90000 && Max[Abs[x]] <= 1.]");
}

int main(void) {
    symtab_init();
    core_init();
    test_mute_messages();
    setup();
    test_direct();
    test_messages();
    test_methods();
    test_det_inverse();
    test_lsf();
    test_matrixexp();
    if (failures) {
        fprintf(stderr, "%d failure(s)\n", failures);
        return 1;
    }
    printf("sparse_la_solve_tests: all passed\n");
    return 0;
}
