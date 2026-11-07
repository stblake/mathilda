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

/* ---- the band solver and the refined Krylov fallback ---------------------- */

/* Without SuiteSparse the sparse direct solver is the RCM-ordered LAPACK band
 * factor (sparse_la_band.c); with it the same systems go through UMFPACK /
 * CHOLMOD. Either way: the accuracy of a direct solve, and fast. True values
 * from the SIAM 100-digit holdout (benchmarks/202, truth.tsv: 75 digits). */
static void test_band(void) {
    chk("pm[n_, offs_, d_] := SparseArray[Join[{Band[{1, 1}] -> N[d]}, "
        "Flatten[Table[{Band[{1, 1 + o}] -> 1., Band[{1 + o, 1}] -> 1.}, {o, offs}]]], {n, n}]; Null", "Null");
    /* tridiagonal 20000 x 20000 (beyond the densify cap: was Krylov) */
    chk_true("t = AbsoluteTiming[x = First[LinearSolve[pm[20000, {1}, Prime[Range[20000]]], "
             "UnitVector[20000, 1]]]]; Abs[x - 0.608978172457445740351] < 3*^-16 && t[[1]] < 1");
    /* pentadiagonal 10000 x 10000 (was a 5 s dense LU), indefinite near the top */
    chk_true("t = AbsoluteTiming[x = First[LinearSolve[pm[10000, {1, 2}, Prime[Range[10000]]], "
             "UnitVector[10000, 1]]]]; Abs[x - 0.648357732531959834025] < 3*^-16 && t[[1]] < 1");
    /* the 5-point Laplacian on a 60 x 60 grid (band Cholesky, bandwidth 60) */
    chk("lap2[m_] := With[{n = m^2, r = Table[If[Mod[i, m] == 0, 0., -1.], {i, m^2 - 1}]}, "
        "SparseArray[{Band[{1, 1}] -> 4., Band[{1, 2}] -> r, Band[{2, 1}] -> r, "
        "Band[{1, 1 + m}] -> -1., Band[{1 + m, 1}] -> -1.}, {n, n}]]; Null", "Null");
    chk_true("x = First[LinearSolve[lap2[60], UnitVector[3600, 1]]]; "
             "Abs[x - 0.302347219326549719479] < 1*^-15");
    chk_true("With[{a = lap2[40]}, Abs[Det[a]/Det[Normal[a]] - 1] < 10^-9]");
    /* a banded matrix hidden by a symmetric permutation: RCM recovers the band */
    chk_true("SeedRandom[3]; pp = RandomSample[Range[2000]]; "
             "a = pm[2000, {1, 2, 3}, Range[2000] + 10.][[pp, pp]]; b = N[Range[2000]]; "
             "res[a, LinearSolve[a, b], b] < 10^-14");
    /* symmetric with a positive diagonal but indefinite (eigenvalues
     * 1 + 2 cos(k pi / 401)): the Cholesky fails,
     * the band LU answers */
    chk_true("With[{a = pm[400, {1}, ConstantArray[1, 400]], b = N[Range[400]]}, "
             "res[a, LinearSolve[a, b], b] < 10^-10 && "
             "Max[Abs[LinearSolve[a, b] - LinearSolve[Normal[a], b]]] < 10^-8]");
    chk_msg("LinearSolve[pm[400, {1}, ConstantArray[1, 400]], N[Range[400]], Method -> \"Cholesky\"]",
            "LinearSolve::npdef");
    /* complex: Hermitian (band Cholesky) and general (band LU); transposed
     * solves through a LinearSolveFunction */
    chk_true("With[{a = SparseArray[{Band[{1,1}] -> 4., Band[{1,2}] -> I, Band[{2,1}] -> -I}, {400, 400}], "
             "b = N[Range[400]] + I}, res[a, LinearSolve[a, b], b] < 10^-13]");
    chk_true("With[{a = SparseArray[{Band[{1,1}] -> 4. + I, Band[{1,2}] -> 2., Band[{3,1}] -> -1. + 2 I}, {400, 400}], "
             "b = N[Range[400]]}, res[a, LinearSolve[a, b], b] < 10^-13 && "
             "Abs[Det[a]/Det[Normal[a]] - 1] < 10^-9]");
    chk_true("With[{a = nonsym[400], b = N[Range[400]]}, "
             "Max[Abs[Inverse[a] . b - LinearSolve[a, b]]] < 10^-12]");
    /* transposed solves (LinearSolveFunction "T") through the narrow in-house
     * band LU (pivoting: a weak diagonal) and the LAPACK band LU (wide) */
    chk_true("With[{a = SparseArray[{Band[{1,1}] -> 0.1, Band[{1,2}] -> 3., Band[{2,1}] -> -2., "
             "Band[{1,3}] -> 1.}, {300, 300}], b = N[Range[300]]}, f = LinearSolve[a]; "
             "Max[Abs[Transpose[a] . f[b, \"T\"] - b]] < 10^-9 && res[a, f[b], b] < 10^-12 && "
             "Abs[Det[a]/Det[Normal[a]] - 1] < 10^-9]");
    chk_true("With[{a = SparseArray[{Band[{1,1}] -> 1.5, Band[{1,2}] -> 3., Band[{12,1}] -> -2., "
             "Band[{1,14}] -> 1.}, {300, 300}], b = N[Range[300]]}, f = LinearSolve[a]; "
             "Max[Abs[Transpose[a] . f[b, \"T\"] - b]] < 10^-9 && res[a, f[b], b] < 10^-12 && "
             "Abs[Det[a]/Det[Normal[a]] - 1] < 10^-9]");
    chk_true("With[{a = lap2[30], b = N[Range[900]]}, f = LinearSolve[a]; "
             "Max[Abs[f[b, \"T\"] - f[b]]] < 10^-12]");
    /* the 2^k-offset matrix of SIAM problem 7 (no band: Krylov), refined to
     * the accuracy of a direct solve -- 12 correct digits before */
    chk_true("x = First[LinearSolve[pm[20000, 2^Range[0, 14], Prime[Range[20000]]], "
             "UnitVector[20000, 1]]]; Abs[x - 0.725078346268401167469] < 5*^-16");
    chk("Clear[pm, lap2, t, x, a, b, pp, f]; Null", "Null");
}

/* ---- arbitrary precision: mixed-precision iterative refinement ------------- */

/* sparse_la_mpsolve.c: a machine factor of the SparseArray, residuals in
 * MPFR. Each answer is checked against an exact rational solve or a
 * 75-digit truth (benchmarks/202 truth.tsv). */
static void test_mpsolve(void) {
    chk("pe[n_, offs_, d_] := SparseArray[Join[{Band[{1, 1}] -> d}, "
        "Flatten[Table[{Band[{1, 1 + o}] -> 1, Band[{1 + o, 1}] -> 1}, {o, offs}]]], {n, n}]; Null", "Null");
    /* 20000 x 20000 tridiagonal at 60 digits: was SparseArray::dense, unevaluated */
    chk_true("x = First[LinearSolve[N[pe[20000, {1}, Prime[Range[20000]]], 60], UnitVector[20000, 1]]]; "
             "Precision[x] >= 59 && Abs[x - 0.608978172457445740351364767768399653994707493429771615746667998835697`75] < 10^-57");
    /* exact references: pentadiagonal, an indefinite 2^k-offset matrix, a
     * 2-D Laplacian; MPFR, exact and SparseArray right-hand sides */
    chk_true("a = pe[300, {1, 2}, Prime[Range[300]]]; b = Range[300]; xe = LinearSolve[a, b]; "
             "x = LinearSolve[N[a, 50], b]; Max[Abs[x - xe]] < 10^-45 && Min[Precision /@ x] >= 49");
    chk_true("a = pe[400, 2^Range[0, 8], Range[400]]; b = UnitVector[400, 1]; xe = LinearSolve[a, b]; "
             "x = LinearSolve[N[a, 40], N[b, 40]]; Max[Abs[x - xe]] < 10^-34");
    chk_true("a = pe[400, 2^Range[0, 8], Range[400]]; xe = LinearSolve[a, UnitVector[400, 7]]; "
             "x = LinearSolve[N[a, 40], SparseArray[{7 -> 1}, 400]]; Max[Abs[x - xe]] < 10^-34");
    chk_true("a = SparseArray[{Band[{1,1}] -> 4, Band[{1,2}] -> -1, Band[{2,1}] -> -1, "
             "Band[{1,21}] -> -1, Band[{21,1}] -> -1}, {400, 400}]; xe = LinearSolve[a, Range[400]]; "
             "x = LinearSolve[N[a, 100], Range[400]]; Max[Abs[x - xe]/Abs[xe]] < 10^-95");
    /* a nonsymmetric real matrix with fractional 30-digit entries */
    chk_true("a = SparseArray[{Band[{1,1}] -> 3, Band[{1,2}] -> 1/3, Band[{3,1}] -> -2/7}, {300, 300}]; "
             "xe = LinearSolve[a, Range[300]]; x = LinearSolve[N[a, 30], N[Range[300], 30]]; "
             "Max[Abs[x - xe]] < 10^-25");
    /* declines -- the dense arbitrary-precision path answers as before:
     * singular, too ill-conditioned for machine corrections, complex */
    chk_true("a = N[SparseArray[{{1,1} -> 1, {2,2} -> 1, {1,2} -> 1, {2,1} -> 1, {3,3} -> 2}, {3, 3}], 30]; "
             "Quiet[LinearSolve[a, N[{1, 1, 1}, 30]] == LinearSolve[Normal[a], N[{1, 1, 1}, 30]]]");
    chk_true("a = SparseArray[HilbertMatrix[16]]; xe = LinearSolve[HilbertMatrix[16], Range[16]]; "
             "x = Quiet[LinearSolve[N[a, 60], Range[16]]]; Max[Abs[x - xe]/Max[Abs[xe]]] < 10^-30");
    chk_true("a = N[SparseArray[{Band[{1,1}] -> 4 + I, Band[{1,2}] -> 1}, {200, 200}], 30]; "
             "x = LinearSolve[a, N[Range[200], 30]]; Max[Abs[a . x - Range[200]]] < 10^-25");
    chk("Clear[pe, a, b, x, xe]; Null", "Null");
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

/* ---- exact rational sparse elimination (sparse_la_exact.c) ------------------ */

static void test_exact(void) {
    chk("ea = SparseArray[{Band[{1,1}] -> 3, Band[{1,2}] -> -1, Band[{2,1}] -> 2, Band[{1, 9}] -> 1, "
        "Band[{5, 1}] -> -2}, {100, 100}]; "
        "eq = SparseArray[{Band[{1,1}] -> 1/3, Band[{1,2}] -> -2/7, Band[{2,1}] -> 5, Band[{1,20}] -> 1/2}, {100, 100}]; "
        "es = SparseArray[{Band[{1,1}] -> 1, Band[{1,2}] -> 1}, {100, 100}] + SparseArray[{{100, 100} -> -1}, {100, 100}]; Null",
        "Null");
    /* against the dense exact path: solve, Det, MatrixRank, rational entries, matrix rhs */
    chk_true("LinearSolve[ea, Range[100]] === LinearSolve[Normal[ea], Range[100]]");
    chk_true("Det[ea] === Det[Normal[ea]] && MatrixRank[ea] === MatrixRank[Normal[ea]]");
    chk_true("LinearSolve[eq, Range[100]] === LinearSolve[Normal[eq], Range[100]] && Det[eq] === Det[Normal[eq]]");
    chk_true("LinearSolve[ea, Transpose[{Range[100], 2 Range[100]}]] === "
             "LinearSolve[Normal[ea], Transpose[{Range[100], 2 Range[100]}]]");
    chk_true("With[{p = SparseArray[Table[{i, Mod[3 i, 100] + 1} -> i, {i, 100}], {100, 100}]}, "
             "Det[p] === Det[Normal[p]] && LinearSolve[p, Range[100]] === LinearSolve[Normal[p], Range[100]]]");
    /* fill with exact cancellation: a fill column must be indexed even when
     * the row comes out no longer than before (the 400 x 400 grid-like
     * matrix answered a wrong solution, A . x != b), and random sparsity */
    chk_true("With[{g = SparseArray[{Band[{1,1}] -> 4, Band[{1,2}] -> -1, Band[{2,1}] -> -1, "
             "Band[{1,21}] -> -1, Band[{21,1}] -> -1}, {400, 400}]}, "
             "x = LinearSolve[g, Range[400]]; g . x == Range[400] && x === LinearSolve[Normal[g], Range[400]] && "
             "Det[g] === Det[Normal[g]]]");
    chk_true("SeedRandom[11]; With[{r = SparseArray[Join[Table[{i, i} -> RandomInteger[{1, 9}], {i, 300}], "
             "Table[RandomInteger[{1, 300}, 2] -> RandomInteger[{-3, 3}], {900}]], {300, 300}]}, "
             "x = LinearSolve[r, Range[300]]; r . x == Range[300] && Det[r] === Det[Normal[r]] && "
             "MatrixRank[r] === MatrixRank[Normal[r]]]");
    /* singular (an explicit stored zero on the diagonal): rank, Det 0, a
     * particular solution of a consistent system, ::nosol otherwise */
    chk("{MatrixRank[es], Det[es]}", "{99, 0}");
    chk_true("With[{c = es . ConstantArray[1, 100]}, es . LinearSolve[es, c] == c]");
    chk_msg("LinearSolve[es, UnitVector[100, 100]]", "LinearSolve::nosol");
    /* past the densify cap */
    chk_true("With[{L = SparseArray[{Band[{1,1}] -> 2, Band[{1,2}] -> -1, Band[{2,1}] -> -1}, {20000, 20000}]}, "
             "x = LinearSolve[L, ConstantArray[1, 20000]]; x[[1;;3]] == {10000, 19999, 29997} && "
             "L . x == ConstantArray[1, 20000] && Det[L] == 20001 && MatrixRank[L] == 20000]");
}

/* ---- MatrixExp of a SparseArray (sparse_la_func.c) ------------------------- */

static void test_matrixexp(void) {
    /* a machine SparseArray answers a SparseArray (reference head); exact a List */
    chk("Head[MatrixExp[SparseArray[{{1,1}->1.,{2,2}->2.,{1,2}->0.5}, {2,2}]]]", "SparseArray");
    chk("Head[MatrixExp[SparseArray[{{1,1}->1,{2,2}->2}]]]", "List");
    chk("Head[MatrixExp[SparseArray[{{1,1}->1.,{2,2}->2.}, {2,2}, 0.5]]]", "List");
    /* an arbitrary-precision SparseArray: SparseArray (reference); the action
     * keeps the precision (dense path) */
    chk("Head[MatrixExp[N[SparseArray[{{1,1}->1,{2,2}->2}],30]]]", "SparseArray");
    /* Precision of a vector is one number, the least over its entries */
    chk("Precision[MatrixExp[N[SparseArray[{{1,1}->1,{2,2}->2}],30], {1,1}]]", "30.0");
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
    test_band();
    test_mpsolve();
    test_messages();
    test_methods();
    test_det_inverse();
    test_lsf();
    test_matrixexp();
    test_exact();
    if (failures) {
        fprintf(stderr, "%d failure(s)\n", failures);
        return 1;
    }
    printf("sparse_la_solve_tests: all passed\n");
    return 0;
}
