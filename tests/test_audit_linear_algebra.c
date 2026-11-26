/* Builtin-audit regressions: linear algebra (review/issues/<Head>.md;
 * tests/test_audit.h). */
#include "test_audit.h"

/* Fit#1 (critical, crash): "Fit with a NormFunction other than Norm[#,1]/
 * Norm[#,2] crashes".  fit_solve_findmin handed its coefficient vector to the
 * Dot (which adopts it) and freed it again.  The minimax line through
 * (0,1),(1,0),(3,2),(5,4) is 0.2 + 0.6 x: residuals -0.8, 0.8, 0, -0.8
 * equioscillate, so the max residual 0.8 is optimal. */
static void test_fit_general_norm(void) {
    audit_run("fd = {{0, 1}, {1, 0}, {3, 2}, {5, 4}}; ClearAll[fx]");
    audit_run("fm = Fit[fd, {1, fx}, fx, NormFunction -> Function[Norm[#, Infinity]]]");
    audit_true("Fit#1", "Abs[Coefficient[fm, fx, 0] - 0.2] < 10^-6 && Abs[Coefficient[fm, fx] - 0.6] < 10^-6");
    audit_true("Fit#1", "Abs[Max[Abs[(fm /. fx -> #[[1]]) - #[[2]] & /@ fd]] - 0.8] < 10^-6");
    /* no line on a fine grid of coefficients does better than 0.8 */
    audit_true("Fit#1",
        "Min[Table[Max[Abs[a + b fd[[All, 1]] - fd[[All, 2]]]], {a, -1, 1, 1/50}, {b, 0, 1, 1/50}]] >= 0.8 - 10^-9");
    audit_true("Fit#1", "MatchQ[Fit[fd, {1, fx}, fx, NormFunction -> (Norm[#, 3] &)], _Real + _Real fx]");
    audit_true("Fit#1", "MatchQ[Fit[fd, {1, fx}, fx, NormFunction -> (Total[Abs[#]^2] &)], _Real + _Real fx]");
    audit_true("Fit#1", "MatchQ[Fit[fd, {1, fx}, fx, NormFunction -> (Norm[#, 1] &), "
                        "FitRegularization -> {\"Tikhonov\", 1}], _Real + _Real fx]");
    /* sanity: the L2 fit is unchanged */
    audit_true("Fit#1", "Abs[Coefficient[Fit[fd, {1, fx}, fx], fx] - 41/59] < 10^-12");
}

/* Builtin-audit regressions: linear algebra (review/issues/; tests/test_audit.h). */
#include "test_audit.h"

/* An impossible size aborts with General::nomem: CheckAbort sees the abort and
 * Check sees the (routed) message. */
static void nomem_chk(const char* id, const char* expr) {
    char buf[512];
    snprintf(buf, sizeof buf, "Quiet[CheckAbort[%s, aborted]]", expr);
    audit_chk(id, buf, "aborted");
    snprintf(buf, sizeof buf, "Quiet[Check[CheckAbort[%s, 0], fired, General::nomem]]", expr);
    audit_chk(id, buf, "fired");
}

/* IdentityMatrix#1 (high): IdentityMatrix[{m, n}] with m*n overflowing size_t
 * built an undersized buffer (heap overflow, garbage entries, SIGBUS). */
static void test_identity_matrix_size_overflow(void) {
    nomem_chk("IdentityMatrix#1", "IdentityMatrix[{4611686018427387905, 4}][[100000000, 1]]");
    nomem_chk("IdentityMatrix#1", "IdentityMatrix[{2147483648, 1073741824}]");
    nomem_chk("IdentityMatrix#1", "IdentityMatrix[2^40]");
    audit_chk("IdentityMatrix#1", "IdentityMatrix[{2, 3}]", "{{1, 0, 0}, {0, 1, 0}}");
    audit_chk("IdentityMatrix#1", "Total[IdentityMatrix[300], 2]", "300");
}

/* DiagonalMatrix#1 (high): DiagonalMatrix with a huge {m, n} overflowed the
 * m*n buffer size and crashed. */
static void test_diagonal_matrix_size_overflow(void) {
    nomem_chk("DiagonalMatrix#1", "DiagonalMatrix[{1}, 0, {2147483648, 1073741824}]");
    nomem_chk("DiagonalMatrix#1", "DiagonalMatrix[{7}, 0, {3037000500, 3037000500}]");
    audit_chk("DiagonalMatrix#1", "DiagonalMatrix[{1}, 0, {3, 2}]", "{{1, 0}, {0, 0}, {0, 0}}");
    audit_chk("DiagonalMatrix#1", "DiagonalMatrix[{1, 2}, 1]", "{{0, 1, 0}, {0, 0, 2}, {0, 0, 0}}");
}

/* HankelMatrix#1 (high): HankelMatrix[2^31] crashed (n*n*8 wrapped to 0). */
static void test_hankel_matrix_size_overflow(void) {
    nomem_chk("HankelMatrix#1", "Dimensions[HankelMatrix[2147483648]]");
    nomem_chk("HankelMatrix#1", "HankelMatrix[2^62]");
    audit_chk("HankelMatrix#1", "HankelMatrix[3]", "{{1, 2, 3}, {2, 3, 0}, {3, 0, 0}}");
}

/* ToeplitzMatrix#1 (medium): ToeplitzMatrix[2^31] crashed. */
static void test_toeplitz_matrix_size_overflow(void) {
    nomem_chk("ToeplitzMatrix#1", "ToeplitzMatrix[2^31]");
    nomem_chk("ToeplitzMatrix#1", "ToeplitzMatrix[2^62]");
    audit_chk("ToeplitzMatrix#1", "ToeplitzMatrix[3]", "{{1, 2, 3}, {2, 1, 2}, {3, 2, 1}}");
}

/* VandermondeMatrix#1 (high): a huge column count k crashed (unchecked
 * malloc), on both the symbolic and the integer-node path. */
static void test_vandermonde_matrix_huge_k(void) {
    nomem_chk("VandermondeMatrix#1", "VandermondeMatrix[{x}, 10^15]");
    nomem_chk("VandermondeMatrix#1", "VandermondeMatrix[{1}, 10^15]");
    nomem_chk("VandermondeMatrix#1", "VandermondeMatrix[{1.5}, 2^62]");
    audit_chk("VandermondeMatrix#1", "VandermondeMatrix[{x, y}, 3]", "{{1, x, x^2}, {1, y, y^2}}");
}

/* PositiveDefiniteMatrixQ#1 (high): "Exact integer matrix with det -1 reported
 * positive definite (rounding to double)".  Exact input is now decided
 * exactly (Sylvester's criterion by fraction-free elimination); the same
 * root cause gave #2 (10^-20 off singular), #4 (10^-400) and #5 (10^400). */
static void test_posdef_exact(void) {
    audit_chk("PositiveDefiniteMatrixQ#1",
              "PositiveDefiniteMatrixQ[{{10^16, 10^16 + 1}, {10^16 + 1, 10^16 + 2}}]", "False");
    audit_true("PositiveDefiniteMatrixQ#1",
               "Det[{{10^16, 10^16 + 1}, {10^16 + 1, 10^16 + 2}}] == -1");
    audit_chk("PositiveDefiniteMatrixQ#1", "PositiveDefiniteMatrixQ[{{1, 1}, {1, 1 + 10^-20}}]", "True");
    audit_chk("PositiveDefiniteMatrixQ#1", "PositiveDefiniteMatrixQ[{{10^-400}}]", "True");
    audit_chk("PositiveDefiniteMatrixQ#1",
              "PositiveDefiniteMatrixQ[{{2*10^400, 10^400}, {10^400, 2*10^400}}]", "True");
    audit_chk("PositiveDefiniteMatrixQ#1", "PositiveDefiniteMatrixQ[HilbertMatrix[14]]", "True");
    audit_chk("PositiveDefiniteMatrixQ#1", "PositiveDefiniteMatrixQ[HilbertMatrix[40]]", "True");
    /* Hermitian part of an exact complex matrix; a non-symmetric real one */
    audit_chk("PositiveDefiniteMatrixQ#1", "PositiveDefiniteMatrixQ[{{2, I}, {-I, 2}}]", "True");
    audit_chk("PositiveDefiniteMatrixQ#1", "PositiveDefiniteMatrixQ[{{1, 2 I}, {-2 I, 1}}]", "False");
    audit_chk("PositiveDefiniteMatrixQ#1", "PositiveDefiniteMatrixQ[{{1, 5}, {-5, 1}}]", "True");
    audit_chk("PositiveDefiniteMatrixQ#1", "PositiveDefiniteMatrixQ[{{1, 2}, {2, 1}}]", "False");
    /* every leading minor of this 8 x 8 Gram matrix is checked, so a single
     * nonpositive one late in the elimination decides False */
    audit_chk("PositiveDefiniteMatrixQ#1",
              "PositiveDefiniteMatrixQ[Table[Min[i, j], {i, 8}, {j, 8}] - "
              "DiagonalMatrix[UnitVector[8, 8]]]", "False");
    audit_chk("PositiveDefiniteMatrixQ#1",
              "PositiveDefiniteMatrixQ[Table[Min[i, j], {i, 8}, {j, 8}]]", "True");
    /* sanity: machine input keeps its LAPACK path */
    audit_chk("PositiveDefiniteMatrixQ#1", "PositiveDefiniteMatrixQ[{{2., 1.}, {1., 2.}}]", "True");
    audit_chk("PositiveDefiniteMatrixQ#1", "PositiveDefiniteMatrixQ[{{1., 2.}, {2., 1.}}]", "False");
}

/* NegativeDefiniteMatrixQ#1 (high): "SparseArray input always gives False". */
static void test_negdef_sparse(void) {
    audit_chk("NegativeDefiniteMatrixQ#1",
              "NegativeDefiniteMatrixQ[SparseArray[{{-2, 0}, {0, -3}}]]", "True");
    audit_chk("NegativeDefiniteMatrixQ#1",
              "NegativeDefiniteMatrixQ[SparseArray[{{-2, 0}, {0, 3}}]]", "False");
    audit_true("NegativeDefiniteMatrixQ#1",
               "NegativeDefiniteMatrixQ[SparseArray[{{-2., 1.}, {1., -2.}}]] === "
               "NegativeDefiniteMatrixQ[{{-2., 1.}, {1., -2.}}] === True");
}

/* NegativeDefiniteMatrixQ#2 (high): "Exact integer matrix with det -1 reported
 * negative definite (rounding to double)". */
static void test_negdef_exact(void) {
    audit_chk("NegativeDefiniteMatrixQ#2",
              "NegativeDefiniteMatrixQ[{{-10^16, -10^16 - 1}, {-10^16 - 1, -10^16 - 2}}]", "False");
    audit_chk("NegativeDefiniteMatrixQ#2", "NegativeDefiniteMatrixQ[{{-1, -1}, {-1, -1 - 10^-20}}]", "True");
    audit_chk("NegativeDefiniteMatrixQ#2", "NegativeDefiniteMatrixQ[{{-10^-400}}]", "True");
    audit_chk("NegativeDefiniteMatrixQ#2", "NegativeDefiniteMatrixQ[-HilbertMatrix[14]]", "True");
    audit_chk("NegativeDefiniteMatrixQ#2", "NegativeDefiniteMatrixQ[HilbertMatrix[3]]", "False");
    audit_chk("NegativeDefiniteMatrixQ#2", "NegativeDefiniteMatrixQ[{{-2., 1.}, {1., -2.}}]", "True");
}

/* MatrixRank#1 (high): "MatrixRank undercounts rank for machine-real matrices
 * with large or tiny entries (naive complex division overflows/underflows)".
 * Rank is invariant under a nonzero scaling. */
static void test_matrixrank_scaled(void) {
    audit_chk("MatrixRank#1", "MatrixRank[{{1.*^-170, 0}, {0, 1.*^-170}}]", "2");
    audit_chk("MatrixRank#1", "MatrixRank[{{1.*^160, 2.*^160}, {3.*^160, 4.*^160}}]", "2");
    audit_chk("MatrixRank#1",
              "MatrixRank[{{1.*^-300, 2.*^-300}, {3.*^-300, 4.*^-300}}, Tolerance -> 0]", "2");
    audit_chk("MatrixRank#1", "MatrixRank[{{1.*^160 I, 2.*^160}, {3.*^160, 4.*^160 I}}]", "2");
    /* sanity */
    audit_chk("MatrixRank#1", "MatrixRank[{{1.*^160, 2.*^160}, {2.*^160, 4.*^160}}]", "1");
    audit_chk("MatrixRank#1", "MatrixRank[{{1., 2.}, {2., 4.}}]", "1");
}

/* Eigensystem#1 (high): "Arbitrary-precision Eigensystem of a defective matrix
 * returns vectors that are not eigenvectors (identity basis)".  Every returned
 * nonzero vector must satisfy m.v == lambda v; a defective eigenvalue's
 * missing vectors are zero, after the true ones. */
static void test_eigensystem_mpfr_defective(void) {
    audit_run("esm = N[{{1, 1}, {0, 1}}, 20]; ess = Eigensystem[esm]");
    audit_true("Eigensystem#1",
               "And @@ Table[Norm[esm . ess[[2, i]] - ess[[1, i]] ess[[2, i]]] < 10^-15, {i, 2}]");
    audit_true("Eigensystem#1", "Norm[ess[[2, 2]]] == 0 && Abs[ess[[2, 1, 2]]] < 10^-15");
    audit_run("esm3 = N[{{2, 1, 0}, {0, 2, 1}, {0, 0, 2}}, 20]; ess3 = Eigensystem[esm3]");
    audit_true("Eigensystem#1",
               "And @@ Table[Norm[esm3 . ess3[[2, i]] - ess3[[1, i]] ess3[[2, i]]] < 10^-15, {i, 3}]");
    audit_true("Eigensystem#1", "Count[Norm /@ ess3[[2]], z_ /; z == 0] == 2");
    audit_true("Eigensystem#1", "Count[Norm /@ Eigenvectors[N[{{2, 1}, {0, 2}}, 30]], z_ /; z == 0] == 1");
    /* sanity: a diagonalizable 30-digit matrix keeps all its vectors */
    audit_run("SeedRandom[3]; esr = N[RandomInteger[{-9, 9}, {6, 6}], 30]; esq = Eigensystem[esr]");
    audit_true("Eigensystem#1",
               "Max[Table[Norm[esr . esq[[2, i]] - esq[[1, i]] esq[[2, i]]], {i, 6}]] < 10^-20 && "
               "Count[Norm /@ esq[[2]], z_ /; z == 0] == 0");
    audit_run("esd = N[{{2, 0, 0}, {0, 2, 1}, {0, 0, 2}}, 20]; esdd = Eigensystem[esd]");
    audit_true("Eigensystem#1", "Count[Norm /@ esdd[[2]], z_ /; z == 0] == 1 && "
               "Max[Table[Norm[esd . esdd[[2, i]] - 2 esdd[[2, i]]], {i, 3}]] < 10^-15");
}

/* Eigenvalues#1 (high): "Exact singular pencil with shared null space returns
 * {Infinity, Infinity} instead of Indeterminate".  m = a = {{1,2},{2,4}}: the
 * shared null space {-2, 1} gives one Indeterminate, range(m) gives 1. */
static void test_eigenvalues_singular_pencil(void) {
    audit_chk("Eigenvalues#1", "Eigenvalues[{{{1, 2}, {2, 4}}, {{1, 2}, {2, 4}}}]", "{1, Indeterminate}");
    audit_chk("Eigenvalues#1",
              "Eigenvalues[{{{2, 0, 0}, {0, 3, 0}, {0, 0, 0}}, {{1, 0, 0}, {0, 1, 0}, {0, 0, 0}}}]",
              "{3, 2, Indeterminate}");
    /* sanity: a regular pencil with a singular a still has its Infinity */
    audit_chk("Eigenvalues#1", "Eigenvalues[{{{1, 2}, {3, 4}}, {{1, 0}, {0, 0}}}]", "{Infinity, -1/2}");
    audit_chk("Eigenvalues#1", "Eigenvalues[{{{2, 0}, {0, 3}}, {{1, 0}, {0, 1}}}]", "{3, 2}");
}

/* Fit#2 (high): "Machine Fit returns 0.0 when the response values are large
 * (>= ~1e13)".  The data lie exactly on 10^14 (1 + x). */
static void test_fit_large_response(void) {
    audit_run("ClearAll[fx]; fl = Fit[{{0, 10^14}, {1, 2 10^14}, {2, 3 10^14}}, {1, fx}, fx]");
    audit_true("Fit#2", "Abs[Coefficient[fl, fx, 0] - 10^14] < 10^2 && Abs[Coefficient[fl, fx] - 10^14] < 10^2");
    audit_run("fl2 = Fit[{10^20, 2 10^20, 3 10^20}, {1, fx}, fx]");
    audit_true("Fit#2", "Abs[(fl2 /. fx -> 2) - 2 10^20] < 10^7 && Abs[(fl2 /. fx -> 3) - 3 10^20] < 10^7");
    /* sanity: a basis function vanishing over the data still drops out */
    audit_chk("Fit#2", "Chop[Fit[N[{{-Pi,4},{-Pi/2,0},{0,1},{Pi/2,-1},{Pi,-4}}], "
                       "{Sin[fx/2], Sin[fx], Sin[2 fx]}, fx]]",
              "-4.0 Sin[1/2 fx] + 2.32843 Sin[fx]");
}

/* LeastSquares#1 / #2 (high): "Method -> Krylov/LSQR give a wrong
 * least-squares solution for exact complex matrices" and "return {0, 0} for
 * machine-precision complex matrices".  The minimiser solves
 * A^H A x == A^H b; x = {11/18, -1/3 + I/18}. */
static void test_leastsquares_complex_krylov(void) {
    audit_run("lsA = {{1, I}, {2, 3}, {0, 1}}; lsb = {1, 0, I}");
    audit_chk("LeastSquares#1", "LeastSquares[lsA, lsb, Method -> \"Krylov\"]", "{11/18, -1/3 + 1/18*I}");
    audit_chk("LeastSquares#1", "LeastSquares[lsA, lsb, Method -> \"LSQR\"]", "{11/18, -1/3 + 1/18*I}");
    audit_true("LeastSquares#1",
               "ConjugateTranspose[lsA] . (lsA . LeastSquares[lsA, lsb, Method -> \"Krylov\"] - lsb) == {0, 0}");
    audit_true("LeastSquares#2",
               "Norm[LeastSquares[N[lsA], lsb, Method -> \"Krylov\"] - {11/18, -1/3 + I/18}] < 10^-10");
    audit_true("LeastSquares#2",
               "Norm[LeastSquares[N[lsA], lsb, Method -> \"LSQR\"] - {11/18, -1/3 + I/18}] < 10^-10");
    /* sanity: real input unchanged */
    audit_chk("LeastSquares#1", "LeastSquares[{{1, 2}, {2, 3}, {0, 1}}, {1, 0, 1}, Method -> \"Krylov\"]",
              "{-5/3, 7/6}");
}

/* MatrixExp#1 (high): "MatrixExp of an exact diagonal matrix with an
 * irrational entry fails with Inverse::sing".  A triangular matrix's
 * eigenvalues are its diagonal; the quadratic formula gave nested radicals
 * the eigenvector step could not zero-test. */
static void test_matrixexp_irrational_diagonal(void) {
    audit_true("MatrixExp#1", "MatrixExp[{{Sqrt[2], 0}, {0, 1}}] === {{E^Sqrt[2], 0}, {0, E}}");
    audit_true("MatrixExp#1", "MatrixExp[{{1, 0}, {0, Log[2]}}] === {{E, 0}, {0, 2}}");
    audit_true("MatrixExp#1", "MatrixExp[{{Log[2], 0}, {0, Log[3]}}] === {{2, 0}, {0, 3}}");
    audit_true("MatrixExp#1", "Norm[N[MatrixExp[{{Sqrt[2], 1}, {0, 1}}]] - "
                              "MatrixExp[N[{{Sqrt[2], 1}, {0, 1}}]]] < 10^-12");
    audit_run("emj = JordanDecomposition[{{Sqrt[2], 0}, {0, 1}}]");
    audit_true("MatrixExp#1", "emj[[1]] . emj[[2]] . Inverse[emj[[1]]] == {{Sqrt[2], 0}, {0, 1}}");
    audit_run("ems = Eigensystem[{{Sqrt[2], 0}, {0, 1}}]");
    audit_true("MatrixExp#1", "And @@ Table[{{Sqrt[2], 0}, {0, 1}} . ems[[2, i]] == ems[[1, i]] ems[[2, i]], {i, 2}]");
    /* sanity */
    audit_true("MatrixExp#1", "MatrixExp[{{a, 0}, {0, b}}] === {{E^a, 0}, {0, E^b}}");
}

/* Inner#1 (high): "Incommensurate lengths are silently truncated or padded
 * with Null instead of Inner::incom". */
static void test_inner_incommensurate(void) {
    audit_chk("Inner#1", "Check[Inner[f, {a, b}, {x, y, z}, g], fired, Inner::incom]", "fired");
    audit_chk("Inner#1", "Quiet[Inner[f, {a, b}, {x, y, z}, g]]", "Inner[f, {a, b}, {x, y, z}, g]");
    audit_chk("Inner#1", "Check[Inner[f, {a, b}, {{x, y}}, g], fired, Inner::incom]", "fired");
    audit_chk("Inner#1", "Inner[f, {a, b}, {x, y}, g]", "g[f[a, x], f[b, y]]");
    audit_chk("Inner#1", "Inner[f, {a, b}, {{x, y}, {z, w}}, g]", "{g[f[a, x], f[b, z]], g[f[a, y], f[b, w]]}");
    audit_chk("Inner#1", "Inner[Times, {{1, 2}, {3, 4}}, {5, 6}, Plus]", "{17, 39}");
}

/* Inner#2 (high): "Inner[f, A, B, g, n] ignores n for any n other than 1". */
static void test_inner_index_n(void) {
    audit_chk("Inner#2", "Dimensions[Inner[f, Array[a, {2, 3, 4}], Array[b, {3, 5}], g, 2]]", "{2, 4, 5}");
    audit_true("Inner#2", "Inner[f, Array[a, {2, 3, 4}], Array[b, {3, 5}], g, 2] === "
                          "Inner[f, Transpose[Array[a, {2, 3, 4}], {1, 3, 2}], Array[b, {3, 5}], g]");
    audit_true("Inner#2", "Inner[f, Array[a, {2, 3}], Array[b, {2, 5}], g, 1] === "
                          "Inner[f, Transpose[Array[a, {2, 3}]], Array[b, {2, 5}], g]");
    audit_true("Inner#2", "Inner[Times, Array[a, {2, 3, 4}], Array[b, {4, 5}], Plus, 3] === "
                          "Array[a, {2, 3, 4}] . Array[b, {4, 5}]");
    audit_chk("Inner#2", "Quiet[Inner[f, {a, b}, {x, y}, g, 3]]", "Inner[f, {a, b}, {x, y}, g, 3]");
    audit_chk("Inner#2", "Check[Inner[f, Array[a, {2, 3, 4}], Array[b, {4, 5}], g, 2], fired, Inner::incom]", "fired");
}

/* NDArray#1 (high): "NDArray + Infinity discards the array and returns bare
 * Infinity". */
static void test_ndarray_infinity(void) {
    audit_chk("NDArray#1", "NDArray[{1., 2.}] + Infinity", "{Infinity, Infinity}");
    audit_chk("NDArray#1", "NDArray[{1., -2.}] Infinity", "{Infinity, -Infinity}");
    audit_chk("NDArray#1", "NDArray[{1., 2.}] - Infinity", "{-Infinity, -Infinity}");
    audit_true("NDArray#1", "(NDArray[{1., 2.}] + Infinity) === ({1., 2.} + Infinity)");
    /* sanity: a genuine symbol still leaves the sum unevaluated */
    audit_chk("NDArray#1", "Quiet[Head[NDArray[{1., 2.}] + q]]", "Plus");
}

/* NDArray#2 (high): "Elementwise Power/Times overflow and 0/0 put raw inf.0 /
 * nan.0 entries into the array".  0./0. is Indeterminate, 1/0. ComplexInfinity
 * (with their messages), 2.^1000000 a finite big Real. */
static void test_ndarray_nonfinite(void) {
    audit_chk("NDArray#2", "Quiet[NDArray[{0., 2.}]^-1]", "{ComplexInfinity, 0.5}");
    audit_chk("NDArray#2", "Quiet[NDArray[{0.}] / NDArray[{0.}]]", "{Indeterminate}");
    audit_chk("NDArray#2", "Check[NDArray[{0.}] / NDArray[{0.}], fired, Infinity::indet]", "fired");
    audit_chk("NDArray#2", "Quiet[Divide[0., 0.]]", "Indeterminate");
    audit_chk("NDArray#2", "Quiet[{0., 1.} / {0., 2.}]", "{Indeterminate, 0.5}");
    audit_true("NDArray#2", "Quiet[(NDArray[{2.}]^1000000)[[1]] == 2.^1000000]");
    audit_true("NDArray#2", "FreeQ[ToString[Quiet[NDArray[{0., 2.}]^-1]], \"inf\"]");
    /* sanity */
    audit_chk("NDArray#2", "NDArray[{1., 4.}] / NDArray[{2., 8.}]", "NDArray[{0.5, 0.5}]");
}

/* SingularValueDecomposition#1 (high): "Generalized SVD returns
 * upper-triangular, non-diagonal sigma matrices with negative entries".
 * w and wa are diagonal (pseudo-diagonal), non-negative, and reconstruct
 * m == u.w.ConjugateTranspose[v], a == ua.wa.ConjugateTranspose[v]. */
static void test_gsvd_diagonal(void) {
    audit_run("gsvdchk[m_, a_] := Module[{u, ua, w, wa, v}, "
              "{{u, ua}, {w, wa}, v} = SingularValueDecomposition[{m, a}]; "
              "Norm[u . w . ConjugateTranspose[v] - m] < 10^-12 && "
              "Norm[ua . wa . ConjugateTranspose[v] - a] < 10^-12 && "
              "Min[Flatten[{w, wa}]] >= 0 && "
              "Count[Flatten[w], x_ /; x != 0] <= Min[Dimensions[w]] && "
              "Count[Flatten[wa], x_ /; x != 0] <= Min[Dimensions[wa]] && "
              "And @@ Flatten[MapIndexed[#1 == 0 || #2[[1]] == #2[[2]] &, w, {2}]]]");
    audit_true("SingularValueDecomposition#1", "gsvdchk[{{1., 2.}, {3., 4.}}, {{2., 1.}, {1., 3.}}]");
    audit_true("SingularValueDecomposition#1",
               "gsvdchk[{{1., 2., 3.}, {3., 4., 5.}}, {{2., 1., 0.}, {1., 3., 1.}, {0., 1., 1.}}]");
    audit_true("SingularValueDecomposition#1", "gsvdchk[{{1. + I, 2.}, {3., 4. I}}, {{2., 1.}, {I, 3.}}]");
    audit_true("SingularValueDecomposition#1", "gsvdchk[{{1., 2.}, {3., 4.}, {5., 6.}}, {{2., 1.}, {0., 1.}}]");
    /* the generalized singular values are alpha/beta = the square roots of
     * the eigenvalues of the pencil (m^T m, a^T a) */
    audit_true("SingularValueDecomposition#1",
               "Module[{w, wa, s}, {w, wa} = SingularValueDecomposition[{{{1., 2.}, {3., 4.}}, "
               "{{2., 1.}, {1., 3.}}}][[2]]; s = Sort[Diagonal[w]/Diagonal[wa]]; "
               "Norm[s - Sort[Sqrt[Eigenvalues[{Transpose[{{1., 2.}, {3., 4.}}] . {{1., 2.}, {3., 4.}}, "
               "Transpose[{{2., 1.}, {1., 3.}}] . {{2., 1.}, {1., 3.}}}]]]] < 10^-12]");
}

/* Wave 3 S12/S13: declared argument counts and option checks
 * (src/arity.h, src/arity_table.c).  The evaluator checks a call
 * against its head's declared arity before the C code runs and
 * prints argx/argr/argrx/argt(u)/argb(u)/argm(u), nonopt or optx;
 * the call stays unevaluated. */
static void test_wave3_arity_options(void) {
    audit_chk("Det#2", "Quiet[Check[Det[], fired, Det::argx]]", "fired");
    audit_chk("Det#2", "Quiet[Det[]]", "Det[]");
    audit_chk("Eigensystem#3", "Quiet[Check[Eigensystem[], fired, Eigensystem::argt]]", "fired");
    audit_chk("Eigensystem#3", "Quiet[Eigensystem[]]", "Eigensystem[]");
    audit_chk("Norm#5", "Quiet[Check[Norm[{1, 2}, 2, 3], fired, Norm::argt]]", "fired");
    audit_chk("Norm#5", "Quiet[Norm[{1, 2}, 2, 3]]", "Norm[{1, 2}, 2, 3]");
    audit_chk("LUDecomposition#2", "Quiet[Check[LUDecomposition[{{1, 2}, {3, 4}}, x], fired, LUDecomposition::argx]]", "fired");
    audit_chk("QRDecomposition#3", "Quiet[Check[QRDecomposition[], fired, QRDecomposition::argx]]", "fired");
    audit_chk("QRDecomposition#3", "Quiet[QRDecomposition[]]", "QRDecomposition[]");
    audit_chk("Tr#1", "Quiet[Check[Tr[{{1, 2}, {3, 4}}, Plus, 1, 2], fired, Tr::argb]]", "fired");
    audit_chk("VandermondeMatrix#2", "Quiet[Check[VandermondeMatrix[{1, 2}, 2, 3], fired, VandermondeMatrix::argt]]", "fired");
    audit_chk("VandermondeMatrix#2", "Quiet[VandermondeMatrix[{1, 2}, 2, 3]]", "VandermondeMatrix[{1, 2}, 2, 3]");
    audit_chk("Eigensystem#2", "Quiet[Check[Eigensystem[{{1, 2}, {3, 4}}, Foo -> 1], fired, Eigensystem::optx]]", "fired");
    audit_chk("Eigenvalues#5", "Quiet[Check[Eigenvalues[{{1, 2}, {3, 4}}, Foo -> 1], fired, Eigenvalues::optx]]", "fired");
    audit_chk("Fit#4", "Quiet[Check[Fit[{{0, 1}, {1, 2}, {2, 3}}, {1, x}, x, Foo -> 1], fired, Fit::optx]]", "fired");
    audit_chk("Inverse#3", "Quiet[Check[Inverse[{{1}}, {{2}}], fired, Inverse::nonopt]]", "fired");
    audit_chk("Inverse#3", "Quiet[Inverse[{{1}}, {{2}}]]", "Inverse[{{1}}, {{2}}]");
    audit_chk("LinearSolve#3", "Quiet[Check[LinearSolve[{{1, 2}, {3, 4}}, {5, 6}, Foo -> 1], fired, LinearSolve::optx]]", "fired");
    audit_chk("DiagonalMatrixQ#3", "Quiet[Check[DiagonalMatrixQ[{{1}}, Foo -> 1], fired, DiagonalMatrixQ::optx]]", "fired");
    audit_chk("DiagonalMatrixQ#3", "Quiet[DiagonalMatrixQ[{{1}}, Foo -> 1]]", "DiagonalMatrixQ[{{1}}, Foo -> 1]");
    audit_chk("MatrixFunction#3", "Quiet[Check[MatrixFunction[Sin, {{1, 2}, {3, 4}}, Foo -> 1], fired, MatrixFunction::optx]]", "fired");
    audit_chk("NullSpace#3", "Quiet[Check[NullSpace[{{1, 2}, {3, 4}}, 5], fired, NullSpace::nonopt]]", "fired");
    audit_chk("NullSpace#3", "Quiet[NullSpace[{{1, 2}, {3, 4}}, 5]]", "NullSpace[{{1, 2}, {3, 4}}, 5]");
    audit_chk("QRDecomposition#4", "Quiet[Check[QRDecomposition[{{1, 2}, {3, 4}}, Foo -> 1], fired, QRDecomposition::optx]]", "fired");
    audit_chk("SchurDecomposition#4", "Quiet[Check[SchurDecomposition[{{1., 2.}, {3., 4.}}, Foo -> 1], fired, SchurDecomposition::optx]]", "fired");
    audit_chk("Eigenvalues#5", "Eigenvalues[{{2, 0}, {0, 3}}, Cubics -> True]", "{3, 2}");
    audit_chk("DiagonalMatrixQ#2", "{DiagonalMatrixQ[{{1, 0.001}, {0, 1}}, Tolerance :> 0.01], DiagonalMatrixQ[{{1, 0.001}, {0, 1}}, {Tolerance -> 0.01}]}", "{True, True}");
    audit_chk("HilbertMatrix#3", "HilbertMatrix[2, {WorkingPrecision -> MachinePrecision}] === HilbertMatrix[2, WorkingPrecision -> MachinePrecision]", "True");
    audit_chk("QRDecomposition#4", "Module[{n = 0}, Do[Quiet[Check[QRDecomposition[{{1, 2}, {3, 4}}, Foo -> 1], n++, QRDecomposition::optx]], {2}]; n]", "2");
    audit_chk("UpperTriangularMatrixQ#2", "{UpperTriangularMatrixQ[{{1, 2}, {1*^-20, 3}}, 0, {Tolerance -> 1*^-10}], UpperTriangularMatrixQ[{{1, 2}, {1*^-20, 3}}, Tolerance :> 1*^-10]}", "{True, True}");
}

/* ================================================================== */
/* Wave 3, shared fix S11: Indeterminate / infinity / pole values     */
/* (src/infinity_rules.c and the per-head rules it does not cover).   */
/* ================================================================== */

static void test_s11_special_values(void) {
    /* Norm#1: a numeric p below 1 is no norm */
    audit_true("Norm#1", "Check[Norm[{1, 2}, 1/2], \"fired\", Norm::ptype] === \"fired\"");
    audit_true("Norm#1", "Quiet[Norm[{1, 2}, 0]] === Norm[{1, 2}, 0]");
    audit_true("Norm#1", "Quiet[Norm[{1, 2}, -1]] === Norm[{1, 2}, -1]");
    audit_true("Norm#1", "Norm[{3, 4}, 1] === 7");
    /* Norm#4: Norm with an infinite entry */
    audit_true("Norm#4", "Norm[{Infinity, 1}] === Infinity");
    audit_true("Norm#4", "Norm[Infinity] === Infinity");
    /* HermitianMatrixQ#2: conjugates of I Infinity entries */
    audit_true("HermitianMatrixQ#2", "Check[HermitianMatrixQ[{{1, I Infinity}, {-I Infinity, 1}}], \"m\"] === True");
    /* Normalize#1: Normalize with an Indeterminate entry */
    audit_true("Normalize#1", "Normalize[{Indeterminate, 1}] === {Indeterminate, Indeterminate}");
}

/* Wave 4: Modulus -> p, the elimination zero test, numeric constants,
 * argument validation. */
static void test_wave4_linear_algebra(void) {
    /* Modulus -> p (linalg_mod.c) */
    audit_chk("Det#1", "Det[{{2, 3}, {5, 7}}, Modulus -> 5]", "4");
    audit_chk("Det#1", "Options[Det]", "{Modulus -> 0}");
    audit_chk("Inverse#1", "Inverse[{{1, 2}, {3, 4}}, Modulus -> 5]", "{{3, 1}, {4, 2}}");
    audit_chk("NullSpace#2", "NullSpace[{{1, 1}, {1, 1}}, Modulus -> 2]", "{{1, 1}}");
    /* Tolerance: singular values below t times the largest are zero */
    audit_true("NullSpace#1", "With[{ns = Quiet[NullSpace[{{1, 2}, {2, 4.001}}, Tolerance -> 0.01]]}, "
               "Length[ns] == 1 && Abs[Abs[ns[[1, 1]]] - 0.894] < 0.001]");
    /* GoldenRatio is not a free indeterminate */
    audit_chk("RowReduce#1", "RowReduce[{{GoldenRatio, 1}, {1, GoldenRatio - 1}}]",
              "{{1, 1/GoldenRatio}, {0, 0}}");
    audit_chk("MatrixRank#2", "MatrixRank[{{GoldenRatio, 1}, {1, GoldenRatio - 1}}]", "1");
    audit_true("PseudoInverse#1", "Max[Abs[N[PseudoInverse[{{GoldenRatio, 1}, {1, GoldenRatio - 1}}]] - "
               "N[{{GoldenRatio, 1}, {1, GoldenRatio - 1}}/5]]] < 1*^-12");
    /* ZeroTest: a cell judged zero is written as 0 */
    audit_chk("RowReduce#2", "RowReduce[{{GoldenRatio, 1}, {1, GoldenRatio - 1}}, ZeroTest -> PossibleZeroQ]",
              "{{1, 1/GoldenRatio}, {0, 0}}");
    audit_chk("RowReduce#2", "RowReduce[{{Sin[x]^2 + Cos[x]^2 - 1, 1}, {1, 0}}, "
              "ZeroTest -> (Simplify[#] === 0 &)]", "{{1, 0}, {0, 1}}");
    /* non-matrices are not "already reduced" */
    audit_chk("RowReduce#3", "Quiet[RowReduce[5]]", "RowReduce[5]");
    audit_chk("RowReduce#3", "Quiet[RowReduce[{1, 2, 3}]]", "RowReduce[{1, 2, 3}]");
    audit_chk("RowReduce#3", "Quiet[RowReduce[{{1, 2}, {3}}]]", "RowReduce[{{1, 2}, {3}}]");
    audit_true("RowReduce#3", "Check[RowReduce[5], $Failed] === $Failed");
    /* numeric constants decide definiteness exactly */
    audit_chk("PositiveDefiniteMatrixQ#3", "PositiveDefiniteMatrixQ[{{Pi, 0}, {0, Sqrt[2]}}]", "True");
    audit_chk("PositiveDefiniteMatrixQ#3", "PositiveDefiniteMatrixQ[{{Pi, 4}, {4, Sqrt[2]}}]", "False");
    audit_chk("NegativeDefiniteMatrixQ#4", "NegativeDefiniteMatrixQ[{{-Pi, 0}, {0, -Sqrt[2]}}]", "True");
    /* an MPFR zero is zero */
    audit_chk("DiagonalMatrixQ#1", "DiagonalMatrixQ[{{1, 1`30 - 1`30}, {0, 1}}]", "True");
    audit_chk("UpperTriangularMatrixQ#1", "UpperTriangularMatrixQ[{{1, 2}, {0``30, 3}}]", "True");
    /* approximate matrix with exact constants is numericized */
    audit_true("HessenbergDecomposition#1", "With[{h = HessenbergDecomposition[{{1., Pi}, {2., 3.}}]}, "
               "Max[Abs[h[[1]] . h[[2]] . ConjugateTranspose[h[[1]]] - {{1., N[Pi]}, {2., 3.}}]] < 1*^-12]");
    audit_true("HessenbergDecomposition#1",
               "Head[HessenbergDecomposition[{{1., E, 1}, {2., 3., 1}, {1, 2, 3}}]] === List");
    /* Eigenvectors' second argument is a count spec */
    audit_chk("Eigenvectors#1", "Quiet[Eigenvectors[{{1, 2}, {3, 4}}, 1.5]]", "Eigenvectors[{{1, 2}, {3, 4}}, 1.5]");
    audit_chk("Eigenvectors#1", "Quiet[Eigenvectors[{{1, 2}, {3, 4}}, x]]", "Eigenvectors[{{1, 2}, {3, 4}}, x]");
    audit_chk("Eigenvectors#1", "Quiet[Eigenvectors[{{1, 2}, {2, 4}}, {{1, 0}, {0, 0}}]]",
              "Eigenvectors[{{1, 2}, {2, 4}}, {{1, 0}, {0, 0}}]");
    /* a singular pencil direction is Indeterminate */
    audit_true("Eigenvalues#2", "With[{e = Eigenvalues[{{{1., 2.}, {2., 4.}}, {{1., 2.}, {2., 4.}}}]}, "
               "Count[e, Indeterminate] == 1 && Abs[First[Cases[e, _Real]] - 1.] < 1*^-10]");
    /* Inner needs one head */
    audit_chk("Inner#3", "Quiet[Inner[f, h[a, b], k[x, y], g]]", "Inner[f, h[a, b], k[x, y], g]");
    audit_chk("Inner#3", "Inner[f, h[a, b], h[x, y], g]", "g[f[a, x], f[b, y]]");
    /* UpTo counts */
    audit_chk("UpTo#1", "{TakeLargest[{3, 1, 2}, 2], TakeLargest[{3, 1, 2}, UpTo[5]], "
              "TakeSmallest[{3, 1, 2}, UpTo[2]], TakeLargestBy[{3, 1, 2}, Minus, UpTo[5]]}",
              "{{3, 2}, {3, 2, 1}, {1, 2}, {1, 2, 3}}");
    /* machine CharacteristicPolynomial is stable */
    audit_true("CharacteristicPolynomial#1", "With[{m = N[Table[1/(i + j - 1), {i, 6}, {j, 6}]]}, "
               "Abs[CharacteristicPolynomial[m, 0]/Det[Table[1/(i + j - 1), {i, 6}, {j, 6}]] - 1] < 1*^-4]");
    /* a replaced matrix is not solved with the stale factor */
    audit_true("LinearSolveFunction#1", "Module[{f = LinearSolve[{{1., 2.}, {3., 4.}}], g}, "
               "g = f /. {{1., 2.}, {3., 4.}} -> {{2., 0.}, {0., 2.}}; Max[Abs[g[{2., 2.}] - {1., 1.}]] < 1*^-12]");
    /* NDArray beside an exact constant */
    audit_true("NDArray#3", "Max[Abs[Normal[NDArray[{1., 2.}] + Pi] - {1. + N[Pi], 2. + N[Pi]}]] < 1*^-15");
    /* Fit through a basis without the variable is not a fit */
    audit_chk("Fit#3", "Quiet[Fit[{{0, 1}, {1, 3}, {2, 5}}, {1, x}, y]]", "Fit[{{0, 1}, {1, 3}, {2, 5}}, {1, x}, y]");
    /* DesignMatrix: the constant column, a single basis function */
    audit_chk("DesignMatrix#1", "DesignMatrix[{{0, 1}, {1, 0}, {3, 2}, {5, 4}}, {x}, x]",
              "{{1, 0}, {1, 1}, {1, 3}, {1, 5}}");
    audit_chk("DesignMatrix#2", "DesignMatrix[{{0, 1}, {1, 0}, {3, 2}, {5, 4}}, x, x]",
              "{{1, 0}, {1, 1}, {1, 3}, {1, 5}}");
    audit_chk("DesignMatrix#1", "DesignMatrix[{{0, 1}, {1, 0}, {3, 2}, {5, 4}}, {x}, x, IncludeConstantBasis -> False]",
              "{{0}, {1}, {3}, {5}}");
    /* a 40 x 40 HNF in well under the budget, and still a valid decomposition */
    audit_run("SeedRandom[1]; auditH = RandomInteger[{-100, 100}, {40, 40}]");
    audit_true("HermiteDecomposition#1", "With[{d = TimeConstrained[HermiteDecomposition[auditH], 5]}, "
               "ListQ[d] && d[[1]] . auditH == d[[2]] && Abs[Det[d[[1]]]] == 1 && "
               "UpperTriangularMatrixQ[d[[2]]]]");
    audit_run("ClearAll[auditH]");
    /* already fixed by an earlier wave; pinned here */
    audit_chk("HilbertMatrix#1", "Options[HilbertMatrix]", "{WorkingPrecision -> Infinity}");
}

/* Wave 5 (G6, low): linear-algebra message, option and bigint fixes. */
static void test_wave5_linalg(void) {
    /* Norm: a string order other than "Frobenius" is ::ptype; a scalar's
     * p-norm is its Abs. */
    audit_chk("Norm#2", "Quiet[Check[Norm[{1, 2}, \"Foo\"], fired, Norm::ptype]]", "fired");
    audit_chk("Norm#2", "Quiet[Norm[{1, 2}, \"Foo\"]]", "Norm[{1, 2}, \"Foo\"]");
    audit_chk("Norm#2", "Norm[{{1, 2}, {3, 4}}, \"Frobenius\"]", "Sqrt[30]");
    audit_chk("Norm#3", "{Norm[-3, 1], Norm[-3, 2], Norm[-3, Infinity], Norm[3 + 4 I, 3]}", "{3, 3, 3, 5}");
    /* Dot of a ragged operand */
    audit_chk("Dot#1", "Quiet[Check[{{1, 2}, {3}} . {1, 2}, fired, Dot::rect]]", "fired");
    audit_chk("Dot#1", "Quiet[{{1, 2}, {3}} . {1, 2}]", "{{1, 2}, {3}} . {1, 2}");
    audit_chk("Dot#1", "Quiet[Check[a . {1, 2}, fired]]", "a . {1, 2}");
    /* BigInt counts and indices */
    audit_chk("Diagonal#1", "{Diagonal[{{1, 2}, {3, 4}}, 2^70], Diagonal[{{1, 2}, {3, 4}}, -2^70]}", "{{}, {}}");
    audit_chk("Minors#1", "{Minors[{{1, 2}, {3, 4}}, 10^30], Minors[{{1, 2}, {3, 4}}, -10^30]}", "{{}, {}}");
    audit_chk("LUBackSubstitution#1",
              "LUBackSubstitution[{{{1, 2}, {3, 4}}, {2^32 + 1, 2}, 0}, {1, 2}]",
              "LUBackSubstitution[{{{1, 2}, {3, 4}}, {4294967297, 2}, 0}, {1, 2}]");
    audit_chk("LUBackSubstitution#1",
              "LUBackSubstitution[{{{1, 2}, {3, 4}}, {-4294967295, 2}, 0}, {1, 2}]",
              "LUBackSubstitution[{{{1, 2}, {3, 4}}, {-4294967295, 2}, 0}, {1, 2}]");
    audit_chk("MatrixPower#1", "MatrixPower[{{1, 1}, {0, 1}}, 10^30]",
              "{{1, 1000000000000000000000000000000}, {0, 1}}");
    audit_chk("MatrixPower#1", "MatrixPower[{{1, 0}, {0, 1}}, 10^30]", "{{1, 0}, {0, 1}}");
    audit_chk("MatrixPower#2", "MatrixPower[{{1, 2}, {3, 4}}, {2, 3}]", "MatrixPower[{{1, 2}, {3, 4}}, {2, 3}]");
    /* argument counts */
    audit_chk("IdentityMatrix#2", "Quiet[Check[IdentityMatrix[], fired, IdentityMatrix::argb]]", "fired");
    audit_chk("SymmetricMatrixQ#1", "Quiet[Check[SymmetricMatrixQ[], fired, SymmetricMatrixQ::argx]]", "fired");
    audit_chk("SymmetricMatrixQ#1", "SymmetricMatrixQ[{{1, 2}, {2, 1}}]", "True");
    audit_chk("Inner#4", "Quiet[Check[Inner[f, {a, b}, {x, y}, g, h, k], fired, Inner::argb]]", "fired");
    audit_chk("Inner#4", "Quiet[Inner[f, {a, b}, {x, y}, g, h, k]]", "Inner[f, {a, b}, {x, y}, g, h, k]");
    audit_chk("Inner#4", "Quiet[Check[Inner[f], fired, Inner::argbu]]", "fired");
    audit_chk("Inner#4", "Inner[f, {a, b}, {x, y}, g]", "g[f[a, x], f[b, y]]");
    /* eigen count specs and methods */
    audit_chk("Eigenvalues#4", "Quiet[Check[Eigenvalues[{{1, 2}, {3, 4}}, x], fired, Eigenvalues::seqs]]", "fired");
    audit_chk("Eigenvalues#4", "Quiet[Check[Eigenvalues[{{1, 2}, {3, 4}}, 1.5], fired, Eigenvalues::seqs]]", "fired");
    audit_chk("Eigenvalues#4", "Quiet[Check[Eigenvalues[{{1, 2}, {3, 4}}, UpTo[-1]], fired, Eigenvalues::seqs]]", "fired");
    audit_chk("Eigenvalues#4", "Quiet[Check[Eigenvalues[{{1, 2}, {2, 4}}, {{1, 2}, {2, 4}}], fired, Eigenvalues::seqs]]", "fired");
    audit_chk("Eigenvalues#4", "Quiet[Eigenvalues[{{1, 2}, {3, 4}}, x]]", "Eigenvalues[{{1, 2}, {3, 4}}, x]");
    audit_chk("Eigenvalues#4", "Eigenvalues[{{2, 0}, {0, 1}}, 1]", "{2}");
    audit_chk("Eigenvalues#4", "Eigenvalues[{{2, 0}, {0, 1}}, UpTo[5]]", "{2, 1}");
    audit_chk("Eigenvectors#2",
              "Quiet[Check[Eigenvectors[{{1., 2.}, {3., 4.}}, Method -> \"Foo\"], fired, Eigenvectors::method]]",
              "fired");
    audit_chk("CharacteristicPolynomial#2",
              "Quiet[Check[CharacteristicPolynomial[{{1, 2, 3}, {4, 5, 6}}, x], fired, CharacteristicPolynomial::matsq]]",
              "fired");
    audit_chk("CharacteristicPolynomial#2",
              "Quiet[Check[CharacteristicPolynomial[{}, x], fired, CharacteristicPolynomial::matsq]]", "fired");
    audit_chk("CharacteristicPolynomial#2", "Quiet[Check[CharacteristicPolynomial[m, x], fired]]",
              "CharacteristicPolynomial[m, x]");
    /* Cross: a non-List argument is silent in either order */
    audit_chk("Cross#1", "Check[Cross[b[x], {1, 2}], fired]", "Cross[b[x], {1, 2}]");
    audit_chk("Cross#1", "Check[Cross[{1, 2}, b[x]], fired]", "Cross[{1, 2}, b[x]]");
    audit_chk("Cross#1", "Quiet[Check[Cross[{1, 2}, {3, 4}], fired, Cross::nonn1]]", "fired");
    audit_chk("DesignMatrix#3",
              "Quiet[Check[DesignMatrix[{{1, 2}, {2, 3}}, {x}, x, WorkingPrecision -> -3], fired, DesignMatrix::wprec]]",
              "fired");
    audit_chk("HilbertMatrix#2",
              "Precision[HilbertMatrix[2, WorkingPrecision -> 30, WorkingPrecision -> Infinity]]", "30.0");
    /* options */
    audit_chk("MatrixRank#3", "MatrixRank[{{1, 1}, {1, 3}}, Modulus -> 2]", "1");
    audit_chk("MatrixRank#3", "MatrixRank[{{1, 1}, {1, 3}}, Modulus -> 3]", "2");
    audit_chk("MatrixRank#4", "MatrixRank[{{1, 1}, {1, 1}}, ZeroTest -> (PossibleZeroQ[#1] &)]", "1");
    audit_chk("MatrixRank#4", "MatrixRank[{{1, 2}, {3, 4}}, ZeroTest -> (PossibleZeroQ[#1] &)]", "2");
    audit_chk("Inverse#2", "Inverse[{{1, 2}, {3, 4}}, ZeroTest -> (#1 == 0 &)]", "{{-2, 1}, {3/2, -1/2}}");
    audit_chk("Inverse#2", "Quiet[Check[Inverse[{{1, 2}, {2, 4}}, ZeroTest -> (#1 == 0 &)], fired, Inverse::sing]]", "fired");
    /* SingularValueDecomposition messages */
    audit_chk("SingularValueDecomposition#2",
              "Quiet[Check[SingularValueDecomposition[{1, 2}], fired, SingularValueDecomposition::matrix]]", "fired");
    audit_chk("SingularValueDecomposition#3",
              "Quiet[Check[SingularValueDecomposition[{{1., 2.}, {3., 4.}}, 2^40], fired, SingularValueDecomposition::sval]]",
              "fired");
    audit_chk("SingularValueDecomposition#3",
              "Quiet[Check[SingularValueDecomposition[{{1., 2.}, {3., 4.}}, 2^70], fired, SingularValueDecomposition::sval]]",
              "fired");
    audit_chk("SingularValueDecomposition#3",
              "Quiet[Check[SingularValueDecomposition[{{1, 2}, {3, 4}}, UpTo[-1]], fired, SingularValueDecomposition::sval]]",
              "fired");
    audit_chk("SingularValueDecomposition#3",
              "Dimensions /@ SingularValueDecomposition[{{1., 2.}, {3., 4.}}, 1]", "{{2, 1}, {1, 1}, {2, 1}}");
    audit_chk("LatticeReduce#1", "Quiet[Check[LatticeReduce[{{0, 0, 0}}], fired, LatticeReduce::dep]]", "fired");
    audit_chk("LatticeReduce#1", "LatticeReduce[{{1, 0}, {0, 1}}]", "{{1, 0}, {0, 1}}");
    /* LinearSolve: one machine operand makes the system machine (contagion),
     * so the mixed form gives the same minimum-norm answer */
    audit_chk("LinearSolve#1",
              "Max[Abs[LinearSolve[{{1, 2.}, {2, 4}}, {3, 6}] - LinearSolve[{{1., 2.}, {2., 4.}}, {3., 6.}]]] < 10^-14",
              "True");
    audit_chk("LinearSolve#1", "Max[Abs[LinearSolve[{{1., 2.}, {2., 4.}}, {3., 6.}] - {0.6, 1.2}]] < 10^-14", "True");
    /* RowReduce options */
    audit_chk("RowReduce#4", "Options[RowReduce]",
              "{Method -> Automatic, Modulus -> 0, Tolerance -> Automatic, ZeroTest -> Automatic}");
    audit_chk("RowReduce#4", "RowReduce[{{1, 2}, {3, 1}}, Modulus -> 5]", "{{1, 2}, {0, 0}}");
    audit_chk("RowReduce#4", "Quiet[Check[RowReduce[{{1, 2}, {3, 4}}, Foo -> 1], fired, RowReduce::optx]]", "fired");
    audit_chk("RowReduce#4", "RowReduce[{{1, 2}, {3, 4}}, ZeroTest -> Automatic, Tolerance -> Automatic]", "{{1, 0}, {0, 1}}");
    audit_chk("RowReduce#5", "RowReduce[{{1., 2, 3}, {4, 5, 6}, {7, 8, 9}}]",
              "{{1.0, 0.0, -1.0}, {0.0, 1.0, 2.0}, {0.0, 0.0, 0.0}}");
    audit_chk("FindIntegerNullVector#1", "FindIntegerNullVector[{1, 2}, Sqrt[5]]", "{-2, 1}");
    audit_chk("MatrixFunction#2",
              "Quiet[Check[MatrixFunction[Sqrt, {{0., 1.}, {0., 0.}}], fired, MatrixFunction::nosol]]", "fired");
    audit_chk("MatrixFunction#2",
              "Quiet[Check[MatrixFunction[Sqrt, {{0., 1.}, {0., 0.}}], fired, MatrixPower::nosol]]",
              "MatrixFunction[Sqrt, {{0.0, 1.0}, {0.0, 0.0}}]");
    audit_chk("PositiveDefiniteMatrixQ#6", "PositiveDefiniteMatrixQ[N[{{1, 1}, {1, 1 + 10^-30}}, 50]]", "True");
    audit_chk("PositiveDefiniteMatrixQ#6", "PositiveDefiniteMatrixQ[N[{{1, 1}, {1, 1 - 10^-30}}, 50]]", "False");
    audit_chk("NegativeDefiniteMatrixQ#7", "NegativeDefiniteMatrixQ[N[{{-1, -1}, {-1, -1 - 10^-30}}, 50]]", "True");
    audit_chk("TensorProduct#1",
              "{TensorProduct[2, 3], TensorProduct[2, 3, {1}], TensorProduct[Pi, {1, 2}], TensorProduct[Sqrt[2], {1, 2}]}",
              "{6, {6}, {Pi, 2 Pi}, {Sqrt[2], 2 Sqrt[2]}}");
    audit_chk("TensorProduct#1", "TensorProduct[a, {1, 2}]", "TensorProduct[a, {1, 2}]");
    audit_chk("UpTo#3", "Quiet[Check[Take[{1, 2, 3}, UpTo[-2]], fired, Take::seqs]]", "fired");
    audit_chk("UpTo#3", "Quiet[Check[Drop[{1, 2, 3}, UpTo[-1]], fired, Drop::seqs]]", "fired");
    audit_chk("UpTo#3", "Quiet[Take[{1, 2, 3}, UpTo[-2]]]", "Take[{1, 2, 3}, UpTo[-2]]");
    audit_chk("UpTo#3", "{Take[{1, 2, 3}, UpTo[5]], Drop[{1, 2, 3}, UpTo[2]]}", "{{1, 2, 3}, {3}}");
    audit_chk("SingularValueList#1", "SingularValueList[{{Indeterminate, 0}, {0, 1.}}]", "{Indeterminate, Indeterminate}");
    audit_chk("SingularValueList#1", "SingularValueList[{{2., 0}, {0, 1.}}]", "{2.0, 1.0}");
    audit_chk("Eigensystem#4", "FreeQ[Eigensystem[{{Indeterminate, 1}, {0, 1}}], Infinity]", "True");
    audit_chk("Eigensystem#4", "MemberQ[First[Eigensystem[{{Indeterminate, 1}, {0, 1}}]], Indeterminate]", "True");
    audit_chk("Eigenvalues#6", "Eigenvalues[{{I, 1}, {0, 2}}]", "{2, I}");
    audit_chk("Outer#1", "Quiet[Check[Outer[f, {a, b}, g[c, d]], fired, Outer::heads]]", "fired");
    audit_chk("Outer#1", "Quiet[Outer[f, <|1 -> a|>, {x}]]", "Outer[f, <|1 -> a|>, {x}]");
    audit_chk("Outer#1", "Outer[f, {a, b}, {x}]", "{{f[a, x]}, {f[b, x]}}");
    audit_chk("Outer#1", "Outer[f, g[a, b], g[x]]", "g[g[f[a, x]], g[f[b, x]]]");
    audit_chk("CholeskyDecomposition#1",
              "c = CholeskyDecomposition[{{4, 0}, {0, 9.`20}}]; {Precision[c[[1, 1]]] >= 19, Abs[c[[1, 1]] - 2] < 10^-18}",
              "{True, True}");
    audit_chk("CholeskyDecomposition#1",
              "FreeQ[CholeskyDecomposition[{{4, 1.`20}, {1.`20, 3}}][[1, 1]], _Integer]", "True");
    audit_chk("SparseArray#1", "Quiet[Check[SparseArray[{{0, 1} -> 1}], fired, SparseArray::pos]]", "fired");
    audit_chk("SparseArray#1", "Quiet[Check[SparseArray[{{0, 1} -> 1}], fired, SparseArray::exdims]]",
              "SparseArray[{{0, 1} -> 1}]");
    audit_chk("SparseArray#1", "Normal[SparseArray[{{1, 2} -> 1}]]", "{{0, 1}}");
    audit_chk("SchurDecomposition#1",
              "Quiet[Check[SchurDecomposition[{{1, 2}, {3, 4}}], fired, SchurDecomposition::schurf]]", "fired");
    audit_chk("SchurDecomposition#3",
              "Quiet[Check[SchurDecomposition[{{1., Infinity}, {0., 1.}}], fired, SchurDecomposition::schurn]]", "fired");
    /* docstrings that now describe the implementation */
    audit_chk("Det#3", "StringContainsQ[Information[Det], \"FLINT\"] && !StringContainsQ[Information[Det], \"Bareiss\"]", "True");
    audit_chk("Adjugate#1", "StringContainsQ[Information[Adjugate], \"one LU factorisation\"]", "True");
    audit_chk("PseudoInverse#2", "StringContainsQ[Information[PseudoInverse], \"singular value\"] && !StringContainsQ[Information[PseudoInverse], \"rationalised\"]", "True");
    audit_chk("SchurDecomposition#2", "StringContainsQ[Information[SchurDecomposition], \"real or complex\"]", "True");
}

/* ---- Wave 6 G3 BEGIN ---- */
static void test_wave6_g3(void) {
    audit_chk("QRDecomposition#1", "QRDecomposition[{{qra, qrb}, {qrc, qrd}}][[1, 1]]", "{Conjugate[qra]/Sqrt[Abs[qra]^2 + Abs[qrc]^2], Conjugate[qrc]/Sqrt[Abs[qra]^2 + Abs[qrc]^2]}");
    audit_chk("QRDecomposition#1", "Block[{a, b, c, d}, Module[{q, r, q2, r2}, {q, r} = QRDecomposition[{{a, b}, {c, d}}]; q2 = q /. {a -> I, b -> 1, c -> 2, d -> 0}; r2 = r /. {a -> I, b -> 1, c -> 2, d -> 0}; Simplify[{q2 . ConjugateTranspose[q2], ConjugateTranspose[q2] . r2}]]]", "{{{1, 0}, {0, 1}}, {{I, 1}, {2, 0}}}");
    audit_chk("QRDecomposition#2", "QRDecomposition[{{0, 0}, {0, 0}}]", "{{{1, 0}}, {{0, 0}}}");
    audit_chk("QRDecomposition#2", "QRDecomposition[{{0, 0}, {0, 0}}, Pivoting -> True]", "{{{1, 0}}, {{0, 0}}, {{1, 0}, {0, 1}}}");
    audit_chk("QRDecomposition#2", "QRDecomposition[{{0, 0, 0}, {0, 0, 0}}, Pivoting -> True]", "{{{1, 0}}, {{0, 0, 0}}, {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}}");
    audit_chk("SparseArray#2", "Block[{$MessageList = {}}, SparseArray[{1 -> 1}, 3] + SparseArray[{2 -> 5}, 4]; Length[$MessageList]]", "1");
    audit_chk("SparseArray#2", "Check[SparseArray[{1 -> 1}, 3] + SparseArray[{2 -> 5}, 4], fired, Thread::tdlen]", "fired");
    audit_chk("NDArray#5", "Head[NDArray[{1., 2.}] + {10, 20}]", "NDArray");
    audit_chk("NDArray#5", "NDArray[{1., 2.}] + {10, 20} == NDArray[{11., 22.}]", "True");
    audit_chk("MatrixExp#2", "MatrixExp[{{0, -Pi}, {Pi, 0}}]", "{{-1, 0}, {0, -1}}");
    audit_chk("MatrixExp#2", "MatrixExp[{{0, -2}, {2, 0}}]", "{{Cos[2], -Sin[2]}, {Sin[2], Cos[2]}}");
    audit_chk("MatrixExp#2", "{Sqrt[-Pi^2], Sqrt[-4 Pi], Sqrt[-E^2]}", "{I Pi, (2*I) Sqrt[Pi], I E}");
    audit_chk("MatrixFunction#1", "Head[TimeConstrained[MatrixFunction[Sin, Table[1/(i+j-1), {i,4},{j,4}]], 10]]", "List");
    audit_chk("MatrixFunction#1", "Module[{h = Table[1/(i+j-1), {i,4},{j,4}]}, Max[Abs[N[MatrixFunction[Sin, h]] - MatrixFunction[Sin, N[h]]]] < 10^-12]", "True");
    audit_chk("MatrixFunction#1", "Module[{h = Table[1/(i+j-1), {i,4},{j,4}]}, Max[Abs[N[TimeConstrained[MatrixExp[h], 10]] - MatrixExp[N[h]]]] < 10^-12]", "True");
}
/* ---- Wave 6 G3 END ---- */

int main(void) {
    audit_init();
    TEST(test_wave4_linear_algebra);
    TEST(test_fit_general_norm);
    TEST(test_identity_matrix_size_overflow);
    TEST(test_diagonal_matrix_size_overflow);
    TEST(test_hankel_matrix_size_overflow);
    TEST(test_toeplitz_matrix_size_overflow);
    TEST(test_vandermonde_matrix_huge_k);
    TEST(test_posdef_exact);
    TEST(test_negdef_sparse);
    TEST(test_negdef_exact);
    TEST(test_matrixrank_scaled);
    TEST(test_eigensystem_mpfr_defective);
    TEST(test_eigenvalues_singular_pencil);
    TEST(test_fit_large_response);
    TEST(test_leastsquares_complex_krylov);
    TEST(test_matrixexp_irrational_diagonal);
    TEST(test_inner_incommensurate);
    TEST(test_inner_index_n);
    TEST(test_ndarray_infinity);
    TEST(test_ndarray_nonfinite);
    TEST(test_gsvd_diagonal);
    TEST(test_wave3_arity_options);
    TEST(test_s11_special_values);
    TEST(test_wave5_linalg);
    TEST(test_wave6_g3);
    return audit_finish("audit_linear_algebra");
}
