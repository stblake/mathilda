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

int main(void) {
    audit_init();
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
    return audit_finish("audit_linear_algebra");
}
