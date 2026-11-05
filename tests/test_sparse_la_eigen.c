/* test_sparse_la_eigen.c -- a few eigenpairs / singular values of machine
 * SparseArrays and the dense Method -> "Arnoldi" (src/sparse/sparse_la_eigen.c,
 * src/linalg/krylov_schur.c).
 *
 * Every check runs twice: with the default backend (ARPACK when the build
 * links it) and with the native Krylov-Schur forced (ks_set_backend --
 * the same switch as MATHILDA_NO_ARPACK=1 in the environment). Unique answers
 * are compared with the dense LAPACK answer or with values printed by the
 * reference system (wolframscript); eigenvectors, which are unique only up to
 * sign / phase, through their residuals ||A v - lambda v|| and (Hermitian)
 * orthonormality. Messages are the reference's: arh, arhm, take, moptx,
 * eslct, arbsym, ssing, arbss, wrgopt, svecnm, maxit2.
 *
 * Failures are counted rather than asserted (Release builds pass -DNDEBUG). */

#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include "krylov_schur.h"

static int failures = 0;
static const char* g_backend = "default";

static void chk(const char* input, const char* expected) {
    struct Expr* p = parse_expression(input);
    if (!p) { fprintf(stderr, "FAIL(parse): %s\n", input); failures++; return; }
    struct Expr* e = evaluate(p);
    expr_free(p);
    char* s = expr_to_string(e);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL [%s]: %s\n  Expected: %s\n  Actual:   %s\n", g_backend, input, expected, s);
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
/* ... and no message at all */
static void chk_nomsg(const char* input) {
    char buf[4096];
    snprintf(buf, sizeof buf, "Quiet[Check[%s; \"none\", \"msg\"]]", input);
    chk(buf, "\"none\"");
}

static void setup(void) {
    chk("lap[n_] := SparseArray[{Band[{1,1}] -> 2., Band[{1,2}] -> -1., Band[{2,1}] -> -1.}, {n, n}]; "
        "nsy = SparseArray[{Band[{1,1}] -> Range[60] - 30.5, Band[{1,2}] -> 1., Band[{2,1}] -> -0.5, "
        "  Band[{1, 5}] -> 2.}, {60, 60}]; "
        "herm = SparseArray[{Band[{1,1}] -> 2. + 0. I, Band[{1,2}] -> -1. I, Band[{2,1}] -> 1. I}, {40, 40}]; "
        "cgen = SparseArray[{Band[{1,1}] -> N[Range[50]] (1 + 0.5 I), Band[{1,2}] -> 1. + 2. I, "
        "  Band[{3,1}] -> -0.7 I}, {50, 50}]; "
        "rect = SparseArray[{Band[{1,1}] -> N[Range[80]], Band[{1,2}] -> 1., Band[{3,1}] -> -0.5}, {100, 80}]; "
        "close[a_, b_, t_] := Max[Abs[a - b]] <= t Max[1, Max[Abs[b]]]; "
        "resid[a_, vals_, vecs_] := Max[Table[Norm[a . vecs[[i]] - vals[[i]] vecs[[i]]], {i, Length[vals]}]]/"
        "  Max[Abs[Flatten[Normal[a]]]]; "
        "orthoQ[v_] := Max[Abs[Flatten[Conjugate[v] . Transpose[v] - IdentityMatrix[Length[v]]]]] < 10^-10; "
        "Null", "Null");
}

/* ---- real symmetric --------------------------------------------------------- */

static void test_symmetric(void) {
    /* reference: {3.99620665747409, 3.984841019343869, 3.965946199367808} */
    chk_true("close[Eigenvalues[lap[50], 3], {3.99620665747409, 3.984841019343869, 3.965946199367808}, 10^-12]");
    chk_true("close[Eigenvalues[lap[50], 3], Take[Eigenvalues[Normal[lap[50]]], 3], 10^-12]");
    /* reference: {0.0340538006321964, 0.015158980656128472, 0.0037933425259118595} */
    chk_true("close[Eigenvalues[lap[50], -3], {0.0340538006321964, 0.015158980656128472, 0.0037933425259118595}, 10^-11]");
    chk_true("Eigenvalues[lap[50], UpTo[3]] == Eigenvalues[lap[50], 3]");
    chk("Head[Eigenvalues[lap[50], 3]]", "List");
    chk("Length[Eigenvalues[lap[50], UpTo[60]]]", "50");
    /* vectors: residual, unit norm, orthonormal; Eigensystem agrees */
    chk_true("With[{v = Eigenvectors[lap[50], 4]}, resid[lap[50], Eigenvalues[lap[50], 4], v] < 10^-12 && orthoQ[v]]");
    chk_true("With[{v = Eigenvectors[lap[50], -4]}, resid[lap[50], Eigenvalues[lap[50], -4], v] < 10^-12 && orthoQ[v]]");
    chk_true("With[{es = Eigensystem[lap[60], 3]}, resid[lap[60], es[[1]], es[[2]]] < 10^-12]");
    chk("Length /@ Eigensystem[lap[60], -2]", "{2, 2}");
    /* the sign convention: the largest component is positive */
    chk_true("With[{v = Eigenvectors[lap[50], 2]}, And @@ (Function[x, Re[x[[First[Select[Range[Length[x]], Abs[x[[#]]] >= (1 - 10^-8) Max[Abs[x]] &]]]]] > 0] /@ v)]");
    /* shift: the eigenvalues nearest 1.05, by decreasing magnitude
     * (reference: {1.1521712185802786, 1.0600465139453599, 0.9704143969803384}) */
    chk_true("close[Eigenvalues[lap[60], 3, Method -> {\"Arnoldi\", \"Shift\" -> 1.05}], "
             "{1.1521712185802786, 1.0600465139453599, 0.9704143969803384}, 10^-12]");
    chk_true("Eigenvalues[lap[60], -3, Method -> {\"Arnoldi\", \"Shift\" -> 1.05}] == "
             "Eigenvalues[lap[60], 3, Method -> {\"Arnoldi\", \"Shift\" -> 1.05}]");
    /* BothEnds: 2 from the top, 2 from the bottom (reference values) */
    chk_true("close[Eigenvalues[lap[60], 4, Method -> {\"Arnoldi\", \"Criteria\" -> \"BothEnds\"}], "
             "{3.997348179769651, 3.9893997512291817, 0.010600248770819672, 0.0026518202303410336}, 10^-11]");
    chk_true("close[Eigenvalues[lap[50], 3, Method -> {\"Arnoldi\", \"Criteria\" -> \"BothEnds\"}], "
             "{3.9962066574740884, 3.9848410193438775, 0.0037933425259108243}, 10^-11]");
    /* RealPart on a symmetric matrix: largest algebraic; with -k the
     * values just above 0 (shift-invert about 0, reference) */
    chk_true("close[Eigenvalues[lap[60], 3, Method -> {\"Arnoldi\", \"Shift\" -> 1.05, \"Criteria\" -> \"RealPart\"}], "
             "{1.2465442127296298, 1.1521712185802784, 1.0600465139453599}, 10^-12]");
    /* the other sub-options are accepted */
    chk_true("close[Eigenvalues[lap[80], 3, Method -> {\"Arnoldi\", \"BasisSize\" -> 30, \"MaxIterations\" -> 5000, "
             "\"Tolerance\" -> 10^-14, \"StartingVector\" -> N[Range[80]]}], Take[Eigenvalues[Normal[lap[80]]], 3], 10^-11]");
    /* a singular matrix: -k shifts off 0 by itself (reference answers too) */
    chk_true("With[{s = lap[60] - SparseArray[{{1,1} -> 1., {60, 60} -> 1.}, {60, 60}]}, "
             "close[Eigenvalues[s, -3], Take[Eigenvalues[Normal[s]], -3], 10^-10]]");
    chk_nomsg("Eigenvalues[lap[60] - SparseArray[{{1,1} -> 1., {60, 60} -> 1.}, {60, 60}], -3]");
    /* a machine zero matrix */
    chk("Eigenvalues[SparseArray[{}, {149, 149}, 0.], 2]", "{0.0, 0.0}");
    chk_true("Max[Abs[Eigenvalues[SparseArray[{}, {149, 149}, 0.], -2]]] < 10^-12");
}

/* ---- real nonsymmetric ------------------------------------------------------ */

static void test_nonsymmetric(void) {
    /* k = 3 splits the second conjugate pair, as the dense answer does */
    chk_true("close[Eigenvalues[nsy, 3], Take[Eigenvalues[Normal[nsy]], 3], 10^-10]");
    chk_true("With[{e = Eigenvalues[nsy, 3]}, Im[e[[1]]] > 0 && e[[2]] == Conjugate[e[[1]]]]");
    chk_true("close[Abs[Eigenvalues[nsy, -4]], Abs[Take[Eigenvalues[Normal[nsy]], -4]], 10^-10]");
    /* RealPart: the largest real parts (reference: {28.749 + 0.26 I, 28.749 - 0.26 I, 27.51}) */
    chk_true("close[Eigenvalues[nsy, 3, Method -> {\"Arnoldi\", \"Criteria\" -> \"RealPart\"}], "
             "{28.749086167060778 + 0.2597731662585135 I, 28.749086167060778 - 0.2597731662585135 I, 27.509848307452515}, 10^-8]");
    /* -k RealPart: shift-invert about 0, largest Re(1/lambda) (reference: {2.5, 1.5, 0.5}) */
    chk_true("close[Eigenvalues[nsy, -3, Method -> {\"Arnoldi\", \"Criteria\" -> \"RealPart\"}], {2.5, 1.5, 0.5}, 10^-8]");
    chk_true("close[Eigenvalues[nsy, 2, Method -> {\"Arnoldi\", \"Shift\" -> 3.}], {3.5, 2.5}, 10^-8]");
    chk_true("With[{es = Eigensystem[nsy, 4]}, resid[nsy, es[[1]], es[[2]]] < 10^-11]");
    chk_true("With[{es = Eigensystem[nsy, -4]}, resid[nsy, es[[1]], es[[2]]] < 10^-11]");
}

/* ---- complex ------------------------------------------------------------- */

static void test_complex(void) {
    /* Hermitian: real eigenvalues (reference: {3.9941316023674958, 3.9765608475607257, 3.947390847755591}) */
    chk_true("close[Eigenvalues[herm, 3], {3.9941316023674958, 3.9765608475607257, 3.947390847755591}, 10^-12]");
    chk_true("And @@ (Head[#] === Real & /@ Eigenvalues[herm, 3])");
    chk_true("With[{v = Eigenvectors[herm, 3]}, resid[herm, Eigenvalues[herm, 3], v] < 10^-12 && orthoQ[v]]");
    chk_true("close[Eigenvalues[herm, -3], Take[Eigenvalues[Normal[herm]], -3], 10^-10]");
    /* general complex */
    chk_true("close[Eigenvalues[cgen, 3], Take[Eigenvalues[Normal[cgen]], 3], 10^-10]");
    chk_true("close[Eigenvalues[cgen, -3], Take[Eigenvalues[Normal[cgen]], -3], 10^-10]");
    chk_true("With[{es = Eigensystem[cgen, 3]}, resid[cgen, es[[1]], es[[2]]] < 10^-12]");
    /* the phase convention: the largest component real and positive */
    chk_true("With[{v = Eigenvectors[cgen, 2]}, And @@ (Function[x, With[{c = x[[First[Select[Range[Length[x]], Abs[x[[#]]] >= (1 - 10^-8) Max[Abs[x]] &]]]]}, "
             "Re[c] > 0 && Im[c] == 0]] /@ v)]");
    /* a complex shift on a real matrix */
    chk_true("close[Eigenvalues[lap[40], 3, Method -> {\"Arnoldi\", \"Shift\" -> 1 + I}], "
             "{1.18186272565732, 1.0445603629754738, 0.9128648999975573}, 10^-10]");
}

/* ---- the generalized problem --------------------------------------------- */

static void test_generalized(void) {
    chk("b60 = SparseArray[DiagonalMatrix[N[Range[60]]]]; Null", "Null");
    /* reference: {2.388015445668369, 1.089974742672494, 0.7056602132332803} */
    chk_true("close[Eigenvalues[{lap[60], b60}, 3, Method -> \"Arnoldi\"], "
             "{2.388015445668369, 1.089974742672494, 0.7056602132332803}, 10^-10]");
    chk_true("close[Eigenvalues[{lap[60], b60}, -3, Method -> \"Arnoldi\"], "
             "Take[Eigenvalues[{Normal[lap[60]], Normal[b60]}], -3], 10^-10]");
    chk_true("With[{es = Eigensystem[{nsy, b60}, 3, Method -> \"Arnoldi\"]}, "
             "Max[Table[Norm[nsy . es[[2, i]] - es[[1, i]] b60 . es[[2, i]]], {i, 3}]] < 10^-9]");
    /* Automatic: the dense QZ under the cap -- the same answer */
    chk_true("close[Eigenvalues[{lap[60], b60}, 3], Eigenvalues[{Normal[lap[60]], Normal[b60]}, 3], 10^-12]");
}

/* ---- singular values -------------------------------------------------- */

static void test_singular(void) {
    /* reference: {80.21590722726376, 79.09089936024506, 78.01843742833465} */
    chk_true("close[SingularValueList[rect, 3], {80.21590722726376, 79.09089936024506, 78.01843742833465}, 10^-12]");
    chk_true("close[SingularValueList[rect, 3], Take[SingularValueList[Normal[rect]], 3], 10^-12]");
    chk_true("close[SingularValueList[rect, -3], Take[SingularValueList[Normal[rect]], -3], 10^-10]");
    chk_true("close[SingularValueList[Transpose[rect], 3], Take[SingularValueList[Normal[rect]], 3], 10^-12]");
    chk_true("close[SingularValueList[Transpose[rect], -2], Take[SingularValueList[Normal[rect]], -2], 10^-10]");
    chk_true("close[SingularValueList[cgen, 2], Take[SingularValueList[Normal[cgen]], 2], 10^-12]");
    chk_true("close[SingularValueList[cgen, -2], Take[SingularValueList[Normal[cgen]], -2], 10^-10]");
    /* {List, SparseArray, List}, u . w . v^H the rank-k part */
    chk("Head /@ SingularValueDecomposition[rect, 2]", "{List, SparseArray, List}");
    chk("Dimensions /@ SingularValueDecomposition[rect, 2]", "{{100, 2}, {2, 2}, {80, 2}}");
    chk_true("With[{d = SingularValueDecomposition[rect, 3]}, Max[Abs[rect . d[[3]] - d[[1]] . d[[2]]]] < 10^-11 && "
             "Max[Abs[ConjugateTranspose[rect] . d[[1]] - d[[3]] . d[[2]]]] < 10^-11]");
    chk_true("With[{d = SingularValueDecomposition[cgen, -2]}, Max[Abs[cgen . d[[3]] - d[[1]] . d[[2]]]] < 10^-10]");
    chk("Head /@ SingularValueDecomposition[rect]", "{List, SparseArray, List}");
    chk("Head /@ SingularValueDecomposition[SparseArray[{{1., 0.}, {0., 2.}}]]", "{List, SparseArray, List}");
    chk_true("With[{d = SingularValueDecomposition[rect]}, Max[Abs[d[[1]] . d[[2]] . ConjugateTranspose[d[[3]]] - rect]] < 10^-11]");
    /* a value under the tolerance is reported as 0. */
    chk("SingularValueList[SparseArray[{}, {40, 30}, 0.], 2]", "{0.0, 0.0}");
}

/* ---- the reference's messages ------------------------------------------- */

static void test_messages(void) {
    chk_msg("Eigenvalues[lap[50], 3, Method -> {\"Arnoldi\", \"Foo\" -> 1}]", "Eigenvalues::moptx");
    chk_msg("Eigensystem[lap[40], 1, Method -> {\"Arnoldi\", \"Foo\" -> 1}]", "Eigensystem::moptx");
    chk_msg("SingularValueList[rect, 3, Method -> {\"Arnoldi\", \"Foo\" -> 1}]", "SingularValueList::moptx");
    chk("Head[Quiet[Eigenvalues[lap[50], 3, Method -> {\"Arnoldi\", \"Foo\" -> 1}]]]", "Eigenvalues");
    chk_msg("Eigenvalues[lap[50], 3, Method -> {\"Arnoldi\", \"Criteria\" -> \"Foo\"}]", "Eigenvalues::eslct");
    chk_msg("Eigenvalues[nsy, 3, Method -> {\"Arnoldi\", \"Criteria\" -> \"BothEnds\"}]", "Eigenvalues::arbsym");
    chk_msg("Eigenvalues[herm, 3, Method -> {\"Arnoldi\", \"Criteria\" -> \"BothEnds\"}]", "Eigenvalues::arbsym");
    chk_msg("Eigenvalues[lap[50], 3, Method -> {\"Arnoldi\", \"Shift\" -> 1.}]", "Eigenvalues::ssing");
    chk("Head[Quiet[Eigenvalues[lap[50], 3, Method -> {\"Arnoldi\", \"Shift\" -> 1.}]]]", "Eigenvalues");
    chk_msg("Eigenvalues[lap[40], 3, Method -> {\"Arnoldi\", \"BasisSize\" -> 2}]", "Eigenvalues::arbss");
    chk_msg("Eigenvalues[lap[40], 3, Method -> {\"Arnoldi\", \"BasisSize\" -> 41}]", "Eigenvalues::arbss");
    chk_msg("Eigenvalues[lap[40], 3, Method -> {\"Arnoldi\", \"MaxIterations\" -> -1}]", "Eigenvalues::wrgopt");
    chk_msg("Eigenvalues[lap[40], 3, Method -> {\"Arnoldi\", \"Tolerance\" -> -1}]", "Eigenvalues::wrgopt");
    chk_msg("Eigenvalues[lap[40], 3, Method -> {\"Arnoldi\", \"Shift\" -> x}]", "Eigenvalues::wrgopt");
    chk_msg("Eigenvalues[lap[40], 3, Method -> {\"Arnoldi\", \"StartingVector\" -> {1, 2}}]", "Eigenvalues::svecnm");
    chk_msg("Eigenvalues[lap[400], 401]", "Eigenvalues::take");
    chk_msg("Eigenvalues[lap[40], -41]", "Eigenvalues::take");
    chk_msg("Eigenvectors[lap[40], 41]", "Eigenvectors::take");
    chk_msg("Eigensystem[lap[40], 41]", "Eigensystem::take");
    chk_msg("SingularValueList[lap[40], -41]", "SingularValueList::take");
    chk_msg("SingularValueDecomposition[lap[40], 41]", "SingularValueDecomposition::take");
    /* ::arh: a dense answer for a large, sparse matrix (n > 100, at most
     * a fifth stored), with a count over n/5 or the whole spectrum */
    chk_msg("Eigenvalues[lap[400], 100]", "Eigenvalues::arh");
    chk_msg("Eigenvectors[lap[400], 100]", "Eigenvectors::arh");
    chk_msg("Eigensystem[lap[400], 100]", "Eigensystem::arh");
    chk_msg("Eigenvalues[lap[101]]", "Eigenvalues::arh");
    chk_nomsg("Eigenvalues[lap[100]]");
    chk_nomsg("Eigenvalues[lap[400], 80]");
    chk_msg("SingularValueList[lap[150]]", "SingularValueList::arh");
    chk_msg("SingularValueDecomposition[lap[150]]", "SingularValueDecomposition::arh");
    /* 300 x 80: the count is measured against min(m, n), the size test
     * against max(m, n) (reference) */
    chk("tall = SparseArray[{Band[{1,1}] -> N[Range[80]], Band[{1,2}] -> 1., Band[{3,1}] -> -0.5}, {300, 80}]; Null", "Null");
    chk_msg("SingularValueList[tall, 17]", "SingularValueList::arh");
    chk_nomsg("SingularValueList[tall, 16]");
    chk_nomsg("SingularValueDecomposition[lap[150], 2]");
    /* exact input: ::arhm, and the zero matrix's -k notice */
    chk_msg("Eigenvalues[SparseArray[{{1,1}->1}, {120, 120}]]", "Eigenvalues::arhm");
    chk_msg("Eigensystem[SparseArray[{}, {149, 149}]]", "Eigensystem::arhm");
    chk_msg("SingularValueList[SparseArray[{}, {149, 149}]]", "SingularValueList::arhm");
    chk_msg("Eigenvalues[SparseArray[{}, {184, 184}], -2]", "Eigenvalues::ssing");
    chk("Quiet[Eigenvalues[SparseArray[{}, {20, 20}], -2]]", "{0, 0}");
    /* the exact zero matrix with a count (reference answers) */
    chk("Eigenvectors[SparseArray[{}, {6, 6}], 2]", "{{0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 1, 0}}");
    chk("Quiet[Eigensystem[SparseArray[{}, {6, 6}], -2]]", "{{0, 0}, {{0, 1, 0, 0, 0, 0}, {1, 0, 0, 0, 0, 0}}}");
    chk("Length[Eigensystem[SparseArray[{}, {400, 400}], 2][[2, 1]]]", "400");
    chk("SingularValueList[SparseArray[{}, {283, 338}], 2]", "{0, 0}");
    /* the iteration limit: the current values come back */
    chk_msg("Eigenvalues[lap[2000], 3, Method -> {\"Arnoldi\", \"MaxIterations\" -> 1, \"BasisSize\" -> 6}]",
            "Eigenvalues::maxit2");
    chk("Length[Quiet[Eigenvalues[lap[2000], 3, Method -> {\"Arnoldi\", \"MaxIterations\" -> 1, \"BasisSize\" -> 6}]]]", "3");
}

/* ---- far above the densify cap ------------------------------------------- */

static void test_large(void) {
    /* the 2-D Laplacian on a 250 x 250 grid: 62500 unknowns;
     * lambda_ij = 4 sin^2(i pi / 502) + 4 sin^2(j pi / 502) */
    chk("t250 = SparseArray[{Band[{1,1}] -> 2., Band[{1,2}] -> -1., Band[{2,1}] -> -1.}, {250, 250}]; "
        "i250 = SparseArray[{Band[{1,1}] -> 1.}, {250, 250}]; "
        "L = KroneckerProduct[t250, i250] + KroneckerProduct[i250, t250]; "
        "ex = Sort[N[Flatten[Table[4 Sin[i Pi/502]^2 + 4 Sin[j Pi/502]^2, {i, 3}, {j, 3}]]]]; Null", "Null");
    chk_true("close[Eigenvalues[L, -4], Reverse[Take[ex, 4]], 10^-9]");
    chk_true("With[{es = Eigensystem[L, -2]}, Max[Table[Norm[L . es[[2, i]] - es[[1, i]] es[[2, i]]], {i, 2}]] < 10^-10]");
    /* a 50000 x 50000 nonsymmetric tridiagonal with separated eigenvalues */
    chk("nb = SparseArray[{Band[{1,1}] -> N[Range[50000]], Band[{1,2}] -> 0.3, Band[{2,1}] -> -0.2}, {50000, 50000}]; Null", "Null");
    chk_true("With[{es = Eigensystem[nb, -2]}, Max[Table[Norm[nb . es[[2, i]] - es[[1, i]] es[[2, i]]], {i, 2}]] < 10^-9]");
    chk_true("With[{s = SingularValueList[nb, -2]}, Length[s] == 2 && s[[1]] > s[[2]] > 0.5 && s[[1]] < 2.5]");
    chk_true("With[{d = SingularValueDecomposition[nb, -2]}, Max[Abs[nb . d[[3]] - d[[1]] . d[[2]]]] < 10^-9]");
    chk_true("With[{es = Eigensystem[nb, 2, Method -> {\"Arnoldi\", \"Shift\" -> 25000.3}]}, Max[Table[Norm[nb . es[[2, i]] - es[[1, i]] es[[2, i]]], {i, 2}]] < 10^-7]");
}

/* ---- the dense Method -> "Arnoldi" ------------------------------------------ */

static void test_dense_arnoldi(void) {
    chk("d60 = Normal[lap[60]]; dn = Normal[nsy]; Null", "Null");
    chk_true("close[Eigenvalues[d60, 3, Method -> \"Arnoldi\"], Eigenvalues[d60, 3, Method -> \"Direct\"], 10^-12]");
    chk_true("close[Eigenvalues[d60, -3, Method -> \"Arnoldi\"], Eigenvalues[d60, -3, Method -> \"Direct\"], 10^-11]");
    chk_true("close[Eigenvalues[dn, 3, Method -> \"Arnoldi\"], Eigenvalues[dn, 3], 10^-10]");
    chk_true("With[{v = Eigenvectors[d60, 3, Method -> \"Arnoldi\"]}, "
             "resid[d60, Eigenvalues[d60, 3], v] < 10^-12 && orthoQ[v]]");
    chk_true("With[{es = Eigensystem[dn, 4, Method -> {\"Arnoldi\", \"Criteria\" -> \"RealPart\"}]}, resid[dn, es[[1]], es[[2]]] < 10^-11]");
    /* Automatic with a count is Direct, exactly the first k of the spectrum */
    chk_true("Eigenvalues[d60, 3] === Take[Eigenvalues[d60], 3]");
    /* a clustered symmetric 120 x 120: residuals were 2e-2 with the single pass */
    chk("cl = N[Table[If[i == j, 1 + i/1000, 1/(1 + Abs[i - j])^3], {i, 120}, {j, 120}]]; Null", "Null");
    chk_true("With[{v = Eigenvectors[cl, 3, Method -> \"Arnoldi\"]}, resid[cl, Eigenvalues[cl, 3], v] < 10^-12]");
}

#include <time.h>
static void timed(const char* what, void (*f)(void)) {
    clock_t t0 = clock();
    f();
    if (getenv("SLA_EIGEN_TIMES"))
        fprintf(stderr, "  [%s] %s: %.2f s cpu\n", g_backend, what, (double)(clock() - t0) / CLOCKS_PER_SEC);
}

static void run_all(const char* name) {
    g_backend = name;
    timed("symmetric", test_symmetric);
    timed("nonsymmetric", test_nonsymmetric);
    timed("complex", test_complex);
    timed("generalized", test_generalized);
    timed("singular", test_singular);
    timed("messages", test_messages);
    timed("large", test_large);
    timed("dense_arnoldi", test_dense_arnoldi);
}

int main(void) {
    symtab_init();
    core_init();
    test_mute_messages();
    setup();
    run_all(ks_arpack_active() ? "ARPACK" : "native");
    if (ks_arpack_active()) {
        ks_set_backend(KS_BACKEND_NATIVE);
        run_all("native");
        ks_set_backend(KS_BACKEND_AUTO);
    }
    if (failures) {
        fprintf(stderr, "%d failure(s)\n", failures);
        return 1;
    }
    printf("sparse_la_eigen_tests: all passed\n");
    return 0;
}
