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

int main(void) {
    audit_init();
    TEST(test_fit_general_norm);
    TEST(test_identity_matrix_size_overflow);
    TEST(test_diagonal_matrix_size_overflow);
    TEST(test_hankel_matrix_size_overflow);
    TEST(test_toeplitz_matrix_size_overflow);
    TEST(test_vandermonde_matrix_huge_k);
    return audit_finish("audit_linear_algebra");
}
