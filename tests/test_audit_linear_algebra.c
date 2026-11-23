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

int main(void) {
    audit_init();
    TEST(test_fit_general_norm);
    return audit_finish("audit_linear_algebra");
}
