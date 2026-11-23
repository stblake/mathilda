/* Builtin-audit regressions: graphics (review/issues/; tests/test_audit.h). */
#include "test_audit.h"

/* Plot3D#1 (high): a huge PlotPoints crashed the kernel -- build_grid's
 * n*n GridPt malloc wrapped (2^32 squared is 0 mod 2^64) and was written
 * through.  The grid is checked first: past memory aborts with General::nomem. */
static void test_plot3d_huge_plot_points(void) {
    audit_chk("Plot3D#1",
              "Quiet[CheckAbort[Head[Plot3D[x + y, {x, 0, 1}, {y, 0, 1}, PlotPoints -> 2^32]], aborted]]",
              "aborted");
    audit_chk("Plot3D#1",
              "Quiet[Check[CheckAbort[Plot3D[x + y, {x, 0, 1}, {y, 0, 1}, PlotPoints -> 2^31], 0], "
              "fired, General::nomem]]",
              "fired");
    audit_chk("Plot3D#1", "Head[Plot3D[x + y, {x, 0, 1}, {y, 0, 1}, PlotPoints -> 10]]", "Graphics3D");
}

int main(void) {
    audit_init();
    TEST(test_plot3d_huge_plot_points);
    return audit_finish("audit_graphics");
}
