/*
 * tests/test_lapack.c
 *
 * Phase-1 smoke test for the BLAS/LAPACK linkage wired in by
 * the QR performance plan (tasks/qr_perf_plan.md, Section 5).
 *
 * These tests do not exercise any QR code yet; they prove that
 *   - mathilda_lapack_probe() reports the correct value for the
 *     current build mode (USE_LAPACK on/off),
 *   - on platforms that should have BLAS, cblas_ddot returns
 *     the right answer for a hand-checkable inner product,
 *   - LAPACK Fortran-ABI symbols are linkable on every supported
 *     platform (we call dgeqrf_ on a 1x1 system as a dynamic-link
 *     smoke test; we do not validate the QR result yet — that
 *     comes in Phase 3).
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "test_utils.h"
#include "lapack.h"

#include "core.h"
#include "symtab.h"

static int g_failures = 0;

#define EXPECT(cond, msg) do { \
    if (!(cond)) { \
        fprintf(stderr, "FAIL: %s (%s:%d)\n", (msg), __FILE__, __LINE__); \
        g_failures++; \
    } else { \
        printf("PASS: %s\n", (msg)); \
    } \
} while (0)

/*
 * The probe must agree with the build configuration: 1 if USE_LAPACK
 * is defined at compile time, 0 otherwise.  Anything else means the
 * wrapper or the build wiring is broken.
 */
static void test_probe_matches_build_flag(void) {
    int p = mathilda_lapack_probe();
#ifdef USE_LAPACK
    EXPECT(p == 1, "mathilda_lapack_probe() returns 1 under USE_LAPACK");
#else
    EXPECT(p == 0, "mathilda_lapack_probe() returns 0 without USE_LAPACK");
#endif
}

/*
 * cblas_ddot smoke tests.  Only meaningful when USE_LAPACK is on; on
 * the OFF build we just record a PASS so the test count stays
 * stable for CI dashboards.
 */
static void test_cblas_ddot_orthogonal(void) {
#ifdef USE_LAPACK
    const double e1[3] = { 1.0, 0.0, 0.0 };
    const double e2[3] = { 0.0, 1.0, 0.0 };
    double r = cblas_ddot(3, e1, 1, e2, 1);
    EXPECT(r == 0.0, "cblas_ddot of orthogonal e1.e2 == 0");
#else
    EXPECT(1, "cblas_ddot orthogonal (skipped: USE_LAPACK off)");
#endif
}

static void test_cblas_ddot_parallel(void) {
#ifdef USE_LAPACK
    const double v[4] = { 1.0, 2.0, 3.0, 4.0 };
    double r = cblas_ddot(4, v, 1, v, 1);
    /* 1+4+9+16 == 30 */
    EXPECT(r == 30.0, "cblas_ddot of v.v == sum of squares");
#else
    EXPECT(1, "cblas_ddot parallel (skipped: USE_LAPACK off)");
#endif
}

static void test_cblas_ddot_strided(void) {
#ifdef USE_LAPACK
    /* Interleaved layout: pick every other element. */
    const double x[6] = { 1.0, 99.0, 2.0, 99.0, 3.0, 99.0 };
    const double y[6] = { 1.0, 99.0, 2.0, 99.0, 3.0, 99.0 };
    double r = cblas_ddot(3, x, 2, y, 2);
    EXPECT(r == 14.0, "cblas_ddot with stride=2 picks correct entries");
#else
    EXPECT(1, "cblas_ddot strided (skipped: USE_LAPACK off)");
#endif
}

/*
 * LAPACK Fortran-ABI link test.  Call dgeqrf_ on a tiny 2x2 matrix
 * and only check that `info == 0` — we are not verifying the QR
 * factorisation here, only that the symbol is resolvable and the
 * Fortran calling convention matches our declarations.
 *
 * On macOS the prototype comes from <Accelerate/Accelerate.h>; on
 * Linux from our own extern declaration in lapack.h.
 */
static void test_lapack_dgeqrf_link(void) {
#ifdef USE_LAPACK
    /* Column-major 2x2: [[1, 0], [0, 1]] -> stored as {1, 0, 0, 1} */
    double A[4] = { 1.0, 0.0, 0.0, 1.0 };
    double tau[2] = { 0.0, 0.0 };
    /* Workspace query: lwork=-1 makes dgeqrf write the optimal size
     * into work[0]. */
    double work_query = 0.0;
    int    m = 2, n = 2, lda = 2, lwork_q = -1, info = 0;
    dgeqrf_(&m, &n, A, &lda, tau, &work_query, &lwork_q, &info);
    EXPECT(info == 0, "dgeqrf_ workspace query returns info==0");

    int lwork = (int)work_query;
    if (lwork < 2) lwork = 2;
    double* work = (double*)malloc((size_t)lwork * sizeof(double));
    info = 0;
    dgeqrf_(&m, &n, A, &lda, tau, work, &lwork, &info);
    EXPECT(info == 0, "dgeqrf_ on 2x2 identity returns info==0");
    /* QR of identity: R should be diagonal with |r_ii|==1.  In LAPACK
     * the upper triangle of A becomes R, so A[0]==R[0,0]. */
    EXPECT(fabs(fabs(A[0]) - 1.0) < 1e-14, "dgeqrf_ R[0,0] has unit magnitude for I");
    EXPECT(fabs(fabs(A[3]) - 1.0) < 1e-14, "dgeqrf_ R[1,1] has unit magnitude for I");
    free(work);
#else
    EXPECT(1, "dgeqrf link (skipped: USE_LAPACK off)");
#endif
}

int main(void) {
    /* The smoke test only needs the symbol table for ASSERT_* infra;
     * full core_init keeps test_utils.h happy for future tests in
     * this file. */
    symtab_init();
    core_init();

    TEST(test_probe_matches_build_flag);
    TEST(test_cblas_ddot_orthogonal);
    TEST(test_cblas_ddot_parallel);
    TEST(test_cblas_ddot_strided);
    TEST(test_lapack_dgeqrf_link);

    if (g_failures == 0) {
        printf("All LAPACK smoke tests passed.\n");
        return 0;
    } else {
        printf("%d LAPACK smoke test(s) failed.\n", g_failures);
        return 1;
    }
}
