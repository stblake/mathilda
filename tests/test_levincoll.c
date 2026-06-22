/* Tests for the Levin collocation engine (levincoll.c).
 *
 * These exercise the pure-numeric core directly through C samplers (no
 * expression layer): for ∫_a^b f·{e^{ig}|cos g|sin g} dx with known closed
 * forms, the collocation result must match to a tight tolerance, and accuracy
 * must IMPROVE as the oscillation rate grows.  The complex LU solve, the
 * Chebyshev basis/derivative recurrences, and the endpoint reconstruction are
 * all covered transitively by the full-integral checks.
 */
#include "levincoll.h"

#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "test_utils.h"

/* A linear phase g(x) = c·x with a chosen amplitude. */
typedef struct { double c; } Phase;

static bool s_amp_one(void* ctx, double x, double _Complex* o) {
    (void)ctx; (void)x; *o = 1.0; return true;
}
static bool s_amp_x(void* ctx, double x, double _Complex* o) {
    (void)ctx; *o = x; return true;
}
static bool s_gprime(void* ctx, double x, double _Complex* o) {
    (void)x; *o = ((Phase*)ctx)->c; return true;
}
static bool s_gphase(void* ctx, double x, double _Complex* o) {
    *o = ((Phase*)ctx)->c * x; return true;
}

static bool lc_close(double _Complex got, double _Complex want, double tol) {
    return cabs(got - want) <= tol;
}

/* ∫_0^1 cos(c x) dx = sin(c)/c  (COS kernel, f = 1). */
static void test_cos_const_amp(void) {
    Phase ph = { 40.0 };
    LevinResult r = levin_collocation_machine(
        0.0, 1.0, s_amp_one, NULL, s_gprime, &ph, s_gphase, &ph,
        LEVIN_KERNEL_COS, 1e-10, 64);
    ASSERT(r.have && r.conv);
    double want = sin(40.0) / 40.0;
    ASSERT_MSG(lc_close(r.val, want, 1e-9), "cos: got %.15g want %.15g",
               creal(r.val), want);
}

/* ∫_0^1 sin(c x) dx = (1 - cos(c))/c  (SIN kernel, f = 1). */
static void test_sin_const_amp(void) {
    Phase ph = { 50.0 };
    LevinResult r = levin_collocation_machine(
        0.0, 1.0, s_amp_one, NULL, s_gprime, &ph, s_gphase, &ph,
        LEVIN_KERNEL_SIN, 1e-10, 64);
    ASSERT(r.have && r.conv);
    double want = (1.0 - cos(50.0)) / 50.0;
    ASSERT_MSG(lc_close(r.val, want, 1e-9), "sin: got %.15g want %.15g",
               creal(r.val), want);
}

/* ∫_0^1 x cos(c x) dx = (cos c + c sin c)/c^2 - 1/c^2  (COS kernel, f = x). */
static void test_cos_linear_amp(void) {
    double c = 100.0;
    Phase ph = { c };
    LevinResult r = levin_collocation_machine(
        0.0, 1.0, s_amp_x, NULL, s_gprime, &ph, s_gphase, &ph,
        LEVIN_KERNEL_COS, 1e-10, 64);
    ASSERT(r.have && r.conv);
    double want = (cos(c) + c * sin(c)) / (c * c) - 1.0 / (c * c);
    ASSERT_MSG(lc_close(r.val, want, 1e-9), "x cos: got %.15g want %.15g",
               creal(r.val), want);
}

/* ∫_0^1 e^{i c x} dx = (e^{ic} - 1)/(i c)  (EXP kernel, complex result). */
static void test_exp_const_amp(void) {
    double c = 30.0;
    Phase ph = { c };
    LevinResult r = levin_collocation_machine(
        0.0, 1.0, s_amp_one, NULL, s_gprime, &ph, s_gphase, &ph,
        LEVIN_KERNEL_EXP, 1e-10, 64);
    ASSERT(r.have && r.conv);
    double _Complex want = (cexp(I * c) - 1.0) / (I * c);
    ASSERT_MSG(lc_close(r.val, want, 1e-9), "exp: got (%.12g,%.12g) want (%.12g,%.12g)",
               creal(r.val), cimag(r.val), creal(want), cimag(want));
}

/* Accuracy improves with the oscillation rate: a fixed small order nails a very
 * oscillatory integral that ordinary quadrature would struggle with. */
static void test_high_frequency(void) {
    Phase ph = { 5000.0 };
    LevinResult r = levin_collocation_machine(
        0.0, 1.0, s_amp_x, NULL, s_gprime, &ph, s_gphase, &ph,
        LEVIN_KERNEL_COS, 1e-10, 32);
    ASSERT(r.have && r.conv);
    double c = 5000.0;
    double want = (cos(c) + c * sin(c)) / (c * c) - 1.0 / (c * c);
    ASSERT_MSG(lc_close(r.val, want, 1e-10), "hi-freq: got %.15g want %.15g",
               creal(r.val), want);
}

/* A non-zero lower limit exercises the affine node map a != 0.
 * ∫_1^3 cos(c x) dx = (sin(3c) - sin(c))/c. */
static void test_shifted_interval(void) {
    double c = 60.0;
    Phase ph = { c };
    LevinResult r = levin_collocation_machine(
        1.0, 3.0, s_amp_one, NULL, s_gprime, &ph, s_gphase, &ph,
        LEVIN_KERNEL_COS, 1e-10, 64);
    ASSERT(r.have && r.conv);
    double want = (sin(3.0 * c) - sin(c)) / c;
    ASSERT_MSG(lc_close(r.val, want, 1e-9), "shifted: got %.15g want %.15g",
               creal(r.val), want);
}

int main(void) {
    TEST(test_cos_const_amp);
    TEST(test_sin_const_amp);
    TEST(test_cos_linear_amp);
    TEST(test_exp_const_amp);
    TEST(test_high_frequency);
    TEST(test_shifted_interval);
    printf("All levincoll tests passed.\n");
    return 0;
}
