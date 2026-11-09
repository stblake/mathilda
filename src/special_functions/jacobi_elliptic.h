/* Mathilda -- the Jacobi elliptic functions JacobiSN ... JacobiSD (parameter m
 * convention), JacobiAmplitude, and InverseJacobiSN / CN / DN. See
 * jacobi_elliptic.c. */
#ifndef MATHILDA_JACOBI_ELLIPTIC_H
#define MATHILDA_JACOBI_ELLIPTIC_H

#include <stdbool.h>
#include "expr.h"

/* Function codes: the twelve pq = p/q (s, c, d, n with n = 1), then the
 * amplitude, then the three inverses. */
enum {
    JAC_SN, JAC_CN, JAC_DN, JAC_CD, JAC_CS, JAC_DC, JAC_DS,
    JAC_NC, JAC_ND, JAC_NS, JAC_SC, JAC_SD,
    JAC_AM, JAC_ISN, JAC_ICN, JAC_IDN, JAC_COUNT
};

/* The machine kernel for function `which` at (u, m), real or complex parts. The
 * same helper behind the scalar builtins, the NDArray / packed map and the
 * Compile[] VM. Declines only for a pole, a value outside the double range, or
 * (inverses) a complex value the interpreter's EllipticF form answers. */
bool sf_machine_jacobi(int which, double ur, double ui, double mr, double mi,
                       double* rr, double* ri);

void jacobi_elliptic_init(void);

#endif /* MATHILDA_JACOBI_ELLIPTIC_H */
