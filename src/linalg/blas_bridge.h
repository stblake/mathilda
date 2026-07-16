/*
 * src/linalg/blas_bridge.h
 *
 * Direct REPL access to the BLAS (Basic Linear Algebra Subprograms) kernels
 * through the `BLAS`` context: `BLAS`ddot`, `BLAS`dgemv`, `BLAS`dgemm`, and
 * their complex `z*` counterparts, using the canonical netlib routine names.
 *
 * The bindings are ergonomic: dimensions and leading dimensions are inferred
 * from the array shapes, so only the mathematically meaningful arguments
 * (`alpha`, `beta`, and the operand arrays) are exposed. Arrays may be passed
 * as `NDArray[...]` values or as ordinary nested `List`s (see numarray.h).
 *
 * Everything is behind `USE_LAPACK` (the same flag that gates the CBLAS/LAPACK
 * link); the no-BLAS build provides a stub `blas_bridge_init()` so the module
 * links either way, matching the flint_bridge.c graceful-degrade policy.
 */
#ifndef MATHILDA_LINALG_BLAS_BRIDGE_H
#define MATHILDA_LINALG_BLAS_BRIDGE_H

#ifdef __cplusplus
extern "C" {
#endif

/* Register the `BLAS`` context builtins, their attributes, and docstrings.
 * Called from core_init(). No-op when built without BLAS/LAPACK. */
void blas_bridge_init(void);

#ifdef __cplusplus
}
#endif

#endif /* MATHILDA_LINALG_BLAS_BRIDGE_H */
