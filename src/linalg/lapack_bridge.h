/*
 * src/linalg/lapack_bridge.h
 *
 * Direct REPL access to the LAPACK driver routines through the `LAPACK`
 * context: `LAPACK`dgesv`, `LAPACK`dgeqrf`, `LAPACK`dgesdd`, `LAPACK`dgeev`,
 * and their complex `z*` counterparts, using the canonical netlib names.
 *
 * The bindings are ergonomic: dimensions and leading dimensions are inferred
 * from the array shapes. Arrays may be passed as `NDArray[...]` values or as
 * nested `List`s (see numarray.h). Multi-output routines return a `List` of
 * the outputs (e.g. `LAPACK`dgeqrf[A]` -> {Q, R}, `LAPACK`dgesdd[A]` ->
 * {U, S, VT}, `LAPACK`dgeev[A]` -> {values, vectors}).
 *
 * Everything is behind `USE_LAPACK`; the no-LAPACK build provides a stub
 * `lapack_bridge_init()` so the module links either way.
 */
#ifndef MATHILDA_LINALG_LAPACK_BRIDGE_H
#define MATHILDA_LINALG_LAPACK_BRIDGE_H

#ifdef __cplusplus
extern "C" {
#endif

/* Register the `LAPACK`` context builtins, attributes, and docstrings.
 * Called from core_init(). No-op when built without LAPACK. */
void lapack_bridge_init(void);

#ifdef __cplusplus
}
#endif

#endif /* MATHILDA_LINALG_LAPACK_BRIDGE_H */
