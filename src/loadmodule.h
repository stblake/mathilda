#ifndef LOADMODULE_H
#define LOADMODULE_H

#include <stddef.h>   /* size_t */
#include "expr.h"

/*
 * loadmodule.{c,h} -- runtime loading of Mathilda (.m) source modules.
 *
 * The kernel keeps performance-sensitive primitives in C and higher-level
 * mathematical identities as Mathilda-language rewrite rules in .m files.
 * This module is the bridge: it reads a .m file, parses and evaluates every
 * expression in it, and -- crucially -- locates internal modules robustly no
 * matter what working directory the process was started from (the REPL runs
 * from the repo root; the CMake test binaries run from tests/build/).
 *
 * It generalises the ad-hoc three-path fallback previously hard-coded for the
 * CRC integral tables (see src/calculus/integrate.c) into a single reusable
 * mechanism that both the C side (core_init) and the Mathilda side
 * (LoadModule["..."]) can call.
 */

/*
 * Read, parse, and evaluate every expression in the file at `path`.
 *
 * On success returns the value of the last evaluated expression (caller owns
 * it; it is `Null` for an empty file) and sets *opened to 1.
 * If the file cannot be opened, returns NULL and sets *opened to 0.
 * `opened` may be NULL if the caller does not care to distinguish.
 */
Expr* mathilda_run_file(const char* path, int* opened);

/*
 * Resolve the internal module named by `relpath` (e.g. "simp/FullSimplify.m",
 * relative to the source tree's src/internal directory) to a readable file on
 * disk, writing the located path into `out` (a buffer of `outsz` bytes).
 * Returns 1 and fills `out` if found, 0 (leaving `out` untouched) otherwise.
 *
 * Resolution is independent of the current working directory so a relocated or
 * installed binary still finds its bundled modules. It tries, in order:
 *   1. $MATHILDA_HOME/<relpath>                     (env var -> src/internal)
 *   2. <exe_dir>/src/internal/<relpath>             (binary carries the tree)
 *   3. <exe_dir>/../share/mathilda/internal/<...>   (installed layout)
 *   4. MATHILDA_PREFIX/share/mathilda/internal/<..> (compile-time, if defined)
 *   5. {., .., ../.., ../../..}/src/internal/<...>  (CWD ladder: repo root,
 *                                                    tests/, tests/build/, …)
 * The first candidate that opens wins; its base directory is cached so later
 * lookups skip the search.
 */
int mathilda_resolve_internal(const char* relpath, char* out, size_t outsz);

/*
 * Locate (via mathilda_resolve_internal) and load the internal module named by
 * `relpath` exactly once. Repeat calls with the same `relpath` are no-ops (the
 * module is loaded only once, preventing duplicate rule registration). Returns
 * 1 if the module is loaded (now or previously), 0 if it could not be located
 * (a LoadModule::nofile diagnostic is printed to stderr in that case).
 */
int mathilda_load_module(const char* relpath);

/* Builtin: LoadModule["relpath"] -> True (loaded) or False (not found). */
Expr* builtin_loadmodule(Expr* res);

/* Register the LoadModule builtin and its attributes. */
void loadmodule_init(void);

#endif /* LOADMODULE_H */
