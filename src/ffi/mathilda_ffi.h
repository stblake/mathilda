/* mathilda_ffi.h — in-process C API for embedding the Mathilda kernel.
 *
 * The desktop notebook talks to a spawned `mathilda` sidecar over stdio. On
 * mobile platforms (iOS / Android) spawning child processes is forbidden by
 * the OS sandbox, so the kernel must run *in process*. This header exposes a
 * tiny, stable C ABI that a host (the Rust/Tauri app, a test harness, any
 * FFI caller) links against directly.
 *
 * Threading: the Mathilda evaluator uses global state (the symbol table), so
 * these functions are NOT reentrant. Callers must serialize access — e.g.
 * behind a mutex. Initialize once, evaluate on a single logical thread.
 */
#ifndef MATHILDA_FFI_H
#define MATHILDA_FFI_H

#ifdef __cplusplus
extern "C" {
#endif

/* Initialize the kernel: symbol table, builtins, and the internal `init.m`
 * bootstrap. Idempotent — the second and later calls are no-ops. Safe to call
 * before the first `mathilda_ffi_eval`.
 *
 * The `init.m` tree is located via mathilda_resolve_internal(); on a bundled
 * app set the MATHILDA_HOME environment variable to the directory containing
 * the `internal/` tree before calling this (see mathilda_ffi_set_home). */
void mathilda_ffi_init(void);

/* Point the loader at the bundled `internal/` module tree. Equivalent to
 * setting the MATHILDA_HOME environment variable, but callable from hosts that
 * cannot easily set the process environment before init. `dir` is copied.
 * Call BEFORE mathilda_ffi_init(). No-op if `dir` is NULL/empty. */
void mathilda_ffi_set_home(const char* dir);

/* Parse, evaluate to a fixed point, and format `input` as Mathilda output
 * (the same text the REPL prints after `Out[n]= `). Returns a newly allocated,
 * NUL-terminated UTF-8 string owned by the CALLER — release it with
 * mathilda_ffi_free(). Never returns NULL: parse failures yield the string
 * "$Failed (parse error)" and expressions that evaluate to nothing yield "".
 *
 * Calls mathilda_ffi_init() implicitly if the kernel is not yet initialized. */
char* mathilda_ffi_eval(const char* input);

/* Like mathilda_ffi_eval, but formats the result as a LaTeX string suitable
 * for a math renderer (KaTeX/MathJax). Returns caller-owned memory; free with
 * mathilda_ffi_free(). Never returns NULL. */
char* mathilda_ffi_eval_latex(const char* input);

/* Release a string returned by mathilda_ffi_eval / _eval_latex. NULL-safe. */
void mathilda_ffi_free(char* s);

/* Version string of the linked kernel (static storage; do not free). */
const char* mathilda_ffi_version(void);

#ifdef __cplusplus
}
#endif

#endif /* MATHILDA_FFI_H */
