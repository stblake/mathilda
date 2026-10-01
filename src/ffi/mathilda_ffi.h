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

/* Parse and evaluate `input` ONCE, returning a single JSON object describing
 * the result — the mobile counterpart to the sidecar's per-line NDJSON (see
 * repl.c: pipe_process_input), minus the "id"/"done" framing the host adds.
 * This is the preferred entry point for the notebook UI: a lone eval keeps
 * side effects (assignments, counters, RandomReal) from running twice, and it
 * carries the plot payload that _eval / _eval_latex cannot.
 *
 * The returned string is one of:
 *   {"type":"error","message":"Parse error: ..."}       parse failure
 *   {"type":"plot","payload":<plotly-json>}              Graphics[...] / Graphics3D[...]
 *   {"type":"expr","payload":"<text>","latex":"<latex>"} normal result (latex omitted if empty)
 *   {"type":"expr","payload":""}                         evaluated to nothing
 *
 * Caller-owned memory; free with mathilda_ffi_free(). Never returns NULL. */
char* mathilda_ffi_eval_json(const char* input);

/* --- Notebook-cell evaluation and editor services ------------------------
 *
 * These exist for a front end that drives the kernel as a notebook or console
 * (e.g. a xeus Jupyter kernel): a cell may hold several statements, Print
 * output and messages must reach the cell, and an editor wants completion and
 * a continuation check. They are the in-process counterpart of the sidecar's
 * `cell:true` NDJSON pipe mode (repl.c: pipe_process_input), and emit the SAME
 * protocol event lines, so one event mapping serves both transports. */

/* Receives one NDJSON event line (no trailing newline), exactly as the sidecar
 * would write it to stdout: {"id":N,"type":"line"|"stream"|"message"|"expr"|
 * "plot"|"image"|"usage"|"names"|"error"|"done", ...}. The embedder parses it
 * and publishes the matching front-end message. Called in order, on the
 * calling thread, before mathilda_ffi_eval_cell returns. */
typedef void (*mathilda_ffi_sink)(void* ctx, const char* json_line);

/* Evaluate `code` for request `id`, delivering each event to `sink`.
 *
 * cell != 0 — notebook-cell semantics: the cell is split into statements like a
 *   Mathematica input cell (newline- or ';'-separated); the whole cell is
 *   parsed before anything runs, so a syntax error anywhere evaluates nothing;
 *   each statement's Print output arrives as a "stream" event and its messages
 *   as "message" events BEFORE its result; a statement ending in ';' or
 *   evaluating to Null sends no result; a per-statement "line" event carries
 *   the session $Line (so `%`, In[n], Out[n] work across cells).
 * cell == 0 — one expression, no history, no capture (the batch-tool path).
 *
 * A terminating {"type":"done"} is always the last event. Rich results route
 * as in mathilda_ffi_eval_json (Graphics -> "plot", Image -> "image", ?name ->
 * "usage"/"names", otherwise "expr" with an optional "latex"). Implicitly
 * initializes the kernel; not reentrant (serialize calls — see file header). */
void mathilda_ffi_eval_cell(const char* code, int id, int cell,
                            mathilda_ffi_sink sink, void* ctx);

/* 1 if `code` forms a complete expression (every (), [], {} closed, no open
 * "..." string or (* ... *) comment), 0 otherwise. Drives a console/notebook
 * continuation prompt — the same rule the REPL's smart-Return uses. */
int mathilda_ffi_is_complete(const char* code);

/* A JSON array (as a string) of defined symbol names beginning with `prefix`,
 * sorted, for tab completion — e.g. prefix "Sin" -> "[\"Sin\",\"Sinh\",
 * \"SinIntegral\"]". Caller-owned; free with mathilda_ffi_free(). Never NULL
 * ("[]" when nothing matches). A NULL/empty prefix lists every symbol. */
char* mathilda_ffi_complete(const char* prefix);

/* Release a string returned by mathilda_ffi_eval / _eval_latex. NULL-safe. */
void mathilda_ffi_free(char* s);

/* Version string of the linked kernel (static storage; do not free). */
const char* mathilda_ffi_version(void);

#ifdef __cplusplus
}
#endif

#endif /* MATHILDA_FFI_H */
