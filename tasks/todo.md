# Auto-compile Plot / NIntegrate / FindRoot / Table (M4 wiring)

Wire the numeric compiler as a transparent fast path into the numeric builtins.
Each already binds its variable as an OwnValue and runs full evaluate() per
sample point (thousands of times). Compile the held body ONCE and branch to the
bytecode; fall back to the interpreter when the body is uncompilable, and — for
the integral/root/table callers — per-point when the compiled result is
non-finite (the interpreter may produce a complex/singular value the real
program can't).

## Shared helper — src/compile/autocompile.{h,c}
- `AutoCompiled* autocompile_new(body, const Expr* const* vars, nvars)` → NULL if
  body uncompilable. Interns each var symbol name; compile_expr with all-CT_REAL
  inputs.
- `bool autocompiled_eval_real(ac, const double* xs, double* out)` — false on
  non-finite OR non-real result (caller excludes/falls back).
- `bool autocompiled_eval_complex(ac, const double* xs, double _Complex* out)` —
  handles CT_REAL/INT/COMPLEX results; false on non-finite.
- `autocompiled_free`.

## Callers
- [ ] Plot (plot.c plot_eval_fn) + Plot3D (plot3d.c): compile once in the setup;
      compiled real path, false ⇒ exclude point (matches Plot excluding non-real).
      No per-point fallback (trust the compiler contract; speed is the point).
- [ ] Table (table.c): ONLY the `is_real` (inexact) iterator branch. Compile body
      as f(iterator), CT_REAL. Per element: compiled_eval_real → expr_new_real;
      on false, fall back to evaluate() (complex/singular/Infinity). Exact
      (Integer/Rational/BigInt) iterators UNTOUCHED — exactness preserved.
- [ ] NIntegrate (nint.c ni_eval_at, 1D machine-double): compile body as f(var).
      If the abscissa value is real, autocompiled_eval_complex; box result for
      ni_to_complex. On false OR complex abscissa (contour) ⇒ interpreter. Leave
      ni_sample_mpfr and the multi-D ni_mc_sample untouched (defer multi-D).
- [ ] FindRoot (findroot.c, FR_PREC_MACHINE only): compile f (and symbolic df
      when available) as f(var). Use in fr_eval_with_bindings-equivalent hot path;
      fall back on non-finite. Leave MPFR drivers untouched.

## Guardrails
- Never touch MPFR paths.
- Table: never compile an exact iterator.
- Parity tests: compiled vs interpreter to machine precision for each caller.
- Each caller commits separately with its own verification + micro-benchmark.

## Deferred
- Multi-D NIntegrate (ni_mc_sample), FindRoot systems, complex-contour NIntegrate
  compilation, MPFR fast paths.

## Review — DONE (2026-07-27)

Shared adapter `src/compile/autocompile.{c,h}` (real+complex eval, interpreter
fallback, self-contained header). Wired into all four builtins, each its own
commit with parity/fallback tests + micro-benchmark, all leaks-clean, MPFR
untouched:
- Plot/Plot3D (27f5ccd): ~215× / ~11×; non-real sample excludes the point.
- Table (4595cf7): machine-real iterator only; exact iterators bit-for-bit
  unchanged; ~128×.
- NIntegrate (ee4eb3a): 1-D machine ni_eval_at; complex per-sample fallback;
  ~353×. Fixed ni_ctx_rebody bug (context copy carried stale compiled program).
- FindRoot (c1edf83): scalar machine real via pointer-identity guard in
  fr_eval_with_bindings; ~19×.
tests/test_autocompile.c covers all four end-to-end. Docs: changelog
2026-07-27.md, design compile.md status, control-flow.md Compile section.

Next: M3 arrays/NDArray; or multi-D NIntegrate / FindRoot systems.
