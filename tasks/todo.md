# NLimit — Implementation Plan

`NLimit[expr, z -> z0, opts]` numerically finds the limiting value of `expr` as
`z -> z0` by constructing a sequence of sample points approaching `z0` and
applying sequence-acceleration / extrapolation. Lives in
`src/numerical_calculus/nlimit.c`, mirroring the existing NResidue / ND / NSeries
trio in that directory.

## Research summary (state of the art)

- **Two methods**, faithful to Mathematica's `NumericalCalculus`NLimit`:
  - **`EulerSum`** (default): telescope the approximant sequence into a series
    `S_n = S_0 + Σ (S_k − S_{k−1})` and apply the **Euler transform** (forward
    differences with `1/2^{k+1}` weighting). Best for alternating / slowly
    decaying tails. Depth controlled by `Terms` (default 7).
  - **`SequenceLimit`**: **Wynn's epsilon algorithm** (iterated Shanks
    transform) on the constructed sequence. Recurrence
    `ε_{-1}=0, ε_0=S_n, ε_{k+1}^{(n)} = ε_{k-1}^{(n+1)} + 1/(ε_k^{(n+1)} − ε_k^{(n)})`;
    even columns are the limit estimates. `WynnDegree=n` ⇒ n iterations,
    requiring `≥ 2(n+1)` terms. Exact in one step for geometric/exponential
    tails — "exponential approaches need fewer terms".
- **Sample sequence**: finite `z0`: `z_k = z0 + d·Scale·(1/2)^k`; infinite `z0`:
  `z_k = rayDir·Scale·2^k` (grow outward on a ray from the origin). `d` is the
  `Direction` (Automatic ≡ −1, i.e. approach from larger values; complex `d`
  selects an arbitrary ray).
- **Robustness**: track the largest extrapolation-table magnitude as a
  precision-loss gauge; declare `NLimit::noise` when successive estimates stop
  converging; do NOT auto-Chop residuals (Mathematica leaves them — user applies
  Chop). Pair `WorkingPrecision` with `Terms` (more precision alone doesn't add
  correct digits). Sources: MathWorld (Wynn/Euler), mpmath `shanks`/`richardson`,
  GSL Levin-u, Wolfram NLimit Function Repository.

## Key conventions confirmed from the codebase

- **Closest analog** is `src/numerical_calculus/nderiv.c` (ND) — it already has
  the machine + MPFR two-path structure, `Method->EulerSum`, `Scale`, `Terms`,
  `WorkingPrecision` option handling, Block-style variable binding, and complex
  coercion helpers. Copy its idioms (file-local `nl_*` statics; the `*_bind_*`
  snapshot/set/eval/restore quartet must be replicated since nderiv's are
  `static`).
- **Argument form differs**: NLimit takes a **`Rule[var, z0]`** (`z->z0`), NOT
  the `{z, z0}` list used by NResidue/ND/NSeries. Detect `head == SYM_Rule ||
  SYM_RuleDelayed`, arg_count 2, `args[0]` an `EXPR_SYMBOL`.
- **Option parsing**: mirror `nr_is_known_option` / `nr_is_option_arg` /
  `nr_apply_option` (nresidue.c:226-257). Peel trailing options, then require
  exactly one positional spec arg after the expr.
- **Numeric coercion helpers** to replicate: `*_to_double_real`,
  `*_to_complex`, `*_from_complex_d` (machine); `get_approx_mpfr` +
  `numeric_mpfr_make_complex` (MPFR); `numericalize(expr, spec)`;
  `numeric_digits_to_bits`, `NUMERIC_MACHINE_PRECISION_DIGITS`.
- **Symbols still missing** (must be added to sym_names.{h,c}): `SYM_NLimit`,
  `SYM_Direction`, `SYM_WynnDegree`, `SYM_SequenceLimit`. Already present:
  EulerSum, Scale, Terms, Method, WorkingPrecision, Automatic, Infinity,
  DirectedInfinity, ComplexInfinity, Conjugate, Chop.

## Files

| Action | File |
|--------|------|
| NEW | `src/numerical_calculus/nlimit.h` — `Expr* builtin_nlimit(Expr*); void nlimit_init(void);` |
| NEW | `src/numerical_calculus/nlimit.c` — implementation |
| NEW | `tests/test_nlimit.c` — unit suite |
| EDIT | `src/sym_names.h` / `src/sym_names.c` — add NLimit, Direction, WynnDegree, SequenceLimit |
| EDIT | `src/core.c` — declare + call `nlimit_init()` (next to nseries_init, ~line 550) |
| EDIT | `src/info.c` — terse `NLimit` docstring (no examples) |
| EDIT | `tests/CMakeLists.txt` — add `nlimit.c` to COMMON_SRC + `nlimit_tests` target |
| EDIT | `docs/spec/builtins/calculus.md` — NLimit reference entry |
| EDIT | `docs/spec/changelog/2026-06-08.md` — change summary (Monday of current ISO week) |

## Implementation steps

1. **sym_names**: add the four interned symbols (decl in `.h`, `intern_symbol`
   in `.c`).
2. **nlimit.h**: public prototypes.
3. **nlimit.c — scaffolding**: file-local `NlOpts` struct (method, mode/bits,
   direction `double _Complex` or MPFR pair, scale, terms, wynn_degree); the
   `nl_warn` helper; option predicates/applier; Block-binding quartet; numeric
   coercion helpers (machine + MPFR, all guarded by `#ifdef USE_MPFR`).
4. **nlimit.c — target classification**: detect infinite targets
   (`SYM_Infinity`; `Times[-1, Infinity]`; `DirectedInfinity[dir]`;
   `Times[I, Infinity]` etc.; `ComplexInfinity`) → ray direction; else
   numericalize the finite `z0`. Resolve `Direction->Automatic` to −1 for finite
   points, and to the target's own ray for infinite points.
5. **nlimit.c — sampling**: build `S[0..terms-1]` by evaluating `expr` at the
   geometric sequence; bail with `nnum` if a sample is non-numeric (matches
   Mathematica's Power::infy/indeterminate failures, e.g. `Tanh[x]→∞` at default
   Terms).
6. **nlimit.c — extrapolation cores** (machine + MPFR each):
   - `nl_euler_*`: telescope to differences, Euler transform.
   - `nl_wynn_*`: epsilon table to `WynnDegree` even columns; enforce
     `2(WynnDegree+1) ≤ Terms` (else `NLimit::ndterm` / fall back).
   - guard near-zero epsilon denominators; track max-magnitude gauge; converge
     test on last two estimates → `NLimit::noise` + return NULL when it fails.
7. **nlimit.c — result build**: `nl_from_complex_d` / `numeric_mpfr_make_complex`;
   keep tiny residuals (no auto-Chop). Restore the binding; free scratch.
8. **nlimit.c — `builtin_nlimit`**: glue (validate, parse rule, peel options,
   dispatch machine vs MPFR, choose method). **Ownership**: do not free `res`;
   return new Expr on success / NULL on can't-eval.
9. **nlimit_init**: `symtab_add_builtin("NLimit", builtin_nlimit)`;
   `ATTR_PROTECTED` only (no ATTR_LISTABLE — threading not meaningful here);
   register docstring via `symtab_set_docstring` (or in info.c).
10. **Wire build**: core.c call; CMakeLists COMMON_SRC + test target.
11. **Docs**: calculus.md entry + changelog note.

## Tests (`tests/test_nlimit.c`) — compare in-language, never parse printed reals

- `Sin[x]/x, x->0` → 1
- `(1+1/n)^n, n->Infinity` → E
- `(1+I/x)^x, n->Infinity` → Exp[I] (Cos[1]+I Sin[1])
- `(10^x-1)/x, x->0, Terms->10, Method->SequenceLimit` → Log[10]
- `Tanh[x], x->Infinity, Terms->3` → 1 (and default Terms → fails: `nnum`)
- complex limit point: `Tanh[Pi x]/(1+x^2), x->I` → ≈ −1.5708 I
- `Method->SequenceLimit` with `WynnDegree->3` improves over `->1`
  (e.g. on `SinIntegral` if available, else a constructed geometric sequence)
- `Direction->1` one-sided (use a real one-sided example that needs no Conjugate
  if Conjugate is unreliable; otherwise `z+Conjugate[z]/z, z->0, Direction->1`→1)
- MPFR: `WorkingPrecision->35` on `Sin[x]/x, x->0` → within 35-digit tol of 1
- diagnostic: `1/x, x->0` → `NLimit::noise`, returns unevaluated
- assertion helper: evaluate `Abs[NLimit[...] - expected] < tol` inside Mathilda
- **valgrind**: run `nlimit_tests` under valgrind, diff against `Sin[1.0]`
  baseline noise (macOS dyld/Accelerate ~12.8KB is not a leak).

## Verification

- `make -j` clean under `-std=c99 -Wall -Wextra`.
- Build + run `nlimit_tests` (add to CMakeLists COMMON_SRC or it won't link).
- Run only the scoped `nlimit_tests` binary (not the full suite).
- valgrind clean (modulo known macOS baseline).
- Spot-check each spec example in the REPL.

## Review (completed 2026-06-13)

**Implemented** `NLimit[expr, z -> z0, opts]` in `src/numerical_calculus/nlimit.{c,h}`
(~850 lines), modeled on `nderiv.c`. Machine + MPFR paths for both methods.

- **EulerSum** (default) = Richardson/Romberg extrapolation (denom `2^j-1`,
  consistent with `ND`). Excellent on smooth cases (e.g. `(10^x-1)/x → Log[10]`
  to 1.6e-8 at default Terms; 1e-13 at Terms 10).
- **SequenceLimit** = Wynn epsilon; estimate read from the best-agreeing entry
  of the `ε_{2·WynnDegree}` column (not the roundoff-amplified corner).
  `WynnDegree` scaling verified: 3.1e-5 → 5.5e-7 → 5.6e-9 (deg 1→2→4).
- Sample points: finite `z_k = z0 - d·Scale·2^-k` (the **minus** matches all 3
  documented `Direction` examples); infinite `z_k = u·Scale·2^k` on the target
  ray (handles `Infinity`, `-Infinity`, `I Infinity`, `DirectedInfinity[d]`,
  `ComplexInfinity`, and `Times[..,Infinity]`).
- **Noise gate** measured against the *sample scale* (`max|S_k|`), not the
  result — so genuine near-zero limits aren't flagged, while divergent `1/x`
  is correctly rejected (`NLimit::noise`).

**Wiring:** 4 new symbols (NLimit, Direction, WynnDegree, SequenceLimit) in
`sym_names.{c,h}`; `nlimit_init` in `core.c`; docstring in `info.c`;
`nlimit.c` + `nlimit_tests` in `tests/CMakeLists.txt`; docs in
`docs/spec/builtins/calculus.md` + changelog `2026-06-08.md`.

**Verification:**
- `make -j` clean (`-std=c99 -Wall -Wextra`, no new warnings).
- `nlimit_tests`: all 15 groups pass (finite/infinite/complex limits,
  SequenceLimit + WynnDegree, all Direction cases, Scale/Terms, MPFR to ~1e-27,
  noise diagnostic, unevaluated edge cases, Protected, correctness battery,
  Terms/WynnDegree sweep, memory loop).
- valgrind: only the known macOS baseline (12,800 B / 400 blocks), no `nlimit`
  frames in any leak context.
- Sibling `nderiv`/`nseries`/`nresidue` suites still pass.

**Known limitations (documented):** `Sin[x]/x → Infinity` (oscillatory) returns
a bounded-but-rough value rather than Mathematica's ~3e-7 — genuine Euler
summation would do better there; Richardson does not. Logarithmic approaches
(`Log[x]/x → 0`) converge slowly.
