# Mathilda Stress-Test & Code Review

**Date:** 2026-06-25
**Reviewer role:** Expert software tester / C auditor
**Scope:** (1) REPL-level behavioural stress testing, (2) source-code audit, (3) documentation review.
**Build:** `make` clean, `./Mathilda` (commit `6b8fd0f`).

This report is organised by phase. Each finding has a severity
(Critical / High / Medium / Low) and, where possible, a verified
reproduction or a `file:line` citation. Findings marked **[REPRO]** were
reproduced live at the REPL; findings marked **[CODE]** were verified by
reading source but not (yet) reproduced end-to-end.

> **Update (2026-06-25): fixed so far** — A (Integrate Listable → definite
> integrals now return unevaluated), E (LCM/GCD int64 overflow), F
> (rbd_factor_mpz unbounded loop), G (mpoly_normalize leak), H (HypergeometricPFQ
> MPFR leak), I (C99 `_GNU_SOURCE`/`strdup`/`__attribute__` cleanup). See the
> 2026-06-22 changelog. Remaining: B, C, D and the consistency/doc items.

---

## Phase 1 — REPL Stress Testing

### Summary

The parser and evaluator are notably **robust against malformed input** —
no segfaults from unbalanced brackets, empty input, junk operators, or
moderately deep nesting; errors are reported cleanly. However, one
**severe functional defect** (definite integration is broken) and several
genuine defects plus a cluster of WL-consistency gaps were found.

### 1.0 — Definite integration `Integrate[f,{x,a,b}]` is completely broken — **Critical [REPRO]**

```
Integrate[x^2,{x,0,1}]   →  {1/3 x^3, Integrate[x^2, 0], Integrate[x^2, 1]}    (should be 1/3)
Integrate[Sin[x],{x,0,Pi}] →  {-Cos[x], Integrate[Sin[x], 0], Pi Sin[x]}        (should be 2)
Integrate[x^2,{x,y,z}]   →  {1/3 x^3, x^2 y, x^2 z}                             (proves the cause)
```

`Integrate` is registered with `ATTR_LISTABLE` (`src/calculus/integrate.c:446`),
so the range spec `{x,a,b}` is threaded element-wise as an ordinary list:
the evaluator computes `{Integrate[f,x], Integrate[f,a], Integrate[f,b]}`
instead of evaluating the antiderivative at the bounds. **No definite
integral works.** The docstring only advertises `Integrate[f, x]`
(indefinite), so definite integration may simply be unimplemented — but the
spurious `ATTR_LISTABLE` turns the standard syntax into silent garbage rather
than a clean "unevaluated". **Fix:** remove `ATTR_LISTABLE` from `Integrate`
and implement the `{x,a,b}` range (antiderivative + bound substitution, with
`Integrate::idiv` for divergent cases like `Integrate[1/x,{x,-1,1}]`).

### 1.1 — `[[ ... ]]` with >64 indices aborts (stack-smashing) — **High [REPRO]**

```
x[[1,1,...,1]]   (200 indices)   →   abort, exit code 134
```

The `Part` index parser uses a fixed `Expr* args[64]` stack buffer with an
unbounded append loop (`src/parse.c:~1143`), so ≥64 indices smash the stack.
The stack-protector catches it (`rc=134`), but it is a hard crash on
attacker/user-controllable input. `parse_list`/`parse_function` already use
growable buffers; this path was missed. **Cross-referenced by the core-engine
audit (Phase 2).**

### 1.2 — Recursion-limit unwind leaks internal `Hold[...]` artifacts — **Medium [REPRO]**

```
g[0]:=1; g[n_]:=n*g[n-1]; g[10000]
→  ...big number... (8982 + Hold[Times][Hold[-1], Hold[1]]) g[8982 + Hold[Times][Hold[-1], Hold[1]] - 1]

{{{...5000 deep...}}}  →  ... Hold[List][Hold[{{...}}]] ...
```

When `$RecursionLimit` (1024) is hit, the partially-unwound expression
exposes internal `Hold[Head]` / `Hold[args]` wrappers that are an
implementation detail and must never be user-visible. WL returns a clean
held form. Closely related: `f[n_]:=f[n-1]; f[5]` returns
`f[-8185 - 1]` — the inner `-8185 - 1` is left unevaluated, so the
returned tree is malformed/non-canonical.

### 1.3 — Many structural ops silently return unevaluated with no message — **Medium [REPRO]**

WL emits a diagnostic and (mostly) returns the input unchanged; Mathilda
returns it unchanged but **emits nothing**, so a user cannot tell a real
error from an intentional no-op. Observed:

| Input | Mathilda | WL behaviour |
|-------|----------|--------------|
| `{1,2,3}[[10]]` | `Part[{1,2,3},10]` silently | `Part::partw` message |
| `{1,2,3}[[-5]]` | unevaluated silently | `Part::partw` |
| `{1,2,3}[[1.5]]` | unevaluated silently | `Part::pspec` |
| `First[{}]` / `Rest[{}]` | unevaluated silently | `First::nofirst` / `Rest::norest` |
| `Take[{1,2,3},10]` | unevaluated silently | `Take::take` |
| `StringTake["abc",10]` | unevaluated silently | `StringTake::take` |
| `D[x^2, 5]` | returns `0` silently | `General::ivar` (5 not a valid variable) |

(Note: `Prime[0]` stays unevaluated because **`Prime` is not a registered
builtin at all** — see Phase 3 doc finding #7 — not merely a missing message.)

This mirrors the project's own `*Q`-predicate rule ("NULL return is only OK
with an emitted diagnostic"). The structural builtins violate the spirit of
that rule. `{1,2,3}[[0]]` correctly returns the head `List`.

### 1.4 — In/Out line counter increments on parse error — **Low [REPRO]**

```
1+1      → Out[1]= 2
(((      → Parse error
2+2      → Out[3]= 4      (Out[2] is skipped)
```

A syntax error consumes a line index, leaving gaps in `In[n]`/`Out[n]`. In
WL a syntax error does not advance the line number. Cosmetic but affects
back-reference (`%2`) reliability.

### 1.5 — Minor WL-semantics gaps (Low, [REPRO])

- `DirectedInfinity[0]` stays unevaluated; WL → `ComplexInfinity`.
- `Quotient[5,0]`, `Mod[5,0]` stay unevaluated (defensible, but WL has
  defined behaviour for `Quotient`).
- `0x1F` parses as `0 * x1F → 0` (no hex literal support); `1.5e` parses as
  `1.5 * e`. Silent surprises rather than bugs.
- `Set::wrsym` is emitted for `Sin=5` (good — assignment correctly blocked),
  but WL's tag is `Set::write`. Message-name cosmetic.

### 1.6 — Linear-algebra / calculus REPL gaps (Low–Medium, [REPRO])

- **`Eigenvalues` of a non-square matrix returns unevaluated with no message**
  (`Eigenvalues[{{1,2,3}}]`), whereas `Inverse`/`Det`/`Dot`/`LinearSolve` all
  emit proper diagnostics (`Inverse::matsq`, `Det::matsq`, `Dot::dotsh`,
  `LinearSolve::nosol`). Inconsistent.
- **Diagnostic message formatting is inconsistent:** `Inverse::matsq` prints
  the pretty form `{{1,2,3},{4,5,6}}`, but `Det::matsq` and `Dot::dotsh` print
  raw FullForm `List[List[1,2,3],...]`. Should use one form (pretty).
- `Transpose[{{1,2},{3}}]` (ragged) returns unevaluated silently; WL →
  `Transpose::nmtx`.
- Linear algebra error handling is otherwise **excellent** — singular,
  non-square, dimension-mismatch, and inconsistent-system cases all caught.

### Things that worked well (no defect)

- Unbalanced `(`, `{`, `[`; empty input; `]]]`; `1++2`; deep parens (5000)
  — all handled without crash.
- `0^0`→Indeterminate, `1/0`→ComplexInfinity (+`Power::infy`),
  `(-1)!`→ComplexInfinity, `Sqrt[-1]`→`I`, `Log[0]`→`-Infinity`,
  `2^1000000`, `PartitionsP[100000]`, `1000!/999!`→1000 — all correct.
- `Solve[x==x,x]`→`{{}}`, `Solve[1==2,x]`→`{}`, `FactorInteger[0]`→`{{0,1}}`,
  `Binomial[-1,2]`→1, `Divisors[-12]` → positive divisors — all WL-correct.

---

## Phase 2 — Source-Code Audit

Five focused auditors swept the source tree against the project's memory-safety
contract (builtin ownership of `res`; GMP/MPFR init/clear pairing), strict-C99
compliance, integer overflow, and bounds. Highlights are promoted here; the
full per-area findings follow.

### Top-priority defects (cross-phase)

| # | Severity | Location | Defect |
|---|----------|----------|--------|
| A | **Critical** | `integrate.c:446` | `Integrate` is `ATTR_LISTABLE` → **definite integration `Integrate[f,{x,a,b}]` is completely broken** (threads over the range list) **[REPRO]** |
| B | **Critical** | `limit.c:2740-2812` | Heap **and** stack buffer overflow in multivariate `Limit` for n≥4 (fixed `[3]` strides indexed to `n`) |
| C | **Critical** | `facint.c:1294/1300` | Use-after-free of `method_opt` in rational `FactorInteger` |
| D | **High** | `parse.c:~1143` | Fixed 64-slot stack buffer in `[[…]]` → overflow/abort **[REPRO §1.1]** |
| E | **High** | `numbertheory.c:257,173` | `LCM`/`GCD` int64 fast path overflows silently **[REPRO]** |
| F | **High** | `facint.c:773-833` | `rbd_factor_mpz` unbounded loop → REPL hang |
| G | **High** | `mpoly.c:234-266` | `mpoly_normalize` leaks `mpz_t` slots on every term-combine (hot path) |
| H | **High** | `hypergeopfq.c:416-420` | MPFR `cpx_t` leak on partial parameter-extraction failure |
| I | **High** | `modular.c`, `trig.c` | `_GNU_SOURCE` / `_USE_MATH_DEFINES` + `strdup` — banned C99 violations |
| J | **Medium** | `print.c:729,748` | `open_memstream` + reassigning `stdout` — non-C99/non-portable |
| K | **Medium** | multiple | `strdup` unguarded (expr/symtab/match.c); `__attribute__((unused))` unguarded (intrischnorman/intrat) |

**[REPRO] for A:** `Integrate[x^2,{x,0,1}]` → `{1/3 x^3, Integrate[x^2,0], Integrate[x^2,1]}`
(should be `1/3`); `Integrate[x^2,{x,y,z}]` → `{1/3 x^3, x^2 y, x^2 z}` proves the
Listable threading. Even trivial definite integrals fail.
**[REPRO] for E:** `LCM[4000000007, 4000000009]` → `-2446744009709551553`
(wrapped negative) instead of `16000000064000000063`.


### 2A — Core Engine (expr / parse / eval / symtab / match / replace / attr / print)

**Memory safety**
- **`match.c:436-441` — Medium.** Named-`Pattern` sub-match failure path does not `env_rollback(env, saved_env_count)` before `return false`; sub-bindings leak and pollute the env for later attempts. Every sibling branch rolls back. Fix: add the rollback.
- **`parse.c:~1143` — High.** `OP_PART` parses `[[…]]` indices into `Expr* args[64]` with an unbounded append → stack buffer overflow (**[REPRO §1.1]**). Convert to a growable buffer like `parse_function`.
- **`parse.c:615,654` — Medium.** `elements = realloc(elements, …)` / `args = realloc(args, …)` with **no NULL check** → leak + NULL-deref under OOM.
- **`eval.c:1125-1129` — Low.** Unchecked `malloc` in the recursion-overflow wrap path (runs exactly when resources are scarce).
- Pervasive unchecked `malloc`/`realloc` in match.c/replace.c/symtab.c/eval.c (codebase convention; crash-on-OOM).

**C99 violations**
- **`expr.c:70,284`, `symtab.c:95`, `match.c:80` — Medium.** `strdup` unguarded (papered over with `#define _GNU_SOURCE`), exactly the pattern SPEC §10 forbids. Add a C99 `mathilda_strdup`.
- **`print.c:729,748` — Medium.** `expr_to_string` relies on `open_memstream` (POSIX.1-2008, not C99) **and reassigns the `stdout` global** (`stdout = stream;`) — not guaranteed a modifiable lvalue, not thread-safe. Build into a growable buffer via `vsnprintf` sizing.

**Overflow / logic**
- **`match.c:165-178` — Low.** `(int)spec->data.integer` truncates `Repeated[p, 5e9]` bounds.
- **`replace.c:377-388` — Low.** `All`/`Infinity` levels encoded as magic `1000000`; depth >10^6 silently stops replacing.
- **`parse.c:884-893` — Low.** Missing `)` only warns to stderr but still returns `inner` as success → structurally wrong tree inside larger expressions (relates to §1.1 lenient error handling).
- **`parse.c:570-589 / 101-115` — Low.** `parse_string`/`parse_symbol` use fixed `char buffer[256]`; over-long strings/identifiers silently truncate (two distinct long names can collide to one interned symbol), and string escapes are not decoded (`\n` → literal `n`).
- **`eval.c:1138-1202` — Low.** Fixed-point cap `MAX_ITERATIONS` (4096) prints "$IterationLimit exceeded" but is a hardcoded constant, not the user-settable `$IterationLimit` (inconsistent with the configurable `$RecursionLimit`).
- **`attr.c:143-163` — Low (perf).** `get_attributes` does an O(n) `strcmp` scan of the ~120-entry table on *every* function node, redundant with the symtab attribute copy; the rest of the engine uses interned-pointer compares.

**Verified non-issues:** ownership contract respected in audited builtins; `expr_copy` refcount-bump is by design; `expr_compare` intentionally in sort.c.

### 2B — Arithmetic & Number Theory (plus / times / power / arithmetic / numbertheory / rat / parfrac / expand / modular / facint)

**Critical**
- **`facint.c:1294/1300`.** Use-after-free: `expr_free(method_opt)` at 1294, then `method_opt` passed to `factorize_mpz(...)` at 1300 when denominator > 1. Triggers on `FactorInteger[1/n, Method->{"BlakeRationalBaseDescent", ...}]`. Move the free after the denominator block.

**High**
- **`modular.c` — `_GNU_SOURCE`+`strdup`** (lines 1-3, 208, 291). Banned. Replace with `malloc+strcpy`.
- **`numbertheory.c:257` (LCM) & `:173` (GCD).** int64 fast path `(a/gcd)*b` overflows silently; bigint path only triggers when an arg is *already* bigint. **[REPRO]** `LCM[4000000007, 4000000009]` → `-2446744009709551553`. Use `__builtin_mul_overflow`/GMP.
- **`facint.c:773-833` — `rbd_factor_mpz`** `for(j=3;;j++)` has no iteration cap; the `"BlakeRationalBaseDescent"` method hangs the REPL on a hard composite. Add a bound.

**Medium**
- **`arithmetic.c:221` — `builtin_divide`** computes `make_rational(n1*d2, d1*n2)` in plain int64 (no `__int128` guard, unlike plus.c/times.c) → silent wrong rational on overflow.
- **`times.c:82` — mixed-BigInt fast path** keeps `int64_t den = d1*d2` while the numerator uses GMP → wrapped denominator seeds the GCD/divide.
- **`power.c:1097` — unguarded `head->data.symbol`** read without `EXPR_SYMBOL` type check (UB for composite-headed base); every sibling check guards first.
- **`numbertheory.c:32,169,250` — `llabs(INT64_MIN)`** UB in abs/GCD/LCM int64 paths (also `arithmetic.c` gcd/lcm). Route INT64_MIN through GMP.

**Low**
- `power.c:947/925/927/907`, `expand.c:318/360` — negating a possibly-`INT64_MIN` exponent (UB, practically unreachable).
- `modular.c:206-208` — fixed 256-byte `snprintf("%s$%lld")` can truncate long names → temp/real symbol collision.
- `facint.c:494-538` — `cfrac_factor_mpz` `while(1)` has no hard iteration cap (terminates in practice).
- `parfrac.c:245` — unguarded `head->data.symbol` read (same pattern as power.c:1097; currently unreachable).
- **`numbertheory.c` — missing docstrings.** Only `FactorialPower` gets `symtab_set_docstring`; `GCD, LCM, ExtendedGCD, PowerMod, PrimitiveRoot, MultiplicativeOrder, Factorial, Binomial, Divisible, CoprimeQ, Divisors` lack one — violates SPEC §10.
- Pervasive unchecked `malloc`/`calloc`.

**Verified clean:** plus.c overflow handling (proper `__int128` + BigInt promotion); power.c edge cases (0^0, x^0, 0^−n, 0^+n) WL-faithful; rat.c/parfrac.c are int64-only (no GMP leak surface); division/modulus-by-zero guarded across PowerMod/Divisible/MultiplicativeOrder/divide.


### 2C — Math & Special Functions (trig / hyperbolic / logexp / complex / special_functions / calculus)

Ownership contract is clean throughout (no `expr_free(res)` violations, no double-frees); MPFR/GMP pairing clean on normal paths.

**Critical**
- **`limit.c:2740,2744,2757,2760` — heap buffer overflow.** `int (*dirs)[3] = calloc(max_dirs, sizeof(*dirs))` (3-int rows) but filled `dirs[nd][j]` for `j<n`; the `n>=4` branch (line 2987) reaches it, so `Limit[f,{w,x,y,z}->{0,0,0,0}]` writes past each row. Fix: stride `n` — `calloc(max_dirs*n, sizeof(int))`, index `dirs[nd*n+j]`.
- **`limit.c:2767,2784,2810-2812` — stack buffer overflow.** `Expr* values[3]` indexed `values[j]` for `j<n` (write/read/free) — overruns for n≥4 with `Expr*` pointers. Same entry path. Fix: heap-size `values` from `n`. **[REPRO note]** live `Limit[w+x+y+z,{w,x,y,z}->{0,0,0,0}]` returned *no output at all* (no result, no error) — consistent with silent corruption; would be a hard crash under ASAN.

**High**
- **`hypergeopfq.c:416-420` — MPFR leak on partial parameter-extraction failure.** The `!ok` branch only `free(ac)/free(bc)`; the already-init'd `ac[0..i)`/`bc[0..j)` `cpx_t` pairs are never `cpx_clear`'d. Fix: clear all prior elements first.

**Medium**
- **`intrischnorman.c:796,1256,1264` & `intrat.c:979` — unguarded `__attribute__((unused))`** — forbidden by SPEC §10; the same file uses `(void)x;` elsewhere. Replace with `(void)` casts.
- **`intrischnorman.c:3270,3279,3307,1906` — unchecked matrix-sized `calloc`/`malloc`** dereferenced in the solver path (NULL-deref on OOM).
- **`harmonicnumber.c:179,199` — unchecked `malloc(sizeof(Expr*)*n)`** (n up to 100000) then dereferenced; inconsistent with the rest of the module.
- **`qpochhammer.c:61-64` — uncapped `n` + unchecked malloc** in the 3-arg finite product: `QPochhammer[a,q,1000000000]` attempts a billion-element alloc/`Times` (Hyperfactorial/BarnesG cap; this doesn't).
- **`complex.c:648,654` — lossy `(double) ==` in `builtin_arg` exact-quadrant path** for Rational components; two distinct rationals rounding to one double wrongly yield exact `Pi/4`. Fix: compare via `expr_eq(re,im)`.

**Low**
- **`trig.c:1-3` — `#define _USE_MATH_DEFINES` and `#define _GNU_SOURCE`** (both forbidden; the `#ifndef M_PI` guard already makes M_PI safe — delete both).
- **`erf.c:329, erfc.c:301, erfi.c:351, expintegralei.c:457` — no double-overflow promotion** in the machine-complex (`prec<=53`) path: a large finite MPFR value silently becomes `±Inf`. Airy/Bessel guard this. Mirror `ai_complex_result`.
- **`fibonacci.c:178, lucas.c:186` — no degree cap on polynomial forms** (`Fibonacci[1000000,x]` runs ~1M Expand iterations); siblings cap (`LEG_POLY_CAP`, etc.).
- **`productlog.c:436` — tautological range check** (`k<LONG_MIN||k>LONG_MAX` always false on LP64 — dead).
- **`zeta.c:431-433` — empty `if` block** (dead code, harmless).

**Verified clean:** gamma/loggamma/polygamma/hurwitzzeta/lerchphi/polylog/beta/pochhammer, bessel/airy/inverf/logintegral, legendre, and the integrator cascade (`integrate*.c`, `intsimp.c`) — ownership honored, recursion depth-capped, no non-terminating loops.

### 2D — Polynomials, Linear Algebra & Simplification (poly/ linalg/ simp/ part / sort / funcprog / purefunc / patterns / stats / iter / piecewise)

Ownership contract clean (no `expr_free(res)` violations). Note: there is no `src/list.c` (list ops live in `src/list/`, out of this agent's named scope — recommend a follow-up sweep there).

**High**
- **`mpoly.c:234-266` — `mpoly_normalize` leaks `mpz_t` slots on every term-combine (hot path).** `n_terms`/`cap` are overwritten with `out` before the trim loop, making the trailing `mpz_clear` loop dead and the trim guard `out < out` always false. Any cancel/combine in `mpoly_mul`/`mpoly_add`/substitution permanently leaks `new_coefs[out..N-1]`. Fix: capture `old_n` before the overwrite, clear `[out..old_n)`, then realloc to `out`.

**Medium**
- **`piecewise.c:123-133` — `Floor`/`Ceiling`/`Round`/`IntegerPart` of large `EXPR_REAL` truncates via `(int64_t)`** (UB for `Floor[1.0e30]`). The MPFR branch routes through `mpz_t` correctly; the plain-double branch doesn't. Fix: route `|res|>=9e18` through `mpfr_get_z`.
- **`poly/poly.c:4374-4397 (apply_floor_to_coeffs)` — unguarded `head->data.symbol` read + raw `.data.integer`** on Rational/Complex components (wrong field for bigint components).
- **`poly/minpoly.c:253-271` — Rational/Complex handled with no `arg_count==2` check** → reads past `args` on under-arity input.
- **`poly/subresultants.c:179-186` & `subresultantpoly.c:175-182` — unchecked `realloc`** (leak + NULL-deref on OOM; also `qaupoly.c:90`).
- **`simp/simp_cuberoot.c:260-264` — signed int64 overflow in the cube-root denesting grid search** → can yield a mathematically *wrong* simplification. Fix: compute in `mpz_t` or gate on magnitude.
- **`stats.c` exact-rational paths (`builtin_mean` 127-134, `builtin_variance` 338-370) — all-int64 rational arithmetic overflows silently** for large values/long lists. Fix: use `mpq_t`.

**Low**
- `linalg/construct.c:33-47` — `IdentityMatrix` truncates int64 dim to `int`, no size cap (signed UB for `m>INT_MAX`); `DiagonalMatrix` already uses int64 counters.
- `part.c:55,97,217,335,583,718` — `idx->...->head->data.symbol` read without `EXPR_SYMBOL` guard (inconsistent with the guarded `part.c:124`).
- `part.c:569-570 (compare_paths)` — reads `args[i]->data.integer` assuming integer path elements (no type check in the qsort comparator).
- `sort.c:440` — `builtin_orderedq` returns `expr_new_symbol("True"/"False")` literals instead of interned `SYM_True/False`.
- `patterns.c:70,383`, `funcprog.c:616,1017,1171` — `realloc(*results,…)` assigned directly (OOM leak/NULL-deref); pervasive unchecked-alloc class.
- `poly/qa.c:63,164,241`, `minpoly.c:257,269` — `mpq_set_si(...,(unsigned long)den)` corrupts a negative denominator (latent; all callers pass positive).
- `poly/qa.c:104` — `qa_alpha` returns 0 (not α) for degree-1 extension (degenerate, unreachable).
- `linalg/dot.c:135-136` — `Dot` in-place compaction skips re-checking index 0 after a merge (logic-only; masked by fixed-point loop).

**Verified clean:** `iter.c`/`purefunc.c` (per-iteration frees balanced, `ITER_SAFETY_CAP` bounds `Do`/`NestWhile`/`FixedPoint`); `funcprog.c` ownership (`Distribute`/`Thread`/`Tuples`/`Outer`); all `simp/` drivers have explicit round caps (no non-terminating rewrites); `groebner*.c` have step caps + overflow-checked weights.

---

## Phase 3 — Documentation Review

SPEC.md has drifted far from the code; the docs/spec builtin reference has a few real gaps. The good news: **SPEC §4's builtin ownership contract is correct** (evaluator frees `res`; the contradicting auto-memory note is itself stale), and all cross-reference links resolve.

**High**
1. **SPEC.md §1 — size grossly understated.** Claims "~22 kLoC across ~44 source modules"; actual is **271 `.c` files / ~195,000 lines** (excl. `src/external/`). README.md ("159,000 lines / 176 modules / ~365 builtins") is also stale. Fix: update to ~195 kLoC / ~271 files / 438 registered builtins and reconcile README.
2. **SPEC.md §1 & §9 — MPFR omitted.** MPFR is default-on (`USE_MPFR ?= 1`, `-lmpfr`, used in 128 files, powers `N[expr,prec]`/Precision/numeric_calculus) but is absent from §1's dependency list and §9's link line (`-lreadline -lgmp -lm`). Fix: add MPFR.
3. **SPEC.md §5 — module init list incomplete + names a nonexistent init.** Lists ~28 `*_init`; `core.c` invokes ~120. Names `load` init, but the real symbol is `loadmodule_init()` (no `load_init`). Fix: regenerate/representative-list it, correct `load`→`loadmodule`.

**Medium**
4. **SPEC.md §2 — layout lists `load.c` which does not exist** (it's `loadmodule.c`).
5. **SPEC.md §2 — omits four populated subdirs:** `src/numerical_calculus/`, `src/numerical_roots/`, `src/product/`, `src/sum/` (internally inconsistent with Mathilda_spec.md, which documents their builtins).
6. **SPEC.md §2 — `numeric.c` and `sym_names.c` missing from layout** (both exist; `sym_names.c` is the CLAUDE-mandated home for internal symbols).
7. **docs/spec — `Prime` and `Exponent` documented but not callable.** Used in examples (`calculus.md:936`, `number-theory.md:200`; `algebra.md:156,189-197`) but neither is registered nor defined in any `.m`. (`PrimeQ`/`PrimePi`/`NextPrime`, `IntegerExponent`/`RealExponent` exist; bare forms don't.) Cross-confirmed at the REPL (`Prime[0]` stays unevaluated). Fix: implement, or rewrite examples.

**Low**
8. **16 registered builtins are undocumented:** `AddTo, Subtract, SubtractFrom, Increment, Decrement, PreIncrement, PreDecrement, Divide, Factorial2, Begin, End, EndPackage, SetAccuracy, SetPrecision, ReplaceList, Quartiles`. Add entries.

**Verified correct:** SPEC §4 ownership contract matches `eval.c:884-886`; all Mathilda_spec.md / docs/extending.md cross-links resolve; the changelog table (10 weekly files through 2026-06-22) is complete; §9 test-binary naming matches the CMake convention.

---

## Appendix — Methodology & Caveats

- **Build:** `make` (already current); `./Mathilda` (commit `6b8fd0f`, `-std=c99 -O3`, GMP+MPFR+readline). The binary is a **normal release build** — no ASAN/UBSan — so several memory-corruption findings (B, the `[[…]]` overflow, the limit overflow) manifest as silent corruption / `__stack_chk_fail` abort rather than clean diagnostics. **Recommend re-running the test suite under `-fsanitize=address,undefined` to surface these deterministically.**
- **Source audit** was performed by five focused agents (core; arithmetic/number-theory; math/special; poly/linalg/simp; documentation). Every reported `file:line` was read by an agent; items not reproduced end-to-end are unmarked, REPL-reproduced items are tagged **[REPRO]**.
- **Coverage gaps / suggested follow-ups:** `src/list/` (not in any agent's named scope), `src/graphics/`, `src/repl.c`, the `.m` bootstrap files, and `tests/` were not deeply audited. A valgrind leak run over a representative session and a full `*_tests` pass under sanitizers would complement this review.

### Recommended fix order

1. **A — `Integrate` `ATTR_LISTABLE`** (one-line attribute fix; restores a core CAS feature).
2. **B, C, D — memory corruption** (`limit.c` overflows, `[[…]]` buffer, `facint.c` UAF) — crash/exploitable.
3. **E, G, H — silent wrong answers / hot-path leaks** (`LCM` overflow, `mpoly_normalize` leak, hypergeo leak).
4. **F, I — non-terminating loops & C99 violations** (`rbd_factor_mpz` cap, `_GNU_SOURCE`/`strdup`/`__attribute__` cleanup).
5. **Consistency pass** — emit WL diagnostics for out-of-range structural ops (§1.3, §1.6); recursion-unwind `Hold[]` leak (§1.2); SPEC.md refresh (Phase 3).
