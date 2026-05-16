# Integrate (Rational Function Integration) — Implementation Plan

> **Note:** Plan mode restricts edits to this file. The final committed plan should be copied to `INTEGRATE_PLAN.md` in the repo root after approval.

---

## Context

Mathilda currently has no `Integrate` builtin. The user wants a faithful C99 port of `IntegrateRational.m` (1989 lines of Mathematica), which implements the modern textbook pipeline for rational function integration:

1. **Mack's linear Hermite Reduction** — splits `f = A/D` into `D[g, x] + h` where `h` has a square-free denominator (Bronstein, *Symbolic Integration I*, p. 44).
2. **Lazard–Rioboo–Trager log-part algorithm** — computes the logarithmic part of the integral using a subresultant pseudo-remainder sequence (PRS) over `K(x)[t]` (Bronstein p. 51, Trager 1976).
3. **Rioboo's `LogToReal` / `LogToAtan`** — converts sums of complex logarithms back to real `Log + ArcTan` form when the resultant `Q(t)` factors over a tractable algebraic extension (Bronstein p. 63).
4. **Pattern-based `LogToArcTan` / `LogToArcTanh`** — final simplification combining `Log[a]−Log[b]` into `ArcTanh[(b−a)/(b+a)]`-style closed forms.

References (in order of importance to the implementation):
- Bronstein, *Symbolic Integration I — Transcendental Functions*, Springer 2nd ed., 2004 — §2.2-2.8.
- Trager, *Algebraic Factoring and Rational Function Integration*, ACM SYMSAC 1976. (PDF in repo root.)
- Geddes et al., *Algorithms for Computer Algebra*, 1992 — Example 11.12 used for the `lc(d)·polymod` corrective shift inside `IntRationalLogPart`.

The Mathilda baseline has surprisingly strong infrastructure already (subresultant-PRS-based multivariate `PolynomialGCD`, `PolynomialExtendedGCD`, `Resultant`, `FactorSquareFree`, `Apart`, `Together`, `Cancel`, `D`, and Trager-algorithm algebraic factoring via `Factor[poly, Extension -> α]`). The main gaps relative to Mathematica are: `Solve`, `RootSum`, `RootReduce`, `ToRadicals`, `Apart` does not accept `Extension`, and `PolynomialQuotientRemainder` is not exposed (only the separated pair).

**Scope decisions confirmed by user:**
- Cover the full pipeline including `LogToArcTan` / `LogToArcTanh`, recreating the missing pieces (Solve-quadratic, basic RootReduce, etc.) where needed.
- User-visible names: `Integrate[f, x]` in `System\`` (auto-dispatch when input is a rational function in `x`); `Integrate\`BronsteinRational[f, x]` is the explicit form (matches the Mathematica file). `Integrate\`` becomes the umbrella context; `Integrate[f, x]` calls it internally.
- Helpers (`HermiteReduce`, `IntRationalLogPart`, `LogToAtan`, `LogToReal`, `LogToArcTan`, `LogToArcTanh`, `IntegratePolynomial`) are exposed under the `Integrate\`` context so each stage can be unit-tested at the REPL.

**Out of scope (explicitly deferred):**
- Inexact / machine-precision integrands (Mathematica's `Rationalize` step) — return unevaluated.
- `Refine`-style assumption tracking on parameters — instead we use `Together`+`Cancel` and trust user input.
- `ConditionalExpression` post-checking — emit final answer without it.
- `Solve` over reals beyond cases reducible to degree ≤ 2, biquadratic (`a t^4 + b t^2 + c`), or `n`-th-root (`a t^n + b`) factors. Higher-degree real conversion returns a symbolic `RootSum` form (added as a passive head).

---

## Architecture

### New files
| File | Purpose | Approx. LOC |
|------|---------|-------------|
| `src/intrat.c` | Full algorithm implementation | ~1500-1800 |
| `src/intrat.h` | Public + helper-builtin entry points | ~80 |
| `src/integrate.c` | Tiny dispatch shim (`Integrate[]` → `Integrate\`BronsteinRational`) | ~120 |
| `src/integrate.h` | Public entry point | ~20 |
| `tests/test_intrat.c` | ~60 unit tests covering each phase | ~600 |

### Modified files
| File | Change |
|------|--------|
| `src/core.c` | Call `integrate_init()` (after `expand_init()`, before `deriv_init()`) |
| `src/sym_names.c` / `src/sym_names.h` | Cache `SYM_Integrate`, `SYM_BronsteinRational`, `SYM_RootSum`, `SYM_HermiteReduce` and a few internals |
| `src/poly.c` / `src/poly.h` | Add `SubresultantPolynomialRemainders[a, b, x]` builtin (the algorithm needs the *full chain*, not just the resultant) and `PolynomialQuotientRemainder[a, b, x]` (returns `{q, r}` pair) |
| `tests/CMakeLists.txt` | Register `test_intrat` |
| `Mathilda_spec.md` | Document new builtins per CLAUDE.md instructions |

### Dependency graph for `intrat.c`
```
intrat.c
├── poly.c           PolynomialGCD, PolynomialExtendedGCD, PolynomialQuotient,
│                    PolynomialRemainder, Resultant, CoefficientList,
│                    Coefficient, Variables, Collect, [new] SubresultantPRS,
│                    [new] PolynomialQuotientRemainder
├── facpoly.c        FactorSquareFree, builtin_factor (with Extension)
├── qafactor.c       qa_factor_with_extension, qaupoly_norm (used to drive
│                    Lazard-Rioboo-Trager when input has algebraic coefficients)
├── parfrac.c        Apart (rational-coeff fast path), [extended] Apart-with-Extension
├── rat.c            Together, Cancel
├── expand.c         Expand, expr_expand
├── deriv.c          builtin_d (used by HermiteReduce: D[d, x])
└── core/arith       PolynomialQ, integer arithmetic, GMP via expr_to_mpz
```

---

## Phased Implementation

Each phase ends with green tests in `tests/test_intrat.c`. Per CLAUDE.md ("If tests pass, move to the next phase immediately") we proceed continuously after each phase.

### Phase 0 — `Extension` option for core polynomial routines (prerequisite)
**Goal:** `Together`, `Cancel`, `PolynomialGCD`, `PolynomialLCM`, and `Apart` all accept an `Extension -> α` option, defaulting to `Extension -> None` (current behavior). This is required because `BronsteinRational` calls `Cancel[Together[..., Extension -> α], Extension -> α]` as its `canonic` step on every iteration — without it, every algebraic-coefficient test case fails before the integrator even starts.

**Why this comes first:** the integrator's correctness depends on these primitives being closed over `Q(α)`. We can't unit-test Phase 1's `HermiteReduce` on inputs with `Sqrt[2]` coefficients otherwise.

**Files modified:**
- `src/poly.c` / `src/poly.h` — `PolynomialGCD`, `PolynomialLCM` accept trailing `Rule[Extension, α]`.
- `src/rat.c` / `src/rat.h` — `Together`, `Cancel` accept the option and propagate it to the internal GCD call.
- `src/parfrac.c` / `src/parfrac.h` — `Apart` (option already mentioned in the original Phase 2 plumbing; bring it forward to Phase 0 for symmetry).
- `src/sym_names.c` / `.h` — `SYM_Extension` and `SYM_None` already exist; reuse.
- `Mathilda_spec.md` — extend the docstring for each builtin.
- `tests/test_extension_options.c` (NEW) — option-parsing + correctness coverage.
- `tests/CMakeLists.txt` — register the new test binary.

**Implementation strategy** (avoids touching `qa.h` / `qafactor.c` internals):
1. Add a shared static helper `expr_extract_extension_option(Expr *args[], size_t *argc, Expr **alpha_out)` (place in `src/options.c`/`.h` if multiple modules end up needing it; otherwise duplicate the existing pattern from `facpoly.c:2782-2800`). Strips a trailing `Rule[Extension, α]` from the argument list, normalizes `Extension -> None` and `Extension -> Automatic` (treat Automatic as None for now — auto-detection is deferred).
2. **`PolynomialGCD[a, b, ..., Extension -> α]`**: when an extension is present, route through factorization rather than the Z[x]-PRS path. For a pair `(a, b)`:
   - `fa = Factor[a, Extension -> α]`, `fb = Factor[b, Extension -> α]` (these already work and dispatch through `qaupoly_gcd` internally — `qafactor.h:154`).
   - Decompose each into a `(factor, multiplicity)` multiset.
   - Intersect the multisets with min-multiplicity, multiply back, expand. Handles the multi-arg form by left-fold.
   - Constant content is treated as a degree-0 factor so `PolynomialGCD[2 (x^2-2), 4 (x-Sqrt[2]), Extension -> Sqrt[2]] == 2 (x - Sqrt[2])` works.
   - Restriction: requires univariate input when the extension is present (raise `PolynomialGCD::ext` and fall back to the no-extension path on multivariate inputs).
3. **`PolynomialLCM[a, b, ..., Extension -> α]`**: `lcm(a, b) = Expand[a * b / gcd(a, b)]` using the extended GCD from step 2.
4. **`Cancel[expr, Extension -> α]`**: extract `n = Numerator[expr]`, `d = Denominator[expr]`, compute `g = PolynomialGCD[n, d, Extension -> α]`, return `PolynomialQuotient[n, g, x] / PolynomialQuotient[d, g, x]` where `x` is auto-detected via `Variables`.
5. **`Together[expr, Extension -> α]`**: existing `Together` logic combines into a single fraction (no extension needed for this step), then apply step-4 `Cancel` with the supplied extension to the result.

**Attributes:** unchanged — these all stay `Protected`. The option just enters via trailing `Rule[Extension, _]` arguments, parsed by the builtin itself (matching the precedent set by `Factor`).

**Docstrings (`symtab_set_docstring`):** every modified builtin's docstring gets a new line:
```
  Option Extension -> α (default None) factors over Q(α) when computing
  the polynomial GCD.  Extension -> {α_1, ..., α_n} builds the compositum.
```

**Tests** (`tests/test_extension_options.c`, ~25 cases):
- *Parsing:* `Together[1/x + 1/x, Extension -> Sqrt[2]]` is recognised and dispatches; `Together[1/x + 1/x, Extension -> None]` matches the no-option output bit-for-bit; `Together[1/x + 1/x, Foo -> Bar]` is rejected (returns unevaluated, mirrors `Factor`).
- *Correctness:*
  - `PolynomialGCD[x^2 - 2, x - Sqrt[2], Extension -> Sqrt[2]] === x - Sqrt[2]`.
  - `PolynomialGCD[x^2 - 2, x - Sqrt[2]] === 1` (default).
  - `PolynomialLCM[x^2 - 2, x - Sqrt[2], Extension -> Sqrt[2]] === x^2 - 2` (up to sign).
  - `Cancel[(x^2 - 2)/(x - Sqrt[2]), Extension -> Sqrt[2]] === x + Sqrt[2]`.
  - `Cancel[(x^2 - 2)/(x - Sqrt[2])] === (x^2 - 2)/(x - Sqrt[2])` (default unchanged).
  - `Together[1/(x - Sqrt[2]) + 1/(x + Sqrt[2]), Extension -> Sqrt[2]] === 2x/(x^2 - 2)`.
  - Tower form: `Cancel[(x^2 - 6)/(x - Sqrt[2] Sqrt[3]), Extension -> {Sqrt[2], Sqrt[3]}] === x + Sqrt[6]`.
- *Backwards compatibility:* run a representative subset of the existing `test_polygcd`, `test_apart`, `test_rat` tests through this code path to confirm zero regression on `Extension -> None` / no-option.

**Phase 0 exit criteria:**
- All five builtins parse the option without error.
- `tests/test_extension_options.c` is green.
- `make` clean, `valgrind --leak-check=full ./test_extension_options` shows zero leaks.
- No regressions in `test_polygcd`, `test_apart`, `test_rat`, `test_qafactor`.

---

### Phase 1 — Scaffolding, polynomial integration, Hermite reduction
**Goal:** `Integrate[poly, x]` and `Integrate[A/D, x]` for any squarefree `D` (covers cases without repeated roots).

**New builtins** (all in `Integrate\`` context except the System dispatcher):
- `Integrate[f, x]` — System dispatcher. Validates `PolynomialQ[f,x] || rationalQ[f,x]`; calls `Integrate\`BronsteinRational[f, x]`; otherwise returns unevaluated.
- `Integrate\`BronsteinRational[f, x]` — top-level. Calls the helpers below.
- `Integrate\`IntegratePolynomial[f, x]` — `Expand` then term-by-term: `a x^n -> a x^(n+1)/(n+1)`, `1/x -> Log[x]`, constant -> `c x`.
- `Integrate\`HermiteReduce[f, x]` — Mack's linear version (Bronstein p. 44). Returns `{g, h}` such that `f = D[g,x] + h` and `Denominator[h]` is squarefree.

**Pre-Hermite fast path — derivative recognition:**
Before invoking `HermiteReduce`, check whether the integrand has the form `c · D'/D^k` for some constant `c` and integer `k ≥ 1`:
- Let `n = Numerator[f]`, `d = Denominator[f]`. Compute `dprime = D[d_base, x]` where `d = d_base^k`.
- If `Cancel[Together[n - c · dprime · d_base^(k-1)]] === 0` for some constant `c` (extract `c` by polynomial division), short-circuit:
  - `k == 1` → emit `c · Log[d_base]`.
  - `k > 1` → emit `−c/((k−1) · d_base^(k−1))`.
- Catches `(2x)/(x²+1)`, `1/(x · Log[x])`-style integrands and the entire numerator-is-derivative-of-denominator family without ever entering the LRT machinery. ~30-50 LOC, expected to short-circuit ~10-15% of real-world inputs.

**Debug/trace scaffolding (used by every later phase):**
- Global `Integrate\`$Verbose = False` (OwnValue, default off). When `True`, every helper prints an `IN <args>` line on entry and `OUT <result>` on exit, mirroring the Mathematica baseline's `debugprint`.
- Single C helper `intrat_trace(const char *fn, const char *direction, Expr *e)` — gated on the `$Verbose` lookup. Cheap when off (single symbol-table lookup); invaluable when debugging Phase 5+.
- Add the entry/exit calls to every public helper from Phase 1 onward — retrofitting later means re-touching every file.
- `Integrate\`Helpers\`Content[p, x]`, `Integrate\`Helpers\`Primitive[p, x]`, `Integrate\`Helpers\`Monic[p, x]`, `Integrate\`Helpers\`LeadingCoefficient[p, x]` — internal utilities. Exposed for testing only.
- New System builtin: `PolynomialQuotientRemainder[a, b, x]` — returns `{q, r}`. Trivial wrapper over the existing `poly_div_rem` C helper (already used by `PolynomialQuotient` / `PolynomialRemainder` separately).

**Algorithms to implement:**
- `content[p, x]` = `PolynomialGCD @@ CoefficientList[p, x]`.
- `primitive[p, x]` = `p / content[p, x]`, expanded.
- `monic[p, x]` = divide all coefficients by `lc[p, x]`.
- `lc[p, x]` = `Coefficient[p, x, Exponent[p, x]]`.
- `exquo[a, b, x]` = `PolynomialQuotient[a, b, x]` after asserting remainder is zero (raise message + return unevaluated quotient on failure).
- `gcd[a, b, x]` = first element of `PolynomialExtendedGCD[a, b, x]`.
- `ExtendedEuclidean[a, b, c, x]` (Diophantine) — exactly mirrors the Mathematica function in `IntegrateRational.m:1203-1210`, using `PolynomialExtendedGCD` then `PolynomialQuotientRemainder` of `q*t` mod `b`.
- `HermiteReduce[f, x]` — direct port of `IntegrateRational.m:1303-1323`.

**Attributes:** All new builtins get `ATTR_PROTECTED | ATTR_READPROTECTED`. Set docstrings via `symtab_set_docstring`.

**Tests:** ~12 cases — pure polynomial, simple `1/(x-a)`, `1/(x-a)^2`, the four `HermiteReduce` cases lifted from `IntegrateRational.m:1331-1358`. Verification: `D[Integrate[f,x], x] - f // Together // Cancel === 0`.

### Phase 2 — Lazard-Rioboo-Trager log part (structural)
**Goal:** Compute the logarithmic part of `∫ A/D dx` where `D` is squarefree, returning a symbolic `RootSum` for each squarefree factor.

**New builtins:**
- `Integrate\`IntRationalLogPart[A/D, x, t]` — direct port of `IntegrateRational.m:751-801`. Returns either a list of `{Q_i[t], S_i[t,x]}` pairs (when `RootSum -> False`) or a sum of `RootSum[Function[t, Q_i[t]], Function[t, t Log[S_i[t,x]]]]` terms.
- `Integrate\`Helpers\`SquareFree[p]` — Mathilda wrapper around `FactorSquareFree` that buckets factors of equal multiplicity (mirrors `IntegrateRational.m:1477`). Returns the list `{ {m_i, factor_of_multiplicity_i} ... }`.
- `Integrate\`Helpers\`ApartList[expr, x, Extension -> opt]` — returns `Apart`'s output already split into a list (mirrors `IntegrateRational.m:955`).
- `Integrate\`Helpers\`ExtractConstants[f, x]` — pulls scalar prefactors out of numerator/denominator so the LRT subresultant chain doesn't blow up; mirrors `IntegrateRational.m:1013-1023`.
- New System builtin: `SubresultantPolynomialRemainders[a, b, x]` — exposes the multivariate subresultant PRS already used internally by `PolynomialGCD` (see `src/poly.c:943-1411`). Returns the chain `{S_n, S_{n-1}, ..., S_0}` as a list. Without this, the algorithm cannot extract the degree-`i` element of the chain (see `Extract[prs, Position[degs, i] // First]` in the Mathematica source).
- New passive head: `RootSum` — registered with no builtin (just a symbol so it can appear in output). Set `ATTR_PROTECTED`.

**Algorithm:**
- Direct port of `IntegrateRational.m:751-801`. The Geddes-Czapor-Labahn corrective branch (lines 764-768) handles the case where the squarefree exponent equals `Exponent[d, x]`.
- `apartList` over Q via existing `Apart`. When `Extension -> Sqrt[c]` is supplied, factor the denominator with `qa_factor_with_extension` first, then call `Apart` on the rewritten expression.

**Cyclotomic / `n`-th-root structural pre-detection:**
Before launching the (expensive) subresultant PRS, inspect the squarefree denominator `D`:
- If `D` matches `a · x^n + b` with `a, b` constant (the `nthrootQ` predicate, `IntegrateRational.m:1093`): immediately split `D` into linear and quadratic real factors using `n`-th roots of `(−b/a)`, then recurse `IntRationalLogPart` on each piece. The `Q[t]` resultant is then linear or quadratic by construction and Phase 4's `LogToReal` handles it cleanly.
- If `D` is cyclotomic (`x^n − 1`, `x^n + 1`, or more generally a product of cyclotomic polynomials): use the explicit cyclotomic factorization (`Φ_d(x)` for `d | n`) before PRS.
- Catches the `1/(x^8 + 1)`, `1/(x^6 − 1)`, `x/(x^8 + 1)` family from `IntegrateRationalTests.m` without running the full LRT pipeline. The Mathematica baseline only does this inside `LogToReal`; promoting it to `IntRationalLogPart` is faster and produces cleaner output. ~80-120 LOC.

**Tests:** ~10 cases targeting just `IntRationalLogPart` (Bronstein examples 2.4.1, 2.5; tests file lines 822-882). Verification: `D[result, x] - integrand` simplifies to zero.

### Phase 3 — `LogToAtan` (Rioboo recursive)
**Goal:** Given `A, B ∈ K[x]`, return a sum of `ArcTan` of polynomials such that the derivative equals `D[I·Log[(A+I·B)/(A−I·B)], x]`. This is the inner workhorse.

**New builtin:**
- `Integrate\`LogToAtan[a, b, x]` — direct port of `IntegrateRational.m:1529-1561`. Recursive: when `Exponent[A,x] < Exponent[B,x]`, recurse with `LogToAtan[-B, A, x]`. Otherwise compute `g, {d, c} = PolynomialExtendedGCD[B, -A, x]`, set `v = (A d + B c)/g`, and return `2 ArcTan[v] + LogToAtan[d, c, x]`.

Uses `PolynomialExtendedGCD`, `PolynomialQuotient` (for `exquo`), and a `polymod`/`PolynomialRemainder` helper (already exposed). All required pieces exist.

**Tests:** ~6 cases lifted from `IntegrateRational.m:1574-1597`.

### Phase 4 — `LogToReal` (Rioboo) with bounded `Solve`
**Goal:** Convert `RootSum[Q[t]==0, α Log[S[α, x]]]` into real `Log + ArcTan` form when `Q[t]` factors over Q (or a supplied algebraic extension) into pieces of degree ≤ 2 / biquadratic / `n`-th-root.

**New builtins:**
- `Integrate\`LogToReal[Q, S, x, t]` — port of `IntegrateRational.m:1618-1658`. Algorithmic structure:
  1. Factor `Q[t]` using `Factor[Q, Extension -> opt]`.
  2. For each linear factor `t − r` (with `r ∈ K`): emit `r Log[S /. t -> r]`.
  3. For each quadratic factor `t² + p t + q`:
     - Compute `disc = p² − 4q`. If `disc` is a non-negative rational, return the two real roots `(−p ± Sqrt[disc])/2` and treat as two linear cases.
     - If `disc < 0`, write the complex conjugate pair `u ± I v` with `u = −p/2`, `v = Sqrt[−disc]/2`. Compute `{A, B} = ReIm[S /. t -> u + I v]` via `ComplexExpand`-style splitting (we'll port `polynomialIm` / `polynomialRe` from `IntegrateRational.m:1827-1847`). Emit `u Log[A² + B²] + v · LogToAtan[A, B, x]`.
     - If `disc` is symbolic and we cannot decide the sign, fall back: emit a `Piecewise` with both branches OR (simpler) leave as `RootSum` and let the user supply assumptions.
  4. For biquadratic factors (`a t^(2m) + b t^m + c`): substitute `t^m -> s`, recurse on the quadratic, back-substitute `m`-th roots of the result. Mirrors `biquadraticQ` at line 1060.
  5. For `n`-th-root factors (`a t^n + b`): roots are `(−b/a)^(1/n) · ω_n^k` for `k = 0…n−1`. Emit accordingly. Mirrors `nthrootQ` at line 1093.
  6. Otherwise: return the original `RootSum` unchanged.
- `Integrate\`Helpers\`PolynomialReal[p, x]` and `Integrate\`Helpers\`PolynomialImag[p, x]` — split a polynomial whose coefficients live in `Q[I, ...]` into real and imaginary parts via coefficient-wise `ComplexExpand`. Tiny helpers.
- `Integrate\`Helpers\`ComplexExpand[e]` — sufficient subset of Mathematica's: distribute `Re[]` / `Im[]` over `Plus`, handle `Re[a + b I] = a`, `Im[a + b I] = b`, `Re[I^n]`, `Re[a*b]`, `Re[a^n]` for integer `n`. Several hundred lines but a clean, well-bounded port.

**Tests:** ~10 cases lifted from `IntegrateRationalTests.m` covering the Bronstein examples and the `1/(x^4+1)` family.

### Phase 5 — `IntegrateRealRationalFunction` glue
Port `IntegrateRational.m:534-587`. Composes Phase 1-4: Hermite-reduce, polynomial-quotient-remainder split, partial-fraction expansion, per-term `IntRationalLogPart` followed by `LogToReal`, with `NaiveLogPart` fallback when `LogToReal` returns the symbolic form for an unresolvable factor.

**New builtin:**
- `Integrate\`IntegrateRealRationalFunction[f, x]` — top-level glue.
- `Integrate\`NaiveLogPart[f, x]` — port of `IntegrateRational.m:1116-1124`. Used as fallback: emit `RootSum[d, a Log[T-x]/d'] /. T -> x`.

### Phase 6 — `LogToArcTan` and `LogToArcTanh` post-processing
**Goal:** Pattern-based simplification combining `c·Log[A] + c·Log[B] -> c·Log[A B]`, `c·Log[A] − c·Log[B] -> c·Log[A/B]`, and `c·Log[A] − c·Log[B] -> 2c·ArcTanh[(A−B)/(A+B)]` when the resulting argument is rational in `x`.

**Implementation choice:** Rather than pattern-rewriting at the Mathilda-rule layer (slow and fragile), implement these as direct C transformations operating on `Plus[...]` of `Log[...]` terms.

**New builtins:**
- `Integrate\`LogToArcTan[e, x]` — port of `IntegrateRational.m:1896-1958`.
- `Integrate\`LogToArcTanh[e, x]` — port of `IntegrateRational.m:1722-1761`.
- `Integrate\`Helpers\`ZeroQ[e]` — `TrueQ[Cancel[Together[e]] === 0]`.

These are pure post-processing on the Phase 5 output; if disabled (option `LogToArcTan -> False`) the result is still correct.

### Phase 7 — `BronsteinRational` top-level wrapper + Options + output cleanup
Port `IntegrateRational.m:67-121`. Implements:
- `Options[BronsteinRational] = {"PFD" -> True, Extension -> Automatic, "LogToArcTan" -> True, "Radicals" -> False}`.
- Stitches Phase 1-6 together: canonicalize, integrate, optionally `LogToArcTan` + `LogToArcTanh`, optionally `ToRadicals` (we'll skip this — leave as identity until `ToRadicals` exists).
- The bare `Integrate[f, x]` shim (in `integrate.c`) just calls `Integrate\`BronsteinRational[f, x]` and returns its result, falling through to unevaluated if the input is non-rational.

**Final-pass output cleanup (applied as the last step before returning):**

1. **`Collect` by transcendental head.** Apply `Collect[result // Expand, {_Log, _ArcTan, _ArcTanh}, simproot]` (mirrors `IntegrateRational.m:99` and `:578`). This merges duplicated transcendental terms — without it, the output of LRT often contains five separate copies of `Sqrt[2]` floating around as coefficients of the same `Log`. The Mathematica baseline does this and it's responsible for most of the "looks-like-a-textbook-answer" output quality.
   - `simproot` here is our local helper: `Cancel[Together[#]]` followed by a future `RootReduce` slot when implemented (no-op for now).
   - Implement as a dedicated `intrat_collect_transcendentals(Expr *e, const char *var)` C function rather than a generic `Collect` invocation, so we can canonicalize coefficients in-place.

2. **Output canonicalization for `Log` / `ArcTan` / `ArcTanh` arguments.** A single normalization pass that runs after `Collect`:
   - **`Log[arg]`**: if the leading coefficient (in `var`) of `arg` is negative, rewrite as `Log[−arg]` and absorb the implicit `+ I π` into a constant (which then drops, since the integration constant is implicit). For `Log[c · p]` with `c` free of `var`, drop `c` (already done by `IntegrateRational.m:1746`); make this systematic.
   - **`ArcTan[arg]`**: normalize so the leading coefficient of `arg` (or of the numerator if `arg` is a quotient) is positive — using `ArcTan[−x] = −ArcTan[x]`. Halves the equivalence class of "structurally different but mathematically identical" outputs.
   - **`ArcTanh[arg]`**: same treatment as `ArcTan` (odd function).
   - The point of canonicalization is reproducible test output: structurally-equivalent results compare equal under string match. Without it, swapping the order of two roots of a quadratic gives a different-looking (but correct) answer, which makes test assertions brittle.

3. **Sort the top-level `Plus`.** Apply Mathilda's existing `expr_compare`-based canonical ordering (already automatic for `ATTR_ORDERLESS`, but `Integrate`'s output goes through user-side `Expand` which can reorder terms). Ensure final form has terms grouped: rational-part first, then `Log` terms, then `ArcTan`, then `ArcTanh` — matches Mathematica's output convention and makes diffing test expectations easy.

---

## Critical files to modify or create

| Path | Action | Notes |
|------|--------|-------|
| `src/intrat.c` | **NEW** | Bulk of the implementation (Phases 1-7) |
| `src/intrat.h` | **NEW** | Public functions invoked by `integrate.c` and direct C callers |
| `src/integrate.c` | **NEW** | Tiny dispatcher; only thing in `System\``  |
| `src/integrate.h` | **NEW** | Single public `integrate_init` |
| `src/poly.c` | **MODIFY** | Add `SubresultantPolynomialRemainders` and `PolynomialQuotientRemainder` builtins. The internal C helpers exist (the multivariate subresultant PRS at lines 943-1411 and `poly_div_rem`). We just need to expose them. |
| `src/poly.h` | **MODIFY** | Declare new builtins |
| `src/parfrac.c` | **MODIFY (small)** | Plumb `Extension` option through to the Factor call (currently hardcoded) so `Apart[expr, x, Extension -> Sqrt[2]]` works. Required for the `apartList` step in Phase 5. |
| `src/parfrac.h` | **MODIFY (small)** | Updated signature for option parsing |
| `src/sym_names.c` / `.h` | **MODIFY** | Add cached pointers for `Integrate`, `BronsteinRational`, `RootSum`, `HermiteReduce`, `IntRationalLogPart`, `LogToAtan`, `LogToReal`, `LogToArcTan`, `LogToArcTanh` |
| `src/core.c` | **MODIFY** | One-line `integrate_init()` call right after `expand_init()` |
| `tests/test_intrat.c` | **NEW** | ~60 cases (per phase + 20 end-to-end from `IntegrateRationalTests.m`) |
| `tests/CMakeLists.txt` | **MODIFY** | Add `test_intrat` |
| `Mathilda_spec.md` | **MODIFY** | Document each new builtin per CLAUDE.md |

---

## Existing functions to reuse (no reimplementation)

| Mathilda builtin / C function | File:line | Used for |
|------------------------------|-----------|----------|
| `PolynomialGCD` | `src/poly.c:3273` | `gcd[a,b,x]` in HermiteReduce, IntRationalLogPart |
| `PolynomialExtendedGCD` | `src/poly.c:3142` | `ExtendedEuclidean`, `LogToAtan` |
| `PolynomialQuotient`, `PolynomialRemainder` | `src/poly.c:3277, 3279` | `exquo`, `polymod` |
| `Resultant` | `src/poly.c:2737` | LRT resultant computation |
| `CoefficientList`, `Coefficient` | `src/poly.c:2122, 422` | `content`, `primitive`, `monic` |
| `FactorSquareFree` | `src/facpoly.c:372` | `SquareFree` wrapper |
| `Factor` (with Extension) | `src/facpoly.c:2779`, dispatches to `qafactor.c` | Real-form conversion in `LogToReal` |
| `qa_factor_with_extension` | `src/qafactor.h:154` | Algebraic-coefficient factoring inside `LogToReal` |
| `Apart` | `src/parfrac.c:18` | `apartList` |
| `Together`, `Cancel` | `src/rat.c:444, 446` | `canonic` and `zeroQ` |
| `Expand` / `expr_expand` | `src/expand.c` | All polynomial normalization |
| `D` (`builtin_d`) | `src/deriv.c:786` | HermiteReduce uses `D[d, x]` |
| `Variables` | `src/sym_names.h` | `extractConstants`, `canonic` |
| Subresultant PRS | `src/poly.c:943-1411` (internal) | Already used for multivariate GCD; we expose it. |
| `poly_div_rem` | `src/poly.c` (internal) | Backs the new `PolynomialQuotientRemainder` builtin |

---

## Verification plan

1. **Per-phase unit tests** in `tests/test_intrat.c`:
   - Phase 1: 12 cases. Polynomial integration, simple rational, Hermite reduction.
   - Phase 2: 10 cases. `IntRationalLogPart` returns the correct `{Q,S}` pair list (compare via `Sort` on `FullForm`).
   - Phase 3: 6 cases. `LogToAtan` returns the correct ArcTan sum (verify via differentiation modulo `Together`+`Cancel`).
   - Phase 4: 10 cases. `LogToReal` for quadratic/biquadratic/`n`-th-root resultants.
   - Phase 6: 8 cases. `LogToArcTan` / `LogToArcTanh` simplifications.
   - End-to-end: 20 cases drawn from `IntegrateRationalTests.m` (lines 15-1005).

2. **Universal correctness check** — the gold-standard test for any integrator:
   ```c
   /* For each test (integrand, var):
        result = evaluate(parse("Integrate[<integrand>, <var>]"));
        diff   = evaluate(parse("D[<result>, <var>] - <integrand>"));
        check  = evaluate(parse("Cancel[Together[<diff>]]"));
        ASSERT(check is 0);                                       */
   ```
   Implemented as a helper `assert_integral_correct(const char *integrand, const char *var)` in `tests/test_intrat.c`. This catches algorithmic bugs even when the printed form differs from the Mathematica reference.

3. **Reference-form check (lighter, optional)** — for a curated subset of ~15 cases where we expect output to match the Mathematica reference up to `Together[expand[diff]] == 0`, compare via `assert_eval_eq`.

4. **Memory** — run `make && cd tests/build && valgrind --leak-check=full --error-exitcode=1 ./test_intrat` per CLAUDE.md.

5. **Strict-C99 build check** — `gcc -std=c99 -Wall -Wextra -O3` on Linux (no POSIX-only deps; follow the `M_PI` guarded-fallback pattern from `src/trig.c`).

---

## Risk register

| Risk | Mitigation |
|------|-----------|
| Subresultant chain extraction differs subtly from Mathematica's `SubresultantPolynomialRemainders` | The `IntRationalLogPart` algorithm only uses two facts about the chain: degrees `Exponent[s, x]` for each `s` and the unique element of degree `i`. Hand-roll the chain by re-running pseudo-division if the existing internal function isn't directly suitable. |
| Algebraic-coefficient `Apart` is hard | Phase 2 implements `apartList` only over `Q`; algebraic-extension partial fractions are achieved by *first* `Factor[d, Extension->...]` and *then* `Apart` on the rewritten expression. This sidesteps modifying the `Apart` solver. |
| `Solve`-replacement is brittle for symbolic discriminants | Restrict `LogToReal`'s analytical solver to numeric (rational) discriminants; for symbolic, fall back to `RootSum`. The user's confirmed scope explicitly accepts this. |
| Output not bit-identical to Mathematica | Universal correctness check (`D[result] − integrand == 0` after `Cancel∘Together`) is the contract. Form matching is best-effort and does not gate phase completion. |
| Context-qualified builtin registration is novel | No existing builtin is registered with a backtick path. The context system (`src/context.c:177`) is set up to support it (absolute names with backticks). Light validation needed: register one helper in Phase 1, run the REPL to confirm `Integrate\`HermiteReduce[...]` resolves correctly, then proceed. |
| Performance on degree-32+ test cases | The Mathematica file has tests with degree up to 48 in the denominator. The PRS backbone is already O(n²) in degree; bigint coefficients via GMP keep arithmetic exact. If a specific case is slow, profile and consider porting the `extractConstants` optimization (`IntegrateRational.m:1013`) which keeps PRS coefficients small. |

---

## Open questions to revisit during execution

1. Should we expose `RootSum` as a more capable head (with `RootSum[fn, summand]` actually expanding when the polynomial factors)? Phase 2 makes it passive; we may upgrade it in Phase 4 when implementing `LogToReal`.
2. Should `Integrate[poly, x]` (pure polynomial) bypass the full pipeline for speed? The Mathematica code does. Easy optimization in Phase 1.
3. Does the Mathematica `polymod[A, B, x]` (using `PolynomialReduce`) match our `PolynomialRemainder[A, B, x]` for multivariate `B`? Yes — for univariate `B` they agree. The algorithm only uses univariate `B`. Confirmed.

---

## Approval criteria

- Plan is reviewed and approved by user.
- After approval: copy this plan to `INTEGRATE_PLAN.md` in the repo root, then begin Phase 1.
- No code changes happen during plan mode.
