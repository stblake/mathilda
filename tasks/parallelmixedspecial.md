# Port `ParallelMixedSpecial` into Mathilda as an `Integrate` method

Plan: `/Users/user/.claude/plans/let-s-port-users-user-documents-research-unified-robin.md`

Kept separate from `tasks/todo.md`, which another session is using for the M61
DSolve corpus wave.

## Target

The three existing ports of the special-function stage agree case-for-case on the
312-case corpus (`<research>/special/stress_special.py`):

| port | PASS | HONEST | FAIL | WEAK | summed case s | median | max |
|---|---|---|---|---|---|---|---|
| Python (SymPy, reference) | 305 | 7 | 0 | 0 | 374.6 | 0.49 | 45.0 |
| Mathematica 14 | 305 | 7 | 0 | 0 | 121.9 | 0.14 | 20.0 |
| Maxima | 305 | 7 | 0 | 0 | 845.6 | 1.97 | 42.9 |

Wall times are NOT comparable (Python serial, WL 6 kernels, Maxima 12 workers,
and three different per-case clocks). Report summed per-case kernel seconds with
WORKERS/TIMEOUT pinned.

## Phase 0 — core gaps  ✅ DONE

- [x] **`Block` now restores `DownValues`** (and UpValues/attributes) of its
      locals — `src/modular.c`. Previously only `own_values` was saved, so a
      `Block[{f}, f[x_] := ...; body]` left the temporary rule installed for the
      rest of the session. That is the mechanism `ExtendedBounds` is built on, so
      without it the first extended call would silently and permanently repoint
      four of Part II's bound-decision hooks.
      - Frames are also threaded on a global stack drained by
        `mth_block_depth_unwind`, called from `tc_run_guarded`
        (`src/core.c`), so a `TimeConstrained` timeout's `siglongjmp` cannot
        leave a Block binding installed. Same pattern as the async-defer count.
      - Attributes are deliberately *not cleared* (Mathematica's Block does not
        clear them either — verified).
- [x] **`Block` excluded from capture-avoiding substitution** —
      `expr_is_binding_scope` in `src/modular.c`. `Block` is *dynamic* scope, so
      a symbol arriving inside a caller's value is exactly what it means to
      rebind; alpha-renaming the local defeated the construct. Verified against
      Mathematica: `g[v_] := Block[{e = 1}, v + e]; g[e + 1]` is `3` for Block
      and `2 + e` for Module, and Mathilda now gives both.
- [x] **A18 closed for the whole iterator family** — `Do`, `Sum` and `Product`
      were all capture-prone; only `Table` was listed as a scoping construct.
      `is_iterator_scope_head` in `src/modular.c` now covers all four.
      `g[v_] := Module[{s = 0}, Do[s += v, {k, 2}]; s]; g[k]` was `3`, is now
      `2 k`; `Sum` was `3`, now `2 k`; `Product` was `2`, now `k^2`.
      Measured no performance cost (min-of-5, alternating runs, indistinguishable
      from HEAD).
- [x] **`CoefficientArrays[polys, vars]`** — new builtin in
      `src/poly/monomials.c`, reusing `build_monomials`. Matches Mathematica on
      every probe including the `x y` upper-triangular placement and
      `Modulus ->`. Returns dense Lists rather than SparseArrays (documented);
      `Normal` of a List is the identity so `Normal[CoefficientArrays[...]]`
      — the spelling `LinSolveZero` uses — works unchanged.
- [x] **Nested Association assignment auto-vivifies** — `src/part.c`.
      `a = <||>; a["k","s"] = 7` was a SILENT no-op (the Set was left
      unevaluated and `;` discarded it); now gives `<|k -> <|s -> 7|>|>` as in
      Wolfram. Sibling keys are not clobbered and value semantics hold.
- [x] **`BooleanQ`**, **`SymbolName`** (`src/core.c`),
      **`Internal`SyntacticNegativeQ`** (`src/minus.c`) — with docstrings and
      `ATTR_PROTECTED`. The last is what SymPy's `could_extract_minus_sign` is
      built on in the package; a `Plus` is deliberately never syntactically
      negative, which is why the package composes it with a `First[]` test.
- [x] **`RowReduce[..., ZeroTest -> f]` needed no work** — it already works. The
      plan listed it as a gap on the strength of a probe whose matrix
      auto-simplified before `RowReduce` saw it, so the option had nothing to
      decide. Re-probed with `Sin[x]^2 + Cos[x]^2 - 1`, which stays symbolic:
      the default reduces to the identity and the zero test correctly collapses
      the row.

## Phase 1 — the elliptic family  ✅ DONE

Needed by corpus group F, cases 100–114 (15 of 312). Mathilda had none of
`EllipticF`/`EllipticE`/`EllipticPi`/`EllipticK`.

- [x] `src/special_functions/elliptic.{c,h}` — the four heads, PARAMETER
      convention `m = k^2` (not the modulus), arity-overloaded `EllipticE` and
      `EllipticPi` as in Wolfram. Exact reductions, then numerics, then
      symbolic.
- [x] Numerics via **FLINT/Arb** `acb_elliptic_*` in `src/flint_num_bridge.c` —
      Arb already uses Mathematica's parameter convention and branch placement,
      and supplies the three hard parts a hand-rolled kernel gets wrong: the
      quasi-periodic extension off the principal strip, complex `phi` (which
      `EllipticF[ArcSin[z], m]` with `|z| > 1` produces routinely), and the
      Cauchy principal value for `EllipticPi` with `n > 1`.
      **Verified against mpmath to 25 digits** on all of: complete K/E,
      incomplete F/E, complete and incomplete Pi, a complex-`phi` case, and the
      `n = 3/2` principal value — including the sign of the imaginary part.
      (The reference values recalled from memory into the first draft of that
      test were themselves wrong; independent quadrature settled it.)
- [x] Machine-precision Carlson `R_F`/`R_D` kernels in `double`, for the
      packed/ND paths and as the no-FLINT real fallback.
- [x] Derivative rules in `src/calculus/deriv.c` — `phi`-derivatives of F, E, Pi
      (these are the integrands, which is what makes a numeric verify of an
      antiderivative close) and `m`-derivatives of K, E(m), E(phi,m).
      **Every shipped formula verified against a central difference to ~1e-16.**
      `D[EllipticF[phi,m], m]` and the `n`/`m`-derivatives of `EllipticPi` are
      left as inert `Derivative[...]` forms rather than guessed — four-term
      expressions whose signs are easy to get wrong, and an inert derivative is
      honest where a wrong formula corrupts every caller silently (BesselJ does
      the same for its order).
- [x] `SYM_EllipticE/F/K/Pi` in `src/sym_names.{c,h}`.
- [x] **End-to-end: every Mathematica reference answer in group F
      differentiates back to its integrand** at the corpus sample points —
      `#100, #101, #102, #103, #105, #113, #114` all at ~1e-29 or exactly 0,
      including the two `EllipticPi` cases with `n = 3/2` and `n = 2`.
- [x] ND/packed kernels — `EllipticK`/`EllipticE` unary and `EllipticF`/`EllipticE`
      binary in `src/ndkernels.c`, declining outside the real principal domain so
      the buffer is abandoned and the exact path answers. All three
      representations (List, packed List, visible NDArray) agree.
      **`EllipticPi` deliberately has no machine kernel** (its principal value
      needs `R_J` with the `p < 0` transformation; a wrong PV is a wrong answer).
- [x] **`Compile[]` came for free** — registering the ND kernels made
      `CompileDiagnostics` report `Compiled -> True` at scalar and rank-1 shapes
      with no lowering written, and the compiled values match the interpreter
      (out-of-domain correctly bails to the interpreter). `LogIntegral` behaves
      the same way, so for a new numeric head: do the ND kernel first and
      re-check compile coverage before assuming a lowering is needed.
- [x] `docs/spec/builtins/special-functions.md`, `tests/test_elliptic.c`
      (12 groups, all passing), changelog, `src/version.h` → **v0.238**.

## Phase 4 — the C method and its three surfaces  ✅ DONE

`src/calculus/integrate.c`, mirroring the PMT block:

- [x] `pms_lazy_load` / `builtin_integrate_pms` / `try_parallelmixedspecial`.
      The worker's contract differs from PMT's: `IntegrateSurfaceSpecial` answers
      with the PAIR `{answer, verified}`, so the builtin unwraps a **verified**
      pair to the bare antiderivative and otherwise hands back exactly what the
      package said. Nothing is invented, so the qualified surface still shows the
      package's own words and `Head[r] === List` stays the caller's test for
      "no answer".
- [x] `METHOD_PARALLEL_MIXED_SPECIAL` + `method_from_string` + dispatch; the
      definite path passes an unknown method through by name, so
      `Integrate[f, {x,a,b}, Method -> ...]` works for free.
- [x] Registration, docstring, `ATTR_PROTECTED`, and the row in `Integrate`'s own
      docstring method table (which the refpages scrape).
- [x] **Both lazy loaders now go through `mathilda_load_module`**, so the 217 KB
      Part II file is no longer parsed and evaluated twice (~1 s each) in a
      session that reaches both methods.
- [x] `Integrate::method`'s valid-name list now includes **both** parallel methods
      — it had been missing `ParallelMixedTower` as well.
- [x] One clean diagnostic when the package is absent: the loader's own
      `LoadModule::nofile` is muted so only the method-level message appears, and
      `Quiet` suppresses it.
- [x] **Cascade position: LAST, after Goursat** — a deliberate change from the
      plan, which said "before Goursat". This stage is the only one that may answer
      with a non-elementary function, and an elementary antiderivative is always
      the better answer, so every stage that can produce one gets first crack.
      It is also **not** gated on `has_pseudoelliptic_radical`, unlike PMT: that
      gate keeps the elementary-only stage off a genus>0 `F/R^p` curve where it
      grinds and closes nothing, and an elliptic pencil is exactly what this stage
      closes (`1/Sqrt[x^3-x]`, by that test pseudo-elliptic, is corpus case #100).
      Complete answers only, so plain `Integrate` can never return a partial
      `answer + Inactive[Integrate][...]`.
- [x] Verified with the package still absent: the qualified symbol and the Method
      surface both decline cleanly and plain `Integrate` is unchanged on Ei, erf,
      Si, elementary and still-open integrands. 32 integrate/calculus/risch/cherry
      suites pass.

## Phase 5 — the stress harness  ✅ WRITTEN (not yet exercised: needs the package)

- [x] `<research>/special/stress_mathilda.py`, a clone of `stress_wl.py` — the
      kernel computes the whole `rec` and `stress_judge.judge` scores it, so the
      verdict rule is not re-implemented and agreement across ports stays
      meaningful. One process per case (Print to a pipe is fully buffered, so a
      process the OS kills loses its record); `$ParallelMixedTimeBudget` raised to
      the cap, since its 45 s default would otherwise silently measure a different
      budget from the other three ports; WORKERS and TIMEOUT recorded in the
      report header, which none of the existing three runs does consistently.
- [x] `stress_compare.py` gained the `mathilda` port for a four-way table, and
      degrades cleanly while `stress_mathilda.json` does not exist.
- [x] The corpus translation was checked to produce valid Mathilda input,
      including the explicit `Tower[...]` of groups E and F.

## Unexpected: one integral changed answer, for the better

The A18 fix made `DerivativeDivides`' substitution search find the fold
`u = Log[x]` it had been losing, after which the base-field Cherry erf engine
closes the reduced integral:

    Integrate[E^(Log[x]^2), x]  ->  I Sqrt[Pi] Erf[-I (1 + 2 Log[x])/2] / (2 E^(1/4))

Derivative matches to ~1e-30; Mathematica returns the same closure.
`tests/test_risch_rde_tower.c` carried `FreeQ[Integrate[E^(Log[x]^2), x],
Integrate] === False` — written to stop the RDE *tower* over-solving, but
forbidding any answer from the whole cascade. **Restated rather than relaxed:** the
RDE tower path must still decline, and an answer from the wider cascade must be
non-elementary *and* differentiate back. An elementary claim here would still be
the defect the guard was written against.

## Pre-existing test failures, verified NOT caused by this work

Settled by building pristine `HEAD` (75cc3960) in a `git worktree` — additive
and safe while another session shares the tree, unlike `git stash`.

| suite | state |
|---|---|
| `crc_corpus_tests` | **red on HEAD** (4 diff-nonzero vs baseline 3). Its count is load-sensitive: a case over the 5 s child alarm counts as TIMEOUT, not DIFF NONZERO, so the total is not comparable across runs. My tree shows 6/2 where HEAD shows 4/5 — the two extra are `1/Sqrt[5 + 3 Tan[2 x]^2]` and `1/Sqrt[5 + 2 Tan[3 x]^2]`, two of the three the test's own comment names as acceptable, which HEAD simply never reached. My tree closes 1313 vs HEAD's 1312. **Compare the `-> DIFF NONZERO (integrand: ...)` lines, never the totals.** |
| `dsolve_m34_stress_tests` | **red on HEAD**, identically, on an idle machine. |
| `dsolve_corpus_tests` | not a failure — exit 124 was my own 400 s harness timeout on a 1204-case corpus. |
| `dsolve_stress_tests` | **red on HEAD**, identically. `t_stress_undetcoeff` asserts `Head[DSolve`UndeterminedCoefficients[y'' - 2y' + y == Cos[2x], y, x]] === List` and the head stays `DSolve`UndeterminedCoefficients` on both trees — TRIG forcing fails while `Exp[3x]` forcing solves, and full `DSolve` still answers via another path. Bisected the hard way first (modular.c, integrate.c, deriv.c each reverted in turn, all still failing) before building pristine HEAD settled it. Note it PASSED in this session's earlier full-suite run, so it is also load- or order-sensitive — another reason to compare against a HEAD build rather than against an earlier run of the same suite. |
| `dsolve_tests` | exit 142 = SIGALRM: the known `alarm(120)` cutoff, which stops the suite partway (see the standing note that units after ~48/266 never run). |

Everything else in the 512-binary suite passes.

## Phases 2 and 3 — DONE, and the first Mathilda corpus numbers

### Phase 2 — the Part II merge

All 11 items of the 276-line additive delta hand-merged into
`src/internal/mixed/ParallelMixed.m`, preserving the 794 lines of Mathilda work.
**Gate met: Part II is behaviourally unchanged** — `parallelmixedtower_tests`
11/11, and a 30-case differential (`Log[x]`, radicals, `Sqrt[Tan[x]]`, the
nonelem certificate, the declines, …) is **byte-identical** pre- and post-merge.
`Exp[2x]/(1 + Exp[x])` now closes as `E^x - Log[1 + E^x]` through the new
structure-theorem path, and `Exp[x] Exp[x^2]` merges to `E^(x + x^2)` on input.

Two splice bugs worth recording, both caught by the differential rather than by
eye: the head splice **overlapped** the reference range (duplicating
`Y = Unique`, `nc`, `f = TPad` and swallowing `unitsBase = {};`), and the iPIM
tail splice stopped one line short, dropping `units = Join[sunitsAll, unitsBase]`.
The second is the instructive one — it lost the *units* silently, so the ansatz
went from 9 unknowns to 8 and Part II answered `{"not elementary", …}` for
`Sqrt[x + Sqrt[x]]`, **a false certificate**. A line-range splice needs an
assertion per seam, not a visual check.

### Phase 3 — the package port

`src/internal/ParallelMixedSpecial.m`, 1868 lines, near-verbatim from the
1825-line `.wl`. Loaded **into** `ParallelMixed`Private`` (the `logrewrite.m`
pattern) so Part II's 41 internals resolve by short name; the original's
`AppendTo[$ContextPath, "ParallelMixed`Private`"]` cannot work here. Only forced
change in the body: five `FreeQ[…, Complex]` → `_Complex` (A19). **The A18
iterator renames that Part II needs were NOT needed** — fixing A18 at the root in
Phase 0 means this file keeps the original's iterator names throughout.

### The numbers (312 cases, TIMEOUT=120 s, WORKERS=6)

| port | PASS | HONEST | WEAK | FAIL | summed case s | median | max |
|---|---|---|---|---|---|---|---|
| Python (SymPy, reference) | 305 | 7 | 0 | 0 | 374.6 | 0.49 | 45.0 |
| Mathematica 14 | 305 | 7 | 0 | 0 | 121.9 | 0.14 | 20.0 |
| Maxima | 305 | 7 | 0 | 0 | 845.6 | 1.97 | 42.9 |
| **Mathilda (first run)** | **247** | **51** | **11** | **3** | 313.8 | 0.08 | 120.0 |

Seven groups are at parity: G not-in-class 14/14, H partial 10/10, Q 20/20,
R elementary 15/15, R radicals 5/5, S elementary 15/15, E radical (same 2 HONEST
as every port). **No false certificate anywhere** — the soundness bar the corpus
exists to police.

**57 of the 58-case gap share ONE root cause.** 46 of the 51 HONEST are
`{"failed", "special answer found, but the integrand is not certified
non-elementary", answer}` — the answer is found and *matches the Mathematica
reference exactly* (`exp(-x^2)` → `-Sqrt[Pi] Erfc[x]/2`, `exp(-x^3)` →
`-Gamma[1/3, x^3]/3`, `log(x) sin(x)` → `CosIntegral[x] - Log[x]/(1 + Tan[x/2]^2)
+ …`, character for character) and strict mode withholds it for want of a
certificate. All 11 WEAK are the same thing one level down: the partial mode left
exactly those sub-integrands (`E^x^2`, `E^(-x^2)`, `E^(-x^3)`, `Sqrt[Log[x]]`,
`ArcTan[x]/x`) in the remainder because it could not certify them.

The cause is in **Part II, not the port, and predates this work**:
`src/internal/mixed/ParallelMixed.m` guards the holomorphic-remainder certificate
with `q =!= None`, which the research reference does not have. The changelog for
v0.161 records it as "the T2/T10 false-certificate protection, which withholds the
new (K4) certificates for `Exp[x^2]`, `Exp[-x^2]`, `Sqrt[x] Exp[x]`" — it was
added because Mathilda emitted **false** certificates on curve-free towers.
Lifting it would hand back ~57 cases and risk the worst failure mode, so the next
step is to find why Mathilda's `VerifiedResidueFree` / `SecondKindAtInfinityQ`
were unsound there, not to delete the guard.

The remaining 5 HONEST are genuine coverage (`"no solution with the special
kernels"`). The 3 FAIL, each attributable and none to the port itself:
**#34** `asin(sqrt(x+1))/sqrt(x)` comes back `verified = False` with a residual
that is purely *imaginary* (`0.0 - 1.29 I`) — the real part is exact and the
`ArcSin` branch in the reconstructed surface is placed differently from the
integrand's. **Part II alone fails the same integrand** (`{failed, "verification
failed"}`), so this is a pre-existing branch-fidelity gap, and the stage does the
right thing by reporting the answer unverified rather than claiming it. And
**#109 / #110 timeouts at 120 s** — the two third-kind elliptic cases needing the
mod-p non-torsion certificate, which take 45 s / 43 s in Python and 42.9 s in
Maxima. That is the GF(p) Jacobian performance gap already in
`MATHILDA_DIVERGENCES.md` §C (512 s against Python's 0.8 s on a genus-4
certificate).

**Speed, like for like.** Mathilda's summed 313.8 s is not comparable to the
others', because 57 cases bail out early at the certificate step instead of doing
the full work. On the 247 cases Mathilda actually PASSes:

| port | summed | median | max |
|---|---|---|---|
| **Mathilda** | **41.2 s** | **0.070 s** | 3.8 s |
| Mathematica | 56.0 s | 0.110 s | 4.3 s |
| Python | 188.2 s | 0.390 s | 11.6 s |
| Maxima | 569.7 s | 1.930 s | 22.8 s |

Median per-case ratio, restricted to cases the other port takes over 0.2 s:
Mathilda is **0.62×** Mathematica, **0.16×** Python, **0.04×** Maxima. So on the
work it does complete, Mathilda is the fastest of the four — but that is a
statement about 247 cases, not 312.

### Two harness bugs found on the way, both silent

- **`AbsoluteTime[]` has integer-second resolution in Mathilda.** The research
  runners time each case as an `AbsoluteTime[]` delta; here that reads `0.0` for
  anything under a second, so every per-case time was 0 until the runner switched
  to `AbsoluteTiming`. Recorded as divergence F7.
- **The two package contexts.** Part II exports seven *public* symbols and keeps
  ~200 private; the special stage is loaded into the private context. Addressing
  either with the wrong context leaves the call unevaluated, which reads
  downstream as a wrong *answer* ("differs from Part II"), not an error — it made
  all 41 group-A cases FAIL. Recorded as F8.

## Remaining phases

2. Hand-merge the 276-line additive Part II delta into
   `src/internal/mixed/ParallelMixed.m` (which is 794 diff lines AHEAD of the
   reference, not behind — it carries the field-aware `Can`, the memos and ~20
   divergence workarounds, so this is a merge, not a copy).
3. Port `ParallelMixedSpecial.wl` (1825 lines) to
   `src/internal/ParallelMixedSpecial.m`, as a context-injected file loaded from
   inside `ParallelMixed`Private`` (the `logrewrite.m` pattern) — Mathilda
   resolves a short private *variable* through `$ContextPath` but not a short
   private *function*, so the package's own `AppendTo[$ContextPath, ...]` trick
   cannot carry its 41 short names.
4. The C method and its three surfaces in `src/calculus/integrate.c`.
5. `special/stress_mathilda.py` + the four-port `stress_compare.py`.
6. Tests, docs, `src/version.h` bump and tag.
