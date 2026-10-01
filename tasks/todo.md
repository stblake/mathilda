# Task — execution-speed comparison of the four ParallelMixedSpecial ports

Goal: compare the **execution speed** of the Python (SymPy), Mathematica and Maxima
implementations of the parallel mixed tower *special* algorithm on the 312-case stress
corpus of `MATHILDA_PARALLEL_MIXED_SPECIAL_PLAN.md`, with Mathilda v0.244 as a fourth
column.

Research tree: `~/Documents/Research/post_phd_research/algebraic_integration/risch_norman_radicals/special`
(referred to as `<research>/special` in the plan document).

## Why the recorded runs cannot answer the question

The four `stress_*.json` files on disk each carry a per-case `time`, and
`stress_compare.py`'s own docstring already says not to quote them as a speed result.
Measured/read confirmations:

- [x] **Different worker counts.** Python ran **sequentially** (`stress_special.main`
      has no executor), WL at **6** kernels, Maxima at **8**, Mathilda at **6** — on an
      8-physical-core i9-9880H. The parallel ports' per-case clocks are inflated by
      contention the sequential one never paid.
- [x] **Different timed regions.** 136 of 312 cases have *harness* work inside at least
      one port's clock:
      - WL and Mathilda bracket the harness's Part II cross-check
        (`ParallelIntegrateMixed`) together with the call — 61 `elem` cases.
      - Maxima prints `elapsed_real_time() - t0` *after* the whole post-processing
        block (`pms_rt_t1` is captured and never used), so its clock carries the
        remainder `ratsimp`/`radcan`/`bfloat` probes (75 partial cases) **and** a second
        full Part II integration (61 `elem` cases).
      - Python brackets `_run`, i.e. the call plus `sp.simplify(rem - expect)` and the
        Part II cross-check.
- [x] **Python ran all 312 cases in one process**, so SymPy's global caches were warm
      across cases; the other three start a fresh process per case.
- [x] Measured instance: case 1 (`log(x)`) reads **2.0 s** in the recorded Maxima run
      and **0.6 s** re-run at 2 workers.

## Plan

- [ ] 1. `speed_bench.py` in `<research>/special`: one uniform clock across the four
      ports. Timed region = **exactly the entry-point call**
      (`integrate_surface_special` / `integrate_surface_partial` and each port's
      spelling of them) — no package load, no process start-up, no harness
      verification, no Part II cross-check. One case per OS process for **every** port,
      Python included. Reuses the existing corpus (`stress_special.C`) and the existing
      expression translators (`stress_wl.to_wl/tower_wl`, `stress_maxima.mx/tower_code`)
      so nothing about the mathematics is re-implemented.
- [ ] 2. `WORKERS=1` (sequential) for the measurement run — no contention, and it also
      sidesteps the Mathematica licence stall seen at 2 concurrent kernels. `CAP=120` s,
      the corpus cap, pinned in-kernel where the port has an inner budget (Mathilda:
      both `$ParallelMixedTimeBudget` and `$SpecialTimeBudget`, *after* the lazy-load
      line).
- [ ] 3. Measure each port's **fixed per-case floor** separately (process start-up +
      package load), 5 reps: it is not part of the algorithm's speed but it is what a
      user pays, so it is reported as its own column rather than smuggled into the clock.
- [ ] 4. Run all four ports, 312 cases each. Cross-check every case's `kind`/`status`
      against the recorded verdict run — a case that is fast because it now errors is
      not a speed result.
- [ ] 5. `speed_report.py`: verdict-agreement check, then per-case ratios, medians,
      totals on the common subset, per-group breakdown, slowest cases per port, and the
      cases where the ranking inverts.
- [ ] 6. Write the findings up. Research-tree measurement of research-tree ports, so no
      Mathilda source change, no version bump, no changelog entry unless a Mathilda
      defect is found.

## Review

All six items done. Two runs were needed, because the first uniform-clock run was
itself unfair.

**Deliverables** (all in `<research>/special`, nothing in the Mathilda tree changed
except this file and the plan-document correction below):

- `speed_bench.py` — the uniform-clock harness, four ports, `--floor`, `--port=all`.
- `speed_report.py` — console report: integrity, three bases, ratio matrix, groups.
- `speed_tex.py` → `PARALLEL_MIXED_SPECIAL_PERF_COMPARISON.{tex,pdf}` (10 pages).
  Every figure is computed from the JSONs by the generator; none is transcribed.
- `speed_{python,wl,maxima,mathilda}.json`, `speed_floor.json`, `speed_cases.tsv`,
  and the discarded first run as `speed_*_nowarm.*`.

**Result** — 282 cases every port answers, algorithm time only, sequential:

| port | common sum | median | p90 | max | floor |
|---|---|---|---|---|---|
| Mathilda 0.244 | 50.8 s | 0.061 | 0.293 | cap | 0.05 s |
| Maxima 5.49 | 63.7 s | 0.043 | 0.429 | 12.46 | 0.69 s |
| Mathematica 14.0 | 69.1 s | 0.083 | 0.318 | 11.95 | 2.98 s |
| Python/SymPy | 154.6 s | 0.262 | 1.392 | 26.26 | 0.53 s |

Integrity: **0 of 312 outcomes changed** for any port against the verdict runs.

**The measurement finding that mattered more than the table.** A first pass put Maxima
at 227 s and 6.45x slower than Mathematica. Wrong: Maxima autoloads/compiles library
code on the first *real* call (0.483 s, then 0.013 s), not at `load()`, so excluding the
package load did not exclude it, and with a fresh process per case it landed on all 312
— 147 s of the 227 s. Python, Mathematica and Mathilda have no such cost, and
Mathilda's runner *already* warmed itself via its lazy-load line, so the port that
looked best was the only one getting the treatment. Warming every process inverted the
ranking: Maxima 227.3 → 79.3 s, now the lowest median of the four. Verified no residual
per-branch autoload (case vs same-group sibling in one warm process: erf 0.062/0.058,
elliptic 0.084/0.084, Ei 0.012/0.012, partial 0.030/0.023). Lesson saved as
`memory/feedback_fresh_process_benchmark_charges_lazy_init.md`.

**Also corrected**: `MATHILDA_PARALLEL_MIXED_SPECIAL_PLAN.md` §4 quoted
`Python 45.0 s / Maxima 42.9 s / ~3 s` for #109/#110 from the contention-inflated
verdict runs, and omitted Mathematica. Clean-clocked it is Python 26.3/22.3, Maxima
10.2/1.04, **Mathematica 0.83/1.65**, Mathilda >120 both. Item 3 is still worth only 2
cases, but the 0.83 s shows the 120 s wall is not intrinsic to the certificate.

No Mathilda source change, so no version bump and no changelog entry.

---

## Notebook: syntax highlighting + bottom-up structural selection (2026-10-01)

Frontend-only (Tauri/Svelte/CodeMirror). No kernel/Rust changes (confirmed with
user: bracket/token-level selection, capitalization heuristic for symbols).

- [x] `frontend/src/lib/mathildaLex.ts` — pure lexer (scanToken, tokenize) mirroring parse.c lexical rules
- [x] `frontend/src/lib/mathildaLang.ts` — CodeMirror StreamLanguage + HighlightStyle (CSS-var colours)
- [x] `frontend/src/lib/structSelect.ts` — buildSpans / ladderAt / chooseExpand (pure)
- [x] `frontend/src/lib/CellShell.svelte` — wire language+highlight, multi-click handler, Alt-Up/Down, cache invalidation
- [x] `frontend/src/App.svelte` — `--cm-*` colour vars (dark `:root`, light `html.light`)
- [x] `frontend/package.json` — add `@lezer/highlight`, `check:selection` script
- [x] `frontend/scripts/check-selection.mjs` — node test of tokenizer + span ladder
- [x] Verify: `npm run check` (tsc/svelte-check), `npm run build`, `npm run check:selection`

### Review

Done, frontend-only (no kernel/Rust change → no version bump, tag, spec, or
changelog). Shared pure lexer `mathildaLex.ts` (mirrors parse.c lexical rules)
feeds both: a CodeMirror `StreamLanguage` + themed `HighlightStyle`
(`mathildaLang.ts`) for highlighting, and bracket/token nesting spans
(`structSelect.ts`) for bottom-up selection. In `CellShell.svelte`: highlight
bundle swapped in for the no-op `defaultHighlightStyle`; a left multi-click
(`event.detail`) selects the token then one enclosing bracket level per further
click (Mathematica's model); Alt-Up/Down do the same from the keyboard; the span
cache invalidates on edit. `--cm-*` colour vars added to both themes in
`App.svelte`.

Verified: `check:selection` (new, 21 assertions) PASS incl. the exact chosen
behaviour (`a + b*c` → `b` → whole, no `b*c` step); `npm run check`
(svelte-check + tsc) 0 errors; `npm run build` OK; existing check:prose/search/
snippets/notebook all still PASS. **Not** verified here: the live GUI (colours in
light/dark, click feel) — needs `cd frontend && npm run tauri dev` with the
sidecar built, which is an interactive desktop app.

Known MVP limits (documented in the plan as future, kernel-backed upgrades):
selection is delimiter/precedence-free (no `b*c` sub-select); symbol colouring is
the capitalisation heuristic (not the real defined/undefined split); chained
calls `f[x][y]` don't merge heads; genuine `[[ ]]` Part may show one extra level.

### Follow-up (same day): user ran the bundled app

Two issues surfaced when the user launched via MathildaNotebook.command:

1. **Clickable selection "completely broken."** Root cause: the first handler keyed
   expansion on `MouseEvent.detail` (the browser's RAPID multi-click counter), so
   slow deliberate clicks never expanded. Rewrote to a pace-independent gesture:
   a 2nd click on the same atom, or any click inside the current selection, expands
   one level (fast double-click still starts it); wrapped in `Prec.highest(...)` to
   preempt CM's word/line select and in try/catch. `CellShell.svelte` only.

2. **`*::nofile` errors loading `.m` files.** NOT the frontend — the bundled kernel
   runs from an arbitrary cwd and the resolver's search paths (`$MATHILDA_HOME` →
   exe-rel → prefix → cwd ladder) found nothing, so EVERY `.m` load failed (init.m
   silently; CRCTable/ParallelMixed loudly). Pre-existing in every bundle; the
   Integrate example just surfaced it. Fix: `tauri.conf.json` bundles `src/internal`
   as resource `internal`; `kernel.rs` `spawn_inner` sets
   `.env("MATHILDA_HOME", resource_dir()/internal)` when present (dev falls back to
   the cwd ladder). Self-contained → won't recur on rebuilds.

Verified: frontend `npm run check` 0 errors; `check:selection` PASS; `cargo check`
clean; rebuilt bundle contains `Resources/internal/` (13 .m files, subdirs intact);
and a kernel pipe-mode test from cwd=/ shows the four errors WITHOUT MATHILDA_HOME
and a clean run WITH it (integral result identical). Lessons added to memory
[[project_notebook_frontend_highlight_and_selection]].

### Follow-up 2: precedence-aware selection (user's explicit spec)

User's test: clicking x in `Integrate[(2 x + 1) Sin[x^2+x+1] - x Log[x] + 1, x]`
must ladder x -> 2 x -> 2 x + 1 -> (2 x + 1) Sin[...] -> whole. That is
PRECEDENCE-aware, which the frontend-only bracket/token model (the earlier chosen
scope) cannot produce. Reversed that decision and built the kernel-authoritative
approach:

- C: `mth_parse_spans` + opt-in `MthSpanSink*` on `ParserState` (NULL on the
  normal path); one `record_span` at the Pratt-loop top captures each precedence
  stage; `flat_continuation` suppresses partial Plus/Times sums; `record_span`
  trims trailing whitespace. Exposed in `src/parse.h`. (src/parse.c)
- Pipe: `{"id":N,"expr":"...","spans":true}` -> `{"type":"spans","payload":[[s,e],..]}`
  via `pipe_process_spans`/`pipe_emit_spans`; one-line route in pipe_mode_loop.
  (src/repl.c)
- Rust: `MathildaKernel::fetch_spans` + `syntax_spans` command + registration.
- Frontend: `fetchSpans` (ipc.ts, UTF-8->UTF-16, raw text); gesture rewritten to
  store {anchor, level} and recompute the ladder each click so the first-click
  local fallback upgrades to kernel spans; first click selects the token, further
  clicks climb; `CellShell` caches kernel spans per doc + prefetch on focus/edit.

Also fixed the earlier "second attempt selects whole expression" bug (anchor the
gesture on the token, not "is the click inside the selection").

Verified: kernel pipe test gives the exact ladder for `2 x + 1`
([[0,1],[2,3],[0,3],[6,7],[0,7]]) and for the integrand; `npm run check` 0 errors;
`check:selection` PASS (incl. the first-click-selects, different-token-restart,
and local->kernel upgrade cases); `cargo check` clean. Pending: version.h bump +
tag + changelog when committed; optional C unit test for mth_parse_spans in the
CMake suite.

### Follow-up 3: output "Convert To" right-click menu

Right-click an output expression -> Convert To -> StandardForm / InputForm /
FullForm / TeXForm / MathML, switching the display form in place.
- Rust: generic `eval_once(expr) -> {payload,latex,error}` command (non-cell eval,
  no $Line pollution) — kernel.rs/commands.rs/lib.rs; `evalOnce` in ipc.ts.
- Output.svelte: contextmenu on `.out-expr`, per-item form + fetched-string cache,
  `exprHtml()` form-aware render, `.convert-menu` popup. InputForm uses `item.text`;
  StandardForm/TeXForm/MathML derive from `item.latex`; FullForm (and TeX/MathML
  without latex) fetched via evalOnce. No C change.
Verified: kernel non-cell eval gives FullForm `Plus[1,Power[x,2]]`, TeXForm
`1+x^{2}`, InputForm `1 + x^2`; `npm run check` 0 errors; `cargo check` clean.

---

## No slow cases in `Integrate`ParallelMixedSpecial`, and no segfaults (2026-10-01)

Goal: eliminate the slow cases in the special-function stage — the centerpiece of
Mathilda's integration suite. Two constraints set by the user: the algorithm must be
**correct by construction with no verification step** inside it (verification belongs in
tests), and there must be **no segfaults**. Plan: `~/.claude/plans/curious-munching-key.md`.
Corpus baseline at v0.244: **303 PASS / 9 HONEST / 0 WEAK / 0 FAIL**.

### What the profiling found — four unrelated mechanisms, not one

- **(A) Cost is multiplicative in the number of independent generator families appearing
  additively.** `Log[x]/x` 0.021 s, `Sin[x]/x` 0.216 s, `Exp[-x^2]` 0.063 s; the two-term
  sum 5.75 s (24x the parts), the three-term sum **20.12 s (67x)**. One ansatz spans the
  *product* of the families.
- **(B) A live segfault**, repro `D[Inactive[Integrate][f[u], Sqrt[x]] + 1, x]`, firing
  seven times in one 21-second window from a test suite.
- **(C/D) #109 and #110 fail for two different reasons, neither the recorded one.** §4 of
  the plan document blames GF(p) Jacobian speed; that arithmetic already landed and is
  60–100x faster than the figures §4 reasons from. #109 is a missing quartic Legendre
  branch in `ThirdInt` plus the crash; #110 is `RealisePoints`/`TorsionRealise` on
  un-normalised nested-radical coordinates (`ToNumberField` gives a degree-**4** field in
  9.8 ms, so "degree 8" is a spelling artefact).
- **(E) The trig/Si–Ci family is 4–5.7x Maxima across ~16 cases.** `Sin[x]/x` costs
  0.194 s where the same mathematics as exponentials costs 0.082 s: a `Tan[x/2]`
  Weierstrass tower is built and `e^{ix}` then recovered back out of it.

### Phase 1 — the segfault, fixed at root  [done]

`deriv_of` (`src/calculus/deriv.c:545`) was a pass-through forwarding `compute_deriv`'s
"cannot differentiate" NULL into callers that store it straight into an argument slot;
`evaluate_step` then dereferences it (`src/eval.c:1471`). **82 call sites, none checking.**
`deriv_of` now returns the inert `D[g, x]` — the behaviour its own docstring already
promised — so the class is impossible rather than merely absent. Answers improve too:
`D[INT + g[x], x]` gives `Derivative[1][g][x] + D[INT, x]`. The wrap stays `D` and not
`Dt` on purpose (`deriv.c:2850`: re-emitting `Dt` re-enters `builtin_dt` unboundedly).

### Phase 2 — additive decomposition  [done]

**The soundness correction that mattered most.** My first equivalence — "terms that cancel
must share a generator" — is **false**, and it survived every probe because I probed
exponentials, which do share a generator. The counterexample is a pair of *primitives*:

    Log[1+x]/x + Log[x]/(1+x)  ->  Log[x] Log[1+x]   elementary
    Log[1+x]/x   alone         ->  "not in class"
    Log[x]/(1+x) alone         ->  "not in class"

Disjoint families, elementary sum, neither half even in the class, because for primitives
with base derivatives `D(c t_A t_B) = c (D t_A) t_B + c t_A (D t_B)` is a two-term sum with
*disjoint* support. A naive split turns this PASS into a decline.

The guard: let `N` be the blocks that did not close and `P` the **coupling-capable** ones
(a `Log` or an `InvRational` head at an argument rational in `x` alone). Compose when
`|N and P| <= 1`, else decline to the joint path. The degree argument is written out at the
call site with the one residual assumption (the *class* boundary is not proved).

| case | before | after | |
|---|---|---|---|
| `Log[x]/x + Sin[x]/x + Exp[-x^2]` | 20.12 s | **0.287 s** | 70x |
| #252 | 6.38 s | **0.28 s** | 23x |
| #251 | 2.62 s | **0.12 s** | 22x |
| #129 | 3.57 s | **0.38 s** | 9.4x |
| #199 | 3.30 s | **0.41 s** | 8.0x |
| #250 | 0.61 s | **0.12 s** | 5.1x |
| #233 | 1.16 s | **0.28 s** | 4.1x |

Corpus kernel time excluding #109/#110: **64.5 s -> 43.9 s (-32%)**. Cost on cases that do
not split is 4–10%. `#199` composing despite two coupling-capable blocks is corollary C2
doing real work: its `ArcTan` block *closes*, so it leaves `N`.

### Verification

- **Corpus: 303 PASS / 9 HONEST / 0 WEAK / 0 FAIL, `verdict changes: NONE`**, HONEST the
  same nine ids `[59, 97, 98, 109, 110, 274, 275, 276, 277]`.
- **PMT 114: 94 solved / 17 nonelem / 3 declined / 0 wrong**; **Charlwood 50: 49 verified**
  — both baseline. (`charlwood_record.py` rewrote section E of `MATHILDA_DIVERGENCES.md`
  with stale template prose contradicting its own measured run; restored from backup. The
  table did not change.)
- **Unit suite: 506 pass / 9 fail / 0 crash** of 515, all nine accounted for:
  `crc_corpus_tests` pre-existing (compared *case identities*, not totals — the documented
  trio plus the `1/Sqrt[a + b Tan[c x]^2]` family its own baseline names acceptable; two
  verified to `Simplify` back to 0); `dsolve_m34_stress` red on pristine HEAD;
  `dsolve_stress`/`dsolve_tests`/`dsolve_corpus` load-sensitive, `alarm()`, or my 300 s cap;
  `parallelmixedspecial_tests` a stale mid-sweep binary; plus the two below.
- New tests: `test_undifferentiable_subexpression` (one case per affected head) and four
  groups in `tests/test_parallelmixedspecial.c` — both counterexamples, split-vs-joint
  agreement *including remainder equality*, the complete-answer property, and faithfulness.

### Two cleanups the sweep forced

- **Orphaned test binaries removed.** `intrischnorman_tests` and `int_rnb_tests` are
  binaries from the C RischNorman/RischNormanBlake implementations deleted in `1a20f7f9` —
  not CMake targets (`make`: "No rule to make target"), sources gone, dated 18 Sep. A
  `ls *_tests` sweep of `tests/build` picks such orphans up and their failures are noise.
- **`moebiusmu_tests` / `primenu_tests` strengthened.** They asserted oracles depending on
  ECM splitting a large composite cofactor; when it fails, `FactorInteger` returns it
  "unfactored with exponent 1" and both functions compute a **confidently wrong answer**.
  Three consecutive runs of the old `PrimeNu[...]` assertion on one unmodified build gave
  **7, 8, 8**. Replaced with products of explicit four-digit primes — above `2^63` so the
  bigint path still runs, every factor inside trial division's reach, expected value known
  by *construction* — plus the `mu = (-1)^nu` / `nu <= omega` identities. 4/4 deterministic.
  The hazard itself is left as a finding: these heads return a wrong number rather than
  staying unevaluated when the factorisation is incomplete.

### Not yet done

- **#231 `Sqrt[Log[x]] + ArcTan[x]/x` is not improved** (3.9 s). Both blocks are
  coupling-capable and neither closes, so the guard declines — conservatively, since no
  coupling actually occurs. Recovering it needs the sharper test of whether block A's
  integrand really has a component along `g_A * D g_B`; that is exactly where a mistake is
  a wrong answer, so the sound-but-conservative rule went in first. #227 is a 7%
  regression, as predicted at design time.
- Phases 3–8 (remove verification / construct the branch and embedding; the quartic
  third-kind hole; the trig front end; #110; the small items; docs).
- **No version bump or commit.** This tree is shared with another active session
  (`frontend/*`, `src/parse.c`, `src/parse.h`, `src/repl.c` are theirs), so `src/version.h`
  is left alone to avoid a collision — and the binary measured here also compiled their
  in-progress parser/repl changes. The corpus matching baseline exactly argues against any
  interaction, but it is not a clean-room measurement.

### Follow-up 4: converted output must be selectable (structural like input)

User: converting output to InputForm/FullForm must allow the same depth-first
clickable selection as input; TeXForm just needs ordinary text selection.
- Extracted the input's click-selection into a shared CM extension
  `structuralSelection.ts` (ViewPlugin + Prec.highest mousedown + Alt-Up/Down);
  CellShell refactored to use it (removed its inline duplicate).
- New read-only `CodeView.svelte` (EditorState.readOnly) renders converted output:
  InputForm/FullForm with highlight + structural selection, TeXForm/MathML as plain
  selectable text. Replaces the non-selectable static {@html <code>}.
Verified: `npm run check` 0 errors (150 files), `check:selection` PASS. Input
selection logic unchanged (same shared gesture) — reverify in-app after relaunch.

### Follow-up 5: fix InputForm/FullForm output selection (focus-steal)

Converted-output structural selection failed because CellShell's `.cell-content`
on:click focused the INPUT editor on any click (incl. output), wiping the output
selection. Fixed: skip input-focus when the click is in `.output-pane`; focus the
clicked view in the shared onMousedown; add drawSelection() to CodeView.
`npm run check` 0 errors. (Save-on-close prompt: still pending — see below.)

### Follow-up 6: save-on-close prompt

Added unsaved-changes tracking + a Save/Don't Save/Cancel prompt on window close.
- notebook.ts: global `dirty` store + markDirty() in the 8 content-mutating store
  methods (not output/status/exec, since serialize() saves only {type,source}).
- App.svelte: markClean() on save success + .lb open + startup; saveFile/saveFileAs/
  doSave now return boolean; `getCurrentWindow().onCloseRequested` → if dirty,
  preventDefault + in-app 3-button modal → Save (saveFile then win.destroy) /
  Don't Save (destroy) / Cancel (stay). Enter=Save, Esc=Cancel.
- capabilities: + core:window:allow-destroy.
Verified: `npm run check` 0 errors; selection test PASS. Caveat: macOS Cmd+Q may
bypass onCloseRequested (red-X close is covered).

### Follow-up 7: 2D render all trig/hyperbolic functions (v0.246)

Integrate[1/(-1+x^2),x] → -ArcTanh[x] didn't typeset (fell back to text) because
print_latex.c's FUNC_MAP lacked ArcTanh et al. Added ArcCot/ArcSec/ArcCsc, Coth/
Sech/Csch, and all six inverse hyperbolics (\operatorname{...} where no KaTeX
builtin; \coth is builtin). Verified via kernel: all now emit proper LaTeX.
version.h 0.245→0.246; changelog note added. (Pending commit → tag v0.246.)

---

## Review — "why are we slower than Maxima/Mathematica", and 109/110 (v0.248)

### What the question turned out to be

Two different questions, with two different answers, and the published table
invites confusing them:

- **all-312 total** (Mathilda 277.6 s vs Maxima 79.3 / MMA 74.6) is **86.5% two
  cases**: #109 and #110 at the 120 s cap = 240.0 s. Excluding them Mathilda is
  37.5 s against Maxima 68.0 and MMA 72.1. Not a breadth problem.
- **per-case ratio to the best other port** is the real answer: median **1.25x**,
  geometric mean 1.21x, **194 of 282 cases slower** than the fastest other port.
  Mathilda wins the common total only because it has no catastrophic tail (max
  4.40 s vs Maxima 12.46 / MMA 11.95). Worst clusters: D trigonometric 4.92x
  median (worst #84 sin(x)/x at 7.10x), B Ei/li 2.15x, P families 1.95x.

### Where the time is (measured, not guessed)

Added `$SpecialProfile` to ParallelMixedSpecial.m, because a `sample` profile
cannot attribute a `.m` pipeline — it names `evaluate_step`/matcher/malloc. The
stage attribution sums to 99.7% of the entry-point clock over 310 cases:

- **48% is Part II re-entries** — 748 calls, ~2.4 complete elementary
  integrations per case. E1 27% (the architecture), **E4 14%**, E2 6%.
- 33% the special-stage algebra proper, 8% certificates, 4% BuildTower.

**E4 was pure waste**: 344 Part II runs, **0 replacements**. Not luck —
a single gamma kernel's column is non-elementary by construction (s is never a
positive integer, so Gamma[s,.] never degenerates; Liouville/Rosenlicht). Added
`SingleKernelNonelementaryQ`; E4 now probes only the whole special part and only
with >= 2 terms. Elliptic kinds deliberately not claimed. **-10.6% corpus-wide**,
verdict run byte-identical (303/9/0/0, zero verdict changes, identical `why`).

### #109 / #110 — not the special stage at all

Both record **zero completed stages**: the whole 120 s is in the FIRST one,
Part II. The branch is `TorsionRealise` — when the mod-p non-torsion certificate
fires Part II declines in 0.1-0.5 s; when it does not, the symbolic group law
runs. **6 of 18 sibling quadratics hang the same way**, and the number field is
not the discriminator (x^2-x-1 fast / x^2-3x+1 hangs, both Q(sqrt5)).

- #109: the group law is never normalised over `Root[]` coordinates — add[P,P]
  8 ms, add[2P,P] 41 ms and a page of `Root[...]^12`, squaring per doubling,
  with a `Root::conv` flood to 4928 bits. `RootReduce` of its defining relation
  is 1.3 ms.
- #110: group law instant (TorsionOrder = 4 in 0.5 ms), logands built in 1.3 s;
  the sink is after that, in poly-GCD/degree churn on growing Q(sqrt2)
  coefficients.

Both point at one fix: normalise `PointsOver`'s coordinates into a single number
field before the group law (`ToNumberField` gives degree 4 in 9.8 ms).

### Two core bugs found and fixed on the way

- **String literals truncated at 255 chars** (`parse_string`, `parse_symbol`).
  Silent. Found because 6 of 312 answers would not re-parse. Buffers now grow;
  `test_long_string_literal_not_truncated` checks the tail, not just the length,
  and that two 301-char symbol names stay distinct.
- **TeXForm dropped parentheses on juxtaposed factors in a fraction**:
  `a/(b (c+d))` -> `\frac{a}{b c+d}`, a wrong answer in typeset form, and the
  shape of every third-kind elliptic antiderivative. Multi-factor slots now
  render at `*` precedence. `test_texform_fraction_parenthesisation`.

### Report

`PARALLEL_MIXED_SPECIAL_PERF_COMPARISON.tex` gained **Appendix A: every case in
the suite** — all 312 in corpus order with the four clocks and Mathilda's own
answer, typeset by Mathilda's `TeXForm` via the new `answers_tex.py`. 34 pages,
0 errors, 0 overfull hboxes.

### Open / next

- The torsion-realisation normalisation (closes #109/#110 and 4 more).
- The trig front end: on the Tan[x/2] tower the residue is `2t/(1+t^2)` where the
  exponential tower has a bare generator, so EiCandidates costs 33.6 ms vs 6.3
  and DecompResidue 8.7 vs 0.9. This is the 4.92x group-D median.
- Part II re-entries at 48% — E1 is the architecture, but worth asking whether
  the E1 result can be reused rather than recomputed.
- `Hash[]` is not a builtin, so `Part2`'s memo key `Hash[{...}]` stays
  unevaluated and compares structurally. It WORKS (exact, no collisions) but is
  O(size). Implementing `Hash` would silently turn that exact memo into a
  collision-prone one — fix the memo first if `Hash` ever lands.

---

## Review — the torsion realisation for #109/#110 (v0.249)

### Outcome

**#109: 120 s cap -> 0.88 s, and HONEST -> PASS** (it now returns a special-function
answer rather than declining). Corpus **304 PASS / 8 HONEST / 0 WEAK / 0 FAIL**,
one verdict change, that one. Verdict-run total 282.9 -> 166.6 s.
**#110 is not fixed** and is blocked on a separate, documented core defect.

Four more of the swept `y^2 = x^3-x` denominators that hit the budget now resolve
in 0.40-0.66 s: `x^2-6`, `x^2-2x+2`, `x^2-5x+1`, `x^2-3x+1` (#109).

### What it actually was — not what I proposed

My proposal was "put the place coordinates into one number field before the group
law". That turned out to be the right *area* and the wrong *fix*: the coordinates
are already fine (`PointsOver` RootReduces them, and `Solve` returns radicals).
Three stacked defects, each silent:

1. **`MinimalPolynomial` refused on a `Root[]` whose numerics do not converge.**
   It picks its irreducible factor by an 80-digit NUMERIC test with a 1e-3
   ceiling — run even with only ONE candidate, where `G(s)=0` holds by
   construction and there is nothing to choose. `N[Root[<dense deg 8, 40-digit
   coefficients>, 3], 80]` returns the Root unevaluated (Root::conv x3), so a
   provably correct answer was refused. Fixed in `src/poly/minpoly.c`; also
   10-20x faster on the single-candidate case, which is the common one.
   The silent chain: Root refinement declines -> MinimalPolynomial unevaluated ->
   `NontorsionDivisor` CoefficientLists it into a one-element list of garbage ->
   the mod-p certificate grinds past the whole budget. No link reported anything.

2. **`DivisionPolyOrder` mixed the two coordinates in its zero test.** Even `n`
   (which carry the `yv` factor) cost 0.003 ... 0.249 s and then hung; odd `n`,
   same degrees up to 264 in X but no `Y0`, stayed under 0.035 s. A product is
   zero iff a factor is.

3. **`TorsionOrder` ran the chord-and-tangent fallback after the division
   polynomials had already decided.** The psi_n criterion is a theorem, not a
   hint; the fallback belongs to the undepressed model only. It was the
   expensive half: 0.005, 0.019, 11.2 s, then nothing.

### #110 — pinned, not fixed

Its divisor is genuinely 4-torsion, so the certificate correctly declines and the
Miller logands get built over Q(sqrt2). `TDiv[T, TowerD[T,u], u]` then costs 5 s
per logand and `Can` GROWS the remainder (81 -> 163 -> 703 -> 879 leaves), because

    Together[x (1/(1 + Sqrt[2] t) + 1/(1 - t))]   BAILS, returns its input
    Together[x (1/(1 + 2 t)       + 1/(1 - t))]   combines
    Together[x (1/(1 + I t)       + 1/(1 - t))]   combines

i.e. the open `exact_poly_div` defect at `src/poly/poly.c:1450`, whose unit test
accepts Integer/Rational/Complex[q,q] and rejects everything else. It is also
exactly why #111 (roots ±I) is fast and #110 (roots 1±Sqrt[2]) is not. Blast
radius is PolynomialGCD/LCM/Mod, Decompose, Cancel and the exact linear solver,
so it is its own change with its own corpus run.

### Regression gates (all clean)

| gate | before | after |
|---|---|---|
| 312-case corpus | 303/9/0/0 | **304/8/0/0**, only #109 changed |
| PMT 114 | 94 solved / 3 declined / 17 nonelem | identical, 94 verified |
| Charlwood 50 | 49 of 50 verified | identical |
| unit suite | 5 red of 513 | identical set, crc_corpus byte-identical |

---

## Review — the bugs this session uncovered (v0.250, v0.251)

### Verified-fresh inventory first

Six candidates from earlier phases; **one had already been fixed** by other work
(`TrigToExp[Sin[x^2]]` -> ComplexInfinity and `TrigToExp[Log[x] Sin[x]]` ->
Indeterminate both now answer correctly). Re-checking before acting is the
cheapest step in the whole exercise.

### Fixed — three silent wrong answers (v0.250)

The class that matters: a confident wrong value rather than a decline.

1. **Number theory on a hard composite.** `FactorInteger`'s Automatic method is
   bounded and returns an unfactored cofactor with exponent 1. `df_factor_mpz`
   (`src/numbertheory/nt_gaussian.c`) read that back as a prime. For an 82-digit
   semiprime: `PrimeNu` = 1 and `PrimeOmega` = 1 (true 2), `MoebiusMu` = -1
   (true 1) -- while `PrimeQ` on the same number correctly said False.
   `facint.c` already stated the rule its own API follows; this helper never saw
   the flag. 40 Miller-Rabin rounds per base; all five consumers already
   declined on failure.
2. **`Factor` / `FactorList` ignored `Modulus`.** `Factor[x^2+1, Modulus -> 5]`
   gave `1 + x^2`; over GF(5) it is `(x+2)(x+3)`. New
   `flint_nmod_poly_factor_list` over FLINT's `nmod_poly_factor` -- already
   linked and already used inside `flint_bridge.c`.
3. **`PolynomialGCD` ignored `Modulus`**, while `PolynomialExtendedGCD` with the
   same option was already correct. Routed to the same FLINT nmod path.

All three **decline** where the fast path cannot go rather than falling back to
the answer over Q.

### Fixed — the last capped case (v0.251), and a wrong diagnosis on the way

I was confident #110 was the `Together` bail over Q(Sqrt[2]). That bug is real
and now fixed -- `exact_poly_div` accepts an `AlgebraicNumber` divisor, whose
inverse is COMPUTED rather than asserted, unlike the symbolic radical the
surrounding comment rules out -- and **it changed #110's time not at all.**

The real cause, found by timing the one call that did not return: `Can`'s
number-field detour. On #110's logand norm (164 leaves, atoms {Sqrt[2]}) the
direct path is Together 1.3 ms -> 56 leaves; the detour is Together 2.7 s ->
6732 leaves and then Cancel never finishes. The detour's only value is sparing
`CanRaw` the `Extension -> Automatic` option, which `CanRaw` asks for only when
a Complex atom is present. Gated on that.

**Corpus 305 PASS / 7 HONEST / 0 WEAK / 0 FAIL -- the parity ceiling.** #110 the
only verdict change; verdict-run total 282.9 -> 41.6 s.

### Gates

| gate | result |
|---|---|
| 312 corpus | **305/7/0/0**, one verdict change (#110), `why` unchanged on 311 |
| PMT 114 | 94 solved / 3 declined / 17 nonelem, all 94 verified -- baseline |
| Charlwood 50 | 49 of 50 verified -- baseline |
| unit suite | running at time of writing; `solve_corpus_tests` red in the batch run was LOAD (130/130 in isolation) |

### Still open

- `x^2 + 2x - 1` on the same curve still reaches the budget (not a corpus case).
- `Together` over a *symbolic* radical denominator still declines, by design.
- `Root` numeric refinement does not always converge (`src/root_numeric.c:699`).
- A branch error in `Integrate[Sqrt[1+x^3], x]` (residual -1.178 at x = 1/3),
  pre-existing and not in the corpus.
- `Hash[]` is still not a builtin. When it lands, fix `Part2`'s memo FIRST: its
  key is `Hash[{...}]`, which currently stays unevaluated and therefore compares
  structurally (exact). An integer hash would make that memo collision-prone.
