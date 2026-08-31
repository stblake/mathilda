---
ticket: RG-3
created: 2026-08-31T15:03:08Z
researcher: Michael Sollami
source_sha: 854419997e4b0d8b35fb401783a17b13db2495c9
branch: find-vertex-coloring
repository: mathilda
topic: "Which builtins must poll tc_check_deadline and do not"
tags: [research, codebase, timeconstrained, abort, calculus, integrate, poly, eval]
subsystems: [core, eval, calculus, poly, graph]
type: research
lifecycle: active
status: complete
last_updated: 2026-08-31
last_updated_by: Michael Sollami
---

# Research: Which builtins must poll `tc_check_deadline` and do not

**Date**: 2026-08-31T15:03:08Z
**Researcher**: Michael Sollami
**Git Commit**: `854419997e4b0d8b35fb401783a17b13db2495c9`
**Branch**: `find-vertex-coloring`
**Repository**: mathilda

## TL;DR

Asked which builtins should poll `tc_check_deadline()` and don't. Found that just 2 of 832
registered heads poll — but that the polling gap is the *smaller* problem. Aborting a long
computation strands static recursion counters, and `Integrate` thereafter silently returns
unevaluated for every hard integrand for the rest of the session. Reproduced on an
unmodified build. Fix the corruption before widening the polling.

## Summary

`TimeConstrained` aborts via `SIGPROF` + `siglongjmp`, with `tc_check_deadline()` as a
cooperative backstop that only helps where the signal is unreliable (WSL 1, which the source
says under-counts CPU time and delivers SIGPROF ~15× late, `src/core.c:3778-3786`). The
evaluator polls it once per *rewrite step* (`src/eval.c:2306`), which a long-running C
builtin never reaches — hence the four hand-placed polls in `poly/groebner.c`,
`poly/gbmod.c`, `poly/groebnerwalk.c` and `graph/vertexcoloring.c`.

So the literal answer: **any builtin with a long C-level loop and no poll is uninterruptible
on WSL 1**. On healthy hosts the signal catches them, which I verified — `FullSimplify` and
`Integrate` both abort correctly here despite no calculus file polling anything.

But that verification surfaced something worse, **now split out as
[`RG-4`](../RG-4/research.md)** at the human's direction: aborting strands static recursion
counters in `src/calculus/`, after which `Integrate` silently returns unevaluated for
structurally-related hard integrands for the rest of the session. A live correctness bug on
every platform, independent of the polling question — the abort that triggers it arrives via
SIGPROF, and no calculus file polls anything. RG-4 carries the reproduction, mechanism and
fixes; only its consequences *for this audit* are kept below.

## Open Questions

### Unresolved

- [ ] Should the `src/calculus/` polling work wait on RG-4, or is excluding that subsystem
      from the first pass sufficient? (Narrowed from a broader sequencing question after
      review — see `## Research Review`.)
- [ ] Do the top-ranked polling targets (`src/solve/`, `src/strings/regex/`,
      `src/poly/mvfactor.c`) hold stranded statics of their own? Unchecked, and worth knowing
      before the first poll lands in any of them.

### Resolved

- [x] Should the state-corruption bug and the polling audit be one ticket or two? — **Two.**
      Split on 2026-08-31 at the human's direction: the corruption is a live correctness bug on
      every platform, the audit is a WSL-1 latency question. Filed as
      `thoughts/shared/tickets/RG-4/research.md`, which carries the reproduction, the
      mechanism and the candidate fixes. RG-3 keeps the audit. _(human decision)_
- [x] Reset counters to 0 or save/restore around the `sigsetjmp`? — Moved to RG-4, where the
      fix lives. Still open there. _(scope move)_
- [x] Does `TimeConstrained` actually interrupt pure-C builtins? — Yes, on this host, via the
      signal path. `FactorInteger` on a hard 50-digit semiprime, `FullSimplify`, and
      `Integrate` all return `$Aborted`/failexpr with the session intact. _(measured)_
- [x] Is the polling gap detectable by testing? — **No**, and this is the methodological
      finding: on any host with working `ITIMER_PROF` the signal catches everything, so an
      empirical abortability test passes for polling and non-polling heads alike. The audit
      must be static (source-based). _(derived from the measurements plus `src/core.c:3778`)_
- [x] Do the heads that already poll recover cleanly? — Yes. `GroebnerBasis` aborted and then
      answered correctly. It leaks heap but holds no static state. _(measured)_
- [x] Is the corruption global or scoped? — Scoped to the subsystem executing at abort time.
      An aborted `Integrate` does not damage `Limit`, `FullSimplify`, `D`, or `Solve`.
      _(measured)_

## Research Review

Reviewed 2026-08-31 by `plan-reviewer` (coverage-gap + contradiction lenses), which verified
the source claims independently. Four Blocking, five Worth Flagging. All resolved before this
document was presented; the reviewer read a pre-correction draft, so two of its findings had
already been fixed in flight and are recorded as such.

### Blocking

_None._

### Worth Flagging

**[WORTH FLAGGING] Shortlist rows 2, 4 and 6-12 are unverified survey output presented in the same voice as the three I checked**
- Where: the ranked shortlist table
- Why: line numbers and constant values give every row equal apparent authority; a reader
  cannot tell which nine of twelve to re-check.
- Addressed: ✓/~ provenance markers per row plus a paragraph stating what each means. Left as
  Worth Flagging rather than Resolved because the rows are still unverified — the marking
  makes that legible, it does not fix it.

**[WORTH FLAGGING] "Reset to 0 is provably correct" is asserted in the summary while the research doc lists it as an open question**
- Where: `research-summary.md` Decision Criteria vs this doc's `### Unresolved`
- Why: no proof is given, and the nested-`TimeConstrained` counterexample is a reachable state
  this document's own measurement table exercises.
- Addressed: "provably" dropped from the summary; both artifacts now defer to the open
  question. Carried into RG-4, where the fix actually lives.

### Resolved

**[BLOCKING] The stated mechanism does not account for the symptom — a stranded `iu_depth` cannot refuse unrelated integrals**
- Why: the depth cap at `integrate_unknown.c:711` only refuses once `iu_depth >= 64`, which a
  0.15 s abort will not reach; and the draft claimed an *unrelated* integrand was broken.
- Resolved: two ways, and the finding improved the result. (a) The "unrelated integrand"
  evidence was already withdrawn in flight — a control showed that integrand returns
  unevaluated on a first call in a clean session, so it proved nothing. (b) The reviewer
  identified the **cycle guard** at `:708-709` as the real mechanism: it needs only one stale
  `iu_stack` key, not depth 64, and it predicts exactly the blast radius measured. Both now
  documented in `RG-4`.

**[BLOCKING] The alternative explanation is never raised and the confidence language overstates the evidence**
- Why: written as established mechanism with no hedge; the competing `intg_fail_tab`
  explanation unmentioned.
- Resolved: confidence now stated explicitly ("symptom reproduced; mechanism inferred, not
  instrumented"), and `intg_fail_tab` is recorded as **ruled out** with the reason —
  `intg_fail_sync_epoch` (`integrate.c:725-728`) zeroes the table per top-level eval, so it
  cannot persist across commands. Verified.

**[BLOCKING] "4 of 754 builtins" is 2 of 832**
- Why: the four `poly/` sites are all internals of one registered head, `GroebnerBasis`; the
  table counted files and the prose relabelled them builtins. The denominator grep's character
  class dropped names.
- Resolved: corrected to 2 heads / 5 call sites / 832 builtins throughout, including the
  summary's ratio.

**[BLOCKING] "69 unbalanced sites across 12 files, one counter per file" is wrong three ways**
- Why: 69 is increments *plus* decrements (21 and 48); they are correctly paired in source and
  unbalanced only on the abort path; and there are 13 counters, not 12 —`integrate.c` holds
  both `crc_depth` and `g_integrate_depth`.
- Resolved: recounted in RG-4, wording fixed, and the 48-vs-21 asymmetry now called out as
  evidence the fix is *less* mechanical than the sizing claimed.

**[WORTH FLAGGING] "~55 local cap constants" is not reproducible from the stated method**
- Resolved: the sweep yields exactly **41**; the number and the exact command are now in the
  text, described as a floor rather than an estimate.

**[WORTH FLAGGING] The reorder recommendation overreaches**
- Why: the corruption surface is `src/calculus/`, but none of the top four polling candidates
  lives there — the sequencing constraint applies to 1 of 12 rows, not the audit as a whole.
- Resolved: narrowed to "do not add polls inside `src/calculus/` until the counters are
  fixed"; the rest of the audit can proceed in parallel.

**[WORTH FLAGGING] The two documents' recommendations differ — split is an open question here and a delivered recommendation in the summary**
- Resolved: the human decided the split on 2026-08-31. Both artifacts now record it as
  decided, and RG-4 exists.

## Requires Approval

Two scope calls. **This research found a live bug outside its own question** — the
`Integrate` state corruption — and recorded it rather than fixing it, per the read-only
posture. Someone must decide whether that becomes RG-3's primary scope, a separate ticket, or
a hotfix ahead of both. **And the polling audit's value is host-dependent**: on Linux/macOS
CI it changes nothing observable, so the case for it rests on WSL 1 users and on the
`Abort[]` finding below, not on a measurable defect here.

---

## Research Question

> which builtins must poll tc_check_deadline and do not

**Process note.** The `grill-me` `research-open` interview was **skipped**, deliberately:
scope, constraints and the reason for the ticket were pre-specified by the human in the
conversation that opened it (the RG-2 Phase 1 acceptance message, which named the question,
its origin in `docs/design/timeconstrained-abort-channel.md`, and that it should be its own
ticket). Recorded rather than left silent — a gate that did not apply is otherwise
indistinguishable from one that did not run. The `research-close` pass was not skipped.

## Detailed Findings

### The mechanism, and why polling exists at all

`TimeConstrained` (`src/core.c:3984-4131`) races two layers into one `siglongjmp` target
(`tc_run_guarded`, `src/core.c:3924-3930`):

- **`SIGPROF` + `ITIMER_PROF`** — asynchronous, CPU-time-based, can interrupt any
  instruction including mid-builtin.
- **`tc_check_deadline()`** (`src/core.c:3889-3902`, declared `src/core.h:88`) — cooperative,
  wall-clock, fires *only* where called by hand. A no-op single branch when no
  `TimeConstrained` is active (`tc_deadline_active`, `src/core.c:3829`).

The source states the division plainly (`src/core.c:3778-3786`):

> "on real Linux and macOS, ITIMER_PROF is precise and the cooperative check is a cheap
> no-op. On hosts that mishandle ITIMER_PROF -- notably WSL 1, whose syscall-translation
> layer under-counts CPU time and delivers SIGPROF ~15x late -- the cooperative wall-clock
> backstop catches the deadline at the next rewrite step. The only case that escapes both
> layers is a single long-running C builtin (e.g. FactorInteger on a huge composite), which
> cannot be aborted cooperatively and must wait for the late SIGPROF on broken hosts."

`src/core.h:84-87` repeats it: "Granularity: limited to between rewrite steps."

### Who polls: 2 registered heads, 5 call sites, of 832 builtins

`grep -rln tc_check_deadline src/` → 7 files, of which `core.c` defines it, `core.h` declares
it, and `eval.c` is the rewrite loop. **Four builtin-side pollers:**

| Site | Function | Frequency |
|---|---|---|
| `src/poly/groebner.c:1025` | `gb_buchberger` | per pair dequeued |
| `src/poly/groebner.c:1234` | `gb_strong_buchberger` | per pair |
| `src/poly/groebnerwalk.c:366` | `gb_groebner_walk` | per walk step |
| `src/poly/gbmod.c:429` | `gfp_buchberger` | per pair |
| `src/graph/vertexcoloring.c:270` | `fvc_bb` | every 4096 nodes |

**Corrected after review — the first draft said "4 of 754" and both numbers were wrong.**
Those five call sites belong to just **two registered heads**: all four `poly/` sites are
internals of `GroebnerBasis` (`src/poly/groebnerbasis.c:920`), and `fvc_bb` serves
`FindVertexColoring`. The table counts *files*; the prose had relabelled them "builtins". And
the denominator is **832**, not 754 — my grep's character class dropped names the reviewer's
`[^"]*` catches. So: **2 heads / 5 call sites / 832 registered builtins.**

`src/eval.c:2306` is the evaluator's poll, at the top of the per-rewrite-step loop — so it
sees a builtin only *between* invocations, never inside one.

### The polling gap is real but host-dependent — and untestable

Measured on this machine (macOS, working `ITIMER_PROF`), all via `-file` scripts:

| Expression | Bound | Result |
|---|---|---|
| `FactorInteger[p*q]`, p,q ~25 digits | 3 s | `$Aborted`, session intact |
| `FullSimplify[Sum[Sin[k x]Cos[k y]/(k^2+1),{k,1,12}]]` | 5 s | aborted (failexpr) |
| `Integrate[Sqrt[1+Sqrt[1+Sqrt[1+x]]],x]` | 5 s | aborted (failexpr) |
| `Do[x=i^2,{i,10^9}]` | 2 s | `$Aborted` |
| nested `TimeConstrained[…, 30]` inside `[…, 3]` | 3 s | outer fired correctly |

**None of `src/calculus/` polls `tc_check_deadline`** (verified by grep), so every one of
those aborts came from the signal alone. Consequence for the audit: *you cannot find the gap
by testing.* On a healthy host, polling and non-polling heads behave identically. Only a
source-level audit — or a WSL 1 machine — distinguishes them.

### The abort-safety bug found while verifying this — now RG-4

Split to [`thoughts/shared/tickets/RG-4/research.md`](../RG-4/research.md) on 2026-08-31 so
the numbers live in exactly one place. In brief: a `siglongjmp` skips the decrement of a
static recursion counter in `src/calculus/`, stranding stale keys in the cycle-guard stack at
`src/calculus/integrate_unknown.c:708-709`, after which any integrand matching one is refused
instantly. Reproduced deterministically. `GroebnerBasis` — which does poll — recovers fully,
because it leaks heap but holds no static counters.

**What it means for THIS audit**, and the only reason it stays in this document:

- The corruption surface is `src/calculus/`, so **do not add polls there until RG-4 lands**.
  A new poll is a new deterministic way to reach the bug.
- It does **not** block the rest of the audit: none of the shortlist's top four candidates is
  in `src/calculus/`. Only row 5 (`intrischnorman.c`) is.
- Whether the other high-ranked targets hold stranded statics of their own is **unchecked**,
  and is the thing to establish before the first poll lands in any of them.

### The non-polling shortlist, ranked

A survey of the subsystems where unbounded search is most likely, ranked by how far the work
can run inside one builtin invocation. The absence of `tc_check_deadline` is verified
tree-wide (only 5 call sites exist, all listed above), so that part holds for every row.

**Provenance — rows differ in how well checked they are.** Rows 1, 3 and 5 I verified against
source myself (marked ✓). Rows 2, 4 and 6-12 come from a survey agent and are **unverified**
(marked ~): their line numbers and constant values are reported as given, and should be
re-checked before any of them is acted on.

| # | Head(s) | Unbounded construct | Existing local guard |
|---|---|---|---|
| 1 ✓ | `Solve[…, Integers]` leaf search | `src/solve/solveint_leaf.c:179`, `:486` — recursion branching on the search-box size | `st->max_visits`, normally `SI_MAX_NODES` = 2×10⁸, **raised to 3×10⁹** at `:457` |
| 2 ~ | `Solve[…, Integers]` meet-in-the-middle | `src/solve/solveint_mitm.c:145`, `:167` — two `for(;;)` odometers | `MITM_HASH_CAP` = 5×10⁶ (`:46`) |
| 3 ✓ | `StringMatchQ`/`StringCases`/`StringReplace`/… on `RegularExpression` | `src/strings/regex/regex_engine.c:98` — the single `pcre2_match` site | **none** — zero `match_limit`/`depth_limit` calls in `src/strings/regex/` (verified) |
| 4 ~ | `Factor`/`FactorList` multivariate | `src/poly/mvfactor.c:634-676` — `while (mask < end)` over `2^r` subsets, Hensel lift per iteration | none beyond the combinatorial `2^r` |
| 5 ✓ | `Integrate` Risch-Norman fallback | `src/calculus/intrischnorman.c` search loops (`risch_coupled.c:190,252,362`, …) | `PMINT_BUDGET_SEC` 4.0 — see below |
| 6 ~ | `Limit` (Gruntz) | `src/calculus/gruntz.c` recursive scale construction | `GRUNTZ_MAX_DEPTH` 80, `GRUNTZ_MAX_WORK` 6000 |
| 7 ~ | `NIntegrate` | `src/numerical_calculus/gkadapt.c:174`, `cubature.c:203` — adaptive subdivision, symbolic integrand per node | `MaxRecursion` (default 800) |
| 8 ~ | `NMinimize`/`FindMinimum` global methods | `nm_dual_annealing.c:134`, `nm_basin_hopping.c:114`, `nm_shgo.c:87` | `NM_DA_MAXFUN` 10⁷, `SHGO_MAX_VERTICES` 8000 |
| 9 ~ | `FindRoot`/`NRoots` | `src/numerical_roots/findroot.c` — seven MPFR iteration loops; `nroots.c:107` | `max_iter` 100, or `100 + 20d` |
| 10 ~ | `Solve[…, Integers]` Thue | `src/solve/solvethue.c:616`, `:640` — LLL rounds | 30 rounds / 2000 iterations |
| 11 ~ | `FindClusters` | `src/list/find_clusters.c` | `FC_MAX_ITER` 100, `FC_LLOYD_MAX_WORK` 2×10⁷ |
| 12 ~ | `GaussianMixture` | `src/ml/gmm.c:171` | `ML_GMM_MAX_ITER` 200 |

**The ratio that sizes the problem: at least 41 local cap constants against 5 poll sites in 2 heads.** A
`#define` sweep for `MAX_ITER|MAXITER|MAX_STEPS|MAX_DEPTH|MAX_NODES|MAX_WORK|MAXFUN|…` finds
exactly **41** distinct constants (corrected from a hand-waved "~55"; the sweep is
`grep -rhoE '#define [A-Z_]*(MAX_ITER|MAXITER|MAX_STEPS|MAX_DEPTH|MAX_NODES|MAX_WORK|MAXFUN)[A-Z_]*' src/ | awk '{print $2}' | sort -u | wc -l`).
It genuinely undercounts — `SI_MAX_NODES`, `MITM_HASH_CAP`, `PMINT_BUDGET_SEC` and
`GB_WALK_MAX_STEPS` use non-matching suffixes and were found separately — but by an unmeasured
amount, so treat 41 as a floor rather than an estimate. Each is
an author who knew their code could run away and reached for a private ceiling instead of the
shared abort channel. That is the pattern `docs/design/timeconstrained-abort-channel.md`
rules against, and it is already the tree's dominant habit by roughly twenty to one by head.

**Two entries deserve singling out.**

`src/calculus/intrischnorman.c` already has a **wall-clock** budget — `PMINT_BUDGET_SEC 4.0`
(`:79`), armed with `clock()` at `:3637` — which is precisely the machine-dependent design
the ADR rejects: the same integral succeeds on a fast host and fails on a slow one. Worse, it
is consulted **once**, as a pre-flight gate before `RowReduce` (`:3087`), not polled inside
the expensive loops, so it cannot stop a long run that has already started. Verified.

The regex path (#3) is the one case where cooperation may be **impossible in principle**:
`pcre2_match` (`regex_engine.c:98`) runs entirely inside PCRE2's own C stack with no callback
reachable from Mathilda. There are no `pcre2_set_match_limit`/`set_depth_limit` calls at all
(verified: zero hits). For this head the remedy is PCRE2's own limits, not a poll.

### Why this inverts the ticket's priority

Widening `tc_check_deadline` polling adds *deterministic* abort points to code that today is
only interrupted by an asynchronous signal. Every new poll in a function holding a static
counter is a new, reliably-reachable corruption site. The corruption is already live on
healthy hosts — my reproduction used the signal path, not a poll — so this is not a
regression the audit would introduce, but the audit would make it easier to hit.

**Narrowed after review.** The strong form — "make abort safe before making anything else
abortable" — overreaches. The corruption surface is `src/calculus/`, and *none* of the
shortlist's top four polling candidates lives there (`src/solve/`, `src/strings/regex/`,
`src/poly/mvfactor.c`). Only row 5, `intrischnorman.c`, is in the corrupting subsystem. So the
defensible constraint is the cheaper one: **do not add polls inside `src/calculus/` until the
counters are fixed.** The rest of the audit can proceed in parallel. Whether the top-ranked
targets hold stranded statics of their own is unchecked, and is worth checking before the
first poll lands in any of them.

### No cleanup mechanism exists

Confirmed absent tree-wide: no arena, registry, or `TC_CLEANUP`-style hook. The codebase
states it (`src/core.c:3792-3796`): "the longjmp unwind cannot run destructors, so any Expr
nodes the in-flight evaluator allocated leak … the leak is bounded by the size of the
abandoned computation and the next session-level GC is `Quit[]`." The one mitigation is
narrower: `tc_gmp_alloc_busy` / `tc_alloc_safepoint` (`src/core.c:3802-3880`) defers a jump
that lands inside a guarded GMP allocation, preventing a libmalloc zone-lock **deadlock** —
it recovers no memory and no state.

Note the documented rationale covers *memory* ("bounded … `Quit[]`"). It does not cover
corrupted static state, which is unbounded in time and invisible.

### Incidental: `Abort[]` is unwired

`SYM_DollarAborted` is produced at exactly two sites, both in `builtin_time_constrained`
(`src/core.c:4020`, `:4123`) — `TimeConstrained` is the *sole* producer of `$Aborted`.
`Abort[]` (`SYM_Abort`, `src/sym_names.c:910`) exists only as a flow-control marker
recognised in `src/iter.c:267`, `src/funcprog.c:3081`, `src/core.c:1032-1036`. There is **no**
`symtab_add_builtin("Abort", …)`, no docstring, and nothing that rewrites a surviving
`Abort[]` to `$Aborted` — unlike the explicit uncaught-`Throw` reporting at
`src/eval.c:2323-2328`. So `Abort[]` typed by a user propagates as an inert expression.

## Code References

- `src/core.c:3889-3902` — `tc_check_deadline`, the cooperative poll
- `src/core.c:3778-3786` — why both layers exist; the WSL 1 quote
- `src/core.c:3792-3796` — the documented leak-on-abort rationale (memory only)
- `src/core.c:3924-3930` — `tc_run_guarded`, the single `sigsetjmp` target
- `src/core.c:3802-3880` — GMP allocator lock-safety deferral
- `src/core.c:4020`, `:4123` — the only two `$Aborted` producers
- `src/core.h:84-87` — "limited to between rewrite steps"
- `src/eval.c:2306` — the evaluator's once-per-rewrite-step poll
- `src/calculus/integrate_unknown.c:698-717` — the unbalanced counter, canonical instance
- `src/calculus/integrate.c:318-360` — `crc_depth` / `MAX_CRC_DEPTH`, same shape
- `src/poly/groebner.c:1025`, `:1234`; `src/poly/groebnerwalk.c:366`; `src/poly/gbmod.c:429`
  — the existing polls, and the precedent
- `src/graph/vertexcoloring.c:257-270` — RG-2's poll and its leak comment

## Architecture Insights

- **Two-layer abort with asymmetric coverage.** The signal is universal but host-fragile; the
  poll is reliable but only where hand-placed. Neither alone is sufficient, and the codebase
  knows this — the design is deliberate and documented, not accidental.
- **"Bounded by `Quit[]`" reasons about memory and silently generalises to state.** The
  comment at `src/core.c:3792` is correct about `Expr` nodes and wrong about static counters,
  which no `Quit[]` short of process exit resets. A correct-sounding local rationale that does
  not cover an adjacent case is exactly the failure mode SPEC.md §9 describes for audits.
- **Local guard constants are the tree's habit** (`IU_MAX_DEPTH`, `MAX_CRC_DEPTH`,
  `FM_MAX_CON`, `FVC_MAX_STEPS`). They bound work but are orthogonal to interruptibility, and
  the depth-counter variety is precisely what abort corrupts.
- **RG-2's `FVC_MAX_STEPS` design is validated by this.** `fvc_bb` keeps its counter in a
  heap `FvcBB` struct passed by pointer, not a file static, so an abort cannot strand it — the
  next call starts clean. Unintentional at the time, but the right pattern, and worth stating
  as the convention.

## Historical Context (from thoughts/)

- `thoughts/shared/tickets/RG-2/plan.md` — the ticket that surfaced this. Its `## Overview`
  carries the correction of the original false premise ("there is no timeout … no abort"), and
  its `## Risks and Rollback` records the abort-path leak for `FindVertexColoring`
  specifically. This research generalises that risk and finds it understated.
- `docs/design/timeconstrained-abort-channel.md` — the decision record that named this audit
  as its open follow-up. Its "Open follow-up" paragraph is what RG-3 exists to close.
- No prior research or plan in `thoughts/` mentions `tc_check_deadline` (grep: only RG-2).

## Related Research

- `thoughts/shared/tickets/RG-2/research.md` — graph subsystem conformance survey
