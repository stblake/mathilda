---
ticket: RG-4
created: 2026-08-31T15:03:08Z
researcher: Michael Sollami
source_sha: 854419997e4b0d8b35fb401783a17b13db2495c9
branch: find-vertex-coloring
repository: mathilda
topic: "Aborting a computation strands static depth counters; Integrate then silently returns unevaluated"
tags: [research, bug, correctness, timeconstrained, abort, calculus, integrate]
subsystems: [core, calculus]
type: research
lifecycle: active
status: complete
severity: correctness — silent wrong answer
last_updated: 2026-08-31
last_updated_by: Michael Sollami
---

# Research: Abort strands static depth counters, silently breaking `Integrate`

**Date**: 2026-08-31T15:03:08Z · **Commit**: `85441999` · **Branch**: `find-vertex-coloring`

> **Split from RG-3** on 2026-08-31 at the human's direction. RG-3 asks *which builtins should
> poll `tc_check_deadline` and do not* — a latency question that only bites on WSL 1. This is
> a **live correctness bug on every platform**, found while verifying RG-3's premise, and it
> outranks the audit. Different fix, different urgency, so: different ticket.

## TL;DR

After any `TimeConstrained` abort that lands inside `Integrate`, every later hard integral in
that session returns unevaluated instantly instead of being attempted. Silent wrong answer —
indistinguishable from "no closed form exists". Cause is a `siglongjmp` skipping the decrement
of a static recursion counter. Reproduced deterministically on an unmodified build.

## Summary

`TimeConstrained` aborts by `siglongjmp` (`src/core.c:3924-3930`) with no unwind cleanup — the
codebase documents this for *memory* (`src/core.c:3792-3796`: "the leak is bounded … the next
session-level GC is `Quit[]`"). That rationale does not extend to **static state**, which no
`Quit[]` short of process exit resets, and `src/calculus/` carries 69 unbalanced
`_depth++`/`_depth--` pairs across 12 files. An abort between an increment and its matching
decrement strands the counter above its cap for the life of the process, after which the guard
that reads it refuses all subsequent work.

The user-visible result is the worst available failure mode: not a crash, not a slowdown, but
`Integrate[...]` coming back unevaluated when the engine could have integrated it.

Contrast that isolates the pattern: `GroebnerBasis` aborts and then answers correctly. It
leaks heap on abort but holds no static counters. **The heads that leak memory stay correct;
the head that stays memory-bounded goes silently wrong.**

## Reproduction

Deterministic, on an unmodified `make -j8` build, macOS. Let
A = `Integrate[Sqrt[1+Sqrt[1+Sqrt[1+x]]], x]` and B = the same with `Sqrt[2+x]`.

```
FRESH session:  TimeConstrained[B, 12]              -> aborts after the full 12 s (it tries)
SAME session:   TimeConstrained[A, 2]  -> $Aborted
                TimeConstrained[B, 12]              -> Integrate, instantly (run total 2.00 s)
```

Timing of the abort is irrelevant — 0.15 s, 0.5 s, 1 s and 3 s all break it identically.
Simple integrals keep working (`Integrate[x^2,x]`, `Integrate[x^7 Sin[x],x]` are still
correct), because they are answered before the guarded recursion is entered.

**Controls run:**

| Control | Result |
|---|---|
| `TimeConstrained` on a *fast* integral (no abort), then B | B still broken? **No** — see caveat |
| Hard integral twice, no `TimeConstrained` at all | establishes the baseline for each integrand |
| `Limit`, `FullSimplify`, `D`, `Solve`, `Expand` after the abort | all unaffected |
| `GroebnerBasis` aborted, then re-run | fully recovers |

Caveat on scope, stated because an earlier draft got it wrong: the blast radius is
demonstrated for B, which shares A's structure. A third integrand initially cited as
"unrelated and also broken" turned out to return unevaluated on a *first* call in a clean
session — Mathilda cannot integrate it at all, so it proved nothing. Whether a structurally
unrelated integrand is affected is **not established**; that control is missing, not negative.

## Suspected mechanism

`src/calculus/integrate_unknown.c:705-718`, in full — the elided middle matters:

```c
Expr* key = canon(expr_copy(f));                                  /* :706        */
for (int i = 0; i < iu_depth; i++) {                              /* :708 CYCLE  */
    if (expr_eq(iu_stack[i], key)) { expr_free(key); return NULL; }/* :709 GUARD */
}
if (iu_depth >= IU_MAX_DEPTH) { expr_free(key); return NULL; }    /* :711 depth  */
iu_stack[iu_depth++] = key;                                       /* :712 enter  */
    ... iu_integrate_core(f, x) ...
iu_depth--;                                                       /* :716 exit   */
expr_free(iu_stack[iu_depth]);                                    /* :717        */
```

An abort between `:712` and `:716` never decrements, so `iu_depth` stays at *d* > 0 **and
`iu_stack[0..d-1]` keeps stale keys from the abandoned computation**.

**It is the cycle guard at `:708-709`, not the depth cap at `:711`, that does the damage.**
This distinction was supplied by review and it matters: the depth cap needs `iu_depth >= 64`
before it refuses anything, which a 0.15 s abort is unlikely to reach — the observation that
made the first draft's explanation implausible. The cycle guard needs only **one** stranded
key that `expr_eq`-matches a subproblem of the new integrand. Since it is designed to detect
recursion cycles, a stale key silently becomes a permanent "already tried this, give up".

This predicts exactly the blast radius measured: integrand B shares A's structure, so its
recursion regenerates a subproblem matching a stranded key and is refused; a structurally
unrelated integrand would *not* match and should be unaffected — which is why the missing
control noted above is the one worth running first.

**Ruled out: the failure memo.** `intg_fail_tab` / `intg_fail_count`
(`src/calculus/integrate.c:722`) is a stranded static of the same family and the obvious
competing explanation, but `intg_fail_sync_epoch` (`:725-728`) zeroes the table whenever
`eval_toplevel_id` changes, so it cannot persist across top-level commands. Verified.

**Confidence: symptom reproduced deterministically; mechanism inferred from code, not
instrumented.** Printing `iu_depth` and dumping `iu_stack` after an abort would settle it in
minutes. That is the first task of this ticket, not a conclusion of it.

## Scope of the pattern

Counts corrected after review; the first draft reported "69 unbalanced sites across 12 files,
one counter per file", which was wrong three ways — see the note below the table.

| Area | `_depth++` (enter) | `_depth--` (exit) | Distinct counters |
|---|---|---|---|
| `src/calculus/` | 21 | 48 | **13** across 12 files |
| `src/simp/` | — | — | 1 (`simp_log.c`, 2 sites) |
| root `src/*.c` | — | — | 13 sites in `eval.c`, `match.c`, `parse.c`, `print.c`, `rat.c`, `random.c`, `message.c` |

The 12 calculus files are `integrate.c`, `integrate_unknown.c`, `integrate_chebychev.c`,
`integrate_derivdivides.c`, `integrate_goursat.c`, `integrate_diffunderint.c`,
`integrate_linrad.c`, `integrate_jeffrey.c`, `integrate_linratiorad.c`,
`integrate_quadrad.c`, `gruntz.c`, `dsolve.c`.

**Three corrections, because the fix sizing rested on the wrong ones.** (1) "69" was the total
of increment *and* decrement lines, not distinct sites: it is 21 enter points and 48 exits.
(2) They are **not** unbalanced in source — they are correctly paired, and become unbalanced
only on the abort path; the original wording read like a pre-existing coding defect. (3) It is
**13** counters, not 12, and not one per file: `integrate.c` carries both `crc_depth` (`:319`)
and the non-static `g_integrate_depth` (`:696`).

The 48-vs-21 asymmetry is itself the important signal: multiple early-return decrement paths
per counter means the fix is **less mechanical** than "reset 12 things", and any save/restore
approach has more exit paths to get right than the first draft assumed.

`gruntz.c` (`Limit`) and `dsolve.c` (`DSolve`) carry the same pattern and are therefore
suspected-vulnerable by inspection. Measured: `Limit` is **unaffected by an `Integrate`
abort** — expected, since the counters are per-file and the abort landed in `Integrate`'s
recursion. Whether aborting *inside* `Limit` breaks `Limit` is untested.

## Candidate fixes (for the plan phase, not decided here)

1. **Reset the counters on the abort path.** After the `siglongjmp` lands in `tc_run_guarded`
   (`src/core.c:3924-3930`), nothing is in flight, so `0` is provably the correct value. One
   reset hook per file, called from one place. Smallest change.
2. **Save/restore around the `sigsetjmp`.** Safer under a nested `TimeConstrained` inside an
   already-running integration, where the correct restore value is nonzero rather than 0.
   Touches 12 files instead of one.
3. **A general cleanup registry.** The complete answer, and what
   `docs/design/timeconstrained-abort-channel.md` gestures at. Much larger — an unwind-safe
   allocation and state discipline tree-wide. Not justified by the instances known today.

Option 1 is the obvious starting point; option 2 is the honest one if nested
`TimeConstrained` inside `Integrate` is reachable. Determining that is an open question.

## Open Questions

### Unresolved

- [ ] Which static is actually responsible — `iu_depth`, another `*_depth`, or
      `intg_fail_count`? Needs instrumentation. First task of the fix.
- [ ] Reset-to-zero or save/restore? Depends on whether a nested `TimeConstrained` can be
      live inside a running integration.
- [ ] Is a structurally unrelated integrand affected? The missing control above.
- [ ] Does aborting inside `Limit`/`DSolve` break `Limit`/`DSolve` the same way? Same pattern
      by inspection, untested.

### Resolved

- [x] Should this be fixed ahead of RG-2 Phase 2, or in parallel? — **Moot, 2026-08-31.**
      RG-2 is stopping at Phase 2 and submitting, so there is no RG-2 work left for this to
      be ahead of or in parallel with. RG-4 is sequenced on its own merits from here.

- [x] Is it a memo/cache rather than an abort effect? — No. B works for a full 12 s in a fresh
      session and dies instantly after an unrelated abort. _(measured)_
- [x] Is it caused by the head *polling* `tc_check_deadline`? — No. No file in
      `src/calculus/` polls at all; the abort arrived via SIGPROF. This bug is entirely
      independent of RG-3's polling question and is live on healthy hosts today. _(verified)_
- [x] Is the damage global? — No, scoped to the subsystem executing at abort time. _(measured)_
- [x] Do heads that already poll recover? — Yes, `GroebnerBasis` does. _(measured)_

## Requires Approval

**Whether this jumps the queue.** It is a silent-wrong-answer bug with a known trigger and no
workaround short of restarting the session, which argues for fixing it before more feature
work. It is also not a regression — it has presumably been present since `TimeConstrained`
landed — which argues it can wait. That is a judgement call about risk appetite, not a
technical question, and it is why this ticket exists separately from RG-3.

## References

- Split from `thoughts/shared/tickets/RG-3/research.md`, which holds the full
  `TimeConstrained` mechanism analysis and the polling audit this was separated from.
- `docs/design/timeconstrained-abort-channel.md` — the decision record; its
  "the leak is bounded by … `Quit[]`" reasoning is what this finding shows to be incomplete.
- `src/core.c:3792-3796`, `:3924-3930` — the documented abort contract and the jump target.
