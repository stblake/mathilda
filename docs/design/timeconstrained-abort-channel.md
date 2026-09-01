# Long-running builtins cooperate with `TimeConstrained` rather than inventing local guards

*A decision record. Filed here, not under a numbered `docs/adr/`, because
`.claude/GUIDANCE_ROLES.md` resolves `architecture-guidance` to `docs/design` and this
directory's convention is descriptive filenames — a lone numbered file among them would be
a convention of one.*

- **Status:** accepted (narrow decision); the tree-wide audit in Consequences is **not done**
- **Date:** 2026-08-31
- **Participants:** Michael Sollami
- **Supersedes:** none
- **Superseded by:** —

## Context

`FindVertexColoring[g]` (ticket RG-2) computes a minimal vertex colouring, which is
NP-hard. Its plan needed a way to give up, and reasoned as follows:

> The subsystem's only error channel is `return NULL` (unevaluated) — there is no timeout,
> no iteration budget, no abort. Rather than invent one, the head refuses outright above a
> named vertex cap.

That premise was **false**, and it was the sole justification for the resulting
128-vertex cap. Established during Phase 1 implementation:

- `TimeConstrained` exists (`src/core.c:4038`) and is built on `sigaction(SIGPROF)` +
  `sigsetjmp`/`siglongjmp` + `setitimer(ITIMER_PROF)`, with a cooperative wall-clock
  deadline as a portability backstop for hosts where `ITIMER_PROF` is unreliable (WSL 1).
- **It interrupts pure-C builtins, not just the rewrite loop.** Verified empirically:
  `TimeConstrained[FactorInteger[<hard 50-digit semiprime>], 3]` returns `$Aborted` after
  3 s with the session intact. `FactorInteger` spends that time inside a C loop that never
  returns to the evaluator, so the signal path — not the cooperative one — is what stopped
  it.
- `tc_check_deadline()` is **exported** in `src/core.h:88`, precisely so a builtin that does
  not return to the evaluator's rewrite loop can poll the deadline itself.

Two further facts shaped the decision:

- **A size cap does not bound cost.** Measured on RG-2's search at edge density ≈0.24: 0.33 s
  at n=80, 30 s at n=100, and *unfinished after 60 s* at n=128 — a graph **under** the
  proposed cap. Sparse graphs at the same n=128 answer in 0.00 s. Cost tracks density, not
  vertex count, so a vertex cap is close to unrelated to the thing it was meant to guard.
- **There is no cleanup registry.** A `siglongjmp` out of a builtin leaks whatever that
  builtin had allocated. This is pre-existing and tree-wide (`FactorInteger` included), not
  specific to any one head.

## Decision

**`TimeConstrained` is the abort channel for long-running operations, and a long-running
builtin's obligation is to cooperate with it — by polling `tc_check_deadline()` on a cheap
interval — not to invent a local guard in its place.**

Where an operation can run unboundedly with nobody present to interrupt it (an unattended
script), it may carry a local budget *in addition*, on two conditions:

1. **It is a backstop, not the responsiveness mechanism.** Size it to bound the pathological
   case, not to make the typical case feel fast. A guard tightened until it converts correct
   answers into refusals has made the operation worse, not safer: RG-2's budget at 2M nodes
   refused a dense n=100 graph after 14 s that the search *answers* in 29 s.
2. **It is deterministic — a work count, not a wall clock.** A time-based cutoff makes the
   *answer* machine-dependent: a fast host proves minimality where a slow host refuses the
   same input. A node/iteration count gives every machine the same answer.

On exceeding either guard, return unevaluated. Never return a plausible-but-unproven result
in place of the guaranteed one — that converts a loud failure into a silent wrong answer.

## Consequences

**Easier:**
- Callers get one idiomatic, Wolfram-compatible way to bound any slow call, rather than a
  different per-head constant to discover for each one.
- New long-running builtins need three lines (`#include "core.h"`, a counter, a periodic
  `tc_check_deadline()`), not a bespoke guard design.
- Local caps can be set generously, because they are no longer the only thing standing
  between a user and a hang.

**Harder:**
- Every long-running builtin now carries a real obligation, and nothing enforces it. A head
  that fails to poll is silently uninterruptible on hosts where `SIGPROF` does not fire.
- Cooperating means accepting that a `siglongjmp` may unwind through the operation, leaking
  its scratch allocations. Correct behaviour here requires either a cleanup registry (does
  not exist) or accepting the leak (what the tree does today).

**Now foreclosed:**
- "This subsystem has no abort facility" is no longer a valid justification for a hard
  refusal. Any future cap has to argue against `TimeConstrained` on its merits.

**Open follow-up, deliberately not closed by this ADR:** *the tree has never been audited
for heads that should poll `tc_check_deadline()` and do not.* Against ~765 registered
builtins, `grep -rn tc_check_deadline src/` names six files: `core.c` (defines it), `eval.c`
(the rewrite loop), `poly/groebner.c`, `poly/gbmod.c`, `poly/groebnerwalk.c`, and
`graph/vertexcoloring.c`. So exactly **four** builtin-side pollers, three of them the
Gröbner family — which is worth noting as precedent: this decision codifies what that
subsystem already does rather than inventing a convention.

Everything else long-running is uninterruptible wherever `ITIMER_PROF` is unreliable. Note
what that does *not* mean: `FactorInteger` does not poll, yet the empirical test above shows
it aborting correctly, because the signal path caught it. Polling is the portability
backstop, so a non-polling head is a host-dependent gap, not a guaranteed hang — which is
precisely the kind of gap that stays invisible until someone runs the code on WSL.

Sizing the gap warrants its own ticket, and probably a `make` audit in the family of
`check-packed-aware`: that shape is known to work here, and SPEC.md §9's argument transfers
verbatim — audit silence means "never looked", not "checked and exempt".

## Alternatives considered

### A hard size cap alone (what RG-2 originally planned)
**Rejected because:** measurement showed it does not bound cost. `CompleteGraph[128]` and
`CycleGraph[128]` are instant while a dense 128-vertex random graph is unbounded — all three
identical under a vertex cap. The cap survives in RG-2 as a bound on *scratch state*, which
is what it actually guarantees, but it is no longer claimed to keep the head responsive.

### A per-operation wall-clock timeout
**Rejected because:** it makes the answer machine-dependent — the same input yields a proven
result on one machine and a refusal on another, which is indefensible for a CAS where
reproducibility is the product. This is also the objection RG-2's own plan raised, correctly,
before it knew `TimeConstrained` existed.

### Building a general abort facility for the graph subsystem
**Rejected because:** it already exists and is general. Building a second one would be
duplicate machinery with a narrower surface, and would not compose with nested
`TimeConstrained` the way the existing implementation deliberately does (it saves and
restores handler, timer, and `jmp_buf`).

### Returning the best-known result on abort instead of unevaluated
**Rejected because** — and this is the nearly-chosen option, since it looks strictly more
useful — a valid-but-unproven colouring is exactly the failure the exact search exists to
prevent. It contradicts the documented semantics silently, returning a plausible list of
integers that no caller can distinguish from a minimal one. A refusal is loud and correct;
this is quiet and wrong. The same reasoning caught a real bug during Phase 1: the
allocation-failure path was returning DSATUR's upper bound as though it were proven.

## Notes

- Ticket and measurements: `thoughts/shared/tickets/RG-2/plan.md` — `## Decisions`, the
  `## Overview` correction paragraph, and `## Performance Considerations`.
- Implementation: `src/graph/vertexcoloring.c` — `FVC_MAX_STEPS` (the reasoning above is
  restated at the constant) and the `tc_check_deadline()` poll in `fvc_bb`.
- Evidence rows: `tests/test_graph_slow.c` (`EXCLUDE_FROM_ALL` target `graph_slow_tests`,
  excluded from CI because the rows must spend the whole budget to prove anything).
- `TimeConstrained` implementation and its nesting/restore contract: `src/core.c:3774-4110`.
