---
ticket: RG-3
created: 2026-08-31T15:03:08Z
researcher: Michael Sollami
topic: "Which builtins must poll tc_check_deadline and do not"
type: research
lifecycle: active
full_research: thoughts/shared/tickets/RG-3/research.md
---

# Research Summary: Which builtins must poll `tc_check_deadline` and do not

**Full research (appendix)**: `thoughts/shared/tickets/RG-3/research.md`

## Recommendation

**Split — decided by the human on 2026-08-31, and done.** The state-corruption bug is now
`thoughts/shared/tickets/RG-4/research.md`; RG-3 keeps the polling audit.

RG-3's answer: only **2 of 832** registered heads poll `tc_check_deadline` (`GroebnerBasis`
and `FindVertexColoring`, across 5 call sites), against **at least 41** ad-hoc local cap
constants — roughly twenty to one. The gap is real but **host-dependent**: on any host with
working `ITIMER_PROF` the SIGPROF path catches non-polling heads anyway, so this changes
nothing observable on Linux/macOS and matters on WSL 1. It is also **not detectable by
testing**, which is the finding that shapes the work: the audit must be static, in the
`make check-packed-aware` family.

One sequencing constraint, narrowed after review: **do not add polls inside `src/calculus/`
until RG-4 lands.** That subsystem is the corruption surface, and a new poll there is a new
reliable way to trigger it. None of the top four polling candidates lives there, so the rest
of the audit can proceed in parallel.

## Options Considered

For how the audit itself should work:

1. **A static `make check-abortable` audit (recommended)** — enumerate heads with unbounded
   C-level loops and diff against the poll sites, in the family of `check-packed-aware`. The
   only option that works, since the gap is invisible to runtime testing on a healthy host.
   Costs: like the other audits, it needs a curated exemption list to stay honest.
2. **Test-based detection** — ruled out, not merely dispreferred: on Linux/macOS the signal
   catches polling and non-polling heads identically, so every test passes either way.
3. **Add polls opportunistically as slow heads are reported** — cheapest, and how the four
   existing sites came to exist. Leaves the gap unmeasured, which is the state that produced
   a twenty-to-one ratio nobody had noticed.

## Decision Criteria

- **The gap is invisible to testing**, so a static audit is the only mechanism that can find
  it. This is the single fact that determines the shape of the work.
- **Audit silence must not mean "never looked"** — SPEC.md §9's argument for the packed-array
  audits transfers verbatim, including the need for an explicit exemption list. Some heads
  genuinely cannot cooperate: `pcre2_match` (`src/strings/regex/regex_engine.c:98`) runs
  inside PCRE2's stack with no reachable callback, so its remedy is PCRE2's own
  `match_limit`/`depth_limit` — currently unset — not a poll.
- **A wall-clock local budget is an anti-pattern to fix, not to copy.**
  `PMINT_BUDGET_SEC 4.0` (`src/calculus/intrischnorman.c:79`) makes the answer
  machine-dependent *and* is consulted once, before `RowReduce` (`:3087`), so it cannot stop a
  run already under way.
- **Do not poll inside `src/calculus/` until RG-4 lands** (see Recommendation).

## Open Questions

### Unresolved

- [ ] Does excluding `src/calculus/` from the first polling pass suffice, or should the whole
      audit wait on RG-4?
- [ ] Do the top-ranked targets (`src/solve/`, `src/strings/regex/`, `src/poly/mvfactor.c`)
      hold stranded statics of their own? Worth checking before a poll lands in any of them.

### Resolved

- [x] One ticket or two? — **Two**, decided by the human 2026-08-31. The corruption bug is
      `thoughts/shared/tickets/RG-4/`.
- [x] Does `TimeConstrained` interrupt pure-C builtins? — Yes, via SIGPROF, verified on four
      heads. The plan-level premise that it does not was already corrected in RG-2.
- [x] Can the polling gap be found by testing? — No. On any host with working `ITIMER_PROF`
      the signal catches polling and non-polling heads alike. The audit must read source.
- [x] Do the existing pollers recover cleanly from abort? — Yes; `GroebnerBasis` leaks heap but
      holds no static state and answers correctly afterwards.

## Requires Approval

One left. **The polling audit's value is host-dependent**: on Linux/macOS it changes nothing
observable, so the case for spending on it rests on WSL 1 users and on not wanting a
twenty-to-one ratio of private caps to shared-channel cooperation — not on any defect
measurable on this machine. Worth an explicit yes before the work is scheduled.

The other scope call is settled: the live bug this research found is split out as RG-4.
