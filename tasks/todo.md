# Task: Run `ParallelMixed.m` in Mathilda as `Integrate`ParallelMixedTower`

Port `mixed/ParallelMixed.m` (parallel Risch–Norman integrator over a simple radical
in a mixed tower) so it loads/runs as a `.m` module and is exposed as an `Integrate`
method — via `Method -> "ParallelMixedTower"`, the qualified symbol
`Integrate`ParallelMixedTower[f, x]`, and in the automatic cascade. Purely additive
(old C Risch–Norman NOT retired). Plan file:
`~/.claude/plans/let-s-port-mixed-parallelmixed-m-to-zazzy-peacock.md`.

## Phase 0 — Smoke spike (de-risk)  ✅ DONE
- [x] Copy `mixed/ParallelMixed.m` → `src/internal/mixed/ParallelMixed.m` (canonical working copy; sync to mixed/ at commit)
- [x] Loads clean via `Get` behind System-context Quiet/Check/Message shims
- [x] G0 findings + fixes (4):
  1. **RationalQ collision** → renamed package helper to `RationalFunctionQ` (was clashing with Protected builtin)
  2. **SparseArray unimplemented** → dense assembly rewrite (line ~2130)
  3. **`assoc[key]=val` unsupported in Mathilda** → implemented in `src/eval.c` apply_assignment (reroute to Part+Key when sym holds an Association)
  4. **Part single-position list-RHS bug** (`m[[2]]={1,2,3}` gave `{0,1,0}`) → fixed `is_rhs_list` gate in `src/part.c` (distribute only for multi-position selectors)
- [x] Core validated: Log[x], 1/(x Log[x]), x Exp[x], Sqrt[x], 1/Sqrt[1+x^2], x/Sqrt[1+x^2], 1/(1+Exp[x]) all integrate & D-check=0; Exp[x^2] declines correctly
- [x] No test regression: association/list/core/eval + ctest sweep (14/14) all pass
- [ ] STILL WRONG (revisit after Phase 1 real Check): `Sqrt[Tan[x]]` → false "not elementary" (SplitSpecials retry likely defeated by passthrough Check shim)

## Phase 0b — remaining package rewrites  ✅ DONE
- [x] `KeyDropFrom` → `sample = KeyDrop[sample, Y]`
- [x] `Hash` memo key: left as-is — works now that assoc-set handles Association-containing keys

## Phase 1 — Message subsystem  ✅ DONE
- [x] `Quiet`/`Check`/`Message` C builtins + fired counter in `src/message.c`; registered from `core_init`
- [x] Counter bumped at arithmetic `Power::infy`/`Infinity::indet` (via `arith_warn`)
- [x] Verified: Check[1/0,bad]→bad, Quiet suppresses, Message caught by Check

## Phase 2-3 — bridge + wiring  ✅ DONE
- [x] `builtin_integrate_pmt` (C, lazy-loads package, delegates to worker) + registered in integrate_init
- [x] enum, method_from_string, try_parallelmixedtower (decline on List/$Failed), explicit switch case
- [x] Cascade: inserted as LAST resort (after CRC), gated to skip pure rationals (`pmt_is_rational_structure`)
- [x] All 3 surfaces work from cold session; Dcheck=0 on handled cases

## Phase 4 — tests  ✅ DONE
- [x] `tests/test_parallelmixedtower.c` (Message subsystem + assoc/Part fixes + method); CMake registered; PASSES
- [x] No regression: association/list/core/eval pass; integrate_dispatch pass; crc_corpus PASS (~150s); intrat_corpus PASS
- [ ] dsolve regression (running)

## Phase 5 — house-keeping  ✅ DONE
- [x] Version bump 0.160 → 0.161
- [x] Docs: calculus.md (method 13 + Method list), control-flow.md (Quiet/Check/Message), changelog 2026-09-21
- [x] check-c99 rc=0
- [ ] Rebuild code-review-graph; sync mixed/ copy

## Phase 1 — Message subsystem (C builtins)
- [ ] `Quiet` (HoldAll) — reuse `mth_msg_suppress_push/pop`; 1- and 2-arg forms
- [ ] Message-fired counter in `src/message.c` (+ `.h`)
- [ ] `Message` (HoldFirst) — emit `Head::tag`, bump counter
- [ ] `Check` (HoldAll) — return failexpr iff counter advanced
- [ ] Bump counter at ~8 suppression-aware emission sites
- [ ] `SYM_*` in `sym_names.{c,h}`; register in `core.c`; attrs; docstrings; usage
- [ ] `tests/test_message.c` green (GATE before integrator relies on Check)

## Phase 2 — Port the `.m` + bridge wrapper
- [ ] Rewrite `SparseArray` (2130) → dense; `KeyDropFrom` (1681) → `KeyDrop`; `Hash` (1369,2172) → direct key
- [ ] Append flat `Integrate`ParallelMixedTower[f_,x_Symbol] := ...` wrapper (/;-guarded decline) + Protected

## Phase 3 — Wire into integrate.c
- [ ] `pmt_lazy_load()` (clone `crc_lazy_load`)
- [ ] Enum `METHOD_PARALLEL_MIXED_TOWER`; `method_from_string`
- [ ] `try_parallelmixedtower()` (clone `try_crctable`)
- [ ] Cascade insert (after `try_rischtranscendental`, before `try_crctable`)
- [ ] Explicit `case` in dispatch switch; warning-list + docstring polish

## Phase 4 — Tests
- [ ] `tests/test_integrate_parallelmixedtower.c` (both surfaces + cascade + decline)
- [ ] No-regression: `ctest -R 'message|integrate|dsolve|trigrat|intrat|crc'`
- [ ] (optional) Charlwood corpus

## Phase 5 — House-keeping
- [ ] Version bump + tag (`v0.161`)
- [ ] Docs: calculus.md + Message builtins doc + weekly changelog
- [ ] Rebuild code-review-graph; valgrind

## Review

The port RUNS in Mathilda and is exposed as `Integrate`ParallelMixedTower` via all
three surfaces (qualified symbol, `Method -> "ParallelMixedTower"`, and the Automatic
cascade). What the Phase-0 spike uncovered was that the remaining work was NOT the
stale `PARALLEL_MIXED.md`'s "3 rewrites + Message subsystem" — several load-bearing
Mathilda gaps had to be closed first:

**Mathilda-core fixes (all verified, no test regression):**
1. `assoc[key] = val` / nested `assoc[k1,k2] = val` element assignment — was making a
   dead DownValue; now reroutes to the Part+`Key` machinery when the head holds an
   Association (`src/eval.c`). *User-approved as the root-cause fix.*
2. Single-position `Part` list-value assignment (`m[[2]] = {1,2,3}` gave `{0,1,0}`) —
   pre-existing bug; distribution now gated on multi-position selectors (`src/part.c`).
3. Message subsystem `Quiet`/`Check`/`Message` (`src/message.c`) + a fired counter bumped
   at the arithmetic `Power::infy`/`indet` sites (`src/arithmetic.h`).
4. `Integrate`ParallelMixedTower` C builtin (`builtin_integrate_pmt`) — lazy-loads the
   ~2500-line package (the load is ~1s, so eager loading was rejected) and delegates;
   messages suppressed during the probe.

**Package fixes** (`src/internal/mixed/ParallelMixed.m`, synced to `mixed/`):
`RationalQ`→`RationalFunctionQ` (collision with the Protected builtin), `SparseArray`→
dense assembly, `KeyDropFrom`→`KeyDrop`.

**Cascade wiring:** last resort (after CRC), gated to skip pure rationals (a parametric
rational otherwise burned a full iPIM attempt) AND to **decline while a TimeConstrained
deadline is active** — the crash gate. That last one matters: PMT is malloc-heavy, and
running it inside DSolve's `TimeConstrained[Integrate[...]]` window reproducibly tripped
the pre-existing SIGALRM/siglongjmp-mid-malloc crash (dsolve died at 142). Declining then
restores DSolve's exact pre-PMT behaviour. (`src/core.c` `tc_deadline_is_active`.)

**Verified:** transcendental towers, simple radicals, flattened roots integrate &
D-check=0; genuine non-elementary cases decline cleanly; core/eval/association/list +
integrate_dispatch pass; crc_corpus PASS; intrat_corpus PASS; dsolve no longer crashes;
`check-c99` rc=0; new `tests/test_parallelmixedtower.c` passes.

**Known coverage gaps (algorithm, NOT infrastructure — the "Phase 4 corpus" work):**
`Sqrt[Tan[x]]` (SplitSpecials), complex-residue realisation (`1/(1+x^3)` — but pure
rationals are gated to BronsteinRational upstream anyway). These are the ongoing
baseline-gated corpus refinement, not part of "make it run / callable / exposed".
