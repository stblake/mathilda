# FindInstance stress round 2 — 6 harder cases

Plan: `/Users/user/.claude/plans/cheeky-jingling-dream.md`. Sound, verify-gated,
algorithmic methods + quiet internal probes.

## Cases
- **5** `Log[z^2]==2Log[z]+2πI` ℂ — WRONG `{{z->0}}` (Log[0] undefined). Soundness.
- **4** exp eqn + sign bounds ℝ — declined; solve-one-var.
- **8** rational-eqn tiny region ℝ — declined + spam; solve-one-var.
- **3** two rational eqns inexact ℝ — declined + Power::infy ×47; feasibility objective.
- **9** Rastrigin<0.1 && x!=0 && y!=0 ℝ — NMinimize::nimpl; feasibility objective.
- **7** two poly eqns + x!=0,y!=0 ℂ — Solve::nsdim; ideal saturation.

## Tasks
- [x] M0. `fi_verify` definedness gate (`fi_has_nonfinite`, `fi_defined_truth`). Fix #5.
- [x] M4. `src/message.{c,h}` counter; routed into `arith_warnings_muted`, `warn_nsdim`, `fm_warn`, `emit_ifun`; wrap FindInstance search. Kills spam.
- [x] M1. `fi_solve_one_sample` (solve one eqn for one var, per-var-compacted grid). #4, #8.
- [x] M2. `fi_collect_penalty` (all-constraint penalty) + `fi_fold_aux` (pin aux consts) + broader seeds. #3, #9.
- [x] M3. `fi_saturate_solve` (Rabinowitsch slack → zero-dim Solve). #7.
- [x] Reorder: exact grid sampler BEFORE numeric feasibility (restores exact In6).
- [x] Verify: all 9 correct + no spam; In[1/2/6] unchanged; each witness verifies True; 1.94s.
- [x] Tests: reduce/comparisons/reduce_corpus(0/158)/solve_corpus(0/99)/nroots/findroot/ratcanon_reduce/rootreduce green; one stale `reduce_tests` assertion updated. check-c99 clean; valgrind no Mathilda frames in leak stacks.
- [x] Docs: solutions-of-equations.md, changelog 2026-08-24.md (v0.108), version 0.107→0.108.

## Review (v0.108, 2026-08-27)

All 6 new stress cases solved; In[1]/In[2]/In[6] unchanged; total 1.94s.
- **#5 (soundness)** — `fi_verify` now requires each relation operand to be DEFINED
  (no infinity/indeterminate sentinel), so the wrong `{{z->0}}` (where `Log[0]=-∞`
  folded `-∞==-∞` to True) becomes the correct `{{z->-I}}`. Acceptance-only ⇒ sound.
- **#4, #8** — `fi_solve_one_sample`; exact witnesses (`c1->E`; `y->10⁻⁶`).
- **#7** — `fi_saturate_solve` (Rabinowitsch slack → `Solve` returns the 10 roots).
- **#3, #9** — `fi_collect_penalty` all-constraint objective + `fi_fold_aux`. Numeric
  feasibility is now the true last resort (after the exact sampler), so In6 keeps
  its exact `{t->1, s->-1}`.
- **Spam** — `src/message.{c,h}` depth routed into the existing arith mute + the
  nsdim/fm_warn/ifun emitters; FindInstance wraps its search in it.

**Lesson:** a verify-gated search is only as sound as its verifier. `expr /. pt ===
True` is insufficient when the evaluator folds arithmetic on infinities; definedness
must be checked per relation operand.
