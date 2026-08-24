# Fix: Integrate::nonelem emitted 4× instead of once

## Bug
`Integrate[Log[x] Log[Log[x]]/x^2, x]` emits the `Integrate::nonelem` message
four times:
- twice for internal recursive sub-integrands named `Integrate`DerivativeDivides`u$N`
- twice for the original integrand `(Log[x] Log[Log[x]])/x^2`

Desired: emit **once**, for the original integrand only.

## Root cause (verified by instrumentation)
Two independent defects:

1. **Leaky recursive messages.** `DerivativeDivides` substitutes `u = Log[x]`
   and recurses via `integrate_in -> Integrate[u Log[u] E^(-u), u]`. That
   nested Integrate reaches the same `RischTranscendental` decision stage and
   emits a `nonelem` message naming an internal gensym (`u$3`). A user-facing
   diagnostic must never name an internal recursion variable.

2. **Double cascade.** `builtin_integrate` is invoked twice at top level with
   the byte-identical integrand (same `expr_hash`). Cause: the integrand
   `Times[Log[x], Log[Log[x]], Power[x,-2]]` reorders under Orderless on the
   first `evaluate_step`, so `next != current` and the fixed-point loop runs a
   second step — re-running the whole (expensive) method search + Eliminate/Solve.
   Affects any failed symbolic integral whose integrand reorders on first eval.

## Plan
- [ ] `eval.c`/`eval.h`: add `eval_toplevel_id()` — a counter bumped once per
      outermost `evaluate()` call. General infra to scope per-command state.
- [ ] `integrate.h`/`integrate.c`: add `g_integrate_depth` counter around the
      method cascade (outermost user call = depth 1; recursion = depth >= 2).
- [ ] `integrate.c`: add a fail-memo keyed on `(toplevel_id, hash(f), hash(x),
      method)`. A matching re-entry in the same top-level evaluation returns
      NULL immediately — skipping the redundant second cascade AND its duplicate
      message. Self-invalidates next command via the id tag.
- [ ] `integrate_risch_transcendental.c`: gate the `nonelem` message on
      `g_integrate_depth <= 1` so only the top-level user integrand is named.
- [ ] Remove debug instrumentation.
- [ ] Verify: reported integrand -> one message; test C -> one message; a
      genuinely-solvable integral still solves; re-run in a new command re-warns.
- [ ] Run integration/calculus tests; docs + changelog.

## Review — DONE

Implemented both halves:
- `eval.c`/`eval.h`: `eval_toplevel_id()` — counter bumped once per outermost
  `evaluate()`. Zero semantic effect (only Integrate reads it).
- `integrate.h`/`integrate.c`: `g_integrate_depth` around the cascade; a 32-slot
  per-command fail-memo (`intg_fail_*`) keyed on `(eval_toplevel_id, hash(f),
  hash(x), method)`, auto-expiring each command, no dynamic allocation.
- `integrate_risch_transcendental.c`: `nonelem` message gated on
  `g_integrate_depth <= 1`.

Verified:
- Reported integral -> **one** message naming `(Log[x] Log[Log[x]])/x^2`; the
  cascade runs **once** (was twice — proven with `[CASCADE-RUN]`/`[MEMO-HIT]`
  instrumentation, since removed). ~0.078 s.
- Single-pass nonelem (`1/Log[Log[x]]`), `Log[Log[x]]`, `E^(x^2) Log[x]` each
  still warn once. Re-run in a later command re-warns (memo self-invalidates).
- Same nonelem twice in one `List` command -> one message; two *different*
  nonelem in one command -> one message each (fixed 32-slot table, not a single
  slot).
- Solvable / definite / list-threaded integrals unchanged.
- Suites green: integrals, integrate_dispatch, derivdivides, deriv,
  intrischnorman, limit, series, ramanujan, eval, eval_timestamps,
  eval_eager_exit, core, match. `make check-c99` clean.

## Follow-up correctness bug — FIXED
`Integrate[Sin[x]/Log[x], x]` returned a wrong **`0`** (also `Cos[x]/Log[x]`,
`Tan[x]/Log[x]`, `Sin[x]/Log[x]^2`, `ArcTan[x]/Log[x]`, `Gamma[x]/Log[x]`,
`Sin[x]/(1+E^x)`, `Cos[x]/(1+E^x)`, `BesselJ[0,x]/Log[x]`, …).

Root cause: the single-extension Risch cases in `risch_singleext.c`
(`rt_frac_try`, `rt_hermite_try`, `rt_hyperexp_case`) kernelize at `t=Log[x]`/`E^x`
and solve a Rothstein-Trager identity via `SolveAlways[..,{t,x}]`, but their gate
never verified the coefficients were rational in x. A `Sin[x]` coefficient passed
as a degree-0 poly in t; `SolveAlways[Sin[x]-k/x==0,{t,x}]` -> `{{Sin[x]->0,k->0}}`,
and the residues zeroed out to `0`. The module trusts the SolveAlways certificate
with NO diff-back, so an under-restricting gate = wrong answers.

Fix: `rt_is_ratl_in_xt(e,x,t)` — a `C(x)(t)` field-membership predicate — now
gates all three cases (`&& rt_is_ratl_in_xt(num,x,tsym) && rt_is_ratl_in_xt(den,...)`).
Out-of-field integrands decline instead of mis-certifying.

Verified: wrong-`0` set now empty across the whole class; `1/(x Log[x])`,
`1/(1+E^x)`, `1/(E^x(1+E^x)^2)`, special-fn recognizers, and elementary
`Sin[x] E^x` all unchanged; suites green incl. integrate_risch_transcendental,
risch_hermite, risch_field, risch_elementaryq, cherry_ei/li, knowles_erf,
intrat, intrat_corpus; `make check-c99` clean.
