# solverad — fragility write-up and strengthening

`src/solve/solverad.c` (`solverad_solve_radicals_equality`) is the leaf solver
for single-variable **radical** equations — anything carrying `Sqrt[...]` or
`x^(p/q)` that the polynomial and inverse-function specialists decline. Its
pipeline is:

```
lhs == rhs
  └─ e_orig = lhs − rhs
  └─ Together + Numerator            ← clear the (radical) denominator
  └─ substitution loop: u_i = base_i^(1/L_i)    (each radical → a fresh gen)
  └─ resultant chain: eliminate every u_i        → a polynomial in `var`
  └─ solvepoly                                     → candidate roots
  └─ per-candidate verification (N[] then Simplify) + dedupe + Solve::nongen
```

This document records the fragile / hung cases found while extending
`Solve`/`Reduce` to logarithmic & exponential equations (the peeled inner
equation of `1/2 Log[t] − Log[x − Sqrt[t]] == C[1]` lands here), the fix that
landed (2026-09-03), a larger fix that was **tried and reverted** with the
lesson it taught, and the limitations that remain.

---

## 1. The failure matrix

Probing `Solve[Sqrt[t]/(x − Sqrt[t]) == RHS, t]` across RHS shapes exposes three
distinct behaviours that all *look* like "solverad can't do it":

| RHS | status | mechanism |
|-----|--------|-----------|
| `q` (symbol), `a + b`, `a b` | SOLVED | — |
| `q^2` (integer power) | **HANG → FIXED** | `Simplify` blow-up in `verify_candidate` (§2) |
| `E^C[1]`, `2^q` (symbolic-exp power) | DECLINE (open) | `Together` won't clear (§3) |
| `Sqrt[5]` (rational-exp power) | DECLINE (open) | `Together` won't clear (§3) |
| same RHS but constant denom `1/(1 − Sqrt[t])` | SOLVED | — |

Two callers sit on top of solverad and shape what "correct" means here:

- `Solve[1/2 Log[t] − Log[x − Sqrt[t]] == C[1], t]` — after log-fusion the
  inverse-function solver peels a single `Log` and hands solverad
  `Sqrt[t]/(x − Sqrt[t]) == E^C[1]`, i.e. exactly the DECLINE row. **Still
  declines** (open, §3/§5).
- `DSolve\`Homogeneous[y'[x] == (2 x + y[x])/(x + 2 y[x]), y, x]` — the `y = v x`
  reduction integrates to and tries to invert
  `1/((v−1)^(3/4) (v+1)^(1/4)) == E^(C[1] + Log[x])`. This case is why the
  denominator-clearing fix was reverted (§4): **DSolve depends on solverad
  *declining* it** so it can return an implicit / `Root[]` first integral.

---

## 2. FIXED — HANG: `Simplify` in `verify_candidate` on a free-parameter residual

`verify_candidate` verifies each candidate root `N[]`-first, `Simplify`-fallback.
For `Sqrt[t]/(x − Sqrt[t]) == q^2` the candidate is `t = q^4 x^2/(1+q^2)^2`; the
residual carries the free parameters `q, x`, so `N[]` cannot numericise and the
code fell through to
`Simplify[ Sqrt[q^4 x^2/(1+q^2)^2]/(x − Sqrt[…]) − q^2 ]`. Simplify's radical
denesting on `Sqrt[q^4 x^2/(1+q^2)^2]` with the `(1+q^2)^2` denominator runs
**unboundedly** — the hang. (`q` differs only in that its `(1+q)^2` denominator
denests trivially, so it terminated by luck. `PossibleZeroQ` on the same residual
returns fast but says `False` — a *branch* artefact of random-complex sampling;
the candidate is valid on the principal positive-real branch.)

The deeper point: when a **free parameter survives**, rejection is *unsound*
anyway. A squared-radical candidate can satisfy the equation on one parameter
regime and fail on another — which is exactly why the parametric case
`Solve[Sqrt[a x + c] + 3 x == 5, x]` correctly keeps **both** branches with
`Solve::nongen`. So the only correct verdict there is `VERIFY_UNKNOWN` (keep +
nongen), and `Simplify` was never going to help.

**Fix (landed).** After `N[]` fails, `has_free_parameter` tests whether the
residual still contains a free parameter (any symbol other than the
`Indeterminate`/`ComplexInfinity`/`Infinity`/`Undefined`/`DirectedInfinity`
markers; `C[k]`'s head `C` counts). If so, return `VERIFY_UNKNOWN` immediately
**without** calling `Simplify`. `Simplify` is now reached only for a
*parameter-free* residual (a removable singularity), where it terminates. This
is both correct (matches the nongen convention) and hang-proof. Regression test:
`test_hang_parametric_rational_radical` in `tests/test_solveradicals.c`.

---

## 3. OPEN — DECLINE: `Together` won't clear a radical denominator against a `Power`

`Together[Sqrt[t]/(x − Sqrt[t]) − E^C[1]]` returns `−E^C[1] + Sqrt[t]/(−Sqrt[t]
+ x)` — the denominator is **not cleared**. (With a plain rational fraction,
`Together[t/(x − t) − E^C[1]]`, it clears fine; with an integer power `q^2` it
also clears. It balks specifically on a `Sqrt[var]` fraction summed with a
`Power` whose exponent is non-integer or symbolic — `E^C[1]`, `Sqrt[5]`, `2^q`.)
`Numerator` of the un-combined form still contains the variable in a
denominator, so `find_first_radical` + the substitution/resultant chain produce
nothing usable and solverad returns `NULL` (declines). This is why
`Sqrt[t]/(x − Sqrt[t]) == E^C[1]`, `== Sqrt[5]`, `== 2^q`, and the peeled case
`1/2 Log[t] − Log[x − Sqrt[t]] == C[1]` all decline.

These declines are **safe** (unevaluated, never wrong, never hang). Closing them
is future work — see §5.

---

## 4. TRIED AND REVERTED — abstract opaque constants before `Together`

The natural fix for §3 is to abstract every var-free non-numeric compound
subterm (`E^C[1]`, `Sqrt[5]`, `2^q`, `Log[x]`, …; **not** `Rational`/`Complex`,
which appear as radical exponents) to a fresh symbol, run `Together + Numerator`
(which then clears), and substitute the constants back. This was implemented and
**did** solve every DECLINE row and case #14 — but it was reverted, for a reason
worth recording:

- Clearing more denominators makes solverad *attempt* equations it used to
  decline. The DSolve homogeneous inversion
  `1/((v−1)^(3/4) (v+1)^(1/4)) == E^(C[1] + Log[x])` — **two distinct
  quarter-power radicals** — then clears into a degree-~16 resultant polynomial
  whose coefficients are the transcendental constant `E^(C[1]+Log[x])`.
  `solvepoly`'s factor/root search over that extension runs unboundedly: a new
  hang.
- A `SOLVERAD_MAX_DEGREE` cap stopped that hang, but did **not** restore prior
  behaviour: `DSolve\`Homogeneous` *relies on solverad declining* that equation
  so it can fall back to an implicit first integral. With the abstraction it
  instead produced a `Root[]` solution whose downstream `PossibleZeroQ`
  verification was itself pathologically slow, pushing `dsolve_tests` past two
  minutes at `t_homogeneous_algebraic`.

**Lesson.** In solverad, *declining is load-bearing.* Several higher-level
callers (DSolve homogeneous / reduction-of-order inversions) are built assuming
solverad gives up on transcendental-coefficient radical equations, and fall back
to implicit/`Root` forms. Making solverad "try harder" is not a local change —
it silently reroutes those callers. Any future denominator-clearing work must be
gated so it does **not** alter the decline set those callers depend on.

---

## 5. Remaining limitations / future work

1. **Power-valued-RHS declines (§3), including case #14.** To close them without
   the §4 regression, the abstraction must be paired with a guard that preserves
   solverad's decline on the equations DSolve relies on. Two viable routes:
   (a) keep the opaque constants **abstracted through the resultant *and*
   solvepoly** — solving a polynomial in `var` with symbolic (not transcendental)
   coefficients and restoring only at the very end, so no transcendental-extension
   factor search happens; combine with a degree cap so multi-radical blow-ups
   still decline. (b) only abstract-and-retry when the equation has a **single**
   radical generator with small `L` (the E^C[1]/case-#14 shape), leaving
   multi-radical inversions on the current decline path. Either way, re-run
   `dsolve_tests` (especially `t_homogeneous_algebraic`) as the gate.

2. **Branch fidelity of parametric solutions is advisory, not decided.** With a
   free parameter present, solverad keeps every squared-radical branch and
   raises `Solve::nongen`; it does not compute the regions on which each branch
   is valid. A principal-branch-aware sampler (random *positive-real* parameter
   probes, only flagging, never rejecting) could separate "always valid" from
   "regime-dependent", but positive-real sampling alone cannot soundly *reject*.

3. **No global work budget.** solverad bounds generator count
   (`SOLVERAD_MAX_GENS = 12`) and exponent-denominator LCM
   (`SOLVERAD_MAX_LCM = 120`) but not the final polynomial degree or a wall-clock
   budget on the resultant/solve/verify tail. A degree or effort cap would make
   solverad hang-proof against *any* future caller that feeds it a blown-up
   equation — the right home for the §4 `SOLVERAD_MAX_DEGREE` idea, decoupled
   from the abstraction.

4. **Verification cost.** `verify_candidate` runs `N[]` (and, for parameter-free
   residuals, `Simplify`) per candidate; no caching across candidates.
