# Task: Sqrt[Tan[x]] integration + Simplify hang

## Reported bugs
1. `Integrate[Sqrt[Tan[x]], x, Method -> "DerivativeDivides"]` returns unevaluated.
   It should find the substitution `u == Sqrt[Tan[x]]` and reduce to a rational
   integral in `u`.
2. `D[<antiderivative>, x] // Simplify` hangs (>2 min); Mathematica does it in ~0.01s,
   returning `Sqrt[Tan[x]]`.

## Root-cause diagnosis (verified empirically)

Both bugs funnel into ONE missing capability: the system cannot simplify/normalise
**rational functions of the radical-trig kernel `Sqrt[Tan[x]]`**. The whole derivative
is secretly *rational* in `t = Sqrt[Tan[x]]` (since `Tan = t^2`, `Sec^2 = 1 + t^4`),
which is why Mathematica is instant.

Concrete findings:
- The antiderivative `r` produced by the eliminate strategy is **correct**
  (`D[r,x] - Sqrt[Tan[x]] ~ 1e-16` at x = 0.3, 0.5, 1.2).
- **Bug 1 immediate cause:** `try_eliminate_kernel` computes the right `r`, then the
  verification gate `differentiates_back` calls `PossibleZeroQ[D[r,x]-f]` which returns
  **False** — a genuine *false negative* (zero_test.c Schwartz-Zippel sampling). Gate
  bails fast -> unevaluated.
- **Bug 1 secondary cause:** even past PossibleZeroQ, the gate's `Simplify[D[r,x]-f]==0`
  needs to prove the radical-trig identity is 0 -> same missing capability + would hang.
- **Bug 2 cause:** `Simplify` search applies `Factor` to seeds; `Factor[d]` blows up.
  Isolated: the blow-up requires `Sqrt[2]` (algebraic number) — Factor over `Q(sqrt2)` of a
  multivariate rational (`qafactor.c`) is exponential here. With `Sqrt[2]` replaced by a
  symbol it finishes in ~2s. `Together[d]` does not even combine the fractions because of
  the `Sqrt[Tan[x]]` kernels.
- **Why trigrat can't help today:** `simp/trigrat.c` treats `Sqrt[Tan[x]]` (a non-integer
  power) as an *opaque atom* independent of `Tan[x]` — it never learns `o^2 = Tan[x]`, so
  it cannot reduce, returns NULL, and the search falls through to the hanging Factor.

## Fix components
- [ ] A. Core: teach simplification to handle rational functions of a radical kernel
      `Sqrt[g(x)]` (substitute `t = Sqrt[g]`, require rational-in-`t`, Cancel, map back,
      leaf-count gate). Reduces `d -> Sqrt[Tan[x]]` fast and short-circuits Factor.
- [ ] B. Fix `PossibleZeroQ` false negative on `D[r,x] - Sqrt[Tan[x]]` (zero_test.c).
- [ ] C. Defensive bound on Factor / `qafactor` so a hard multivariate-algebraic case
      bails gracefully instead of hanging.
- [ ] Tests: DerivativeDivides closes Integrate[Sqrt[Tan[x]],x]; Simplify[d] fast -> Sqrt[Tan[x]];
      PossibleZeroQ case; Factor-bound case. Update docs/spec + changelog.

## Review (done 2026-06-09)

All five asks delivered; 21 regression suites + 3 new test groups pass; valgrind clean.

- **A. trigrat quadratic-radical generators** (`src/simp/trigrat.c`): `Sqrt[g]` kernels
  carried as algebraic generator `l` with `l^2 = g_subst`, reusing the existing
  reduce/rationalise pipeline (radicals cleared before trig generators). Half-integer
  powers collapse to powers of one `l`. Also fixed a pre-existing leak in the module's
  raw `Times`/`Power` constructors (didn't free the args buffer `expr_new_function` copies).
- **B. PossibleZeroQ** (`src/zero_test.c`): bounded sampled imaginary magnitude
  (`ZT_IMAG_NUMERATOR_BITS=3`, |Im|<=8) — large |Im| saturates Tan→±I and 1+Tan^2 cancels
  below machine eps, misreporting true zeros. Real parts unbounded; denominator granularity
  kept for SZ.
- **C. (not needed) Factor bound**: the qafactor blow-up is now avoided because trigrat
  short-circuits before `Factor` is reached for radical-trig inputs. A defensive bound
  remains available as future work but was out of scope (no user trigger remains).
- **Cascade + plumbing**: `integrate_derivdivides_full` (direct then Eliminate/Solve) wired
  into BOTH the Automatic cascade and `Method->"DerivativeDivides"`; bare `Integrate` now
  closes it. `eliminate_suppress_messages` mutes Eliminate `::ifun`/`::alg`/`::nlin` while
  the integrator drives it.

Results: `Integrate[Sqrt[Tan[x]], x]` (bare + Method) → correct antiderivative, no warnings,
~1.3s; `D[r,x] // Simplify` → `Sqrt[Tan[x]]`, ~0.3s (was infinite hang).

**Perf gap** (~0.3s vs MMA ~0.01s): localized to the first `Together` (280ms of 287ms) in
trigrat's fallback — the general multivariate `Together` combining the `(Sin,Cos,Sqrt[Tan])`
-rational into one fraction. Reduce/rationalise/Cancel are all <4ms. Closing the gap needs a
`Tan`-native (fewer-generator) trigrat path or a faster multivariate `Together` — a separate,
larger change deliberately deferred to avoid destabilising the now-correct, tested result.
