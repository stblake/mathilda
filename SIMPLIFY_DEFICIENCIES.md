# Simplify deficiencies (noted 2026-07-14)

Surfaced while making the transcendental Risch integrator return real closed forms
for rational trigonometric integrands (roadmap P2). Each item is a concrete case
where `Simplify` (and friends) cannot prove a true identity or reduce a form to a
clean equivalent. Correctness note: several "non-reductions" are actually *correct*
(the two sides differ off a principal branch); those are marked NOT-A-BUG.

Directive (user): these should be **noted and fixed** — do not paper over them with
numeric verification or Weierstrass substitutions. This file is the "noted" half;
each fix is a scoped Simplify/TrigReduce improvement.

## What already works (baseline — do not regress)
- `Simplify[D[Tan[x], x] - Sec[x]^2]` -> `0`.
- `Simplify[D[(Sec[x] Tan[x] + Log[Sec[x] + Tan[x]])/2, x] - Sec[x]^3]` -> `0`
  (the CLEAN Sec^3 antiderivative diff-backs fine; only the ugly multiple-angle
  form below does not).
- `Simplify[Log[2 - 2 Cos[x]] - Log[2 + 2 Cos[x]] - 2 Log[Tan[x/2]]]` -> `0`.

## D1 — exp-log real collapse  (REAL deficiency)
`Simplify[I x - Log[1 + E^(2 I x)] + Log[Cos[x]]]` returns unchanged; the expression
is `0` (the standard `I x - Log[1 + E^(2 I x)] == -Log[Cos[x]]`, up to a constant).
No path (`Simplify`/`FullSimplify`/`ComplexExpand`/`TrigReduce`) collapses an
`a x + Sum c_k Log[1 + b_k E^(m_k I x)]` combination to its real trig closed form.
- Impact: the transcendental Risch trig front-end's raw output is I-laden here; the
  new `rt_realify` (real/imag decomposition) sidesteps it for integrator output, but
  the underlying Simplify identity is still missing for user-facing expressions.
- Fix sketch: a `Simplify`/`ComplexExpand` rule that rewrites `Log[1 + r E^(I u)]`
  (|r| suitably bounded) via `Log[...]= (1/2)Log[1 + 2 r Cos u + r^2] + I ArcTan[...]`
  and folds the resulting `I`-linear terms — i.e. implement `ComplexExpand` for
  `Log`/`ArcTan` over real arguments (currently `ComplexExpand` is an unimplemented
  bare symbol; see D4).

## D2 — trig-power form reduction  (form-quality deficiency)
`Integrate\`RischTranscendental[Sec[x]^3, x]` closes (correct, real) but as a
multiple-angle mess:
`(8 Sin x + 4 Sin 5x + 12 Sin 3x + (...)Log[2-2 Sin x] + (...)Log[2+2 Sin x])/(40 + ...)`
rather than the clean `(Sec x Tan x + Log[Sec x + Tan x])/2`. `Simplify` does not
reduce the multiple-angle form to the compact one, and (crucially) `Simplify` cannot
prove the diff-back `D[uglyform] - Sec[x]^3 == 0` — which is why the integrator
currently accepts it via a numeric diff-back fallback (provisional; a Simplify fix
retires it).
- Fix sketch: strengthen `TrigReduce`/`Simplify` to (a) reduce `Sin[k x]`/`Cos[k x]`
  polynomial-over-polynomial forms to `Sec`/`Tan` powers, and (b) prove multiple-angle
  identities zero (expand to a common `Sin[x]`/`Cos[x]` basis and cancel).
- `TrigReduce[Sec[x]^3]` returns `Sec[x]^3` unchanged (does not expand secant powers).

## D3 — ComplexExpand is unimplemented  (missing capability, blocks D1)
`ComplexExpand[3 + 4 I]`, `ComplexExpand[Log[1 + I x]]`, `ComplexExpand[Abs[Cos[x] +
I Sin[x]]]` all return unevaluated — `ComplexExpand` is a registered symbol with no
implementation. Implementing it (real/imag decomposition assuming symbols real, for
`Plus/Times/Power/Log/ArcTan/ArcTanh/Exp/trig`) is the general capability that also
fixes D1 and lets user-facing I-laden trig expressions be realified (the integrator
now carries a private `cx_reim` doing exactly this for its own output; promoting it
to a real `ComplexExpand` builtin is the general fix).

## NOT-A-BUG (Simplify is correct to leave these)
- `Simplify[ArcTan[Cos[x], Sin[x]] - x]` stays unreduced: `ArcTan[Cos x, Sin x]`
  (= `Arg[E^(I x)]`) equals `x` only on `(-Pi, Pi)`; they differ by `2 Pi k`
  elsewhere, so a global rewrite would be WRONG. The integrator's atan2 output is
  the globally-correct form; do not add a rule collapsing it to `x`.
