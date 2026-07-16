# DiffUnderInt — Blockers for the remaining examples

> **Update 2026-07-08 round.** Three blockers are now (partly) resolved inside
> `integrate_diffunderint.c`, no engine changes required:
> - **B1 Gaussian (partial):** a `gaussian_halfline` moment family
>   (`∫₀^∞ xⁿ e^{-p x²}{1,cos} dx`) plus a direct `Erf` parameter back-integration
>   now close the differentiate-to-cosine-moment cases, e.g.
>   `∫₀^∞ e^{-x²} Sin[a x]/x dx = (π/2)Erf[a/2]`. (The self-similar `e^{-x²−a²/x²}`
>   substitution and the Sin-Gaussian/Dawson moment remain open.)
> - **Decaying sinc / non-even half-line (new):** `rational_halfline_general`
>   integrates a non-even real rational over `[0,∞)` to real `ArcTan`/`Log`,
>   closing e.g. `∫₀^∞ e^{-c x}(1−Cos[a x])/x² dx`. This is the capability B4
>   named for the decay-sinc case, obtained *without* complex-Log reduction.
> - **B4 output (partial):** a `PowerExpand`-based `diui_finalize` now crushes the
>   `Log[1/Sqrt[u]]+Log[Sqrt[v]]` outputs to clean `Log` differences (#5/#7/#21
>   style). The genuine complex-Log→ArcTan reduction and the piecewise (B5),
>   trig-period (B2b), and log-rational-base (B6, #2/#16) items remain open.
> - The spurious `Power::infy` message flood on declined piecewise cases (#5/#19)
>   is fixed (search wrapped in `arith_warnings_mute`).


`Integrate\`DiffUnderInt` (differentiation under the integral sign / Feynman's
trick) is implemented, wired into the definite cascade, tested, documented, and
leak-clean. It closes **12 of the 24 target examples** correctly and fast
(≤ 6 s, zero hangs): #1, 2, 3, 5, 6, 7, 13, 14, 17, 20, 21, 22.

The method itself is sound — the differentiate → inner-integrate → integrate-back
→ fix-constant loop works and is verified symbolically (`Simplify[D[I,p]−J]===0`
plus an exact base). The remaining 12 examples are blocked **not by the Feynman
logic but by capabilities missing from the surrounding engine**, which the
inner-integral families cannot yet substitute for. Each blocker below names its
root cause, the examples it blocks, and a concrete fix.

Source: `src/calculus/integrate_diffunderint.c`.

---

## B1 — The engine cannot integrate Gaussians (blocks #4, #10, #12, #24)

**Symptoms.**
- `Integrate[Exp[-x^2], {x, 0, Infinity}]` → returns unevaluated (should be `Sqrt[Pi]/2`).
- `Integrate[Exp[-x^2] Cos[a x], x]` (indefinite) → **hangs ~14–25 s**.
- `Limit[Erf[x], x -> Infinity]` → unevaluated (should be `1`).

**Consequence.** DiffUnderInt pre-screens Gaussian integrands
(`contains_gaussian_exp`: an `Exp` of a nonlinear/singular power of `x`) and
declines them up front, so #4/#10/#12/#24 return unevaluated fast rather than
hanging.

**Fix.** A Gaussian moment-integral evaluator for the half-line forms
`∫₀^∞ xⁿ e^{-p x²} {1, cos(q x), sin(q x)} dx` → closed forms in
`Sqrt[Pi]`, `Erf`, `e^{-q²/4p}` (analogous to the existing Laplace/rational
families), plus teaching `Limit` that `Erf[∞]=1`, `Erfc[∞]=0`. This directly
unlocks #4, #12, #24. #10 (`Exp[-x² − a²/x²]`) additionally needs the
self-similar substitution `x → a/x` (Glasser's master theorem) to reduce to a
Gaussian — a small dedicated recognizer.

**Effort.** Medium. One new family file mirroring `laplace_halfline`, plus a
couple of `Erf`/`Erfc` limit rules. Highest value: 4 examples.

---

## B2 — The engine hangs on elementary Laplace/Fourier and trig-period definite integrals

**Symptoms.**
- `Integrate[Exp[-a x] Cos[b x], {x, 0, Infinity}, Assumptions -> a > 0]`
  → **hangs ~30 s** (the indefinite `Integrate[Exp[-a x] Cos[b x], x]` alone
  takes ~14 s, then the symbolic `Limit` at ∞ adds more).
- `Integrate[<rational-in-trig>, {x, 0, Pi}]` → **hangs** (the inner integrals of
  #8, #9, #15, #19).

**Root cause.** The indefinite cascade (Risch–Norman / table search) is
exponentially slow on `exp × trig`, and the definite path compounds it with slow
symbolic limits.

**Consequence.** DiffUnderInt never hands half-line-exp or finite-period-trig
forms to the general engine (see B3); it uses its own families instead. Anything
outside the implemented families is declined. This is why trig-period cases
(#9, #15) are not yet covered.

**Fix (two independent options).**
- (a) Speed up the general engine's `exp × trig` indefinite integration — a
  system-wide win, larger and riskier.
- (b) Build a **trig-period definite evaluator** in DiffUnderInt: Weierstrass
  `t = tan(x/2)` (or `t = tan(x)`) substitution → rational integrand over a
  finite interval → residue / partial fractions, for `{0, π}` and `{0, π/2}`.
  This is the natural family for #9 and #15 (and the finite base of #18).

**Effort.** Option (b): Medium–High. Unlocks 2 (plus enables #18's base).

---

## B3 — `TimeConstrained` does not bound a nested `evaluate` (infrastructure)

**Symptom.** `TimeConstrained[Integrate[Sin[b x]/x, {x,0,Infinity}], 6]` returns
`$Aborted` after 6 s at the REPL, but the **same call issued from inside a
running builtin** ran for **72 s** — the nested alarm / `longjmp` does not
re-arm (the documented async-safety limitation of the `TimeConstrained` handler).

**Consequence.** DiffUnderInt cannot safely use the general engine as a
time-bounded fallback, which forced the **families-only** architecture: it must
be able to compute each inner integral itself, or decline. This is a robustness
*win* (the method never hangs), but it means coverage is bounded by the families
rather than by the (slow, incomplete) general engine.

**Fix.** Make `TimeConstrained`'s alarm/`longjmp` re-arm correctly under a nested
`evaluate` (repl/eval infrastructure), **or** add a cooperative deadline flag
that the hot loops (Risch, Gröbner, factoring) poll. Either would let DiffUnderInt
fall back to the engine safely for forms no family covers.

**Effort.** Medium (infrastructure), but broadly useful beyond this feature.

---

## B4 — Simplify lacks complex-Log → real (ArcTan/Log) reduction (blocks #23, #16; ugly forms for #5, #21)

**Symptoms.**
- `Simplify[Log[-I a] - Log[I a], a > 0]` → does **not** reduce to `-I Pi`.
- `ComplexExpand[Re[1/(a - I b)]]`, `ComplexExpand[Log[c - I a]]` → unevaluated.
- `Simplify[ArcTan[b/a] + ArcTan[a/b], a>0 && b>0]` → not reduced to `Pi/2`.
- `Simplify[Log[a] - Log[b] - Log[a/b]]` → not reduced to `0`.

**Consequence.**
- The k=1 sinc-**with-decay** result `∫₀^∞ e^{-c x} Sin[a x]/x dx = -Σ cⱼ Log(-αⱼ)`
  (e.g. #23 → `ArcTan[a/c]`) stays as an unsimplifiable complex-Log form, so only
  the pure-**oscillatory** sinc (whose `M(s)` is even and closes via
  `rational_halfline`) is covered. #23 is therefore declined.
- #5 and #21 come out **correct** but in `Sqrt[b²]` / split-`Log` forms that
  Simplify cannot crush to the clean `ArcTan[b/a]` / combined `Log`.

**Fix.** Implement `ComplexExpand[Log[…]]` / `Re` / `Im` reduction for complex
arguments, plus the identities `ArcTan[x] + ArcTan[1/x] = π/2 · Sign[x]` and
`Log[a] − Log[b] − Log[a/b] = 0`. Then the direct Frullani formula
`-Σ cⱼ Log(-αⱼ)` yields clean real ArcTan/Log output and #23 closes.

**Effort.** Medium (Simplify/ComplexExpand work), reusable across the system.

---

## B5 — No piecewise / case-split output (blocks #8, #19)

**Root cause.** These integrals genuinely have **piecewise** results — #8 is `0`
for `|a|<1` and `2π Log|a|` for `|a|>1`; #19 is `(π/2) Min[a,b]`. The
differentiated inner integral `J` is *conditional* (`Sign`/`Boole`/`Piecewise`),
which the current families and the parameter-integration treat as opaque and
decline.

**Fix (Stage C, self-contained in DiffUnderInt).** Detect a conditional `J`;
find the branch-boundary locus (roots of the condition arguments); integrate each
branch over the parameter separately, matching the additive constant across
boundaries by continuity of `I`; emit `Piecewise`, collapsing to `Min`/`Max`
when recognizable.

**Effort.** High (branch analysis + continuity matching).

---

## B6 — No log-rational or finite-radical family (blocks #16, #18)

- **#16** `Log[x²+a²]/(x²+b²)` on `{0,∞}`: its base is
  `∫₀^∞ Log[x²]/(x²+b²) dx = (π/b) Log[b]` — a **log × rational** half-line
  integral no family covers. Needs a log-rational family (differentiate a
  parameter into the log, or a Mellin evaluator). (Also depends on B4 to clean
  the result.)
- **#18** `ArcTan[a x]/(x Sqrt[1−x²])` on `{0,1}`: its inner integral
  `∫₀^1 dx/((1+a²x²)Sqrt(1−x²))` is a **finite-interval radical** form (engine
  slow ~17 s; no family). Fix via the trig substitution `x = sin θ` →
  rational-in-trig over `{0, π/2}` — i.e. it rides on the trig-period evaluator
  from B2(b).

**Effort.** Medium each; #18 largely falls out of B2(b).

---

## B7 — `#11` needs a non-Feynman path (minor)

`Integrate[1/(x²+a²)³, {x,0,∞}] = 3π/(16 a⁵)` is a plain **even-rational
half-line** integral — `rational_halfline` already computes exactly this form.
But DiffUnderInt *differentiates first*, which raises the power and makes it
harder. The clean fix is for the definite cascade to try the direct rational
half-line (residue) *before* DiffUnderInt for such integrands (they don't need
Feynman at all), or for DiffUnderInt to short-circuit when the integrand is
already directly integrable by a family.

**Effort.** Low.

---

## Blocker → example map

| Example | Blocked by |
|--------|------------|
| #4, #12, #24 | B1 (Gaussian) |
| #10 | B1 + self-similar substitution |
| #9, #15 | B2(b) (trig-period) |
| #18 | B2(b) + finite-radical (B6) |
| #8, #19 | B5 (piecewise) |
| #23 | B4 (complex-Log reduction) |
| #16 | B6 (log-rational) + B4 |
| #11 | B7 (direct rational before Feynman) |

## Recommended order (value / effort)

1. **B1 Gaussian family** — unlocks 3–4 examples, self-contained, medium effort.
2. **B4 complex-Log reduction** — unlocks #23, cleans #5/#21/#16, reusable.
3. **B2(b) trig-period evaluator** — unlocks #9/#15 and #18's base.
4. **B5 piecewise (Stage C)** — #8/#19, higher effort.
5. **B7** — cheap cascade-ordering tweak for #11.
6. **B3 TimeConstrained nesting** — infrastructure; makes an engine fallback safe
   and is broadly useful, but not required for the family-based fixes above.
