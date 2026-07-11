# `Integrate`RischTranscendental` — Implementation Status

The **recursive** transcendental Risch integration algorithm (the
Bronstein/Roach decision-procedure lineage), grounded in Mathilda's
Expr/poly/rat machinery.

- Source: `src/calculus/integrate_risch_transcendental.c` (+ `.h`)
- Tests: `tests/test_integrate_risch_transcendental.c`
- Dispatch: wired into `src/calculus/integrate.c`
- Reference docs: `docs/spec/builtins/calculus.md`,
  `docs/spec/changelog/2026-07-06.md`

Reachable as `Integrate`RischTranscendental[f, x]` (explicit backtick form) and
`Integrate[f, x, Method -> "RischTranscendental"]`, and inserted into the `Integrate`
Automatic cascade **after** `Integrate`RischNorman`.

---

## 0. Build state (2026-07-11)

**Compiles clean** under the project's strict flags (`gcc -std=c99 -Wall -Wextra`,
no warnings from the module).  **Tests green:** the **extensive** suite in
`tests/test_integrate_risch_transcendental.c` — 20 `TEST` functions, ~290 assertions, one
per case / sub-case / sub-sub-case of the tower (see §7), plus a white-box unit test
of the `rt_rde_var_bound` degree arithmetic — passes in ~15 s; the
broader `integrals_tests` and `intrat_tests` suites are unaffected (the resultant
LRT reuses `intrat.c`'s log-part core additively, so the pure-rational path is
byte-identical); `valgrind` shows **no leak stack through the module** (the LRT
path adds zero blocks over the `Sin[1.0]` baseline: identical 13,440 B / 420
blocks, all macOS dyld/Accelerate).  Every branch is either an exact
decision procedure (SolveAlways / structural certificate) or a diff-back-verified
bounded search, so a mis-reduction can only *decline*, never emit a wrong closed
form — and **every** assertion in the suite was empirically classified against the
built integrator (no `WRONG` result across ~230 probed integrands) before being
pinned.

**Implemented (closes; else declines cleanly):**

| Area | Status |
|------|--------|
| Rational base case (delegated Bronstein rational Risch) | ✅ |
| Single-kernel logarithmic polynomial (`P(x, Log u)`) | ✅ |
| Single-kernel exponential Laurent (`Σ p_i E^(i u)`), **incl. rational exponents `E^(1/x)` and rational coefficients** (Phase C) | ✅ |
| Fractional Rothstein–Trager log-part (log & exp kernels) | ✅ |
| Hermite reduction for repeated poles (log & exp kernels) | ✅ |
| Coupled hyperexponential (Laurent + Hermite + log, `θ=0` pole) | ✅ |
| Multi-kernel sum-of-exponentials (`E^x Sin[x]`, …) | ✅ |
| Nested log / nested exp towers (flat unified ansatz) | ✅ |
| **Genuine one-extension recursion** (mixed exp/log towers; rational lower-field coefficients) | ✅ |
| **Tower proper part** — Hermite + Rothstein–Trager, **log top** (`rt_field_ratint`) and **exp top** (`rt_field_hyperexp_coupled`) | ✅ |
| **Field Risch DE** in the lower field, **arbitrary (non-monomial) denominator** via the denominator theorem (`rt_field_rde`) | ✅ |
| Trig / hyperbolic front-end (`TrigToExp` → exp machinery → `ExpToTrig`) | ✅ |
| **Evaluator-merged exponential monomials** (`E^x E^(E^x) → E^(x+E^x)` re-split into an independent tower basis) | ✅ |
| **Multiplicatively commensurate merged kernels** (`E^(2 E^x) = (E^(E^x))^2` reduced to a primitive + integer powers in `rt_tower_build`) | ✅ |
| **Single-kernel pure resultant LRT** (algebraic residues → `ArcTan`/`Log` via `Res_t` + Rioboo `LogToReal`, `rt_frac_lrt`) | ✅ |
| **Tower-level resultant LRT**, **log top** (`rt_field_lrt_logpart` → `rt_field_ratint`; nested-log algebraic residues, e.g. `1/(x Log[x](Log[Log[x]]^2+1)) → ArcTan[Log[Log[x]]]`) | ✅ |
| **Tower-level resultant LRT**, **exp top** (`rt_field_lrt_logpart` → `rt_field_hyperexp_coupled`; e.g. `E^x E^(E^x)/(1+E^(2 E^x)) → ArcTan[E^(E^x)]`) | ✅ |
| Special functions — `Erf`/`Erfi`, `ExpIntegralEi`, `LogIntegral`, `PolyLog` | ✅ |

**Not yet implemented** (integrands needing these decline; see §6 for detail):
The full Bronstein **SPDE** degree machinery is now **complete** — the **RDE-solver degree
bounds and the flat-tower / proper-part Hermite ansatz bounds are principled and cap-free**
(Bronstein `RdeBoundDegree` via `rt_rde_var_bound`; exact top-kernel log/exp Laurent
bounds; derived inner-kernel windows via `rt_var_mult_at_zero`; all degrees — §3.12,
§6.1 item 1), **including the leading-coefficient *cancellation/resonance* sub-case**
(Bronstein's recursive degree reduction), now folded into `rt_rde_var_bound` as a monotone
resonance-aware widening with live caller-side detection (`rt_resonance_int`) in the
exponential-top field RDE (§3.12).  Remaining: **Phase D** algebraic extensions
(`Sqrt`, `RootSum`); **Phase E** per-call flags to force/suppress the special
functions; and **Phase F** `Simplify` post-processing for the I-laden /
unfactored outputs.

**Multiplicatively commensurate** merged kernels (`E^(2 E^x) = (E^(E^x))^2`) are
now handled (tenth increment, §3.11): the tower builder reduces each
commensurability class of exp exponents to a single primitive and aliases the
rest to integer powers of its tower variable, so both the *additive* merge
(`E^x E^(E^x) → E^(x+E^x)`, ninth increment) and the *multiplicative* one close.

The resultant LRT is **lifted into the tower recursion** at a *logarithmic* top
(`rt_field_ratint`, §3.4) — a nested-log proper part with algebraic residues
such as `1/(x Log[x] (Log[Log[x]]^2+1)) → ArcTan[Log[Log[x]]]` closes — and, with
the commensurate reduction in place, at an *exponential* top too
(`rt_field_hyperexp_coupled`): `E^x E^(E^x)/(1+E^(2 E^x)) → ArcTan[E^(E^x)]`.

---

## 1. Design principles (non-negotiable)

1. **Recursive Risch, NOT parallel Risch.** This is the recursive
   decision-procedure algorithm (the Bronstein/Roach lineage). It is
   deliberately distinct from `Integrate`RischNorman`, which is
   the *parallel* Risch (Bronstein pmint / Norman) heuristic. RischTranscendental
   **never** falls back on the parallel-Risch engine.

2. **Correct by construction — no differentiation check.** A Risch integrator
   is a decision procedure. Every branch emits a closed form only behind an
   exact structural certificate (a logarithmic derivative that `Cancel`s to a
   polynomial; an exactly linear/coprime denominator; a `SolveAlways`-certified
   polynomial identity; ...) that already proves the answer. Results are **not**
   verified by differentiating them.

3. **Completeness over aesthetics.** A correct-by-construction branch is never
   suppressed because its output is not in the prettiest form. If the answer is
   correct, it ships; residual un-simplified shapes are flagged as `Simplify`
   improvement opportunities (see §5), not dropped.

4. **Grounded in Mathilda's `Expr`/poly/rat machinery**, not a separate CRE
   polynomial representation. The transcendental extension is handled by aliasing each
   kernel `theta` to a fresh polynomial variable `t` and running poly/rat
   primitives (`Together`, `Factor`, `PolynomialGCD`, `Cancel`, `SolveAlways`,
   `D`) over `{x, t}`, with the transcendental derivation applied explicitly as
   `D = d/dx + (D theta) d/dtheta`.

5. **Never wrong; decline instead.** Anything out of scope, or where a heuristic
   degree bound is insufficient, returns the integral **unevaluated** rather
   than a wrong or garbage result.

---

## 2. Dispatch order

`builtin_rischtranscendental` -> `rt_integrate`:

```
1. rt_rational_case         rational in x            -> delegate Integrate`BronsteinRational
2. rt_transcendental_case:
     a. rt_log_poly_case     polynomial in Log[u]
     b. rt_exp_poly_case     Laurent polynomial in E^u   (Risch DE per power)
     c. rt_frac_case         squarefree fractional log-part (Rothstein-Trager), log & exp kernels
                             — rational residues via SolveAlways, then rt_frac_lrt
                               (resultant LRT) for algebraic residues -> ArcTan/Log
     d. rt_hermite_case      repeated poles (Hermite),  log & exp kernels
     e. rt_hyperexp_case     coupled hyperexponential (mixed poly + log, Laurent-coupled Hermite)
     f. rt_expsum_case       sum of non-commensurate exponentials E^(W_k)  (decoupled Risch DE)
     g. rt_log_tower_case    nested logarithmic tower (depth > 1) over triangular D_tower
     h. rt_exp_tower_case    nested exponential tower (depth > 1), Laurent ansatz
     i. rt_recursive_tower_case  genuine one-extension recursion (mixed towers,
                                 rational lower-field coefficients)
     j. rt_trig_frontend     TrigToExp -> exp machinery (+ expsum) -> ExpToTrig
3. rt_special_case          Erf/Erfi, ExpIntegralEi, LogIntegral, PolyLog (special-function outputs)
```

Cheaper/more specific cases run first; the general `SolveAlways`-ansatz cases
(`hermite`, `hyperexp`) run later so the easy integrands never pay for them.

---

## 3. Implemented cases

### 3.1 Rational base case  (`rt_rational_case`)
Delegated to `Integrate`BronsteinRational` (Mathilda's recursive rational Risch:
Hermite + Lazard–Rioboo–Trager). This is the rational Risch base case.

### 3.2 Logarithmic-polynomial  (`rt_log_poly_case`)
`P(x, Log[u])` polynomial in `theta = Log[u]`, `u` rational in `x`, integrated by
the recursive primitive-polynomial **coefficient matching**
`q_i' + (i+1) q_{i+1} eta = p_i` (top-down), with a **limited-integration
oracle** (`rt_limited_integrate`) that folds a would-be new logarithm back into
the tower by bumping `q_{i+1} += c_i/(i+1)`.

- `Integrate[Log[x], x] = x Log[x] - x`
- `Integrate[Log[x]^2, x]`, `Integrate[x Log[x], x]`, `Integrate[Log[x]^3, x]`
- `Integrate[Log[x^2+1], x] = x Log[1+x^2] - 2x + 2 ArcTan[x]`
- `Integrate[Log[2x+3], x]`, `Integrate[Log[x]/x, x] = Log[x]^2/2`

### 3.3 Exponential Laurent-polynomial  (`rt_exp_poly_case`)
`sum_i p_i(x) E^(i u)`, single primitive kernel `E^u`, `i`
positive/negative/multiple. Powers decouple: the `i = 0` term integrates in
`K`, each `i != 0` term solves the **Risch differential equation**
`q_i' + i u' q_i = p_i`.  When `u` and `p_i` are polynomial in `x`, a bounded
polynomial ansatz solved exactly with `SolveAlways` (`rt_solve_rde`); when `u` (the
exponent) OR `p_i` is **rational** in `x` — Phase C, e.g. `E^(1/x)` — the
denominator-theorem ansatz `q_i = h/Denominator[p_i]` with a bounded polynomial `h`
(`rt_solve_rde_rational`), gated so `p_i`/`u'` must be genuine rational functions of
`x` (a transcendental coefficient like the `Sin[x]` of `E^x Sin[x]` is rejected and
left to the expsum / trig front-end).

- `Integrate[x E^x, x] = (x-1) E^x`, `Integrate[x^2 E^x, x]`
- `Integrate[x E^(x^2), x] = E^(x^2)/2`, `Integrate[x E^(-x^2), x] = -E^(-x^2)/2`
- `Integrate[E^(2x), x]`, and Laurent `Integrate[(E^x+E^(-x))/2, x] = Sinh[x]`,
  `Integrate[E^x + E^(-x) + E^(2x), x]`
- **Rational exponent / coefficient (Phase C):** `Integrate[-E^(1/x)/x^2, x] =
  E^(1/x)`, `Integrate[E^x/x - E^x/x^2, x] = E^x/x`.
- When the RDE has no polynomial solution (`N < 0`) the case declines
  (so `E^(-x^2)` itself falls through to the Erf recognizer); `E^(1/x)` (no rational
  RDE solution) likewise declines.

### 3.4 Fractional Rothstein–Trager log-part  (`rt_frac_case` / `rt_frac_try`)
Proper rational in `theta` with **squarefree** denominator `d = prod g_i`,
integrates to `sum_i c_i D(g_i)/g_i` with **constant** residues found from the
exact identity `num = sum_i c_i D(g_i) (d/g_i)` solved for all `theta` and `x`
via `SolveAlways[..., {t, x}]`. Both log and exp kernels; `D = d/dx +
(D theta) d/dtheta`.

- `Integrate[1/(x (1+Log[x])), x] = Log[1+Log[x]]`
- `Integrate[E^x/(1+E^x), x] = Log[1+E^x]`
- `Integrate[1/(x Log[x] (1+Log[x])), x] = Log[Log[x]] - Log[1+Log[x]]`
- `Integrate[1/(x (Log[x]^2-1)), x]` (two residues)

**Pure resultant Lazard–Rioboo–Trager (`rt_frac_lrt`, 2026-07-11).**  The
single-constant-per-factor `SolveAlways` ansatz above can only express
**rational** residues; an irreducible-over-Q denominator factor in `theta` whose
Rothstein–Trager residues are **algebraic** (e.g. the `±I/2` residues that split
`Log[x]^2+1` into a conjugate log pair) has no such solution, so `rt_frac_try`
declines.  When it does, `rt_frac_lrt` runs the exact resultant LRT: the residues
are the roots of `Res_t(a - z D(d), d)` (the monomial derivation `D(d)`, residue
variable `z`), and each squarefree factor `Q_i(z)` with its log argument
`S_i(z,t)` is collapsed to real `Log + ArcTan` form by Rioboo's `LogToReal`.  The
exact resultant + subresultant-PRS log-argument extraction + `LogToReal` are
delegated to the new builtin `Integrate`TranscendentalLogPart` in `intrat.c`,
which **reuses the tested rational-LRT machinery** (only `D[d,x]` → `D(d)` and the
elimination variable change, plus an x-content strip so the residues are gated
free of `x` — a residue that depends on `x` means non-elementary → decline).
The exp kernel is `rt_exp_kernelize`d first so `E^(2x)` collapses onto `t^2`.
Because the reuse crosses content/denominator boundaries, this path is **diff-back
VERIFIED** (`rt_verify_antideriv`) — a mis-reduction declines, never ships wrong.

- `Integrate[1/(x (Log[x]^2+1)), x] = ArcTan[Log[x]]`
- `Integrate[1/(x (Log[x]^2+4)), x] = ArcTan[Log[x]/2]/2`
- `Integrate[E^x/(E^(2x)+1), x] = ArcTan[E^x]`
- `Integrate[(1+Log[x])/(x (Log[x]^2+1)), x] = ArcTan[Log[x]] + Log[1+Log[x]^2]/2`
  (mixed numerator: LogToReal returns the full real form)
- `Integrate[1/(x (Log[x]^3+1)), x]` (linear + irreducible-quadratic factors →
  `Log` + `Log`/`ArcTan` mix)
- Non-elementary siblings (`1/(Log[x]^2+1)`, no `1/x`; `1/(x^2 (Log[x]^2+1))`,
  x-dependent residues) decline via the x-content gate.

### 3.5 Hermite reduction (repeated poles)  (`rt_hermite_case` / `rt_hermite_try`)
Proper rational in `theta` with a **repeated** pole reduces to
`Q = H(theta)/Hden(theta) + sum_j c_j Log(g_j)`, where
`Hden = PolynomialGCD[D, dD/dtheta]` (`= prod D_m^(m-1)`, robustly obtained —
avoids parsing multiplicities), `g_j` = squarefree factors of `D/Hden`, `H` a
polynomial in `theta` of degree `< deg(Hden)` with bounded-degree
polynomial-in-`x` coefficients. All constants solved by
`SolveAlways[Q'-f==0, {t,x}]`; result `Cancel`ed.

- Log kernel: `Integrate[1/(x (1+Log[x])^2), x] = -1/(1+Log[x])`,
  `Integrate[1/(x Log[x]^2), x] = -1/Log[x]`,
  `Integrate[(1+Log[x])/(x Log[x]^2), x] = Log[Log[x]] - 1/Log[x]` (Hermite + log)
- Exp kernel (when `D` coprime to `theta`):
  `Integrate[E^x/(1+E^x)^2, x] = -1/(1+E^x)`, `Integrate[E^x/(2+E^x)^2, x]`

### 3.6 Coupled hyperexponential  (`rt_hyperexp_case`)
General rational function of `theta = E^u` whose integral **mixes** a
Laurent-polynomial part with a log part (where `D(g)/g` is itself improper in
`theta`, so the parts do not separate). Unified ansatz
`Q = sum_i w_i(x) theta^i + H(theta)/Hden(theta) + sum_j c_j Log(g_j)`
(bounded-degree polynomial `w_i`, Hermite numerator `H` proper in `theta` with
polynomial-in-`x` coefficients, constant `c_j`) solved at once by
`SolveAlways[Q'-f==0, {t,x}]`. The `theta`-coprime denominator is split into its
repeated part `Hden = PolynomialGCD[Dtil, dDtil/dtheta]` (absorbed by the
Hermite term `H/Hden`) and its squarefree radical `rad = Dtil/Hden` (whose
distinct factors give the log terms); when `Dtil` is squarefree `Hden = 1` and
the ansatz collapses to the plain Laurent + log form.

- `Integrate[1/(1+E^x), x] = x - Log[1+E^x]`, `Integrate[1/(2+E^x), x]`
- `Integrate[E^(2x)/(1+E^x), x]`
- **Repeated / `theta = 0` poles (Phase A, Laurent-coupled Hermite):**
  `Integrate[1/(1+E^x)^2, x] = x + 1/(1+E^x) - Log[1+E^x]`,
  `Integrate[1/(E^x (1+E^x)^2), x]` (a `theta = 0` pole coupled with a repeated
  pole at `theta = -1`), `Integrate[1/(2+E^x)^2, x]`,
  `Integrate[1/(1+E^x)^3, x]` (higher repeated pole)

### 3.9 Multi-kernel exponential sums  (`rt_expsum_case`) — Phase B, first increment
A sum `f = sum_k p_k(x) E^(W_k)` of exponentials whose exponents `W_k` are **not**
integer multiples of a single primitive (so `rt_exp_kernelize` cannot reduce them
to one `t`) — e.g. the `(1 ± I) x` pair from `E^x Sin[x]` after `TrigToExp`. The
distinct exponentials are independent transcendental extensions; because `d/dx`
maps each `E^(W_k)` to a multiple of itself and never mixes distinct exponents,
the terms **decouple**: `INT p_k E^(W_k) dx = q_k E^(W_k)` with `q_k` the exact
solution of the Risch DE `q_k' + W_k' q_k = p_k` (`rt_solve_rde`, `W_k` polynomial
in `x`), plus any `W_k = 0` term integrated in `K`. The term splitter descends
into nested (un-flattened) `Times` so `p_k` is the full exponential-free cofactor.
Correct by construction (each `q_k` is `SolveAlways`-certified, and
`d/dx sum_k q_k E^(W_k) = f` by linearity); declines whole if any term is
non-elementary (`E^(x^2)`), non-polynomial exponent (Phase C), or carries a
residual log. Reached directly and inside the trig front-end.

- `Integrate[E^x Sin[x], x]`, `Integrate[E^x Cos[x], x]`,
  `Integrate[x E^x Sin[x], x]`, `Integrate[E^(2x) Sin[3x], x]`,
  `Integrate[x^2 E^x Cos[x], x]`
- **Known `Simplify` gap** (same family as §5.1): through the complex
  exponentials `E^((a ± b I) x)`, `ExpToTrig` leaves the answer in a correct but
  I-laden `Cosh`/`Sinh`-of-complex form; the diff-back is exactly `0`.

### 3.10 Nested logarithmic tower  (`rt_log_tower_case`) — Phase B, second increment
A rational function of a **chain of nested logarithms**
`t_1 = Log[u_1](x)`, `t_2 = Log[u_2](x, t_1)`, ..., `t_n = Log[u_n](x, t_1..t_{n-1})`
(`n >= 2`), which the single-extension cases (one kernel over `C(x)`) cannot
model. The single-kernel derivation `D = d/dx + Dt d/dt` generalizes to the whole
tower:
```
D_tower = d/dx + sum_i Dt_i d/dt_i,   Dt_i = Cancel[D[u_i,x]/u_i] |_{kernels->t}
```
(triangular: each `Dt_i` lies in `C(x, t_1..t_{i-1})`, checked explicitly). One
unified ansatz
```
Q = sum_{k=0}^{Ntop} P_k(x, t_1..t_{n-1}) t_n^k  +  sum_j c_j Log(g_j)
```
(`P_k` bounded-degree polynomials with unknown constant coefficients; `g_j` the
squarefree `t_n`-denominator factors; `c_j` constants) is solved at once by
`SolveAlways[D_tower[Q] - f == 0, {t_n,...,t_1,x}]`. Correct by construction — the
`t_i` are algebraically independent transcendentals, so a certified polynomial
identity in them proves the functional identity.

- `Integrate[1/(x Log[x] Log[Log[x]]), x] = Log[Log[Log[x]]]`
- `Integrate[Log[Log[x]]/(x Log[x]), x] = Log[Log[x]]^2/2`
- `Integrate[Log[Log[x]]/x, x] = Log[x] Log[Log[x]] - Log[x]` (lower-field `t_1`
  coefficient), depth-3 `1/(x Log[x] Log[Log[x]] Log[Log[Log[x]]])`, and
  independent (non-nested) chains such as `Log[x] + Log[x+1]`.

**Safety gate (never wrong; decline instead).** After substituting all kernels to
`t`, the integrand must be a genuine rational function of the *whole* tower
`{x, t_1..t_n}` — enforced by `PolynomialQ[num, {x,t_1..t_n}]` and the same for
`den`. Without it a residual non-rational kernel of an inner argument (e.g.
`Sin[Log[x]] -> Sin[t_1]`) would pass the per-top-variable test and let the ansatz
certify a wrong closed form (it once emitted `Integrate[Sin[Log[x]] Log[Log[x]]] =
0`). Genuinely non-elementary nested integrands (`Log[x] Log[Log[x]]`, which needs
`ExpIntegralEi`) also decline. Repeated `t_n` poles (tower Hermite) and rational
lower-field coefficients (the full recursion) are out of scope for this increment
and decline cleanly.

### 3.11 Nested exponential tower  (`rt_exp_tower_case`) — Phase B, third increment
Dual of §3.10 for a chain of nested **exponentials** `t_i = E^(u_i)` with `u_i` in
`C(x, t_1..t_{i-1})` (e.g. `t_1 = E^x`, `t_2 = E^(E^x)`). An exponential kernel's
derivative is a multiple of itself, so `Dt_i = (D[u_i,x]|_{ker->t}) t_i`, and the
ansatz is **Laurent** — the inner exp kernels are invertible, so their coefficient
exponents range over negatives too (crucial: `INT E^(x+E^x) dx = E^(E^x) = t_2/t_1`
needs a `t_1^{-1}`). `Q = sum_{i=ilo}^{ihi} P_i(x, t_1..t_{n-1}) t_n^i + sum_j
c_j Log(g_j)` solved over `{t_n,...,t_1,x}`. Closes
`Integrate[E^x E^(E^x), x] = E^(E^x)`, `Integrate[E^(2 E^x) E^x, x]`.

**Both tower cases are diff-back VERIFIED** (`rt_verify_antideriv`): unlike the
single-kernel decision-procedure cases, the tower cases are bounded-ansatz
searches over a multi-variable `SolveAlways`, whose certificate can be spurious
for a non-elementary integrand — so `Simplify[D[result,x]-f]` must be exactly `0`
or the case declines (e.g. the non-elementary `E^(E^x)/(1+E^(E^x))` is rejected).

**Single-kernel nesting gate (`rt_kernel_simple`).** A pre-existing correctness
bug: `rt_frac_case` (and the other single-kernel cases) substituted only ONE
exponential kernel, leaving a nested inner kernel (e.g. `E^x` inside `E^(E^x)`) in
the derivation `Dt = u' theta`; `SolveAlways` then treated it as a free parameter
and certified the WRONG residue `Integrate[E^(E^x)/(1+E^(E^x))] = Log[1+E^(E^x)]/E^x`.
Fixed: every single-kernel case now requires its kernel's defining `u` to be
rational in `x` alone (`rt_find_exp_of_x(u,x)==NULL && rt_find_log_of_x(u,x)==NULL`);
nested integrands fall through to the tower cases.

### 3.12 Genuine one-extension recursion  (`rt_recursive_tower_case`) — Phase B, fourth increment
The flat tower cases (§3.10, §3.11) solve ONE `SolveAlways` over every tower
variable at once, with bounded-degree **polynomial** lower-field coefficients.
That structurally cannot express two whole families: **mixed exp/log towers** (each
flat case is single-kind) and **rational lower-field coefficients** (a `t_n`
coefficient that is `1/x`, a nonlinear unknown for the flat ansatz).  This engine
is the genuine recursive Risch (Bronstein ch. 5): it peels
**one extension at a time**.

- `rt_tower_build` constructs the ordered differential tower (`RtTower`): collect
  every x-dependent `Log`/`E^` kernel, order innermost-first by structural
  containment (tie-break EXP-deeper so the primitive recursion sits on top and the
  exponential Risch DEs bottom out in `C(x)`), assign tower variables `t_i`, and
  compute each derivation coefficient `Dcoef_i` (log: `u_i'/u_i`; exp: `w_i'`).  The
  **structure-theorem soundness check** requires every `Dcoef_i` to lie in
  `K_{i-1} = C(x, t_1..t_{i-1})` (triangular: free of `t_i..t_n` and of any residual
  foreign kernel), else the tower is rejected.
- `rt_field_integrate(F, L)` integrates `F` in `K_L`: `L < 0` is the rational base
  case (`Integrate`BronsteinRational`).  Otherwise it splits `F` in the top kernel
  `t_L` and integrates the polynomial/Laurent part **coefficient by coefficient**,
  each coefficient integral being an integration in the LOWER field `K_{L-1}` (the
  recursion).  A **logarithmic** top uses the primitive-polynomial coefficient
  recursion `q_i' + (i+1) q_{i+1} Dt = p_i` (`rt_int_primitive_poly`), with each
  residual integrated by `rt_limited_field_integrate` (recursion + new-log fold-back
  of the single would-be new logarithm `= t_L`).  An **exponential** top uses the
  hyperexponential Laurent split (`rt_int_hyperexp_poly`): the `i = 0` term recurses,
  each `i != 0` term solves the Risch DE `q_i' + i w_L' q_i = p_i` (`rt_field_rde`,
  base case only this increment).
- Because each coefficient is integrated *in its own field* it may be **rational**,
  and because each level dispatches on its own kernel kind the tower may **mix** logs
  and exponentials.
- Correctness: a bounded search (not a decision procedure), so **diff-back verified**
  (`rt_verify_antideriv`); a spurious multi-variable certificate cannot ship.

Closes `Integrate[E^x/x + E^x Log[x], x] = E^x Log[x]` (mixed independent),
`Integrate[Log[1+E^x] + x E^x/(1+E^x), x] = x Log[1+E^x]` (nested log-over-exp),
`Integrate[1/(x^2 Log[x]) - Log[Log[x]]/x^2, x] = Log[Log[x]]/x` (rational
coefficient).  Non-elementary mixed integrands (`E^x Log[x]`, `Log[x] Log[Log[x]]`)
decline.

**Proper rational part (fifth increment, 2026-07-11) — `rt_field_ratint`.** A
logarithmic top level with a genuine `t_n`-pole (proper fraction) is integrated
instead of declining: `rt_field_integrate` splits `F` into its `t_n`-polynomial
part (the primitive recursion) and its proper part, and `rt_field_ratint`
integrates the proper part to `Q = H(t_n)/Hden + sum_j c_j Log(g_j)` where
`Hden = gcd(den, d den/dt_n)` is the repeated part (**tower Hermite**), `g_j` the
distinct `t_n`-factors of the squarefree radical `den/Hden` (**Rothstein–Trager**
log arguments), `H` a `t_n`-polynomial of degree `< deg(Hden)` with bounded
lower-field coefficients, and `c_j` constants — all solved at once by
`D_tower[Q] = num/den` (`rt_tower_deriv`) over every tower variable via
`SolveAlways`.  The constant-residue requirement is automatic (a non-constant
residue has no constant-symbol solution ⇒ declines).  The Hermite fraction is
`Cancel`ed for a clean form.  Closes `1/(x Log[x] (1+Log[Log[x]])^2)` →
`-1/(1+Log[Log[x]])` (repeated tower pole), `1/(x(1+Log[x])) + E^x` →
`E^x + Log[1+Log[x]]` (poly-plus-proper split), `Log[Log[x]]/(x Log[x]
(1+Log[Log[x]]))` → `Log[Log[x]] - Log[1+Log[Log[x]]]` (split with new-log
fold-back).  `1/(Log[x](1+Log[Log[x]]))` (non-constant residue) declines.

**Tower-level resultant LRT (2026-07-11) — `rt_field_lrt_logpart`.** The
`SolveAlways` ansatz above expresses exactly one *constant* residue per
squarefree `t_n`-factor, so a factor whose Rothstein–Trager residues are
*algebraic* (e.g. the `±I/2` that split `t_n^2+1` into a conjugate log pair =
`ArcTan`) has no constant-symbol solution and declines.  When it does,
`rt_field_ratint` (and `rt_field_hyperexp_coupled`) fall back to the exact
resultant LRT — the same `Integrate`TranscendentalLogPart` used by the
single-kernel `rt_frac_lrt` (§3.4), but now with the *tower* derivation
`D_tower[rad]` and a **gate list** `{x, t_0, …, t_{L-1}}` requiring the residues
to be constants of the whole tower derivation (free of every lower-field
variable), after stripping the lower-field content that clearing the tower
denominators leaves in the resultant.  Only the squarefree case (no repeated
`t_n`-pole) is handled; the top-level caller diff-back verifies.  Closes the
nested-log algebraic-residue class: `1/(x Log[x] (Log[Log[x]]^2+1)) →
ArcTan[Log[Log[x]]]`, `1/(x Log[x] (Log[Log[x]]^2+4))`, mixed `ArcTan`+`Log`
numerators, and cubic denominators (`Log` + `Log`/`ArcTan` mix).
`1/(Log[x](Log[Log[x]]^2+1))` (residues depend on `x`) declines via the gate.
The exponential-top fallback is wired identically but its degree-≥2 targets need
the commensurate-exponent kernelization of §6.1 item 3 first.

**Coupled hyperexponential proper part (sixth increment, 2026-07-11) —
`rt_field_hyperexp_coupled`.** The exponential-top counterpart of the fifth
increment.  For an exponential top kernel the Laurent and log parts do NOT separate
(`D_tower[Log(1+t_n)] = w_n' t_n/(1+t_n)` feeds both the `t_n^0` Laurent coefficient
and the proper part), so the log/poly split used for the logarithmic case is
invalid.  When `rt_int_hyperexp_poly` declines because a genuine `t_n`-pole is
present, this solves ONE unified ansatz
`Q = sum_i w_i t_n^i + H(t_n)/Hden + sum_j c_j Log(g_j)` (`w_i`/`H` bounded
lower-field polynomials, `c_j` constants) against `D_tower[Q] = num/den` over every
tower variable — the single-kernel coupled hyperexponential case lifted to the tower
derivation.  Closes `(2 Log[x]/x) E^(Log[x]^2)/(1+E^(Log[x]^2))` →
`Log[1+E^(Log[x]^2)]` (nested exp-over-log); the non-elementary `1/(1+E^(Log[x]^2))`
declines.

**General field Risch DE (seventh increment, 2026-07-11) — `rt_field_rde` general
branch.** The exponential Laurent step (`rt_int_hyperexp_poly`) solves, per power
`i != 0`, the Risch DE `q_i' + i w_n' q_i = p_i` for `q_i` in the lower field.  Beyond
the base case (`C(x)`, polynomial `q_i` via `rt_solve_rde`), a `q_i` that is
**rational** in the lower field (monomial denominator, e.g. `1/Log[x]`) — which the
coupled-hyperexponential *polynomial* coefficients cannot express — is now found by a
bounded ansatz solved by `SolveAlways` of `D_tower[q_i] + i w_n' q_i = p_i`.  Closes
`(2/x - 1/(x Log[x]^2)) E^(Log[x]^2)` → `E^(Log[x]^2)/Log[x]`; the non-elementary
`E^(Log[x]^2)` declines.

**Field RDE denominator theorem (eighth increment, 2026-07-11).**  Generalizes the
field RDE from a monomial-Laurent ansatz to an arbitrary denominator: by the RDE
denominator theorem `denom(q_i) | Denominator[p_i]` (a pole of `q_i` gives a
higher-order pole in `q_i' + i w_n' q_i = p_i` that nothing cancels), so
`q_i = h/Denominator[p_i]` with `h` a bounded *polynomial* numerator over the lower
field — subsuming the monomial case and adding non-monomial denominators.  Closes
`(2 Log[x]/(x(1+Log[x])) - 1/(x(1+Log[x])^2)) E^(Log[x]^2)` →
`E^(Log[x]^2)/(1+Log[x])` (`q = 1/(1+Log[x])`).  The `h/pd` result is `Cancel`ed.

**Principled RDE degree bounds — no arbitrary caps (2026-07-11) — SPDE §6.1 item 1.**
The RDE solvers previously carried *arbitrary degree caps* — `rt_field_rde` used a loose
`deg(pd)+deg(pn)+1` proxy **capped at 5** plus a `nmono ≤ 60` monomial ceiling, and
`rt_solve_rde_rational` used `N = deg(pd)+deg(pn)+2` **capped at 10** — so a
higher-degree but perfectly solvable coefficient declined.  Those caps and ceilings are
removed and replaced by Bronstein's exact leading-degree bound (`RdeBoundDegree`), shared
in `rt_rde_var_bound(deg_v(p), deg_v(f), deriv_lowers)`.  For a Risch DE `D[q] + f q = p`,
matching the top `v`-degree:
- **deriv-lowering `v`** (base `x` or a *logarithmic* monomial, `D` lowers `deg_v`):
  `deg_v(f) ≥ 0` → `deg_v(q) = deg_v(p) − deg_v(f)` (exact, `f q` dominates); else
  `deg_v(q) = max(deg_v(p)+1, deg_v(p)−deg_v(f))` (integration rise vs a pole in `f`).
- **deriv-preserving `v`** (an *exponential* monomial, `D[t^k] = k w' t^k`):
  `deg_v(q) = deg_v(p) − deg_v(f)`.

Then `deg_v(h) = deg_v(q) + deg_v(pd)` for `q = h/pd`.  The bound is a function of the
equation's degrees alone — no magic constants — so the RDE solver works for **all
degrees**.  Correctness is unaffected: every solution is `SolveAlways`-certified and the
enclosing tower case diff-back verified, so a bound can only ever *decline*, never ship a
wrong closed form.  Closes exponential-Laurent coefficients of arbitrary degree
(`(6 Log[x]^5 + 2 Log[x]^7)/x · E^(Log[x]^2) → Log[x]^6 E^(Log[x]^2)`, and the deg-8 /
-12 / -20 analogues), and keeps `-E^(1/x)/x^2 → E^(1/x)` exact; the non-elementary
`E^(Log[x]^2)`, `E^(1/x)` still decline.

**Leading-coefficient cancellation / resonance (2026-07-11) — `rt_resonance_int` +
`rt_rde_var_bound` widening.**  The leading-degree balance above is exact except where the
two leading coefficients of `D[q]` and `f q` *cancel*, letting `deg_v(q)` exceed the naive
value (Bronstein's recursive degree reduction).  This is now handled at the bound level in
two configurations, keyed by the integer `m_res` (the Bronstein resonance integer, or `-1`
when none): an **exponential** monomial at `deg_v(f) = 0` cancels when
`n = -(i·Dcoef_L)/Dcoef_v` is a nonnegative integer (`rt_resonance_int` computes it live in
the exp-top field RDE), and a **primitive** monomial at `deg_v(f) = -1` cancels at the
homogeneous-solution degree; in each case the bound is widened to `max(naive, m_res)`.  The
widening is **monotone** (only ever *raises* the bound), so it can never ship a wrong result
— every solution stays `SolveAlways`-certified and the enclosing case diff-back verified —
and a spurious `m_res` at worst adds unused ansatz terms.

*Reachability.*  In the current tower both cancellation configurations are provably
pre-empted, so `m_res = -1` on every reachable call: (a) an exponential resonance
`n = -(i w_L')/w_j'` being an integer means the top and lower exp exponents are
**commensurate**, and the commensurate-exponent reduction in `rt_tower_build` (§3.11)
collapses such kernels to one primitive *before* any RDE solve; (b) `deg_v(f) = -1` needs a
*simple* pole in `f = i·Dcoef`, which a rational tower element's derivative cannot have (it
would integrate to a `Log`) and the only kernel that could — a log top with `Dcoef = u'/u`
— never routes through the field RDE (it uses the primitive-polynomial recursion).  The
detection is computed live regardless, so the degree bound is exact per Bronstein should a
future kernel type ever expose either configuration.  The pure bound arithmetic — every
branch plus the resonance widening — is unit-tested directly (`test_rde_var_bound`,
`rt_rde_var_bound` exposed in the header) without needing an integrand to reach each config.

**Cap-free Hermite / flat-tower ansätze (2026-07-11) — SPDE §6.1 item 1, companion.**
The seven remaining flat-ansatz sites that truncated a *derived* bound with a magic
constant (`if(d>N)d=N`, `ihi>4→4`, `ilo<−4→−4`, the hardcoded inner-exp window
`t_i ∈ [−2,2]`) or gated `SolveAlways` with a magic completeness ceiling
(`nunk ≤ 60/80`) are now cap-free — `rt_log_tower_case`, `rt_exp_tower_case`, the
single-kernel Hermite and hyperexponential proper parts, `rt_field_hyperexp_coupled`,
and the hyperexponential Hermite reduction.  Top-kernel bounds are exact
(log: `Ntop = deg_top(f)+1`; exp: Laurent `[−mult_top(den), deg_top(num)−deg_top(den)]`);
inner exp-kernel Laurent windows are each kernel's own extent widened by the reach of the
top derivation coefficient `w' = D[u_n]` (derived via a new `rt_var_mult_at_zero` helper,
the coupling behind `∫ E^(x+E^x) dx = E^(E^x)`); lower-field degree proxies keep their
derived `deg_v(num)+deg_v(den)+1` form uncapped; ceilings become the overflow-only
`nunk > 0` guard.  Same diff-back + `SolveAlways` certification, so decline-only.  New
anchors the old caps declined now close: `Log[Log[x]]^5/(x Log[x]) → Log[Log[x]]^6/6`
(`Ntop=6`), the `Ntop=8` analogue, and `E^x E^(6 E^x)/(1+E^(E^x))` (top Laurent deg 5);
the `E^(E^x)` / `E^(E^x)/(1+E^(E^x))` decline guards still hold.

**Evaluator-merged exponential monomials (ninth increment, 2026-07-11) —
`rt_expand_exp_sums` + `rt_subst_kernels`.**  The evaluator automatically folds a
product of exponentials into a single power with a summed exponent
(`E^x E^(E^x) → E^(x+E^x)`).  That merged exponent is not a valid tower monomial —
its argument `x + E^x` carries the foreign kernel `E^x`, so `rt_tower_build`'s
structure-theorem check rejects the tower and the case declines.  A structural
pre-pass `rt_expand_exp_sums` re-splits `E^(a+b) → E^a E^b` (built directly, *not*
evaluated, so the evaluator cannot re-merge it), restoring the independent basis
`{E^x, E^(E^x)}`.  The kernel→tower-variable aliasing is likewise done structurally
(`rt_subst_kernels`) instead of by an evaluated `ReplaceAll`, which would re-merge
the split product before substitution.  Closes
`Integrate[E^x E^(E^x)/(1+E^(E^x)), x] = Log[1+E^(E^x)]` (through the coupled
hyperexponential proper part).  The non-elementary sibling `E^(E^x)/(1+E^(E^x))`
still declines (no `E^x` factor to supply the derivation).

**Multiplicatively commensurate merged kernels (tenth increment, 2026-07-11) —
`rt_tower_build` commensurate reduction.**  Two collected exponents `w_i`, `w_j`
whose ratio `w_i/w_j` is a nonzero rational define *algebraically dependent*
kernels — `E^w = (E^p)^k` for a class primitive `p` and integer `k` (e.g.
`E^(2 E^x) = (E^(E^x))^2`).  Left unreduced, the dependent kernel spuriously adds a
tower extension, breaking independence, and the case declines.  The tower builder
now partitions the exp exponents into commensurability classes, keeps only a
primitive (whose integer multiples cover its class) as a tower variable, and
aliases each remaining member `E^w` to the integer power `t_p^k` of the
primitive's variable — new `RtTower.marg/mprim/mmult` fields populated in
`rt_tower_build` and consumed by `rt_subst_kernels`.  A class with no
integer-ratio primitive (e.g. exponents `E^x/2` and `E^x/3`) is out of scope and
declines the whole tower (never wrong).  Closes the merged-kernel class and, since
its targets carry a degree-≥2 exp denominator `E^(2 u) = t^2`, **unblocks the
exponential-top algebraic-residue LRT** (`rt_field_lrt_logpart`):
`Integrate[E^x E^(2 E^x)/(1+E^(E^x)), x] = E^(E^x) - Log[1+E^(E^x)]`,
`Integrate[E^x E^(E^x)/(1+E^(2 E^x)), x] = ArcTan[E^(E^x)]` (exp-top LRT),
`E^x E^(3 E^x)/(1+E^(E^x))`, `E^x E^(2 E^x)/(1+E^(2 E^x))`,
`E^x E^(E^x)/(4+E^(2 E^x))`.  Non-elementary siblings the reduction now lets the
tower *build* but which lack an `E^x` derivation factor
(`E^(E^x)/(1+E^(2 E^x))`, `E^(2 E^x)/(1+E^(E^x))`) still decline via the diff-back
gate.

**Still deferred (declines cleanly):** a genuinely *rational* Hermite numerator
coefficient; the recursive degree-reduction half of the Bronstein SPDE (bounded
polynomial numerator ansatz here); and **non-integer** multiplicatively
commensurate exponents (a class whose members have no common integer-ratio
primitive), which the reduction declines.

### 3.7 Trig / hyperbolic front-end  (`rt_trig_frontend`)
`TrigToExp` -> the exponential machinery -> `ExpToTrig`. Both rewrites exact.

- `Sin`, `Cos`, `Sinh`, `Cosh`, `Sin[x]^2`, `Cos[x]^2`, `Sin[x] Cos[x]`,
  `Sin[2x]`, `Cosh[x]^2`, `Tan[x]`, `Tanh[x]`
- **Known `Simplify` gap:** through the complex substitution `u = I x`,
  `Tan`/`Tanh` close to a correct but I-laden form such as
  `I x - Log[1 + E^(2 I x)]` (`= -Log[Cos[x]]`) that no current simplifier
  collapses to real closed form. The correct answer is returned (see §5).

### 3.8 Special functions  (`rt_special_case`)  — Erf / dilog / Ei / li outputs
Default-ON (WL-faithful). Each recognizer is exact given its structural
certificate.

- Gaussian `K E^(a x^2 + b x + c)` (`a != 0`) -> `Erf`/`Erfi`
  (certificate: `Cancel[f'/f]` is a degree-1 polynomial; `K' = f|_{x=0}`;
  includes the completing-the-square constant `E^(-b^2/(4a))`)
- `(M E^(a x + b))/(c x + d)` -> `ExpIntegralEi` — the exponential kernel `E^v`
  is extracted **directly** (not via `Together`), so a **negative** leading
  coefficient (`E^(-x)/x`, whose `E^(-x)` `Together` would push into the
  denominator, hiding the linear form) and a nonzero exponent constant `b`
  close uniformly. Certificate: `v` exactly linear, and the rational cofactor
  `R = Together[f E^(-v)]` has a constant numerator over a linear denominator.
- `c w^(p-1) w'/Log[w]` -> `c LogIntegral[w^p]` — a `Log[w]` kernel, integer
  `p >= 1`, constant `c` (certificate: `Together[f Log[w]/(w^(p-1) w')]` is a
  nonzero constant). Subsumes the bare `K/Log[x]` (p=1, w=x) and adds a
  **scaled/affine argument** (`1/Log[2x] -> LogIntegral[2x]/2`, p=1) and a
  **monomial numerator** (`x/Log[x] -> LogIntegral[x^2]`, p=2).
- `K Log[1 + p x] / x` -> `-K PolyLog[2, -p x]`

Examples: `Integrate[E^(-x^2), x] = (Sqrt[Pi]/2) Erf[x]`,
`Integrate[E^x/x, x] = ExpIntegralEi[x]`, `Integrate[E^(-x)/x, x] =
ExpIntegralEi[-x]`, `Integrate[1/Log[x], x] = LogIntegral[x]`,
`Integrate[1/Log[2 x], x] = LogIntegral[2 x]/2`, `Integrate[x/Log[x], x] =
LogIntegral[x^2]`, `Integrate[Log[1-x]/x, x] = -PolyLog[2, x]`.

**Rothstein–Trager constant-residue guard (frac case).**  Surfaced while
widening the li recognizer: `1/Log[c x + d]` (`d != 0`) is `LogIntegral[c x+d]`,
but `rt_frac_case` was wrongly closing it with a *polynomial-coefficient*
logarithm `(1/c)(c x+d) Log[Log[c x+d]]`.  Root cause: `SolveAlways` can return
an **x-dependent** pseudo-solution for a residual with a constant-term structure
(`SolveAlways[1 - k/(x+c) == 0, {t,x}]` returns `k -> x + c`, not `False`), and
`rt_frac_try` used it without checking the residues are constant.  Fixed by
enforcing that every Rothstein–Trager residue is free of `x` and `t`, else the
case declines (letting the correct li recognizer close it).  This is a genuine
never-wrong (§1.5) fix, not just a completeness win.

---

## 4. Key internal machinery

- `rt_eval_call` / `rt_eval_own` / `rt_eval1/2/3` — build-and-evaluate helpers
  (evaluate does NOT consume its argument; these own the constructed node).
- `rt_free_of_x`, `rt_free_of_head`, `rt_is_poly`, `rt_degree`
  (`= Length[CoefficientList] - 1`; Mathilda has no `Exponent` builtin),
  `rt_coeff`, `rt_is_zero` (`Together`-based, exact and cheap).
- `rt_template` — parse a template, `ReplaceAll` placeholder symbols, evaluate.
- `rt_find_log_of_x`, `rt_find_exp_of_x`, `rt_collect_exp_exponents`,
  `rt_exp_kernelize` (pick a primitive `u` so every `E^w` is `E^(k u)`,
  substitute `E^(k u) -> t^k`).
- `SolveAlways[eqn, {t, x}]` is the workhorse certificate: it returns constant
  solutions to a polynomial identity that must hold for all `t` and `x`, or **no
  solution** (which safely declines). Used by frac / hermite / hyperexp.

**Memory discipline:** the builtin borrows `res`/`f`/`x` (never frees them;
splices `expr_copy`). Every ansatz array, `SolveAlways` call, and intermediate
is freed on all paths. Verified valgrind-clean: no leak stack passes through the
module (the residual "definitely lost" is the documented macOS dyld/Accelerate
baseline plus pre-existing `Solve`/poly-machinery leaks triggered by legitimate
use).

---

## 5. Known gaps / `Simplify` improvement opportunities

1. **`Tan`/`Tanh` (and similar) real closed form.** Through `u = I x` they close
   to `I x - Log[1 + E^(2 I x)]` etc., which is `-Log[Cos[x]]` but is not
   reduced by `Simplify` / `FullSimplify` / `TrigReduce` / `PowerExpand`.
   *Fix belongs in `Simplify`*: teach it the log-of-product / half-angle
   collapse `I x - Log[1 + E^(2 I x)] -> -Log[Cos[x]]`. The integrator already
   returns the correct answer.

2. **Occasional unreduced denominators**, e.g. `Integrate[1/(x (1+Log[x])^3), x]`
   returns `-1/(2 + 4 Log[x] + 2 Log[x]^2)` rather than `-1/(2 (1+Log[x])^2)`.
   Correct; a `Factor`-of-denominator pass in `Simplify` would tidy it.

3. ~~**Narrow special-function recognizers.**~~ **DONE (2026-07-11).** The three
   Ei/li integrands that previously declined now close: `E^(-x)/x ->
   ExpIntegralEi[-x]` (Ei with `a < 0`, via direct `E^v` extraction), `1/Log[2x]
   -> LogIntegral[2x]/2` (scaled argument) and `x/Log[x] -> LogIntegral[x^2]`
   (monomial numerator, both via `c LogIntegral[w^p]`).  See §3.8.  The widening
   also exposed and fixed a **wrong-answer** bug in the frac case (non-constant
   Rothstein–Trager residue for `1/Log[c x + d]`, `d != 0`) — §3.8.

---

## 6. Next phases

Ordered by value / tractability. All are correct-by-construction extensions;
until done, the corresponding integrands decline cleanly.

### Phase A — Exponential Hermite with a `theta = 0` pole (Laurent-coupled) ✅ DONE (2026-07-10)
Implemented by generalizing `rt_hyperexp_case`: the `theta`-coprime denominator
`Dtil = den/theta^a` is split into its repeated part
`Hden = PolynomialGCD[Dtil, dDtil/dtheta]` and squarefree radical
`rad = Dtil/Hden`, and a Hermite term `H(theta)/Hden(theta)` (proper in `theta`,
polynomial-in-`x` numerator coefficients) is fused into the Laurent + log ansatz
`Q = sum_i w_i theta^i + H(theta)/Hden(theta) + sum_j c_j Log(g_j)`, all
constants solved at once by one `SolveAlways[Q'-f==0, {t,x}]`. When `Dtil` is
squarefree `Hden = 1` and the ansatz collapses to the previous form, so existing
cases are unchanged. Closes the repeated / `theta = 0` exponential shapes
`1/(1+E^x)^2`, `1/(E^x (1+E^x)^2)`, `1/(2+E^x)^2`, `1/(1+E^x)^3`, ... which the
pure Hermite path (cannot supply the `x`/Laurent term) and the squarefree
hyperexp path both declined. Tests: `test_hyperexponential_case`.

### Phase B — Nested / multi-kernel towers (depth > 1, > 1 independent extension)
**First increment DONE (2026-07-10):** the decoupled **sum-of-exponentials** case
(`rt_expsum_case`, §3.9) — integrands that are a sum `sum_k p_k(x) E^(W_k)` of
exponentials with non-commensurate exponents `W_k` (polynomial in `x`), e.g.
`Integrate[E^x Sin[x], x]` (two non-commensurate `(1 ± I) x` after `TrigToExp`).
The distinct exponentials decouple, so each term is closed by its own Risch DE
without building a coupled tower.

**Second increment DONE (2026-07-10):** nested **logarithmic** towers
(`rt_log_tower_case`, §3.10) — a rational function of a chain of nested logarithms
`Log[x]` / `Log[Log[x]]` / ... is integrated over the triangular tower derivation
`D_tower` by one unified `SolveAlways` ansatz. Closes
`Integrate[1/(x Log[x] Log[Log[x]]), x]` and siblings (up to depth 4), guarded by
a whole-tower rationality gate so non-rational inner kernels decline rather than
certify a wrong answer.

**Third increment DONE (2026-07-10):** nested **exponential** towers
(`rt_exp_tower_case`, §3.11) via a Laurent ansatz over the exp tower derivation,
diff-back verified; plus the `rt_kernel_simple` fix closing a pre-existing
wrong-answer bug in the single-kernel cases for nested exponents.

**Fourth increment DONE (2026-07-11):** the **genuine one-extension-at-a-time
recursion** with a structure-theorem kernel basis (`rt_recursive_tower_case`,
§3.12).  Peels one kernel at a time, integrating the top-kernel polynomial/Laurent
part coefficient-by-coefficient with each coefficient integrated recursively in the
lower field.  Closes **mixed exp/log towers** and **rational lower-field
coefficients** — the two families the flat single-kind tower ansätze structurally
cannot express.  Diff-back verified.

**Fifth increment DONE (2026-07-11):** the **proper rational part** of the
recursion for a logarithmic top — tower Hermite reduction + Rothstein–Trager
(`rt_field_ratint`, §3.12): `H/Hden + sum c_j Log(g_j)` (constant residues).

**Sixth increment DONE (2026-07-11):** the **coupled hyperexponential proper part**
for an exponential top — `rt_field_hyperexp_coupled` (§3.12): one unified
Laurent + Hermite + log ansatz over the tower derivation (the exp Laurent and log
parts don't separate).

**Seventh increment DONE (2026-07-11):** the **general field Risch DE** — the
`rt_field_rde` general branch (§3.12): finds a `q_i` rational in the lower field.

**Eighth increment DONE (2026-07-11):** the **field-RDE denominator theorem** —
`q_i = h/Denominator[p_i]` with a bounded polynomial numerator (§3.12), extending the
field RDE to non-monomial denominators (e.g. `q = 1/(1+Log[x])`).

**Phase B is complete** — the one-extension recursion, both proper-part kinds, and
the field RDE with arbitrary denominators are all in place.  What remains within the
recursion is refinement, listed under §6.1.

### 6.1 Remaining work (prioritized)

The recursive architecture is in place; the remaining items are increasingly
specialized and each unlocks progressively rarer integrands.  All decline cleanly
until done.

1. **Full Bronstein SPDE (Ch. 6).** *First increment DONE (2026-07-11) — principled,
   cap-free RDE degree bounds.*  Both RDE solvers (`rt_field_rde`,
   `rt_solve_rde_rational`) had **arbitrary degree caps** (cap-at-5 + `nmono ≤ 60`;
   cap-at-10) that declined higher-degree but solvable coefficients.  Removed and
   replaced by Bronstein's exact leading-degree bound (`RdeBoundDegree`) in a shared
   `rt_rde_var_bound` — a function of the equation's degrees alone, monomial-type-aware
   (a logarithmic monomial's derivative lowers its degree, an exponential's preserves
   it), with **no magic constants and no ceiling**, so the RDE solver works for **all
   degrees** (§3.12).  Closes exponential-Laurent coefficients of arbitrary degree,
   `(6 Log[x]^5 + 2 Log[x]^7)/x · E^(Log[x]^2) → Log[x]^6 E^(Log[x]^2)` and its
   deg-8/-12/-20 analogues.  *Second increment DONE (2026-07-11) — cap-free
   Hermite / flat-tower ansätze.*  The same principled treatment now covers all seven
   flat-ansatz sites (`rt_log_tower_case`, `rt_exp_tower_case`, single-kernel Hermite +
   hyperexp, `rt_field_hyperexp_coupled`, hyperexp Hermite): exact top-kernel log/exp
   Laurent bounds, derived inner-exp windows (`rt_var_mult_at_zero` + top-coefficient
   reach), uncapped lower-field proxies, `nunk > 0` overflow guards replacing the
   `≤ 60/80` ceilings.  Closes `Log[Log[x]]^5/(x Log[x]) → Log[Log[x]]^6/6` (`Ntop=6`)
   and `E^x E^(6 E^x)/(1+E^(E^x))` (top Laurent deg 5).  *Third increment DONE
   (2026-07-11) — leading-coefficient cancellation/resonance.*  The last SPDE degree
   piece: `rt_rde_var_bound` now takes the Bronstein resonance integer `m_res` and, in
   the two cancellation configs (exponential `deg_v(f)=0`; primitive `deg_v(f)=-1`),
   widens the bound monotonically to `max(naive, m_res)`; the exp resonance integer
   `n = -(i·Dcoef_L)/Dcoef_v` is detected live by `rt_resonance_int` in the exp-top field
   RDE (§3.12).  The widening only ever *raises* the bound, so it is never-wrong by the
   same `SolveAlways` + diff-back gate.  Both configs are provably pre-empted in the
   current tower (exp resonance ⟺ commensurability, already reduced; a rational element's
   derivative has no simple pole), so `m_res = -1` on every reachable call — the detection
   hardens the bound to be exact per Bronstein for future kernel types, and the pure bound
   arithmetic (all branches + widening) is unit-tested directly (`test_rde_var_bound`).
   The Bronstein SPDE degree machinery is now **complete**.
2. ~~**Pure resultant Lazard–Rioboo–Trager.**~~ **DONE (2026-07-11).**  The
   single-kernel frac case now runs the exact resultant LRT (`rt_frac_lrt` →
   `Integrate`TranscendentalLogPart`, §3.4) when its `SolveAlways` rational-residue
   ansatz declines, closing integrands whose log part has **algebraic** residues:
   `1/(x (Log[x]^2+1)) → ArcTan[Log[x]]`, `E^x/(E^(2x)+1) → ArcTan[E^x]`, and mixed
   `Log`/`ArcTan` forms.  Reuses the tested rational-LRT resultant + subresultant-PRS
   + Rioboo `LogToReal` (only the derivation and elimination variable change), gated
   x-content-free and diff-back verified.  **Lifted into the tower recursion
   (2026-07-11):** `rt_field_lrt_logpart` runs the same resultant LRT for the
   tower proper part when the bounded `SolveAlways` ansatz declines, gating the
   residues free of *every* lower-field variable (`{x, t_0, …, t_{L-1}}`).  At a
   *logarithmic* top this closes nested-log algebraic residues, e.g.
   `1/(x Log[x] (Log[Log[x]]^2+1)) → ArcTan[Log[Log[x]]]`; with item 3's
   commensurate reduction in place, the *exponential*-top fallback now fires too:
   `E^x E^(E^x)/(1+E^(2 E^x)) → ArcTan[E^(E^x)]`.
3. ~~**Algebraically dependent kernels the evaluator merges.**~~ **DONE
   (2026-07-11).**  The *additive* merge (`E^x E^(E^x) → E^(x+E^x)`) is handled by
   the `rt_expand_exp_sums` pre-pass (ninth increment), and the *multiplicatively
   commensurate* case (`E^(2 E^x) = (E^(E^x))^2`) by the **commensurate-exponent
   reduction in `rt_tower_build`** (tenth increment, §3.11): each commensurability
   class of exp exponents keeps one primitive as a tower variable and aliases the
   rest to integer powers `t^k` of it.  Closes `E^x E^(2 E^x)/(1+E^(E^x))` and
   unblocks the exp-top LRT `E^x E^(E^x)/(1+E^(2 E^x)) → ArcTan[E^(E^x)]`.  Only
   *non-integer* commensurate exponents (no common integer-ratio primitive)
   remain out of scope.
4. **Phase D — algebraic extensions** (`Sqrt`, `n`-th roots, `RootSum`).  A separate,
   large subsystem (Bronstein Ch. 2–4 / Trager); Mathilda's radical/Goursat
   integrators already cover important algebraic families in the main cascade, so
   this is lower marginal value here.
5. **Phase E — per-call special-function flags** (options to force/suppress the
   Erf / Ei / li / dilog outputs per call).  Low effort; option parsing only.
6. **Phase F — `Simplify` post-processing** for the correct-but-unsimplified outputs
   (I-laden `Tan`/`Tanh` complex-log → real trig; unfactored repeated-pole
   denominators — see §5).  Belongs in `Simplify`, not the integrator.
7. ~~**Wider Ei/li special-function recognizers.**~~ **DONE (2026-07-11).**
   `E^(-x)/x`, `1/Log[2x]`, and `x/Log[x]` now close (§3.8), and the widening
   surfaced + fixed a never-wrong bug in the frac case (non-constant residue).

Recommended stopping point: items 1–3 are diminishing-returns refinements of a
complete engine; items 4–6 are separate subsystems / cosmetic.  Item 7 (the
small, self-contained completeness win) is done.  The integrator is in a
consolidated, correct, well-tested state as of Phase C.

### Phase C — Rational-argument exponentials / non-polynomial `u` ✅ DONE (2026-07-11)
The exponential cases assumed `u` polynomial in `x`.  `rt_solve_rde` now defers,
when `u` or `p` is rational in `x`, to `rt_solve_rde_rational`, which solves
`q = h/Denominator[p]` with `h` a bounded polynomial in `x` (the denominator-theorem
ansatz, as in the tower field RDE).  Closes `Integrate[-E^(1/x)/x^2, x] = E^(1/x)`
and `Integrate[E^x/x - E^x/x^2, x] = E^x/x`.  **Correctness gate (critical):** `p`
and `u'` must be genuine rational functions of `x`; without it a transcendental
coefficient (the `Sin[x]` of a raw `E^x Sin[x]`, reached via `rt_exp_poly_case`)
would let `SolveAlways` certify a spurious `q = 0` — such integrands fall through to
the trig front-end instead.  The full `SPDE` (weak normalization, exact degree
bounds, non-monomial-denominator numerators — Bronstein ch. 6) remains the
refinement; the current numerator ansatz is bounded.

### Phase D — Algebraic extensions
`Sqrt`, `n`-th roots, `RootSum` forms are declined (they are not transcendental
monomials). This is the algebraic Risch case (Bronstein ch. 2–4 / Trager) — a
separate, large subsystem; Mathilda's radical/Goursat integrators already cover
important algebraic families in the main cascade.

### Phase E — Per-call special-function flags
Expose options to force or suppress the special-function outputs (Erf / Ei / li
/ dilog) per call (they are currently default-ON). Low effort; mostly option
parsing on the trailing arguments.

### Phase F — `Simplify` post-processing hooks (the §5 gaps)
Teach `Simplify` the complex-log-to-real-trig collapse and denominator
re-factoring so the correct-but-unsimplified outputs render cleanly. Belongs in
`Simplify`, not the integrator.

---

## 7. Testing & verification

- `tests/test_integrate_risch_transcendental.c` — an **extensive** suite (20 `TEST`
  functions, ~290 assertions) with a dedicated function per case, sub-case and
  sub-sub-case of the transcendental tower, each with many representative
  integrands:
  - **Unit — degree bound:** `rde_var_bound` white-boxes the pure `rt_rde_var_bound`
    arithmetic (exposed in the header): every branch of the Bronstein leading-degree
    bound (deriv-lowering `dfv≥0` / `dfv<0`; deriv-preserving; the clamp at 0) plus the
    cancellation/resonance widening — monotonicity (never lowers the bound) and that
    `m_res` fires only in its exact config (exp `dfv=0`, primitive `dfv=-1`).  Pins the
    "no arbitrary caps" contract directly, without needing an integrand per config.
  - **Base:** `rational_case` (polynomials; distinct / repeated real poles;
    irreducible quadratics; higher-degree denominators; repeated quadratics;
    improper fractions), `rational_agreement` (matches `BronsteinRational`).
  - **Single logarithmic extension:** `logarithmic_case` (θ-powers, poly-in-x
    coefficients, `Log[a x+b]`, `Log` of a quadratic, new-log fold-back),
    `fractional_log_case` (squarefree Rothstein–Trager, single & multi-residue,
    **plus resultant-LRT algebraic residues → `ArcTan[Log[x]]`**),
    `hermite_log_case` (repeated poles, Hermite + log).
  - **Single exponential extension:** `exponential_case` (Laurent, quadratic
    exponents, rational exponent/coefficient), `fractional_exp_case` (incl.
    **resultant-LRT `E^x/(E^(2x)+1) → ArcTan[E^x]`**),
    `hermite_exp_case`, `hyperexponential_case` (coupled Laurent + log, repeated
    / `θ=0` poles).
  - **Front-end / multi-kernel:** `trig_frontend` (powers, products, multiple
    angles, `Tan`/`Tanh`), `multikernel_case` (sum-of-exponentials).
  - **Towers / recursion:** `log_tower_case`, `exp_tower_case` (incl. the
    evaluator-merged monomial, the **multiplicatively commensurate** kernels
    `E^(2 E^x) = (E^(E^x))^2`, and the **exp-top resultant LRT**
    `E^x E^(E^x)/(1+E^(2 E^x)) → ArcTan[E^(E^x)]`), `recursive_tower_case`
    (mixed towers, rational coefficients, proper parts, field RDE, and the
    **cap-free RDE degree bound** `Log[x]^{6,8,12,20} E^(Log[x]^2)`).
  - **Special functions:** `special_functions` (Gaussian→Erf/Erfi, Ei, li,
    dilog).
  - **Guards:** `cascade_default` (Automatic cascade still closes them),
    `method_plumbing` (both surface forms agree), `strict_unevaluated` +
    `strict_misc` (every out-of-scope / non-elementary integrand — Fresnel, Si,
    Ci, a quadratic-denominator Ei form, non-constant residues, non-elementary
    field RDEs, and the `E^(E^x)/(1+E^(E^x))` wrong-answer regression — stays
    unevaluated, never garbage).  The widened Ei/li closes (`E^(-x)/x`,
    `1/Log[2x]`, `x/Log[x]`, `1/Log[c x+d]`) live in `special_functions` (numeric
    diff-back) and double as the regression guard for the frac constant-residue
    bug (§3.8).
- Correctness asserted by diff-back `Simplify[D[Integrate[f,x],x] - f] === 0`
  (`assert_rm_diff_zero`, 172 cases + 4 via the `Method` surface form) or a
  three-point numeric check `assert_rm_num` (special-function / I-laden trig
  outputs whose exact `Simplify` is slower, incl. the new widened Ei/li cases);
  declines asserted by `assert_head_unevaluated`.  Every integrand was empirically
  classified against the built integrator before being pinned, so the suite runs
  green in ~15 s (well under the harness `alarm(60)`).

Build & run (foreground; the test harness sets `alarm(60)`):

```bash
make -j$(nproc)
cd tests && cmake -S . -B build && cmake --build build --target integrate_risch_transcendental_tests
./build/integrate_risch_transcendental_tests            # expect: All Integrate RischTranscendental tests passed!
valgrind --leak-check=full ./build/integrate_risch_transcendental_tests   # no leak stack through the module
```
