# RISCH_AUDIT_A2 — Faithfulness audit of `integrate_risch_transcendental.c`

Reference: Bronstein, *Symbolic Integration I: Transcendental Functions*, 2nd ed.
(2004), §§5.3–5.9 (book pp. 139–164). Target file:
`src/calculus/integrate_risch_transcendental.c`; the Lazard–Rioboo–Trager
resultant core lives in `src/calculus/intrat.c`
(`intrat_log_part_core`, exported as `Integrate`TranscendentalLogPart`).

## Executive summary

- **No CRITICAL (unsound / wrong-answer) findings.** Every path that produces a
  closed form is gated by an *exact* certificate before it can ship:
  - single-kernel LRT / Hermite: `SolveAlways[Q'−f==0,{t,x}]` (a polynomial
    identity) plus `rt_verify_antideriv` (`Simplify[D[Q]−f]===0`), lines 1745, 1954;
  - tower path: an **exact rational-identity** check
    `Together[D_tower[Q] − F] === 0` in the tower variables (lines 4342–4360),
    with a defensive `Simplify` diff-back fallback (4371).
  A spurious certificate cannot be produced by a sound `SolveAlways` / `Together`,
  so declines are the only failure mode.
- The code is a **faithful transcription** of Bronstein for: the LRT
  `ResidueReduce` resultant/subresultant-PRS core (`intrat_log_part_core`), the
  primitive-polynomial recursion (`IntegratePrimitivePolynomial`, §5.8), the
  hyperexponential Laurent recursion (`IntegrateHyperexponentialPolynomial`,
  §5.9), and the Risch DE (`rt_solve_rde` → `rde_base`, Chapter 6 SPDE).
- It **substitutes an undetermined-coefficient ansatz** (correct-by-construction,
  never a literal transcription) for `HermiteReduce` (§5.3) and for the tower
  proper/coupled parts. This is where the MAJOR completeness gaps live.

Severity legend: CRITICAL = unsound; MAJOR = misses valid elementary cases;
MINOR = cosmetic; OK-BY-DESIGN = deliberate, correct specialization.

---

## 1. Hermite reduction (§5.3, `HermiteReduce`, p.139)

### Book algorithm (quadratic version, p.139)
```
(fp,fs,fn) ← CanonicalRepresentation(f,D)
(a,d) ← (numerator(fn), denominator(fn))          (d monic)
(d1,...,dm) ← SquareFree(d)
g ← 0
for i ← 2 to m such that deg(di) > 0 do
    v ← di ;  u ← d/v^i
    for j ← i-1 to 1 step -1 do
        (b,c) ← ExtendedEuclidean(u Dv, v, -a/j)
        g ← g + b/v^j
        a ← -jc - u Db
    d ← uv
(q,r) ← PolyDivide(a, uv)
return(g, r/(uv) + fp + fs)
```
The reduced result `H/D_red + Σ c_j Log(g_j)` is built by the squarefree-part
outer loop and an `ExtendedEuclidean` solve **per multiplicity level j**, and the
rational-part coefficients live in `k = C(x)` (arbitrary rational functions).

### Code (`rt_hermite_try`, lines 1776–1974; tower analogue `rt_field_ratint`, 3441–3598)
The code does **not** implement this loop. It:
1. forms `Hden = PolynomialGCD(den, ∂den/∂t)` (the repeated part `∏ D_m^{m-1}`)
   and `rad = den/Hden` (squarefree radical) — lines 1832–1839;
2. reads the distinct t-dependent factors `g_j` of `rad` — 1845–1864;
3. builds a **generic ansatz** `Q = H(t)/Hden + Σ c_j Log(g_j)` with
   `deg_t H < deg(Hden)` and **bounded-degree polynomial-in-x coefficients**
   (lines 1872–1897), then solves all unknowns with
   `SolveAlways[Q'−F==0,{t,x}]` (1916) and back-substitutes.

#### Finding 1.1 — MAJOR (completeness): H-coefficients restricted to polynomials of a heuristic x-degree
- **Book:** the Hermite numerator's coefficients are elements of `k = C(x)`, i.e.
  arbitrary **rational** functions of `x`, and the reduction is *complete*.
- **Code:** the ansatz coefficient of each `x^k t^p` term is an *unknown scalar*
  and `k` ranges only `0 ≤ k ≤ Nx` with
  `Nx = max(deg_x num, deg_x den) + 2` (lines 1866–1868; tower: `bd[j] =
  max(deg_v num, deg_v den)+1`, 3481–3486). So a Hermite numerator whose x-
  coefficients are genuinely rational, or polynomial of degree `> Nx`, is not
  representable and the ansatz **declines**.
- **Severity:** MAJOR — misses valid elementary integrals of repeated-pole
  integrands. It is **not** a soundness bug (`SolveAlways` + diff-back gate).
  The source comment (3437–3439) already acknowledges "a genuinely rational
  Hermite numerator coefficient … remain a later refinement."
- **Suggested fix:** implement the literal `HermiteReduce` (squarefree
  factorization of the t-denominator + the `ExtendedEuclidean(u Dv, v, -a/j)`
  level loop) over `k = C(x)`, which yields the exact rational-coefficient
  numerator and removes the x-degree heuristic entirely. The monomial-extension
  version is valid (Bronstein Exercise 5.2); `D` = the tower derivation
  `rt_tower_deriv`.

#### Finding 1.2 — OK-BY-DESIGN: multiplicity handling, no off-by-one
- The ansatz denominator is `Hden = gcd(den, ∂den/∂t)` and `H` is generated with
  `p` in `0..dH-1` where `dH = deg_t(Hden)` (loop `for (long p = 0; p < dH; …)`,
  1872). This reproduces exactly `deg_t(H) < deg(Hden)`, the structural output of
  Bronstein's loop. **No off-by-one** in the multiplicity range: the book's
  inner loop `j = i-1 … 1` and the collapsed `Hden = ∏ D_m^{m-1}` give the same
  `deg(Hden) = Σ (m-1) deg(D_m)`, matched here.
- The `∂/∂t`-only gcd (rather than the full derivation `D`) is used **only to
  size the ansatz**; correctness is enforced by the full `D_tower` in the
  `SolveAlways` residual, so any mismatch just declines. Correct.

---

## 2. Residue criterion / Rothstein–Trager / LRT (§5.6, pp. 147–154)

Two implementations coexist:
- **`rt_frac_try` (single kernel, lines 1463–1590+)** — the *rational-residue*
  ansatz: solves `num = Σ c_i D(g_i)(d/g_i)` for **constant** `c_i` via
  `SolveAlways`, with an explicit **constant-residue check** (1572–1588): every
  solved residue must be free of `x` and `t`, else decline. This directly
  enforces Lemma 5.6.2 (residues ∈ Const(k)) and the memory note
  `project_solvealways_nonconstant_residue`. Correct.
- **`rt_frac_lrt` / `rt_field_lrt_logpart` → `intrat_log_part_core`** — the
  genuine **Lazard–Rioboo–Trager** resultant reduction for *algebraic* residues.

### Book `ResidueReduce` (LRT version, p.153)
```
d ← denominator(f) ;  (p,a) ← PolyDivide(numerator(f), d)
z ← new indeterminate
if deg(Dd) ≤ deg(d) then (r,(R0,…,Rq)) ← SubResultant_x(d, a - zDd)
                    else (r,(R0,…,Rq)) ← SubResultant_x(a - zDd, d)
((n1,…,nn),(s1,…,sn)) ← SplitSquarefreeFactor(r, κ_D)
for i ← 1 to n such that deg(si) > 0 do
    if i = deg(d) then Si ← d
    else  Si ← Rm where deg_t(Rm)=i ; (A1,…,As)←SquareFree(lc_t(Si)) ; …
if ∏ ni_i ∈ k then b ← 1 else b ← 0
return(Σ_i Σ_{α|si(α)=0} α log(Si(α,t)),  b)
```
Residue criterion (Thm 5.6.1): `r = r_s r_n`; **`r_n ∈ k` ⇔ elementary**;
`r_n ∉ k` ⇒ `f − Dg ∉ k[t]` and, via Liouville, no elementary integral.

### Code (`intrat_log_part_core`, intrat.c 1685–1911)
- **Resultant + subresultant PRS:** computes both
  `SubresultantPolynomialRemainders(d, a−t·D(d), x)` (1713–1716) for the `S_i`
  lookups and `primitive_t(Resultant(d, a−t·Dd, x))` (1735–1743) for the residue
  polynomial. The `i == deg_d` special case uses `d` itself (1804–1808), and the
  `A_j`-squarefree division of leading coefficients (1815–1845) — a **faithful**
  transcription of the p.153 inner loop.
- **Constant-residue decision (`b`)** is implemented as the `xgate` gate
  (1746–1776): the residue polynomial, after stripping x-only content, must be
  `FreeQ` of `x` **and** every lower tower variable, i.e. *all* roots are
  constants of the derivation. This is exactly `r_n ∈ k`.

#### Finding 2.1 — OK-BY-DESIGN: `SubResultant` argument order never swapped
- **Book:** swaps to `SubResultant(a−zDd, d)` when `deg(Dd) > deg(d)`.
- **Code:** always `SubresultantPolynomialRemainders(d, a−t·Dd, x)` (1713),
  never swaps.
- **Assessment:** correct here. The PRS main variable is the **monomial** `t`.
  For a log kernel `D` *lowers* deg_t; for an exp kernel `D(t^n)=n w' t^n`
  *preserves* deg_t. Hence `deg_t(Dd) ≤ deg_t(d)` **always**, so the swap branch
  is unreachable. No fix needed. (Would become a bug only if a non-monomial main
  variable were ever passed — it never is.)

#### Finding 2.2 — MAJOR (completeness): `b = 0` partial elementary part is discarded
- **Book:** even when `r_n ∉ k` (b=0), `ResidueReduce` **returns `g`** — the
  logarithms coming from the *constant* roots `r_s` — and only reports the `r_n`
  part as non-elementary. Bronstein's own Example 5.6.2
  (`∫(2log²x−logx−x²)/(log³x−x²logx) dx`) returns
  `½ log((logx+x)/(logx−x)) + Li(x)`, where the `½ log(…)` comes from the
  constant residues `z=±½` even though `z=x` is a non-constant residue.
- **Code:** the `xgate` gate (intrat.c 1759–1775) requires the **entire**
  resultant free of the lower-field variables; if *any* residue is non-constant,
  `free_all=false` and it returns `NULL` — the whole log part, **including the
  elementary constant-residue logs**, is dropped.
- **Severity:** MAJOR completeness gap (not soundness — the overall integrand in
  such a case is itself non-elementary, so Mathilda correctly returns it
  unevaluated; it simply cannot surface the partial closed form). Consistent with
  Mathilda's whole-or-nothing integral policy, but it means Bronstein's headline
  §5.6 example cannot be reproduced even partially.
- **Suggested fix:** split the resultant by `κ_D` into `r_s` (constant roots) and
  `r_n` (non-constant), build `g` from `r_s` only, and return `g` with an
  unevaluated remainder for the `r_n` part rather than declining outright.

#### Finding 2.3 — MAJOR / OK-BY-DESIGN: decision procedure reduced to a semi-decision
- **Book:** `ResidueReduce` returns a Boolean `b` that, combined with Liouville
  (Thm 5.5.1) and the primitive/hyperexponential polynomial cases, makes the
  whole chapter a **decision procedure** — it *proves* non-elementarity.
- **Code:** `b=0` is only ever consumed as "decline (return NULL/unevaluated)".
  Mathilda never emits a positive "provably non-elementary" verdict.
- **Severity:** structural MAJOR, but reasonable OK-BY-DESIGN for a CAS whose
  contract is "closed form or unevaluated." Worth documenting that the recursive
  Risch here is a *semi*-decision procedure, not the full Bronstein decision
  procedure.

---

## 3. Primitive polynomial case (§5.8, `IntegratePrimitivePolynomial`, p.158)

### Book
```
if p ∈ k then return(0,1)
a ← lc(p)
(b,c) ← LimitedIntegrate(a, Dt, D)         (* a = Db + cDt, c constant *)
if (b,c) = "no solution" then return(0,0)
m ← deg(p)
q0 ← c t^{m+1}/(m+1) + b t^m
(q,β) ← IntegratePrimitivePolynomial(p − Dq0, D)
return(q + q0, β)
```
i.e. `Q=Σq_i t^i` with `D Q = p` solved by `Dq_i + (i+1)q_{i+1}Dt = p_i`, top
coefficient first; `LimitedIntegrate` returns `∫a = b + c t`, the `c t^{m+1}/(m+1)`
term feeding the next-higher coefficient.

### Code (`rt_log_poly_case`, 622–715; tower `rt_int_primitive_poly`, 3867–3924)
Single-pass unrolling of the recursion (i = m … 0):
```
r_i = p_i − (i+1) q[i+1] Dt                       (664–671 / 3884–3891)
(s,c) = LimitedIntegrate(r_i)                     (673 / 3893) → s + c·t
q[i+1] += c/(i+1)                                 (677–685 / 3897–3903)
q[i]   = s
```
- **Correct mapping to the book.** `q[i]=s=b` and `q[i+1]+=c/(i+1)` reproduce
  exactly `q0 = c t^{m+1}/(m+1) + b t^m`. Subtracting `Dq0` from `p` modifies
  `p_{m-1} ← p_{m-1} − m b Dt`, which the loop reproduces as
  `r_{m-1}=p_{m-1} − m·q_m·Dt`. **Verified equivalent.** No off-by-one; `q`
  array `0..m+1` correctly holds the degree-raised top log. `m<1` decline
  (643) is fine — the `p∈k` (m=0) base is handled by the caller integrating in
  K directly.

#### Finding 3.1 — OK-BY-DESIGN: `LimitedIntegrate` realized by "integrate-then-shape"
- **Book:** `LimitedIntegrate(a, Dt)` is a dedicated Chapter-7 primitive.
- **Code:** `rt_limited_integrate` (581–608) / `rt_limited_field_integrate`
  (3932–3958) instead fully integrate `r_i` in the lower field and *check* the
  result has the form `s + c·t` with `c` a constant of the derivation (free of
  `x` and all tower vars, 3946–3948), else decline. This is equivalent to
  `LimitedIntegrate` whenever the lower-field integrator is complete, because
  `a = Db + cDt ⟺ ∫a = b + ct`. The "recursive InFieldDerivative-style check for
  the top coefficient" the audit asked about is present and correct: the
  constant-`c` gate is the elementarity certificate for the top term.
- **Severity:** OK-BY-DESIGN. A genuine `LimitedIntegrate` would be marginally
  more general only in pathological lower fields; not an issue for `k=C(x)` towers.

---

## 4. Hyperexponential polynomial case (§5.9, `IntegrateHyperexponentialPolynomial`, p.161)

### Book
```
q ← 0, β ← 1
for i ← ν_t(p) to −ν_∞(p) such that i ≠ 0 do
    a ← coefficient(p, t^i)
    v ← RischDE(i Dt/t, a)          (* a = Dv + i(Dt/t) v *)
    if v = "no solution" then β ← 0 else q ← q + v t^i
return(q, β)
```
`p ∈ k[t,t^{-1}]`; the Laurent range is `[ν_t(p), −ν_∞(p)]`; each nonzero power
solves a **Risch DE**; `i=0` is deferred to the caller (remains `∈ k`).

### Code (`rt_exp_poly_case`, 1364–1449; tower `rt_int_hyperexp_poly`, 3966–4022)
- **Kernelization** (`rt_exp_kernelize`, 1320–1362): collapses multiplicatively
  commensurate exponents onto one primitive `t=E^u` (`E^(2u)→t²`, `E^{x/2},E^{x/3}
  →E^{x/6}`), declining genuinely independent kernels — sound and beyond the raw
  book statement.
- **Monomial-denominator gate:** requires `den = c·t^M` (1382–1386 / 3980), i.e.
  a genuine Laurent polynomial; a real proper fraction is declined and left to
  the fractional/coupled case. Matches §5.9's `p ∈ k[t,t^{-1}]` precondition.
- **Laurent range:** `i = j − M`, `j = 0..deg_t(num)` ⇒
  `i ∈ [−M, deg_t(num) − M]` = `[ν_t(p), −ν_∞(p)]`. **Exact, no cap, no
  off-by-one** (M = mult of t at 0 in den = `−ν_t`; top = deg num − deg den =
  `−ν_∞`). Faithful.
- **Per-power solve:** `i≠0 → rt_solve_rde(p_i, i, u, x)` solving
  `q_i' + i u' q_i = p_i` with `f = i·u' = i·Dt/t` (since `Dt = u' t`). Exactly
  `RischDE(i Dt/t, a)`. `rt_solve_rde → rde_base` is the genuine Chapter-6 SPDE
  rational solver (`WeakNormalizer → RdeNormalDenominator → RdeBoundDegree →
  SPDE → PolyRischDENoCancel`), degree bound `rt_rde_var_bound` derived, no
  arbitrary cap. Faithful.

#### Finding 4.1 — OK-BY-DESIGN (deviation, more complete): `i=0` term is integrated here
- **Book:** the loop condition `i ≠ 0` **skips** the `t^0` coefficient, deferring
  the remaining `p_0 ∈ k` to the caller.
- **Code:** integrates it in-line — `i==0 → rt_integrate_in_K_with_logs` (1408) /
  `rt_field_integrate(…,L−1,…)` (3996). This is **more** complete than the book
  fragment (the total antiderivative is finished, not left as `g + (p_0∈k)`), and
  remains correct because `p_0` is free of `t`. Note as a deliberate deviation,
  not a defect.

#### Finding 4.2 — OK-BY-DESIGN: coupled proper part uses ansatz, not the split algorithm
- The exp Hermite/log parts do not separate (`D Log(g)=D_tower(g)/g` mixes into
  `t^0`); `rt_field_hyperexp_coupled` (3614–3797) solves a **unified** Laurent +
  Hermite + log ansatz via `SolveAlways`. Same completeness caveat as Finding 1.1
  (bounded polynomial lower-field coefficients, `ihi=dnum−dden`, `ilo=−a`). The
  pure-Laurent recursion (`rt_int_hyperexp_poly`) is tried first (faithful), and
  the coupled ansatz is a correct-by-construction fallback; the pure-resultant LRT
  (`rt_field_lrt_logpart`) closes algebraic residues when `a==0` (3788). Sound;
  MAJOR completeness ceiling for non-polynomial coefficients as in §1.

---

## 5. Cross-cutting notes

- **RDE degree bound (`rt_rde_var_bound`, 770–787; `rt_resonance_int`, 798–809):**
  a genuine Bronstein `RdeBoundDegree` derivation (leading-degree balance +
  exponential/primitive resonance widening), monotone, no arbitrary cap. Matches
  the memory note `feedback_no_arbitrary_caps_decision_procedures`. Sound —
  a too-small bound only declines; a too-large one wastes ansatz terms.
- **Verification net:** single-kernel LRT/Hermite → `rt_verify_antideriv`
  (`Simplify===0`, 2372–2382, uses `Simplify` not `PossibleZeroQ`, per memory
  note `project_possiblezeroq_decay_false_positive`); tower → exact
  `Together[D_tower[Q]−F]===0` (4342–4360). Both are sound certificates.

## 6. Verdict table

| # | Location | Book ref | Severity | One-line |
|---|----------|----------|----------|----------|
| 1.1 | `rt_hermite_try` 1866–1897; `rt_field_ratint` 3481–3534 | §5.3 p.139 | MAJOR | Hermite done by bounded-**polynomial**-coeff ansatz, not the exact `k`-rational `HermiteReduce`; declines rational/high-degree numerators |
| 1.2 | `rt_hermite_try` 1872, 1835–1839 | §5.3 | OK-BY-DESIGN | Multiplicity range `deg_t H < deg(Hden)` correct; no off-by-one |
| 2.1 | `intrat_log_part_core` intrat.c 1713 | §5.6 p.153 | OK-BY-DESIGN | No `SubResultant` swap needed — `t` monomial ⇒ `deg_t(Dd) ≤ deg_t(d)` always |
| 2.2 | `intrat_log_part_core` intrat.c 1759–1775 | Thm 5.6.1 | MAJOR | `b=0` case discards the elementary constant-residue (`r_s`) logs instead of returning them |
| 2.3 | LRT `b` consumers | §5.6 | MAJOR / OK-BY-DESIGN | Decision procedure reduced to semi-decision; never proves non-elementarity |
| 3.1 | `rt_limited_(field_)integrate` 581–608 / 3932–3958 | §5.8 p.158 | OK-BY-DESIGN | `LimitedIntegrate` realized by integrate-then-shape + constant-`c` gate; equivalent |
| — | `rt_log_poly_case`/`rt_int_primitive_poly` recursion | §5.8 | OK (faithful) | `q_i+=c/(i+1)`, `r_i=p_i−(i+1)q_{i+1}Dt` matches `q0`; verified equivalent |
| 4.1 | `rt_exp_poly_case` 1408 / `rt_int_hyperexp_poly` 3996 | §5.9 p.161 | OK-BY-DESIGN | Integrates `i=0` term (book defers it); more complete, still correct |
| 4.2 | `rt_field_hyperexp_coupled` 3614–3797 | §5.9 | OK-BY-DESIGN | Coupled part = unified ansatz; same poly-coeff ceiling as 1.1 |
| — | Laurent range `[−M, degnum−M]` 1405 / 3663–3665 | §5.9 | OK (faithful) | Equals `[ν_t(p),−ν_∞(p)]`; no off-by-one |
| — | `rt_solve_rde`→`rde_base` | Chap. 6 | OK (faithful) | Genuine SPDE rational RDE, derived degree bound |

**Bottom line:** the file is *sound* (no CRITICAL). The two MAJOR items are both
completeness gaps flowing from one design choice — an undetermined-coefficient
ansatz with polynomial (not rational) lower-field coefficients in place of the
literal `HermiteReduce`, and the whole-or-nothing treatment of `ResidueReduce`'s
`b=0` split. Both are safe (gated by exact certificates) but keep Mathilda from
reproducing Bronstein's own §5.3/§5.6 examples in full.
