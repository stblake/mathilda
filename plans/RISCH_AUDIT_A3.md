# RISCH_AUDIT_A3 — Faithfulness audit vs Bronstein, *Symbolic Integration I* (2nd ed., 2004)

Scope: five modules against the cited algorithm boxes / theorems. Book page N =
PDF page N+17. All references verified against the extracted PDF pages.

**Bottom line: no CRITICAL and no MAJOR inconsistencies found.** All five
modules are faithful ports of the corresponding Bronstein algorithm boxes. The
findings below are OK-BY-DESIGN specializations (deliberate, sound) and a few
MINOR notes. Each deviation was checked for soundness (does it ever produce a
wrong antiderivative / a false classification / a false "reducible").

---

## 1. `intrat.c` — Hermite reduction (§2.2) and LRT log part (§2.3/§2.5)

### 1.1 HermiteReduce — Mack's linear version (book p.44)

Book box:
```
g ← 0 ; D⁻ ← gcd(D, D'); D* ← D/D⁻
while deg(D⁻)>0:
    D⁻² ← gcd(D⁻, D⁻'); D⁻* ← D⁻/D⁻²
    (B,C) ← ExtendedEuclidean(−D* D⁻'/D⁻, D⁻*, A)
    A ← C − B' · D*/D⁻*
    g ← g + B/D⁻
    D⁻ ← D⁻²
return (g, A/D*)
```
Code `intrat_hermite_reduce` (lines 437–532), naming dbar=D⁻, dstar=D*,
dbartwo=D⁻², dbarhat=D⁻*:
- `dbar = monic-gcd(d, d')`, `dstar = exquo(d, dbar)` ✓
- `dbartwo = monic-gcd(dbar, dbar')`, `dbarhat = exquo(dbar, dbartwo)` ✓
- `dtil = −exquo(dstar·dbar', dbar)` = −D* D⁻'/D⁻ ✓ (first ExtEuc arg)
- `(b_coef,c_coef) = ExtendedEuclidean(dtil, dbarhat, a)` ✓
- `a := c − exquo(b'·dstar, dbarhat)` = C − B'·D*/D⁻* ✓
- `g += canonic(b/dbar)`; `dbar := dbartwo`; final `h = canonic(a/dstar)` ✓

`intrat_extended_euclidean` (327–386) verified to return (r, sp) with
r·a + sp·b = c and deg(r) < deg(b) — exactly Bronstein's `ExtendedEuclidean`
3-arg semantics. **Faithful, verdict OK.**

MINOR — `intrat_exquo` (164–173) emits a warning but returns the quotient even
when the remainder is non-zero. Mack's identities guarantee exactness at every
call site, so this never fires on valid input; it is a robustness escape hatch,
not an algorithmic deviation. (Matches the reference `.m`.) Severity MINOR.

### 1.2 IntRationalLogPart — Lazard–Rioboo–Trager (book p.51)

Book box (§2.5):
```
(R,(R0,…,Rk,0)) ← SubResultant_x(D, A − t D')
(Q1,…,Qn) ← SquareFree(R)
for i=1..n with deg_t(Qi)>0:
    if i = deg(D): Si ← D
    else:
        Si ← Rm with deg_x(Rm)=i
        (A1,…,Aq) ← SquareFree(lc_x(Si))
        for j=1..q: Si ← Si / gcd(Aj,Qi)^j       (exact quotient)
return Σ_i Σ_{a|Qi(a)=0} a·log(Si(a,x))
```
Code `intrat_log_part_core` (1685–1911):
- PRS via `SubresultantPolynomialRemainders[d, a − t·D(d), x]` ✓
- `resultant = primitive(Resultant[d, a−t·D(d), x], t)`, then
  `Q = SquareFree(resultant)` ✓
- loop over Q densely indexed by multiplicity i; skip `deg_t(Qi) ≤ 0` ✓
- `i == deg_d` branch ✓ (see note below)
- else branch: `s = find_prs_at_degree(prs, x, i)`; `s = primitive(s, x)`;
  `A = SquareFree(primitive(lc_x(s), t))`; divide out `gcd(Aj,Qi)^(j+1)` via
  `exquo`; `Si = primitive_part_mod(s, Qi, x, t)` ✓
- output builds `RootSum[Function[t,Qi], Function[t, t·Log[Si]]]` ✓

Verified subtleties:
- **Subresultant PRS vs pseudo-remainder chain.** Mathilda's
  `SubresultantPolynomialRemainders` is a *pseudo-remainder* chain, not the
  Lazard-scaled subresultant chain (self-documented at 1722–1734). This is
  sound: (a) the true resultant is computed separately via `Resultant[]` for the
  SquareFree step, and (b) the S_i lookups take `primitive(·, x)`, and
  pp_x(pseudo-Rm) = pp_x(subresultant-Rm) up to a unit because the two chains
  differ only by factors in K[t] (constant in x). The degree sequence of any PRS
  variant is identical, so `find_prs_at_degree` never misses a defective-case
  entry. **No completeness gap.**
- **`i == deg_d` case (1804–1808).** Book says `Si ← D`; code uses
  `primitive_part_mod(d, Qi, x, t)`. Since d is t-free, this reduces to
  `primitive_x(d/lc_x(d))`, i.e. D scaled by a constant. log(Si) differs from
  log(D) by log(constant); derivative unchanged. Bronstein explicitly notes
  (p.51) the monic/pp normalization "is optional". **Sound — OK-BY-DESIGN**
  (documented GCL correction).
- **Squarefree exponent indexing.** `intrat_squarefree_list` returns entries
  densely indexed so position j (0-based) is the multiplicity-(j+1) factor;
  `jpow = j+1`; `gcd(Aj,Qi)^(j+1)` matches book's `gcd(Aj,Qi)^j` with j = the
  multiplicity. ✓
- **Constant field / algebraic residues.** Residues are kept as formal roots of
  Q_i (RootSum) — the minimal-extension "formal sum over roots" form of
  Thm 2.5.1; no splitting field is computed. ✓

**Faithful, verdict OK.** The palindromic-quartic / simple-rootsum radical
expanders (1188–1508) and LogToAtan/LogToReal (Phase 3/4) are Mathilda
extensions beyond the book's IntRationalLogPart (which stops at RootSum); they
are sign-gated and diff-back verified and out of scope for §2.3 faithfulness.

---

## 2. `risch_field.c` — normal/special (Def 3.4.2) and field arithmetic

Def 3.4.2: p **normal** iff gcd(p,Dp)=1; p **special** iff gcd(p,Dp)=p (p|Dp).

- `risch_field_is_normal` (190–196): `deg_t(gcd(p,Dp)) == 0`. ✓ matches gcd=1.
- `risch_field_is_special` (198–206): guards `deg_t(p) ≤ 0 → false`, else
  `deg_t(gcd(p,Dp)) == deg_t(p)`. For positive-degree p, deg(gcd)=deg(p) ⟺
  gcd = c·p ⟺ p|Dp. ✓ matches p|Dp.

Field arithmetic (k = C(other vars) treated as a field):
- `risch_field_gcd_t` (153–171): clears pure-lower-field denominators
  (`Numerator[Together]`), `PolynomialGCD` over Q[all vars], returns 1 when the
  result has t-degree ≤ 0 (shares only a unit of k), else divides by lc_t to make
  monic. Verified correct by Gauss's lemma: for G = gcd over the UFD
  Q[x,t₁,…][t], cont_t(G) | lc_t(G), so G/lc_t(G) = monic·pp_t(G) = the monic
  field gcd. Pure-x/lower-tower common factors (t-degree 0) collapse to the unit
  1. ✓
- `risch_field_divexact_t` (177–184): `Cancel[a/b]` then `PolynomialQ[·,t]`
  gate; NULL when not in k[t]. ✓
- `risch_field_deriv` (111–121): D[p] = Σ dvars[i]·∂p/∂vars[i] — the monomial
  derivation (Lemma 3.2.2). ✓
- `risch_field_polynomial_reduce` (302–330): Bronstein PolynomialReduce (§5.4,
  book p.141). Requires `delta = deg_t(Dt) ≥ 2` (nonlinear, Thm 5.4.1); m =
  dr−δ+1; c = lc_t(r)/(m·λ), λ = lc_t(Dt); q += c·t^m; r −= D[q0]. Matches box
  exactly. ✓

MINOR — units (elements of k, t-degree 0) are classified normal and **not**
special, whereas Bronstein remarks that p ∈ k* is *both* normal and special
((p)=(1)). The binary convention here is the algorithmically useful one and is
consistent with SplitFactor (deg(S)=0 ⇒ p_n=p, p_s=1). No downstream effect.
Severity MINOR (documented in code comments).

**Faithful, verdict OK.**

---

## 3. `risch_canonical.c` — SplitFactor / SplitSquarefreeFactor / CanonicalRepresentation (§3.5)

- `risch_split_factor` (70–102) vs book box p.100:
  `S = divexact(gcd(p,Dp), gcd(p, dp/dt))`; `deg(S) ≤ 0 → (p,1)`; else recurse
  on p/S and return `(qn, Expand(S·qs))`. Exact match to
  `S ← gcd(p,Dp)/gcd(p,dp/dt); if deg(S)=0 return(p,1); (qn,qs)←SplitFactor(p/S);
  return(qn, S·qs)`. Note the second gcd uses the *ordinary* d/dt (`rc_ddt`),
  matching the book. ✓
- `builtin_risch_splitsquarefree` (231–255) vs box p.102: `(p1..pm)=Squarefree(p)`
  (Yun over d/dt, `risch_squarefree_t`), then `Si=gcd(pi,Dpi)`, `Ni=pi/Si`. ✓
  Squarefree factorization uses d/dt; the S_i split uses the monomial D — exactly
  Bronstein's split. ✓
- `risch_canonical_representation` (113–145) vs box p.103:
  `(a,d)=NumDen(f)` (d monic); `(q,r)=PolyDivide(a,d)`; `(dn,ds)=SplitFactor(d)`;
  `(b,c)=Diophantine(dn,ds,r)` with deg_t(b)<deg_t(ds); return
  `(q, b/ds, c/dn)`. Matches `return(q, b/d_s, c/d_n)`. ✓

Properness verified: `risch_field_diophantine_t` (332–350) guarantees
deg_t(b) < deg_t(ds) (reduces b0 mod ds). deg_t(c) < deg_t(dn) then follows
algebraically from c·ds = r − b·dn with deg(r) < deg(ds)+deg(dn) — so f_s = b/ds
and f_n = c/dn are both proper fractions, and f = q + b/ds + c/dn holds. ✓
Diophantine identity b·dn + c·ds = r verified by expansion. ✓

**Faithful, verdict OK.**

---

## 4. `risch_structure.c` — RationalSpan / LogReducible / ExpReducible (Cor 9.3.1)

Cor 9.3.1 (book p.284):
- (i) Da/a is a derivative of an element of K ⟺ ∃ rᵢ ∈ **Q**:
  Σ_{i∈L} rᵢ·Dtᵢ + Σ_{i∈E} rᵢ·(Dtᵢ/tᵢ) = Da/a.
- (ii) Db is the log-derivative of a K-radical ⟺ ∃ rᵢ ∈ **Q**:
  Σ_{i∈L} rᵢ·Dtᵢ + Σ_{i∈E} rᵢ·(Dtᵢ/tᵢ) = Db.

Where L = log monomials (generator Dtᵢ), E = exp monomials (generator Dtᵢ/tᵢ).

Code:
- `rs_decode_tower` (173–229): generator for `"Log"` = Dtᵢ ✓, for `"Exp"` =
  Dtᵢ/tᵢ ✓ — exactly the L/E generators of (9.8)/(9.9). (`"Tan"` = Dtᵢ/(tᵢ²+1)
  extends to the arctangent index set A of the real theorem Cor 9.3.2.)
- `rs_reducible` (233–258): target θ = Cancel(Da/a) for the log test (i) ✓,
  θ = Db for the exp test (ii) ✓, where Da/Db use the full-tower derivation.
- `risch_rational_span` (70–163): sets up θ − Σ cᵢ gᵢ = 0 as a cleared-denominator
  polynomial identity `Numerator[Together[·]] == 0`, solves via `SolveAlways`
  over vars = {x, t₁,…,tₙ}, pins free parameters to 0, then **requires every
  coefficient to be a rational number** (`rs_is_rational`). ✓

Soundness analysis (audit item 4 — "any false 'reducible'?"):
- The rationality gate (149) is the key: a solution with irrational/symbolic/
  non-constant coefficients is rejected → returns False. A genuine rational
  solution to the identity-in-vars certifies (9.8)/(9.9) by Cor 9.3.1, so a
  reported "reducible" is always genuine. **No false positive.**
- θ = 0 short-circuit (73–82) returns all-zero coefficients (0 = Σ 0·gᵢ,
  reducible) — correct. m = 0 with θ ≠ 0 → NULL → False — correct.
- Free-parameter pinning to 0 is valid (free params keep the identity true for
  any value, including 0), so the pinned particular solution is a legitimate
  rational witness. **Sound.**
- Treating x and the tower monomials as independent indeterminates for
  coefficient-matching is exactly Bronstein's "same constant solutions by
  Lemma 7.1.2" reduction — the tower monomials are algebraically independent. ✓

Absence of a rational solution → False is a completeness question, not soundness;
`SolveAlways` is exact (RowReduce/FLINT), so no gap on realizable fields.

**Faithful, verdict OK.**

---

## 5. `risch_hypertangent.c` — IntegrateHypertangentPolynomial (Thm 5.10.2, book p.167)

Book box:
```
(q,r) ← PolynomialReduce(p, D)         (deg r ≤ 1)
α ← Dt/(t²+1)
c ← coefficient(r,t)/(2α)
return (q,c)
```
Code `risch_integrate_hypertangent_poly` (44–70):
- `a = Cancel(Dt/(t²+1))`; gate `FreeQ[a, t]` (this *is* the hypertangent test,
  Dt/(t²+1) ∈ k) ✓
- `(q,r) = PolynomialReduce(p,t,d)` (delta=2 for a hypertangent ⇒ deg_t(r) ≤ 1) ✓
- `r1 = Coefficient(r, t, 1)` = coefficient of t¹ in r ✓
- `c = Cancel(r1 / (2·a))` = coefficient(r,t)/(2α) ✓
- returns {q, c} ✓

**Faithful, verdict OK.** Note (not a defect): the routine returns (q,c) without
computing/checking Dc — this exactly mirrors the book's algorithm box, where the
"Dc≠0 ⇒ not elementary" clause is a *specification comment*, not a computed step;
the non-elementarity certificate is the caller's responsibility (Thm 5.10.2). The
global preconditions √−1 ∉ k and "√−1·Dt/(t²+1) not a log-derivative of a
k(√−1)-radical" are likewise caller obligations, as in the book. The docstring's
"If D[c] is nonzero the integral is not elementary" describes this property
correctly without claiming the function verifies it.

---

## Summary table

| # | Location | Deviation | Severity |
|---|----------|-----------|----------|
| 1.1 | `intrat_exquo` 164–173 | proceeds on non-zero remainder (warns) | MINOR |
| 1.2 | log-part `i==deg_d` 1804–1808 | `primitive_part_mod(d)` instead of literal `D` | OK-BY-DESIGN |
| 1.2 | log-part 1713–1744 | separate `Resultant[]` + pseudo-PRS instead of one `SubResultant` | OK-BY-DESIGN |
| 2 | `is_special`/`is_normal` | units classified normal, not "both" | MINOR |
| 4 | `rs_decode_tower` "Tan" | arctangent generator (extends Cor 9.3.1 → 9.3.2) | OK-BY-DESIGN |
| 5 | hypertangent | Dc≠0 certificate left to caller (as in book box) | OK-BY-DESIGN |

No CRITICAL, no MAJOR. No suggested fixes required for soundness or completeness.
