# Risch Audit A1 — Chapter 6 RDE block faithfulness to Bronstein (2nd ed., 2004)

Target: `src/calculus/integrate_risch_transcendental.c`, lines 828–1272
Reference: Bronstein, *Symbolic Integration I*, §6.1 (pp.183,185), §6.3 (p.199),
§6.4 (p.203), §6.5 (pp.208–209).

**Bottom line: the implementation is a faithful, correct base-field (C(x),
D=d/dx) specialization of Bronstein's RischDE pipeline. No CRITICAL and no
MAJOR faithfulness defects were found.** Every deviation from the printed
algorithm is either a behavior-equivalent simplification or a deliberate,
provably-safe base-field specialization. Details per algorithm below.

---

## 1. WeakNormalizer — `rde_weak_normalizer` (lines 1050–1128)

Book (p.183):
```
(d_n,d_s) ← SplitFactor(denominator(f),D)
g ← gcd(d_n, d_n')
d* ← d_n/g
d1 ← d*/gcd(d*,g)
(a,b) ← ExtendedEuclidean(denominator(f)/d1, d1, numerator(f))
r ← resultant_t(a - z·Dd1, d1)
(n1,…,ns) ← positive integer roots of r
return ∏ gcd(a - n_i·Dd1, d1)^{n_i}
```

Code mapping (all correct):
- SplitFactor ≡ identity in the base field (`d_n = den`, `d_s = 1`) — correct:
  over C(x) every irreducible p has `Dp = p'` of lower degree, so `p ∤ p'`,
  hence every irreducible is normal and S^irr = ∅.
- `g=gcd(den,den')`, `dstar=den/g`, `d1=dstar/gcd(dstar,g)` — exact match.
- `A=den/d1`; `PolynomialExtendedGCD(A,d1)` gives cofactor `s` of A; then
  `a = s·num mod d1` (deg<deg d1) — this is exactly `ExtendedEuclidean(den/d1,
  d1, num)` reduced (valid because `gg==1`).
- `azd = a - z·Dd1`, `res = Resultant(azd,d1,x)`, positive-integer roots via
  `Solve`, `w *= gcd(a - nv·Dd1, d1)^nv` — exact match of the resultant/residue
  step and the product.

Findings:
- **FINDING 1 (OK-BY-DESIGN).** Lines 1079–1122 gate the residue computation on
  the extended-gcd returning `gg == 1`, else leave `w = 1`. This is *not* a
  completeness gap: `gcd(den/d1, d1) = 1` holds for **every** input, not just
  "the common case". Proof: writing `den = ∏ p_i^{e_i}`, one gets
  `d1 = ∏_{e_i=1} p_i` and `den/d1 = ∏_{e_i≥2} p_i^{e_i}`, which are supported on
  disjoint irreducibles, so their gcd is 1 unconditionally. The gate therefore
  never spuriously declines. The docstring wording ("Only the common gg==1
  case") understates this — it is in fact the *only* case. Cosmetic only.

Verdict: faithful.

---

## 2. RdeNormalDenominator — inlined in `rde_base` step 2 (lines 1173–1207)

Book (p.185):
```
(d_n,d_s) ← SplitFactor(denominator(f),D)
(e_n,e_s) ← SplitFactor(denominator(g),D)
p ← gcd(d_n,e_n)
h ← gcd(e_n, e_n') / gcd(p, p')
if e_n ∤ d_n·h^2 then return "no solution"
return (d_n·h,  d_n·h·f − d_n·Dh,  d_n·h^2·g,  h)
```

Code mapping (all correct, base field ⇒ `d_n=den(f)`, `e_n=den(g)`):
- `dn=den(fbar)`, `en=den(gbar)`, `p=gcd(dn,en)` — match.
- `h = gcd(en,en') / gcd(p,p')` (lines 1178–1182) — exact match of the `h`
  formula.
- Guard `en | dn·h^2` (lines 1185–1188) — exact match of Corollary 6.1.1(ii);
  declining (result stays NULL) is the correct "no solution".
- `a = dn·h` (line 1199) — match.
- `b = dn·h·f − dn·Dh` (lines 1200–1203): `hnf = h·num_f = h·(dn·fbar) = dn·h·f`,
  `dndh = dn·Dh`, `bb = hnf − dndh` — exact match, correct sign.
- `c = dn·h^2·g` (lines 1204–1205): `(dn·h^2/en)·num_g = (dn·h^2/en)·(en·gbar)
  = dn·h^2·g` — exact match.

Note the comment at 1192–1196 justifies computing a,b,c as *exact polynomials*
(via `Numerator`/exact `PolynomialQuotient`) instead of `Cancel` of a rational,
to avoid an unreduced factored/expanded mismatch. This is a robustness
improvement, not a deviation; the results equal the book's a,b,c.

Findings:
- **FINDING 2 (MINOR / defensive).** Line 1208 gates the rest on
  `rt_is_poly(aa,x) && rt_is_poly(bb,x) && rt_is_poly(cc,x)`. By Corollary 6.1.1
  a,b,c are *guaranteed* polynomials, so this guard is defensive; if it ever
  triggered on a valid input (Cancel/Together failing to reduce) the code would
  silently decline a solvable problem. In practice the exact-polynomial
  construction above makes it unreachable. Keep as a safety net.

Verdict: faithful.

---

## 3. RdeBoundDegreeBase — inlined in `rde_base` step 3 (lines 1209–1226)

Book (p.199):
```
d_a←deg(a), d_b←deg(b), d_c←deg(c)
n ← max(0, d_c − max(d_b, d_a−1))
if d_b = d_a−1 then
    m ← −lc(b)/lc(a)
    if m ∈ Z then n ← max(0, m, d_c − d_b)
return n
```

Code mapping:
- `mx = max(db, da−1)`, `n = max(0, dc − mx)` (lines 1211–1212) — exact match.
- When `db==da−1`: `mx=db`, so `n=max(0,dc−db)`; block adds `m=−lc(b)/lc(a)` and
  `cand=max(0,dc−db)`, taking `n=max(n,m,cand)=max(0,m,dc−db)` — exact match.
- Non-integer `m` (Rational/rational-fn after `Cancel`) is not `EXPR_INTEGER`,
  so skipped — matches `m ∈ Z` test. Negative integer `m` contributes nothing
  via `if (mv>n)` — matches `max(0,m,…)`.

Findings:
- **FINDING 3 (OK-BY-DESIGN).** Line 1213 adds `&& da >= 1` to the book's plain
  `d_b = d_a − 1` test. The only case this excludes is `da=0, db=−1` (i.e.
  `a∈C*, b=0`, equation `Dq=c/a`). There the book's special branch yields
  `m=0∈Z ⇒ n=max(0,0,dc+1)=dc+1`, and the *base* formula already gives
  `n=max(0,dc−max(−1,−1))=dc+1`. Identical result, so the gate is harmless.

Verdict: faithful.

---

## 4. SPDE — `rde_spde` (lines 904–979)

Book (p.203):
```
if n<0 then { if c=0 return (0,0,0,0,0) else return "no solution" }
g ← gcd(a,b)
if g ∤ c then return "no solution"
a←a/g, b←b/g, c←c/g
if deg(a)=0 then return (b/a, c/a, n, 1, 0)
(r,z) ← ExtendedEuclidean(b,a,c)          (* b·r + a·z = c, deg(r)<deg(a) *)
u ← SPDE(a, b+Da, z−Dr, D, n−deg(a))
if u="no solution" return "no solution"
(b̄,c̄,m,α,β) ← u
return (b̄, c̄, m, a·α, a·β + r)
```

Code mapping (all correct):
- `n<0` base (lines 911–918): `c=0` ⇒ success tuple `(0,0,0,0,0)`; else decline.
  This is exactly right: at the recursive entry the passed `c` is `z−Dr` from
  the parent, so `c=0` is precisely the condition that `h=0` solves (6.16). ✓
- `g=gcd(a,b)`, `g∤c` decline via remainder (lines 919–923) — match.
- `a1=a/g, b1=b/g, c1=c/g` — match.
- `deg(a1)=0` base (lines 929–937): returns `(b1/a1, c1/a1, n, 1, 0)`; dividing
  by the constant `a1` correctly makes the residual equation monic
  (`DH + b̄H = c̄`), which is exactly what the downstream NoCancel solvers
  expect (leading Dq-coefficient 1). ✓
- ExtendedEuclidean (lines 940–951): `PolynomialExtendedGCD(b1,a1)` returns
  cofactor `s` of `b1` (Mathematica order: first cofactor ↔ first poly);
  `r = s·c1 mod a1` gives `b1·r ≡ c1 (mod a1)` with `deg(r)<deg(a1)`. Since
  `gcd(a1,b1)=1` (we divided out `g`), the returned gcd is the constant 1, so
  the scaling is exact. ✓ Matches `(r,z)=ExtendedEuclidean(b,a,c)`.
- `z=(c1−b1·r)/a1` exact quotient (lines 952–957) — match; exactness is
  guaranteed since `c1−b1·r ≡ 0 (mod a1)`.
- Recursion `SPDE(a1, b1+Da1, z−Dr, n−da)` (lines 958–966) — exact match,
  including the `n − deg(a/g)` degree reduction.
- Reassembly (lines 971–977): `α = a1·α_sub`, `β = a1·β_sub + r`, adopt inner
  `(b̄,c̄,m)` — exact match of `(b̄, c̄, m, a·α, a·β+r)`.

Findings: none. This is a textbook-exact transcription of Rothstein's SPDE,
including the tricky `(r,z)` computation, the `n−deg(a)` reduction, the base
cases, and the `α,β` reassembly. The termination argument (deg reduced by
`deg(a)≥1` per level) holds.

Verdict: faithful.

---

## 5. PolyRischDENoCancel1 (b≠0) — `rde_polyrischde_nocancel1` (lines 981–1015)

Book (p.208):
```
q←0
while c≠0 do
    m ← deg(c)−deg(b)
    if n<0 or m<0 or m>n then return "no solution"
    p ← (lc(c)/lc(b))·t^m
    q ← q+p
    n ← m−1
    c ← c − Dp − bp
return q
```

Code mapping: `m=dc−db`, the `n<0||m<0||m>n` decline, `coeff=lc(c)/lc(b)`,
`p=coeff·x^m`, `q+=p`, `c←c−Dp−bp`, `n=m−1` — exact match. (The `n=m−1` update
happens after the `c` update in code, but since it does not depend on `c` this
is irrelevant.)

Dispatch: `rde_base` routes here when `sp.b ≠ 0` (line 1235). Correct: for the
base field `D=d/dx`, NoCancel1's applicability precondition "`b≠0` and
(`D=d/dt` or …)" is satisfied for *any* `b≠0`, and Lemma 6.5.1(i) guarantees no
leading-term cancellation.

Verdict: faithful.

---

## 6. PolyRischDENoCancel2 (b=0) — `rde_polyrischde_integrate` (lines 1017–1048)

Book (p.209), specialized to base field (`b=0`, `δ(t)=0`, `λ(t)=1`):
```
q←0
while c≠0 do
    if n=0 then m←0 else m←deg(c)−δ+1        (= deg(c)+1)
    if n<0 or m<0 or m>n then return "no solution"
    if m>0 then p ← (lc(c)/(m·λ))·t^m
    else { if deg(b)≠deg(c) return "no solution"; … }   (b=0 ⇒ always declines)
    q←q+p ; n←m−1 ; c←c−Dp−bp
return q
```

Code mapping: `m=dc+1`, decline on `n<0||m<0||m>n`, `coeff=lc(c)/m`,
`p=coeff·x^m`, `q+=p`, `c←c−Dp` (b=0), `n=m−1` — matches the `m>0` branch
exactly (antidifferentiation).

Findings:
- **FINDING 4 (OK-BY-DESIGN / behavior-equivalent).** The code does not
  replicate the book's `if n=0 then m←0` branch. In the base field this branch
  is unreachable-as-a-success: with `b=0` and `c≠0`, `m←0` falls into the `else`
  clause where `deg(b)≠deg(c)` (since `deg(0) ≠ deg(c≥0)`) forces "no solution".
  The code instead reaches the same decline through `m=dc+1 ≥ 1 > 0 = n ⇒
  m>n ⇒ NULL`. Both decline; when `c=0` both return `q=0`. Equivalent.

Dispatch: `rde_base` routes here when `sp.b = 0` (line 1233). Correct — in the
base field the only `b=0` sub-case is NoCancel2, and there is no
δ(t)≥2 cancellation case (NoCancel3) to worry about.

Verdict: faithful.

---

## 7. Driver — `rde_base` (lines 1130–1255) and the substitution chain

- Weak normalization (lines 1156–1171): `fbar = f − Dw/w`, `gbar = w·g`, and the
  final `y = q/(h·w)`. This matches Bronstein's p.183 substitution `z = w·y`
  (⇒ `y = z/w`) composed with Corollary 6.1.1's `q = y·h` (⇒ `y = q/h`), giving
  `y = q/(h·w)`. Correct (lines 1242–1243).
- Zero-`g` short-circuit (line 1137): `Dy+fy=0` ⇒ `y=0`. Correct.
- Pipeline order WeakNormalizer → RdeNormalDenominator → RdeBoundDegreeBase →
  SPDE → PolyRischDENoCancel1/2 matches Bronstein exactly.

**Base-field completeness (task question 5).** Dropping the special part and
using SplitFactor=identity is *correct and complete* for `D=d/dx` over `C(x)`:
S^irr = ∅ (every irreducible is normal), so there is no RdeSpecialDenom step to
perform. δ(t)=0 means the only non-cancellation branches are NoCancel1 (`b≠0`)
and NoCancel2 (`b=0`); the cancellation case NoCancel3 requires δ(t)≥2 and never
arises. The two-way `sp.b == 0 ? integrate : nocancel1` dispatch therefore
covers every base-field sub-equation. No valid base-field case is silently
dropped.

---

## 8. Cross-file note (outside the audited line range)

- **FINDING 5 (verify downstream, not a defect in these lines).** Lines
  1139–1154 short-circuit through `flint_rde_base_solve_fg` when `x` is a
  symbol and the FLINT verdict `nr>=0` is treated as authoritative
  (`nr==0` = decline without falling through to the Expr path). The audited
  Expr path is only reached when FLINT returns `-1` ("out of scope"). The
  correctness of the fast path depends entirely on `flint_rde_base_solve_fg`'s
  scope claim (f polynomial over Q(x), univariate) being exact — i.e. it must
  return `-1`, never `0`, for any f/g it does not fully handle (in particular
  the f-rational-with-positive-integer-residue case where weak normalization is
  non-trivial). That routine lives in another file and should be audited
  separately to confirm it never returns a false `0`.

---

## Summary table

| # | Location | Book step | Code | Severity |
|---|----------|-----------|------|----------|
| 1 | rde_weak_normalizer 1079–1122 | ExtendedEuclidean + residue product | `gg==1` gate; provably always 1 | OK-BY-DESIGN (cosmetic docstring) |
| 2 | rde_base 1208 | a,b,c ∈ k[t] (Cor 6.1.1) | `rt_is_poly` defensive guard | MINOR (defensive) |
| 3 | rde_base 1213 | `if d_b=d_a−1` | extra `&& da>=1`; excluded case identical | OK-BY-DESIGN |
| 4 | rde_polyrischde_integrate 1028 | `if n=0 then m←0` | folded into `m>n` decline | OK-BY-DESIGN |
| 5 | rde_base 1139–1154 | (n/a) | FLINT fast-path authoritative `nr>=0` | verify downstream |

No CRITICAL or MAJOR findings. SPDE, RdeNormalDenominator, RdeBoundDegreeBase,
NoCancel1, and the WeakNormalizer resultant/residue product are all
textbook-exact for the base field.
