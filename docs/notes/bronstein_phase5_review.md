# Adversarial review of phase5_theory.tex (E2, E1, Phase 5 plan)

Reviewed against [SI] = Bronstein, *Symbolic Integration I*, 2nd ed. (2004); [B] = Bronstein 1990,
*Integration of Elementary Functions* (JSC), Part 3 (the Risch equation); [B91] = Bronstein, *The Risch
Differential Equation on an Algebraic Curve*; [T] = Trager 1984 thesis. Page numbers are the printed ones.

Verdict: **one FATAL finding** (a false "no" from the plan as written), **seven GAPs**, and several
MINOR findings. The core admission theorem (Thm E2struct) is applied within its hypotheses. I found no
case where two monomials are wrongly declared independent as long as Const(k_{i-1}) = Q really holds.
The false-certificate risk is in what happens *after* a dependence is found (A3).

---

## FATAL

### F1. (A3) elimination puts non-rational constants into the field, but E2/E1 fix Const = Q everywhere

**Where.** E2.2 (A3); Prop E2field(iii); Lemma E2const(i) ("R/lc(R) ∈ Q[z]"); E2.3 (PRDE/LI with
c ∈ Q^m); E1.3 (`brBCLimited` tests "(a − Dv)/η ∈ Q"); Prop E1ld(iii) (c ∈ Q*); Plan item 2 ("a rational
vector eliminates").

**The problem.** (A3) says a dependent logarithm becomes Σ r_g t_g + Σ r_g w_g + κ, and that Prop T10
"holds verbatim". By T10 the offset is (Log κ + 2πi m)/d, where κ ∈ Const(k_{i-1}) = Q*. So the offset is
log 2, iπ, and the like. A dependent exponential with integer r_g gives the factor e^c with c ∈ Q, so
E = e^1 can appear. After the substitution the integrand lies in C(x)(y)(t…) with C = Q(log 2), or
Q(i)(π), or Q(e). Every oracle in E2.3 and E1 still tests membership in **Q**:
- Lemma E2const(i) rejects residues in Q(τ)\Q.
- LI_0 (`brBCLimited`) and PRDE return "no solution" when the needed constant is 1/log 2.
- ROOT_0 asks for c ∈ Q*.

The original T10 section does carry C = Q(α)(τ), and it relies on Prop T1 for genericity. But E2 never
imports that. OL1 says Prop T1 is unproved on the curve, and OL2 says number-field constants are
unchecked. The iπ case also puts i into the constants, which is exactly the situation that (A0)
declines.

**Counterexample.** Let z = √(1+x³), t1 = log(x+z), t2 = log(2x+2z), and take

  g3 = D(log(x+z)) · log(x+z) / (log(2x+2z) − log(x+z)).

Then t2 − t1 = log 2 on principal branches. The antiderivative is log(x+z)²/(2 log 2), which is elementary.
- Admission finds Dt2 = Dt1 (r = 1). (A3) eliminates t2 = t1 + log 2.
- The integrand becomes (1/log 2)·t1·Dt1, a polynomial of degree 1 in the primitive t1.
- Thm E2level step 3 calls LI_0((1/log 2)·Dt1; Dt1). The answer needs c = 1/log 2 ∉ Q.
- With c restricted to Q this gives "no solution", and by Thm E2cert(b) that is a certificate. It is
  false.

Checks with the binary:
- The genus-0 analogue (y = √(1+x²)) is integrated correctly today: `Method -> "Bronstein"` returns
  log(x+y)²/(2 log 2), and `PossibleZeroQ` confirms it.
- The genus-1 g3 is currently *declined* with "transcendental Risch over an algebraic base". So Phase 5,
  by lifting `Tower.m:205`, is precisely what exposes it.
- The same happens with log(x − z) beside log(x + z) (product −1 when z = √(1+x²); offset iπ on the real
  line), and with e^{xz+1} beside e^{xz} (factor e).

**Fix.** Either (a) decline every (A3) elimination whose constant offset is nonzero, or (b) carry
C = Q(α)(τ) through every oracle. Option (b) touches Lemma E2const(i), `brDFSpan`'s projection, PRDE/LI
constants, `brBCLimited`'s "∈ Q" test and ROOT_0's c ∈ C*, and it needs OL1 and OL2 discharged first. A
unit test on g3 should be added to the M1 gates.

---

## GAP

### G1. Prop E1ld(ii): "projecting on a Q-basis containing 1" is wrong if read literally, and so is the source

The note cites [SI] §7.3 last paragraph (p. 255). That paragraph says that projecting the C-coefficient
system "on 1 yields a system with … the same rational solutions". This is false. Take C = Q(√2) and the
equation √2·r1 − √2·r2 = 0. Its rational solutions are r1 = r2. Its 1-projection is 0 = 0, which every
pair satisfies.

The correct reduction projects on **every** basis element and keeps all the projections.

If only the 1-projection is used, (ii) can report a spurious (n0, m0, U). (iii) then returns a witness z,
and LD_j answers "yes" falsely.
- In a degree bound that is harmless: by remark (2) of E1.1, an extra candidate only weakens the bound.
- In the primitive cancellation case, PolyRischDECancelPrim rewrites the equation using z. A wrong z gives
  a wrong equation, and so a possible false "no solution".

**Required.** Project on all basis elements, *and* check n0·α = DU/U exactly in normal form before using
any positive auxiliary answer. This is cheap and closes the case. E2.4 already demands soundness of
auxiliary answers in both directions; this is the place where it can actually fail.

### G2. Prop E1ld(ii) needs Thm E2struct and Prop E2span over Const = Q(c); only the Q case is proved or planned

F = k_j(c)(admitted log v_l) has constants Q(c). The note only asserts "constants of k_j⊗Q(c) are Q(c)".
That claim is true: [SI] Cor 3.3.1 (p. 86), plus the fact that Q is algebraically closed in k_j. But the
admission test over F and the normal-form span test with Q(c)-coefficients are never stated, and
`brDFSpan` as specified works over Q only.

[SI] Thm 9.3.1 (p. 283) does allow any constant field C, so the mathematics goes through. What is missing
is the algorithm (normal forms over Q(c), with the projection done as in G1) and its proof.

### G3. ROOT_0 is not given as an algorithm

The principal-divisor test with a generator has three missing pieces:
- It needs places over Q̄, including the places at infinity.
- It needs an integral basis that is normal at infinity ([T] Ch. 2 §3; [T] Ch. 6 §1, pp. 64–65).
- The Möbius frame x ↦ 1/(x−a) must also transform b, the integral basis, and the divisor, with a chosen
  so that x = a is not in the support and is not a branch point.

None of this is written down, and OL1 (specialisation) is open. Thm E2decide assumes ROOT_0 as a
hypothesis, so the "decision procedure" is conditional on an unwritten oracle.

### G4. RDE_0 bounds: per-place values, and [B] Prop 3.7's hypothesis on the integral basis

[B] Prop 3.3(iii) and Prop 3.5 (pp. 160–161) give a candidate order −r·μ(f) *per place*. The plan's G1
reads τ from a resultant over all places above P. That is sound only if the code takes the most negative
candidate (a weaker bound). If it picks one place's value, a pole at another place can be missed, which
is a false "no".

[B] Prop 3.7 (p. 162) also requires w_1, …, w_n to be an integral basis **normal at infinity**. The note
does not say that `brAlgBounds` uses such a basis. Every "no" from RDE_0, and through E1 every "no" above
it, rests on these two points, and they are recorded only as "implementation" (G1/OL).

### G5. "Decision procedure for K" is overstated: in-scope towers are declined

- (A3) declines a dependent exponential with a fractional r. Example with y = √x: t1 = log x, then
  t2 = e^{t1/2}. Here r = 1/2, yet t2 = y ∈ k_0.
- Another example: e^{xy} admitted first, then e^{xy/2}.
- T10 declines when more than one transcendental constant appears.

Risch's algorithm handles both by re-choosing generators (normalising to e^{w/N}) or by accepting the
algebraic element. Thm E2decide is a decision procedure only for towers that are already presented in
admitted, normalised form. Either the front end normalises generators, or the claim must be restricted.
This is a completeness gap, not a soundness one.

### G6. Front end: the spelling of y and stray algebraic constants

Prop E2analytic and Prop T0 assume ι(F) = f on U. Kernel arguments can spell the radical differently, and
these spellings agree with ±y only up to sign, region by region:
- √(4+4x³) = 2y;
- (1+x³)^{3/2} = (1+x³)·y only on principal branches;
- √x·√(1+x) versus √(x²+x), which differ by sign for x < −1.

(A1)'s "argument in normal form" must be backed by a branch check of the substitution. The plan does not
mention one.

Separately, (A0) checks only Y^n − b. An integrand containing I or √2 (e.g. log(x + I·y)) has constant
field Q(i) even when (A0) holds. It must be declined (OL2) or handled. The plan has no such check.

### G7. LI in the primitive degree bound needs an integer, not a rational

The degree bound of [SI] §6.3 for d_b = d_a − 1, d_b = d_a needs β = Dv + m·Dt with **m ∈ Z**. LI returns
c ∈ Q, so the step "keep the candidate only when c ∈ Z" must be explicit.

On its own a non-integer c is harmless if it is discarded (the bound stays valid). It is unsound only if
the c is rounded. The table and remark (1) in E1.1 omit this.

---

## MINOR

1. **E2.1.** The proof of Prop E2field(iii) cites "[T] Ch. 2" for "Q is algebraically closed in k_0 iff
   Y^n − b is absolutely irreducible". That material is [T] Ch. 3 (absolute irreducibility, exact
   coefficient fields, p. 32ff).
2. **E2.1 derivation formula.** "Dt_l ∈ k_{l−1}" is false for an exponential, where Dt_l = t_l·Dw_l ∈ k_l.
   The formula itself is fine.
3. **E1.1 remark (2), sign.** The correct identity is (n′m − nm′)·η = D(u′^n/u^{n′}) / (u′^n/u^{n′}), not
   the logarithmic derivative of u^{n′}/u′^n. The conclusion m/n = m′/n′ is unaffected.
4. **E1.4 "Where Bronstein stops".** It claims the parametric RDE on a curve is "not in [B]". In fact
   [B] Part 3 is stated for the parametric equation (R1) z′ + fz = Σ c_i g_i with c_i ∈ C, and
   γ_p = min_i ν_p(g_i) (p. 160). Lemma E1uniform is essentially [B] Props 3.3–3.7 and should cite them;
   that strengthens the lemma.
5. **[B] Part 3 hypothesis.** It assumes k contains x with x′ = 1, and k = Q does not (E1.3: "k = Q,
   θ = x"). This looks harmless: Lemmas 1.7/1.8 (pp. 122–123) only use the Puiseux coefficients being
   constants. The note should still say so.
6. **Lemma E1uniform proof.** The quoted bound min(0, γ_P + r, −τ_P) omits [B] Prop 3.3(i)
   (β_p < −r ⇒ α_p ≥ min(0, γ_p − β_p)). Monotonicity in γ still holds, so the lemma stands.
7. **E1.5 Example 2.** "The two ramified places over 1+x³ = 0": there are three (plus the one at infinity;
   genus 1 has 2g+2 = 4 branch points). The argument is unaffected.
8. **E2.5 Example 1.** "LI_1(1; η2) gives (b, c) = (x, 0) (unique)": b is unique only up to an additive
   constant. This is harmless, because κ′·η2 is a logarithmic derivative.
9. **Lemma E2const(ii) and the paragraph after it.** "Squarefreeness is unaffected" is true for d itself,
   but the *norm* denominator can fail to be squarefree even when the true denominator is squarefree.
   Example: (t − y)(t − x) gives (t² − b)(t − x)². Hermite reduction is still correct, but the stated
   invariant holds only if `brBCSplitDen` removes y-free factors before taking norms. State that.
10. **Thm E1rde proof.** It describes [SI] Ch. 6 as "undetermined-coefficient solving in the bounded
    space". It is SPDE recursion followed by the base RDE. Only the level-0 step is an ansatz.
11. **Thm E2level step 3 (exponential).** The residue-reduction log terms contribute deg(g_α)·α·Dw to the
    t⁰ coefficient, which feeds Bronstein's IntegrateHyperexponential. This is implicit in "the Laurent
    part" but worth stating.

---

## Checked and correct

- **Thm E2struct (a)–(c).** [SI] Thm 9.3.1 / Cor 9.3.1 (pp. 283–285) hold for any field C and any
  elementary extension of (C(x), d/dx) with Const = C. Algebraic generators are allowed and excluded from
  E and L ((9.6)–(9.7), p. 282). [SI] Thm 5.1.1 and 5.1.2 (p. 133ff) give the monomial criteria, and the
  proof of (c) is right.
- **Probes of the admission test.** These are decided correctly given Const = Q:
  - log(x+√(1+x²)) beside log(x−√(1+x²)): dependent, r = −1 (constants: see F1).
  - log(x+√(1+x²)) beside log(1+x²): independent, since the y¹ and y⁰ coordinates are disjoint.
  - e^{√x} beside e^{2√x}: r = 2, κ = 1.
  - log y beside log b: r = 2, offset 0 for n = 2 on principal branches (for n = 3 it is region-dependent,
    and T10 handles it).
  - e^{(log x)/2 + √x} beside log x: independent, which is correct.
- **Prop E2field.** (i), (ii) and (iv) are correct. (iii) at i = 0 is correct, with the citation fixed
  (MINOR 1).
- **Prop E2span.** Exact; it is the rational case.
- **Lemma E2const.** (i) is the [SI] Thm 5.6.1 criterion, and the constants are correct when
  Const(k_{s−1}) = Q. (ii) is checked root by root. (iii) is valid for monic d.
- **Prop E2analytic.** Uses [SI] Thm 5.1.1/5.1.2 in M(U) correctly.
- **Thm E2cert (a)–(e).** The strong Liouville theorem ([SI] Thm 5.5.3, p. 145) puts v in K, so
  c ∈ Const(k) for degree ≥ 1. Step 4 (remainder h) is valid by definition of an elementary integral.
- **Prop E1ld.** (i) is correct. In (ii), the norm argument uses [SI] Thm 3.2.4(ii) correctly. (iii) is
  correct at both levels. LD from a non-minimal n0 is correct.
- **E1.3 LI_0 completeness.** The residue-by-Q-component argument is valid.
- **E1.1 remark (2).** A spurious PLD candidate only weakens both the special-denominator bound (min) and
  the degree bound (max).

---

## Numerical and symbolic checks (./Mathilda -file)

All of the following returned `True` under `PossibleZeroQ`:
- D(xy) = (2+5x³)/(2y);
- f1 = log(x+y) + a·e^{xy};
- the η2 normal form;
- the RDE solution q = 1/x;
- x·η2 = (3x³+2xy)/(2xy+2(1+x³));
- D[e^{xy}/x + x·log(x+y)] − f1 = x·η2;
- the f2 identity, and D(y + log x) = (3x³+2y)/(2xy).

Also:
- lc(R)⁻¹·R(0) for 1/(y + log x) equals −2xy/(3x³+2y): the difference simplifies to 0.
- Trager declares ∫x·η2 not elementary. This is consistent with the note but not independent of it.
- `Bronstein` certifies f1 and 1/(y + log x) as non-elementary. It currently declines f2, which agrees
  with the note.
- The F1 analogue with y = √(1+x²) has antiderivative log(x+y)²/(2 log 2); `PossibleZeroQ` confirms it,
  and a numerical residual check gives 1e-16. The genus-1 g3 is currently declined.

Hand checks of the valuations (agree):
- E1.5 Example 2: ν∞(q) = 1, so q = 0.
- E1.5 Example 1: L(p₊ + p₋) = span{1, 1/x}, giving q = 1/x.

Scripts are in /private/tmp/claude-501/-Users-user-Mathilda-dev-Mathilda/f0aa3cf6-c83f-4fb9-a39f-4fc388b010e6/scratchpad/p5/.

---

# Second pass (revision 2 of phase5_theory.tex; v1 kept as phase5_theory.v1.tex)

Verdict: **no FATAL finding.** F1 is closed in substance: under (H_Q) with the zero-offset
rule, no level gets a constant field larger than Q. Three GAPs remain:
- the specification does not say that the κ = 1 test is exact;
- the accepted spelling (1/b)^{k/n} ↦ y^{-k} is not an identity on the real interval where
  b < 0;
- the M1 decline list misses the exponential cancellation case, which needs a PLD witness.

The new Lemmas E1res, E1choice and E1inf are correct. Two of them need a hypothesis stated
(uniform ramification) or have a notation clash.

## 1. Is F1 closed?

I looked for every path by which a non-rational constant could enter a level.

| Path | Status |
|---|---|
| Dependent log, nonzero offset (log 2, iπ) | Declined by (A3). Sound. |
| Dependent exp, offset e^c with c ≠ 0, or fractional r | Declined. c is read exactly from the normal form. Sound. |
| Dependent log, κ = 1, m = 0, rational r | Substitutes Σ r_g t_g ∈ k_{i-1}. Constants stay Q. Sound. |
| Constants in the integrand or in kernel arguments (I, √2, π, Root) | Excluded syntactically by (H_Q). Sound. |
| Constants in the radicand (b = 2+2x³; b = 2x² is not absolutely irreducible) | (A0) catches the bad cases. Sound. |
| Residues α ∈ Q̄ of the Rothstein–Trager step at level s | The remainder h − g₂ is in k_{s−1}(t) ([SI] Thm 5.6.1(i), p. 149). Sound. |
| Exponential correction −Σ_α α·deg(g_α)·Dw_s (added in revision 2) | Rational, because deg g_α is constant on a conjugacy class, so the sum is a trace. It must be **computed** as a rational. A RootSum left unevaluated in the t_s⁰ coefficient would hand a non-Q constant to a Q-only oracle at level s−1. **MINOR:** state it. |
| Constants in INT/LI/PRDE solutions (c ∈ Q) | Complete because Const = Q. Sound. |

### GAP S1: the κ = 1 test must be exact; only m may be numerical

Plan item 2 says "κ = 1 and m = 0 at every T10 sample point". Read literally, that makes κ a
numerical quantity too.

κ = v^d / ∏u_j^{d r_j} · ∏ t_k^{−d s_k} is an element of k_{i−1} with Dκ = 0, so by Prop E2field(iv)
it is a rational number. That number must be computed exactly in the brDField normal form and
compared with 1 exactly.

A numerical test at T10's tolerance of 1e-8 would accept κ = (10¹²+1)/10¹², from log((10¹²+1)(x+z))
beside log(10¹²(x+z)). That admits a nonzero offset log κ ≈ 1e-12, which is exactly the F1 situation.

The existing C rule (`rt_dependent_log_rule`, risch_tower.c:1017) does compute κ symbolically
(`Together`, then a free-of-x check). Only the integer m comes from `N[...]` at the sample points. The
note should say the same for the curve path, and that κ must be found constant by `ConstantQ` in normal
form. `Together` on Sqrt-spelled kernels leaves (x−y)(x+y) unexpanded, which is safe (it declines) but
useless.

The m test itself is not exact. It is integer separation: δ is an integer wherever the branches are
analytic, and the check only needs an error below 1/2. At machine precision that is reliable but not
certified. At a sample point that lies on a principal-Log cut (for z = √(1+x³), every real x < −1, e.g.
x = −13/8, sits on the cut of y), rounding picks one adjacent side. The value read is then the δ of an
adjacent open region. So m = 0 still yields a region on which the relation holds, and the certificate
(a statement about the germ, Prop T0) stays valid. **MINOR:** state the integer-separation argument, and
consider interval evaluation.

### GAP S2: (1/b)^{k/n} ↦ y^{−k} is not an identity of functions where b < 0; the declined example is mislabelled

(A1) accepts (1/b)^{k/n} as y^{−k} and calls the conversion "an identity of functions on U". Off the cut
it is not identical: Log(1/b) = −Log b + 2πi wherever b is negative real, so (1/b)^{1/n} = e^{2πi/n}·y^{−1}
there. For b = 1+x³ that set contains the whole real interval x < −1.

Checked with the binary at x = −13/8: √(1/b) = +0.551233 i, while 1/√b = −0.551233 i. Just above the
axis (x = −13/8 + i/1000) they agree.

Consequence: h = √(1/(1+x³)) + 1/√(1+x³) is identically 0 on x < −1, as a real function. Yet
`Method -> "Bronstein"` certifies it "no antiderivative elementary in x" today, because on U off the cut
it equals 2/y. This is correct only as a statement about the germ on a U that avoids b⁻¹((−∞,0)).

Either (A1) must say that U avoids the preimage of the negative axis of b (and T10's real samples
x = −7/19, −13/8, where b < 0 for −13/8, then lie outside every admissible U for this spelling), or
(1/b)^{k/n} must be declined like the other respellings.

Conversely, the note's declined example √(4+4x³) "agrees with ±y only region by region". That is false:
√(4c) = 2√c holds for every c under the principal branch, so the respelling is a global identity.
Likewise b^{k/n} = y^k holds globally, so (1+x³)^{3/2} in OL6 needs no branch check. Declining them is
safe; only the stated reasons are wrong (**MINOR**).

## 2. The new lemmas

**Lemma E1res — correct.** dw is exact, so its residues are zero; res_P(dU/U) = ord_P U ∈ Z; divide by n.

The use in Prop E1ld(ii), first case, is correct and needs less than stated: it is enough that each c_l
lies in the Q-span of the residues, not that the c_l form a basis. A c_l ∉ Q then forces some residue
∉ Q, so there is no radical solution. The lemma covers places at infinity and is restricted to level 0,
consistent with OL3.

**Lemma E1choice — correct**, with one hypothesis left unstated.

(a) holds when the exponent e used to clear P is exactly the pole order at P. Then
P^e·∏_{pole branches}(z − u_σ) tends to a nonzero constant, because the coefficient with the largest pole
has exactly one term of full order, so nothing cancels. The limit polynomial is
C·∏_{regular branches}(z − u_σ(ξ)). Clearing with a larger power makes T ≡ 0, which the code reports as
`Unknown` (I1): safe.

The derivation α_p = −u(p) with u = r·f·P/P' matches [B] Prop 3.3(iii) (p. 160), and only positive
integer values matter. In (b), the case-(i) inequality γ − β > γ + r is right, and taking the most
negative candidate over the places above P is right. γ_P ≤ γ_p holds because the w_l are integral at
finite places.

Unstated: χ_P uses a single r for all places above P. That is true here because Q̄(x)(y)/Q̄(x) is
Kummer, hence Galois, so ramification is uniform above each point, and conjugate roots of P behave
alike. The same holds at infinity. **MINOR:** state it.

The infinity rule (minimum over rational roots a of r·a, with [B] Prop 3.5 at k = Q, where μ(u) = u for
u ∈ Q) is valid.

**Lemma E1inf — correct.** The x^e coefficient Σ_{E_i=e} lc(h_i)·lc(d_i)⁻¹·v^i is a nonzero polynomial
in v of degree < n. The v range over the n distinct roots of v^n = lc(b), so on some branch it is
nonzero. That branch lies in some place above infinity, where ν = −r·e ≥ b_Z. Lower-order terms of the
(1+O(1/x))^{1/n} expansion do not interfere.

The lemma does need curveBasis to be an integral basis of the integral closure, so that
Z = Σ h_i w_i with h_i ∈ Q[x]. It also needs uniform r at infinity, which holds as above.

**MINOR, notation:** the letter b means the radicand, the infinity bound in Lemma E1choice, and the base
of b_Z = b − r·deg 𝒟 in the remark after Lemma E1inf, and d means both deg b and the d_i. Rename.

Lemma E1uniform now rests on these two lemmas, plus monotonicity in γ. Correct.

## 3. Is the M1 scope consistent with the theorems?

I traced every "no" M1 can emit:
- the residue criterion at levels 1 and 2 (Lemma E2const under (H_Q));
- LI_1 via PRDE_1 ([SI] §7.1, whose parametric cancellation needs only PRDE_0 (p. 241), plus
  Lemma E1uniform);
- RDE_1 ([SI] Ch. 6, using RDE_0 via Lemmas E1choice and E1inf, LI_0 via `brBCLimited`, LD_0/PLD_0 via
  Prop E1ld at j = 0);
- INT_0 (Trager, assumed).

Each is covered under its stated hypotheses, with these exceptions.

### GAP S3: the exponential cancellation case needs an n = 1 PLD witness, and the decline list does not cover it

[SI] PolyRischDECancelExp (p. 213) branches on "b = Dz/z + m·Dt/t for some z ∈ k*, m ∈ Z", and then
**uses z** (it integrates c·z·t^m). When that test fails, it runs the top-down RischDE(b + mη, ·) pass.
The top-down pass relies on uniqueness of each coefficient's solution (Lemma 5.9.1), and uniqueness fails
exactly when a homogeneous solution exists.

Prop E1ld(ii) yields a radical solution (n0, m0, U). When n0 > 1 and m0/n0 ∈ Z, an n = 1 solution
exists, but recovering z needs ROOT_0(U, n0), which is OL5. If the implementation reports "not of the
form" and falls through to the top-down pass, it can return a **false "no solution"**, and through
Thm E2cert(c) a false certificate.

This case arises easily. Trager's choice of basis scaling makes n0 > 1 whenever its logand is a power,
e.g. c = 1/2 with logand u² gives n0 = 2 although u ∈ k_0.

The M1 decline list says only "an LD call whose radical solution has n0 > 1". It must also cover every
PLD call whose n = 1 witness is needed, i.e. the exponential cancellation branch. That branch in
practice needs the witness with m = m0/n0. The PLD uses in bounds and in the special denominator only
need m0/n0 and are fine.

Related **MINOR** points:
- The table in Thm E1rde omits that the z-branch of both cancellation cases needs a "derivative in k"
  test: Dp = c·z·t^m with p ∈ k(t) ([SI] §5.12 variant), i.e. INT_{j−1} with an empty log part, plus
  RDE_{j−1} for the Laurent coefficients.
- The M1 decline list should also name "rational log coefficients from INT_0 whose logands are not given
  in k_0" (RootSum output). The third case of Prop E1ld(ii) requires v_l ∈ k_0, which is I2.

### Other M1 items: consistent

- Thm E2cert now assumes (H_Q), admission, and zero-offset elimination, matching M1's "decides" list.
- PLD is called only at level 0, consistent with OL3.
- Every positive auxiliary answer is verified exactly, which closes first-pass G1.
- An LI c enters a bound only if c ∈ Z_{>0}, which closes G7.
- ROOT_0 and `Unknown` are declined, settling G3 and G4 as declines.
- OL6 makes G5 and G6 declines.
- The worked examples are unchanged; Example 2's LD_0(−Dy) needs no witness and no ROOT. They were
  verified in the first pass.

## Second-pass summary

| ID | Class | Item |
|---|---|---|
| S1 | GAP | κ = 1 must be an exact normal-form test; only the integer m may be numerical (integer separation; cut points are safe under germ semantics). |
| S2 | GAP | (1/b)^{k/n} ↦ y^{−k} differs by e^{2πi/n} on b⁻¹(negative reals), e.g. real x < −1 for 1+x³. The certificate holds only for U off that set (h = √(1/b) + 1/√b is ≡ 0 on x < −1 yet certified). Restrict U or decline the spelling. √(4+4x³) and (1+x³)^{3/2} are global identities, so the stated reasons for declining them are wrong (MINOR). |
| S3 | GAP | Exponential cancellation needs an n = 1 PLD witness. With n0 > 1, ROOT_0 is required, so M1 must decline there, not fall through to the top-down pass (otherwise a possible false "no"). |
| — | MINOR | Exponential log-part correction must be computed as a rational trace. Uniform ramification is an unstated hypothesis of Lemmas E1choice and E1inf. b/d notation clash. The table in Thm E1rde omits the cancellation z-branch's dependence on a "derivative in k" test. Decline list should name RootSum logands with rational coefficients. |

# Third pass (2026-11-22): Lemma E1irr, which closes I2 (Phase L1)

Reviewer: an independent agent given only the note's text and the sources. Verdict: **no
FATAL**, one GAP, five MINOR. All are addressed in the note's revised "Use" paragraph.

| ID | Class | Item | Resolution |
|---|---|---|---|
| T1 | GAP | The certificate relies on the Q-independence of 1, beta_1..beta_r. With a non-canonical basis (sqrt(3+2 sqrt 2) and -sqrt 2 taken as independent, though they sum to 1), or a reducible "minimal" polynomial, a false "no" results. | Use step 2: one primitive element with an irreducible, factored minimal polynomial; unique coordinates; every logand rational in x, y over K, with no other radical. |
| T2 | MINOR | RootSum coefficients live in the splitting field, not Q(z). | RootSum declines (step 1). |
| T3 | MINOR | The W0 recipe was circular. | Divide the canonical normal form (monic d) by a nonzero coefficient; W/c = W0/r is in k0 (step 4). |
| T4 | MINOR | Test sum N s_lk Dv_l/v_l = 0 instead of forming V_k. | Step 3. |
| T5 | MINOR | An ArcTan in Trager's surface output gives coefficients +-i/2, so K must contain i. | Step 1 rewrites ArcTan/ArcTanh as logs. |
| T6 | MINOR | (a)'s LI_0 sentence is valid only for eta_i with rational residues; Lemma E1res's places are over Qbar. | Both stated. |

Checked and sound:
- (a): the residue formula at every place, including infinite and ramified places, and the zero
  of a nonconstant V_k on the complete curve;
- (a), the LI_0 claim, within the stated restriction;
- (c): Gal(Kk0/k0) = Gal(K/Q) by (A0); Const(Kk0) = K; the cocycle identity; Hilbert 90.

The remaining false-certificate surface is exactness only: an inexact zero test, an unverified
Liouville identity, or (a) applied to an eta with irrational residues.

# Fourth pass (2026-11-22): Corollary E1hom, the n = 1 witness without ROOT_0 (Phase L2)

Reviewer: an independent agent. It re-derived the local analysis of z' + f z = 0 from scratch
and checked it against brAlgBounds, brAlgRDEHom, brDTHomRDE0, brDTDegreeBound and
brDTPolyRDE. Verdict: **no FATAL, no GAP**.

Checked:
- **Regular places, ramified places included:** no pole and no zero is possible.
- **Poles of f:** beta < -r leaves no solution; beta = -r gives alpha = -tau, with tau a
  positive integer exactly as brAlgTauCandidates keeps; -r < beta < 0 gives alpha = 0.
- **Infinity:** alpha = r a.
- **Bound data:** Lemma E1inf holds; the pole factors are complete with g = 0, on the Kummer
  integral basis.
- **Solution space:** dimension at most 1 by (A0); the null-space correspondence is exact.
- **Use paragraph:** in the exponential cancellation case only i = m0/n0 can occur, so the
  top-down pass is valid with no witness; in the primitive case m = 0; in the degree bound,
  no witness means the beta branch does not apply.

MINOR findings, all addressed in the note:
- stale "declined" text in remark (3) and in the ROOT_0 paragraph;
- the corollary must name brAlgRDEHom, since brAlgRDE short-circuits on g = 0;
- the (H_Q) dependence of (i) is now stated;
- the I1 path looks unreachable for the tau candidates (kept as Unknown);
- soundness still rests on Prop E1ld (and Lemma E1irr) for the "PLD says none" steps.

Found while implementing: a bare Return inside With in the cancellation branch dropped the
witness's answer, and the top-down pass then said "no solution". That was a false certificate,
caught by a unit test before release. A source scan then found three more such Returns (inside
Do in AlgebraicPart.m, inside With in TorsionModP.m), all fixed in v0.6692.

# Fifth pass (2026-11-22): Lemma E2off, one transcendental offset (Phase L3, theory only)

Verdict: **no FATAL**. Lemma (a)-(c) is correct: the P^D argument holds, the isomorphism is
differential, and g3's integral is t1^2/(2 tau0). Four GAPs, all in the scope and the code
plan, are resolved in the revised note:

| ID | Class | Item | Resolution |
|---|---|---|---|
| L3-1 | GAP | Prime factorisation ignores the sign: Log(-2) = log 2 + i pi needs two constants | kappa_i > 0 with m = 0 (tau0 = log kappa0), or kappa_i = +-1 with any m (tau0 = i pi) |
| L3-2 | GAP | "m = 0" excluded the i pi case | stated as above, with T10's reading of m |
| L3-3 | GAP | Exponential offsets are multiplicative factors e^c | all c_i in (1/N)Z, tau0 = e^(1/N) |
| L3-4 | GAP | After elimination, later kernels can bring algebraic constants (e^(t2/2) = Sqrt[2] e^(t1/2)) | declined until (A3) reads c in Q(tau) |
| — | MINOR | Const(K(tau)) = Q(tau) unproved; the embedding for RootSum constants; the admitted-tower claim | stated |
| — | code | "in Q" (integers, exponents, torsion) vs "in Const(K)" were conflated; RootSum residues live over the algebraic closure of Q(tau); missing: brBCLimited, ROOT_0, E1ld(ii), ConstantQ, the (A3) tests, torsion over Q(tau), zero tests, T10 | split into separate lists |

Status: theory reviewed; the code is a staged follow-up. No corpus case of 256-260 reaches
the offset declines.
