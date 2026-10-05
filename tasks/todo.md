# Residue method: second round of general algorithmic extensions

Plan: `~/.claude/plans/pasted-content-id-e7e1-the-following-humming-bachman.md`

Policy: correct-by-construction tests (exact string or `{FreeQ[r,Integrate],
Chop[N[r-ref]]}`→`{True,0}`; no committed NIntegrate). Verify each case against
NIntegrate during dev. Bump `src/version.h` + tag per substantive phase. Route
warnings via `mth_message`/Quiet.

All reference values twice-confirmed vs NIntegrate (see plan triage table).

## Phase 1 — Fourier conjugation correctness (In[14] WRONG=0, highest severity)  [DONE v0.277]
- [x] TWO bugs found. (a) Fourier closing: `ReplaceAll[I->-I]` dropped stored -I atoms
      with a mixed-sign (axis+enclosed) residue sum -> switched Cos/Sin extraction to
      `ComplexExpand[Re/Im[J]]`, removed res_conjugate. (b) CORE: `FactorTerms` returned
      content 0 from a nonzero Laurent input (b^-1 left by a failed Together), poisoning
      Simplify (Simplify[(I-I Exp[-a b])/b^2]->0). Fixed facpoly_factorterms.inc: content
      of nonzero is a unit (1), never 0.
- [x] Full Fourier regression green (In[2], x Sin/(x²+b²), Sin[x]/x, Sin[x]/(x(x²+1)),
      Cos/(x²+b²)^2); FactorTerms/FullSimplify/Together/ratcanon/limit suites green.
- [x] Pins: In[14] (Pi(1-E^(-a b)))/b^2 (exact + non-vacuous); In[15] numeric pin;
      FactorTerms no-spurious-zero regression. version 0.277 + docs + changelog + tag

## Phase 2 — Symbolic-param plumbing robustness (In[5])  [DONE v0.278]
- [x] Refine-based strict-sign convergence gate: res_region_neg + g_assume (raw
      assumptions mirror g_inst); replaced the Gaussian param_interval(A).hi<0 gate.
      Proves A=-a<0 under a>0, AND works on the Element[_,Reals]/FindInstance path.
- [x] CORE Exponent fix (exponent.c): Exponent[-a x^2 + I x, x] gave 1 (base_exp_in_monomial
      didn't recurse into the non-flat Times[-1,Times[a,x^2]] a complex sibling produces).
- [x] build_instantiation per-param: NOT NEEDED (Refine gate sidesteps it). Left unchanged.
- [x] Pins: In[5] whole/half + Element[b,Reals] (test_gaussian); Exponent nonflat_times
      (test_exponent). version 0.278 + docs + changelog + tag. Suites green.

## Phase 3 — Assumption-aware closing simplification  [DONE v0.279, partial]
- [x] res_close_positive helper (Simplify[Refine[Simplify[PowerExpand[e]],g_assume]])
      under g_inst&&g_all_pos; wired into close_algebraic AFTER RootReduce.
- [x] In[7] 1/((x²+a²)(x²+b²)) on {0,∞} -> Pi/(2ab(a+b)) CLEAN (test_rational_symbolic).
- [~] In[17] Log[x]/(x²+a²) and In[3] x^a Log[x]/(1+x²)²: DOCUMENTED correct-but-unsimplified.
      Blockers are Simplify/Refine gaps (Arg[I a] doesn't reduce to Pi/2 for symbolic a>0;
      (1-E^(2πia)) doesn't collapse), NOT residue-method gaps. Not forced per plan. a0/general_a
      closers left unchanged (res_close_positive wouldn't help -- Arg[I a] survives PowerExpand).

## Phase 4 — Periodic-strip hyperbolic rectangle + scale normalization  [DONE v0.280]
- [x] residue_family_hyperbolic_strip: N(x)/Cosh[b x], quasi-period fold -> single pole,
      per-term (Pi/b)Sec[alpha Pi/2b]; TrigToExp numerator; Refine convergence. Wired
      after rectangular. In[4] Sec[a/2], In[21] (Pi/b)Sec, In[18] (Pi/b)Sech CLEAN.
- [x] scale-normalization in rectangular (x->x/c for common positive c): In[12]
      x/Sinh[a x] {0,Inf} EVALUATES = Pi^2/(4a^2) (form messy for symbolic a -- keyhole
      Arg surface, same as In[17]; numeric pin). Concrete x/Sinh[2x]=Pi^2/16 clean.
- [x] BONUS (sound, general): Refine[Arg[pos]]->0, Arg[neg]->Pi, Arg[i*pos]->+-Pi/2
      (simp_assume_rewrite.c) -- prov_pos descends Times/Power. NOT applied to keyhole
      closers (PowerExpand there is unsound: Arg[-4a^2]=Pi -> Arg[a^2]=0 drops a factor).
- [x] Pins in test_hyperbolic_strip + refine Arg pins. version 0.280 + docs + changelog + tag.

## Phase 5 — Mellin after power substitution (In[9], In[10])  [DONE v0.281]
- [x] residue_family_mellin_power: C x^(mu-1) G(kappa x^nu), G in {Exp,Sin,Cos}, u=x^nu
      -> (C/nu) M[G](mu/nu). Exp->Gamma c^-s (Re s>0); Sin/Cos->Gamma k^-s {Sin,Cos}[Pi s/2]
      (0<Re s<1). res_region_* convergence gates.
- [x] In[9] x^(s-1)e^(-a x^2) = (1/2)a^(-s/2)Gamma[s/2] (needs s>0 not Re[s]>0 -- Refine
      can't combine a Re[...] conjunct with a sibling); In[10] x^p Sin[x^2] clean; Cos sibling;
      Sin[x^3]=Gamma[4/3]/2. test_mellin_power. version 0.281 + docs + changelog + tag.

## Phase 6 — Sector symbolic powers (In[19]) + generalized Beta (In[16])  [DONE v0.282]
- [x] monomial_split_sym (symbolic numerator exponent); sector gates via res_region_*
      (c>0, n>0, 0<s<n). In[19] x^(2m)/(1+x^(2n)) = (Pi/2n)Csc[Pi(2m+1)/2n] with n>=m+1.
      NOTE: strict n>m needs integer-gap (n>m & ints => n>=m+1) Refine lacks -> use n>=m+1.
      Removed now-unused integer monomial_split.
- [x] residue_family_beta: x^a/(x+b)^c (non-integer c, branch pt at -b) via x=b t ->
      b^(a+1-c)Gamma[a+1]Gamma[c-a-1]/Gamma[c]. Tried last on half-line (integer-c -> mellin).
      In[16] clean. test_beta + In[19] in test_sector. version 0.282 + docs + changelog + tag.

## Phase 7 — Trig combined b Cos+c Sin + warning-leak fix (In[11])  [DONE v0.283]
- [x] Quiet apply_inst (+ res_reim_direct/numeric_double N-probes): the degenerate-point
      1/0 (pole denom vanishes at b=c=0) no longer leaks Power::infy/Infinity::indet.
- [x] residue_family_trig wrapper: on undecidable classification (coupled FindInstance
      mode only) retry at an all-nonzero representative (build_nonzero_inst from g_assume).
      Sound: full-period rational is analytic off circle-poles (core flags idiv). Confined
      to trig -> Fourier under-constrained declines (Cos[kx]/(x^2+a^2), k>0) unchanged.
      NOTE: a global nonzero-FindInstance was UNSOUND (broke that control) -> reverted;
      localized to trig instead.
- [x] In[11] = 2πa/(a²-b²-c²)^(3/2) (form 16πa/(4a²-4b²-4c²)^(3/2), correct). Pin in
      test_trig_symbolic (non-vacuous). version 0.283 + docs + changelog + tag; check-messages OK.

## Phase 8 — Hard tail  [DONE v0.284]
- [x] In[13] (Cos[a x]-Cos[b x])/x²: SOLVED (not a decline). New residue_family_fourier_sum
      (+fsum_term) handles R(x)·Σ C_i K[ω_i x] (K all Cos or all Sin, ω_i same sign); lifts
      to Σ C_i e^(iω_i x), closes half-plane. Axis-pole-by-cancellation gate: admit a real
      pole of R only when Limit[f,x->z0] of the ORIGINAL integrand is finite (the numerator
      cancelled the singularity), half-residue weight 1 vs enclosed weight 2. Half-line
      π(b-a)/2, whole-line π(b-a). Single Cos[a x]/x² (no cancellation) DECLINES (divergent).
- [x] In[20] arctan: HONEST DECLINE. Branch-cut (ArcTan) integrand; value (π/2)Log[1+a/b]
      is Feynman parametric-diff, not residue. Pinned as unevaluated-form control.
- [x] Pins: In[13] exact half+whole + non-vacuous numeric + divergence-decline control +
      (Sin[a x]-Sin[b x])/x->0 sibling (test_fourier_symbolic); In[20] unevaluated pin
      (test_honest_declines). docs multi-freq-Fourier bullet + In[20] note; changelog.
      version 0.284 + tag. Residue suite + integrate/series/limit/fresnel/beta/exponent/
      factor/simplify suites green.

## Review

All 21 probe cases (In[1]–In[21]) from the v0.276 transcript are accounted for.
In[6] is a parser error (user typo `Method -> ]`), not a math case.

OK on arrival (no change needed): In[1] πCsc[πa], In[2] πe^(-ab)/b, In[8] (π/n)Csc[π/n].

Fixed (correctness / new families), one commit + tag each:
- In[14] WRONG=0 → (π/b²)(1-e^(-ab))  [v0.277] — TWO bugs: Fourier conjugation via
  ReplaceAll[I->-I] (now ComplexExpand[Re/Im]); CORE FactorTerms returning 0 content
  for a nonzero input (facpoly_factorterms.inc: content of nonzero is a unit).
- In[15] rode the conjugation fix (cleaner form).  [v0.277]
- In[5] ½√(π/a)e^(-b²/4a)  [v0.278] — Refine-based sign gate (res_region_neg + g_assume)
  replacing the brittle interval test. CORE Exponent fix (recurse into non-flat Times).
- In[7] π/(2ab(a+b))  [v0.279] — res_close_positive (PowerExpand+Refine) in close_algebraic.
- In[4] Sec[a/2], In[18] (π/b)Sech[πa/2b], In[21] (π/b)Sec[πa/2b], In[12] π²/(4a²)  [v0.280]
  — new residue_family_hyperbolic_strip (quasi-period fold) + scale normalization in
  rectangular. BONUS sound Refine[Arg[pos/neg/i·pos]] rewrite.
- In[9] ½a^(-s/2)Γ(s/2) (s>0), In[10] ½Γ((p+1)/2)Sin[π(p+1)/4]  [v0.281] — new
  residue_family_mellin_power (u=x^ν → (C/ν)M[G](μ/ν)).
- In[19] (π/2n)Csc[π(2m+1)/2n] (n≥m+1), In[16] b^(a+1-c)Γ(a+1)Γ(c-a-1)/Γ(c)  [v0.282] —
  sector symbolic-exponent (monomial_split_sym) + new residue_family_beta.
- In[11] 2πa/(a²-b²-c²)^(3/2)  [v0.283] — unit-circle b Cos+c Sin degenerate-instance retry
  + Quiet the classification probe (no Power::infy leak; message-routing correctness).
- In[13] π(b-a)/2  [v0.284] — new residue_family_fourier_sum (multi-frequency diff).

Documented correct-but-unsimplified (Refine/Simplify gaps, NOT residue-method gaps; not
forced per plan):
- In[3] x^a Log[x]/(1+x²)²: (1-E^(2πia)) factor does not collapse.
- In[17] Log[x]/(x²+a²) = πLog[a]/(2a): Arg[I a] does not reduce for symbolic a>0.

Honest decline (documented + pinned as a negative control):
- In[20] ArcTan[a x]/(x(1+b²x²)) = (π/2)Log[1+a/b]: branch-cut/Feynman, not residue.

Core bugs fixed in passing (each with its own regression test):
- FactorTerms fabricated content 0 from a nonzero input (poisoned Simplify).  [v0.277]
- Exponent miscounted a non-flat Times[-1, Times[a, x^2]] (gave 1).  [v0.278]
- Refine now reduces Arg[positive/negative/i·positive] under assumptions.  [v0.280]

No regressions: residue, integrate (newton_leibniz/ramanujan/symmetry/fresnel/beta),
series, limit, factor_terms, exponent, simplify, beta suites all green. check-messages OK.
