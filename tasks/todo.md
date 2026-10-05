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

## Phase 2 — Symbolic-param plumbing robustness (In[5]; enables 9/18)
- [ ] build_instantiation per-param bounds (Element[_,Reals]/Re[s]>0 don't nuke others)
- [ ] Refine-based strict-sign convergence gate (shared helper); replace Gaussian interval gate
- [ ] Pin In[5] ½√(π/a)e^(-b²/4a) on {0,∞}. version + docs + changelog + tag

## Phase 3 — Assumption-aware closing simplification (In[3], In[7], In[17])
- [ ] res_close_positive helper (PowerExpand + Arg/Abs/Log-of-I·pos + Sqrt[-c²]) under g_all_pos
- [ ] wire into a0_keyhole_log, keyhole_log_general_a, close_algebraic/half-line; res_powerclean poles
- [ ] Pins: In[17] π Log[a]/(2a); In[7] π/(2ab(a+b)); In[3] clean. version + docs + changelog + tag

## Phase 4 — Periodic-strip hyperbolic rectangle + scale normalization (In[4],18,21,12)
- [ ] residue_family_hyperbolic_strip (quasi-period shift, geometric 1/(1-λ))
- [ ] scale-normalization pre-step (u = c x) → In[12] via existing a=1 route
- [ ] Pins: In[4] Sec[a/2]; In[18] (π/b)Sech[πa/2b]; In[21] (π/b)Sec[πa/2b]; In[12] π²/4a². version+docs+tag

## Phase 5 — Mellin after power substitution (In[9], In[10])
- [ ] residue_family_mellin_power: u=x^ν → (1/ν)M[f](μ/ν), f∈{Exp,Sin,Cos}
- [ ] Pins: In[9] ½a^(-s/2)Γ(s/2); In[10] ½Γ((p+1)/2)Sin[π(p+1)/4]. version+docs+tag

## Phase 6 — Sector symbolic powers (In[19]) + generalized Beta (In[16])
- [ ] extend residue_family_sector (symbolic num/den exponents, Refine convergence)
- [ ] residue_family_beta: x^a/(x+b)^c → b^(a+1-c)Γ(a+1)Γ(c-a-1)/Γ(c)
- [ ] Pins: In[19], In[16]. version+docs+tag

## Phase 7 — Trig combined b Cos+c Sin + warning-leak fix (In[11])
- [ ] amplitude-phase pre-normalization (β Cos+γ Sin → R Cos[θ-φ]); shift θ→θ+φ
- [ ] Quiet the internal N-probe in res_reim/apply_inst (no Power::infy leak)
- [ ] Pin In[11] 2πa/(a²-b²-c²)^(3/2). version+docs+tag; make check-messages

## Phase 8 — Hard tail
- [ ] In[13] (Cos[a x]-Cos[b x])/x²: best-effort multi-frequency Fourier diff; else decline
- [ ] In[20] arctan: honest decline (parametric-diff, not residue); pin negative control
- [ ] docs note + tests. (bump only if src behavior changed)

## Review
(filled at end)
