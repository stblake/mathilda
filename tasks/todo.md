# Task: General algorithmic extensions to the Mellin-transform Integrate method

Plan: `~/.claude/plans/pasted-content-id-9f48-after-our-frolicking-aurora.md`
Target file: `src/calculus/integrate_ramanujan.c` (+ tests, spec, changelog, version).

Engine is correct-by-construction (NO runtime NIntegrate — project rule, file
header L36). Numeric verification lives in `tests/test_integrate_ramanujan.c`.

## Phase G1 — single functions → recognized kernel (low risk)
- [ ] EllipticK/EllipticE → 2F1 reduce rules (closes In[16])
- [ ] BesselY → BesselJ combination reduce rule (closes In[17])
- [ ] rec_gamma_upper: Gamma[nu, a x] → Gamma(s+nu)/s (closes In[9])
- [ ] Extend reduce_to_hypergeometric cheap-guard head list
- [ ] Tests + numeric verification for In[9],[16],[17]

## Phase G2 — hyperbolic Dirichlet recognizer (medium risk)
- [ ] rec_hypergeom_dirichlet: Csch[a x] → 2 Gamma(s)(1-2^-s)Zeta(s) (In[13])
- [ ] Sech[a x] → 2^(1-s) Gamma(s) LerchPhi(-1,s,1/2) (In[12])
- [ ] Tests + numeric verification for In[12],[13]

## Phase G3 — Mellin-convolution product engine (headline; high risk)
- [ ] G3a: convolution scaffold + equal-scale Bessel products (In[1-3])
- [ ] G3b: unequal-scale Exp×Bessel → 2F1/1F1 (In[4-5])
- [ ] G3c: AiryAi^2 reduce rule (fallback convolution) (In[14])
- [ ] Update file-header scope comment (L30-37) — products no longer all NULL
- [ ] Tests + numeric verification for In[1-5],[14]

## Cross-cutting
- [ ] docs/spec/builtins/calculus.md — new transforms + fix scope note L2072-2074
- [ ] docs/spec/changelog/2026-10-05.md — changelog sections
- [ ] src/version.h bump per substantive phase
- [ ] make clean build + make check-c99 + full tests/build suite
- [ ] Rebuild code-review-graph after the batch

## Review (v0.297, complete)

Status of the 17-input stress set: **13/17 close** (was 3/17). 10 new closures.

- **G1 (single → recognized kernel):** EllipticK/EllipticE → 2F1 reduce rules;
  BesselY → BesselJ combination reduce rule; `rec_gamma_upper` for upper
  incomplete Gamma[ν,·]. Closes In9, In16, In17.
- **G2 (hyperbolic Dirichlet):** `rec_hypergeom_dirichlet` for Csch (→Zeta) and
  Sech (→LerchPhi/β). Closes In12, In13.
- **G3 (Mellin convolution engine):** new `rec_convolution` branch in
  `dispatch_term` + helpers `conv_bessel/exp_rate/gauss_a2/equal` and
  `conv_JJ/KK/JK/expJ/gaussJ`. J·J (Weber-Schafheitlin via 2F3→rec_pfq), K·K
  (Barnes first lemma), J·K equal order (Kummer), exp·J (2F1), Gaussian·J (1F1).
  Closes In1-5.

All closed forms numerically cross-checked vs NIntegrate during development;
pinned in `tests/test_integrate_ramanujan.c` (`test_mellin_single_extensions`,
`test_mellin_convolution`) — full suite green. Build clean under
`-Wall -Wextra`; `make check-c99` PASS; `make check-messages` PASS (no new
stderr). valgrind: no leak traces through any new function; residual leaks are
the pre-existing `fullsimp2`/evaluate-on-Simplify pattern (unchanged in
character, affects all Mellin cases).

**Deferred (documented in spec scope note + changelog):** In6/In7/In8 need
WhittakerW / WhittakerM / ParabolicCylinderD as new special-function builtins
(separate campaigns, per user); In14 AiryAi² is a genuine Meijer-G with 1/3-step
Gamma coefficients that does not reduce to a simple pFq/Gamma ratio.

**Known minor gap (not in user's set):** a *squared* Bessel product written as
`Power[K,2]` (e.g. `BesselK[0,x]^2`) is one kernel, not `Times[K,K]`, so
`conv_KK` does not see it (mirrors how only `BesselJ[ν,z]^2` has a dedicated
square reduce rule). The user's In2/In3 use distinct-order products and work.

Not committed/tagged (awaiting user): version.h bumped 0.296 → 0.297; tag
`v0.297` to be applied on commit per release-tagging rule.
