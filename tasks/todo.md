# Task: RischMacsyma — remove arbitrary caps from the proper-part Hermite / flat-tower ansätze

## Goal
Extend the "no arbitrary caps" principle (already applied to both RDE solvers via
`rm_rde_var_bound`) to the remaining ansatz sites that still truncate a *derived*
degree bound with a magic constant (`if(d>N)d=N`) or gate the SolveAlways system
with a magic completeness ceiling (`nunk <= 60/80`). Every degree/Laurent bound
must be a **function of the input's degrees**, not a hardcoded constant.

## Why correct-by-construction (never wrong regardless of the bound)
Each of these ansätze is SolveAlways-certified AND the enclosing case is diff-back
verified (`rm_verify_antideriv`: `Simplify[D[result,x]-f]===0`). So a bound can only
DECLINE (too small) or waste SolveAlways time (too large) — never ship a wrong
closed form. Truncating a derived bound with a magic constant is therefore pure
incompleteness for no principled reason: a hack. Removing the truncation strictly
WIDENS the search → never regresses a currently-closing integrand.

## The capped sites (grep)
1. `rm_field_hyperexp_coupled` (~1662): `Nx = …; if(Nx>8)Nx=8` + `… <= 80` ceiling.
2. `rm_log_tower_case` (~2107): `Ntop>4→4`, lower-var `d = a+b+1; if(d>3)d=3`, `nunk<=80`.
3. `rm_exp_tower_case` (~2467): `ihi>4→4`, `ilo<-4→-4`, x-deg `d>2→2`,
   **inner exp kernels hardcoded `lo[j]=-2; bd[j]=4` (range -2..2)**, `nunk<=80`.
4. single-kernel Hermite (~2958): lower-var `d = max(a,b)+1; if(d>3)d=3`, `nunk<=60`.
5. single-kernel hyperexp (~3135): `ihi>4→4`, `ilo<-4→-4`, `d>3→3`, `nunk<=60`.

## Plan
- [ ] Add helper `rm_var_mult_at_zero(p, v)` (lowest power of v with nonzero coeff;
      factors out the inline `a`-computation already in the exp cases).
- [ ] Site 1: drop `if(Nx>8)Nx=8`; replace `<=80` ceiling with the `nunk>0`
      overflow guard (matching the `rm_field_rde` precedent).
- [ ] Site 2: drop `if(Ntop>4)Ntop=4` (Ntop=deg+1 is the exact log bound) and
      `if(d>3)d=3` (keep derived `a+b+1`); `<=80` → `nunk>0`.
- [ ] Site 3: drop `ihi/ilo` truncations (ihi=dnum-dden, ilo=-a are the exact top
      Laurent structure) and x-deg cap; **derive the inner exp-kernel Laurent window
      from f's per-kernel degree/pole structure widened by the top Laurent span**
      (a function of the input, replacing the magic -2..2); `<=80` → `nunk>0`.
- [ ] Sites 4,5: drop `d>N→N` and `ihi/ilo` truncations; `<=60` → `nunk>0`.
- [ ] Keep the benign structural clamps (`if(d<0)d=0`, `if(Ntop<0)Ntop=0`,
      `if(ihi<0)ihi=0`, `ng>=16`, `nunk>0`) — these are correctness/overflow, not
      completeness caps.

## Verification (empirical is the arbiter — diff-back protects correctness)
- [ ] `integrate_risch_macsyma_tests` green — CRUCIAL guards:
      - must STILL CLOSE: coupled exp towers `E^x E^(E^x)` … `E^x E^(E^x)/(4+E^(2E^x))`
        (lines 379-392); log towers (362-368); Hermite/hyperexp coupled cases.
      - must STILL DECLINE: non-elementary `E^(E^x)`, `E^(E^x)/(1+E^(E^x))`,
        `E^(E^x)/(1+E^(2E^x))`, `E^(2E^x)/(1+E^(E^x))` (496-505) — the diff-back
        regression guards.
- [ ] Add cap-free anchor(s): a lower-var / Laurent degree beyond the old cap.
- [ ] `integrals_tests`, `intrat_tests` green; no timing blowup.
- [ ] valgrind baseline-identical (13,440 B / 420 blocks, no module frame).

## Review

Done (2026-07-11).  All seven capped flat-ansatz / Hermite sites in
`src/calculus/integrate_risch_macsyma.c` are now cap-free — the companion to the
RDE-solver change, extending the same "no arbitrary caps" principle.

- **New shared helper** `rm_var_mult_at_zero(p, v)` (lowest power of `v` with a
  nonzero coeff = exact negative Laurent extent for an exp kernel); also factors out
  the inline `a`-computation the exp cases already had.
- **Top-kernel bounds now exact** (were truncated at 4): LOG top `Ntop = deg_top(f)+1`
  (D lowers degree); EXP top Laurent `[-mult_top(den), deg_top(num)-deg_top(den)]`
  (D preserves degree).
- **Inner exp-kernel Laurent windows derived** (was hardcoded `t_i ∈ [-2,2]`): each
  kernel's own extent widened by the reach of the top derivation coefficient
  `w' = D[u_n]` — a function of the input, the coupling behind `∫E^(x+E^x)=E^(E^x)`.
- **Lower-field polynomial proxies** keep their derived `deg_v(num)+deg_v(den)+1` /
  `max(...)+1` form with the `if(d>N)d=N` truncation dropped; the `+1`/`+2` additive
  slack is the derived-bound family (untouched), not a rejection threshold.
- **Completeness ceilings removed**: every `nunk ≤ 60/80` / `nh+ng ≤ 80` → the
  overflow-only `nunk > 0` guard (matching the RDE-solver precedent).
- **Structural boundaries kept** (not degree caps): tower depth `nl ∈ [2,4]` (deeper
  towers are the recursive engine's job), `ng < 16`, `if(d<0)d=0`, `Ntop≥0`, `ihi≥0`.

Verified:
- `integrate_risch_macsyma_tests` green (all 19 tests), incl. new cap-free anchors
  `Log[Log[x]]^{5,7}/(x Log[x])` (Ntop 6/8) and `E^x E^(6 E^x)/(1+E^(E^x))`
  (top Laurent deg 5); coupled towers `E^x E^(E^x)` … `E^x E^(E^x)/(4+E^(2E^x))` and
  the `E^(E^x)` / `E^(E^x)/(1+E^(E^x))` decline guards all hold.
- `integrals_tests`, `intrat_tests` green; no timing blowup (3.65s / 0.66s).
- valgrind baseline-identical (13,440 B / 420 blocks, **zero** module frames in any
  lost stack).
- Docs: `RISCH_STATUS.md` (§0, §3.12, §6.1 item 1), `calculus.md`,
  `docs/spec/changelog/2026-07-06.md`.

**Remaining follow-up:** only the deep leading-coefficient *cancellation/resonance*
sub-case (Bronstein's recursive degree reduction) — decline-only, rare.
