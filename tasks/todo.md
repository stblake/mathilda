# Task: Lift the pure resultant LRT into the tower-level proper part

## Goal
The single-kernel pure resultant Lazard–Rioboo–Trager log part (`rm_frac_lrt`,
prior increment) closes frac integrands with **algebraic (non-rational)
residues** (`1/(x(Log[x]^2+1)) → ArcTan[Log[x]]`).  The tower proper part
(`rm_field_ratint` log top, `rm_field_hyperexp_coupled` exp top) still used the
bounded single-constant `SolveAlways` ansatz, so a *nested* algebraic-residue log
part declined.  Lift the resultant LRT into the tower recursion.

Target (was declining; full `Integrate[]` closes it):
- `Integrate`RischMacsyma[1/(x Log[x] (Log[Log[x]]^2+1)), x]` → `ArcTan[Log[Log[x]]]`

## Key insight
The tower log part is the SAME `Integrate`TranscendentalLogPart` as the
single-kernel case, with two changes: the monomial derivation becomes the **tower**
derivation `D_tower[rad]` (`rm_tower_deriv`), and the residue-constant gate must be
free of **every** lower-field variable `{x, t_0, …, t_{L-1}}`, not just `x`.  The
resultant's content strip (over the residue variable) already removes the
multi-variable content that clearing the tower denominators introduces, so only
the final `FreeQ` gate needed generalizing to a List.

## Plan
- [x] `intrat.c`: generalize the `xgate` in `intrat_log_part_core` to accept a
      `List` of gate variables (require freeness of each); relax the
      `Integrate`TranscendentalLogPart` builtin to accept a symbol OR a List.
- [x] `integrate_risch_macsyma.c`: add `rm_field_lrt_logpart(Rr, den, T, L, x)` —
      strip the t-free content (squarefree case, dH==0), build `D_tower[rad]`,
      gate `{x, t_0..t_{L-1}}`, call the builtin; wire it as the fallback after
      the `SolveAlways` ansatz declines in `rm_field_ratint` and (gated `a==0`)
      `rm_field_hyperexp_coupled`.
- [x] Tests: log-tower algebraic-residue closures + a non-elementary decline pin.
- [x] Docs: RISCH_STATUS.md (table, §3.12, §6.1 item 2, "not yet implemented"),
      changelog 2026-07-06.md, calculus.md (frac bullet + builtin ref).

## Review (DONE 2026-07-11)
- **Closes (log tower, diff-back verified):**
  `1/(x Log[x] (Log[Log[x]]^2+1)) → ArcTan[Log[Log[x]]]`,
  `1/(x Log[x] (Log[Log[x]]^2+4))`, mixed `ArcTan`+`Log` numerators,
  cubic denominators (`Log` + `Log`/`ArcTan` mix).
- **Declines (verified):** `1/(Log[x](Log[Log[x]]^2+1))` — residues depend on `x`,
  caught by the gate list.
- **Exp tower:** the fallback is wired into `rm_field_hyperexp_coupled` and the
  `Integrate`TranscendentalLogPart` builtin handles the exp derivation directly
  (verified: `[a b, 1+b^2, b, z, 2 a b^2, {x,a}] → ArcTan[b]`), but the concrete
  targets (`E^x E^(E^x)/(E^(2 E^x)+1)`) build a depth-3 tower with `E^(2 E^x)` as a
  SEPARATE commensurate kernel rather than `(E^(E^x))^2`.  That is the
  commensurate-exponent tower-builder gap (RISCH_STATUS §6.1 item 3), upstream of
  and orthogonal to the LRT — left as the next refinement.
- Tests: `integrate_risch_macsyma_tests` + `intrat_tests` + `integrals_tests` all
  green; single-kernel LRT + rational-residue tower cases unchanged.
- Build clean (`-std=c99 -Wall -Wextra`); valgrind byte-identical to the `Sin[1.0]`
  baseline (13,440 B / 420 blocks, zero module frames) on both the closing and the
  declining paths.

## Remaining refinement (out of scope)
- Commensurate-exponent kernelization inside `rm_tower_build` (`E^(k u) → t^k`),
  which would let the exp-tower LRT fallback fire (§6.1 item 3).
- Full Bronstein SPDE for the tower Hermite numerator (repeated-pole + algebraic
  residue combined; `rm_field_lrt_logpart` handles only the squarefree case).
