---
source: src/simp/trigsimp.c
---
**Algorithm.** `builtin_exptotrig` is the partial inverse of TrigToExp. It runs
a four-stage pipeline on the argument: (1) `ReplaceRepeated` with
`exp_to_trig_rules` — `E^(I x) :> Cos[x] + I Sin[x]`, `E^(-x) :> Cosh[x] -
Sinh[x]`, `E^x :> Cosh[x] + Sinh[x]`, plus `Log`-combination patterns that fold
back to `ArcTan`/`ArcTanh`/`ArcSin`/`ArcCsch`/... ; (2) `Together` to combine
over a common denominator; (3) `Cancel` to reduce the rational; (4) a final
`ReplaceRepeated` with `exp_to_trig_simp` that recovers reciprocal heads
(`Sin/Cos :> Tan`, `Cos^-1 :> Sec`, etc.). The trig canonicalizer is suppressed
across the whole pipeline (`trig_canon_suppress_inc`/`dec`) so the intermediate
`Sin/Cos` forms survive long enough for `Together`/`Cancel` to act.

**Data structures.** Three static rule lists (`exp_to_trig_rules`,
`exp_to_trig_simp`, plus the shared `trig_factor_*` machinery) parsed once in
`trigsimp_init`. Results are routed through the active `FactorMemo` by the public
`builtin_exptotrig` wrapper path via `trig_memo_call` semantics shared with the
other trig builtins.
