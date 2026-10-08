# Reduce / CAD / QE deficiencies — seed for a future dev campaign

Surfaced while landing the `Minimize` campaign II (v0.313–v0.318). Every item
below is **verified at the REPL** and is the reason a specific `Minimize` example
stays declined or slow. Fixing these in `src/solve/` (`reduce_cad.c`,
`reduce_qe.c`, `reduce_realdiag.c`, `reduce_realfn.c`, `reduce_sys.c`) would turn
sound declines into solves with no change to `Minimize` itself.

Impact legend: which 21–40 Minimize example each blocks.

## D1 — `Reduce[expr, {vars}, Reals]` (3-arg) grinds where 2-arg fails fast
The explicit-variable-list form attempts full elimination and consumes the whole
`TimeConstraint` (often overrunning it, see D6) on a hard/high-degree system,
whereas `Reduce[expr, Reals]` returns unevaluated in milliseconds.
```
Reduce[x+y+z==1 && x>=0&&y>=0&&z>=0 && x^6+y^6+z^6-x y z < -8/243, Reals]   (* ~2ms unevaluated *)
Reduce[x+y+z==1 && x>=0&&y>=0&&z>=0 && x^6+y^6+z^6-x y z < -8/243, {x,y,z}, Reals]  (* grinds *)
```
The `Minimize` lower-bound certificate (`mz_entails`) uses the 3-arg form, so it
cannot fail fast; the compact-region shortcut had to be ordered *before* the
certificate as a workaround. **Direction:** the 3-arg form should recognise an
undecidable/too-costly instance and return unevaluated quickly (as the 2-arg form
does) rather than grind.

## D2 — multi-variable emptiness: `Reduce` gives up where `Resolve[Exists[…]]` decides
For the identical question "is this region non-empty", the direct `Reduce` leaves
a 3+-variable equality-constrained region unevaluated, but the existential-QE path
decides it.
```
Reduce[x+y+z==1 && x>=0 && y>=0 && z>=0 && x^2+y^2+z^2>4, Reals]            (* unevaluated *)
Resolve[Exists[{x,y,z}, x+y+z==1 && x>=0&&y>=0&&z>=0 && x^2+y^2+z^2>4], Reals]  (* -> False, ~80ms *)
```
`Reduce`'s direct emptiness decision is strictly weaker than its own `Exists` QE
for the same statement. **Direction:** route `Reduce[cons, vars, Reals] === False`
emptiness checks through the `Exists` QE internally (or share the machinery).
(`Minimize`'s boundedness probe already works around this by calling
`Resolve[Exists[…]]`.)

## D3 — universal QE over a constrained region returns unevaluated — THE big one
`Reduce[ForAll[{vars}, Implies[cons, f >= b]], {b}, Reals]` is the exact,
sound way to read a constrained infimum. It is decided for small/easy cases but
left **unevaluated** for genuinely positive-dimensional constrained minima:
```
(* #28: determinant over the 4-ball, true answer b <= -1/2 *)
Reduce[ForAll[{x,y,z,w}, Implies[x^2+y^2+z^2+w^2<=1, x w - y z >= b]], {b}, Reals]  (* unevaluated *)
(* #36 *)
Reduce[ForAll[{x,y,z}, Implies[x^3+y^3+z^3-3 x y z==1 && x>=0 && y>=0, z >= b]], {b}, Reals]  (* unevaluated *)
(* works for easy cases: *)
Reduce[ForAll[{x,y}, Implies[x^2+y^2<=1, x+y>=b]], {b}, Reals]              (* -> b <= -Sqrt[2] *)
```
Blocks **#28, #36** (and, via the epigraph route, #33). The `Minimize` constrained
QE scaffold (v0.318) is in place and will solve these automatically once the QE
engine can eliminate them. **Direction:** strengthen `ForAll` elimination over a
semialgebraic region for n>=3 / degree>=3.

## D4 — `Max` / `Abs` under a quantifier aborts
```
Reduce[ForAll[{x,y,z}, Max[Abs[x-y],Abs[y-z],Abs[z-x]]+x^2+y^2+z^2 >= b], {b}, Reals]  (* $Aborted ~25s *)
```
Even the epigraph polynomialisation (`Max[…]->t` with `t>=each`, `Abs[g]->s` with
`s>=±g`) produces a system the QE can't settle in time. Blocks **#33** (true min 0
at the origin). **Direction:** piecewise/`Max`/`Abs` preprocessing into the CAD, or
a cheaper epigraph-aware path.

## D5 — high-degree CAD wall even with rational coefficients
The algebraic-coefficient wall is known; separately, the CAD cannot settle a
**rational-coefficient** degree-6-in-3-vars (or degree-3 on the open simplex)
decision:
```
Reduce[x+y+z==1 && x>=0&&y>=0&&z>=0 && x^6+y^6+z^6-x y z < -8/243, Reals]   (* unevaluated (#34 cert) *)
Reduce[x+y+z==1 && x>0&&y>0&&z>0 && x y + x z + y z < 9 x y z, Reals]       (* unevaluated (#40 cert) *)
```
`Minimize` dodges #34 via the compact-region shortcut, but #40 (rational objective
with boundary poles) stays declined because its certificate is exactly this.
**Direction:** raise the CAD's practical degree/dimension ceiling, or add a
positivity/SOS certificate path for these inequalities.

## D6 — `TimeConstrained` does not interrupt a `Reduce`/`Solve` C-loop
A single CAD/Solve probe overruns its `TimeConstrained[..., t]` budget because the
inner C loop never checks the async abort flag (see the existing memory note
`project_timeconstrained_async_malloc_lock`). Consequences: `Minimize` declines
are not tightly bounded (e.g. **#26** went from a 9 ms fast-decline to a ~60 s
slow-decline once elimination exposed a harder residual), and D1's grind can't be
clamped. **Direction:** add periodic abort-flag checks at CAD/Solve recursion
points (projection, lifting, Gröbner steps) so `TimeConstrained` can preempt them.

---

### Suggested campaign order (ROI × tractability)
1. **D2** (route emptiness through `Exists` QE) — small, removes a whole class of
   spurious "unevaluated" and would let `Minimize` drop its `Resolve` workaround.
2. **D6** (abort-flag checks) — cross-cutting robustness; makes every other probe
   bounded and D1 clampable.
3. **D1** (3-arg fast-give-up) — pairs with D6.
4. **D3** (stronger constrained `ForAll` elimination) — biggest capability win
   (#28, #36), hardest.
5. **D5** (SOS / higher CAD ceiling), **D4** (`Max`/`Abs` under quantifiers) — specialised.
