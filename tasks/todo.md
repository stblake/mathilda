# Reduce: multivariate Max/Min, unbounded polynomial-in-Floor, and general piecewise functions

## Goal (approved)
1. Solve multivariate `Max`/`Min`: `Reduce[Max[x,y]>2,{x,y}] -> x>2 || y>2`.
2. Solve unbounded polynomial-in-`Floor`: `Reduce[Floor[x]^2>5,x,Reals] -> x<-2 || x>=3`.
3. Same treatment for other piecewise functions: `Piecewise`, `Sign`, `UnitStep`,
   `Ramp`, `Clip`, `HeavisideTheta`, `Boole`, `UnitBox`, `IntegerPart`,
   `FractionalPart`.

## Validated design (manual splits all reduce correctly; CAD/FM soundly decline
   on residual non-polynomial heads)

### Phase 1 — unified piecewise case-split + multivariate dispatch
- `reduce_realfn.c`: one `eliminate_piecewise` driven by a `piecewise_clauses`
  table that maps each head to `{(guard_i, value_i)}` + optional default:
  - Piecewise[{{v,c}..},def]: clauses (c_i,v_i), default def(or 0)
  - Sign[u]: (u<0,-1),(u>0,1), default 0
  - UnitStep[u]: (u<0,0), default 1     Ramp[u]: (u<0,0), default u
  - Boole[c]: (c,1), default 0          UnitBox[u]: (u<-1/2,0),(u>1/2,0), default 1
  - HeavisideTheta[u]: (u<0,0),(u>0,1), NO default (u==0 excluded)
  - Clip[u]/Clip[u,{a,b}]/Clip[u,{a,b},{va,vb}]: below/above clauses, default u
  - IntegerPart[u]: (u<0,Ceiling[u]), default Floor[u]
  - FractionalPart[u]: (u<0,u-Ceiling[u]), default u-Floor[u]
  Split: `Or_i[ And[ !g_0..!g_{i-1}, g_i, stmt|_{node->v_i} ] ]` (+ default branch);
  no-default heads exclude the uncovered region. Recurses; folded into the
  fixpoint driver alongside Abs/Max-Min/Mod/IP.
- Factor `apply_selector_splits` (Abs, Max/Min, piecewise) — nv/domain-agnostic —
  used by both the univariate `reduce_realfn_preprocess` and a new
  `reduce_piecewise_preprocess` (any nv).
- `reduce.c`: for nv>=2 with a selector head, run `reduce_piecewise_preprocess`,
  force Reals, then FM/CAD.
- Add `SYM_Ramp`, `SYM_UnitBox`.

### Phase 2 — unbounded integers (fixes unbounded Floor, generalizes Reduce/Integers)
- `reduce_univar.c` `reduce_univar_integers`: keep bounded-run enumeration
  (`k==n`, preserves all pinned tests) but add UNBOUNDED-TAIL rays: sample the
  left/right tail integer beyond the extreme root; emit `x<=b1` / `x>=a2`.
- `reduce_realfn.c` `translate_ksol`: accept all six `k REL int` + Inequality, map
  via `ip_defining`. Extend `ip_defining` for Ceiling (all six); Floor already
  full; Round stays Equal-only (round-half-even -> inequalities decline, sound).
- Corpus: `dec-int-unbounded`(x>0)/`dec-int-unbounded2`(x<5) become solved
  (`x>=1`/`x<=4`).

## Verify
- [ ] Targets 1-3 return expected answers; existing reduce_tests + corpus green.
- [ ] New corpus cases (multivar max/min, unbounded floor, each piecewise head).
- [ ] check-c99, valgrind, docs + changelog + memory.

## Review
All three goals delivered. Implemented in two phases, both green.

Phase 1 — general piecewise case-split + multivariate dispatch
- `reduce_realfn.c`: `piecewise_clauses` table + `eliminate_piecewise` (Piecewise/
  Sign/UnitStep/Ramp/Clip/HeavisideTheta/Boole/UnitBox/IntegerPart/FractionalPart);
  `apply_selector_splits` factored; `reduce_piecewise_preprocess` (nv>=2);
  `reduce_stmt_has_piecewise`; piecewise heads added to `node_is_realfn`.
- `reduce.c`: nv>=2 selector dispatch (force Reals -> FM/CAD).
- `sym_names.{c,h}`: `SYM_Ramp`, `SYM_UnitBox`.

Phase 2 — unbounded integers
- `reduce_univar.c` `reduce_univar_integers`: one-sided rays for satisfied tails
  (kept bounded enumeration -> all pinned tests preserved).
- `reduce_realfn.c`: `translate_ksol` accepts all six `k REL int` + `flip_rel`;
  `ip_defining` full Ceiling relations (Round stays Equal-only, half-even).

Verification
- reduce_tests: all pass (added `test_piecewise_functions`, moved 3 integer
  inequalities decline->solved). reduce_corpus: 141/141 (+23 cases).
- solve_tests, piecewise_tests: pass. check-c99: clean. valgrind: 0 additional
  errors vs baseline; no leak frames in new code.
- Docs (`solutions-of-equations.md`), changelog (`2026-08-24.md`), memory updated.

Known sound limitations (decline, never wrong): `FractionalPart[x]<1/2` and a
`Floor` whose inner value also appears outside it (periodic/mixed sets a finite
interval list can't express); unbounded `Round` inequalities (round-half-even);
multivariate integer-part (`Floor[x]+y>2`).
