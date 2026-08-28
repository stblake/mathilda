# Reduce — Implementation Plan

*A design for implementing `Reduce` (and its quantifier/logic family: `LogicalExpand`,
`Exists`, `ForAll`, `Resolve`, `FindInstance`, `CylindricalDecomposition`) in Mathilda,
reusing the `Solve` infrastructure and adding a CAD-based real-inequality engine.*

---

## Status (as of 2026-08-26)

Phases **0–5 are implemented, tested (`tests/test_reduce.c`), and leak-clean**;
**Phase 6a–6c (two-variable CAD over the Reals) has landed** (`reduce_cad.{c,h}`,
`reduce_real_util.{c,h}`); and **Phase 6d (n-variable recursive CAD over the Reals) has landed in full**
— Stage A (`reduce_cad_nvar` + the recursive lift) and Stage B (the n-D boundary
merge that closes outer ranges for closed regions), all in `reduce_cad.c`, v0.088.
The **Phase-8 options polish has landed** (all seven `Options[Reduce]` registered
and honored, `reduce_opts.{c,h}`, 2026-08-25). **Phase 7 (quantifier elimination)
has landed in v1** (`Exists`/`ForAll`/`Resolve`, `reduce_qe.{c,h}` +
`reduce_cad_qe` seam, v0.095, 2026-08-26): the fully-quantified decision procedure
and single-free-variable parametric elimination over the rational-fibre regime.
**A zero-dimensional nonlinear-system engine (`reduce_zerodim`) has landed**
(v0.097, 2026-08-26, issue #69): when the equations pin the variety to finitely
many points, both `Reduce` and `Solve` solve them exactly and filter each branch
by the inequalities and realness via the `qqbar` oracle — covering the
irrational-fibre zero-dimensional case that CAD declines, but NOT the
positive-dimensional one (see *Known limitations*).
**Phase 8 companions are now mostly landed**: `LogicalExpand` + `NotElement`
(v0.101), `Xor`/`Implies` evaluation on Booleans (v0.102), and **`FindInstance`
over Complexes / Reals / Integers / Rationals / Booleans (v0.104, 2026-08-26;
extended v0.105, 2026-08-27 with parametric-family instantiation, a solve-the-
parameter step for periodic instances, indexed variables `c[i]`, and a bounded
integer-box search)**, and **`CylindricalDecomposition` (v0.111, 2026-08-28)** —
completing the Phase-8 companion family.
**Phase 6b (real-algebraic-coefficient fibre isolation) has landed** (v0.113,
2026-08-28, `reduce_algfiber.{c,h}` + `reduce_cad.c`): the positive-dimensional
irrational case is now solved by iterated-resultant tower projection + exact
`qqbar` filtering, so `x^2==2 && y<x`, the radius-√2 ball, the sphere positive
octant and cubic/hyperbola fibres all decompose.
The remaining pieces are **6e (McCallum well-orientedness augmentation)** and
**7-extended (≥2 free vars / alternating quantifiers / algebraic free-variable
boundaries — the algebraic discriminant-variety case still needs a nested-QE
emitter beyond 6b fibres)**.
Implementation order followed was
`0 → 1 → 2 → 3 → 5 → 4 → 6(2-var) → 6d(n-var A) → 6d(n-var B) → 8-opts → 7(v1)`.

| Phase | What | Status |
|---|---|---|
| 0 | Front-end + DNF normal-form layer | ✅ done |
| 1 | Complete univariate equations (Complexes) | ✅ done |
| 2 | Univariate real sign diagram (Reals) | ✅ done |
| 3 | Linear real systems (Fourier–Motzkin) | ✅ done |
| 4 | Parametric linear systems (Complexes) | ✅ done |
| 5 | Integers / Rationals | ✅ done |
| 6 | Multivariate nonlinear CAD (Reals) | ◧ 2-var done (6a–6c); n-var done (6d Stage A + Stage B n-D boundary merge); 6b (algebraic-coeff fibres, tower) ✅ done (v0.113); 6e (well-orientedness) pending |
| 7 | Quantifier elimination (`Exists`/`ForAll`/`Resolve`) | ◧ v1 done (fully-quantified decision procedure + single-free-var parametric QE, rational-fibre regime); ≥2 free vars / alternating / algebraic-boundary parametric deferred (blocked on 6b) |
| 8 | Companion builtins + polish | ✅ LogicalExpand + NotElement + FindInstance (C/R/Z/Q/Booleans) + CylindricalDecomposition done |
| 9 | Elementary real functions (radicals, `Abs`, `Log`, inverse-trig, `Floor`/`Mod`) over the Reals | ✅ done (+ multivariate `Sqrt` rationalization, 2026-08-24) |
| Opt | Options: `Backsubstitution`, `Cubics`, `GeneratedParameters`, `Method`, `Modulus`, `Quartics`, `WorkingPrecision` | ✅ done (2026-08-25) |

> **2026-08-28 — Phase 6b real-algebraic-coefficient fibre isolation (tower), v0.113.**
> The CAD Reals engine no longer requires a non-innermost breakpoint to be rational.
> A new primitive `rru_algebraic_fiber_roots` (`src/solve/reduce_algfiber.{c,h}`)
> isolates the real roots of a fibre polynomial with real-algebraic-number
> coefficients by **iterated-resultant tower projection** — the outer assignment is
> a tower ℚ⊆ℚ(α₀)⊆ℚ(α₀,α₁)⊆…, each irrational αᵢ a root of a known factor over ℚ
> (carried by the CAD's root provenance); substitute the rational levels, eliminate
> each algebraic level by `Resultant`, isolate the resulting rational-coefficient
> univariate's real roots, then keep exactly the candidates β with
> `factor(vals…,β)==0` decided by the exact `qqbar` oracle (discarding the
> conjugate-spurious roots the resultant introduces). Wired into both drivers
> (`reduce_cad.c`): the rational-fibre gates (former `reduce_cad.c:960-962` / `:1448`)
> are removed and a section's defining factor is threaded through
> `cad_build`/`cad_leaf`/`lift_fiber`/`reduce_cad_qe` via a new borrowed `asgdef[]`;
> the all-rational assignment keeps its unchanged `rru_collect_roots` fast path.
> **Deviation from the §4 sketch:** no primitive-element / number-field construction
> is built — the tower is projected by successive **public** `Resultant` calls
> (reusing exactly what the projection code already uses), and the conjugate-spurious
> roots are removed by an exact `qqbar` zero-test rather than by primitive-element
> re-expression. Full tower height is supported (single- and multi-extension alike);
> a generous var-degree / node-count budget declines rather than blows up. **Now
> solved** (all declined before): `x^2==2 && y<x`, `x^2+y^2+z^2<=2`, the sphere
> positive octant, `x^3+y^3==1 && x>0 && y>0`, hyperbola branches, closed 3-/4-balls
> — flowing through `CylindricalDecomposition` too. **Soundness preserved:** a
> transcendental breakpoint (`Sqrt[Pi]`), an undecidable `qqbar`, an identically-zero
> resultant, a budget overrun, or FLINT absent all decline. Tests:
> `test_cad_algebraic_fibre`; seven `alg-*` corpus rows (166/166); stale
> rational-fibre decline assertions updated. `check-c99` clean; reduce+solve
> suites/corpora green; valgrind at macOS baseline (13,440/420, zero per-call leak).
> **Still pending: 6e (well-orientedness) and 7-extended (≥2 free vars /
> alternating / algebraic discriminant-variety boundaries — needs a nested-QE
> emitter, not just 6b fibres).**
>
> **2026-08-28 — Phase 8 `CylindricalDecomposition` (last companion), v0.111.**
> `CylindricalDecomposition[expr, vars]` returns the **merged cylindrical formula** of the
> real solution set — a quantifier-free `And`/`Or` in which each variable is bounded
> cylindrically in terms of the earlier ones (`x^2+y^2<=1, {x,y}` →
> `-1<=x<=1 && -Sqrt[1-x^2]<=y<=Sqrt[1-x^2]`; `x y>1` → `(x<0&&y<1/x)||(x>0&&y>1/x)`;
> `x^2==-1` → `False`). **Deviation from this sketch:** it did NOT need a `cad_extract`
> cell-list exposure "before DNF merging". Mathematica's `CylindricalDecomposition` returns
> the *merged* cylindrical formula (verified against 14.0), which is exactly what our Reals
> engine (`reduce_fm`/`reduce_univar`/`reduce_cad`) already emits — the only semantic
> difference from `Reduce` is that CD is **Reals-only**. So it shipped as a thin front-end
> (`src/solve/reduce_companions.c`, `builtin_cylindrical_decomposition`): validate arity
> (`[expr,vars]`, or a redundant `[expr,vars,Reals]`; any other 3rd positional declines) and
> `vars` (symbol or `List` of symbols), then build+evaluate `Reduce[expr, vars, Reals,
> <trailing option Rules…>]` under message suppression, forwarding all of `Reduce`'s options,
> and decline (return NULL) iff the result is still headed by `Reduce`. Zero engine
> duplication; the whole `Reduce` front-end preprocessing (Abs sign-split, real radicals,
> `Mod`→`Floor`), DNF build and per-arity Reals dispatch are reused. Sound declines carry over
> from `Reduce` — notably the positive-dimensional irrational-fibre case (Phase 6b limit).
> New `SYM_CylindricalDecomposition`; `ATTR_PROTECTED`; docstring. Tests:
> `test_cylindrical_decomposition` (14 asserts); reduce suite + corpus 160/160; `check-c99`
> clean; valgrind at macOS baseline (13,440/420, no new leak). This completes the Phase-8
> companion family; 6b / 6e / 7-extended remain.
>
> **2026-08-27 — `FindInstance` round 2, v0.106.** Three more verification-gated
> capabilities (`src/solve/reduce_companions.c`, `src/boolean.c`): (1) **`Equivalent`
> evaluates** (`builtin_equivalent`) and the DNF engine rewrites it to a cyclic
> `And[Implies[...]]`, so `LogicalExpand` / `FindInstance[…,Booleans]` handle
> `Xor[p,q] && Implies[q,r] && Not[Equivalent[p,r]]` → `{{p→T,q→F,r→F}}`; (2) **numerical
> witness** for transcendental/inexact Real systems Reduce decides unsoundly — its
> `False` is not trusted (`fi_is_transc_inexact`), `NMinimize[{0,expr},vars]` supplies a
> verified inexact point for `0<x<0.001 && Sin[1/x]>0.999` → `{{x→0.000903}}`; (3)
> **Gröbner emptiness** (`fi_groebner_empty`) — Rabinowitsch `t·(∏dᵢ)−1` + equalities,
> basis `{1}` certifies emptiness over ℂ⊇ℝ⊇ℤ, so the 2×2 nilpotent `M^2==0 && det≠0` →
> `{}`. Guards kept: `a^5+b^5==c^5`/`p·q==prime` → `{}` (sound Reduce False), and all
> round-1 cases. suite + corpus 158/158; new code at macOS valgrind baseline.
>
> **2026-08-27 — `FindInstance` extended, v0.105.** Four verification-gated
> extensions reach witnesses the raw Reduce/Solve outputs do not surface
> (`src/solve/reduce_companions.c`): (1) **generated-parameter instantiation** — a
> parametric family `x -> ConditionalExpression[value(C[1]), C[1]>=1]` is tried over
> a small integer grid, reaching the fundamental Pell solution of `x^2-61 y^2==1` at
> `C[1]==1`; (2) **solve-the-parameter** — a single leftover parameter is solved
> against the remaining constraints over the Reals and rounded (`Ceiling` of e.g.
> `50000/Pi`), reaching the periodic `Sin[1/x]==0 && 0<x<10^-5` at `C[1]==15916`;
> (3) **indexed variables** `c[i]` matched with `expr_eq` (Reduce, which rejects
> them, is skipped as the oracle to avoid a spurious `Reduce::ivar`); (4) **bounded
> integer search** (`fi_integer_search`) over the Integers — a finite box is
> decidable (Diophantine witness, or `{}` on exhaustion / an instant linear
> reach-range `{}` for the 0/1 knapsack `Σ pᵢ cᵢ==500, 0<=cᵢ<=1`), an unbounded
> domain is best-effort (witness or sound decline, never `{}` — so `a^4+b^4+c^4==d^4`
> declines within budget). `test_find_instance` extended; suite 267 PASS; corpus
> 158/158; new code at macOS valgrind baseline.
>
> **2026-08-26 — Phase 8 `FindInstance` (Complexes / Reals / Integers /
> Rationals / Booleans), v0.104.** `FindInstance[expr, vars, dom, n]` returns up
> to `n` witness points in `Solve`'s form, `{}` when the set is provably empty,
> and unevaluated otherwise (`src/solve/reduce_companions.c`, `builtin_find_instance`;
> `SYM_FindInstance`). **Design: witness-by-verification, not CAD surgery.** The
> three CAD engines are `static` across separate `nu==1`/`nu==2`/`nu>=3` paths with
> no witness extractor; instead witnesses are read off the **public**
> already-cylindrical output of `Reduce` (satisfiability oracle: `False` → `{}`,
> formula → walk `Or` clauses, sample free-var intervals via `rru_rational_between`,
> resolve `Equal` pins via internal `Solve`) and off `Solve`'s rule-lists (free
> listed vars instantiated to 0, grid for extra instances — so it succeeds on
> `x^2 - y z == 1` where `Reduce` declines). **Every** candidate is verified
> (`expr /. point === True`), so a returned instance is always true, `{}` only when
> `Reduce` proves emptiness, and unwitnessed/unrefuted ⇒ decline. `Booleans` reuses
> the `LogicalExpand` DNF engine (clause → partial assignment, unconstrained ⇒
> `False`, verified with the now-evaluating `Xor`/`Implies`). Options: `Modulus -> p`
> honored; `Method`/`WorkingPrecision`/`RandomSeeding` accepted-and-ignored; unknown
> ⇒ `FindInstance::optx` + unevaluated. Sound declines (documented): transcendental
> numeric instances, positive-dim irrational-fibre real systems (Reduce declines),
> region constraints, bare-`!=` statements. Tests: `test_find_instance` (21 asserts);
> reduce corpus 158/158; `check-c99` clean; valgrind at macOS baseline (no
> attributable leak); module also cleared of its pre-existing compiler warnings.
> `CylindricalDecomposition` remains the last unshipped Phase-8 companion.
>
> **2026-08-26 — zero-dimensional nonlinear systems (`reduce_zerodim`), v0.097
> (issue #69).** `Reduce` and `Solve` now solve a conjunction whose polynomial
> **equations** pin the variety to finitely many points (a zero-dimensional
> ideal), optionally carrying inequality / disequation side relations — the
> three-circle system `u^2+v^2==9 && u^2+(a+v)^2==36 && (a+u)^2+v^2==25 && u>0 &&
> v>0 && a>0` and simpler cases like `x^2==1 && y^2==4`. New shared engine
> `src/solve/reduce_zerodim.{c,h}`: split each conjunct into equations `E` and
> side relations `O` (`<`,`<=`,`!=`); solve `E` over the complexes with the
> existing polynomial-system solver (a parametric / underdetermined answer ⇒
> positive-dimensional ⇒ **decline**); keep each solution branch only if every
> side relation holds there and — over the Reals — every coordinate is real,
> decided **exactly** by the FLINT `qqbar` oracle (`rru_sign_of`,
> `flint_qqbar_equal`, and the new `flint_qqbar_is_real`), an undecidable test
> ⇒ decline. Complete for zero-dimensional systems (the finite solution set is
> enumerated and filtered exactly), sound everywhere. Wired into `reduce.c`
> (Complexes after `reduce_eq_system`; Reals `nu>=2` after `reduce_fm` +
> `reduce_cad` — CAD declines on exactly the irrational fibres this covers) and
> `solve.c` (an equations-with-constraints pre-pass mirroring the Integers
> pre-pass, gated on the presence of a side constraint so pure-equation `Solve`
> is untouched). This is the **zero-dimensional** complement to Phase 6b — it
> lifts the rational-fibre restriction *only* when the equations already reduce
> the variety to points, so no CAD algebraic-fibre isolation is needed. The
> **positive-dimensional** case is unchanged and still declines: a nonlinear
> system over the Reals with irrational fibres and a *free* dimension — e.g.
> `x^2+y^2+z^2==1 && x>0 && y>0 && z>0`, whose solution set is a 2-D surface
> region — is the n-variable CAD's rational-sample limitation (Phase 6b), a
> separate enhancement (see *Known limitations* below). Tests:
> `zdim-*` corpus rows + `test_reduce.c` updates; reduce corpus 158/158, solve
> corpus 99/99; `check-c99` clean; no new leaks.
>
> **2026-08-26 — Phase 7 quantifier elimination (v1), v0.095.** `Exists`, `ForAll`
> and `Resolve` landed. `Exists`/`ForAll` are inert (`HoldAll`) quantifier wrappers;
> `Resolve[expr, dom]` and a top-level `Exists`/`ForAll` inside `Reduce[...]` both
> route to one engine over the Reals (front-end `reduce_qe.{c,h}`; new symbols
> `SYM_Exists`/`SYM_ForAll`/`SYM_Resolve`; peel in `builtin_reduce`; `reduce_qe_init`
> from `reduce_init`). The engine reuses the existing CAD: `cad_build` already
> computes the fold (`cell->empty` is `!Exists` over the bound fibre, `all_true` is
> the `ForAll` roll-up), so QE is a per-cell read. Three regimes: **fully quantified**
> → a real-closed-field decision procedure (`Exists = Reduce[φ,B,Reals] =!= False`,
> `ForAll = ... === True`, reusing the whole engine); **one free variable** → a new
> public seam `reduce_cad_qe` (`reduce_cad.{c,h}`) that builds the CAD with the free
> variable outermost, labels its cells by the quantifier verdict over the bound
> fibre, and emits via `rru_emit_sign_diagram` (`Exists[y,x^2+y^2<1]→-1<x<1`,
> `ForAll[y,x^2+y^2>=1]→x<=-1||x>=1`, `Exists[x,x^2==a]→a>=0`); **everything else**
> (≥2 free vars, alternating prefix, non-`Reals` domain, non-rational free-variable
> breakpoint — so the two-parameter `Resolve[Exists[x, x^2+bx+c==0]]` — or any
> undecidable sign) declines, preserving soundness. Same-kind nested quantifiers
> flatten; 3-arg `Exists[x,c,g]→c&&g`, `ForAll[x,c,g]→!c||g`; `nbound==0` strips;
> a listed-but-absent variable is left unconstrained (not a decline). Tests:
> `test_quantifiers_{decision,parametric,decline}` in `tests/test_reduce.c`;
> `check-c99`/`check-packed-aware` clean; no QE-attributable leak. **Deviation from
> the original Phase-7 sketch:** `ForAll` is done DIRECTLY via the `all_true`
> roll-up rather than as `Not[Exists[Not φ]]` (no `rform_not_*` needed — those
> helpers never existed), and the fully-quantified case reuses `builtin_reduce`
> wholesale instead of a bespoke CAD fold.
>
> **2026-08-25 — Reduce options (all seven registered and honored).**
> `Options[Reduce]` now reports `{Backsubstitution -> False, Cubics -> False,
> GeneratedParameters -> C, Method -> Automatic, Modulus -> 0, Quartics -> False,
> WorkingPrecision -> Infinity}`, and each is honored. A Solve-style trailing-option
> peeler (`is_reduce_option_name`/`reduce_apply_option`, unknown name → `Reduce::optx`,
> unevaluated) fills a `ReduceOpts` (`src/solve/reduce_opts.{c,h}`) threaded through the
> engines. The four options overlapping `Solve` reuse its machinery: **Cubics/Quartics**
> are forwarded onto the internal `Solve[...]` calls (radicals vs `Root[]`), via a shared
> `reduce_opts_build_solve`; **Modulus** routes the equational statement through Solve's
> `solvemod` residue enumeration as a top-level pre-pass (`reduce_modular`), reformatted
> as an `Or` of `x == r`; **GeneratedParameters** renames the Diophantine engine's `C[k]`
> to `h[k]` (`rename_param_head`). **Backsubstitution** is accept/validate/echo (the
> current linear engine's grafted output is the `False` default and also what `True`
> requests — no fork); **WorkingPrecision** threads a numeric-fallback tolerance into the
> transcendental sign decisions (`Infinity` keeps the exact path); **Method** is reserved
> (`Automatic` only). New `SYM_Backsubstitution`; defaults registered in
> `options_builtin.c`. Tests: nine `test_option_*` groups in `tests/test_reduce.c`. Two
> scoped exclusions, documented not silent: CAD fibre isolation keeps `Root[]` regardless
> of Cubics/Quartics (radical fibres over the multivariate-real path are deferred), and
> the residual radical-path leak under `Cubics/Quartics -> True` is pre-existing in
> `solvepoly.c` (reached identically via `Solve`), not introduced here.

> **2026-08-24 — radicals compose with `Abs`; univariate domain-gate soundness fix.**
> A Phase-9 preprocessing extension rationalizes square-root radicals for the
> **multivariate** path too: `Sqrt[u] REL c` is squared under sign guards into a
> polynomial `And`/`Or` in the same variables, so `Sqrt[Abs[x]]+Abs[y]<1` and
> `Sqrt[x^2+y^2]<1` (previously unevaluated for any radical in ≥2 variables) now
> solve on the existing FM/CAD. Separately, the univariate general sign diagram's
> domain gate is now **per-conjunct** (was global), fixing a wrong `Sqrt[Abs[x]]<1
> -> x==0` / `Log[Abs[x]]<0 -> False`. See `docs/spec/builtins/solutions-of-equations.md`.
>
> This rationalize-by-square approach was chosen over the general **aux-variable
> purification** method (rewrite each non-polynomial atom to a fresh variable with a
> defining relation, then project it out) because that method needs **two features
> still pending here**: **Phase 7** quantifier elimination (`reduce_qe.c` — the
> existential fold over CAD cells, unbuilt) to drop the auxiliary variable from the
> answer, **and** **Phase 6b** algebraic-coefficient CAD fibres (lifting the
> rational-breakpoint restriction) so the *real* variables — pushed to non-innermost
> CAD levels once the auxiliary variable is innermost — may carry irrational
> boundaries. Both remain on the roadmap; purification becomes viable only once they
> land. Rationalize-by-square keeps the problem at its original dimension and needs
> neither.

The **[Deviations from this plan](#deviations-discovered-during-implementation)**
section at the end records where the built code differs from the original design
(soundness fix, extra `RAtom` fields, actual file split, known limitations for the
Phase-8 polish pass). Everything shipped honours the hard invariant: an undecidable
sign/ordering makes `Reduce` return unevaluated, never a wrong formula.

---

## Context — why this, and what "done" means

Mathilda has a mature **`Solve`** subsystem (`src/solve/`) but no `Reduce`. `Solve` returns
the *generic* solution of **equations** as a list of rules (`{{x -> a}, ...}`), silently
dropping the parametric/degenerate cases and offering **no inequality support over the
reals**. `Reduce` is the strictly more powerful primitive: it returns a **complete,
quantifier-free logical description** of the solution set of **equations *and*
inequalities** over a domain, including every degenerate case.

```
Reduce[a x == b, x]              (a != 0 && x == b/a) || (a == 0 && b == 0)
Reduce[x^2 > 1, x, Reals]        x < -1 || x > 1
Reduce[x^2 + y^2 <= 1, {x,y}, Reals]
                                 -1 <= x <= 1 && -Sqrt[1-x^2] <= y <= Sqrt[1-x^2]
Reduce[x^2 == 4, x, Integers]    x == -2 || x == 2
```

The strategic goal (chosen scope: **full CAD roadmap**, **full companion family**):
reuse the Solve infrastructure wherever it already answers a sub-question, and add only
the genuinely-missing machinery — a **logical normal-form layer**, a **real-inequality
engine** culminating in **Cylindrical Algebraic Decomposition (CAD)**, and **quantifier
elimination**. The companion builtins that Reduce's machinery enables — `LogicalExpand`,
`Exists`, `ForAll`, `Resolve`, `FindInstance`, `CylindricalDecomposition` — are included
as later phases.

"Done" for the whole program = the phased roadmap in §7 lands, each phase independently
verified. "Done" for *this* task = `REDUCE_PLAN.md` is written.

**Hard invariant, everywhere:** an undecidable sign/ordering (the real-algebraic oracle
returning "unknown", or FLINT being absent) must make `Reduce` **return unevaluated
(`NULL`)** — never a wrong formula. Soundness over completeness.

---

## What `Reduce` must do (semantics)

- Signatures: `Reduce[expr, vars]`, `Reduce[expr, vars, dom]`,
  `dom ∈ {Complexes (default), Reals, Integers, Rationals}`.
- `expr` is a logical combination (`&&`, `||`, `!`, `Implies`, `Xor`, chained
  `Inequality`) of atoms: equations (`==`, `!=`) and inequalities (`<`, `<=`, `>`, `>=`),
  possibly with free parameters.
- Output: `True`/`False`, or an `And`/`Or`/`Not` tree of relational atoms
  (`==`,`!=`,`<`,`<=`,`>`,`>=`,`Element`) describing the **complete** solution set.
- Over Complexes only equations (and `!=`) are meaningful (no ordering); over Reals the
  full inequality machinery engages; over Integers/Rationals the Diophantine engine drives.
- Quantifiers `Exists[{y..}, φ]` / `ForAll[{y..}, φ]` and `Resolve[...]` eliminate the
  quantified variables (later phase).

---

## Reuse inventory (the point of the exercise)

Mathilda already ships almost every hard *primitive*; what is missing is the *orchestration*.

| Capability | Reuse target | Location |
|---|---|---|
| **Front-end template** (option peel, var validation, dom arg, True/False short-circuit) | `builtin_solve` | `src/solve/solve.c:817` |
| **Univariate equation solving** → roots, emits `Root[]` for irreducibles | `solvepoly_solve_polynomial_equality(eqn, var, dom, opts)` | `src/poly/solvepoly.h:41` |
| **Linear / nonlinear systems** (Gauss-Jordan; lex Gröbner + back-sub) | `solvelinsys_solve_linear_system`, `solvenlsys_solve_nonlinear_system` | `src/solve/solvelinsys.c`, `solvenlsys.c` |
| **Integer/Diophantine engine — already handles inequalities & orderings** | `solveint_solve_integer`, `SICtx`, `classify_conjunct`, `register_inequality`, `flatten_conjuncts`, `build_result`, HNF linear | `src/solve/solveint.c:106`, `solve_common.c`, `solveint_internal.h` |
| **CAD projection primitives** | `flint_polynomial_resultant`, `resultant_subresultant` (Bronstein PRS), `builtin_discriminant`, `Subresultants`, `SubresultantPolynomials` | `src/poly/flint_bridge.c`, `poly.c:4303,4527,4751` |
| **Elimination (lex / elimination order)** | `gb_buchberger`, `GroebnerBasis` (`GB_ORDER_ELIM`), `Eliminate` | `src/poly/groebner*.c`, `eliminate.c` |
| **Algebraic numbers + rigorous real sign oracle** | `Root[Function[t,p],k]` (real-first ascending index); `root_numericalize`; **`flint_qqbar_compare(a,b)` → sign(a−b) ∈ {−1,0,1}, −2=undecided — defined but currently UNUSED**; `flint_qqbar_equal`, `flint_qqbar_is_constant_algebraic` | `src/root.{c,h}`, `src/root_numeric.h`, `src/poly/flint_qqbar.{c,h}` |
| **Scalar sign/zero decision** | `decide_pair(op,a,b)` → 1/0/−1; `compare_numeric`; `zero_test_decide`, `zero_test_decide_assuming` + `AssumeCtx` | `src/comparisons.c:448`, `src/simp/simp_assume.c` |
| **Multivariate poly container** | `MPoly`, `expr_to_mpoly`/`mpoly_to_expr`, `mpoly_deg_var`, `mpoly_coef_of_var`, `mpoly_lc_var`, `mpoly_subst_var_mpz`, `mpoly_total_deg` | `src/poly/mpoly.h` |
| **Real-root count (Sturm certificate)** | `sturm_real_root_count` (static; expose if needed) | `src/root_numeric.c:184` |
| **Parametric case-split precedent** | `SolveAlways` = coefficient-vanishing + recursive `Solve` | `src/solve/solvealways.c` |

**Two structural gaps that force new code** (flagged so scope is honest):
1. **No logical normal-form layer.** `src/boolean.c` only collapses `And`/`Or`/`Not` on
   literal `True`/`False`. No DNF/CNF, `LogicalExpand`, `BooleanConvert`, `Resolve`,
   `Exists`, `ForAll`.
2. **No CAD / Fourier–Motzkin / sign-diagram engine.** `flint_qqbar_compare` is the only
   real-algebraic sign machinery present, and nothing calls it yet.

---

## Architecture

### Module layout (all new, under `src/solve/`)

```
reduce.h / reduce.c        [done] public builtin_reduce + reduce_init; front-end
                           (mirrors builtin_solve): dom, True/False, var validation,
                           dispatch routing
reduce_form.h / .c         [done] INTERNAL NORMAL FORM: RAtom / RConj / RForm, builders,
                           simplifiers, DNF ops, emit-to-Expr
reduce_atom.c              [done] atom canonicalisation (poly REL 0), sign-normalise,
                           denominator-clearing detection, form_from_expr parser
reduce_eq.{c,h}            [done] Phase 1: complete univariate equation solver (Complexes)
reduce_univar.{c,h}        [done] Phase 2: univariate real sign diagram; Phase 5's bounded
                           integer enumeration (reduce_univar_integers)
reduce_fm.{c,h}            [done] Phase 3: Fourier–Motzkin for linear real (in)equalities
reduce_int.{c,h}           [done] Phase 5: Integers/Rationals adapter over Solve
reduce_sys.{c,h}           [done] Phase 4: parametric linear systems (symbolic Gauss)
reduce_real_util.{c,h}     [done] shared real-algebraic primitives (qqbar sign oracle,
                           rational-sample selection, Solve-based real-root isolation)
reduce_cad.{c,h}           [done] Phase 6: McCallum projection + partial-CAD lifting;
                           2-var driver (reduce_cad) + n-var recursive engine
                           (reduce_cad_nvar) with the n-D boundary merge (Stage A+B)
reduce_zerodim.{c,h}       [done] zero-dimensional nonlinear systems (Complexes & Reals),
                           shared by Reduce and Solve: solve the equations exactly, filter
                           branches by side relations + realness via the qqbar oracle
                           (issue #69). Covers the irrational-fibre ZERO-dim case CAD
                           declines; positive-dim stays with 6b.
reduce_qe.{c,h}            [pending] Phase 7: Exists / ForAll / Resolve via CAD cells
reduce_companions.{c,h}    [done] Phase 8: LogicalExpand, NotElement,
                           FindInstance (witness-by-verification off Reduce/Solve
                           outputs; C/R/Z/Q/Booleans), and CylindricalDecomposition
                           (Reals-only front-end delegating to Reduce)
reduce_realfn.{c,h}        [done] Phase 9: preprocessing (Abs sign-split, Mod->Floor,
                           integer-part isolation) + the head->real-domain table
reduce_realdiag.{c,h}      [done] Phase 9: general univariate real sign diagram over
                           radical / pole / bounded-domain-transcendental atoms
```

Wire `reduce_init()` into `src/core.c` next to `solve_init()` (~line 836); add
`SYM_Reduce`, `SYM_Resolve`, `SYM_Exists`, `SYM_ForAll`, `SYM_FindInstance`,
`SYM_CylindricalDecomposition`, `SYM_LogicalExpand` to `src/sym_names.{h,c}` (3 sites
each: `extern` decl, `= NULL` def, `intern_symbol(...)`). `src/solve` is already on the
Makefile `-I` path, so no build-graph change is needed. Attributes: `ATTR_PROTECTED`
(like `Solve` — Reduce does **not** hold its args; `Exists`/`ForAll` get `ATTR_HOLDALL`
so their bound-variable lists are not evaluated).

### Internal normal form — **DNF** (`reduce_form.h`)

Every engine below *produces* a disjunction of cell-conjunctions, and real `Reduce`
output is overwhelmingly a top-level `Or` of `And`s — so DNF is the natural target.

```c
/* As built: >/>= are canonicalised away (operands swapped), so only these
 * five relations ever reach an engine. */
typedef enum { R_EQ, R_NE, R_LT, R_LE, R_ELEM } RRel;

typedef struct {            /* one atom: (poly) REL 0, canonical */
    Expr*  poly;            /* owned; LHS-RHS moved here, RHS forced to 0 */
    RRel   rel;
    Expr*  elem_dom;        /* for R_ELEM */
    Expr*  display;         /* [added] solved-form override, e.g. Equal[x, b/a] */
    bool   nonconst_denom;  /* [added] canonicalisation cleared a variable denom */
    int    main_var, deg_main;   /* cheap cached classification */
    bool   is_linear;            /* total degree <= 1 in all reduce vars */
} RAtom;

typedef struct { RAtom* a; int n, cap; bool is_false; } RConj;  /* an And */
typedef struct { RConj** c; int n, cap; bool is_true; } RForm;  /* an Or  */
```

`RConj.is_false` and `RForm.is_true` are absorbing sentinels (mirroring
`builtin_and`/`builtin_or`). Key helpers:

```c
RAtom  reduce_atom_canonicalize(Expr* rel, Expr** vars, int nv);
void   rconj_push(RConj*, RAtom);              /* dedup + contradiction detect */
bool   ratom_eq(const RAtom*, const RAtom*);
RForm* rform_or(RForm*, RForm*);               /* concat conjunction lists */
RForm* rform_and(RForm*, RForm*);              /* distributive product */
RForm* rform_not_atom(const RAtom*);           /* De Morgan on one atom */
void   rform_simplify(RForm*);                 /* drop false, absorb, dedup, subsume */
RForm* rinterval_emit(int var, const RInterval*, Expr** vars);
Expr*  rform_to_expr(const RForm*, Expr** vars, int nv);
```

**Atom canonicalisation** (reuse existing simplifiers — do not hand-roll poly arithmetic):
1. Head → `RRel`; flip `>`/`>=` to `<`/`<=` by swapping sides (halves the case matrix).
2. `poly = Numerator[Together[lhs - rhs]]` via the existing `Together`/`Numerator`
   builtins (exactly how `solvealways.c` clears denominators). **Denominator caveat:** for
   strict/non-strict relations the cleared denominator can flip the sense — emit an
   auxiliary `den != 0` atom rather than assume positivity; for `EQ`/`NE` the numerator
   alone is sound.
3. Sign-normalise (leading coeff of lex-lowest var positive) so equal atoms are
   structurally identical for dedup.
4. Fill `main_var`/`deg_main`/`is_linear` from `expr_to_mpoly` + `mpoly_deg_var`.

**Simplification inside a conjunction** (`rconj_push`): constant atoms decided by
`decide_pair` / `flint_qqbar_compare` (True→drop, False→`is_false`, unknown→keep as a
parametric condition); same-poly clashes (`p==0` ∧ `p!=0` → false; keep the stronger of
`<`/`<=`); **redundant one-variable bounds** (`x>0 && x>1 → x>1`) collapse into an
`RInterval` per variable, with bound comparison via `flint_qqbar_compare` so
`Root[]`/radical bounds order correctly. `rform_simplify` drops false conjunctions,
promotes an empty conjunction to `is_true`, dedups and does a cheap subsumption pass (no
full BooleanMinimize — Mathematica's own output is not minimal either).

**Emission** builds `Equal/Unequal/Less/LessEqual[poly, 0]` (or `Root`-bounded interval
atoms), then a single `evaluate()` lets the existing `And`/`Or` flatten and optionally
fuses `r1 < x && x < r2` via the existing `Inequality` collapser.

### Dispatch routing (`reduce.c`, mirrors `builtin_solve`)

After the Solve-style front-end (option peel, `is_valid_solve_vars`, dom at `arg[2]`,
`True`→`True`/`False`→`False`), peel any top-level `Exists`/`ForAll` to `reduce_qe`, build
the `RForm` from the input (De Morgan `Implies`→`!a||b`, expand `Xor`, split `Inequality`
chains), then route **cheaply** off the cached atom flags:

```
route(RForm F, vars V, dom D):
  D == Integers | Rationals            -> reduce_integer(F,V,D)     # solveint wrapper
  all atoms EQ/NE (no ordering):
      D == Complexes (default)         -> reduce_equational(F,V)    # complete eq engine
  every atom is_linear && D == Reals   -> reduce_fm(F,V)            # Fourier-Motzkin
  |V| == 1 && D == Reals               -> reduce_univar(F, V[0])    # sign diagram
  D == Reals                           -> reduce_cad(F,V)           # full CAD
  else                                 -> NULL (unevaluated)
```

Ordering matters for **efficiency**: linear-real → Fourier–Motzkin **before** CAD (avoids
exponential lifting for the common `a x + b y <= c` case); univariate-real → sign-diagram
before CAD (skips projection bookkeeping). Detection is O(#atoms), no re-parsing.

---

## The engines

### 1. Complete equation solver, parametric case split (`reduce_eq.c`)

Turns Solve's *generic* rules into Reduce's *complete* tree via recursive leading-
coefficient-vanishing splits (precedent: `SolveAlways`):

```
solve_case(poly p in x):
  lc = mpoly_lc_var(p, x)
  if lc is a nonzero constant:                       # generic terminal
      solvepoly_solve_polynomial_equality(p==0, x, dom, opts)
      -> Or_k (x == r_k)                             # Root[]/radical rules reformatted
  else:
      A = (lc != 0)  &&  solutions_of_degree_n(p)    # generic branch
      B = (lc == 0)  &&  solve_case(p without leading term)   # degenerate branch
      return A || B
  base (p constant a_0):  a_0 == 0 ? True : (a_0 == 0)
```

`Reduce[a x == b, x]` → `(a!=0 && x==b/a) || (a==0 && b==0)`. Over Reals, `Root`
solutions are kept only if real (`sturm_real_root_count`/qqbar). `p != 0` is the De Morgan
complement of `p == 0`. **Systems:** linear via `solvelinsys_solve_linear_system` with a
pivot-by-pivot `d!=0` vs `d==0` split (multivariate analogue of the lc-split); nonlinear
via `solvenlsys_solve_nonlinear_system` for the generic branch — **full comprehensive-
Gröbner case analysis is deferred**; genuinely-nonlinear parametric systems route through
CAD (Reals) or emit the generic branch guarded by the non-vanishing leading coeffs
(sound-on-generic, matching Solve today). Over **Complexes this equation route is the
entire engine** — the common `Reduce[poly system, vars]` path, shipped first.

### 2. Univariate real sign diagram (`reduce_univar.c`)

`p_i(x) REL 0` atoms over one real variable:
1. **Breakpoints:** real roots of each `p_i` via `solvepoly_...(p_i==0, x, Reals)` →
   rationals/radicals/`Root[]`.
2. **Order & dedup** with `flint_qqbar_compare` (works directly on `Root[]`); any `−2`
   (undecided / FLINT off) → **abort to NULL**.
3. **Cells:** `2m+1` alternating open intervals and point cells across `b_1<...<b_m`.
4. **Sample points:** breakpoint for a point cell; a **rational strictly between**
   neighbours for an interval (numericalise ends with `root_numericalize`, pick a rational,
   *certify* with `flint_qqbar_compare`, bump precision on collision); a rational beyond
   the ends for the unbounded cells. Evaluate each atom's sign at the sample via
   `flint_qqbar_compare(p_i(sample), 0)`; combine per `F`'s Boolean tree → cell truth.
5. **Emit:** merge maximal runs of true cells with correct endpoint openness
   (`<` vs `<=`), isolated true points → `x == b_j`, `!=` holes split runs, unbounded true
   ends → `x < b_1` / `x > b_m`. `Reduce[x^2>1,x,Reals]` → `x<-1 || x>1`.

### 3. Fourier–Motzkin for linear real systems (`reduce_fm.c`)

For all-`is_linear` atoms over Reals: eliminate variables one at a time by pairing each
positive and negative bound to generate implied constraints; a residual `c < 0` with `c`
constant → `False`. New code, but small and self-contained; keeps the linear class off the
CAD path entirely.

### 4. Multivariate CAD over Reals (`reduce_cad.c`) — the large phase

**Projection (McCallum), eliminating the *last* listed variable first** (so the output
formula reads in the given variable order):
- `disc_{x_j}(p)` for each `p` (reuse `builtin_discriminant` /
  `flint_polynomial_resultant(p, ∂p/∂x_j, x_j)`),
- leading coeff `mpoly_lc_var(p, x_j)` (degree-drop locus),
- pairwise `Res_{x_j}(p, q)` (reuse `flint_polynomial_resultant`, or `resultant_subresultant`
  when FLINT is off),
- **factor every projection poly** and keep distinct irreducible factors (biggest
  constant-factor win, controls subresultant blow-up).

**Lifting** recursively from the univariate base upward: substitute the fixed sample into
each poly (`mpoly_subst_var_mpz`/`ReplaceAll`), isolate & order the fiber's real roots
(the §2 machinery) → sub-cells; recurse; at the top level evaluate the input formula's
truth via `flint_qqbar_compare(p(sample),0)`.

**Partial CAD (Collins–Hong):** before lifting a cell's children, evaluate the part of
`F` that depends only on already-fixed vars; if truth is already forced, **skip the
subtree**. This is what makes `x^2+y^2<=1` tractable.

**Cell → DNF extraction:** each true leaf contributes one conjunction of section
(`x_i == Root[...]`) / sector (`Root[..lo..] < x_i < Root[..hi..]`) atoms in projection
order; merge adjacent true sectors/sections (the §2 endpoint logic, one level up) so
`x^2+y^2<=1` emerges as `-1<=x<=1 && -Sqrt[1-x^2]<=y<=Sqrt[1-x^2]`.

**Well-orientedness caveat (McCallum):** if a projection poly *nullifies* on a lower cell,
soundness isn't guaranteed — detect it, add the nullified poly's coefficients (Brown/Hong
augmentation) and re-lift; if still nullified, **bail to NULL**. The genuinely-missing
sub-capability is *real-root isolation of a univariate poly with real-algebraic-number
coefficients* (a fiber at a non-rational sample): prefer **rational sample points**
wherever the stack allows (fiber then has rational coeffs and `solvepoly`+qqbar apply
directly); otherwise numericalise high-precision, isolate, and certify count/order with
`sturm_real_root_count` + `flint_qqbar_compare`, bailing to NULL on undecided.

### 5. Quantifier elimination (`reduce_qe.c`)

`Exists[{y..}, φ]` / `ForAll` / `Resolve` run CAD with **quantified vars projected first**
(free vars outermost), then fold truth over each free-variable cell's children:
`Exists` = OR over child sub-cells, `ForAll` = AND (equivalently `!Exists !φ`, reusing
`rform_not_*`). A fully-quantified sentence (`X` empty) returns `True`/`False` — a decision
procedure. Because partial CAD already computes truth bottom-up, QE is just the fold; no
extra projection machinery.

### 6. Companion builtins (`reduce_companions.c`)

- **`LogicalExpand`** — expose the internal DNF distributor (`rform_and`/`rform_or` +
  `rform_simplify`) as a builtin; needed internally anyway.
- **`FindInstance[expr, vars, dom]`** — one witness: return the sample point of any true
  CAD cell (or one Solve solution for the equational case); `{}` when unsatisfiable.
- **`CylindricalDecomposition[expr, vars]`** — *(shipped v0.111 as a Reals-only front-end
  delegating to `Reduce`, NOT as a raw `cad_extract` cell dump — see the Deviations note.)*
  The merged cylindrical formula is what our Reals engine already emits and what Mathematica
  returns, so it forces the Reals domain and reuses the whole `Reduce` pipeline.

### 9. Elementary real functions over the Reals (`reduce_realfn.c`, `reduce_realdiag.c`)

A whole class of ordinary problems — `Reduce[Sqrt[x+3-4Sqrt[x-1]]+Sqrt[x+8-6Sqrt[x-1]]==1, x]`
→ `5 <= x <= 10`, `Reduce[Abs[Abs[x]-2]+Abs[Abs[x]-5]==5, x]`, `Reduce[Floor[2x-1]==3, x]`,
`Reduce[ArcSin[x]+ArcCos[x]==Pi/2, x]` — used to bubble back unevaluated because every
Phase-0 atom is a **polynomial** (`Numerator[Together[lhs-rhs]]`) and `reduce_univar`'s
`collect_breakpoints` declines the moment it sees a non-polynomial. Phase 9 is a **general
univariate real sign diagram** that consumes these, reusing the existing 1-D cell
decomposition and its emission (`rru_emit_sign_diagram`, shared with `reduce_univar`).

The insight: over the reals such a statement's truth is **piecewise-constant**, changing
only at a finite set of breakpoints. Three pieces supply what the polynomial engine lacks:

- **Preprocessing (`reduce_realfn_preprocess`)** — an Expr→Expr rewrite run in `reduce.c`
  *before* `reduce_form_from_expr`: (1) `Mod[u,m]` → `u - m*Floor[u/m]`; (2) a relational
  leaf linear in a single `Floor`/`Ceiling`/`Round` is expanded to its defining
  inequalities (`Floor[v]==n` → `n<=v<n+1`, etc.); (3) every `Abs[u]` is eliminated by
  **sign-splitting** into `Or[And[u>=0, …/.Abs[u]->u], And[u<0, …/.Abs[u]->-u]]` (so `Solve`
  is never asked to invert `Abs`, which it cannot). Detecting an elementary real function of
  the variable also **forces the domain to Reals**, so a pure radical *equation* with no
  ordering routes here instead of the Complexes equational path.

- **The head→real-domain table (`reduce_real_domain_collect`)** — the one place that
  enumerates supported partial-domain functions: even radical `u^(p/q)` → `u>=0`; `Log[u]` →
  `u>0`; `ArcSin`/`ArcCos[u]` → `-1<=u<=1`; `ArcTanh`/`ArcCosh`/`ArcSech` likewise; a
  rational pole → `denom != 0`. Each contributes a per-sample **domain-gate** constraint and
  a **boundary** whose roots are breakpoints.

- **The sign diagram (`reduce_univar_general`)** — breakpoints are the soft `Solve[·==0, x,
  Reals]` roots of every atom (descending polynomial factor structure so an isolated root of
  a squared factor survives even when Solve returns the identity `{{}}`), the poles, and the
  domain boundaries; sorted with an exact (qqbar) compare that **falls back to a numeric N
  sign** for a transcendental breakpoint (a multiple of `Pi`). At one sample per cell the
  truth oracle applies the domain gate (out-of-domain ⇒ excluded — essential: over ℝ,
  `ArcSin[2]+ArcCos[2]==Pi/2` and `Sqrt[x^2-4]==Sqrt[x-2]Sqrt[x+2]` at x=0 both read `True`
  in ℂ but are out of the real domain), the pole gate, then decides equations/`Unequal` by
  `evaluate` (so a radical/transcendental identity resolves via `PossibleZeroQ`) and
  orderings by an exact-then-numeric sign.

A supporting soundness fix lives in `reduce_atom.c`: `canonical_poly`/`canonical_denom`
**skip `Together`** when the difference contains a branch-cut transcendental (`Log`, inverse
trig/hyperbolic), because `Together[Log[x^2]-2Log[-x]]` collapses to the constant `-2 I Pi`
— wrong for real `x<0`, where the difference is `0`.

Soundness invariant preserved throughout: a free parameter, an undecidable sign, or an
unsupported domain node makes the engine return NULL and `Reduce` stays unevaluated. One
documented deviation: at a **removable 0/0 singularity** the engine reports the sound open
boundary (`x>0` for `x/Sqrt[x^2]+Sqrt[x^2]/x==2`), which can differ from Mathematica at that
single point. Verified by `tests/test_reduce.c` (`test_real_functions`) and the sampling
corpus (`tests/reduce_corpus.m`, `rf-*` records).

---

## 7. Phased roadmap (ship order **0→1→2→5→3→4→6→7→8**)

Phases 0–5 cover the bulk of real-world `Reduce` calls before the CAD investment. Each
phase is independently verified against `expr_to_string_fullform` in a new
`tests/test_reduce.c` (mirror `tests/test_solve.c`; add a `reduce_tests` target **with an
`add_test(...)` line** — note `solve_tests` lacks one — and register in
`tests/CMakeLists.txt`). Optionally add a `tests/reduce_corpus.m` + verifier mirroring
`solve_corpus.m`.

| Phase | Deliverable | New/modified files | Reuse | Verify (input → output) |
|---|---|---|---|---|
| **0** Front-end + normal form | `builtin_reduce` (Solve-style parse, True/False, bad-var); `RForm`/`RAtom`, `rform_to_expr`, atom canon; constant-atom eval | `reduce.{c,h}`, `reduce_form.{c,h}`, `reduce_atom.c`; `core.c`; `sym_names.{h,c}` | `builtin_solve`, `decide_pair`, `Together`/`Numerator` | `Reduce[True,x]→True`; `Reduce[1<2,x]→True`; `Reduce[x==x,x]→True`; `Reduce[3<2,x]→False` |
| **1** Complete univariate equations (Complexes) | `reduce_eq_univariate` lc-vanishing split | `reduce_eq.c` | `solvepoly_...`, `mpoly_lc_var`, `SolveAlways` pattern | `Reduce[a x==b,x]→(a!=0&&x==b/a)||(a==0&&b==0)`; `Reduce[x^2==4,x]→x==-2||x==2` |
| **2** Univariate real sign diagram | `reduce_univar`, `RInterval`+`rinterval_emit` | `reduce_univar.c`, `reduce_form.c` | `solvepoly` (Reals), `flint_qqbar_compare`, `root_numericalize` | `x^2>1→x<-1\|\|x>1`; `x^2>=1→x<=-1\|\|x>=1`; `x^2<1→-1<x<1`; `(x-1)(x-2)(x-3)>0→1<x<2\|\|x>3`; `x^2!=1→x!=-1&&x!=1` |
| **3** Linear real systems (Fourier–Motzkin) | `reduce_fm` | `reduce_fm.c` | `RAtom.is_linear`, rational arith | `x+y<1&&x>0&&y>0,{x,y}→x>0&&0<y<1-x`; `x>1&&x<0→False` |
| **4** Parametric linear systems *(built as `reduce_sys.c`, symbolic Gaussian elimination — see Deviations)* | `reduce_eq_system` (const/symbolic pivot split, p==0 substitute-recurse, back-sub) | `reduce_sys.{c,h}` | `Coefficient`/`Together`/`Solve` | `{a x+y==1, x+y==0},{x,y}` → `1-a!=0 && x==1/(a-1) && y==1/(1-a)`; `a x==1 && x==2` → `2a-1==0 && x==2` |
| **5** Integers / Rationals *(built as `reduce_int.c` + bounded enumeration)* | reformat `Solve[..,dom]` output to `\|\|` of `==`/`Element`; univariate inequality → sign-diagram enumeration | `reduce_int.{c,h}`, `reduce_univar_integers` | `Solve[..,Integers\|Rationals]`, `collect_breakpoints` | `Reduce[x^2==4,x,Integers]→x==-2\|\|x==2`; `Reduce[x^2<10&&x>0,x,Integers]→x==1\|\|x==2\|\|x==3` |
| **6** Multivariate CAD (Reals) — land incrementally: 6a 2-var full CAD rational fibers · 6b Root-coeff fibers + certification · 6c partial-CAD pruning · 6d n-var recursion · 6e well-orientedness + augment/bail | `reduce_cad` | `reduce_cad.c` | `flint_polynomial_resultant`/`resultant_subresultant`, `builtin_discriminant`, `Subresultants`, `solvepoly`, `flint_qqbar_compare`, `MPoly` | `x^2+y^2<=1,{x,y}→-1<=x<=1 && -Sqrt[1-x^2]<=y<=Sqrt[1-x^2]`; `x^2+y^2<0→False`; `x y>0→(x>0&&y>0)\|\|(x<0&&y<0)` |
| **7** Quantifier elimination | `Exists`/`ForAll`/`Resolve` on CAD truth tree | `reduce_qe.c`; register symbols | `reduce_cad` | `Resolve[Exists[x, x^2+b x+c==0], Reals]→b^2-4c>=0`; `Reduce[Exists[y, x^2+y^2<1],{x},Reals]→-1<x<1` |
| **8** Companions + polish | `LogicalExpand`, `FindInstance`, `CylindricalDecomposition`; `Element` I/O; `Inequality` fusion; `Assumptions` bridge via `zero_test_decide_assuming`/`AssumeCtx` to prune impossible parametric branches | `reduce_companions.c`, `reduce.c` | DNF ops, CAD cells, `AssumeCtx` | `FindInstance[x^2+y^2<1,{x,y},Reals]` → one witness; `CylindricalDecomposition[...]` cell list |

**Per-phase project hygiene (every phase):** register the builtin + `ATTR_PROTECTED`
(+`ATTR_HOLDALL` for `Exists`/`ForAll`) + `symtab_set_docstring` in `reduce_init`; update
`docs/spec/builtins/solutions-of-equations.md`; add a section to the current week's
`docs/spec/changelog/<Monday-of-ISO-week>.md`. Run `make check-c99` (POSIX/`M_*`/`int64_t`
guards) and `valgrind` on the new corpus.

---

## Constraints (non-negotiable)

- **C99 strictly** (`-std=c99 -Wall -Wextra`): `#ifndef` fallback for any `M_*` constant;
  feature-test macro *before the first include* for any POSIX function; `int64_t` uses the
  `ci_*_i64` checked-int family, never the `long long` one. `make check-c99` gates this.
- **Memory/ownership:** each builtin **borrows `res`**, returns a fresh `Expr*` (evaluator
  frees `res`) or `NULL` (unevaluated, `res` retained). NULL-out reused sub-nodes before
  the wrapper is freed. Every new `RForm`/`RAtom`/`Expr` owned is paired with a free;
  valgrind-clean.
- **FLINT/qqbar guarding:** every `flint_qqbar_*` / FLINT resultant call under
  `#ifdef USE_FLINT`; a `−2`/absent oracle degrades to **NULL (unevaluated)** — never a
  wrong answer.
- **Numeric fast-path surfaces:** `Reduce` is a symbolic/structural head returning logical
  formulas (not machine arrays), so it is legitimately exempt from the packed/NDArray and
  `Compile[]` surfaces — record the exemption with a one-line reason in the audit tool's
  `EXEMPT` list rather than leaving a silent omission.

## Where reuse is insufficient (new code, eyes open)

The logical DNF layer (`reduce_form.*`), Fourier–Motzkin, the sign-diagram cell
construction/emission, the parametric case-split recursion, the CAD driver (projection
orchestration, lifting recursion, cell tree, partial-CAD pruning, section/sector
extraction, well-orientedness handling), and real-root isolation of univariate polys with
real-algebraic-number coefficients (mitigated by preferring rational samples + Sturm/qqbar
certification, bailing to NULL on undecided).

## Critical files to mirror / reuse

- `src/solve/solve.c:817` — front-end template.
- `src/poly/flint_qqbar.{c,h}` (`flint_qqbar_compare`, `.c:580`) — real-algebraic sign oracle.
- `src/poly/solvepoly.c` — univariate real-root isolation (sign-diagram + CAD fibers).
- `src/poly/flint_bridge.c` (`flint_polynomial_resultant`) / `src/poly/poly.c:4303`
  (`resultant_subresultant`), `poly.c:4751` (`builtin_discriminant`) — CAD projection.
- `src/solve/solveint.c:106` (`solveint_solve_integer`) + `solveint_internal.h` (`SICtx`) — Integers.
- `src/solve/solvealways.c` — parametric case-split precedent.
- `src/boolean.c`, `src/comparisons.c:448` (`decide_pair`), `src/poly/mpoly.h`, `src/simp/simp_assume.c`.
- `tests/test_solve.c` + `tests/CMakeLists.txt` — test pattern (add the `add_test` line).

## End-to-end verification

1. `tests/test_reduce.c` — per-phase FullForm assertions from the tables above.
2. `tests/reduce_corpus.m` + back-substitution/containment verifier (mirror
   `solve_corpus.m`) — spot-checks: every reported cell's sample point satisfies the input;
   negation of the output is unsatisfiable on a random real/integer sample grid.
3. Differential check against `Solve` where they overlap (equations, generic branch).
4. `make check-c99`; `valgrind --leak-check=full` on the corpus; docs + changelog updated.

---

## Deviations discovered during implementation

Recorded so the remaining phases (6–8) build on what actually exists, not the
original sketch. Phases 0–5 are otherwise as designed.

- **Soundness fix — rational-function inequalities (the "Denominator caveat",
  §"Atom canonicalisation" step 2).** Clearing a *variable* denominator via
  `Numerator[Together[...]]` drops the pole and can flip an inequality's sense, so
  `1/x < 1` reduced (wrongly) to `x > 1`. As built, `reduce_atom_canonicalize`
  detects a non-constant cleared denominator and sets a new `RAtom.nonconst_denom`
  flag; such an **inequality** is neither constant-decided nor accepted by the real
  engines — Reduce **declines** (stays unevaluated). Equations, where clearing a
  fully-cancelled denominator is sound, are unaffected (`1/x == 0 → False`). The
  plan's alternative (emit an auxiliary `den != 0` atom) was not taken; decline is
  the sound minimum and a candidate for a later, more complete treatment.
- **Extra `RAtom` fields.** `display` (a solved-form emission override, so a
  solution prints `x == b/a` rather than `x - b/a == 0`) and `nonconst_denom` above.
- **`RRel` has no `R_GT`/`R_GE`.** `>`/`>=` are canonicalised to `<`/`<=` by swapping
  operands, halving the case matrix. `RInterval` was not needed — the sign diagram
  emits segments locally.
- **Phase 4 is its own file, `reduce_sys.c`, not an extension of `reduce_eq.c`.**
  It does direct **symbolic Gaussian elimination** over rational-function coefficient
  vectors (nonzero-const pivot direct; symbolic pivot `p` → `p != 0` branch +
  `p == 0` branch solved via `Solve`, substituted, recursed), rather than wrapping
  `solvelinsys`. An LSol "DNF of cases" intermediate carries per-branch conditions +
  assignments; back-substitution ("graft") expresses each variable in the parameters.
- **Phase 5 gained a bounded-enumeration fallback.** `Solve[..., Integers]` only
  engages when an equation is present, so a pure bounded inequality
  (`x^2 < 10 && x > 0`) is enumerated over the Phase-2 sign-diagram breakpoints
  (`reduce_univar_integers`, sharing the factored-out `collect_breakpoints`).
- **Implementation-order note / bug.** A double-free in the Phase-4 elimination loop
  (passing an owned `Expr*` into `expr_new_function`, which takes ownership, then
  freeing it again) corrupted the heap and surfaced as runaway evaluator recursion —
  see the `expr_new_function`-consumes-args memory.

- **Phase 6 landed as two-variable CAD (`reduce_cad.c`) + a shared-primitive
  refactor.** The exact sign oracle, rational-sample selection and Solve-based
  real-root isolation were extracted out of `reduce_univar.c` into
  `reduce_real_util.{c,h}` (with a provenance-carrying `rru_collect_roots`), so the
  CAD and the univariate sign diagram share one implementation. The single
  structural correction from the design review: the *atom* polynomials are factored
  into a distinct-irreducible squarefree basis **before** projection (not just the
  projection polynomials), which makes `disc`/`res` never identically zero and
  confines nullification to 0-dimensional sections — so the "bail on interval
  nullification" rule is provably non-restrictive at nv==2. The dispatcher prunes
  variables absent from every atom, so an effective-1-var call delegates to the sign
  diagram and only genuinely-2-var work runs CAD; nv≥3 declines.
- **Phase-6 v1 boundaries (all sound-declines).** An irrational base breakpoint is
  declined (real-algebraic-coefficient fibre isolation is deferred to the later
  numericalize+Sturm+qqbar fallback); a fibre whose `Solve` does not return a clean,
  1:1-matchable, orderable branch list is declined; the McCallum well-orientedness
  augmentation (6e) is not implemented — an interval nullification bails.
- **Boundary-merge pass (implemented).** Emission builds a structured per-cell
  y-region (`YRegion` of symbolic-bound segments) rather than an Expr, and a merge
  pass fuses a run of consecutive same-template interval cells across a breakpoint
  when the template's limit there equals that section's own fibre — closing the
  x-range for a closed region (`x^2+y^2<=1 -> -1<=x<=1 && ...`) while leaving strict
  regions open and preserving genuine holes (`x^2+y^2<=1 && x!=0` splits at 0). It is
  a cosmetic post-pass: region-equality is decided by sampling both regions, and any
  undecidable comparison leaves the already-correct unmerged form. Non-absorbed
  breakpoints with a non-empty fibre are emitted as standalone `x==b && …` sections.
- **Phase 6d landed as a HYBRID (Stage A), n-var path separate from the 2-var
  driver.** Rather than the unified recursion in §4, the shipped engine keeps the
  2-variable driver `reduce_cad` (with its boundary-merge) byte-for-byte untouched
  and adds a separate recursive engine `reduce_cad_nvar` for `nu>=3`; the dispatcher
  gate now routes `nu==1 → reduce_univar`, `nu==2 → reduce_cad` (unchanged),
  `nu>=3 → reduce_cad_nvar`. The recursive engine is an iterated McCallum projection
  stack (`PolySet pstack[]`, `cad_project_out`) plus a recursive lift
  (`cad_recurse` over the outer levels, bottoming out at `cad_leaf` — a
  parameter-generalized `lift_fiber`), with per-level partial-CAD pruning
  (`cell_dead_n`) and symbolic bounds via `symbolic_branch_lvl`. **v1 scope**: the
  rational-fibre regime — a breakpoint at any non-innermost level must be rational
  (given the rational assignment above it), else decline (Phase 6b, deferred);
  interval nullification bails (6e). **Emission (Stage A)** is a flat DNF of true
  cells with the innermost dimension merged into a `YRegion`: STRICT inequalities
  read as one clean nested conjunction, while CLOSED regions emit a correct but
  verbose union of cells (the boundary sphere/arc/pole cells listed separately).
  **Bugfix during 6d:** the nullification bail must be guarded by
  `contains_symbol(factor, decomposition_var)` — a lower-variable factor vanishing
  at its own section is not a McCallum fibre nullification (without the guard
  `x y z > 0` wrongly declined). Tests: `test_cad_nvar` + ten `cad3-*`/`cad4-*`
  corpus rows (form-agnostic oracle), leak-clean, `check-c99` clean.
- **Phase 6d Stage B — n-D boundary merge (landed, v0.088).** The recursive
  emission was restructured to build a cell TREE (`CADRegion`/`CADCell`,
  `cad_build`) and merge over it (`cad_region_expr`), generalizing the 2-var
  boundary-merge so a CLOSED region closes its outer ranges
  (`x^2+y^2+z^2<=1 -> -1<=x<=1 && -Sqrt[1-x^2]<=y<=… && …`) instead of listing the
  boundary cells; the sphere surface, closed half-ball and the 4-var closed ball
  merge likewise, while strict regions stay open. The 2-var `templates_equal`/
  `breakpoint_absorbable` generalize to `cad_templates_equal` (structural) and
  `cad_absorbable`, which decides cell-equality by **sampling**: it emits both
  cells' sub-formulas and requires them to agree on a grid drawn from the cell
  structure (`cad_sample_cell` walks the tree, `formula_truth_at` evaluates). This
  sampling formulation sidesteps the n-D boundary-degeneracy that a cell-lookup
  comparison hits (coincident breakpoints when an outer value is substituted).
  Cosmetic post-pass: any undecidable comparison leaves the already-correct
  unmerged (verbose) form, and the corpus sample-point oracle certifies the merged
  output equivalent to the input. The flat-DNF `cad_recurse` emission of Stage A
  was replaced by this tree; `xyz>0` consequently reads as a nested `Or` factored
  by the sign of the outer variable rather than a flat 4-octant list (both correct).

### Known limitations for the Phase-8 polish pass

- **Default-domain inequalities.** `Reduce[x^2 > 1, x]` (no explicit `Reals`) stays
  unevaluated; Mathematica reads a bare inequality as real. Route inequality-bearing
  default-domain input to Reals.
- **Unbounded integer sets.** `Reduce[x > 0, x, Integers]` declines; should emit an
  `x >= 1`-style form.
- **Parametric-condition display.** Conditions print as produced (`1 - a != 0`,
  `2 a - 1 == 0`) rather than Mathematica's `a != 1`, `a == 1/2` — correct but not
  minimal; a normalisation pass would tidy them.
- **Comprehensive Gröbner case analysis** for *non-linear* parametric systems remains
  deferred (declines today); Phase 6 CAD covers the Reals case.
- **Positive-dimensional nonlinear systems over the Reals with irrational fibres.**
  A conjunction whose nonlinear constraints leave a *free* dimension and whose fibre
  boundaries are irrational algebraic numbers — e.g.
  `Reduce[x^2+y^2+z^2==1 && x>0 && y>0 && z>0, {x,y,z}, Reals]`, a 2-D surface region —
  still declines. This is the n-variable CAD's **rational-sample limitation** (Phase 6b:
  real-algebraic-coefficient fibre isolation), *not* a regression and *not* covered by the
  zero-dimensional engine above (which fires only when the equations pin the variety to
  finitely many points). Closing it is the Phase 6b enhancement: lift the rational-fibre
  restriction by isolating univariate roots at a real-algebraic sample point.
