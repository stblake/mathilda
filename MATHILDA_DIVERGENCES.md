# Mathilda vs Mathematica: divergences found while porting the non-torsion certificate to `ParallelMixed.m` (2026-09-22)

Every item below was reproduced with the binary at the repo root (`./Mathilda -file probe.m`, build 0.169)
against Mathematica 13.2 running the same Wolfram-language code (`ParallelMixed.wl`, the research copy).
Items are grouped by kind; each gives a one-line repro, what Mathematica returns, what Mathilda returns,
the impact on `src/internal/mixed/ParallelMixed.m`, and the workaround now in the `.m` (search the
tag in the source), so that the workaround can be removed once the core is fixed.

Priority order for a follow-up session: **A1, A2, A4** (silent wrong results on the hot path),
then A3, A5, A6, A7, then B and C.

---

## STATUS (reconciled 2026-09-27, build v0.213)

The inline "Mathilda:" lines below are **pre-fix snapshots from build 0.169** and are kept as the
historical record. Verified against the live binary, most of section A is now fixed:

| Item | State | Fixed by / note |
|---|---|---|
| A1 `Do`/`Table` over a packed list | **FIXED** | commit `3c12c301` (v0.170); iterator materialises NDArray bounds |
| A2 `PolynomialExtendedGCD[…, Modulus->p]` | **FIXED** | `3c12c301` (FLINT `nmod_poly_xgcd`) |
| A3 `PolynomialQuotient/Remainder/QR[…, Modulus->p]` | **FIXED** | `3c12c301` |
| A4 `PadRight[raggedList]` | **FIXED** | `3c12c301` (`pr_is_atomic`) |
| A5 `Lookup[assoc, Key[k], default]` | **FIXED** | `3c12c301` |
| A6 `PolynomialMod` rational coefficients | **FIXED** | `3c12c301` (`rational_mod_int`) |
| A7 `Discriminant` of a degree-1 polynomial | **FIXED** | `3c12c301` |
| A8 `Transpose` of ragged equal-length rows | **FIXED** | `3c12c301` (top-two-levels fallback) |
| A9 `OptionValue` under `OptionsPattern[other]` | **FIXED** (common case) | resolves against `other` |
| A10 nested/stored `Function` closures | **FIXED** (common cases) | |
| A11 `Module`/pattern-var captures an arg-value symbol | **FIXED** | v0.212 — capture-avoiding `replace_bindings` |
| A12 `RootSum` restricted / non-numeric | **FIXED** | v0.214/v0.215 — Rothstein–Trager reduction, general linear denominator, numeric `N[RootSum]`; all five documented cases evaluate |
| A13 `Series` of a large radical quotient | **FIXED** | evaluates on the live binary |
| A14 `PolynomialGCD[…, Extension]` / `ToNumberField[Root, gen]` | **FIXED** | `ToNumberField` half v0.213 (A14a); `PolynomialGCD`/`Extension` Root generator v0.217 (A14b) |
| A15a `{} . {}` segfault | **FIXED** | v0.210 |
| A15b `Coefficient[…, x, i]` symbolic exponent | **FIXED** | v0.211 |
| A16 `PolynomialGCD[p, e, Extension -> Automatic]`, `e` an unexpanded algebraic constant equal to 0 | **OPEN** | returns 1 (found on the v0.216 corpus re-run, 2026-09-27 late); `.m` expands and reduces before every extension gcd |
| A17 `RowReduce[m, Method -> "OneStepRowReduction"]` with parametric entries | **OPEN** | does not finish on a 20 x 15 system linear in two parameters; `.m` uses the default method when parameters are present |
| A18 a `Do` iterator inside a package captures a caller's same-named symbol | **FIXED** (v0.238) | and it was never only `Do` -- `Sum` and `Product` too; only `Table` was registered as a scoping construct. `is_iterator_scope_head` (`src/modular.c`) now covers all four, so the `.m` renames `ma`/`ent`/`ne` are no longer needed. See F3. |
| A19 `FreeQ[e, Complex]` is True on a complex atom | **OPEN** | found 2026-09-28 porting logrewrite; `logrewrite.m` tests `_Complex` everywhere |
| A20 `ComplexExpand` writes a real nested radical as `Cos[Arg[...]]`, `Arg[1 - Sqrt[5]]` unevaluated | **OPEN** | `logrewrite.m` puts constants in rectangular form by rules (`RectPow`); `RRad` of the package inherits the hazard |
| A21 `CountRoots` not implemented | **OPEN** | `SturmCount` in `logrewrite.m` |
| A22 `$InputFileName`, `DirectoryName` not implemented | **OPEN** | a package cannot `Get` a sibling file; `LoadModule["mixed/logrewrite.m"]` instead |
| A23 `NumericQ[Root[...]]` is False; `ToRadicals` gives the Ferrari form for every quartic | **OPEN** | `RootRadicals` in `logrewrite.m` (biquadratic / palindromic forms, root picked numerically) |
| A24 `Can`'s field detour returned a `Dot[{}, Inverse[{}], {}]` coefficient | **OPEN** (state-dependent, no standalone repro) | `logrewrite.m` canonicalises with `CanRaw` |
| A18 (re-checked on v0.221) | **FIXED** (v0.238) | `g[v_] := Module[{s = 0}, Do[s += v, {k, 2}]; s]; g[k]` gave 3, now gives `2 k` |
| B7 an outer `TimeConstrained` cannot interrupt an inner one | **OPEN** (behavioural) | the package budget cannot cut the rewrite's own `TimeConstrained` short |
| B1 `ToNumberField` non-canonical primitive element | **OPEN** (behavioural, by design) | `.m` reads whatever theta comes back |
| B2–B6 | behavioural; see each entry | mostly by-design / hard |

Items A19-A24 and B7 were found on 2026-09-28 (builds 0.221/0.222) while porting the real form of
the logarithmic part (`src/internal/mixed/logrewrite.m`); each has a one-line repro below and a
workaround in that module (section D).

Remaining core work: the three items found when the review corpus was re-run
on v0.216 (2026-09-27, late): **A16** (extension gcd against an unexpanded zero constant, the cause of
three FALSE non-elementary certificates on the raw run), **A17** (`OneStepRowReduction` on a parametric
matrix); **A18** (`Do` iterators capture a caller's symbol; `Table` does not) is FIXED in v0.238. Each has a one-line
repro below and a `.m` workaround (section D). The **A14 `PolynomialGCD` Root generator** half (A14b)
is now **FIXED** in v0.217 — `qa_resolve_extension` recognises a `Root[]` object as an algebraic
generator through the radical-oriented autodetect/tower/extension-gcd pipeline (a `GEN_ROOT`
autodetect kind; reducible Root defining polynomials are declined to the gcd over Q). Everything
else in section A is resolved.

## A. Correctness: wrong or missing results

### A1. `Do[body, {i, list}]` runs zero iterations when `list` is a packed array held in a variable or produced by an expression

```
ps = Select[Range[3, 12], PrimeQ];
Do[Print[p], {p, ps}]                              (* prints nothing *)
Do[Print[p], {p, Select[Range[3, 12], PrimeQ]}]    (* prints nothing *)
Do[Print[p], {p, Prime[Range[2, 5]]}]              (* prints nothing *)
f[l_] := Do[Print[p], {p, l}]; f[ps]               (* prints nothing *)
Do[Print[p], {p, {3, 5, 7}}]                       (* works: literal list *)
Do[Print[p], {p, Range[3, 5]}]                     (* works *)
Do[Print[p], {p, Evaluate[ps]}]                    (* works *)
Do[Print[e], {e, FactorList[x^2 - 1]}]             (* works: not packed *)
```
Mathematica iterates in every case. The failing forms are exactly the NDArray-backed (packed
machine-integer) lists that reach `Do` through a symbol or an unevaluated expression. `Table` with the
same iterator is affected in the other direction: `Table[g[d], {d, Divisors[175]}]` binds `d` to the whole
list (one iteration) and then hangs inside `g` (observed with a `While` on a non-boolean test).

Impact: `NontorsionCertificate` (the Cohen-type non-torsion certificate of Proposition 9.4 for
`[oo+ - oo-]`) iterated over `{p, Select[Range[3, pmax], PrimeQ]}` and therefore never examined a single
prime in Mathilda: the certificate was dead and every such integral was reported "failed" instead of
"not elementary". Any `Do` over a computed integer list is at risk.

Workaround in the `.m`: index loops `ps = Select[...]; Do[p = ps[[ip]]; ..., {ip, Length[ps]}]` in
`NontorsionCertificate` and `NontorsionDivisorCertificate` (tag: "an index loop").

### A2. `PolynomialExtendedGCD[a, b, x, Modulus -> p]` ignores the modulus when `a` and `b` are coprime over Q

```
PolynomialExtendedGCD[x^2 + 5 x + 1, x^2 + 6 x, x, Modulus -> 7]
  Mathematica: {6 + x, {6, 1}}                      (* both vanish at x = 1 mod 7 *)
  Mathilda:    {1, {1 + 6/35 x, 174/35 - 216/35 x}} (* the gcd over Q, rational Bezout coefficients *)
PolynomialExtendedGCD[x + 6, 7 + 7 x, x, Modulus -> 7]   -> {6 + x, {1, 0}}   (* correct: b is 0 mod 7 *)
PolynomialExtendedGCD[x + 3, x + 3, x, Modulus -> 7]     -> {3 + x, {0, 1}}   (* correct *)
```
Only inputs that share a factor over Q are reduced correctly. Impact: Cantor's algorithm in the Jacobian
over GF(p) hit this whenever two divisors shared support mod p; the rational garbage then propagated as
unevaluated `Mod[... PowerMod[35, -1, 7] ...]` expressions (a 12 s "slow path" that was really a wrong
one). Workaround: `PolyExtGCDModP` (a Euclidean extended gcd on reduced polynomials) replaces the builtin
in `CantorAdd`.

### A3. `PolynomialQuotient` / `PolynomialRemainder` / `PolynomialQuotientRemainder` ignore `Modulus -> p` (return unevaluated)

```
PolynomialQuotient[x^3 + 2, x + 1, x, Modulus -> 7]   -> unevaluated   (Mathematica: 1 + 6 x + x^2)
PolynomialRemainder[x^3 + 2, x + 1, x, Modulus -> 7]  -> unevaluated   (Mathematica: 1)
```
Impact: `UnitDegreeModP` (the continued fraction of Sqrt[q] over GF(p)) never terminated a step and returned
`None`, which with A1 kept the unit certificate dead; the Jacobian arithmetic needs the same divisions.
Workaround: `ReduceModP` / `PolyQuoModP` / `PolyRemModP` (division over Q, then coefficient-wise reduction
with `PowerMod`; valid whenever the divisor's leading coefficient is a unit mod p, which holds throughout).

### A4. `PadRight[raggedList]` (one-argument form) does not pad to a rectangular array

```
PadRight[{{1}, {0, 1/8}}]
  Mathematica: {{1, 0}, {0, 1/8}}
  Mathilda:    {{{1, 0}, {0, 0}}, {{0, 0}, 1/8}}
PadRight[{1}, 2]   -> {1, 0}                      (* the two-argument form is fine *)
```
Impact: `QCoords` (the coordinates of algebraic constants over Q, used by the Q-basis decomposition of a
residue divisor in `TorsionRealise`) returned nonsense whenever a rational and an algebraic residue met,
so that decomposition -- and with it the torsion realisation of Algorithm 3(d) beyond the single-place
case -- never worked in Mathilda. Workaround: `PadRows` (`PadRight[#, Max[Length /@ rows]] & /@ rows`).

### A5. `Lookup[assoc, Key[k], default]` returns `default` for a present key

```
a = <||>; a[{1, 0}] = 2;
Lookup[a, Key[{1, 0}], None]   -> None      (Mathematica: 2)
Lookup[a, {{1, 0}}, None]      -> {2}       (* list form works *)
a[{1, 0}]                      -> 2         (* direct part works *)
KeyExistsQ[a, {1, 0}]          -> True
```
Workaround: `If[KeyExistsQ[assoc, k], assoc[k], None]` in the baby-step table of `DivisorOrderModP`.

### A6. `PolynomialMod[poly, p]` leaves rational coefficients unreduced

```
PolynomialMod[1/2 + x, 7]   -> 1/2 + x      (Mathematica: 4 + x)
```
Workaround: never feed rationals to `PolynomialMod`; `ReduceModP` handles them.

### A7. `Discriminant` of a degree-1 polynomial is 0

```
Discriminant[z, z]        -> 0   (Mathematica: 1)
Discriminant[2 z + 3, z]  -> 0   (Mathematica: 1)
```
Harmless in the `.m` (the value is filtered by `DeleteCases[..., 0]`), but a wrong value.

### A8. `Transpose` of a list of equal-length rows whose entries are lists of different shapes is left unevaluated
(Found earlier, 2026-09-22 morning, commit 5440fb16.)
```
Transpose[{{{a, b, c}, {d, e, f}}, {t1, t2}}]   -> unevaluated   (Mathematica: {{{a, b, c}, t1}, {{d, e, f}, t2}})
```
Mathematica transposes the top two levels only. Impact: a false non-elementarity certificate (`1/(x (Log[x]^2+1))`).
Workaround: `Thread[{rows...}]` at every such site (RealisePoints, the curve and m >= 3 paths, TorsionRealise,
NontorsionDivisor).

### A9. `OptionValue["name"]` inside a function declared with `OptionsPattern[other]` resolves an unpassed option's default against the enclosing function's own `Options`, not against `other`
(Found earlier; see the comment above `Options[iPIM] = Options[ParallelIntegrateMixed]`.) An unpassed option
leaked as an unevaluated `OptionValue[...]` into results. Workaround: declare `Options[iPIM]`.

### A10. `Function` does not close over an enclosing `Function`'s parameter in every form
(Found earlier; see `test_method_split_specials` in `tests/test_parallelmixedtower.c`.) The simple case
`Function[x, Function[y, x + y]][1][2]` gives 3 today; the failing shape was a pure function built inside
another and stored for later use (`splittable` stayed False, SplitSpecials never ran). Workaround: the
inner predicate was rewritten without the nested closure.

### A11. A `Module` local or a pattern variable captures a user symbol of the same name inside an argument value

```
g[v_] := Module[{e = 1}, v + e];   g[e + 1]                 (* 3;  Mathematica: 2 + e *)
Integrate`ParallelMixedTower[Log[x + a], x]                 (* {"failed", "no solution within bounds", {1, 2}} *)
Integrate`ParallelMixedTower[Log[x + alpha], x]             (* -x + alpha Log[alpha + x] + x Log[alpha + x] *)
Integrate`ParallelMixedTower[Log[x + b], x]                 (* $RecursionLimit::reclim, then a segmentation fault *)
Integrate`ParallelMixedTower[Log[e + x], x]                 (* the tower prints e$73 + x: the symbol was renamed *)
```
A symbol that reaches a function inside a value (the integrand's parameter `a`, `b`, `e`, ...) is
renamed or bound together with a same-named `Module` local (`a`, `b` in `VInf`/`ResInf`/`RealiseClass`,
`e` in `LaurentPolyTimes`, ...), so every integrand with a free parameter whose name is a local
somewhere in the package fails: `1/(x^2 + a^2)`, `x Exp[a x]`, `x^a`, `1/Sqrt[x^2 + a]`, `Log[x + a]`
("no solution", "verification failed"), `Exp[a x] Sin[b x]` (crash). The same integrands with the
parameter named `alpha` all solve. Mathematica's `Module` renames only the symbols of the body it
holds, never those inside values. Review corpus (`review_corpus.py` of the research directory,
family "parameters"): R348-R352, R356, R366 are this defect.

### A12. `RootSum` evaluates only a restricted rational form, and not numerically

```
RootSum[Function[z, z^4 - 3 z + 3], Function[z, ((z^2 + 1)/(4 z^3 - 3))/(x - z)]]      (* (1 + x^2)/(3 - 3 x + x^4): ok *)
RootSum[Function[z, z^3 - 2], Function[z, (-284/625 + 571/11250 z + 1516/5625 z^2)/(x - z)]]
                                              (* unevaluated (the numerator's rational constants are merged into the denominator) *)
RootSum[Function[z, z^3 + z + 1], Function[z, (a z + 1)/(x - z)]]                      (* unevaluated: a parameter *)
RootSum[Function[z, z^3 + z + 1], Function[z, z (x + z^2)/(x + z)]]                    (* unevaluated *)
N[RootSum[Function[z, z^3 + z + 1], Function[z, z Log[x - z]]] /. x -> 1/3, 30]        (* unevaluated *)
pf = Function[z, z^4 - 3 z + 3]; RootSum[pf, Function[z, (1 + a z)/(x - z)]]            (* RootSum[pf, ...]: the symbol pf is not evaluated *)
```
Mathematica evaluates every rational form symbolically and `N` numerically. `D[RootSum[f, T(z) Log[x - z]], x]`
does evaluate (to the rational form), which the verify gate relies on; the `.m` therefore never asks the
kernel to evaluate a RootSum: the residual is reduced by `A/p`, `A = T p' mod p` (`RootSumLogand`), the
RootSum is built with `Apply` and appears only in the answer, with the coefficient in the form `A(z)/p'(z)`.

### A13. `Series` leaves a large radical quotient unevaluated (and `SeriesCoefficient` returns 0 for it)

```
b = 2 u^18/(-1 + u^4 - 2 u^10 + 2 u^14 - u^20 + u^24);
Series[(2 uu Sqrt[1 + uu^10])/(1 - uu^4), {uu, 0, 2}]                                 (* SeriesData: ok *)
Series[(2 uu Sqrt[1 + uu^10])/(1 - uu^4 + 2 uu^10 - 2 uu^14 + uu^20 - uu^24), {uu, 0, 2}]   (* returned unchanged (Head Times) *)
SeriesCoefficient[(2 uu Sqrt[1 + uu^10])/(1 - uu^4 + 2 uu^10 - 2 uu^14 + uu^20 - uu^24), {uu, 0, 1}]   (* 0; the coefficient is 2 *)
```
`Normal` of the unevaluated call is the input and `Coefficient`/`Exponent` of that are meaningless, so a
valuation or a residue at infinity computed through `Series` was silently wrong (the residues of the genus-4
certificate came out as `{0, 0}` instead of `{-2, 2}`). The `.m` computes every expansion at a place over
infinity by polynomial arithmetic (`LaurentPolyTimes`, `InfLaurent`: truncated inverse modulo `t^(n+1)`,
binomial series for `qsn^alpha`), and no longer calls `Series` there.

### A14. `PolynomialGCD[p, q, Extension -> Automatic]` with a `Root` object returns 1; `ToNumberField[a, theta]` with a `Root` theta stays unevaluated  — **FIXED** (A14a v0.213, A14b v0.217)

```
p = x^5 - x + 1; tau0 = 256/2869 - 625 x/2869 - 500 x^2/2869 - 400 x^3/2869 - 320 x^4/2869;
c = Root[-1 + 15 #1 - 80 #1^2 + 160 #1^3 + 2869 #1^5 &, 1];
PolynomialGCD[p, tau0 - c, Extension -> Automatic]         (* was 1; now the degree-1 gcd over Q(c), Mathematica's answer up to the unit 625 *)
ToNumberField[Root[-1 - #1 + #1^5 &, 1], c]               (* was unevaluated; now evaluates (A14a, v0.213) *)
```
Both halves are fixed. **A14a** (v0.213): `ToNumberField[a, theta]` with a `Root` theta escalates
membership precision at every degree and evaluates. **A14b** (v0.217): `qa_resolve_extension`
(`src/poly/qafactor.c`) recognises a `Root[Function[...], k]` object as an algebraic generator —
its (monic-over-Q) defining polynomial becomes the minimal polynomial of `Q(c)`, rendered back in
terms of the Root; a `GEN_ROOT` autodetect kind drives `Extension -> Automatic`. The fix flows to
`Cancel`/`Together`/`PolynomialLCM`/`Quotient`/`Remainder`/`Factor`/`IrreduciblePolynomialQ` under
`Extension`. A **reducible** Root defining polynomial is declined (falls back to the gcd over Q).
The pre-fix inline comments above are kept as the historical snapshot. Consequently `ResidueClasses`
can now take that gcd over `Q(c)` directly rather than realising the whole logarithmic part as one
RootSum (the `hasRoot` workaround in the `.m` can be revisited).

### A15. `Sum` with a symbolic body, `Dot` with a symbolic vector, and `{} . {}`

```
Sum[Coefficient[t + 2 t^2, t, i] t^i, {i, 0, 1}]        (* -1/(-1 + t) + t^2/(-1 + t);  Mathematica: t *)
Coefficient[t + 2 t^2, t, i]                            (* 1 for a symbolic i;  Mathematica: unevaluated *)
{1, 2, 3} . t^Range[0, 2]                               (* a geometric closed form, not 1 + 2 t + 3 t^2 *)
{} . {}                                                 (* segmentation fault *)
```
`Sum` evaluates the body with the iterator symbolic (`Coefficient[..., i]` -> 1) and then sums `t^i` in
closed form. The `.m`'s truncation helper (`PMTruncate`) is an explicit `Do` loop over `CoefficientList`.

### A16. `PolynomialGCD[p, e, Extension -> Automatic]` returns 1 when `e` is an unexpanded algebraic constant that equals 0

```
c = 1/2 (-1 - I Sqrt[3]); e = c^3 - 1;                (* e is 0: Expand, RootReduce, Together, Simplify all give 0 *)
PolynomialGCD[x, e, Extension -> Automatic]           (* 1;  Mathematica: x *)
PolynomialGCD[x, Expand[e], Extension -> Automatic]   (* x *)
PolynomialGCD[x^2 + 1, e x + (c^2 + c + 1), Extension -> Automatic]   (* 1;  with Expand on the second argument: 1 + x^2 *)
```
`ResidueClasses` substitutes a residue value `c` into the characteristic polynomial of the residue and takes
the gcd with the prime. For the cube-root binomials `x^-1 (1 + x^2)^(1/3)`, `1/(x (x^2 - 1)^(1/3))`, ... the
quadratic factor `z^2 + z + 1` of the residue polynomial gave a trivial gcd, the two classes were dropped,
and -- the package then skipping a class with no places -- the holomorphic-remainder certificate fired on an
incomplete residue divisor: three FALSE "not elementary" verdicts (corpus R060, R131, R137) on the raw
v0.216 run. On the old core the residue polynomial itself came out degenerate ("degenerate residue
polynomial"), so this line was never reached. The `.m` now expands and reduces the coefficients
(`Collect[Expand[..], g, RRad]`) before every extension gcd in `ResidueClasses`, and an empty gcd for a
root of the residue polynomial -- mathematically impossible -- is an honest `{"failed", "residue class
lost: ..."}`; the same guard is in the research WL, Python and Maxima ports.

### A17. `RowReduce[m, Method -> "OneStepRowReduction"]` does not terminate on a small parametric matrix

```
Get["mixed/rowreduce_onestep_repro.m"]   (* aug: 20 x 15, entries linear in a, b -- the ansatz system of Exp[a x] Sin[b x] *)
RowReduce[aug]                                        (* 0.0 s, rank 12, the same reduced matrix as Mathematica *)
RowReduce[aug, Method -> "OneStepRowReduction"]       (* > 90 s, killed;  Mathematica: 1e-5 s *)
RowReduce[aug, Method -> "DivisionFreeRowReduction"]  (* 0.0 s *)
```
`AnsatzSystem` used the one-step method for every system (two orders of magnitude faster than the default
on AlgebraicNumber entries). With parameters present (`! FreeQ[aug, _Symbol?(! NumericQ[#] &)]`) the `.m`
now uses the default method; `Exp[a x] Sin[b x]` went from "time budget exceeded" to 1 s.

### A18. A `Do` iterator inside a package body captures a caller's symbol of the same name (`Table` does not)  — **FIXED** (v0.238; and `Sum`/`Product` had it too, see F3)

```
BeginPackage["P`"]; g::usage = "g"; h::usage = "h"; Begin["`Private`"];
g[v_] := Module[{s = 0}, Do[s += v, {a, 3}]; s];   h[v_] := Table[v, {a, 3}];
End[]; EndPackage[];
g[a]                    (* 6;         Mathematica: 3 a  (its iterator is P`Private`a) *)
h[a]                    (* {a, a, a}: correct *)
g[b]                    (* 3 b:       correct *)
Names["P`Private`*"]    (* {P`Private`s$1}: no P`Private`a is ever created *)
```
The A11 fix (v0.212) covers `Module` locals and pattern variables; a `Do` iterator is still bound by name.
`AnsatzSystem` iterated over the monomials with `{a, Length[monos]}` and over the matrix entries with
`{e, ...}`, so every integrand with a parameter named `a` or `e` was assembled with the parameter replaced
by loop indices: `Log[x + a]` "no solution within bounds", `x Exp[a x]`, `Exp[a x] Sin[b x]`, `Log[e + x]`,
`x^a` a wrong surface, rejected by the verify gate (never a wrong answer). The `.m` renames those iterators
(`ma`, `ent`, `ne`); the other single-letter iterators of the package (`i`, `k`, `c`, `n`, `p`, `g`, ...
about 150 sites) stay exposed to a parameter of that name until `Do` is made capture-avoiding like `Table`.

### A19. `FreeQ[e, Complex]` is True on a complex atom

```
FreeQ[1 + 2 I, Complex]           (* True;  Mathematica: False *)
FreeQ[1 + 2 I, _Complex]          (* False: correct *)
MatchQ[1 + 2 I, Complex[_Integer | _Rational, _Integer | _Rational]]   (* True: correct *)
(3 + 2 I) x /. Complex[a_, b_] :> a + b ii                              (* (3 + 2 ii) x: correct *)
```
A bare symbol as the second argument of `FreeQ` matches compound heads (`FreeQ[Log[x] + 1, Log]` is False)
but not the type of an atom. `logrewrite.m` writes every test for I with `_Complex` (and `_Root`, `_Re`,
... for the heads); found 2026-09-28 when the rule `Power[b_, 1/2] /; ! FreeQ[ComplexExpand[b], Complex]`
of `logrewrite.wl` never fired.

### A20. `ComplexExpand` writes a real nested radical as `Cos[Arg[...]]` terms; `Arg[1 - Sqrt[5]]` stays unevaluated

```
ComplexExpand[Sqrt[2 + Sqrt[3]]]      (* (7 + 4 Sqrt[3])^(1/4) Cos[1/2 Arg[2 + Sqrt[3]]] + I (...) Sin[...];  Mathematica: Sqrt[2 + Sqrt[3]] *)
ComplexExpand[Sqrt[1 - Sqrt[5]]]      (* (6 - 2 Sqrt[5])^(1/4) Cos[1/2 Arg[1 - Sqrt[5]]] + ...;          Mathematica: I Sqrt[-1 + Sqrt[5]] *)
Arg[1 - Sqrt[5]]                      (* unevaluated;  Mathematica: Pi *)
ComplexExpand[(-1)^(1/3)]             (* 1/2 + (I/2) Sqrt[3]: correct *)
```
The sign of a real algebraic number under a radical is not decided. `logrewrite.m` puts a constant in
rectangular form by rules (`RectConst`/`RectPow`: the sign of a real base from 30 digits, the principal
square root of a non-real base written by hand) and applies `ComplexExpand` to `(-1)^k` only. The
package's `RRad` (`ComplexExpand[ToRadicals[RootReduce[e]]]`) inherits the hazard on nested radicals.

### A21. `CountRoots` is not implemented

```
CountRoots[x^4 + 1, x]                (* unevaluated;  Mathematica: 0 *)
```
`logrewrite.m` counts the real roots with Sturm's theorem (`SturmCount`, the squarefree part and a
remainder sequence with `RootReduce`-canonical coefficients).

### A22. `$InputFileName` and `DirectoryName` are not implemented

```
$InputFileName                        (* unevaluated inside a file run with -file;  Mathematica: the file's path *)
DirectoryName["/a/b/c.m"]             (* unevaluated;  Mathematica: "/a/b/" *)
```
A package cannot `Get` a file beside itself the way `ParallelMixed.wl` reads `logrewrite.wl`;
`ParallelMixed.m` uses `LoadModule["mixed/logrewrite.m"]`, which resolves `src/internal/` like the lazy
load and evaluates the file in the current (private) context.

### A23. `NumericQ` is False on a `Root` object; `ToRadicals` gives the Ferrari form for every quartic root

```
NumericQ[Root[1 + #^4 &, 2]]                          (* False;  Mathematica: True *)
ToRadicals[Root[1 + 2 # - 2 #^2 + 2 #^3 + #^4 &, 4]]  (* -1/2 + 1/2 (Sqrt[-4 (-7/4 - 1/6 (-7 + (1/2 (-416 + 240 Sqrt[3]))^(1/3) + ...;
                                                         Mathematica: (-1 + Sqrt[5] + I Sqrt[2 (-1 + Sqrt[5])])/2 *)
ToRadicals[Root[1 + #^4 &, 2]]                        (* (-1 + I)/Sqrt[2]: correct *)
```
`N[Root[...], 30]` and `ToRadicals` of biquadratic roots are correct. `logrewrite.m` replaces every Root
object before a `NumericQ` test (`RootFree`) and writes a quartic root through its quadratic-in-disguise
form when the kernel's radical form is not tame (`RootRadicals`: biquadratic, palindromic and
antipalindromic quartics, the root picked by its 30-digit value). Charlwood P4's constants are the
roots of the palindromic `1 + 2 z - 2 z^2 + 2 z^3 + z^4`.

### A24. The field detour of `Can` returned a `Dot[{}, Inverse[{}], {}]` coefficient

Observed 2026-09-28 in one kernel that had integrated Charlwood P4 (`Log[1 + x Sqrt[1 + x^2]]`) and
then A40 (`ArcTan[x Sqrt[1 - x^2]]`): a logand of the rewritten logarithmic part came out as
`Log[Dot[{}, Inverse[{}], {}] + x^2 + 2 Dot[{}, Inverse[{}], {}] Sqrt[1 - x^2]]`, the back-conversion
`(v . Binv) . basisRad` of `FieldData` applied to an empty coordinate vector; the verify gate rejected
the surface. State-dependent (A40 alone spent its 30 s rewrite budget instead), no standalone repro
yet. `logrewrite.m` canonicalises with `CanRaw` (Cancel over the extension, which is what `Can` is in
`ParallelMixed.wl`) and writes the m = 2 quotient of `YQuot` with the conjugate formula instead of
`Pdiv`.

### A25. `Series` did not resolve a series variable that carries an OwnValue  (FIXED, v0.229)

`Series` and `SeriesCoefficient` are `HoldAll`, so the series variable arrives unevaluated.
Mathematica resolves one whose OwnValue names another SYMBOL; Mathilda treated the held symbol
literally, so an expression that does not contain it is constant in it and the "series" is the
input itself -- returned silently, with no message:

```
Module[{w}, w = Unique["w"]; Normal[Series[Sqrt[1 + w^4], {w, 0, 4}]]]
  Mathematica  1 + w11^4/2
  Mathilda     Sqrt[1 + w11^4]        (* before v0.229 *)
```

A NUMERIC value leaves the call unevaluated in Mathematica (`Series[f, {5, 0, 4}]`) and an
EXPRESSION value is treated as a variable the input is constant in; both are matched.
`ParallelMixed.m`'s `RealiseClass` expands `Sqrt[q(1/w) w^(2h)/lc]` in a `Unique["w"]`, so every
unbalanced configuration at infinity was solved against a matrix carrying a square root of a
polynomial where a rational belonged. Fixing it took the review corpus from 304 to 329 correct.

### A26. Multivariate `PolynomialGCD` over a number field was wrong  (FIXED, v0.230)

```
a = AlgebraicNumber[Sqrt[2], {0, 1}];   d = x^3 + Sqrt[2] x y + y^2 + 1;
                                                    Mathematica        Mathilda (v0.229)
PolynomialGCD[(x+a y)(x+1), (x+a y)(x+2)]              1              (x+a y)(x+2)  NOT A DIVISOR
  ... , Extension -> Automatic]               2y + x AlgebraicNumber  (x+a y)(x+2)  NOT A DIVISOR
PolynomialGCD[x + a y, (x+a y)(x+2)]                   1              (x+a y)(x+2)  degree > operand
PolynomialGCD[d(x^2+Sqrt[2]y+3), d(x^2+2Sqrt[2]y-1), Extension -> Automatic]
                                                       d                     1      factor missed
```

**Three separate defects, only the first of which this entry originally recorded.**

*The wrong answer.* `collect_variables` correctly treats an `AlgebraicNumber` as a CONSTANT
(v0.229), so the coefficient ring is `K`. But `poly_content` bottoms out in `my_number_gcd` →
`get_int_content` (`src/poly/poly.c:1251-1291`), which is INTEGER content and returns 1 for any
`K`-coefficient. Both operands therefore enter the pseudo-remainder sequence non-primitive over
`K`, and on the FIRST iteration `pseudo_rem` computes `lc(B)*A - lc(A)*B`, which vanishes
identically whenever the two share a factor and have equal degree in the main variable. The loop
then does `U = V`, breaks, and returns `contGCD * U` — **the second operand**, which is not a
common divisor at all. `PolynomialGCD[f, g]` returned `g` and `PolynomialGCD[g, f]` returned `f`.

*The missed factor.* `Extension -> Automatic` on RADICAL input reached the Phase D tower path,
which computes the `Q[gamma, x, y]`-GCD rather than the `Q(gamma)[x, y]`-GCD (its own comment
says so at `qafactor.c:2798`). That is right only while the cofactors carry no algebraic
constants, and silently answers 1 when they do.

*The unreached engine.* Every explicit `Extension -> <value>` form (a value, a list, `All`,
`None`) went to the classical path too — only `Extension -> Automatic` reached any field engine.

**The fix** is the native `flint_field_gcd` named here (`src/poly/flint_bridge.c`), hooked into
both `poly_gcd_internal` and `builtin_polynomialgcd` so `Factor`, `SquareFreeQ`, `FactorTerms`,
`PolynomialLCM`, `Cancel`/`Together` and the Risch consumers are covered too. FLINT has no
multivariate GCD over a number field (`gr_mpoly.h` declares no arithmetic at all), so it is the
classical modular (Encarnación) algorithm over `fq_nmod_mpoly_gcd`, CRT and rational
reconstruction; the residue ring has to be SPLIT into `M mod p`'s irreducible factors and the
component gcds CRT'd back, because for a non-cyclic Galois group such as `Q(sqrt2, sqrt3)` no
prime keeps `M` irreducible. Every answer is certified before return — the candidate is monic so
its leading monomial is tau-free while `M`'s is `tau^n`, making `{G, M}` a Gröbner basis by
Buchberger's first criterion, so `fmpq_mpoly_divrem_ideal` decides divisibility over `K` exactly.
Radical input is mapped into one common field and rendered back through the PRODUCT BASIS of the
caller's own atoms, so radicals in gives radicals out rather than qqbar's `Root` spelling of the
primitive element. `MATHILDA_NO_FIELD_GCD=1` A/Bs it.

Anything the engine still declines (a build without FLINT, a generator qqbar does not model, an
input that never certified) is **checked** at the end of the classical multivariate path: if the
PRS answer does not divide both operands it is replaced by **1** — always a common divisor, and
what Mathematica returns for these inputs under its default `Extension -> None`. Checked rather
than pre-empted: refusing every multivariate input that merely carries an algebraic constant is
too blunt, and cost a `DSolve` case whose answer runs through `Q(Sqrt[17])`. Univariate is never
checked (one variable means the un-stripped `K`-content only scales the answer, so the PRS lands
on the true gcd) and neither is `Q(i)`, whose content `get_int_content` does strip.

Remaining, deliberate: the gcd over a field is defined only up to a constant of that field, and
Mathilda returns the MONIC associate where Mathematica returns another (`x + Sqrt[2] y` against
Mathematica's `Sqrt[2] x + 2 y`). Both are correct; tests assert by divisibility, not by printed
form. `logrewrite.m`'s `KGcd` no longer needs its multivariate workaround.

Still declined, and still answering 1: an expression that MIXES spellings of the same field
(`AlgebraicNumber[Sqrt[2], ..]` alongside a bare `Sqrt[3]`), or carries two structurally distinct
`AlgebraicNumber` generators. `field_scan` reports a conflict and the engine declines rather than
building the compositum. Pre-existing and unchanged: before v0.230 the mixed case answered with
an unsimplified **zero** (`Sqrt[3] AlgebraicNumber[Sqrt[2], {0, -4}] + Sqrt[3]
AlgebraicNumber[Sqrt[2], {0, 4}] + ...`), which the divisibility check does not catch because it
is not structurally zero.

#### A26a. What the v0.232 stress pass found  (FIXED, v0.232)

The engine shipped at v0.230 correct on its tests and faster than the baseline on Charlwood, but
its tests were all a handful of terms in a degree 2–4 field. Pushed on coefficient size, field
degree, term count and variable count it turned out to have three silent failures — all of which
**answered**, which is why none had been noticed. `tests/bench_field_gcd.c` is the harness, and it
treats a decline as a failure rather than a slow row precisely because that is the failure mode.

| | v0.230 | v0.232 |
|---|---|---|
| coefficient ceiling | declines past **~831 bits**, so the caller answers 1 | no decline at 13,288 bits |
| `[K:Q] = 6` via radicals | declines for **every** generator | works |
| compositum `{Sqrt[2], Sqrt[3]}` | correct but spelled as a degree-4 `Root` per coefficient | `1 + x^3 + Sqrt[2] x y + y^2` |
| 628-term operands | 63 ms | 30 ms |
| residue-field work | — | 2.4–2.6× less |

*The ceiling* was arithmetic, not subtlety: `FG_MAX_PRIMES 64` × 29-bit primes is 1856 bits of
modulus and rational reconstruction needs about twice the coefficient size. Primes are now 62-bit
and the cap is a backstop rather than a limit, which is safe because the certificate — not the
budget — is what makes the answer correct, so extra primes only ever cost time.

*The degree-6 decline* was upstream, in the qqbar compositum. It chose a primitive element by
TRIAL membership, and `qqbar_express_in_field_esc` deliberately does not escalate past 64 bits of
working precision when the generator has degree `<= 6`. Degrees 2–5 resolve inside 64 bits and 7+
are allowed to escalate, so **degree 6 alone** failed — for every candidate `c`, hence the whole
field, hence `ToNumberField[{2^(1/6), 2^(1/3)}]` unevaluated and every degree-6 radical gcd
answering 1. The primitive element is now chosen by DEGREE: `alpha + c*b` always lies in
`Q(alpha, b)`, so `Q(alpha + c*b)` equals the compositum exactly when their degrees agree, and no
membership test that could fail for an unrelated reason is needed.

*The render-back* gave up whenever the operands mentioned more atoms than a basis needs — and they
routinely do, because `Expand` folds `Sqrt[2] Sqrt[3]` into `Sqrt[6]`, making the atom set of a
`Q(sqrt2, sqrt3)` problem `{sqrt2, sqrt3, sqrt6}` with a product basis of 8 against `[K:Q] = 4`.
The overshoot was read as "these atoms do not span". Atoms are now taken greedily, largest degree
first, until the product reaches `n`; the invertibility of the change-of-basis matrix is the test
that they really span, so nothing is assumed.

The speed came from choosing primes that split LEAST, since `fg_gcd_mod_p` runs one multivariate
gcd per irreducible factor of `M mod p` and those gcds are essentially the entire cost (0.062 s
against 0.008 s for everything else combined, over 200 calls). The direction was measured, not
assumed — a split prime has cheaper coefficient arithmetic but pays FLINT's per-call cost `r`
times, and choosing the most-split prime instead was 1.6× slower at `n = 2` and 2.6× slower at
`n = 8`. `MATHILDA_FIELD_GCD_STATS=1` prints the per-stage profile.

Two further findings, both **pre-existing and outside this engine** (identical at v0.230 and with
`MATHILDA_NO_FIELD_GCD=1`), recorded here rather than fixed: `PolynomialGCD[0, f]` with algebraic
coefficients answers 1 where Mathematica answers `f`; and a coefficient mixing an inexact real
with an algebraic constant returns a garbage near-zero float
(`PolynomialGCD[Expand[(x + 1.5 Sqrt[2] y)(x+1)], Expand[(x + 1.5 Sqrt[2] y)(x+2)]]` →
`3.71618e-16`) rather than declining. The pure-float case is fine, so it is specifically the
mixture.

### A27. `SparseArray` densifies but is still not an array  (PARTLY FIXED, v0.233)

```
Normal[SparseArray[{{1, 2} -> 5, {2, 1} -> 7}, {2, 3}]]                {{0, 5, 0}, {7, 0, 0}}   (* fixed *)
Normal[SparseArray[Automatic, {2, 2}, 0, {1, {{0, 1, 2}, {{1}, {2}}}, {9, 8}}]]  {{9, 0}, {0, 8}}   (* fixed *)

s = SparseArray[{{1, 2} -> 5, {2, 1} -> 7}, {2, 3}];
                  Mathilda                Mathematica
Dimensions[s]     {2}                     {2, 3}
s[[1, 2]]         {2, 1} -> 7             5
MatrixQ[s]        False                   True
ArrayRules[s]     unevaluated             {{1, 2} -> 5, {2, 1} -> 7, {_, _} -> 0}
s . {1, 2, 3}     unevaluated             {10, 7}
```

`Normal` was taught to densify by the worked-example sweep (`src/sparsearray.c`, v0.231, merged at
v0.233) — both the rule-list spelling and the `Automatic` CSR one — so the repro this entry was
written for now matches Mathematica. The head itself is still inert: it registers no builtin, so
every *other* operation reads `SparseArray[rules, dims]` as an ordinary two-argument expression.

`Dimensions[s]` answering `{2}` is the case to watch, because it is the argument count wearing the
shape of an answer rather than an unevaluated form a caller would notice; `s[[1, 2]]` returning the
second rule is the same failure. So the original warning stands for everything except `Normal`: a
caller that assembles a matrix this way still gets a non-matrix with no message and no error, and
the downstream `RowReduce` still "solves" it. Found while rewriting the ansatz assembly, which
builds its augmented matrix row by row instead.

## B. Behavioural differences (not wrong, but code written against Mathematica's output breaks)

### B1. `ToNumberField` chooses a non-canonical primitive element, and a different one per call

```
ToNumberField[{Sqrt[2]/8, 1, -1, Sqrt[2]}]  -> {AlgebraicNumber[-4 Sqrt[2], {0, -1/32}], 1, -1, AlgebraicNumber[-4 Sqrt[2], {0, -1/4}]}
ToNumberField[{Sqrt[2]/8, 1}]               -> {AlgebraicNumber[32 + 4 Sqrt[2], {-1, 1/32}], 1}
ToNumberField[{I, 2}]                       -> {AlgebraicNumber[2 + I, {-2, 1}], 2}
  Mathematica: AlgebraicNumber[Sqrt[2], {0, 1/8}], ..., AlgebraicNumber[I, {0, 1}]
```
The arithmetic is right (theta is a primitive element and the coordinates are correct), but anything that
pattern-matches `AlgebraicNumber[Sqrt[2], _]`, compares thetas across calls, or expects small
discriminants (`MinimalPolynomial[32 + 4 Sqrt[2], z]` = `992 - 64 z + z^2`) breaks. The `.m` makes one
`ToNumberField` call per field and reads the minimal polynomial of whatever theta comes back.

### B2. `ToString[expr, InputForm]` prints `1/(4 Sqrt[2])` as `1/4/Sqrt[2]`
Mathematica: `1/(4*Sqrt[2])`. Re-parses to the same value; cosmetic (it shows in the certificate's divisor string).

### B3. `Print` output to a pipe is fully buffered
Nothing appears until the process exits, and a process killed by a timeout loses all of it (Mathematica's
`wolframscript` flushes per `Print`). Every `timeout -s KILL ./Mathilda -file x.m` diagnostic in this session
had to be split into one process per statement. A line-buffered (or `Print`-flushing) mode would help.

### B4. `Select`/`Range` results are packed and interact with A1
`Select[Range[3], # > 5 &] === {}` is True and `Length` is 0, so empty-list tests are fine; only iteration (A1)
is affected.

### B5. `TimeConstrained` does not interrupt a kernel-level wait
`TimeConstrained[Pause[3]; 1, 1, "TO"]` returns 1 after 3 s (Mathematica: "TO" after 1 s). The budget of
`ParallelIntegrateMixed` is checked between evaluation steps only; a long single kernel operation
(a large `Cancel`, `Solve`) runs to completion.

### B6. `Exp[c Log[x]]` and `Exp[x Log[2]]` evaluate back to `x^c`, `2^x` (as Mathematica does)
The tower builder holds such exponentials as `PMExp[c Log[x]]` until the generator is created (same in
the research `ParallelMixed.wl`).

### B7. An outer `TimeConstrained` cannot interrupt an inner one
```
TimeConstrained[TimeConstrained[While[True, 1], 0.3, "inner"], 5, "outer"]   (* "inner": correct *)
TimeConstrained[TimeConstrained[While[True, 1], 5, "inner"], 0.3, "outer"]   (* "inner" after 5 s;  Mathematica: "outer" after 0.3 s *)
```
The 45 s budget of `ParallelIntegrateMixed` therefore cannot cut short the rewrite of the logarithmic
part, which runs under its own `$RewriteBudget` (30 s) inside it.

## C. Performance of polynomial arithmetic over GF(p) (the reason the new certificate exceeds the 45 s budget)

Measured on the same machine, same algorithm (`CantorAdd`, Mumford pairs, `PolyExtGCDModP`,
`PolyQuoModP`):

| operation | Mathilda | Python (sympy `galoistools`) | Mathematica |
|---|---|---|---|
| one group operation in Jac(y^2 = f) over GF(17), genus 4 | 10-20 ms (builtin gcd) / ~90 ms (own gcd) | 0.1 ms | ~0.5 ms |
| order of a class, genus 2, p = 7 (BSGS, ~80 operations) | 3.7 s | 0.01 s | 0.05 s |
| certificate for the genus-2 form `1/(x (1 - x^2) (1 + x^5)^(3/2))` (p = 7, 17) | 9.7 s (17.9 s with the pipeline; inside the 45 s budget) | 0.7 s | 0.2 s |
| certificate for `Coth[x]/(1 + Sech[x]^5)^(3/2)` (genus 4, p = 17, 41, ~7000 operations) | 512 s (correct: `{True, {{17, 29}, {41, 155}}}`; measured with another job on the machine) | 0.8 s | 3 s |

The per-operation cost is evaluator overhead: every `ReduceModP` is a `CoefficientList` -> `Sum` round trip
through expressions, every `PolynomialQuotient` runs over Q with rationals. The fix that removes both the
workarounds and the cost is in the C core: implement `Modulus -> p` for `PolynomialQuotient`,
`PolynomialRemainder`, `PolynomialQuotientRemainder`, `PolynomialExtendedGCD` (correctly, A2) and
`PolynomialMod` on FLINT's `nmod_poly` (already linked), so that the straightforward Wolfram-language
port (`ParallelMixed.wl`, which uses the `Modulus` forms directly) runs as is.

Consequence today: `Integrate`ParallelMixedTower[Coth[x]/(1 + Sech[x]^5)^(3/2), x]` returns
`{"failed", "time budget exceeded"}` (a clean decline, never a wrong answer) where Mathematica and Python
return the certificate `{"not elementary", "residue divisor not torsion: reduction mod p", {{17, 29}, {41, 155}}, ...}`;
`Sqrt[1 + Sec[x] Tan[x]]`, which now reaches `TorsionRealise` for the first time (A4), also exceeds the budget.

## D. Workarounds present in `src/internal/mixed/ParallelMixed.m` (remove when the core is fixed)

| tag / function | replaces | item |
|---|---|---|
| `PadRows` in `QCoords`, `NontorsionDivisor` | `PadRight[rows]` | A4 |
| `ReduceModP`, `PolyQuoModP`, `PolyRemModP` (also in `UnitDegreeModP`) | `PolynomialQuotient[..., Modulus -> p]` etc. | A3, A6 |
| `PolyExtGCDModP` in `CantorAdd` | `PolynomialExtendedGCD[..., Modulus -> p]` | A2 |
| index loops in `NontorsionCertificate`, `NontorsionDivisorCertificate` | `Do[..., {p, Select[...]}]` | A1 |
| `If[KeyExistsQ[baby, G], baby[G], None]` in `DivisorOrderModP` | `Lookup[baby, Key[G], None]` | A5 |
| `Thread[{...}]` (several sites) | `Transpose[{...}]` | A8 |
| `Options[iPIM] = ...` | -- | A9 |
| one `ToNumberField` call per field in `NontorsionDivisor` | -- | B1 |
| `hasRoot` in `ResidueClasses`, RootSum built with `Apply`, residual by `ModP` (`RootSumLogand`) | the gcd over `Q(c)` and the kernel's RootSum evaluation | A12, A14 |
| `LaurentPolyTimes` / `InfLaurent` (exact expansions at infinity), `PMTruncate` as a `Do` loop | `Series`, `Sum` over a symbolic iterator | A13, A15 |
| `Together` on the reduced matrix and on the residual with symbolic parameters (`AnsatzSystem`) | -- | A11 (partial) |
| parameters substituted by fixed rationals in the verify gate | -- | A11 (partial) |
| `Collect[Expand[..], g, RRad]` before every extension gcd in `ResidueClasses` (and the honest failure for a lost class) | -- | A16 |
| the default `RowReduce` method when the ansatz matrix has parameters (`AnsatzSystem`) | `Method -> "OneStepRowReduction"` for every system | A17 |
| iterators `ma`, `ent` in `AnsatzSystem` and `ne` in the realisation | `a`, `e` | A18 |
| `logrewrite.m` (the real form of the logarithmic part): `_Complex` tests; `RectConst`/`RectPow` by rules; `SturmCount`; `RootFree`/`RootRadicals`; `CanRaw` and the conjugate formula in `YQuot`; lr-prefixed `Do` iterators; `LoadModule["mixed/logrewrite.m"]` from `ParallelMixed.m` with a default `LogToReal` when the file is missing | `FreeQ[.., Complex]`, `ComplexExpand`, `CountRoots`, `ToRadicals`/`NumericQ` on Root objects, `Can`, `Get` of a sibling file | A18-A24 |

The research copy `ParallelMixed.wl` (Mathematica) has none of these and is the reference for the intended
code -- except that it carries the A16 expansion and honest failure (a soundness guard in every port), the
A17 method choice and the A18 renames as well, all harmless there.

<!-- charlwood-300 -->
## E. Charlwood's fifty integrals at the 300 s protocol (2026-09-23, build 0.175)

Run for the Maxima paper (`rn-radicals-maxima.tex` in the research directory) with
`python3.11 charlwood_wl.py mathilda` (a copy of the runner is in `mixed/charlwood_wl.py`; it
needs the research `charlwood.py` next to it and `MATHILDA_BIN`): one `Mathilda -file` process
per integral, `Integrate[f, x, Method -> "ParallelMixedTower"]` after two warm-up integrals,
`$ParallelMixedTimeBudget = 300` (the symbol is `System`$ParallelMixedTimeBudget`; the bare name
reaches it through `$ContextPath`), the answer parsed into SymPy and verified by `charlwood.verify`
(1e-18 at three rational points, 30 digits); the reason of a decline from
`Integrate`ParallelMixedTower[f, x]` afterwards.  Two runner details that matter in Mathilda: the
unevaluated `Integrate[...]` of a decline is at once rewritten `/. Integrate -> PMDeclined`, because
Mathilda keeps no evaluated-flag on expressions and every later reference to it (even `Head[r]`)
runs the method again; and the reason call gets the full budget when the decline came inside the
limit (A35 declines at 232 s with "solution does not verify"; asked under the 68 s that were left it
reported "time budget exceeded").  Result file:
`/Users/user/Documents/Research/post_phd_research/algebraic_integration/risch_norman_radicals/mixed/charlwood_mathilda.json`.  Times are the kernel's `AbsoluteTiming` of the call.
SymPy = `t_tower + t_integrate` of the Charlwood paper's run, Mathematica = `ParallelMixed.wl`
under 14.0, Maxima = the port of the same date.

**Summary (50 of 50 run): 49 verified** (35 in the run of the morning of 2026-09-23 on
0.174; the fixes in between are in E.3).  On the integrals it solves Mathilda is fast: faster than
SymPy on 45 of 49, median ratio Mathilda/SymPy 0.45, 25.5 s against
55.0 s on those; faster than Mathematica on 15, slower than Maxima on all but one.

**This is the practical ceiling.** The single remaining miss is **A39**,
`ArcSin[x Sqrt[1 - x^2]]`, which is genuinely non-elementary — every one of the four
CAS fails it, so 49 of 50 is the score to beat and Mathilda has it. The `Root`-object
rungs that used to time out or over-assemble (P4, P8, A2, A3, A27, A35, A40) all
solve now, and nothing declines at the 300 s budget any more; the slowest case in the
table is well inside it.

Re-recorded 2026-10-01 at v0.244 by `mixed/stress/charlwood_record.py`. The previous
table in this section was from build 0.174 (42 verified), so the 42 → 49 step is the
accumulated work of v0.175–v0.243 and **not** attributable to any single change; the
v0.244 verify-gate change that prompted the re-run was measured separately against a
pristine v0.243 worktree and is per-case identical (50/50 run, 49 verified, no case
changed status either way).

| id | Mathilda | s | SymPy s | Mathematica s | Maxima s | reason |
|---|---|---|---|---|---|---|
| P1 | ok | 0.20 | 0.49 | 0.10 | 0.16 |  |
| P2 | ok | 0.08 | 0.16 | 0.04 | 0.03 |  |
| P3 | ok | 0.39 | 0.61 | 0.18 | 0.11 |  |
| P4 | ok | 2.17 | 5.52 | 1.90 | 1.70 |  |
| P5 | ok | 0.60 | 1.90 | 0.38 | 0.13 |  |
| P6 | ok | 0.16 | 0.40 | 0.13 | 0.06 |  |
| P7 | ok | 0.05 | 0.12 | 0.03 | 0.02 |  |
| P8 | ok | 1.75 | 6.59 | 1.07 | 1.36 |  |
| P9 | ok | 0.11 | 0.32 | 0.06 | 0.05 |  |
| P10 | ok | 0.30 | 0.44 | 0.10 | 0.10 |  |
| A1 | ok | 0.48 | 0.98 | 0.77 | 3.46 |  |
| A2 | ok | 1.46 | 3.28 | 0.93 | 0.73 |  |
| A3 | ok | 1.50 | 2.85 | 0.97 | 0.75 |  |
| A4 | ok | 0.15 | 0.64 | 0.09 | 0.04 |  |
| A5 | ok | 0.35 | 0.71 | 0.30 | 0.15 |  |
| A6 | ok | 0.39 | 0.92 | 0.24 | 0.17 |  |
| A7 | ok | 0.23 | 0.56 | 0.19 | 0.12 |  |
| A8 | ok | 0.11 | 0.27 | 0.08 | 0.06 |  |
| A9 | ok | 0.03 | 0.14 | 0.01 | 0.02 |  |
| A10 | ok | 0.08 | 0.11 | 0.04 | 0.03 |  |
| A11 | ok | 0.40 | 0.79 | 0.25 | 0.16 |  |
| A12 | ok | 0.40 | 1.05 | 0.28 | 0.16 |  |
| A13 | ok | 0.26 | 0.35 | 0.16 | 0.43 |  |
| A14 | ok | 0.12 | 0.29 | 0.10 | 0.07 |  |
| A15 | ok | 0.11 | 0.08 | 0.03 | 0.03 |  |
| A16 | ok | 0.64 | 1.17 | 0.46 | 0.34 |  |
| A17 | ok | 0.12 | 0.24 | 0.08 | 0.07 |  |
| A18 | ok | 0.13 | 0.24 | 0.07 | 0.08 |  |
| A19 | ok | 0.48 | 0.96 | 0.46 | 0.11 |  |
| A20 | ok | 0.15 | 0.52 | 0.18 | 0.09 |  |
| A21 | ok | 0.15 | 0.60 | 0.30 | 0.08 |  |
| A22 | ok | 0.20 | 0.35 | 0.17 | 0.04 |  |
| A23 | ok | 0.12 | 0.24 | 0.09 | 0.06 |  |
| A24 | ok | 0.14 | 0.22 | 0.08 | 0.05 |  |
| A25 | ok | 0.15 | 0.28 | 0.09 | 0.06 |  |
| A26 | ok | 0.11 | 0.18 | 0.08 | 0.05 |  |
| A27 | ok | 1.52 | 3.64 | 0.41 | 0.09 |  |
| A28 | ok | 4.32 | 3.68 | 0.58 | 0.22 |  |
| A29 | ok | 0.37 | 2.83 | 0.31 | 0.20 |  |
| A30 | ok | 0.12 | 0.26 | 0.06 | 0.04 |  |
| A31 | ok | 0.25 | 0.59 | 0.17 | 0.21 |  |
| A32 | ok | 1.25 | 0.45 | 0.17 | 0.03 |  |
| A33 | ok | 0.17 | 0.37 | 0.13 | 0.05 |  |
| A34 | ok | 0.15 | 0.45 | 0.15 | 0.04 |  |
| A35 | ok | 0.66 | 3.07 | 0.37 | 0.15 |  |
| A36 | ok | 0.11 | 0.10 | 0.07 | 0.07 |  |
| A37 | ok | 0.35 | 1.99 | 0.31 | 0.17 |  |
| A38 | ok | 0.14 | 0.33 | 0.08 | 0.06 |  |
| A39 | F | 0.56 | 1.61 | 1.05 | 0.20 | {"no solution within bounds", {2, 2}} |
| A40 | ok | 1.74 | 2.65 | 2.16 | 0.20 |  |

### E.1 Trace-level diagnostics (traces in `mixed/stress/charlwood_traces/`, `mathilda_0175_<id>.txt` against `mathematica_P4_A2_A3_A27_A35_A40.txt`; budget 90 s for the traces)

Verbose traces (`ParallelMixed`ParallelIntegrateMixed[f, x, "Verbose" -> True]`) of Mathilda and of
Mathematica running the same package agree line for line on the classification, the residues,
the S'-units and the bounds, and differ exactly here:

- **A2, A3** (`ArcTan[x + Sqrt[1 - x^2]]`, `x ArcTan[x + Sqrt[1 - x^2]]/Sqrt[1 - x^2]`, 300 s):
  on the conic rung split over `{I + w, -I + w, w - Root[1 - 2# + 2#^2 + 2#^3 + #^4 &, k]}`
  Mathematica assembles **34 equations in 15 unknowns, solved** (1.0 and 1.1 s in all);
  Mathilda assembles **72 and 73 equations**, no solution, and moves on to the split of the base
  rung over `Root[1 - #^2 + #^4 &, k]` with guessed bounds ({4, 2}: 84 equations; {8, 2}: ...)
  until the budget.  Same defect as A1's 81-vs-80 before the AlgAtoms fix, now with Root
  coefficients: a coefficient that is not cancelled to zero.
- **A35** (`Sqrt[Sqrt[Sec[x] + 1] - Sqrt[Sec[x] - 1]]`, 232 s, `{"failed", "internal: solution
  does not verify", {1, 0}}`): on the rung split over `u - Root[2 - 2#^2 + #^4 &, k]` Mathematica
  assembles **8 equations in 10 unknowns, solved** (0.36 s); Mathilda **16**, "solved", and the
  exact post-solve check rejects the solution -- after 230 s spent in that check.
- **P8** (`Sqrt[Tan[x]^2 + 2 Tan[x] + 2]`, 30 s, `{"failed", "no solution within bounds", {1, 0}}`):
  unchanged from the first run: first rung 18 equations against 12; conic rung split over
  `Root[1 - 2# + 2#^3 + #^4 &, k]`: 14 equations against 4, "solved" and rejected; then the
  x-split `{I + t, -I + t}` (12 equations, no solution) and the decline.
- **A40** (`ArcTan[x Sqrt[1 - x^2]]`, 300 s): the FIRST rung -- `18 unknowns, 49 equations` in both
  systems -- is **solved by Mathematica (2.1 s) and has no solution in Mathilda**: the same count,
  so a wrong coefficient rather than a spurious equation (the S'-units of `x^4 - x^2 - 1` carry
  Sqrt[5]).  The ladder then parametrises the conic and splits over
  `{I + w, -I + w, +-Sqrt[Root[1 + 8# - 2#^2 + 8#^3 + #^4 &, k]] + w}`, whose assembly never prints
  its equation count.
- **P4** (`Log[1 + x Sqrt[1 + x^2]]`, 300 s): base rung and conic rung as Mathematica (49/101/149
  and 34 equations, no solution); on the conic rung split over `Root[1 + 2# - 2#^2 + 2#^3 + #^4 &, k]`
  Mathematica assembles 34 equations in 15 unknowns and solves (2.1 s in all); Mathilda prints the
  `bounds:` line of that rung and **never its `ansatz:` line** -- the assembly over the compositum
  of the four roots does not finish in 300 s.
- **A27** (`Sqrt[Sin[x]]/(1 + Sin[x]^2)`, 300 s): the four residue classes over
  `1 + 6 u^4 + u^8` are found as in Mathematica; their realisation then emits
  `NullSpace::method: Method option value is not one of "Automatic", "DivisionFreeRowReduction",
  "OneStepRowReduction", "CofactorExpansion"` (178 times in 90 s) and `Dot::dotsh` (24 700 times):
  the `NullSpace[..., Method -> "OneStepRowReduction"]` of the Q-basis / group-law step is
  rejected, the `Dot` that follows gets mismatched shapes, and the loop never converges.
  Mathematica realises the classes and solves the first rung (7 unknowns, 10 equations) in 0.41 s.
- **A39**: `{"failed", "no solution within bounds", {2, 2}}` in 1.3 s, the honest failure of every
  system (not elementary; the certificate is withheld in this tower).

### E.2 What to fix (leads, in order)

1. The equation count over Root constants (A2, A3, A35, P8): diff the assembled systems against
   Mathematica's on the split rung; the extra equations are monomials whose coefficient is an
   algebraic number in the compositum of the Root objects that `Together`/`Cancel[...,
   Extension -> Automatic]`/`CoefficientRules` does not cancel to zero.  The AlgAtoms fix of E.3
   did exactly this for `I` (A1, A16, A19, A20, A37); the Root case needs the same canonical zero
   test (ToNumberField once per rung, as the .wl relies on, or RootReduce on every coefficient).
2. A40's first rung: same 49 equations, no solution -- solve the two systems side by side; the
   defect is in a coefficient over Q(Sqrt[5]) (the S'-unit logs) or in LinearSolve over it.
3. A27: accept `Method -> "OneStepRowReduction"` in `NullSpace` (or make the .m fall back to
   Automatic in Mathilda) and check the `Dot` shapes in the residue-class realisation.
4. P4 and A40's nested split: the assembly over a compositum of Root objects (degree 8-16) --
   the AlgebraicNumber Plus/Times round trip expr -> qqbar -> expr per operation (todo G4c);
   native nf_elem arithmetic through the assembly.
5. A35's 230 s post-solve check over the Root compositum: the same arithmetic, in the verification.
6. Minor: merge S'-units equal up to sign (A1: 5 listed, 3 distinct; 79 vs 77 unknowns).

### E.3 Fixed between the two runs of 2026-09-23 (0.174 -> 0.175 working tree; see tasks/todo.md)

- **AlgAtoms / FreeQ** (`FreeQ[e, Complex]` is True in Mathilda where WL is False, so `I` was never
  added to the field and the split rung over Q(i) carried one spurious equation): detection by
  `Cases[e, _Complex, ...]`; A1's split rung now `31 unknowns, 80 equations; solved` as in
  Mathematica (5.4 s), A16 `13 unknowns, 28 equations; solved` (0.9 s), and A19, A20, A37 -- whose
  conic rung had failed for the same reason and escalated to a degree-8 compositum -- solve in
  2.2, 0.1 and 0.3 s.
- **Branch-resolved gate**: the raw surface of A11 differentiated to -f everywhere (the other root
  of y^2 = q) and that of A34 to +-f across the branch cut of Sqrt[1 - Sin^6] = |Cos| Sqrt[...]; the
  gate now pins the sign numerically (accepts only exactly +-f at every real sample); A11 0.4 s,
  A34 0.1 s, both verified externally.
- **Guarded FLINT allocator**: `TimeConstrained` preempts inside FLINT work; every budget decline
  is now at 300.0 s (was: the 420 s OS kill, with the output lost in the pipe buffer).
- **Root-index memo** per minimal polynomial in `qqbar_to_expr`: A16 13 s -> 0.9 s, A19 39 s -> 2.2 s.
<!-- /charlwood-300 -->

<!-- special-stage-port -->
## F. Found porting the special-function stage (`ParallelMixedSpecial`), 2026-09-30, build v0.237

All six were reproduced against the binary at the repo root and checked against
Mathematica's documented behaviour. Unlike section A, **every one of these is now
FIXED at the root** (v0.238) rather than worked around in a `.m`, so this section
is a record of what the port exposed, not a list of live hazards.

### F1. `Block` did not restore `DownValues`  (FIXED, v0.238)

```
gg[a_] := "orig";
Block[{gg}, gg[a_] := "patched"; gg[1]]     (* "patched": correct                 *)
gg[1]                                        (* was "patched"; Mathematica: "orig" *)
```

Only `own_values` was saved and restored, so a rule written to a `Block`-ed
*function* symbol stayed installed for the rest of the session. This is the
mechanism the special stage's `ExtendedBounds` is built on — it `Block`s four of
Part II's bound-decision symbols and `SetDelayed`s extended versions inside — so
the **first** extended call would have permanently repointed them, silently
corrupting the existing `ParallelMixedTower` method. Attributes are deliberately
still not cleared: Mathematica's `Block` does not clear them either (verified).

A `TimeConstrained` timeout `siglongjmp`s past the restore, so frames are also
threaded on a stack that `tc_run_guarded` drains — the same treatment
`tc_async_deferred` and the message-suppression depth already had.

### F2. `Block` was alpha-renamed by capture-avoiding substitution  (FIXED, v0.238)

```
g[v_] := Block[{e = 1}, v + e];   g[e + 1]     (* was 2 + e;  Mathematica: 3     *)
g[v_] := Module[{e = 1}, v + e];  g[e + 1]     (* 2 + e: correct, lexical scope  *)
```

The A11 capture-avoidance fix (v0.212) treated every scoping construct alike, but
`Block` is **dynamic** scope: a symbol arriving inside a caller's value is exactly
what it means to rebind. Renaming the local defeated the construct, and the
visible symptom was that a held body handed into a hook-installing `Block` kept
reaching the *original* symbol.

### F3. `Sum` and `Product` iterators captured a caller's symbol  (FIXED, v0.238)

The other half of A18. `Sum[v, {k, 2}]` with `v = k` gave 3 (Mathematica `2 k`)
and `Product[v, {k, 2}]` gave 2 (Mathematica `k^2`). Only `Table` was listed as a
scoping construct; all four iterator heads are structurally identical and now are.

### F4. Nested `Association` element assignment was a silent no-op  (FIXED, v0.238)

```
a = <||>; a["k", "s"] = 7;  a       (* was <||>;  Mathematica: <|k -> <|s -> 7|>|> *)
a = <|"k" -> <||>|>; a["k","s"] = 7 (* worked: the intermediate key existed        *)
```

A deep write through a key that did not exist yet placed nothing: the `Set` was
left unevaluated, `;` discarded it, and the association stayed empty **with no
message**. The memo-table idiom `$cache[key, "field"] = v` — which the special
stage uses for its per-integrand analyses — depends on the auto-vivification.

### F5. `CoefficientArrays` was not implemented  (FIXED, v0.238)

```
CoefficientArrays[{x + 2 y - 3, 3 x - y}, {x, y}]   (* was unevaluated *)
```

`{b, M} = Normal[CoefficientArrays[eqs, vars]]` is how the stage's `LinSolveZero`
(its SymPy-`linsolve` emulation) assembles every linear system. Mathilda's arrays
are dense `List`s rather than `SparseArray`s — a deliberate divergence, and
harmless to that spelling because `Normal` of a `List` is the identity.

### F6. `BooleanQ`, `SymbolName`, `Internal`SyntacticNegativeQ` were missing  (FIXED, v0.238)

`SymbolName` is used by the stage's generator test
(`StringMatchQ[SymbolName[#], "Y" ~~ ___]`); `Internal`SyntacticNegativeQ` backs
its port of SymPy's `could_extract_minus_sign`; `BooleanQ` is used by the stress
harness to tell an `{answer, verified}` pair from a status list. `Hash` is still
not implemented but needs no fix: it stays unevaluated, which is a perfectly
deterministic `Association` key, and that is exactly what the existing Part II
memos rely on.

### F7. `AbsoluteTime[]` has INTEGER-SECOND resolution

```
t0 = AbsoluteTime[]; Do[Integrate[x^3 + 1/x, x], {i, 60}];
AbsoluteTime[] - t0                        (* 0.0;  the work took 4.9 ms      *)
First[AbsoluteTiming[Do[Integrate[x^3 + 1/x, x], {i, 60}]]]   (* 0.004941    *)
```

Mathematica's `AbsoluteTime[]` carries sub-millisecond fractional seconds, so the
idiom `t0 = AbsoluteTime[]; body; AbsoluteTime[] - t0` is a normal way to time a
region. Here it reads **0.0** for anything under a second, which is not an error
and not a message -- just a silently useless number. The research stress runner
times each case exactly that way (`stress_wl.py`'s `wltime`), so the Mathilda
runner had to switch to `AbsoluteTiming`, which is sub-millisecond. Worth fixing
in the kernel, since the failure mode is a plausible-looking zero.

### F8. The two contexts of the mixed-tower packages are easy to confuse

Not a kernel divergence but the same class of silent failure, recorded because it
cost two debugging cycles. Part II exports **seven** public symbols -- `Tower`,
`TowerD`, `ClassifyPrime`, `CertifyNonconstant`, `ParallelIntegrateMixed`,
`BuildTower`, `Undecided` -- and keeps ~200 more in `ParallelMixed`Private``. The
special stage is loaded *into* that private context, so its eight entry points are
private too. Addressing either with the wrong context leaves the call
**unevaluated**, and an unevaluated call reads downstream as a wrong *answer*
(`r[[1]] === b` is False, reported as "differs from Part II") rather than as an
error. When a call into these packages returns something structurally odd, check
the context before the algorithm.

### F9. The Automatic integer-factorisation budget is marginal, and perturbation-sensitive

Not caused by this port, but exposed by it, and worth a entry because the failure
is invisible in the answer.

`FactorInteger`'s Automatic method bounds its search (`factint_warn_incomplete`,
`src/facint.c`), and for a ~29-35 digit semiprime the bound is *marginal*: whether
a given number factors depends on unrelated process state. Two checked-in tests
sit on opposite sides of that line and **have never both been green**:

| build | `MoebiusMu[10^50 + 1]` (wants -1) | `PrimeNu[2491...238]` (wants 8) |
|---|---|---|
| v0.237 (before this port) | **-1** ok | 7 — `nofac` |
| v0.238 | **-1** ok | 7 — `nofac` |
| v0.239 | 1 — `nofac` | **8** ok |

So `moebiusmu_tests` and `primenu_tests` trade places; v0.239 did not break
factorisation, it moved which case lands on the failing side. The state
dependence is direct and reproducible **within one build**:

```
MoebiusMu[10^50 + 1]                                (* 1, with FactorInteger::nofac *)
FactorInteger[10^50 + 1]; MoebiusMu[10^50 + 1]      (* -1: correct, after a warm-up  *)
```

Both go through the same `internal_factorinteger`, so running `FactorInteger`
first changes the outcome of the second call.

**Why the bound is marginal.** In the Automatic path `pollard_rho_brent_mpz` uses
`max_iters = 14`, and `r` doubles per iteration, so one `(y_start, c)` attempt
covers ~2^14 = 16 K inner steps and the 98 attempts ~1.6 M. Pollard-rho needs
about `sqrt(p)` steps, and the unfactored cofactor here,
`23702464296258769770591950101`, is `14103673319201 * 1680588011350901` — its
smaller factor needs ~3.7 M. ECM then gets 7 B1 bounds x 10 curves, which is also
marginal for a 14-digit factor, and its curve parameter is random per call.

**The fix is a budget decision, not a bug fix**, which is why it is recorded here
rather than changed in a commit about integration: raising the rho budget (say
`max_iters` scaled to `size(n)` so a 29-digit semiprime gets its ~4 M steps) buys
completeness at the cost of time on every genuinely hard composite, and the whole
purpose of the present bound is to keep `FactorInteger` from hanging. Whoever
takes it should fix both tests at once and measure the cost on a deliberately hard
input (a 40-digit semiprime), not just on these two.

### F10. A context-injected `.m` prefixes its functions but NOT its variables

The sharp edge behind F8, found while scoping a budget from C and worth its own
entry because the wrong guess is **silent and self-confirming**.

`ParallelMixedSpecial.m` and `mixed/logrewrite.m` are loaded *into*
`ParallelMixed`Private`` (no `BeginPackage`; a `Begin[...]` inside the file) so
the body can name Part II's private functions by short name — the A22
workaround. That injection is asymmetric:

| defined in the file | ends up in |
|---|---|
| functions (`IntegrateSurfaceSpecial`, `BuildTower`, …) | `ParallelMixed`Private`` |
| top-level variables (`$SpecialTimeBudget`, `$StrictSP`, `$ParallelMixedTimeBudget`) | **`Global`** |

Measured at v0.240 after forcing the load:

```
$SpecialTimeBudget                        (* 45                          *)
ParallelMixed`Private`$SpecialTimeBudget  (* itself -- unbound           *)
Context[$SpecialTimeBudget]               (* Global`                     *)
Names["*SpecialTimeBudget*"]
  (* {$SpecialTimeBudget, ParallelMixed`Private`$SpecialTimeBudget} *)
```

Note the second element of that `Names` list: **the qualified spelling creates a
new symbol rather than failing**. So

```
Block[{ParallelMixed`Private`$SpecialTimeBudget = 3}, ...]
```

binds a fresh unbound symbol, leaves the real default in force, raises no
message, and the only symptom is that the budget does not apply — a 3 s budget
that still ran the full 45 s. The fix is to Block the **unqualified** name, which
is why the stress harnesses raise `$ParallelMixedTimeBudget` unqualified too. One
qualified name in the same expression resolving (the *function*) is not evidence
that the other (the *variable*) does.

### F11. The special stage shipped without a time budget (FIXED, v0.240)

Mathilda's own bug, not a divergence from Mathematica, recorded here because it
is the same shape as the rest of section F: a stage that cannot bound itself.

Part II has carried `$ParallelMixedTimeBudget = 45` since v0.161, with the reason
in its own comment — *"interactive Integrate has no timeout of its own"*. The
v0.239 special stage had **no budget at all**, and it is the *last* stage in the
Automatic cascade, so every integrand nothing else closed reached an unbounded
search. Measured: `Integrate[Sin[x^2 + Log[x]] Cos[x], x]` went from 47.6 s on a
pre-v0.238 binary to **97.3 s**, and `integrate_newton_leibniz_tests` from 50.5 s
(rc 0) to **120.06 s / SIGALRM**.

Fixed by `$SpecialTimeBudget = 45` on both entry points (parity with Part II, and
the value for an explicit `Method ->`) plus a 10 s scope for the Automatic
cascade — see the v0.240 changelog entry for how 10 s was chosen from the
measured cost of every close the stage brings to the cascade.

The general lesson, already in `tasks/lessons.md` in its DSolve form: **a
last-resort stage added to a cascade is charged to every input the cascade fails
to close**, so its worst case, not its typical case, is what the cascade pays.

### F12. Four false non-elementarity certificates from dependent tower generators (FIXED, v0.241)

Mathilda's own bug again, and the worst class there is: a certificate that the
integrand has **no** elementary antiderivative, issued for an integrand that has
one. Everything else in this file is a missing or slow answer; this is a wrong one,
and it is unfalsifiable from inside — the user has no way to tell a proof from a
mistake.

`Options[ParallelIntegrateMixed]` defaulted `"StructureTheorem" -> False`, so
`BuildTower` would raise a tower whose generators are algebraically **dependent**:
`E^x` alongside `E^(2x)`, `Log[x]` alongside `Log[x^2]`. Over such a tower every
residue argument is vacuous — a residue that looks non-constant in one generator is
constant once the dependency is used — so the Proposition 9.2(a) residue path
certified four elementary integrands:

| integrand | Mathilda said | the answer |
|---|---|---|
| `Exp[2x]/(1 + Exp[x])` | `{"not elementary", 1 + E^x, -E^(2x)}` | `E^x - Log[1 + E^x]` |
| `Exp[x]/(1 + Exp[x] + Exp[2x])` | `{"not elementary", …}` | `-2 ArcTan[(-1 - 2 E^x)/Sqrt[3]]/Sqrt[3]` |
| `1/(Exp[x] - Exp[-x])` | `{"not elementary", E^x - E^(-x), E^x/2}` | `(Log[-1 + E^x] - Log[1 + E^x])/2` |
| `Log[x^2]/Log[x]` | `{"not elementary", Log[x], x Log[x^2]}` | `2 x` |

Each also emitted `Integrate::nonelem` through the cascade. Present since the port;
confirmed pre-existing by reproducing all four on a pristine `HEAD` in a
`git worktree` rather than by arguing from the payload shape. The special-function
stage was immune throughout because its entry points always passed
`"StructureTheorem" -> True` explicitly — which is exactly the fix: the default is
now `True` for `ParallelIntegrateMixed` too. `Options[BuildTower]` keeps `False`,
since its direct callers pass the option themselves.

Cost of the default being on: the PMT 114-case corpus went 85 → 94 solved with zero
wrong answers, so the structure theorem *buys* coverage as well as soundness.

The lesson, in `tasks/lessons.md`: a soundness-critical option must default to the
safe side, and "the caller that matters always passes it" is not a defence — it
means only that the unsafe default is reached by every *other* caller.

### Not a divergence after all: `RowReduce[..., ZeroTest -> f]`

Listed as a seventh gap on the strength of

```
RowReduce[{{1, 1/(x+1) - 1/(x+1)}, {0, 1}}, ZeroTest -> (False &)]
```

returning the same reduced matrix as the default. The entry `1/(x+1) - 1/(x+1)`
evaluates to `0` **before** `RowReduce` sees it, so the option had nothing to
decide. Re-probed with `Sin[x]^2 + Cos[x]^2 - 1`, which stays symbolic: the
default reduces to the identity and the zero test correctly collapses the row.
`matsol_parse_zerotest_option` has been there all along.
<!-- /special-stage-port -->
