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
| A18 a `Do` iterator inside a package captures a caller's same-named symbol | **OPEN** | `Table` is capture-avoiding, `Do` is not (the A11 fix does not cover it); `.m` iterators on the assembly path renamed |
| A19 `FreeQ[e, Complex]` is True on a complex atom | **OPEN** | found 2026-09-28 porting logrewrite; `logrewrite.m` tests `_Complex` everywhere |
| A20 `ComplexExpand` writes a real nested radical as `Cos[Arg[...]]`, `Arg[1 - Sqrt[5]]` unevaluated | **OPEN** | `logrewrite.m` puts constants in rectangular form by rules (`RectPow`); `RRad` of the package inherits the hazard |
| A21 `CountRoots` not implemented | **OPEN** | `SturmCount` in `logrewrite.m` |
| A22 `$InputFileName`, `DirectoryName` not implemented | **OPEN** | a package cannot `Get` a sibling file; `LoadModule["mixed/logrewrite.m"]` instead |
| A23 `NumericQ[Root[...]]` is False; `ToRadicals` gives the Ferrari form for every quartic | **OPEN** | `RootRadicals` in `logrewrite.m` (biquadratic / palindromic forms, root picked numerically) |
| A24 `Can`'s field detour returned a `Dot[{}, Inverse[{}], {}]` coefficient | **OPEN** (state-dependent, no standalone repro) | `logrewrite.m` canonicalises with `CanRaw` |
| A18 (re-checked on v0.221) | still **OPEN** | `g[v_] := Module[{s = 0}, Do[s += v, {k, 2}]; s]; g[k]` gives 3 |
| B7 an outer `TimeConstrained` cannot interrupt an inner one | **OPEN** (behavioural) | the package budget cannot cut the rewrite's own `TimeConstrained` short |
| B1 `ToNumberField` non-canonical primitive element | **OPEN** (behavioural, by design) | `.m` reads whatever theta comes back |
| B2–B6 | behavioural; see each entry | mostly by-design / hard |

Items A19-A24 and B7 were found on 2026-09-28 (builds 0.221/0.222) while porting the real form of
the logarithmic part (`src/internal/mixed/logrewrite.m`); each has a one-line repro below and a
workaround in that module (section D).

Remaining core work: the three items found when the review corpus was re-run
on v0.216 (2026-09-27, late): **A16** (extension gcd against an unexpanded zero constant, the cause of
three FALSE non-elementary certificates on the raw run), **A17** (`OneStepRowReduction` on a parametric
matrix) and **A18** (`Do` iterators capture a caller's symbol; `Table` does not). Each has a one-line
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

### A18. A `Do` iterator inside a package body captures a caller's symbol of the same name (`Table` does not)

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

### A26. Multivariate `PolynomialGCD` with `AlgebraicNumber` coefficients is wrong

```
a = AlgebraicNumber[Sqrt[2], {0, 1}];
PolynomialGCD[Expand[(x + a y) (x + 1)], Expand[(x + a y) (x + 2)]]
  Mathematica   x + Sqrt[2] y            (on the radical form)
  Mathilda      (x + a y) (x + 2)        -- the second operand, not a common divisor
PolynomialGCD[x + a y, Expand[(x + a y) (x + 2)]]
  Mathilda      (x + a y) (x + 2)        -- a gcd of higher degree than an operand
```

Univariate is correct, and so is `Extension -> Automatic` since v0.229 (`autodetect_walk` no
longer mines the field label `theta` out of an `AlgebraicNumber`, which used to answer `1`).
`collect_variables` no longer enrols an `AlgebraicNumber` as a polynomial VARIABLE (v0.229 --
`Variables[a + x]` is `{x}` now, as in Mathematica), which removes the inconsistent-ring hazard
of two Q-linearly dependent "variables" (`Plus`/`Times` fold a rational scalar into the
coordinate vector, so `2 a` is a structurally distinct atom). What remains is the classical
path's content computation over `K`: with equal degree in the main variable the pseudo-remainder
vanishes and the loop returns the second operand. `Cancel`/`Together` of such a fraction decline
rather than return garbage. The fix is a native `flint_field_gcd` over `K[x_1..x_n]` -- every
component exists (`kx_ctx_init`, `kx_scalar`, `field_subst_tau`, `field_mpoly_to_expr`), none is
wired to `PolynomialGCD`. `logrewrite.m`'s `KGcd` keeps the radical + `Extension` route for
multivariate input because of this.

### A27. `SparseArray` is not implemented

```
Normal[SparseArray[{{1, 2} -> 5, {2, 1} -> 7}, {2, 3}]]
  Mathematica   {{0, 5, 0}, {7, 0, 0}}
  Mathilda      SparseArray[{{1, 2} -> 5, {2, 1} -> 7}, {2, 3}]   (* unevaluated *)
```

`Normal` of the unevaluated head hands the `SparseArray[...]` expression straight back, so a
caller that assembles a matrix this way gets a non-matrix with no message and no error -- the
downstream `RowReduce` then "solves" it. Found while rewriting the ansatz assembly, which builds
its augmented matrix row by row instead.

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

**Summary (50 of 50 run): 42 verified** (35 in the run of the morning of 2026-09-23 on
0.174; the fixes in between are in E.3).  On the integrals it solves Mathilda is fast: faster than
SymPy on 37 of 42, median ratio Mathilda/SymPy 0.33, 19.9 s against
27.4 s on those; faster than Mathematica on 15, slower than Maxima on all but one.  The
budget now interrupts every long computation (P4, A2, A3, A27, A40 decline at 300.0 s; no OS
kill).  Everything still wrong is on a rung whose constants are `Root` objects: the linear system
is assembled with MORE equations than Mathematica assembles from the same package (A2, A3, A35,
P8), or the same system has no solution (A40's first rung), or the assembly over the compositum
never finishes (P4, A40's nested split), or a primitive misbehaves (A27: `NullSpace`'s `Method`).
The method is not at fault: Mathematica running the same `.wl` solves every one of these in
0.4-2.2 s.

| id | Mathilda | s | SymPy s | Mathematica s | Maxima s | reason |
|---|---|---|---|---|---|---|
| P1 | ok | 0.18 | 0.49 | 0.10 | 0.16 |  |
| P2 | ok | 0.04 | 0.16 | 0.04 | 0.03 |  |
| P3 | ok | 0.29 | 0.61 | 0.18 | 0.11 |  |
| P4 | T (budget) | 300.67 | 5.52 | 1.90 | 1.70 | {"time budget exceeded"} |
| P5 | ok | 0.40 | 1.90 | 0.38 | 0.13 |  |
| P6 | ok | 0.10 | 0.40 | 0.13 | 0.06 |  |
| P7 | ok | 0.03 | 0.12 | 0.03 | 0.02 |  |
| P8 | F | 30.38 | 6.59 | 1.07 | 1.36 | {"no solution within bounds", {1, 0}} |
| P9 | ok | 0.07 | 0.32 | 0.06 | 0.05 |  |
| P10 | ok | 0.18 | 0.44 | 0.10 | 0.10 |  |
| A1 | ok | 5.41 | 0.98 | 0.77 | 3.46 |  |
| A2 | T (budget) | 300.00 | 3.28 | 0.93 | 0.73 | {"time budget exceeded"} |
| A3 | T (budget) | 300.00 | 2.85 | 0.97 | 0.75 | {"time budget exceeded"} |
| A4 | ok | 0.09 | 0.64 | 0.09 | 0.04 |  |
| A5 | ok | 0.33 | 0.71 | 0.30 | 0.15 |  |
| A6 | ok | 0.39 | 0.92 | 0.24 | 0.17 |  |
| A7 | ok | 0.17 | 0.56 | 0.19 | 0.12 |  |
| A8 | ok | 0.08 | 0.27 | 0.08 | 0.06 |  |
| A9 | ok | 0.01 | 0.14 | 0.01 | 0.02 |  |
| A10 | ok | 0.04 | 0.11 | 0.04 | 0.03 |  |
| A11 | ok | 0.38 | 0.79 | 0.25 | 0.16 |  |
| A12 | ok | 0.34 | 1.05 | 0.28 | 0.16 |  |
| A13 | ok | 0.58 | 0.35 | 0.16 | 0.43 |  |
| A14 | ok | 0.10 | 0.29 | 0.10 | 0.07 |  |
| A15 | ok | 0.04 | 0.08 | 0.03 | 0.03 |  |
| A16 | ok | 0.91 | 1.17 | 0.46 | 0.34 |  |
| A17 | ok | 0.08 | 0.24 | 0.08 | 0.07 |  |
| A18 | ok | 0.10 | 0.24 | 0.07 | 0.08 |  |
| A19 | ok | 2.24 | 0.96 | 0.46 | 0.11 |  |
| A20 | ok | 0.11 | 0.52 | 0.18 | 0.09 |  |
| A21 | ok | 0.66 | 0.60 | 0.30 | 0.08 |  |
| A22 | ok | 0.11 | 0.35 | 0.17 | 0.04 |  |
| A23 | ok | 0.08 | 0.24 | 0.09 | 0.06 |  |
| A24 | ok | 0.06 | 0.22 | 0.08 | 0.05 |  |
| A25 | ok | 0.10 | 0.28 | 0.09 | 0.06 |  |
| A26 | ok | 0.06 | 0.18 | 0.08 | 0.05 |  |
| A27 | T (budget) | 300.00 | 3.64 | 0.41 | 0.09 | {"time budget exceeded"} |
| A28 | ok | 4.77 | 3.68 | 0.58 | 0.22 |  |
| A29 | ok | 0.29 | 2.83 | 0.31 | 0.20 |  |
| A30 | ok | 0.07 | 0.26 | 0.06 | 0.04 |  |
| A31 | ok | 0.28 | 0.59 | 0.17 | 0.21 |  |
| A32 | ok | 0.11 | 0.45 | 0.17 | 0.03 |  |
| A33 | ok | 0.08 | 0.37 | 0.13 | 0.05 |  |
| A34 | ok | 0.08 | 0.45 | 0.15 | 0.04 |  |
| A35 | F | 231.98 | 3.07 | 0.37 | 0.15 | {"internal: solution does not verify", {1, 0}} |
| A36 | ok | 0.09 | 0.10 | 0.07 | 0.07 |  |
| A37 | ok | 0.26 | 1.99 | 0.31 | 0.17 |  |
| A38 | ok | 0.10 | 0.33 | 0.08 | 0.06 |  |
| A39 | F | 1.34 | 1.61 | 1.05 | 0.20 | {"no solution within bounds", {2, 2}} |
| A40 | T (budget) | 300.04 | 2.65 | 2.16 | 0.20 | {"time budget exceeded"} |

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
