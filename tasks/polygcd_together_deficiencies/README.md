# PolynomialGCD / Together over nested-radical coefficient fields

Found while fixing issue #133 (`Integrate[(x^2 + x + 1)/(x^4 + 1), x]` hung,
v0.6502). The Integrate hang itself is fixed at its root in the rational
integrator (`intrat_primitive` now divides polynomial content exactly, and
`intsimp_pos_sqrt` denests `Sqrt[a + b Sqrt[c]]`), so Integrate no longer
reaches these inputs. The two deficiencies below are still open in the algebra
layer. `repro.m` prints the current behaviour of each.

The field is `K = Q(I, Sqrt[2], q)`, `q = Sqrt[3/4 - 1/Sqrt[2]]`. Note that `q`
denests: `q = 1/Sqrt[2] - 1/2`, so `K = Q(I, Sqrt[2])`. Neither engine sees
that, because both treat `Sqrt[2]`, `q` and `q^3 = (3/4 - 1/Sqrt[2])^(3/2)` as
unrelated atoms, while `Expand` applies only some of their relations
(`q q -> 3/4 - 1/Sqrt[2]`, but `q q^3` stays `(3/4 - 1/Sqrt[2])^2`).

## D1: PolynomialGCD returned 1 for a non-trivial gcd (FIXED, v0.6505)

    d1 = 7 + 5 Sqrt[2] - 4 x q;   d2 = 3 - Sqrt[2] x + q^3;
    PolynomialGCD[Expand[d1 d2], Expand[d1^2]]                         (* -17/2 - 6 Sqrt[2] + x *)
    PolynomialGCD[Expand[d1 d2], Expand[d1^2], Extension -> Automatic] (* -17/2 - 6 Sqrt[2] + x *)

Both are `d1` up to a unit. Before the fix both answered `1`.

**Root cause.** The field-coefficient gcd engine (`flint_field_gcd` in
`src/poly/flint_bridge.c`) can already map a set of algebraic atoms into one
number field with qqbar, compute a modular gcd there, certify it, and render the
answer back in the caller's radicals. But it declined every univariate input
before looking at the atoms, on the grounds that the univariate engines were
already right. They are right only for independent atoms. `flint_numberfield_gcd`,
`flint_tower_gcd` and `flint_cyclotomic_gcd` take `Sqrt[n]` and roots of unity only,
so they declined here. The classical PRS that ran next treats `q`, `q^3`,
`(3/4 - 1/Sqrt[2])^2` and `Sqrt[2]` as independent indeterminates, and under that
model the operands are coprime.

**Fix.** A univariate input now goes to the field engine when its atoms are
algebraically dependent, that is, when the product of the atoms' degrees exceeds
`[K:Q]`. The verdict is cached with the field images, so a repeat costs a cache
hit. Inputs with one atom or with independent atoms keep their old route.
The engine's prime scan also re-factored `M mod p` for every scanned prime on
every call. For small operands that was about 3/4 of the cost, so split counts
are now memoised per minimal polynomial (`fg_split_count_memo`).

## D2: Together (and Simplify) did not finish (FIXED, v0.6507)

`d2_together_input.m` is a 4639-leaf rational function of `x` with coefficients
in `K`. It is the `(A d + B c)/g` quotient that `LogToAtan` handed to
`intrat_canonic` before the #133 fix. `Together` and `Simplify` both ran past
20 s, with or without `Extension -> Automatic`. Its value at `x = 1/3` is
`1.4714045207910316829`. It is in fact `1 + Sqrt[2] x`, and `x` never occurs in a
denominator: it is a polynomial in `x` with nested-radical numeric coefficients.
The copy used by the unit tests is `tests/together_nested_radical_d2.m`.

### Reference comparison (2026-11-11)

Each translation was checked numerically at `x = 1/3` before its result was
trusted. The Wolfram Language reads the file directly. Maxima and FriCAS got a
translation generated from the parsed expression (explicit `*`, `^`, `sqrt`, `%i`,
every subexpression parenthesised): Maxima's input evaluates to
`1.4714045207910162 - 4.0e-16 %i` in doubles, and FriCAS's normalised form to
`1.471404520791031682933896` at 200 digits. LeafCount is the Wolfram Language's
(for Maxima, the same leaf count of its tree; for FriCAS, the length of the
`unparse` string, since its output is not comparable term by term). Times are
wall clock on the same machine.

| System | Call | Result | Size | Time |
|---|---|---|---|---|
| Mathilda v0.6503 | `Together[e]`, `Together[e, Extension -> Automatic]`, `Simplify[e]` | did not finish | | > 20 s each |
| Mathilda v0.6505 (D1 gcd fix only) | `Together[e]`, `Together[e, Extension -> Automatic]` | did not finish | | > 20 s each |
| WL 14.0 | `Together[e]` | `(7 - 5 Sqrt[2] - 10 x + 7 Sqrt[2] x)/(Sqrt[3 - 2 Sqrt[2]] (-3 + 2 Sqrt[2]))` | LC 45 | 0.10 s |
| WL 14.0 | `Together[e, Extension -> Automatic]` | `((1 + Sqrt[2]) (7 - 5 Sqrt[2] - 10 x + 7 Sqrt[2] x))/(-3 + 2 Sqrt[2])` | LC 39 | 0.12 s |
| WL 14.0 | `Simplify[e]` | `(-7 + 5 Sqrt[2] + (10 - 7 Sqrt[2]) x)/(3 - 2 Sqrt[2])^(3/2)` | LC 34 | 0.22 s |
| WL 14.0 | `Together[RootReduce[e]]`, `FullSimplify[e]` | `1 + Sqrt[2] x` | LC 9 | 0.22 s, 0.54 s |
| Maxima 5.49.0 | `ratsimp(e)` | radicals as kernels, ~140-digit integer coefficients | 161 leaves | 0.022 s |
| Maxima 5.49.0 | `fullratsimp(e)` | the same shape | 129 leaves | 0.040 s |
| Maxima 5.49.0 | `ratsimp(e)`, `algebraic: true` | `sqrt(3-2^(3/2))*((sqrt(2)+2)*x+sqrt(2)+1)` | 18 leaves | 0.0045 s |
| Maxima 5.49.0 | `ratsimp(sqrtdenest(e))`, `algebraic: true` | `sqrt(2)*x+1` | 6 leaves | 0.026 s |
| FriCAS 1.3.6 | `normalize(e)` (also `rootSimp`, `simplify`) | ratio of two linear forms in `sqrt((3 sqrt(2) - 4)/sqrt(2))`, `sqrt(2)`, `%i` with 37-digit rational coefficients | 1661 chars | 0.45 s (parse and normalise) |
| Mathilda v0.6507 | `Together[e]` | equal-valued form over `Sqrt[2]`, `q`, `I` | LC 163 | 0.09 s |
| Mathilda v0.6507 | `Together[e, Extension -> Automatic]` | `1 + Sqrt[2] x` | LC 9 | 0.009 s |
| Mathilda v0.6507 | `Cancel[e]`; `Cancel[e, Extension -> Automatic]` | as for Together | LC 163; LC 9 | 0.26 s; 0.007 s |
| Mathilda v0.6507 | `Simplify[e]`, `FullSimplify[e]` | `1 + Sqrt[2] x` | LC 9 | 0.08 s, 0.08 s |

Notes on the references:
- Maxima's `ratsimp` and `fullratsimp` without `algebraic: true` are exact, but
  evaluating them in doubles is meaningless (`1.708 - 1.317 %i` and `-3.0 %i` at
  `x = 1/3`): the coefficients have about 140 digits and cancel.
- FriCAS's answer is correct to 200 digits, but at `x = 1/3` its numerator and
  denominator are each about `5e-75`. FriCAS's own `complexNumeric(.., 60)` fails
  on it with "division by zero". It treats `sqrt((3 sqrt(2) - 4)/sqrt(2))` as a
  degree-2 kernel over `Q(sqrt(2))`, which it is not, since that radical lies in
  `Q(sqrt(2))`.
- None of the three denests by default. The Wolfram Language with
  `Extension -> Automatic` and Maxima with `sqrtdenest` do.

### Root cause

The atoms `Sqrt[3/4 - 1/Sqrt[2]]` and `(3/4 - 1/Sqrt[2])^(3/2)` are radicals
whose radicand is itself algebraic. The generator model in `Together` reads
every radical as an indeterminate tied to its radicand. For a nested radical that
relation has algebraic coefficients and here is reducible (`q = 1/Sqrt[2] - 1/2`),
so the model cannot see the relations among `I`, `Sqrt[2]`, `q` and `q^3`. The
gcds of `cancel_recursive -> PolynomialGCD -> pseudo_rem` then swelled without
bound. Replacing `q` by a plain symbol made the same `Together` take 0.24 s,
which located the defect in the treatment of the nested atom rather than in the
gcd engine.

### Fix (`src/poly/algfield.c`, wired in `src/rat.c` and Simplify)

Gated on one linear walk for a radical with a constant, non-rational algebraic
base (read through Plus / Times / Power only):
- **Default** (reference-consistent: no denesting, no rationalisation). The fast
  exact paths run first. Only when they decline, just before the classical
  fallback that hung, each such radical's base takes its Together form (so a
  radicand that is really a radical monomial goes back to the ordinary path), and
  each remaining nested radical `b^(p/q)` is an opaque generator `g = b^(1/m)`. The
  combination runs with `g` a fresh symbol, and numerator and denominator are
  reduced modulo `g^m - b`. A call with no symbolic variable (numbers and numeric
  kernels only) keeps its former path.
- **`Extension -> Automatic`**: over the number field `K` of the coefficients.
  Every atom becomes an `AlgebraicNumber` over one qqbar primitive element. A
  nested atom whose exact canonical form is a radical expression over `Q` (a
  denestable one) is not a field generator; its atoms are used instead, which
  keeps the primitive element small (`Q(I, Sqrt[2])`, degree 4). Together or
  Cancel runs over `K`. The answer is written over a basis of radical monomials in
  the atoms, simplest atoms first: the basis grows by `B g^j` until a block is
  dependent, so an atom already in the field of simpler ones drops out. The
  denominator is monic over `K` and the rational content is cleared jointly.
- **Simplify** takes that field normal form when it scores lower, then polishes it.

Interaction with D1 (v0.6505, `flint_field_gcd`): the two fixes are independent.
D1's gcd engine makes each univariate gcd over dependent atoms correct and fast,
but on this input `Together` still timed out after it, because the time was in the
pipeline around the ~5300 gcds (coefficient extraction, `Expand`, AlgebraicNumber
products, per-gcd conversions). The nested-radical path here never reaches those
gcds: the default reads `q` as one opaque generator, and `Extension -> Automatic`
moves the whole rational function into `K` once. `Cancel[Expand[d1 d2]/Expand[d1^2],
Extension -> Automatic]` gives `d2/d1` through the field engine.

Two existing expectations changed, both value-equal:
- `extension_auto_builtins_tests`: `Together[1/(5 + 2 Sqrt[6])^(3/2) + (Sqrt[2] -
  Sqrt[3])^3, Extension -> Automatic]` is now `0`, as in the reference (it used
  to skip the extension for nested radicals and equal the default answer).
- `intrat_tests`: `Integrate[1/(x^2 - 3 - 2 Sqrt[2]), x]` is now
  `(1 - Sqrt[2]) ArcTanh[(-1 + Sqrt[2]) x]` (was `(-1 + Sqrt[2]) ArcTanh[(1 -
  Sqrt[2]) x]`), because an inner `Together` now returns the reduced
  `-(4 x)/Sqrt[48 + 32 Sqrt[2]]`. It is the reference's
  `-(ArcTanh[x/Sqrt[3 + 2 Sqrt[2]]]/Sqrt[3 + 2 Sqrt[2]])` with the radical denested.

## D3: open — many unrelated nested radicals at once (stress, timeouts only)

After the D1 (v0.6505) and D2 (v0.6507) fixes I ran a randomised stress check.
It uses 30 random sums of 1 to 3 fractions whose coefficients mix up to three of
`Sqrt[2]`, `Sqrt[3]`, `I`, `Sqrt[3/4 - 1/Sqrt[2]]`, `Sqrt[1 + Sqrt[2]]`,
`Sqrt[3 + 2 Sqrt[2]]`, `Sqrt[5 - 2 Sqrt[6]]`, `(3/4 - 1/Sqrt[2])^(3/2)`,
`(2 + Sqrt[3])^(1/3)` and `2^(1/3)`. The inputs are in `d3_stress_inputs.m`.
Run them with `./Mathilda -file tasks/polygcd_together_deficiencies/d3_stress.m`.

Every result was checked numerically at 3 points (30 digits), and no wrong
answer was found. The table counts calls that finished within 10 s:

| Operation | v0.6506 | v0.6507 |
|---|---|---|
| `Together` | 12 / 30 | 20 / 30 |
| `Together`, `Extension -> Automatic` | 13 / 30 | 15 / 30 |
| `Cancel` | 17 / 30 | 25 / 30 |
| `Cancel`, `Extension -> Automatic` | 17 / 30 | 17 / 30 |
| `Simplify` | 1 / 30 | 1 / 30 |

The remaining timeouts all predate these fixes. Their coefficient fields have
degree in the hundreds, e.g. `Q(I, Sqrt[2], Sqrt[3], Sqrt[1 + Sqrt[2]],
(2 + Sqrt[3])^(1/3))`.

Directions to investigate:
- a field-degree cap with an early, cheap decline, so `Simplify` returns its
  input instead of timing out;
- a tower (per-generator) representation instead of one primitive element.
