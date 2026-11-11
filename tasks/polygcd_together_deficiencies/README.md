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

## D2: Together (and Simplify) do not finish

`d2_together_input.m` is a 4639-leaf rational function of `x` with coefficients
in `K`. It is the `(A d + B c)/g` quotient that `LogToAtan` handed to
`intrat_canonic` before the #133 fix. `Together` and `Simplify` both run past
20 s, with or without `Extension -> Automatic`. Its value at `x = 1/3` is
`1.4714045207910316829`, so the answer is small.

The stack is `Together -> cancel_recursive -> PolynomialGCD ->
poly_gcd_internal -> pseudo_rem`, where the pseudo-remainders over the
multi-atom ring (`I`, `Sqrt[2]`, `q`, `q^3`, nested reciprocals) blow up. Fix
direction: as for D1. Map the algebraic atoms to one primitive element of `K`
first (here degree 4) and run the combination over that number field.

**Status after the D1 fix (v0.6505): still open.** The univariate gcds inside
`Together` now reach the field engine: about 5300 field gcds in the 20 s budget,
against about 950 before. The engine's own stages take about 0.03 s, and both runs
still time out at 20 s. In a profile, the time goes to the rest of the pipeline:
`CoefficientList` / `get_all_coeffs_expanded`, `Expand` and `AlgebraicNumber`
products over the multi-atom ring, `is_zero_poly`, and the per-gcd substitution
into the field and back. The remaining fix is to move the whole rational function
into `K` once, rather than one gcd at a time.
