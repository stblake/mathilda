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

## D1: PolynomialGCD returns 1 for a non-trivial gcd (wrong answer)

    d1 = 7 + 5 Sqrt[2] - 4 x q;   d2 = 3 - Sqrt[2] x + q^3;
    PolynomialGCD[Expand[d1 d2], Expand[d1^2]]                         (* 1 *)
    PolynomialGCD[Expand[d1 d2], Expand[d1^2], Extension -> Automatic] (* 1 *)

Correct: `d1`, up to a unit. After `Expand`, the products mix `q`, `q^3`,
`(3/4 - 1/Sqrt[2])^2` and `Sqrt[2]`, atoms tied by relations the indeterminate
model does not know. The gcd over that model is 1, but over `K` it
is not. A one-level radical works (`PolynomialGCD[Expand[(x - Sqrt[2])(x + 3)],
Expand[(x - Sqrt[2])(x + 5)]]` is `-Sqrt[2] + x`; over `Sqrt[3 - Sqrt[2]]` the
answer is right up to a unit). `Extension -> Automatic` does not detect the
nested tower.

Fix direction: under `Extension -> Automatic` (and as a certificate otherwise),
collect the radical atoms, build the tower with `ToNumberField`/qqbar (where
`q` collapses into `Q(Sqrt[2])`), and run the gcd there. A gcd of 1 over
algebraic atoms should be certified, not assumed.

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
