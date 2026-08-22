# Diophantine equations

The [Solutions of equations tutorial](../07-solutions-of-equations.md) solves
equations over the complex numbers, the reals, and — briefly — the integers. A
**Diophantine equation** asks the integer question in earnest: find *all* the
integer (or all the positive-integer) solutions of a polynomial equation, or
prove there are none. That one extra word — "integer" — changes the character of
the problem completely. `x² = 2` has the two real roots \(\pm\sqrt2\) but **no**
integer solutions; `x + 2y = 5` has infinitely many; `x⁵ + y⁵ + z⁵ + w⁵ = r⁵`
has a first solution that eluded everyone until a 1966 computer search. There is
no single formula: integer solving is a *dispatch* of specialised decision
procedures, one per shape of equation.

Mathilda's engine is the `Integers` domain of `Solve`:

```mathematica
In[1]:= Solve[x^2 + y^2 == z^2 && x + y + z == 12 && 0 < x < y && z > 0, {x, y, z}, Integers]
Out[1]= {{x -> 3, y -> 4, z -> 5}}
```

You write the equation with `==`, add any `&&` constraints (bounds, positivity,
`!=`), list the unknowns, and pass `Integers` as the domain. Behind that uniform
surface, `Solve` classifies the equation and routes it to the right method: the
extended Euclidean algorithm for a linear equation, continued fractions for a
Pell equation, Legendre's theorem for a ternary quadratic, the Tzanakis–de Weger
method for a Thue equation, meet-in-the-middle for a power-sum search, and about
a dozen more. These pages open up that dispatch, one family at a time.

---

## The two guarantees

Integer solving is unusual in that "no solution" is a *mathematical claim* as
strong as any answer — Fermat's Last Theorem is the statement that one particular
equation returns the empty set. Mathilda therefore holds itself to two
invariants, and every method on these pages obeys them:

1. **An empty result `{}` is always a proof.** Mathilda never returns `{}`
   because it gave up or ran out of a search window. A `{}` means *provably no
   solution* — an exhausted finite search, or a divisibility, congruence, or
   descent argument. When it cannot prove emptiness, it declines instead.
2. **A decline is always safe; a wrong answer is not.** For an equation outside
   the engine's reach, `Solve` returns the input **unevaluated** rather than
   guess. That is an honest "I cannot decide this", never a silently incomplete
   list. You will see this happen deliberately in the examples below.

Together these mean you can trust a concrete answer, trust an empty answer, and
recognise a gap on sight (the expression comes back unchanged).

---

## The families

<div class="grid cards" markdown>

-   :material-vector-polyline: __[Linear equations and Pell](linear-and-pell.md)__

    Single linear equations and systems (the extended Euclidean staircase and
    Hermite Normal Form), positivity rays, huge-coefficient lattice search, and
    the Pell family \(x^2 - Dy^2 = N\) — including negative Pell — solved by
    continued fractions.

-   :material-vector-circle: __[Quadratic forms and conics](quadratic-forms.md)__

    Euler's prime-generating conic, factorable (Runge) forms, rotated ellipses,
    the **ternary quadratic** decided by Legendre's theorem, and the Pythagorean
    triangle of fixed perimeter.

-   :material-function-variant: __[Cubics: Mordell, three cubes, and Thue](cubics-and-thue.md)__

    Integral points on Mordell curves \(y^2 = x^3 + k\), Booker's sum-of-three-cubes
    search, and **Thue equations** \(F(x,y) = m\) via the Tzanakis–de Weger
    method over a real number field.

-   :material-function: __[Exponential Diophantine equations](exponential.md)__

    Equations with the unknown in the *exponent*: Catalan's equation and
    Mihăilescu's theorem, and the Ramanujan–Nagell equation \(2^n - 7 = x^2\).

-   :material-pi: __[Famous results and impossibility proofs](famous-power-sums.md)__

    The showpieces: the **Lander–Parkin** counterexample to Euler's conjecture,
    **Frye's** quartic, taxicab numbers, the Euler brick, and Fermat's Last
    Theorem returned as a proof.

-   :material-calculator-variant: __[How Mathilda compares](performance.md)__

    A head-to-head against **sympy** and **PARI/GP** — where Mathilda wins on
    coverage and speed, where it is sound but slower, and the honest gaps that
    remain.

</div>

!!! tip "Following along"
    Start the REPL with `./Mathilda` and type each `In[...]` line yourself
    (without the prompt). Every `Out[...]` in these pages was produced by the
    actual binary. Solution sets are printed in Mathilda's canonical order, which
    may differ from a textbook, but the set is always complete and correct. Large
    sets are shown through `Length[...]` to keep the page readable — drop the
    wrapper to see every tuple.
