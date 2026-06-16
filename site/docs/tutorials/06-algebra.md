# Algebra

This tutorial is a hands-on tour of Mathilda's algebra — the manipulation and
simplification of polynomials and rational expressions. You will multiply
expressions out and factor them back, pull a polynomial apart to inspect its
pieces, reorganise it into Horner or square-free form, run the polynomial
division/GCD/resultant toolkit, put fractions over a common denominator and
split them up again, and simplify with and without assumptions. Everything here
is pure rewriting of exact expressions. Equation solving is a topic in its own
right and lives in the [next tutorial](07-solutions-of-equations.md); the
[calculus tutorial](08-calculus.md) then builds on exactly this machinery.

Every transcript below was produced by the actual Mathilda binary. Type the
`In[...]` lines yourself (without the prompt) and you will see the same
`Out[...]`. Mathilda's output is sometimes ordered differently from a textbook —
it sorts terms into a canonical order, so constants and lower-degree terms tend
to come first — but it is always mathematically correct.

## Expand and Factor: two sides of one coin

`Expand` multiplies everything out into a flat sum of terms. `Factor` runs the
other way, collapsing a polynomial back into a product. They are inverses, and
chaining one into the other is a satisfying way to see both at work.

```mathematica
In[1]:= Expand[(x + 2)(x - 3)]
Out[1]= -6 - x + x^2

In[2]:= Expand[(1 + x)^5]
Out[2]= 1 + 5 x + 10 x^2 + 10 x^3 + 5 x^4 + x^5

In[3]:= Factor[1 + 5 x + 10 x^2 + 10 x^3 + 5 x^4 + x^5]
Out[3]= (1 + x)^5
```

`In[1]` multiplies a product of two binomials; `In[2]` raises a binomial to the
fifth power, producing the row of binomial coefficients `1, 5, 10, 10, 5, 1`; and
`In[3]` factors that same degree-five polynomial straight back into `(1 + x)^5`.

`Expand` handles several variables at once, distributing across every term:

```mathematica
In[1]:= Expand[(x + y)^3]
Out[1]= x^3 + 3 x^2 y + 3 x y^2 + y^3

In[2]:= Expand[(a + b + c)^2]
Out[2]= a^2 + 2 a b + b^2 + 2 a c + 2 b c + c^2
```

`Factor` is the more impressive of the pair, because finding factors is genuinely
hard. It recognises the standard factorisations and more:

```mathematica
In[1]:= Factor[x^2 - 5 x + 6]
Out[1]= (-3 + x) (-2 + x)

In[2]:= Factor[6 x^2 + 11 x + 3]
Out[2]= (1 + 3 x) (3 + 2 x)

In[3]:= Factor[x^4 - 1]
Out[3]= (-1 + x) (1 + x) (1 + x^2)

In[4]:= Factor[x^6 - 1]
Out[4]= (-1 + x) (1 + x) (1 + x + x^2) (1 - x + x^2)
```

The quadratic in `In[1]` splits into two linear factors; `In[2]` shows `Factor`
coping with a leading coefficient. The differences of powers in `In[3]` and
`In[4]` break into all of their cyclotomic pieces — note that `1 + x^2` is left
intact because it has no real roots.

By default `Factor` works over the rationals, so an irreducible-over-the-rationals
polynomial is returned unchanged. To factor over a larger number field, name the
algebraic number to adjoin with the `Extension` option:

```mathematica
In[1]:= Factor[x^2 - 2]
Out[1]= -2 + x^2

In[2]:= Factor[x^2 - 2, Extension -> Sqrt[2]]
Out[2]= (Sqrt[2] + x) (-Sqrt[2] + x)
```

`x^2 - 2` has no rational roots, so `In[1]` leaves it alone; adjoining `Sqrt[2]`
in `In[2]` lets `Factor` split it into `(x - √2)(x + √2)`.

`Factor` has a few specialised cousins. `FactorSquareFree` only separates out
*repeated* factors, leaving the square-free part as one block instead of
factoring it fully — cheaper than `Factor`, and exactly what you want before
square-free–based algorithms like partial fractions or symbolic integration:

```mathematica
In[1]:= FactorSquareFree[(x - 1)^2 (x + 1)^3]
Out[1]= (-1 + x)^2 (1 + x)^3

In[2]:= FactorSquareFree[x^4 + 2 x^3 + 2 x^2 + 2 x + 1]
Out[2]= (1 + x)^2 (1 + x^2)
```

`In[2]` finds the repeated linear factor `(1 + x)^2` but leaves the remaining
square-free `1 + x^2` un-factored. `FactorTerms` goes the other way and pulls out
only the overall *numeric* content, leaving the polynomial part alone, while
`FactorTermsList` returns that content split into integer and content parts
alongside the remaining polynomial:

```mathematica
In[1]:= FactorTerms[2 x^2 + 4 x + 6]
Out[1]= 2 (3 + 2 x + x^2)

In[2]:= FactorTermsList[2 x^2 + 4 x + 6, x]
Out[2]= {2, 1, 3 + 2 x + x^2}
```

`In[1]` factors the common `2` out front; `In[2]` reports the same content as the
first list element and the monic remainder as the last.

## Looking inside a polynomial

A polynomial is just an expression, so you can interrogate its structure. The
most direct questions are "what is the coefficient of this power?" and "what are
all the coefficients?":

```mathematica
In[1]:= Coefficient[x^2 + 3 x + 2, x]
Out[1]= 3

In[2]:= Coefficient[3 x^2 + 5 x + 7, x, 2]
Out[2]= 3

In[3]:= CoefficientList[(x + 1)^3, x]
Out[3]= {1, 3, 3, 1}
```

`Coefficient[poly, x]` returns the coefficient of the first power of `x` (`In[1]`),
while the three-argument form `Coefficient[poly, x, n]` picks out the coefficient
of `x^n` (`In[2]` reads off the leading `3`). `CoefficientList` returns *all* the
coefficients at once, lowest power first, so `In[3]` gives the cube's row of
coefficients as a list. These compose with `Expand` — to read the coefficient of
`x^3` in `(1 + x)^5`, expand first and then ask:

```mathematica
In[1]:= Coefficient[Expand[(1 + x)^5], x, 3]
Out[1]= 10
```

To discover which symbols even appear, use `Variables`:

```mathematica
In[1]:= Variables[x^2 + y z + 3]
Out[1]= {x, y, z}
```

The dual operation to expanding is *grouping*. `Collect[expr, x]` gathers all the
terms that share a power of `x`, folding their coefficients together:

```mathematica
In[1]:= Collect[x y + x z + a x, x]
Out[1]= x (a + y + z)

In[2]:= Collect[a x^2 + b x^2 + x, x]
Out[2]= x + (a + b) x^2
```

`In[1]` factors the common `x` out of three terms; `In[2]` keeps the powers of
`x` separate but merges the two `x^2` coefficients into `a + b`. `Collect` is the
tool of choice when you want a polynomial organised *by* one variable while
treating the others as parameters.

`PolynomialQ` is the predicate that decides whether an expression *is* a
polynomial in the given variable — the gatekeeper before you reach for any of the
tools above:

```mathematica
In[1]:= PolynomialQ[x^2 + 1, x]
Out[1]= True

In[2]:= PolynomialQ[1/x + 1, x]
Out[2]= False
```

`1/x` is `x^(-1)`, a negative power, so `In[2]` is not a polynomial in `x`. Two
final reshaping tools change a polynomial's *form* without changing its value.
`HornerForm` rewrites it in nested Horner form, the arrangement that evaluates in
the fewest multiplications, and `Decompose` finds a *functional* decomposition —
expressing the polynomial as `p(q(x))` for smaller `p` and `q`:

```mathematica
In[1]:= HornerForm[1 + 2 x + 3 x^2 + 4 x^3]
Out[1]= 1 + x (2 + x (3 + 4 x))

In[2]:= Decompose[x^4 + 2 x^2 + 3, x]
Out[2]= {3 + 2 x + x^2, x^2}
```

`In[1]` nests the four coefficients so the cubic costs only three multiplies to
evaluate. `In[2]` discovers that `x^4 + 2 x^2 + 3` is `(u^2 + 2 u + 3)` composed
with `u = x^2`; the returned list `{p, q}` reads outermost-first, so the original
polynomial is `p` applied to `q`.

## The polynomial toolkit

Polynomials support their own arithmetic. `PolynomialGCD` finds the greatest
common divisor of two polynomials, exactly as `GCD` does for integers, and
`PolynomialQuotient`/`PolynomialRemainder` perform division with remainder:

```mathematica
In[1]:= PolynomialGCD[x^2 - 1, x^2 - 3 x + 2]
Out[1]= -1 + x

In[2]:= PolynomialQuotient[x^2 - 1, x - 1, x]
Out[2]= 1 + x

In[3]:= PolynomialRemainder[x^2 + 1, x - 1, x]
Out[3]= 2
```

Both `x^2 - 1` and `x^2 - 3x + 2` share the factor `x - 1`, which `In[1]` finds.
`In[2]` divides `x^2 - 1` by `x - 1` exactly (quotient `x + 1`, no remainder),
while `In[3]` divides `x^2 + 1` by `x - 1` and reports the leftover `2` —
precisely the value of `x^2 + 1` at `x = 1`, as the remainder theorem promises.

Two further tools answer questions about *roots* without computing them.
`Resultant[p, q, x]` is zero exactly when `p` and `q` share a common root, and
`Discriminant[p, x]` is zero exactly when `p` has a repeated root:

```mathematica
In[1]:= Resultant[x^2 - 1, x - 1, x]
Out[1]= 0

In[2]:= Discriminant[x^2 + b x + c, x]
Out[2]= b^2 - 4 c
```

`In[1]` is zero because `x = 1` is a root of both polynomials. `In[2]` recovers
the familiar quadratic discriminant `b² − 4c` symbolically — the very quantity
whose sign tells you how many real roots `x² + bx + c` has.

The GCD comes with the rest of the divisibility family. `PolynomialLCM` is the
least common multiple, and `PolynomialExtendedGCD` returns not just the GCD but
the Bézout cofactors `{a, b}` with `a p + b q = gcd`:

```mathematica
In[1]:= PolynomialLCM[x^2 - 1, x^2 - 3 x + 2]
Out[1]= (1 + x) (2 - 3 x + x^2)

In[2]:= PolynomialExtendedGCD[x^2 - 1, x^2 - 3 x + 2, x]
Out[2]= {-1 + x, {1/3, -1/3}}
```

`In[1]` multiplies the two quadratics' shared factor in only once. In `In[2]` the
first element `-1 + x` is the GCD and `{1/3, -1/3}` are the cofactors: one third
of the first polynomial minus one third of the second is exactly `x - 1`.
`PolynomialMod` reduces every coefficient modulo an integer — the entry point to
arithmetic over finite fields:

```mathematica
In[1]:= PolynomialMod[x^2 + 5 x + 7, 3]
Out[1]= 1 + 2 x + x^2
```

`5 ≡ 2` and `7 ≡ 1` modulo `3`, so the coefficients collapse accordingly.
`IrreduciblePolynomialQ` asks whether a polynomial factors at all over the
rationals — the predicate underlying `Factor`:

```mathematica
In[1]:= IrreduciblePolynomialQ[x^2 + 1]
Out[1]= True

In[2]:= IrreduciblePolynomialQ[x^2 - 1]
Out[2]= False
```

`x^2 + 1` has no rational roots and stays whole; `x^2 - 1` splits into
`(x - 1)(x + 1)`. `Resultant` is the top of a whole ladder of *subresultants*.
`Subresultants` returns their leading coefficients (the principal subresultant
coefficients) and `SubresultantPolynomials` returns the polynomials themselves —
the same remainder sequence a Euclidean GCD walks, but computed without fractions:

```mathematica
In[1]:= Subresultants[x^3 + x + 1, x^2 + x + 1, x]
Out[1]= {3, 1, 1}

In[2]:= SubresultantPolynomials[x^3 + x + 1, x^2 + x + 1, x]
Out[2]= {3, 2 + x, 1 + x + x^2}
```

The first entry of `In[1]`, `3`, is the resultant itself; being non-zero, it
confirms the two cubics share no common root. `In[2]` lists the corresponding
subresultant polynomials, ending in the higher of the two inputs.

For *several* polynomials in *several* variables at once, `GroebnerBasis`
computes a Gröbner basis — a canonical generating set for the ideal the
polynomials define. It is the multivariate generalisation of both `PolynomialGCD`
(one variable) and Gaussian elimination (linear systems), and it is the engine
behind solving polynomial systems:

```mathematica
In[1]:= GroebnerBasis[{x^2 + y^2 - 1, x - y}, {x, y}]
Out[1]= {-1 + 2 y^2, x - y}
```

The basis `In[1]` is *triangular*: the first polynomial `2 y^2 - 1` involves only
`y`, so it can be solved on its own, and then `x - y` pins down `x`. We will lean
on this structure heavily in the applications below.

## Rational expressions

A ratio of polynomials is a *rational expression*, and Mathilda gives you precise
control over its shape without ever changing its value. `Together` puts a sum of
fractions over a single common denominator:

```mathematica
In[1]:= Together[1/x + 1/(x + 1)]
Out[1]= (1 + 2 x)/(x + x^2)

In[2]:= Together[a/b + c/d]
Out[2]= (b c + a d)/(b d)
```

`Apart` does the reverse, splitting one fraction into a sum of simpler partial
fractions — the decomposition you learned for integrating rational functions by
hand:

```mathematica
In[1]:= Apart[(2 x + 3)/((x + 1)(x + 2))]
Out[1]= 1/(1 + x) + 1/(2 + x)

In[2]:= Apart[(x^3 + x + 1)/(x^2 - 1)]
Out[2]= x + 1/2/(1 + x) + 3/2/(-1 + x)
```

`In[1]` breaks a proper fraction into two simple ones. `In[2]` first divides out
a polynomial part (`x`, because the numerator's degree exceeds the
denominator's) and then expresses the remaining proper fraction in partial
fractions.

`Cancel` removes common factors between numerator and denominator, and
`Numerator`/`Denominator` extract the two halves:

```mathematica
In[1]:= Cancel[(x^2 - 1)/(x - 1)]
Out[1]= 1 + x

In[2]:= Numerator[(x + 1)/(x - 1)]
Out[2]= 1 + x

In[3]:= Denominator[(x + 1)/(x - 1)]
Out[3]= -1 + x
```

In `In[1]`, the numerator factors as `(x - 1)(x + 1)` and the `x - 1` cancels
against the denominator, leaving `1 + x`. When you want to expand *one half* of a
fraction while leaving the other in factored form, reach for `ExpandNumerator`
and `ExpandDenominator`:

```mathematica
In[1]:= ExpandNumerator[(x + 1)^2/(x - 1)^2]
Out[1]= (1 + 2 x + x^2)/(-1 + x)^2

In[2]:= ExpandDenominator[(x + 1)^2/(x - 1)^2]
Out[2]= (1 + x)^2/(1 - 2 x + x^2)
```

`In[1]` multiplies out the numerator only; `In[2]` does the same for the
denominator only — finer control than `Expand`, which would flatten both at once.

A different normalisation applies to powers and logarithms. `PowerExpand` rewrites
`(a b)^n` as `a^n b^n` and `Log[a b]` as `Log[a] + Log[b]`, treating every base
as positive — a deliberately *un-conditional* expansion that is only valid under
positivity assumptions, so use it knowingly:

```mathematica
In[1]:= PowerExpand[Sqrt[x^2 y]]
Out[1]= x Sqrt[y]

In[2]:= PowerExpand[Log[a b]]
Out[2]= Log[a] + Log[b]
```

`In[1]` pulls `x` out of the radical (which assumes `x ≥ 0`), and `In[2]` splits
the logarithm of a product into a sum of logarithms.

## Simplify and assumptions

`Simplify` is the all-purpose tool: it tries a battery of transformations —
factoring, cancelling, expanding, applying identities — and keeps whichever
result is smallest. It subsumes much of what the previous sections did manually:

```mathematica
In[1]:= Simplify[(x^2 - 1)/(x - 1)]
Out[1]= 1 + x

In[2]:= Simplify[(x^3 - 1)/(x - 1)]
Out[2]= 1 + x + x^2

In[3]:= Simplify[Sin[x]^2 + Cos[x]^2]
Out[3]= 1
```

`In[1]` and `In[2]` cancel removable factors, and `In[3]` applies the
Pythagorean trigonometric identity — a simplification no amount of factoring
alone would find.

Some simplifications are only valid under extra conditions. For instance
`Sqrt[x^2]` equals `x` only when `x` is non-negative (otherwise it is `-x`).
Supply such facts as a second argument to `Simplify`, or wrap a whole block in
`Assuming`:

```mathematica
In[1]:= Simplify[Sqrt[x^2], x > 0]
Out[1]= x

In[2]:= Assuming[x > 0, Simplify[Sqrt[x^2]]]
Out[2]= x
```

Both say the same thing: *given* that `x > 0`, `Sqrt[x^2]` simplifies to `x`.
Without the assumption Mathilda leaves the expression alone, because the
simplification would be wrong for negative `x`.

## Putting the toolkit to work

The real power of the algebra toolkit shows up when several tools combine to crack
a problem that has no obvious direct command. Here are three.

### Implicitising a parametric curve with `Resultant`

A curve is often given *parametrically*, as `x = f(t)`, `y = g(t)`. To find the
single polynomial equation `F(x, y) = 0` that the points `(x, y)` satisfy — its
*implicit* equation — you must eliminate the parameter `t`. That is exactly what
`Resultant` does: the resultant of `x - f(t)` and `y - g(t)` with respect to `t`
vanishes precisely when both share a common `t`, i.e. for points actually on the
curve. Take the parametrisation `x = t^2 - 1`, `y = t^3 - t`:

```mathematica
In[1]:= Resultant[x - (t^2 - 1), y - (t^3 - t), t]
Out[1]= x^2 + x^3 - y^2
```

So every point of the curve satisfies `y^2 = x^3 + x^2` — the **nodal cubic**, a
classic curve that crosses itself at the origin. We started from two cubics in a
hidden parameter and recovered one clean equation in `x` and `y` alone, with no
trigonometry or substitution by hand.

### Solving a symmetric polynomial system with `GroebnerBasis`

Suppose three unknowns satisfy the power-sum conditions

```
x + y + z = 6,   x² + y² + z² = 14,   x³ + y³ + z³ = 36.
```

This is a genuinely nonlinear system. Computing a Gröbner basis with the
variables ordered so that `z` is eliminated last triangularises it:

```mathematica
In[1]:= GroebnerBasis[{x + y + z - 6, x^2 + y^2 + z^2 - 14, x^3 + y^3 + z^3 - 36}, {x, y, z}]
Out[1]= {-6 + 11 z - 6 z^2 + z^3, 11 - 6 y + y^2 - 6 z + y z + z^2, -6 + x + y + z}
```

The first basis element involves only `z`: it is `z^3 - 6 z^2 + 11 z - 6`. Factor
it and the roots fall out:

```mathematica
In[1]:= Factor[z^3 - 6 z^2 + 11 z - 6]
Out[1]= (-3 + z) (-2 + z) (-1 + z)
```

So `z ∈ {1, 2, 3}`, and by the symmetry of the original equations `{x, y, z}` is
just a permutation of `{1, 2, 3}`. The Gröbner basis turned an intimidating
nonlinear system into a single one-variable cubic — the same elimination strategy
Gaussian elimination uses for linear systems, generalised to polynomials.

### Building a minimal polynomial for a nested radical

How do you find a polynomial equation with integer coefficients satisfied by a
messy nested radical such as `√(2 + √3)`? `MinimalPolynomial` answers directly:

```mathematica
In[1]:= MinimalPolynomial[Sqrt[2 + Sqrt[3]], x]
Out[1]= 1 - 4 x^2 + x^4
```

So `√(2 + √3)` is a root of `x⁴ - 4 x² + 1`. We can confirm it by substituting the
radical back in and simplifying to zero:

```mathematica
In[1]:= Simplify[(Sqrt[2 + Sqrt[3]])^4 - 4 (Sqrt[2 + Sqrt[3]])^2 + 1]
Out[1]= 0
```

The expression collapses to `0`, proving the radical really does satisfy that
quartic — a degree-4 algebraic number captured exactly, with no decimals.

## Where to next

You can now expand and factor polynomials, dissect them with `Coefficient` and
`Collect`, reshape them into Horner and square-free form, run the full
division/GCD/resultant/Gröbner toolkit, reshape rational expressions with
`Together`/`Apart`/`Cancel`, and simplify with and without assumptions. This is
the algebraic foundation the rest of Mathilda's symbolic mathematics rests on.

- **[7. Solutions of equations](07-solutions-of-equations.md)** — turn this
  machinery on equations and systems: `Solve`, `Reduce`, `Roots`, and friends
  build directly on the factoring, resultant, and Gröbner tools above.
- **[8. Calculus](08-calculus.md)** — put this algebra to work: differentiate,
  integrate (leaning on `Apart` for rational functions), expand functions as
  power series, take limits, evaluate sums, and solve problems numerically.
- **[Function reference](../documentation/index.md)** — the full catalogue of
  built-in functions. The
  [arithmetic](../documentation/arithmetic/index.md),
  [algebra](../documentation/algebra/index.md), and
  [simplification](../documentation/simplification/index.md) sections cover every
  function used above, including the complete options for `Factor` and `Simplify`.
