# Algebra

This tutorial is a hands-on tour of Mathilda's algebra — the manipulation of
polynomials, rational expressions, and equations. You will multiply expressions
out and factor them back, pull a polynomial apart to inspect its pieces, put
fractions over a common denominator and split them up again, simplify with and
without assumptions, and finally solve equations and systems. Everything here is
pure rewriting of exact expressions; the [next tutorial](06-calculus.md) builds
calculus on top of exactly this machinery.

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
against the denominator, leaving `1 + x`.

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

## Solving equations

`Solve` finds the values of an unknown that satisfy an equation. Equations use
`==` (double equals); a single `=` is assignment. Start with a quadratic and a
linear equation:

```mathematica
In[1]:= Solve[x^2 == 4, x]
Out[1]= {{x -> -2}, {x -> 2}}

In[2]:= Solve[2 x + 6 == 0, x]
Out[2]= {{x -> -3}}
```

The answer is a list of solutions, each one a list of replacement rules. You can
feed those rules straight into another expression with `/.` to substitute a
solution back in. `Solve` is happy with symbolic coefficients, returning the
formula for the root:

```mathematica
In[1]:= Solve[a x + b == 0, x]
Out[1]= {{x -> -b/a}}
```

Roots need not be rational or even real. Irrational roots come back in radical
form, and complex roots in terms of `I`:

```mathematica
In[1]:= Solve[x^2 - 2 == 0, x]
Out[1]= {{x -> -Sqrt[2]}, {x -> Sqrt[2]}}

In[2]:= Solve[x^2 + 1 == 0, x]
Out[2]= {{x -> -I}, {x -> I}}
```

Higher-degree polynomials are no obstacle when they factor nicely. This cubic has
three integer roots:

```mathematica
In[1]:= Solve[x^3 - 6 x^2 + 11 x - 6 == 0, x]
Out[1]= {{x -> 1}, {x -> 2}, {x -> 3}}
```

Finally, `Solve` handles small *linear systems*. Pass a list of equations and a
list of the variables to solve for:

```mathematica
In[1]:= Solve[{2 x + y == 5, x - y == 1}, {x, y}]
Out[1]= {{x -> 2, y -> 1}}
```

Adding the two equations eliminates `y` and gives `3x = 6`, so `x = 2` and
`y = 1` — exactly what `Solve` reports, as a single solution binding both
variables.

## Where to next

You can now expand and factor polynomials, dissect them with `Coefficient` and
`Collect`, divide them and probe their roots, reshape rational expressions with
`Together`/`Apart`/`Cancel`, simplify with and without assumptions, and solve
equations and linear systems. This is the algebraic foundation the rest of
Mathilda's symbolic mathematics rests on.

- **[6. Calculus](06-calculus.md)** — put this algebra to work: differentiate,
  integrate (leaning on `Apart` for rational functions), expand functions as
  power series, take limits, evaluate sums, and solve problems numerically.
- **[Function reference](../documentation/index.md)** — the full catalogue of
  built-in functions. The
  [arithmetic and algebra](../documentation/arithmetic-and-algebra/index.md) and
  [simplification](../documentation/simplification/index.md) sections cover every
  function used above, including the complete options for `Factor`, `Solve`, and
  `Simplify`.
