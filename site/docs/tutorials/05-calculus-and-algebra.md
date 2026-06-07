# Calculus and algebra

This tutorial is a hands-on tour of the symbolic mathematics at the heart of
Mathilda. We start with algebra — expanding, factoring, simplifying, and solving
— and then climb into calculus: derivatives, integrals, power series, limits,
and sums. The arc is deliberate. Algebra teaches you how Mathilda rewrites
expressions, and calculus builds on exactly the same machinery.

Every transcript below was produced by the actual Mathilda binary. Type the
`In[...]` lines yourself and you will see the same `Out[...]`. Mathilda's output
form is sometimes arranged differently from a textbook (terms may be reordered,
`E^x` is written instead of `Exp[x]`), but it is always mathematically correct.

## Expand and Factor: two sides of one coin

`Expand` multiplies everything out into a flat sum of terms. `Factor` runs the
other way, collapsing a polynomial back into a product. They are inverses, and
chaining one into the other is a satisfying way to see both at work.

```mathematica
In[1]:= Expand[(x+1)^3]
Out[1]= 1 + 3 x + 3 x^2 + x^3

In[2]:= Factor[1 + 3 x + 3 x^2 + x^3]
Out[2]= (1 + x)^3
```

`Expand` turned the cube into the four-term binomial expansion, and `Factor`
recognised that same polynomial and rebuilt the original cube. Notice that
Mathilda writes the factor as `1 + x` rather than `x + 1`; it sorts terms into a
canonical order, so constants come first.

## Simplify: collapsing the obvious

`Simplify` tries a battery of transformations and keeps whichever result is
smallest. A classic case is a rational expression with a removable factor:

```mathematica
In[1]:= Simplify[(x^2-1)/(x-1)]
Out[1]= 1 + x
```

The numerator factors as `(x - 1)(x + 1)`, the `x - 1` cancels against the
denominator, and what remains is the linear polynomial `1 + x`.

## Together and Apart

These two functions manipulate the *shape* of a rational expression without
changing its value. `Together` puts a sum of fractions over a common
denominator:

```mathematica
In[1]:= Together[1/x + 1/y]
Out[1]= (x + y)/(x y)
```

`Apart` does the reverse, splitting a single fraction into a sum of simpler
partial fractions:

```mathematica
In[1]:= Apart[(x+2)/(x (x+1))]
Out[1]= 2/x - 1/(1 + x)
```

`Apart` is the engine behind integrating rational functions by hand — break the
fraction apart, then integrate each piece. We will lean on exactly that pattern
shortly.

## Solve

`Solve` finds the values of an unknown that make an equation true. Equations use
`==` (double equals); a single `=` is assignment. Here is a quadratic:

```mathematica
In[1]:= Solve[x^2 - 5 x + 6 == 0, x]
Out[1]= {{x -> 2}, {x -> 3}}
```

The answer is a list of solutions, each one a list of replacement rules. The
quadratic `x^2 - 5x + 6` factors as `(x - 2)(x - 3)`, so the roots are `2` and
`3` — exactly what `Solve` reports.

`Solve` handles small linear systems too. Pass a list of equations and a list of
the variables to solve for:

```mathematica
In[1]:= Solve[{x + y == 3, x - y == 1}, {x, y}]
Out[1]= {{x -> 2, y -> 1}}
```

Adding the two equations gives `2x = 4`, so `x = 2` and `y = 1`. You can feed
these rules straight into another expression with `/.` (ReplaceAll) to
substitute the solution back in.

## Differentiation with D

`D[expr, x]` differentiates `expr` with respect to `x`. Mathilda knows the
standard rules of calculus and applies them recursively, so you can throw fairly
involved expressions at it.

Start with the power rule:

```mathematica
In[1]:= D[x^3, x]
Out[1]= 3 x^2
```

The product rule appears automatically when you differentiate a product. Here
`x` times `Sin[x]`:

```mathematica
In[1]:= D[x Sin[x], x]
Out[1]= x Cos[x] + Sin[x]
```

Mathilda differentiated each factor in turn — `x` contributes `Sin[x]` and
`Sin[x]` contributes `x Cos[x]` — and summed them. The chain rule shows up when
the argument of a function is itself a function of `x`:

```mathematica
In[1]:= D[Sin[x^2], x]
Out[1]= 2 x Cos[x^2]
```

The outer derivative `Cos[x^2]` is multiplied by the inner derivative `2 x`. The
same logic handles the exponential:

```mathematica
In[1]:= D[Exp[x^2], x]
Out[1]= 2 x E^x^2
```

Note that Mathilda prints `Exp[x^2]` as `E^x^2` — `E` is the symbol for Euler's
number, and `Exp` is just shorthand for a power of `E`.

### Higher-order derivatives

To differentiate more than once, use the `{x, n}` form. This computes the
*n*-th derivative with respect to `x`:

```mathematica
In[1]:= D[x^4, {x, 2}]
Out[1]= 12 x^2
```

The first derivative of `x^4` is `4 x^3`, and differentiating again gives
`12 x^2`.

## Integration with Integrate

`Integrate[expr, x]` computes an indefinite integral (an antiderivative) with
respect to `x`. Mathilda omits the constant of integration, as is conventional
in symbolic systems. Integration is genuinely harder than differentiation —
there is no single mechanical rule — so Mathilda runs a cascade of methods.

A polynomial is the easy case:

```mathematica
In[1]:= Integrate[x^2, x]
Out[1]= 1/3 x^3
```

Two integrals that everyone meets early are worth committing to memory. The
reciprocal integrates to a logarithm:

```mathematica
In[1]:= Integrate[1/x, x]
Out[1]= Log[x]
```

and the reciprocal of `x^2 + 1` integrates to an inverse tangent:

```mathematica
In[1]:= Integrate[1/(x^2+1), x]
Out[1]= ArcTan[x]
```

For something that needs integration by parts, try `x Exp[x]`:

```mathematica
In[1]:= Integrate[x Exp[x], x]
Out[1]= -E^x + x E^x
```

You can check any integral by differentiating it: `D[-E^x + x E^x, x]` returns
`x E^x`, confirming the result. The same trick works on `x Cos[x]`, which
integration by parts turns into a mix of sine and cosine:

```mathematica
In[1]:= Integrate[x Cos[x], x]
Out[1]= 1 + Cos[x] + x Sin[x]
```

A quick reminder: `Integrate` in this build computes *indefinite* integrals.
Differentiate the answer if you want to verify it, and apply your own bounds by
hand when you need a definite value.

## Series expansions

`Series[f, {x, a, n}]` expands `f` as a power series in `x` around the point `a`,
keeping terms up to order `n`. The trailing `O[x]^...` records the order of the
first term that was dropped, so you always know how far the approximation is
trusted.

The exponential is the cleanest example — every coefficient is `1/k!`:

```mathematica
In[1]:= Series[Exp[x], {x, 0, 4}]
Out[1]= 1 + x + 1/2 x^2 + 1/6 x^3 + 1/24 x^4 + O[x]^5
```

Sine has only odd-power terms with alternating signs:

```mathematica
In[1]:= Series[Sin[x], {x, 0, 5}]
Out[1]= x - 1/6 x^3 + 1/120 x^5 + O[x]^6
```

And the geometric series `1/(1-x)` is the textbook `1 + x + x^2 + ...`:

```mathematica
In[1]:= Series[1/(1-x), {x, 0, 4}]
Out[1]= 1 + x + x^2 + x^3 + x^4 + O[x]^5
```

Series expansions are how Mathilda — and a real numerical library — approximates
functions like `Sin` and `Exp` near a point.

## Limits

`Limit[expr, x -> a]` evaluates the limiting value of `expr` as `x` approaches
`a`. The most famous limit in all of calculus is `Sin[x]/x` at the origin, which
is `0/0` if you just substitute but in fact tends to `1`:

```mathematica
In[1]:= Limit[Sin[x]/x, x -> 0]
Out[1]= 1
```

Limits at infinity are equally at home. Use the symbol `Infinity` as the target,
and Mathilda compares the growth rates of numerator and denominator:

```mathematica
In[1]:= Limit[(2 x + 1)/(x + 3), x -> Infinity]
Out[1]= 2
```

As `x` grows without bound the `+1` and `+3` become negligible, leaving the
ratio of leading coefficients, `2/1 = 2`.

## Sums

`Sum[expr, {k, lo, hi}]` adds up `expr` as the index `k` runs from `lo` to `hi`.
When the bounds are concrete numbers you get a number, but the real fun is a
*symbolic* upper bound — Mathilda returns a closed form.

The sum of the first `n` integers is the well-known Gauss formula:

```mathematica
In[1]:= Sum[k, {k, 1, n}]
Out[1]= 1/2 n (1 + n)
```

The sum of squares has its own closed form:

```mathematica
In[1]:= Sum[k^2, {k, 1, n}]
Out[1]= 1/6 n (1 + n) (1 + 2 n)
```

Geometric sums work too. Summing `2^k` from `0` to `n`:

```mathematica
In[1]:= Sum[2^k, {k, 0, n}]
Out[1]= -1 + 2^(1 + n)
```

which is the familiar `2^(n+1) - 1`. Mathilda will even do this symbolically in
the base. Summing `x^k` produces the general geometric-series formula:

```mathematica
In[1]:= Sum[x^k, {k, 0, n}]
Out[1]= -1/(-1 + x) + x^(1 + n)/(-1 + x)
```

Rearranged, this is `(x^(n+1) - 1)/(x - 1)` — the closed form for a finite
geometric series, derived entirely symbolically.

## Where to next

You have now touched every major pillar of Mathilda's symbolic engine: algebraic
manipulation, differentiation, integration, series, limits, and sums. Each of
these has more depth than a single tutorial can cover — extra options, edge
cases, and related functions.

- For the complete reference on any function used here, browse the
  [function documentation](../documentation/index.md). The
  [calculus](../documentation/calculus/index.md) and
  [arithmetic and algebra](../documentation/arithmetic-and-algebra/index.md)
  sections are the natural follow-ups.
- To continue the guided path, head back to the
  [tutorials index](../tutorials/index.md) and pick the next walkthrough.

A good habit while exploring: type `?Name` at the prompt (for example
`?Integrate`) to read a function's built-in help string without leaving the
REPL.
