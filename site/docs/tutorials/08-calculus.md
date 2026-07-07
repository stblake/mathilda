# Calculus

This tutorial climbs from algebra into calculus. We differentiate, integrate,
expand functions as power series, take limits, evaluate sums, and finish with
numerical root-finding and optimisation. All of it is built on the same
expression-rewriting engine you met in the [algebra tutorial](06-algebra.md) —
calculus in Mathilda is mostly a large, well-organised collection of rewrite
rules.

Every transcript below was produced by the actual Mathilda binary. Type the
`In[...]` lines yourself and you will see the same `Out[...]`. Mathilda's output
form is sometimes arranged differently from a textbook (terms may be reordered,
`E^x` is written instead of `Exp[x]`, and indefinite integrals may differ by a
constant), but it is always mathematically correct.

## Differentiation with D

`D[expr, x]` differentiates `expr` with respect to `x`. Mathilda knows the
standard rules of calculus and applies them recursively, so you can throw fairly
involved expressions at it. Start with the power rule, which works for symbolic
exponents too:

```mathematica
In[1]:= D[x^5, x]
Out[1]= 5 x^4

In[2]:= D[x^n, x]
Out[2]= n x^(-1 + n)
```

The product rule appears automatically when you differentiate a product, and the
chain rule when a function's argument is itself a function of `x`:

```mathematica
In[1]:= D[x Sin[x], x]
Out[1]= x Cos[x] + Sin[x]

In[2]:= D[Sin[x^2], x]
Out[2]= 2 x Cos[x^2]

In[3]:= D[Exp[3 x], x]
Out[3]= 3 E^(3 x)
```

In `In[1]` each factor is differentiated in turn — `x` contributes `Sin[x]` and
`Sin[x]` contributes `x Cos[x]` — and the results are summed. In `In[2]` the
outer derivative `Cos[x^2]` is multiplied by the inner derivative `2 x`; the same
logic produces the factor of `3` in `In[3]`.

The elementary functions all carry their derivatives:

```mathematica
In[1]:= D[Log[x], x]
Out[1]= 1/x

In[2]:= D[Tan[x], x]
Out[2]= Sec[x]^2

In[3]:= D[ArcTan[x], x]
Out[3]= 1/(1 + x^2)
```

### Higher-order derivatives

To differentiate more than once, use the `{x, n}` form, which computes the
*n*-th derivative with respect to `x`. Here is the fourth derivative of sine,
which — after four quarter-turns through cosine and back — returns to sine:

```mathematica
In[1]:= D[Sin[x], {x, 4}]
Out[1]= Sin[x]
```

### Partial derivatives

When an expression contains several variables, `D[expr, x]` differentiates with
respect to `x` and treats the others as constants:

```mathematica
In[1]:= D[x^2 y^3, x]
Out[1]= 2 x y^3

In[2]:= D[x^2 y^3, y]
Out[2]= 3 x^2 y^2
```

In `In[1]` the factor `y^3` rides along unchanged; in `In[2]` it is `x^2` that is
held constant. This is exactly the partial-derivative convention.

### Derivatives of unknown functions

Differentiating an *abstract* function — one with no definition — produces
`Derivative` notation rather than a concrete answer, which is how the chain and
product rules express themselves in general:

```mathematica
In[1]:= D[f[g[x]], x]
Out[1]= Derivative[1][g][x] Derivative[1][f][g[x]]

In[2]:= D[f[x], {x, 2}]
Out[2]= Derivative[2][f][x]
```

`Derivative[1][f]` is "the first derivative of `f`", and `Derivative[2][f]` the
second. `In[1]` is precisely the chain rule `f'(g(x)) · g'(x)` written in
Mathilda's notation.

### Total derivatives with Dt

`D` computes *partial* derivatives. Its cousin `Dt` computes *total*
differentials, treating every symbol as potentially depending on the
differentiation variable:

```mathematica
In[1]:= Dt[x y]
Out[1]= Dt[x] y + x Dt[y]

In[2]:= Dt[x^2, x]
Out[2]= 2 x
```

`In[1]` gives the total differential of a product — `y dx + x dy` in classical
notation. `In[2]` takes the total derivative with respect to `x`; with no other
dependencies declared, it agrees with `D`.

## Integration with Integrate

`Integrate[expr, x]` computes an indefinite integral (an antiderivative) with
respect to `x`. Mathilda omits the constant of integration, as is conventional in
symbolic systems — so two correct answers may differ by a constant. Integration
is genuinely harder than differentiation (there is no single mechanical rule), so
Mathilda runs a cascade of methods.

Polynomials and the basic functions are the easy cases:

```mathematica
In[1]:= Integrate[x^4, x]
Out[1]= 1/5 x^5

In[2]:= Integrate[1/x, x]
Out[2]= Log[x]

In[3]:= Integrate[Cos[x], x]
Out[3]= Sin[x]

In[4]:= Integrate[Exp[x], x]
Out[4]= E^x
```

Two more worth committing to memory: the reciprocal of `x^2 + 1` integrates to an
inverse tangent, and `1/(x^2 - 1)` to an inverse hyperbolic tangent:

```mathematica
In[1]:= Integrate[1/(x^2 + 1), x]
Out[1]= ArcTan[x]

In[2]:= Integrate[1/(x^2 - 1), x]
Out[2]= -ArcTanh[x]
```

For integrands that need *integration by parts*, Mathilda applies the technique
automatically. Two classic cases are `x e^x` and `x^2 log x`:

```mathematica
In[1]:= Integrate[x Exp[x], x]
Out[1]= -E^x + x E^x

In[2]:= Integrate[x^2 Log[x], x]
Out[2]= 1/9 (-x^3 + 3 x^3 Log[x])
```

You can always check an integral by differentiating it — antidifferentiation and
differentiation are inverses, so the derivative of the answer must recover the
integrand:

```mathematica
In[1]:= D[-E^x + x E^x, x]
Out[1]= x E^x
```

The derivative of `In[1]` of the previous block returns `x E^x`, confirming the
integration by parts.

A *rational function* is integrated by first splitting it into partial fractions
(the `Apart` of the algebra tutorial) and integrating each piece. Mathilda does
this for you:

```mathematica
In[1]:= Integrate[(2 x + 3)/((x + 1)(x + 2)), x]
Out[1]= Log[2 + 3 x + x^2]
```

Recall that `Apart[(2 x + 3)/((x + 1)(x + 2))]` is `1/(x + 1) + 1/(x + 2)`;
integrating those two pieces gives `Log[x + 1] + Log[x + 2]`, which combines into
the single logarithm Mathilda reports (its argument `2 + 3x + x^2` is just
`(x + 1)(x + 2)` expanded).

!!! note "Indefinite only"
    `Integrate` in this build computes *indefinite* integrals — it does not
    accept integration bounds. To evaluate a definite integral, find the
    antiderivative with `Integrate` and apply the limits yourself.

## Series expansions

`Series[f, {x, a, n}]` expands `f` as a power series in `x` around the point `a`,
keeping terms up to order `n`. The trailing `O[x]^...` records the order of the
first dropped term, so you always know how far the approximation is trusted.

The exponential is the cleanest example — every coefficient is `1/k!` — while
cosine keeps only even powers and `Log[1 + x]` has the alternating harmonic
coefficients:

```mathematica
In[1]:= Series[Exp[x], {x, 0, 5}]
Out[1]= 1 + x + 1/2 x^2 + 1/6 x^3 + 1/24 x^4 + 1/120 x^5 + O[x]^6

In[2]:= Series[Cos[x], {x, 0, 6}]
Out[2]= 1 - 1/2 x^2 + 1/24 x^4 - 1/720 x^6 + O[x]^7

In[3]:= Series[Log[1 + x], {x, 0, 4}]
Out[3]= x - 1/2 x^2 + 1/3 x^3 - 1/4 x^4 + O[x]^5
```

The geometric series and the arctangent series are textbook results, and Mathilda
even expands functions like `Tan` whose coefficients have no simple closed form:

```mathematica
In[1]:= Series[1/(1 - x), {x, 0, 5}]
Out[1]= 1 + x + x^2 + x^3 + x^4 + x^5 + O[x]^6

In[2]:= Series[ArcTan[x], {x, 0, 5}]
Out[2]= x - 1/3 x^3 + 1/5 x^5 + O[x]^6

In[3]:= Series[Tan[x], {x, 0, 5}]
Out[3]= x + 1/3 x^3 + 2/15 x^5 + O[x]^6
```

You can expand around a point other than the origin. The series is then written
in powers of `(x - a)`:

```mathematica
In[1]:= Series[Exp[x], {x, 1, 2}]
Out[1]= E + E (x - 1) + 1/2 E (x - 1)^2 + O[x - 1]^3
```

Every coefficient is a power of `E`, because all the derivatives of `e^x` at
`x = 1` equal `e`. To turn a series back into an ordinary polynomial — dropping
the `O[x]` bookkeeping term — apply `Normal`:

```mathematica
In[1]:= Normal[Series[Exp[x], {x, 0, 3}]]
Out[1]= 1 + x + 1/2 x^2 + 1/6 x^3
```

The result is now a plain expression you can differentiate, integrate, or
evaluate like any other.

## Limits

`Limit[expr, x -> a]` evaluates the limiting value of `expr` as `x` approaches
`a`. Its real value is computing *indeterminate* forms — expressions that are
`0/0` or `∞/∞` on direct substitution. The most famous is `Sin[x]/x` at the
origin:

```mathematica
In[1]:= Limit[Sin[x]/x, x -> 0]
Out[1]= 1

In[2]:= Limit[(1 - Cos[x])/x^2, x -> 0]
Out[2]= 1/2

In[3]:= Limit[(Exp[x] - 1)/x, x -> 0]
Out[3]= 1
```

All three are `0/0` if you just substitute, yet each has a finite limit that
`Limit` finds. Limits at infinity are equally at home — use the symbol `Infinity`
as the target:

```mathematica
In[1]:= Limit[(2 x^2 + 1)/(x^2 + x), x -> Infinity]
Out[1]= 2

In[2]:= Limit[(1 + 1/x)^x, x -> Infinity]
Out[2]= E
```

`In[1]` compares the growth rates of numerator and denominator: as `x` grows, the
lower-order terms become negligible and the ratio tends to that of the leading
coefficients, `2`. `In[2]` is the celebrated limit that *defines* Euler's number.

When a limit is genuinely infinite, Mathilda says so. A real-signed blow-up gives
`Infinity`; a blow-up with no single direction gives `ComplexInfinity`:

```mathematica
In[1]:= Limit[Log[x], x -> 0]
Out[1]= -Infinity

In[2]:= Limit[1/x, x -> 0]
Out[2]= ComplexInfinity
```

`Log[x]` decreases without bound as `x -> 0`, so `In[1]` is `-Infinity`. `1/x`
runs to `+∞` from the right but `-∞` from the left, so it has no ordinary
two-sided limit — Mathilda reports `ComplexInfinity` to signal exactly that.

### Choosing a strategy

Under the hood `Limit` tries a cascade of strategies — direct substitution, a
rational-function shortcut, series expansion, L'Hôpital's rule, asymptotic
reductions, and squeeze arguments — stopping at the first that succeeds. The
`Method` option lets you pin a single one. `Method -> Automatic` (the default)
runs the whole cascade; a named method runs *only* that strategy, and leaves the
`Limit` unevaluated if it does not apply:

```mathematica
In[1]:= Limit[Sin[x]/x, x -> 0, Method -> "Series"]
Out[1]= 1

In[2]:= Limit[(2 x^2 + 1)/(x^2 + x), x -> Infinity, Method -> "RationalFunction"]
Out[2]= 2

In[3]:= Limit[Sin[x]/x, x -> 0, Method -> "RationalFunction"]
Out[3]= Limit[Sin[x]/x, x -> 0, Method -> "RationalFunction"]
```

`In[1]` forces the series-expansion route; `In[2]` uses the leading-degree
comparison for rational functions. `In[3]` asks for that same rational shortcut
on a transcendental ratio it cannot handle, so the limit comes back untouched.
The available methods are `"Substitution"`, `"RationalFunction"`, `"Series"`,
`"LHospital"`, `"Asymptotic"` and `"Bounded"` (see the [`Limit`](../documentation/calculus/Limit.md)
reference for what each covers). Pinning a method is mainly useful for teaching
and for isolating which strategy handles a given limit; for everyday work,
`Automatic` is the right choice.

## Sums

`Sum[expr, {k, lo, hi}]` adds up `expr` as the index `k` runs from `lo` to `hi`.
With concrete numeric bounds you simply get a number:

```mathematica
In[1]:= Sum[k^2, {k, 1, 10}]
Out[1]= 385
```

The real power, though, is a *symbolic* upper bound — Mathilda returns a closed
form. The sums of the first `n` integers, squares, and cubes all have classic
formulas:

```mathematica
In[1]:= Sum[k, {k, 1, n}]
Out[1]= 1/2 n (1 + n)

In[2]:= Sum[k^2, {k, 1, n}]
Out[2]= 1/6 n (1 + n) (1 + 2 n)

In[3]:= Sum[k^3, {k, 1, n}]
Out[3]= 1/4 n^2 (1 + n)^2
```

The last is the lovely identity that the sum of the first `n` cubes equals the
*square* of the sum of the first `n` integers. Geometric sums work too, even with
a symbolic base:

```mathematica
In[1]:= Sum[2^k, {k, 0, n}]
Out[1]= -1 + 2^(1 + n)

In[2]:= Sum[x^k, {k, 0, n}]
Out[2]= -1/(-1 + x) + x^(1 + n)/(-1 + x)
```

`In[1]` is the familiar `2^(n+1) - 1`. Rearranged, `In[2]` is `(x^(n+1) - 1)/(x - 1)`
— the closed form for a finite geometric series, derived entirely symbolically.

`Sum` also evaluates *infinite* sums when the upper bound is `Infinity` —
closing famous series like the Basel problem into constants:

```mathematica
In[1]:= Sum[1/k^2, {k, 1, Infinity}]
Out[1]= 1/6 Pi^2
```

!!! tip "Going further"
    Infinite summation and its companion, infinite products, are large enough
    topics to have their own tutorials: [symbolic summation](11-symbolic-summation.md)
    (telescoping, the zeta family, Euler sums, `π`-machines) and
    [infinite products](12-infinite-products.md) (Wallis, Viète, Euler prime
    products, and the Gamma/Glaisher constants).

## Numerical calculus

Not every equation has a closed-form solution, and not every integral or root can
be written down. For those, Mathilda offers numerical methods that return
machine-precision answers (see the
[precision tutorial](04-arithmetic.md) for what
those decimals mean).

`FindRoot[lhs == rhs, {x, x0}]` searches for a numerical root, starting from your
initial guess `x0`:

```mathematica
In[1]:= FindRoot[Cos[x] == x, {x, 1}]
Out[1]= {x -> 0.739085}

In[2]:= FindRoot[x^2 - 2 == 0, {x, 1}]
Out[2]= {x -> 1.41421}
```

`In[1]` finds the unique solution of `cos x = x`, a transcendental equation with
no algebraic answer. `In[2]` converges to `√2 ≈ 1.41421` — the same number
`Solve[x^2 - 2 == 0, x]` returns exactly as `Sqrt[2]`, here computed numerically
from the guess `x = 1`.

`FindMinimum` and `FindMaximum` locate a local extremum, returning both the
optimal value and the location that achieves it:

```mathematica
In[1]:= FindMinimum[x^2 - 4 x + 7, {x, 0}]
Out[1]= {3.0, {x -> 2.0}}

In[2]:= FindMaximum[-x^2 + 2 x + 1, {x, 0}]
Out[2]= {2.0, {x -> 1.0}}

In[3]:= FindMinimum[Sin[x], {x, 4}]
Out[3]= {-1.0, {x -> 4.71239}}
```

The parabola in `In[1]` has its minimum value `3` at `x = 2` (its vertex), and
the downward parabola in `In[2]` peaks at `2` when `x = 1`. `In[3]` starts near
`x = 4` and slides downhill to the nearest trough of sine at `x ≈ 4.71239`
(that is `3π/2`), where the value is `-1`. Because the search is local, the
starting point decides *which* extremum you find.

## Where to next

You have now toured every major pillar of Mathilda's calculus: differentiation
(partial, higher-order, and of abstract functions), integration by several
methods, power series, limits, symbolic sums, and numerical root-finding and
optimisation. Together with the [algebra tutorial](06-algebra.md) this covers the
heart of the system's symbolic mathematics.

- For the complete reference on any function used here, browse the
  [function documentation](../documentation/index.md). The
  [calculus](../documentation/calculus/index.md) and
  [power series](../documentation/power-series/index.md) sections are the natural
  follow-ups, with the full options for `D`, `Integrate`, `Series`, `Limit`,
  `Sum`, and `FindRoot`.
- To revisit the guided path from the start, head back to the
  [tutorials index](index.md).

A good habit while exploring: type `?Name` at the prompt (for example
`?Integrate`) to read a function's built-in help string without leaving the REPL.
