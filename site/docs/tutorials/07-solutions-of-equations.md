# Solutions of Equations

This tutorial is a hands-on tour of how Mathilda *solves* equations — from a
single linear equation through quadratics, higher-degree polynomials,
simultaneous systems, and on to transcendental equations whose solutions form
infinite families. It picks up where the [algebra tutorial](06-algebra.md)
left off: there you reshaped expressions, here you ask Mathilda to find the
values of an unknown that make an equation true. The [next
tutorial](08-calculus.md) reuses several of these tools — solving a
derivative equal to zero to locate maxima and minima, for instance.

Every transcript below was produced by the actual Mathilda binary. Type the
`In[...]` lines yourself (without the prompt) and you will see the same
`Out[...]`. The single most important thing to remember is the difference
between the two equals signs: equations use `==` (a double equals), while a
single `=` is *assignment*. Writing `x = 4` sets `x` to `4`; writing `x == 4`
is the equation you hand to `Solve`. Output is sometimes ordered differently
from a textbook — Mathilda sorts terms and solutions into a canonical order —
but it is always mathematically correct.

## Solving a single equation

`Solve[eqn, x]` finds the values of `x` that satisfy `eqn`. Start with a linear
equation and a quadratic:

```mathematica
In[1]:= Solve[2 x + 6 == 0, x]
Out[1]= {{x -> -3}}

In[2]:= Solve[x^2 == 4, x]
Out[2]= {{x -> -2}, {x -> 2}}
```

The answer is a **list of solutions**, and each solution is itself a list of
replacement rules of the form `x -> value`. The linear equation has one
solution, so `In[1]` is a list with one entry; the quadratic has two, so `In[2]`
lists both. This nested `{{...}}` shape is deliberate and uniform — even a
single solution comes wrapped in its own list — which makes the result easy to
process programmatically.

Because each solution is a rule, you can feed it straight back into an
expression with `/.` (`ReplaceAll`) to substitute the solution and check or
use it. Assign the solution list to a name and replace into another expression:

```mathematica
In[1]:= sol = Solve[x^2 == 4, x]
Out[1]= {{x -> -2}, {x -> 2}}

In[2]:= x^2 + 1 /. sol
Out[2]= {5, 5}
```

`In[2]` applies *both* solutions, returning a list with one value per solution;
since `(-2)^2 + 1` and `2^2 + 1` are both `5`, you get `{5, 5}`. This
substitute-back idiom is how you verify a solution, evaluate a quantity at a
solution, or chain solving into a larger computation.

`Solve` is perfectly happy with **symbolic coefficients**, returning the
formula for the root rather than a number. The general linear and quadratic
equations give back exactly the formulas you learned by hand:

```mathematica
In[1]:= Solve[a x + b == 0, x]
Out[1]= {{x -> -b/a}}

In[2]:= Solve[a x^2 + b x + c == 0, x]
Out[2]= {{x -> (1/2 (-b + Sqrt[b^2 - 4 a c]))/a}, {x -> (1/2 (-b - Sqrt[b^2 - 4 a c]))/a}}
```

`In[1]` is the solution of `ax + b = 0`, and `In[2]` is the quadratic formula
itself — the two roots `(-b ± √(b² − 4ac)) / 2a`, with the discriminant
`b² − 4ac` appearing under the radical.

Roots need not be rational or even real. Irrational roots come back in radical
form, and complex roots in terms of `I`, the imaginary unit:

```mathematica
In[1]:= Solve[x^2 - 2 == 0, x]
Out[1]= {{x -> -Sqrt[2]}, {x -> Sqrt[2]}}

In[2]:= Solve[x^2 + 1 == 0, x]
Out[2]= {{x -> -I}, {x -> I}}
```

`x² = 2` has the two irrational roots `±√2`, and `x² = −1` has the two
imaginary roots `±i`. Mathilda works with these exactly; nothing is rounded.

## Higher-degree polynomials

A polynomial equation has as many roots as its degree (counted with
multiplicity). When the polynomial factors over the rationals, `Solve` finds
every root exactly. This cubic has three integer roots, and so does the same
polynomial written as a product:

```mathematica
In[1]:= Solve[x^3 - 6 x^2 + 11 x - 6 == 0, x]
Out[1]= {{x -> 1}, {x -> 2}, {x -> 3}}

In[2]:= Solve[(x - 1)(x - 2)(x - 3) == 0, x]
Out[2]= {{x -> 1}, {x -> 2}, {x -> 3}}
```

A biquadratic (a quartic in `x²`) is solved just as readily, and a cubic with
roots of unity returns its complex roots in closed form:

```mathematica
In[1]:= Solve[x^4 - 5 x^2 + 4 == 0, x]
Out[1]= {{x -> -2}, {x -> -1}, {x -> 1}, {x -> 2}}

In[2]:= Solve[x^3 == 1, x]
Out[2]= {{x -> 1}, {x -> -(-1)^(1/3)}, {x -> (-1)^(2/3)}}
```

`In[1]` finds all four real roots of `x⁴ − 5x² + 4 = (x² − 1)(x² − 4)`. `In[2]`
gives the three cube roots of unity: the real root `1` and the two complex
roots written in terms of `(-1)^(1/3)`.

### Root objects

Not every polynomial can be solved in radicals. The general quintic and
beyond have no formula in terms of `+`, `−`, `×`, `÷`, and roots (this is the
Abel–Ruffini theorem). When `Solve` cannot — or by default will not — produce
a radical formula, it returns the roots as exact, symbolic **`Root` objects**:

```mathematica
In[1]:= Solve[x^5 - x + 1 == 0, x]
Out[1]= {{x -> Root[1 - #1 + #1^5 &, 1]}, {x -> Root[1 - #1 + #1^5 &, 2]}, {x -> Root[1 - #1 + #1^5 &, 3]}, {x -> Root[1 - #1 + #1^5 &, 4]}, {x -> Root[1 - #1 + #1^5 &, 5]}}
```

`Root[f &, k]` denotes the `k`-th root of the polynomial `f[#1] == 0` in
Mathilda's canonical ordering. These are genuine exact numbers — they can be
numericised to any precision and they obey the algebra of their defining
polynomial — they simply have no finite radical expression. Even a cubic is
returned as `Root` objects by default, because the radical formulas are
notoriously unwieldy:

```mathematica
In[1]:= Solve[x^3 - 3 x + 1 == 0, x]
Out[1]= {{x -> Root[1 - 3 #1 + #1^3 &, 1]}, {x -> Root[1 - 3 #1 + #1^3 &, 2]}, {x -> Root[1 - 3 #1 + #1^3 &, 3]}}
```

### Forcing radical formulas: Cubics and Quartics

To request the explicit radical (Cardano / Ferrari) formulas for degree-3 and
degree-4 equations, set the `Cubics` or `Quartics` option to `True`. Be warned
that the cubic formula is verbose:

```mathematica
In[1]:= Solve[x^3 - 3 x + 1 == 0, x, Cubics -> True]
Out[1]= {{x -> -1/3 ((1/2 (27 + (27*I) Sqrt[3]))^(1/3) + 9/(1/2 (27 + (27*I) Sqrt[3]))^(1/3))}, {x -> -1/3 (-(-1)^(1/3) (1/2 (27 + (27*I) Sqrt[3]))^(1/3) + 9 (-1)^(2/3)/(1/2 (27 + (27*I) Sqrt[3]))^(1/3))}, {x -> -1/3 ((-1)^(2/3) (1/2 (27 + (27*I) Sqrt[3]))^(1/3) - 9 (-1)^(1/3)/(1/2 (27 + (27*I) Sqrt[3]))^(1/3))}}

In[2]:= Solve[x^4 - 2 == 0, x, Quartics -> True]
Out[2]= {{x -> -2^(1/4)}, {x -> 2^(1/4)}, {x -> -I 2^(1/4)}, {x -> I 2^(1/4)}}
```

`In[1]` is Cardano's formula spelled out for all three roots — correct, but a
good illustration of *why* `Root` objects are the sensible default. `In[2]`,
solving `x⁴ = 2`, is far tidier: the four fourth-roots `±2^(1/4)` and
`±i·2^(1/4)`.

### Converting Root objects to radicals

`ToRadicals[expr]` attempts the reverse direction: it rewrites any `Root`
objects in `expr` as radicals when a closed form exists. It always succeeds for
degree ≤ 4 and for *binomial* roots `a #^n + b` of any degree, and leaves
genuinely unsolvable roots untouched:

```mathematica
In[1]:= ToRadicals[Root[#^4 - 2 &, 1]]
Out[1]= -2^(1/4)

In[2]:= ToRadicals[Root[#^5 - # + 1 &, 1]]
Out[2]= Root[#1^5 - #1 + 1 &, 1]
```

`In[1]` turns the first root of `x⁴ − 2` into `−2^(1/4)`. `In[2]` is a true
quintic with no radical form, so `ToRadicals` returns the `Root` object
unchanged — a faithful answer rather than a false one.

## Simultaneous equations

Pass `Solve` a list of equations and a list of the variables to solve for, and
it solves the system. Linear systems are handled directly, with numeric or
symbolic coefficients:

```mathematica
In[1]:= Solve[{2 x + y == 5, x - y == 1}, {x, y}]
Out[1]= {{x -> 2, y -> 1}}

In[2]:= Solve[{a x + b y == e, c x + d y == f}, {x, y}]
Out[2]= {{x -> (-d e + b f)/(b c - a d), y -> (c e - a f)/(b c - a d)}}
```

`In[1]` solves the concrete `2 × 2` system: adding the equations eliminates `y`
and gives `3x = 6`, so `x = 2`, `y = 1`. `In[2]` returns the general Cramer's-rule
solution, with the determinant `bc − ad` in every denominator.

### Eliminating a variable

Sometimes you do not want to solve outright — you want to *eliminate* a variable
and see the relation that survives among the rest. `Eliminate[eqns, vars]` does
exactly that, removing `vars` and returning the consequence:

```mathematica
In[1]:= Eliminate[{x + y == 2, x - y == 4}, y]
Out[1]= x == 3

In[2]:= Eliminate[{x == a + b, y == a - b}, a]
Out[2]= 2 b + y == x
```

`In[1]` removes `y` from the pair of equations, leaving the single condition
`x == 3`. `In[2]` removes `a`, leaving a relation tying `x`, `y`, and `b`
together. `Eliminate` is the engine behind several of the advanced examples
below — it is how you turn a parametric description into an implicit equation.

## Transcendental equations

Equations involving `Sin`, `Cos`, `Exp`, `Log`, and the like usually have
*infinitely many* solutions. Mathilda represents an entire family with a single
`ConditionalExpression` carrying an integer parameter:

```mathematica
In[1]:= Solve[Sin[x] == 0, x]
Out[1]= {{x -> ConditionalExpression[Pi + 2 C[1] Pi, Element[C[1], Integers]]}, {x -> ConditionalExpression[2 C[1] Pi, Element[C[1], Integers]]}}

In[2]:= Solve[Cos[x] == 1/2, x]
Out[2]= {{x -> ConditionalExpression[2 C[1] Pi - 1/3 Pi, Element[C[1], Integers]]}, {x -> ConditionalExpression[2 C[1] Pi + 1/3 Pi, Element[C[1], Integers]]}}
```

`In[1]` captures *all* solutions of `sin x = 0` — namely `x = kπ`, split into
the even multiples `2C[1]π` and the odd multiples `π + 2C[1]π`, with the
parameter `C[1]` ranging over the integers. `In[2]` gives `x = ±π/3 + 2C[1]π`.
The `Element[C[1], Integers]` clause records that `C[1]` must be an integer.
Exponential equations work the same way through the complex logarithm:

```mathematica
In[1]:= Solve[Exp[x] == 2, x]
Out[1]= {{x -> ConditionalExpression[Log[2] + (2*I) C[1] Pi, Element[C[1], Integers]]}}
```

The real solution `x = log 2` is the `C[1] = 0` member of the family; the others
differ by integer multiples of `2πi`, the period of the complex exponential.

### Naming the parameter: GeneratedParameters

The fresh parameter is called `C` by default, giving `C[1]`, `C[2]`, …. The
`GeneratedParameters` option chooses a different head:

```mathematica
In[1]:= Solve[Sin[x] == 0, x, GeneratedParameters -> n]
Out[1]= {{x -> ConditionalExpression[Pi + 2 Pi n[1], Element[n[1], Integers]]}, {x -> ConditionalExpression[2 Pi n[1], Element[n[1], Integers]]}}
```

The solution is identical to before, but the integer parameter is now reported
as `n[1]` — handy when `C` is already in use, or simply for readability.

### Turning inversion off: InverseFunctions

Peeling `Exp`, `Log`, and the trig functions off an equation relies on the
**inverse-function specialist**, controlled by the `InverseFunctions` option
(enabled by default). Setting it to `False` disables that machinery, so an
equation that can *only* be solved by inverting a transcendental head is
returned unevaluated:

```mathematica
In[1]:= Solve[Exp[x] == 2, x, InverseFunctions -> False]
Out[1]= Solve[E^x == 2, x, InverseFunctions -> False]
```

With inversion turned off Mathilda has no other way to reach the solution, so it
returns the input unchanged — a deliberate "I will not guess" rather than a
wrong answer. Leave the option at its default to get the `Log[2] + 2πi C[1]`
family shown earlier.

## Identities: SolveAlways

`Solve` answers "for which `x` is this true?"; `SolveAlways[eqn, x]` answers a
different question — "for which values of the *other* symbols is this an
**identity**, true for *every* `x`?" It returns conditions on the parameters,
not on `x`:

```mathematica
In[1]:= SolveAlways[a x + b == c x + d, x]
Out[1]= {{d -> b, c -> a}}

In[2]:= SolveAlways[(x + 1)^2 == x^2 + 2 x + a, x]
Out[2]= {{a -> 1}}
```

`In[1]` says `ax + b = cx + d` holds for all `x` precisely when the coefficients
match — `c = a` and `d = b`. `In[2]` finds the value `a = 1` that turns
`(x + 1)² = x² + 2x + a` into a true identity. This is exactly the
matching-coefficients reasoning you use to fix undetermined constants, done
automatically.

## Summing over roots: RootSum

`RootSum[f &, body &]` denotes the sum of `body[α]` over all roots `α` of the
polynomial `f[α] == 0`, *without* computing the individual roots. It stays
symbolic in general, but Mathilda recognises one important closed form — the
Lagrange / partial-fraction collapse that appears when you differentiate the
logarithmic part of a rational integral:

```mathematica
In[1]:= RootSum[#^2 - 3 # + 2 &, 1/(x - #)/(2 # - 3) &]
Out[1]= 1/(2 - 3 x + x^2)

In[2]:= RootSum[#^3 + p # + q &, 1/(x - #)/(3 #^2 + p) &]
Out[2]= 1/(q + p x + x^3)
```

In each case the body has the shape `1 / ((x − α) f'(α))`, and the sum over the
roots `α` collapses to `1 / f(x)` — the Hermite identity
`Σ 1/((x − αᵢ) f'(αᵢ)) = 1/f(x)`. `In[1]` does this for the quadratic
`x² − 3x + 2` and `In[2]` for the general depressed cubic `x³ + px + q`, the sum
folding back into the denominator polynomial. This is why `RootSum` matters even
though it rarely "expands": it lets the integrator carry an exact answer over
the roots of a polynomial it cannot factor.

## Advanced applications

The real payoff comes from combining `Solve`, `D`, and `Eliminate` to answer
genuine problems end to end. Each of the following runs exactly as shown.

### Where a curve has a horizontal tangent

A curve `y = f(x)` has a horizontal tangent wherever its derivative vanishes.
Take `f(x) = x³ − 3x`: differentiate, set the derivative to zero, and solve.

```mathematica
In[1]:= D[x^3 - 3 x, x]
Out[1]= -3 + 3 x^2

In[2]:= Solve[D[x^3 - 3 x, x] == 0, x]
Out[2]= {{x -> -1}, {x -> 1}}
```

`In[1]` gives `f'(x) = 3x² − 3`, and `In[2]` finds the two points `x = ±1` where
it is zero — the locations of the local maximum and minimum. To classify them,
read off the slope just by substituting back, confirming both are genuine
turning points (slope exactly `0`):

```mathematica
In[1]:= 3 - 3 x^2 /. x -> -1
Out[1]= 0

In[2]:= 3 - 3 x^2 /. x -> 1
Out[2]= 0
```

### An optimisation problem

The same machinery finds extrema. Minimise `g(x) = x³/3 − 4x` over the reals:
differentiate, solve `g'(x) = 0` for the critical points, then let Mathilda
classify them with the **second derivative test** — no graphing or guesswork.

```mathematica
In[1]:= D[x^3/3 - 4 x, x]
Out[1]= -4 + x^2

In[2]:= Solve[x^2 - 4 == 0, x]
Out[2]= {{x -> -2}, {x -> 2}}
```

The derivative `g'(x) = x² − 4` vanishes at `x = ±2`. To tell a maximum from a
minimum, take the *second* derivative with `D[g, {x, 2}]` and evaluate its sign
at each critical point: `g''(x) < 0` means the curve bends downward (a local
maximum), `g''(x) > 0` means it bends upward (a local minimum).

```mathematica
In[3]:= D[x^3/3 - 4 x, {x, 2}]
Out[3]= 2 x

In[4]:= D[x^3/3 - 4 x, {x, 2}] /. x -> -2
Out[4]= -4

In[5]:= D[x^3/3 - 4 x, {x, 2}] /. x -> 2
Out[5]= 4
```

`g''(x) = 2x` is negative at `x = −2` (`g'' = −4 < 0`, a local **maximum**) and
positive at `x = 2` (`g'' = 4 > 0`, a local **minimum**). Substituting back into
`g` confirms the values these extrema take:

```mathematica
In[6]:= (x^3/3 - 4 x) /. x -> -2
Out[6]= 16/3

In[7]:= (x^3/3 - 4 x) /. x -> 2
Out[7]= -16/3
```

So `g(−2) = 16/3` is the local maximum and `g(2) = −16/3` the local minimum.
Differentiate, solve `g' = 0`, then read the sign of `g''` at each root — the
entire calculus of optimisation, with the maximum/minimum verdict computed
rather than eyeballed.

### A geometric locus by elimination

`Eliminate` turns a *parametric* curve into its implicit equation. The curve
traced by `x = t²`, `y = t³` is the cuspidal cubic; eliminate the parameter `t`
to see the relation between `x` and `y`:

```mathematica
In[1]:= Eliminate[{x == t^2, y == t^3}, t]
Out[1]= y^2 == x^3
```

Eliminating `t` yields `y² = x³` — the famous semicubical (cuspidal) parabola,
recovered with no algebra on your part. The unit circle comes out the same way
from its trigonometric parametrisation:

```mathematica
In[1]:= Eliminate[{x == Cos[t], y == Sin[t]}, t]
Out[1]= x^2 + y^2 == 1
```

`Eliminate` discharges the parameter using `cos²t + sin²t = 1`, returning the
implicit equation `x² + y² = 1`.

### Intersecting a parabola and a line

To find where two curves meet, eliminate one coordinate to get a single
equation, then `Solve` it. Where does the parabola `y = x²` cross the line
`y = 2x + 3`?

```mathematica
In[1]:= Eliminate[{y == x^2, y == 2 x + 3}, y]
Out[1]= x^2 == 3 + 2 x

In[2]:= Solve[x^2 - 2 x - 3 == 0, x]
Out[2]= {{x -> -1}, {x -> 3}}
```

`In[1]` eliminates `y` to leave `x² = 2x + 3`; `In[2]` solves the resulting
quadratic for `x = −1` and `x = 3`. Substituting each back into `y = 2x + 3`
gives the two intersection points `(−1, 1)` and `(3, 9)` — the classic
eliminate-then-solve recipe for intersecting any two curves.

### Solving a polynomial system with `Resultant`

When both equations are nonlinear, handing the *system* straight to `Solve` can
leave it unevaluated — there is no single variable to isolate by inspection:

```mathematica
In[1]:= Solve[{x^2 + y^2 == 25, y == x^2 - 13}, {x, y}]
Out[1]= Solve[{x^2 + y^2 == 25, y == -13 + x^2}, {x, y}]
```

`Resultant[p, q, x]` is the classical elimination tool that breaks the deadlock.
It returns a polynomial in the *remaining* variable that vanishes exactly when
`p` and `q` share a root in `x` — eliminating `x` from the pair algebraically,
where eyeballing fails. Here the circle `x² + y² = 25` meets the parabola
`y = x² − 13`:

```mathematica
In[1]:= Resultant[x^2 + y^2 - 25, x^2 - y - 13, x]
Out[1]= 144 - 24 y - 23 y^2 + 2 y^3 + y^4

In[2]:= Solve[144 - 24 y - 23 y^2 + 2 y^3 + y^4 == 0, y]
Out[2]= {{y -> -4}, {y -> -4}, {y -> 3}, {y -> 3}}
```

Eliminating `x` collapses the two-variable system to one quartic in `y`, which
factors as `(y² + y − 12)² = ((y + 4)(y − 3))²`. The roots repeat because each
`y` is reached by a symmetric pair `x = ±√(y + 13)`: at `y = 3`, `x = ±4`; at
`y = −4`, `x = ±3` — the four intersection points `(±4, 3)` and `(±3, −4)`.
`Resultant` is the engine underneath `Eliminate`; reaching for it directly pays
off when you want the eliminated polynomial itself, multiplicities and all.

### Cracking a high-degree equation with `Decompose`

A quartic like `x⁴ + x² + 1` looks daunting head-on, but it is secretly a
*quadratic of a quadratic*. `Decompose` exposes that hidden layering, writing a
polynomial as a chain of lower-degree maps:

```mathematica
In[1]:= Decompose[x^4 + x^2 + 1, x]
Out[1]= {1 + x + x^2, x^2}
```

The result reads outermost-first: `x⁴ + x² + 1 = p(q(x))` with the outer
`p(u) = u² + u + 1` and the inner `q(x) = x²`. That immediately halves the work
— solve the outer quadratic for `u`, then the inner one for `x`:

```mathematica
In[2]:= Solve[u^2 + u + 1 == 0, u]
Out[2]= {{u -> 1/2 (-1 - I Sqrt[3])}, {u -> 1/2 (-1 + I Sqrt[3])}}

In[3]:= Solve[x^4 + x^2 + 1 == 0, x]
Out[3]= {{x -> -Sqrt[1/2 (-1 - I Sqrt[3])]}, {x -> Sqrt[1/2 (-1 - I Sqrt[3])]}, {x -> -Sqrt[1/2 (-1 + I Sqrt[3])]}, {x -> Sqrt[1/2 (-1 + I Sqrt[3])]}}
```

Each root of the outer quadratic feeds `x = ±√u`, and `In[3]` confirms it: the
four roots of the quartic are exactly the square roots of the two `u` values.
`Decompose` turns "solve a degree-4 equation" into "solve two quadratics" — and
it nests deeper, peeling a sextic such as `x⁶ + 3x⁴ + 3x² + 2` into
`{2 + 3x + 3x² + x³, x²}`, a cubic composed with a square.

## Where to next

You can now solve linear and polynomial equations, read solutions back as
substitution rules, handle irrational and complex roots, control radical
formulas with `Cubics`/`Quartics`/`ToRadicals`, work with `Root` objects for the
unsolvable cases, solve simultaneous systems and eliminate variables, capture
the infinite solution families of transcendental equations, verify identities
with `SolveAlways`, and combine all of this with `D` to solve real geometric and
optimisation problems.

- **[8. Calculus](08-calculus.md)** — differentiate and integrate, expand
  functions as power series, take limits and sums, and lean on `Solve` to locate
  the maxima, minima, and stationary points of the functions you build.
- **[Function reference](../documentation/index.md)** — the full catalogue of
  built-in functions. The
  [solutions of equations](../documentation/solutions-of-equations/index.md)
  category documents every function used above —
  `Solve`, `SolveAlways`, `Root`, `ToRadicals`, `Eliminate`, `RootSum`,
  and the `Cubics`, `Quartics`, `GeneratedParameters`, `InverseFunctions`, and
  `VerifySolutions` options — in full.
