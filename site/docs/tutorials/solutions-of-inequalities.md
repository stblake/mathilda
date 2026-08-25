# Solutions of Inequalities

An equation pins down a set of *points*; an inequality carves out a *region*.
This tutorial is a hands-on tour of how Mathilda solves inequalities — from a
single quadratic in one variable through rational functions, absolute values and
piecewise functions, and on to two- and three-variable regions of the plane and
space. The one tool at the centre of it is [`Reduce`](07-solutions-of-equations.md),
introduced in the [solutions of equations](07-solutions-of-equations.md) tutorial;
here it does the thing `Solve` cannot, because the answer is almost never a
finite list of values.

Every transcript below was produced by the actual Mathilda binary. Type the
`In[...]` lines yourself (without the prompt) and you will see the same
`Out[...]`. As with equations, inequalities use the ordinary comparison
operators `<`, `<=`, `>`, `>=` (and `!=` for "not equal"); a chained form such as
`-5 < 3 x + 7 <= 22` means both parts at once. The answer comes back as a logical
combination — `&&` ("and"), `||` ("or") — of simpler inequalities describing the
whole solution set.

!!! note "A word on domains"
    Ordering (`<`, `>`, …) is undefined over the complex numbers, so a statement
    that contains an inequality and is given no explicit domain is solved over
    the **reals** — matching Mathematica. That means the third argument `Reals`
    is optional whenever the statement already contains an inequality; the
    examples below often include it for clarity, but `Reduce[x^2 < 4, x]` and
    `Reduce[x^2 < 4, x, Reals]` mean exactly the same thing.

## Sign diagrams: one polynomial inequality

To solve a polynomial inequality, `Reduce` finds the real roots of the
polynomial — the points where it can change sign — and reports the union of the
intervals on which the inequality holds. This is the classic **sign diagram**
you draw by hand, done exactly:

```mathematica
In[1]:= Reduce[x^2 > 1, x, Reals]
Out[1]= x < -1 || x > 1

In[2]:= Reduce[x^2 - 5 x + 6 < 0, x, Reals]
Out[2]= 2 < x < 3
```

`In[1]` splits the line at the roots `±1` and keeps the two outer intervals where
`x² − 1 > 0`. `In[2]` factors as `(x − 2)(x − 3) < 0`, which holds only *between*
the roots. Higher degrees are no harder — the diagram just has more breakpoints:

```mathematica
In[1]:= Reduce[(x - 1) (x - 2) (x - 3) > 0, x, Reals]
Out[1]= 1 < x < 2 || x > 3

In[2]:= Reduce[x^4 - 5 x^2 + 4 < 0, x, Reals]
Out[2]= -2 < x < -1 || 1 < x < 2
```

The cubic alternates sign across its three roots, so the positive set is every
other interval; the quartic `(x² − 1)(x² − 4)` is negative on the two bands
between consecutive roots. When the breakpoints are irrational, `Reduce` keeps
them **exact** rather than rounding — it orders and signs algebraic numbers with
an exact real-algebraic oracle:

```mathematica
In[1]:= Reduce[x^2 < 2, x, Reals]
Out[1]= -Sqrt[2] < x < Sqrt[2]
```

A `<=` or `>=` closes the corresponding endpoint, and `!=` punches out isolated
points:

```mathematica
In[1]:= Reduce[x^2 <= 4, x, Reals]
Out[1]= -2 <= x <= 2

In[2]:= Reduce[x^2 != 1, x, Reals]
Out[2]= x != -1 && x != 1
```

## Combining conditions

Because the result is a logical formula, you can hand `Reduce` any combination of
inequalities joined with `&&` and `||`, and it intersects and unions the
solution sets accordingly. A chained inequality is just two bounds at once:

```mathematica
In[1]:= Reduce[-5 < 3 x + 7 <= 22, x]
Out[1]= -4 < x <= 5

In[2]:= Reduce[x^2 > 1 && x < 3, x, Reals]
Out[2]= x < -1 || 1 < x < 3

In[3]:= Reduce[1 <= x^2 <= 4, x, Reals]
Out[3]= -2 <= x <= -1 || 1 <= x <= 2
```

Note that `In[1]` was given no domain: because it contains an ordering, it is
solved over the reals automatically. `In[3]` — the set where `x²` lies between 1
and 4 — is a good illustration of why the answer is a formula and not a list: it
is the union of two symmetric closed bands.

## Rational functions and their poles

A rational inequality cannot be solved by "multiplying out" the denominator —
that would flip the inequality wherever the denominator is negative, and would
wrongly *include* the points where it vanishes. `Reduce` instead treats the roots
of the denominator as **poles**: extra breakpoints where the expression is
undefined and which are excluded from the solution set.

```mathematica
In[1]:= Reduce[1/x < 1, x, Reals]
Out[1]= x < 0 || x > 1

In[2]:= Reduce[1/x >= 0, x, Reals]
Out[2]= x > 0

In[3]:= Reduce[(x - 1)/(x - 2) > 0, x, Reals]
Out[3]= x < 1 || x > 2
```

`In[1]` is the tell-tale case: `1/x < 1` is satisfied by every negative `x` (where
`1/x` is itself negative) *and* by `x > 1`, but **not** by the interval `0 < x <= 1`.
Naïvely clearing the denominator would have given the wrong `x > 1` alone. In
`In[2]` the pole at `x = 0` is excluded, so the answer is the open ray `x > 0`,
not `x >= 0` — `1/x` is undefined at the origin.

## Deciding a statement: `True` and `False`

Some inequalities hold for *every* value of the variable, or for *none*. There
the "solution set" is the whole line or the empty set, and `Reduce` returns a
bare `True` or `False` — a genuine proof, not a sampled guess:

```mathematica
In[1]:= Reduce[x^2 + 1 > 0, x, Reals]
Out[1]= True

In[2]:= Reduce[x^2 >= 0, x, Reals]
Out[2]= True

In[3]:= Reduce[x^2 < 0, x, Reals]
Out[3]= False
```

`In[1]` proves that `x² + 1` is strictly positive for all real `x`; `In[3]` proves
that a real square is never negative, so the constraint is contradictory.

## Absolute values and piecewise functions

`Reduce` understands the real-valued constructs that make an expression
piecewise-defined: `Abs`, `Min`/`Max`, `Sign`, `Floor`/`Ceiling`/`Round`, `Mod`,
and the piecewise heads (`Piecewise`, `UnitStep`, `Ramp`, `Clip`,
`HeavisideTheta`, …). It eliminates each one by **case-splitting** — an `Abs` on
the sign of its argument, a `Max` on which argument is largest — solving each
branch on its own sign diagram and unioning the results.

```mathematica
In[1]:= Reduce[Abs[x] < 1, x, Reals]
Out[1]= -1 < x < 1

In[2]:= Reduce[Abs[x] >= 2, x, Reals]
Out[2]= x <= -2 || x >= 2

In[3]:= Reduce[Abs[2 x - 1] <= 3, x, Reals]
Out[3]= -1 <= x <= 2
```

The method is not limited to a single absolute value. A sum of them — the kind of
expression that defines a distance — is handled by splitting on every argument at
once:

```mathematica
In[1]:= Reduce[Abs[x - 1] + Abs[x + 1] <= 4, x, Reals]
Out[1]= -2 <= x <= 2

In[2]:= Reduce[Min[x, 1 - x] > 1/4, x, Reals]
Out[2]= 1/4 < x < 3/4

In[3]:= Reduce[Sign[x - 1] < 0, x, Reals]
Out[3]= x < 1
```

## Over the integers

Pass `Integers` as the domain and `Reduce` reports the integer solutions. When
the inequalities bound the variable to a finite range, it enumerates them:

```mathematica
In[1]:= Reduce[x^2 < 10 && x > 0, x, Integers]
Out[1]= x == 1 || x == 2 || x == 3
```

The bounded real set `0 < x < Sqrt[10]` contains exactly the integers 1, 2 and 3,
and `Reduce` lists them explicitly.

## Two variables: regions of the plane

With more than one variable the solution set is a region, and `Reduce` describes
it as a **triangular** cascade: a range for the first variable, then bounds on the
next that may depend on the first, and so on. A system of *linear* inequalities is
handled by Fourier–Motzkin elimination:

```mathematica
In[1]:= Reduce[x + y < 1 && x > 0 && y > 0, {x, y}, Reals]
Out[1]= 0 < x < 1 && 0 < y < 1 - x

In[2]:= Reduce[2 x + y <= 4 && x >= 0 && y >= 0, {x, y}, Reals]
Out[2]= 0 <= x <= 2 && 0 <= y <= 4 - 2 x
```

Read `In[1]` as a sweep: `x` ranges over `(0, 1)`, and for each such `x`, `y` runs
from `0` up to `1 − x` — exactly the open triangle with corners `(0,0)`, `(1,0)`,
`(0,1)`.

When the boundaries curve, the linear method no longer applies and `Reduce`
switches to **cylindrical algebraic decomposition** (CAD): it partitions the plane
into cells on which every defining polynomial keeps a constant sign, then keeps
the cells that satisfy the statement and reads off their bounds.

```mathematica
In[1]:= Reduce[x y > 0, {x, y}, Reals]
Out[1]= x < 0 && y < 0 || x > 0 && y > 0

In[2]:= Reduce[y > x^2, {x, y}, Reals]
Out[2]= y > x^2

In[3]:= Reduce[y >= x^2 && y <= x + 2, {x, y}, Reals]
Out[3]= -1 <= x <= 2 && x^2 <= y <= 2 + x
```

`In[1]` is the pair of open quadrants where `x` and `y` share a sign. `In[3]` is
the lens between the parabola `y = x²` and the line `y = x + 2`: the two curves
meet at `x = −1` and `x = 2`, and between them the region runs from the parabola
up to the line.

A curved boundary brings a radical into the answer. The open unit disk, for
instance, is `x` between `−1` and `1` with `y` between `±√(1 − x²)`:

```mathematica
In[1]:= Reduce[x^2 + y^2 < 1, {x, y}, Reals]
Out[1]= -1 < x < 1 && -1/2 Sqrt[4 - 4 x^2] < y < 1/2 Sqrt[4 - 4 x^2]
```

The bound `1/2 Sqrt[4 - 4 x^2]` is just `Sqrt[1 - x^2]` written in Mathilda's
canonical surface form — a good reminder that `Out[]` shows the system's exact
value, which is not always the tidiest way you would write it by hand. Closing
the inequality to `<=` closes the region's boundary, and the `x`-range along with
it:

```mathematica
In[1]:= Reduce[x^2 + y^2 <= 1, {x, y}, Reals]
Out[1]= -1 <= x <= 1 && -1/2 Sqrt[4 - 4 x^2] <= y <= 1/2 Sqrt[4 - 4 x^2]
```

The same machinery extends to three variables and beyond — an origin-centred
solid ball `x² + y² + z² < 1` comes back as a nested cascade of bounds, each
inner variable limited by the outer ones.

## Soundness: when `Reduce` declines

`Reduce` never returns a wrong or incomplete answer. When it cannot decide a sign
or ordering exactly — most often because a symbolic parameter's sign is unknown —
it leaves the input unevaluated rather than guess:

```mathematica
In[1]:= Reduce[a x + b > 0, x]
Out[1]= Reduce[b + a x > 0, x]
```

Solving `a x + b > 0` for `x` means dividing by `a` — which flips the inequality
when `a < 0` and collapses it when `a = 0`. Without knowing the sign of `a`,
there is no single correct direction, so `Reduce` declines. (Contrast this with
the *equation* `Reduce[a x == b, x]` from the previous tutorial, which *does*
return the full case tree, because an equation has no direction to flip.) A
returned formula, then, is always guaranteed to describe the whole solution set
exactly, and a returned `False` is a genuine proof that no solution exists.

## Where to next

You can now solve polynomial, rational, absolute-value and piecewise
inequalities in one variable, decide always-true and never-true statements,
enumerate integer solutions, and describe multi-variable regions of the plane and
space with `Reduce`.

- **[7. Solutions of equations](07-solutions-of-equations.md)** — the companion
  tutorial: `Solve` and `Reduce` on equations and systems, `Root` objects,
  `Eliminate`, and the transcendental solution families.
- **[8. Calculus](08-calculus.md)** — where inequalities describe the domain of a
  function, the sign of a derivative, or the region of convergence you are about
  to integrate over.
- **[`Reduce` reference page](../documentation/solutions-of-equations/Reduce.md)**
  — every domain and argument form, the engines behind them (sign diagrams,
  Fourier–Motzkin, CAD, the Diophantine reformatting), and their documented
  limits, in full.
