# Reduce

!!! warning "Status: Partial"
    implemented with documented limitations or caveats; some argument forms fall through to symbolic/unevaluated output.

## Description

**`Reduce[expr, vars]`**

Reduces the statement expr -- a logical combination (&&, ||, !, Implies, Xor) of equations (==, !=) and inequalities (\<, \<=, \>, \>=) -- to a complete, quantifier-free description of its solution set for the variables vars.  The default domain is Complexes, or Reals when expr contains an ordering inequality (ordering is undefined over the complexes), so e.g. Reduce\[-5 \< 3x+7 \<= 22, x\] is solved over the Reals.

**`Reduce[expr, vars, dom]`**

Reduces over the domain dom: Complexes, Reals, Integers, or Rationals.

<details>
<summary>Notes</summary>

Where Solve returns the generic solution of a set of equations as a list of replacement rules and silently drops the degenerate cases, Reduce returns an And/Or tree of equations and inequalities that describes the WHOLE solution set -- every parametric and boundary case kept -- and it solves inequalities over the reals.  A statement that decides returns True or False; an out-of-reach input is left unevaluated, never guessed. Handled so far: - Complexes: univariate polynomial equations carrying the full leading-coefficient case tree -- Reduce\[a x == b, x\] -\> (a != 0 && x == b/a) || (a == 0 && b == 0) -- and parametric linear systems by symbolic Gaussian elimination. - Reals, one variable: any Boolean combination of polynomial and rational-function equations and inequalities, solved as a union of intervals and points on an exact real-algebraic sign diagram (denominator roots are breakpoints and are excluded as poles); plus statements built from Abs, real radicals, Log, bounded inverse-trig, Floor/Ceiling/Round/Mod/IntegerPart, Min/Max, and the piecewise heads (Piecewise, Sign, UnitStep, Ramp, Clip, HeavisideTheta, Boole, ...). - Reals, several variables: linear systems by Fourier-Motzkin elimination and nonlinear systems (conics and beyond) by Cylindrical Algebraic Decomposition (McCallum projection), including multivariate Abs/Min/Max/piecewise selectors and square-root radical rationalization. - Integers / Rationals: the Solve Diophantine engine, reformatted as an Or of Ands with Element\[C\[k\], dom\] for a free parameter. Reduce is sound over complete: an undecidable sign, an unsupported construct, or a not-yet-wired engine (nonlinear equations over Complexes, and quantifier elimination) leaves the input unevaluated rather than risk a wrong formula.

</details>

## Examples (18)

Every input below was run against the current Mathilda build and its output recorded.

### Worked examples (1)

```mathematica
In[1]:= Reduce[x^2 < 1 && y^2 < 1 && z^2 < 1, {x, y, z}, Reals]
Out[1]= -1 < x < 1 && -1 < y < 1 && -1 < z < 1
```

### Applications (17)

```mathematica
In[2]:= Reduce[a x == b, x]
Out[2]= a != 0 && x == b/a || a == 0 && b == 0

In[3]:= Solve[a x == b, x]
Out[3]= {{x -> b/a}}

In[4]:= Reduce[x^2 == 4, x]
Out[4]= x == -2 || x == 2

In[5]:= Reduce[x^2 == -1, x]
Out[5]= x == -I || x == I

In[6]:= Reduce[a x^2 + b x + c == 0, x]
Out[6]= a != 0 && x == (1/2 (-b + Sqrt[b^2 - 4 a c]))/a || a != 0 && x == (1/2 (-b - Sqrt[b^2 - 4 a c]))/a || a == 0 && b != 0 && x == -c/b || a == 0 && b == 0 && c == 0

In[7]:= Reduce[x^2 > 1, x, Reals]
Out[7]= x < -1 || x > 1

In[8]:= Reduce[x^2 < 2, x, Reals]
Out[8]= -Sqrt[2] < x < Sqrt[2]

In[9]:= Reduce[(x - 1) (x - 2) (x - 3) > 0, x, Reals]
Out[9]= 1 < x < 2 || x > 3

In[10]:= Reduce[1/x < 1, x, Reals]
Out[10]= x < 0 || x > 1

In[11]:= Reduce[-5 < 3 x + 7 <= 22, x]
Out[11]= -4 < x <= 5

In[12]:= Reduce[Abs[x] < 1, x, Reals]
Out[12]= -1 < x < 1

In[13]:= Reduce[x^2 + 1 > 0, x, Reals]
Out[13]= True

In[14]:= Reduce[x < 0 && x > 1, x, Reals]
Out[14]= False

In[15]:= Reduce[x + y < 1 && x > 0 && y > 0, {x, y}, Reals]
Out[15]= 0 < x < 1 && 0 < y < 1 - x

In[16]:= Reduce[x y > 0, {x, y}, Reals]
Out[16]= x < 0 && y < 0 || x > 0 && y > 0

In[17]:= Reduce[2 x + 3 y == 1, {x, y}, Integers]
Out[17]= Element[C[1], Integers] && x == -1 + 3 C[1] && y == 1 - 2 C[1]

In[18]:= Reduce[x^2 < 10 && x > 0, x, Integers]
Out[18]= x == 1 || x == 2 || x == 3
```

## Algorithm

reduce.c

```text
`Reduce` -- front-end and dispatch skeleton (REDUCE_PLAN.md, Phase 0).
```

Mirrors builtin_solve's argument handling: positional `expr [, vars [, dom]]`, a Solve::ivar-style bad-variable diagnostic, and a True/False short-circuit. The input logical combination is normalised into the internal DNF layer (reduce_form.h); Phase 0 returns True/False when the statement decides and

```text
leaves everything else unevaluated (NULL).  The per-domain solving engines
```

(equational, sign-diagram, Fourier-Motzkin, CAD, integer, quantifier) are wired on top of this skeleton in later phases.

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [Unequal](../../comparisons/Unequal/), [Inequality](../../comparisons/Inequality/), [Solve](../../solutions-of-equations/Solve/), [Abs](../../arithmetic/Abs/), [Log](../../elementary-functions/Log/), [Floor](../../arithmetic/Floor/), [Ceiling](../../arithmetic/Ceiling/), [Round](../../arithmetic/Round/)

- Collins & Hong, "Partial Cylindrical Algebraic Decomposition for Quantifier Elimination", J. Symbolic Computation 12 (1991).
- McCallum, "An Improved Projection Operation for Cylindrical Algebraic Decomposition", in Quantifier Elimination and Cylindrical Algebraic Decomposition (1998).
- Basu, Pollack & Roy, "Algorithms in Real Algebraic Geometry" (2nd ed., 2006), Ch. 5 & 11 (sign determination, CAD).
- Schrijver, "Theory of Linear and Integer Programming" (1986), §12.2 (Fourier-Motzkin elimination).
- Source: [`src/solve/reduce.c`](https://github.com/stblake/mathilda/blob/main/src/solve/reduce.c)
- Specification: [`docs/spec/builtins/solutions-of-equations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/solutions-of-equations.md)
- Tests: [`tests/test_reduce.c`](https://github.com/stblake/mathilda/blob/main/tests/test_reduce.c)
- Tests: [`tests/test_reduce_corpus.c`](https://github.com/stblake/mathilda/blob/main/tests/test_reduce_corpus.c)

## Notes & additional examples

### Notes

`Reduce` and `Solve` answer different questions. `Solve[eqns, vars]` returns the
*generic* solution as a list of replacement rules and quietly discards the
special cases — `Solve[a x == b, x]` gives `{{x -> b/a}}`, silently assuming
`a != 0`. `Reduce[expr, vars]` returns instead a quantifier-free `And`/`Or` tree
of equations and inequalities that describes the **entire** solution set, every
parametric and boundary case included, so it reports both `a != 0 && x == b/a`
and the degenerate branch `a == 0 && b == 0`. Use `Solve` when you want to
substitute a solution back with `/.`; use `Reduce` when you need the complete
condition under which a statement holds — especially with symbolic parameters,
or with inequalities, where there is no finite list of rules to return.

**Domain.** The default domain is `Complexes`. Because an ordering (`<`, `<=`,
`>`, `>=`) is undefined over the complex numbers, a statement that contains one
and is given no explicit domain is solved over the `Reals` — matching
Mathematica — so the `Reals` argument is optional whenever the statement
already contains an inequality. Equations (`==`) and `Unequal` (`!=`) keep the
`Complexes` default. Pass a third argument (`Complexes`, `Reals`, `Integers`, or
`Rationals`) to choose explicitly.

**Soundness over completeness.** `Reduce` never returns a wrong or incomplete
formula: when a sign or ordering cannot be decided exactly, when it meets a
construct it does not yet handle, or when it reaches an engine not yet wired
(nonlinear equations over the complexes, and general quantifier elimination),
it leaves the input unevaluated rather than guess. An empty logical result
(`False`) is therefore a genuine proof that no solution exists, and a returned
formula is guaranteed to describe the whole solution set exactly.
