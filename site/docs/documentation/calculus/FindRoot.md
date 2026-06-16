# FindRoot

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FindRoot[f, {x, x0}]
    searches for a numerical root of f starting from x = x0.
FindRoot[lhs == rhs, {x, x0}]
    searches for a numerical solution to the equation.
FindRoot[f, {x, x0, x1}]
    uses a variant of the secant method with x0 and x1 as the first two approximations.
FindRoot[f, {x, xstart, xmin, xmax}]
    uses Brent's method on the bracket [xmin, xmax].
FindRoot[{f1, f2, ...}, {{x, x0}, {y, y0}, ...}]
    searches for a simultaneous numerical root of the system.

Options: Method ('Newton' | 'Secant' | 'Brent' | Automatic), WorkingPrecision, MaxIterations, AccuracyGoal, PrecisionGoal, DampingFactor, Jacobian, StepMonitor, EvaluationMonitor.  FindRoot has HoldAll and effectively uses Block to localize variables.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= FindRoot[Sin[x] + Exp[x], {x, 0}]
Out[1]= {x -> -0.588533}

In[2]:= FindRoot[Cos[x] == x, {x, 0}]
Out[2]= {x -> 0.739085}

In[3]:= FindRoot[{y == Exp[x], x + y == 2}, {{x, 1}, {y, 1}}]
Out[3]= {x -> 0.442854, y -> 1.55715}

In[4]:= FindRoot[Sin[x], {x, 3}, WorkingPrecision -> 50]
Out[4]= {x -> 3.14159265358979323846264338328757519744320987120168}

In[5]:= FindRoot[(Cos[z + I] - 2) (z + 2), {z, 1.0 + 0.1 I}]
Out[5]= {z -> -1.66935e-13 + 0.316958*I}

In[6]:= FindRoot[Cos[x] - x, {x, 0, 1}, Method -> "Brent"]
Out[6]= {x -> 0.739085}

In[7]:= FindRoot[x^2 - 2, {x, 1.0, 2.0}, Method -> "Secant"]
Out[7]= {x -> 1.41421}

In[8]:= FindRoot[(x - 1)^3, {x, 0.5}, DampingFactor -> 3]
Out[8]= {x -> 1.0}
```

## Implementation notes

**Algorithm.** `FindRoot` (`HoldAll | Protected`) does iterative numerical
root-finding (src/findroot.c). It Block-binds the search variables via temporary
OwnValues so the global symbol table is unperturbed, restoring them on every
exit path. The variable spec selects the method:

- `{x, x0}` → **Newton** from one start. The derivative is obtained
  symbolically by evaluating `D[f, x]` (`fr_compute_derivative`) and re-evaluated
  numerically each step; the update is `x -= damping·f/f'`. Real and complex
  variants exist (`fr_run_newton_real`, `fr_run_newton_complex`).
- `{x, x0, x1}` → **secant** (`fr_run_secant_real`), used as the derivative-free
  fallback.
- `{x, xstart, xmin, xmax}` → **Brent** bracketing (`fr_run_brent_real`).
- `{{x, x0}, {y, y0}, ...}` with a list of equations → **multivariate Newton**
  (`fr_run_newton_system_real`): the Jacobian `J[i][j] = D[f_i, var_j]` is
  precomputed symbolically (or taken from a user `Jacobian` option),
  re-evaluated each step, and the Newton step solves `J·dx = f` by **Gaussian
  elimination with partial pivoting** (inline LU-style forward elimination +
  back substitution), then `x -= damping·dx`.

Equation forms `lhs == rhs` are normalised to `lhs - rhs`. Options include
`Method` (Automatic/"Newton"/"Secant"/"Brent"), `WorkingPrecision`
(MachinePrecision or an MPFR bit count), `MaxIterations` (default 100),
`AccuracyGoal`/`PrecisionGoal`, `DampingFactor`, `Jacobian`, and held
`StepMonitor`/`EvaluationMonitor`. When `WorkingPrecision` requests extended
precision and MPFR is built in, dedicated MPFR Newton paths
(`fr_run_newton_mpfr_real`, …) run the iteration at the requested bit width.

**Data structures.** `double` / C99 `double complex` for machine precision;
`mpfr_t` scalars for extended precision; a flat row-major `double` augmented
matrix for the system solve. Function/derivative evaluation re-enters the
Mathilda evaluator with the current numeric binding installed.

**Complexity / limits.** Per-iteration cost is dominated by symbolic-expression
evaluation (and `O(n^3)` Gaussian elimination for an n-variable system).
Convergence is local — quadratic for Newton near a simple root. Emits
diagnostics (`dsing`, `noconv`, `nlnum`) and returns NULL/unevaluated when the
derivative vanishes, the Jacobian is singular, evaluation is non-numeric, or it
fails to converge.

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- W. H. Press et al., *Numerical Recipes*, 3rd ed. (Cambridge, 2007) — Newton, secant, Brent's method.
- Source: [`src/findroot.c`](https://github.com/stblake/mathilda/blob/main/src/findroot.c)
- Specification: [`docs/spec/builtins/calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/calculus.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= FindRoot[x^2 - 2, {x, 1}]
Out[1]= {x -> 1.41421}
```

```mathematica
In[1]:= FindRoot[Cos[x] == x, {x, 1}, WorkingPrecision -> 40]
Out[1]= {x -> 0.73908513321516064165531208767387340401341}
```

```mathematica
In[1]:= FindRoot[BesselJ[0, x], {x, 2, 3}]
Out[1]= {x -> 2.40483}
```

```mathematica
In[1]:= FindRoot[Sin[x] == 0, {x, 3}, WorkingPrecision -> 40]
Out[1]= {x -> 3.1415926535897932384626433832875751974431}
```

```mathematica
In[1]:= FindRoot[{x^2 + y^2 == 1, x == y}, {{x, 1}, {y, 1}}]
Out[1]= {x -> 0.707107, y -> 0.707107}
```

### Notes

`FindRoot` accepts a function (sought equal to zero) or an explicit
equation, and honours `WorkingPrecision` via MPFR. The second example pins
the Dottie number — the unique real fixed point of cosine — to 40 digits;
the third finds the first positive zero of the Bessel function `J0` by
bracketing on `[2, 3]`; the fourth recovers `π` as a root of `Sin` to full
40-digit precision. The last example solves a nonlinear 2x2 system with a
Newton step using the analytic Jacobian.
