# NMaximize

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`NMaximize[f, x]`**

searches for a global maximum of f with respect to x.

**`NMaximize[f, {x, y, ...}]`**

global maximum with respect to several variables.

**`NMaximize[{f, cons}, vars]`**

global maximum of f subject to the constraints cons.

<details>
<summary>Notes</summary>

NMaximize (Protected, not HoldAll) shares NMinimize's methods, options, and constraint/domain handling.  Internally maximises by minimising -f, then negates the objective value in the result.  Returns {fmax, {x -\> xmax, ...}}.

</details>

## Examples (9)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (6)

```mathematica
In[1]:= NMinimize[x^4 - 3 x^2 - x, x]
Out[1]= {-3.51391, {x -> 1.30084}}

In[2]:= NMinimize[{x + y, x^2 + y^2 <= 9}, {x, y}]
Out[2]= {-4.24264, {x -> -2.12132, y -> -2.12132}}

In[3]:= NMinimize[{x + 2 y, x^2 + 2 y^2 <= 3, x + y == 2, x >= 1}, {x, y}]
Out[3]= {2.33333, {x -> 1.66667, y -> 0.333333}}

In[4]:= NMinimize[{x + y, x + 2 y >= 3, x >= -2}, {Element[x, Integers], Element[y, Integers]}]
Out[4]= {1.0, {x -> -2, y -> 3}}

In[5]:= NMinimize[{x, x > 2 && x < 1}, x]
Out[5]= {Infinity, {x -> Indeterminate}}

In[6]:= NMaximize[{x + y, x^2 + y^2 <= 1}, {x, y}]
Out[6]= {1.41421, {x -> 0.707107, y -> 0.707107}}
```

### Applications (3)

A downward parabola, maximum at the vertex

```mathematica
In[7]:= NMaximize[-x^2 + 4 x, x]
Out[7]= {4.0, {x -> 2.0}}
```

A concave bowl, maximum at its centre

```mathematica
In[8]:= NMaximize[4 - (x - 1)^2 - (y - 2)^2, {x, y}]
Out[8]= {4.0, {x -> 1.0, y -> 2.0}}
```

The arithmetic-geometric-mean extremum on x + y = 10

```mathematica
In[9]:= NMaximize[{x y, x + y == 10}, {x, y}]
Out[9]= {25.0, {x -> 5.0, y -> 5.0}}
```

## Implementation notes

**Algorithm.** `builtin_nmaximize` is a thin wrapper over the `NMinimize`
driver. It constructs a synthetic `NMinimize` call that **minimises `−f`**: for a
bare objective it negates `f`; for a `{f, c1, c2, …}` list it negates only the
objective element (`Times[-1, f]`) and copies every constraint unchanged — a
fix for an earlier bug where negating the whole list threaded `Times[-1, …]` over
the constraints and corrupted them. It then calls `nm_minimize_driver(...,
"NMaximize")` and **negates the reported optimum** back (`mpfr_neg` or `-fmin`,
preserving the numeric type). Everything else — the global-method catalogue
(differential evolution by default, with Deb's feasibility rules), the variable
and constraint grammar, integer domains, the BFGS/augmented-Lagrangian local
polish, auto-compilation at `MachinePrecision`, and fixed-seed determinism — is
`NMinimize`'s, unchanged. See the `NMinimize` implementation notes for the full
account.

**Data structures.** Identical to `NMinimize`: the synthetic objective is
compiled to bytecode and evaluated per trial point; variables are bound and
restored `Block`-style; populations and simplices are flat `double*` buffers; a
`WorkingPrecision > machine` request uses the MPFR BFGS refinement.

**Complexity / limits.** As for `NMinimize`. Returns `{fmax, {x -> xmax, …}}`;
an empty or unmet feasible set returns `{-Infinity, {x -> Indeterminate, …}}`
(the negation of `NMinimize`'s infeasible `Infinity`). `NMaximize` is `Protected`
but **not** `HoldAll`, and shares `NMinimize`'s options
(`Method`, `MaxIterations`, `WorkingPrecision`, `AccuracyGoal`, `PrecisionGoal`,
`EvaluationMonitor`, `StepMonitor`) and message routing.

**Attributes:** `Protected`.

## References

**See also:** [NMinimize](../../numerical-calculus/NMinimize/), [FindMinimum](../../numerical-calculus/FindMinimum/), [Block](../../scoping-constructs/Block/), [HoldAll](../../expression-information/HoldAll/), [Rule](../../assignment-and-rules/Rule/), [Sqrt](../../arithmetic/Sqrt/), [Round](../../arithmetic/Round/), [AccuracyGoal](../../other-advanced/AccuracyGoal/)

- R. Storn and K. Price, *Differential Evolution — a simple and efficient heuristic for global optimization over continuous spaces*, J. Global Optim. **11** (1997) 341–359.
- K. Deb, *An efficient constraint handling method for genetic algorithms*, Comput. Methods Appl. Mech. Engrg. **186** (2000) 311–338.
- Source: [`src/numerical_calculus/nm_driver.c`](https://github.com/stblake/mathilda/blob/main/src/numerical_calculus/nm_driver.c)
- Specification: [`docs/spec/builtins/numerical-calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/numerical-calculus.md)
- Tests: [`tests/test_basin_hopping.c`](https://github.com/stblake/mathilda/blob/main/tests/test_basin_hopping.c)
- Tests: [`tests/test_direct.c`](https://github.com/stblake/mathilda/blob/main/tests/test_direct.c)
- Tests: [`tests/test_dual_annealing.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dual_annealing.c)
- Tests: [`tests/test_nminimize.c`](https://github.com/stblake/mathilda/blob/main/tests/test_nminimize.c)

## Notes & additional examples

### Notes

`NMaximize[f, vars]` searches for a global maximum and returns
`{fmax, {x -> xmax, ...}}`. It is implemented by **minimising `-f` and negating
the objective value**, so it shares `NMinimize`'s methods, options, and
constraint/domain handling in full — see the `NMinimize` page for the method
catalogue (`DifferentialEvolution` by default), the variable and constraint
grammar, integer domains, auto-compilation at `MachinePrecision`, and the
fixed-seed determinism these results rely on. The constrained example attains
`x y = 25` at `x = y = 5`, the classic AM-GM extremum on `x + y = 10`.
`NMaximize` is `Protected` but not `HoldAll`.
