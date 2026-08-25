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

## Examples (6)

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
Out[4]= {1.0, {x -> -1, y -> 2}}

In[5]:= NMinimize[{x, x > 2 && x < 1}, x]
Out[5]= {Infinity, {x -> Indeterminate}}

In[6]:= NMaximize[{x + y, x^2 + y^2 <= 1}, {x, y}]
Out[6]= {1.41421, {x -> 0.707107, y -> 0.707107}}
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [NMinimize](../../numerical-calculus/NMinimize/), [FindMinimum](../../numerical-calculus/FindMinimum/), [Block](../../scoping-constructs/Block/), [HoldAll](../../expression-information/HoldAll/), [Rule](../../assignment-and-rules/Rule/), [Sqrt](../../arithmetic/Sqrt/), [Round](../../arithmetic/Round/), [AccuracyGoal](../../other-advanced/AccuracyGoal/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/numerical-calculus.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/numerical-calculus.md)
- Tests: [`tests/test_basin_hopping.c`](https://github.com/stblake/mathilda/blob/main/tests/test_basin_hopping.c)
- Tests: [`tests/test_direct.c`](https://github.com/stblake/mathilda/blob/main/tests/test_direct.c)
- Tests: [`tests/test_dual_annealing.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dual_annealing.c)
- Tests: [`tests/test_nminimize.c`](https://github.com/stblake/mathilda/blob/main/tests/test_nminimize.c)
