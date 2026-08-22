# SubtractFrom

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`SubtractFrom[x, dx] or x -= dx`**

subtracts dx from x and returns the new value of x. x -= dx is equivalent to x = x - dx.

<details>
<summary>Notes</summary>

SubtractFrom has attribute HoldFirst. The first argument x can be a symbol or a Part expression referring to an existing value; dx may be a number, a symbolic expression, or a list. If x has no assigned value, SubtractFrom::rvalue is emitted and the expression is left unevaluated.

</details>

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= q = 5; q *= 3; q
Out[1]= 15

In[2]:= w = {1., 2., 3.}; w[[3]] /= 4.; w
Out[2]= {1.0, 2.0, 0.75}
```

### Applications (3)

```mathematica
In[3]:= y = 10
Out[3]= 10

In[4]:= y -= 4
Out[4]= 6

In[5]:= y
Out[5]= 6
```

## Implementation notes

`builtin_subtractfrom` (`x -= dx`) delegates to the shared `increment_core(lhs, dx, negate=true, pre=true, "SubtractFrom")`. That helper requires `lhs` to be an lvalue with an existing OwnValue (else it emits `SubtractFrom::rvalue` and returns `NULL`), reads the old value via `evaluate`, builds and evaluates `Plus[old, Times[-1, dx]]`, writes the result back through an evaluated `Set[lhs, new]` (preserving lvalue shape such as `Part[...]` via Set's `HoldFirst`), and returns the new value. Has `ATTR_HOLDFIRST` so the target symbol is not evaluated before mutation.

**Attributes:** `HoldFirst`, `Protected`.

## References

**See also:** [AddTo](../../assignment-and-rules/AddTo/), [TimesBy](../../assignment-and-rules/TimesBy/), [DivideBy](../../assignment-and-rules/DivideBy/), [HoldFirst](../../other-advanced/HoldFirst/), [Part](../../structural-manipulation/Part/), [Plus](../../arithmetic/Plus/), [Times](../../arithmetic/Times/), [Increment](../../assignment-and-rules/Increment/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)
- Tests: [`tests/test_increment.c`](https://github.com/stblake/mathilda/blob/main/tests/test_increment.c)

## Notes & additional examples

### Notes

`y -= dy` (`SubtractFrom`) is equivalent to `y = y - dy`: it updates `y` in place and returns the new value.
