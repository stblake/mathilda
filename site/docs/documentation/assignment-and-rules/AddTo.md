# AddTo

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`AddTo[x, dx] or x += dx`**

adds dx to x and returns the new value of x. x += dx is equivalent to x = x + dx.

<details>
<summary>Notes</summary>

AddTo has attribute HoldFirst. The first argument x can be a symbol or a Part expression referring to an existing value; dx may be a number, a symbolic expression, or a list (combined element-wise via the Listable attribute of Plus). If x has no assigned value, AddTo::rvalue is emitted and the expression is left unevaluated.

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
In[3]:= x = 10
Out[3]= 10

In[4]:= x += 3
Out[4]= 13

In[5]:= x
Out[5]= 13
```

## Implementation notes

`builtin_addto` (`src/core.c`) implements `x += dx` via the shared `increment_core` helper (negate=false, pre=true). `increment_core` requires the lvalue to be a symbol with an existing OwnValue (else `AddTo::rvalue`), evaluates the current value, builds and evaluates `Plus[old, dx]`, then writes the new value back through an evaluated `Set` call (Set's `HoldFirst` preserves complex lvalue shapes like `Part[list, i]`). The "pre" flag means it returns the new value. `AddTo` itself is `ATTR_HOLDFIRST` so the target is not pre-evaluated.

**Attributes:** `HoldFirst`, `Protected`.

## References

**See also:** [SubtractFrom](../../assignment-and-rules/SubtractFrom/), [TimesBy](../../assignment-and-rules/TimesBy/), [DivideBy](../../assignment-and-rules/DivideBy/), [HoldFirst](../../other-advanced/HoldFirst/), [Part](../../structural-manipulation/Part/), [Plus](../../arithmetic/Plus/), [Times](../../arithmetic/Times/), [Increment](../../assignment-and-rules/Increment/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)
- Tests: [`tests/test_increment.c`](https://github.com/stblake/mathilda/blob/main/tests/test_increment.c)

## Notes & additional examples

### Notes

`x += dx` (`AddTo`) is equivalent to `x = x + dx`: it updates `x` in place and returns the new value.
