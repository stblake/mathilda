# TimesBy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`TimesBy[x, dx] or x *= dx`**

multiplies x by dx and returns the new value of x. x \*= dx is equivalent to x = x dx.

<details>
<summary>Notes</summary>

TimesBy has attribute HoldFirst. The first argument x can be a symbol or a Part expression referring to an existing value; dx may be a number, a symbolic expression, or a list (combined element-wise via the Listable attribute of Times). If x has no assigned value, TimesBy::rvalue is emitted and the expression is left unevaluated.

</details>

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= q = 5; q *= 3; q
Out[1]= 15

In[2]:= w = {1., 2., 3.}; w[[3]] /= 4.; w
Out[2]= {1.0, 2.0, 0.75}
```

## Implementation notes

**Attributes:** `HoldFirst`, `Protected`.

## See also

[AddTo](../../assignment-and-rules/AddTo/), [SubtractFrom](../../assignment-and-rules/SubtractFrom/), [DivideBy](../../assignment-and-rules/DivideBy/), [HoldFirst](../../other-advanced/HoldFirst/), [Part](../../structural-manipulation/Part/), [Plus](../../arithmetic/Plus/), [Times](../../arithmetic/Times/), [Increment](../../assignment-and-rules/Increment/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
