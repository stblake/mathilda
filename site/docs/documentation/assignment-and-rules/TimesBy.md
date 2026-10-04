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
In[3]:= x = 5
Out[3]= 5

In[4]:= x *= 4
Out[4]= 20

In[5]:= x
Out[5]= 20
```

## Implementation notes

**Algorithm.** `builtin_timesby` (`src/core.c`) implements `x *= dx` through the
shared `increment_core` helper with `mode = IC_MUL` and `pre = true`. The mode is
kept rather than a `negate` flag precisely because `TimesBy`/`DivideBy` differ
from `AddTo`/`SubtractFrom` in the combining **head**, not the sign: every other
step of the six increment operators is identical. `increment_core` resolves the
lvalue to its backing symbol (`lvalue_symbol_name` handles a plain symbol or a
`Part[sym, ...]` target) and, if that symbol has no OwnValue, emits
`TimesBy::rvalue` and returns `NULL` (the call stays unevaluated).

Otherwise it evaluates the current value, builds and evaluates `Times[old, dx]`,
then writes the result back via an evaluated `Set` (whose `HoldFirst` preserves a
`Part` lvalue so the mutation lands in place). The "pre" flag returns the **new**
value.

**Attributes & limits.** `TimesBy` is `HoldFirst | Protected` so the target is
not pre-evaluated; it requires exactly two arguments.

**Attributes:** `HoldFirst`, `Protected`.

## References

**See also:** [AddTo](../../assignment-and-rules/AddTo/), [SubtractFrom](../../assignment-and-rules/SubtractFrom/), [DivideBy](../../assignment-and-rules/DivideBy/), [HoldFirst](../../other-advanced/HoldFirst/), [Part](../../data-structures/Part/), [Plus](../../arithmetic/Plus/), [Times](../../arithmetic/Times/), [Increment](../../assignment-and-rules/Increment/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)

## Notes & additional examples

### Notes

`x *= dx` (`TimesBy`) is shorthand for `x = x dx`: it multiplies the current
value of `x` in place and returns the **new** value. It shares its whole
machinery with `AddTo`/`SubtractFrom`/`DivideBy`, differing only in that the
combining head is `Times`.

The target must already hold a value; applied to an unassigned symbol the call
is left unevaluated with a `TimesBy::rvalue` message. The left-hand side may be a
part target (`v[[i]] *= c`), and the product is formed and re-evaluated, so
threading over lists and symbolic simplification behave as for `x = x dx`.
