# DivideBy

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DivideBy[x, dx] or x /= dx`**

divides x by dx and returns the new value of x. x /= dx is equivalent to x = x/dx.

<details>
<summary>Notes</summary>

DivideBy has attribute HoldFirst. The first argument x can be a symbol or a Part expression referring to an existing value. If x has no assigned value, DivideBy::rvalue is emitted and the expression is left unevaluated.

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
In[3]:= x = 12
Out[3]= 12

In[4]:= x /= 3
Out[4]= 4

In[5]:= x
Out[5]= 4
```

## Implementation notes

**Algorithm.** `builtin_divideby` (`src/core.c`) implements `x /= dx` through the
shared `increment_core` helper with `mode = IC_DIV` and `pre = true` — the same
worker behind `AddTo`/`SubtractFrom`/`TimesBy`, which differ only in how the
current value and `dx` are combined. `increment_core` first resolves the lvalue
to the symbol that actually holds the value (`lvalue_symbol_name` accepts a plain
symbol or a `Part[sym, ...]` target); if that symbol has no existing OwnValue it
emits `DivideBy::rvalue` and returns `NULL`, leaving the expression unevaluated.

Otherwise it evaluates the lvalue to the current value, builds and evaluates
`Times[old, dx^-1]` (so list threading and symbolic simplification happen exactly
as for the written-out form), and writes the result back through an evaluated
`Set`. `Set`'s `HoldFirst` preserves a compound lvalue shape such as
`Part[list, i]` for the assignment to update in place. The "pre" flag means the
**new** value is returned.

**Attributes & limits.** `DivideBy` is `HoldFirst | Protected` so the target is
not pre-evaluated; it requires exactly two arguments.

**Attributes:** `HoldFirst`, `Protected`.

## References

**See also:** [AddTo](../../assignment-and-rules/AddTo/), [SubtractFrom](../../assignment-and-rules/SubtractFrom/), [TimesBy](../../assignment-and-rules/TimesBy/), [HoldFirst](../../other-advanced/HoldFirst/), [Part](../../data-structures/Part/), [Plus](../../arithmetic/Plus/), [Times](../../arithmetic/Times/), [Increment](../../assignment-and-rules/Increment/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)

## Notes & additional examples

### Notes

`x /= dx` (`DivideBy`) is shorthand for `x = x/dx`: it divides the current value
of `x` in place and returns the **new** value. The target must already be a
variable with a value — applying it to an unassigned symbol leaves the
expression unevaluated and issues `DivideBy::rvalue`.

Because the update runs through `Set`, the left-hand side may be a compound
target such as `v[[i]] /= c`, and because the new value is formed by evaluating
`Times[old, dx^-1]`, list threading and symbolic simplification apply exactly as
they would to the written-out quotient.
