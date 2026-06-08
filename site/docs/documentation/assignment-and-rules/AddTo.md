# AddTo

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
AddTo[x, dx] or x += dx
    adds dx to x and returns the new value of x.
    x += dx is equivalent to x = x + dx.

AddTo has attribute HoldFirst. The first argument x can be a symbol or
a Part expression referring to an existing value; dx may be a number,
a symbolic expression, or a list (combined element-wise via the Listable
attribute of Plus). If x has no assigned value, AddTo::rvalue is emitted
and the expression is left unevaluated.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_addto` (`src/core.c`) implements `x += dx` via the shared `increment_core` helper (negate=false, pre=true). `increment_core` requires the lvalue to be a symbol with an existing OwnValue (else `AddTo::rvalue`), evaluates the current value, builds and evaluates `Plus[old, dx]`, then writes the new value back through an evaluated `Set` call (Set's `HoldFirst` preserves complex lvalue shapes like `Part[list, i]`). The "pre" flag means it returns the new value. `AddTo` itself is `ATTR_HOLDFIRST` so the target is not pre-evaluated.

**Attributes:** `HoldFirst`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= x = 10
Out[1]= 10

In[2]:= x += 3
Out[2]= 13

In[3]:= x
Out[3]= 13
```

### Notes

`x += dx` (`AddTo`) is equivalent to `x = x + dx`: it updates `x` in place and returns the new value.
