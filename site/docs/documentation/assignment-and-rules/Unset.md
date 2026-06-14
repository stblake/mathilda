# Unset

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Unset[lhs] or lhs =.
    removes any rule whose left-hand side is lhs, up to renaming of
    pattern variables. A bare symbol clears its value; a function form
    clears the matching definition on the head symbol.
Unset has attribute HoldFirst; Protected symbols are not affected.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= x = 5; x =.; x
Out[1]= x

In[2]:= f[x_] := x^2; f[x_] =.; f[3]
Out[2]= f[3]

In[3]:= fact[1] = 1; fact[n_] := n fact[n - 1]; fact[1] =.; fact[1]
Out[3]= 0
```

## Implementation notes

- `=.` is a low-precedence postfix operator (precedence 40, like `Set`), so it

**Attributes:** `HoldFirst`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)
