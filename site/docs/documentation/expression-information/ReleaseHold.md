# ReleaseHold

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
ReleaseHold[expr]
    removes Hold, HoldForm, HoldPattern, and HoldComplete in expr.
ReleaseHold removes only one layer of Hold etc.; it does not remove inner occurrences in nested Hold etc. functions.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Hold[1+1]
Out[1]= Hold[1 + 1]

In[2]:= ReleaseHold[Hold[1+1]]
Out[2]= 2

In[3]:= ReleaseHold /@ {Hold[1+2], HoldForm[2+3], HoldComplete[3+4]}
Out[3]= {3, 5, 7}

In[4]:= ReleaseHold[f[Hold[1+2]]]
Out[4]= f[3]

In[5]:= ReleaseHold[f[Hold[1+g[Hold[2+3]]]]]
Out[5]= f[1 + g[Hold[2 + 3]]]

In[6]:= ReleaseHold[Hold[Hold[1+1]]]
Out[6]= Hold[1 + 1]

In[7]:= ReleaseHold[Hold[Sin[Pi/6]]]
Out[7]= 1/2

In[8]:= ReleaseHold[{f[Hold[1+2]], g[HoldForm[3+4]]}]
Out[8]= {f[3], g[7]}
```

## Implementation notes

- `Protected`.
- `ReleaseHold` removes only one layer of `Hold` etc.; it does not remove inner occurrences in nested `Hold` etc. functions.
- `ReleaseHold` traverses into subexpressions of `expr` and strips any hold wrapper it finds, but does not recurse into the contents of the stripped wrapper.
- When `expr` does not contain any hold wrappers, `ReleaseHold` acts as identity.
- `ReleaseHold` removes all standard unevaluated containers: `Hold`, `HoldForm`, `HoldComplete`, and `HoldPattern`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
