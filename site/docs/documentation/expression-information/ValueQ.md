# ValueQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`ValueQ[expr]`**

gives True if a value has been defined for expr, False otherwise.

<details>
<summary>Notes</summary>

HoldAll: inspects the symbol itself, not its evaluated value. A bare symbol tests OwnValues; f\[...\] tests whether head f has any DownValues.

</details>

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= ValueQ[x]
Out[1]= False

In[2]:= x = 5; ValueQ[x]
Out[2]= True

In[3]:= f[x_] := x^2; {ValueQ[f[2]], ValueQ[f[a, b]]}   (* head has DownValues *)
Out[3]= {True, True}

In[4]:= ValueQ[f]                       (* bare symbol, only DownValues *)
Out[4]= False

In[5]:= ValueQ /@ Unevaluated[{x, y}]   (* HoldAll preserved via Unevaluated *)
Out[5]= {True, False}
```

## Implementation notes

**Attributes:** `HoldAll`, `Protected`.

## See also

[HoldAll](../../expression-information/HoldAll/)

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_core.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core.c)
