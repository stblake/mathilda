# Activate

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Activate[expr] reactivates every Inactive[h] in expr (replacing Inactive[h] by h) and re-evaluates, so an inactive integral Inactive[Integrate][g,x] becomes Integrate[g,x] and evaluates.`**

## Examples (1)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (1)

```mathematica
In[1]:= Activate[Inactive[Integrate][2 y, y]]
Out[1]= y^2
```

## Implementation notes

- `Protected`.
- Reverses `Inactive`: `Activate[Inactive[Integrate][g, x]]` becomes `Integrate[g, x]` and evaluates.

**Attributes:** `Protected`.

## References

**See also:** [Inactive](../../expression-information/Inactive/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_inactive.c`](https://github.com/stblake/mathilda/blob/main/tests/test_inactive.c)
