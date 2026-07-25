# Break

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Break[] exits the nearest enclosing Do, For, or While loop.
After Break[], the enclosing loop returns Null.
Break[] takes effect as soon as it is evaluated.
Break has attribute Protected.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= For[i = 1, i <= 10, i++, If[i > 2, Break[]]]; i
Out[1]= 3
```

## Implementation notes

- Has attribute `Protected`.
- Takes effect as soon as it is evaluated (e.g. inside an `If` within the body),

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
