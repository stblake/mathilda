# Continue

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Continue[] proceeds to the next iteration of the nearest enclosing Do, For, or While loop.
Continue[] skips the remainder of the current loop body.
Continue[] takes effect as soon as it is evaluated.
Continue has attribute Protected.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= r = 0; Do[If[EvenQ[i], Continue[]]; r += i, {i, 10}]; r
Out[1]= 25

In[2]:= r = 0; For[i = 1, i <= 10, i++, If[EvenQ[i], Continue[]]; r += i]; r
Out[2]= 25
```

## Implementation notes

- Has attribute `Protected`.
- Takes effect as soon as it is evaluated. In `Do` it advances the iterator and

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/control-flow.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/control-flow.md)
