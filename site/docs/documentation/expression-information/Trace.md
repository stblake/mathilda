# Trace

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Trace[expr]
    Generates a list of the successive top-level expressions
    produced while evaluating expr. Returns {} when expr needs no
    rewriting. Sub-evaluations of arguments are not recorded (flat form).
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Trace[1 + 1]
Out[1]= {1 + 1, 2}

In[2]:= Trace[5]
Out[2]= {}

In[3]:= x = 3; Trace[x + 1]
Out[3]= {x + 1, 4}

In[4]:= Trace[{Trace[1 + 1], 2 + 2}]
Out[4]= {{Trace[1 + 1], 2 + 2}, {{1 + 1, 2}, 4}}
```

## Implementation notes

- `HoldAll`, `Protected`. The argument is held so its rewrite sequence can be

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/trace.c`](https://github.com/stblake/mathilda/blob/main/src/trace.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
