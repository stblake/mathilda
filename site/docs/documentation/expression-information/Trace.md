# Trace

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Trace[expr]
    Generates a nested list of the expressions produced while
    evaluating expr. Each argument sub-evaluation that takes a step appears
    as a sublist, mirroring the structure of the evaluation. Returns {} when
    expr needs no rewriting.
Trace[expr, form]
    Includes only the steps whose expression matches the
    pattern form (e.g. Trace[expr, _Integer]).
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Trace[1 + 1]
Out[1]= {1 + 1, 2}

In[2]:= Trace[5]
Out[2]= {}

In[3]:= Trace[2^3 + 4^2 + 1]
Out[3]= {{2^3, 8}, {4^2, 16}, 8 + 16 + 1, 25}

In[4]:= Trace[x^Range[10]]
Out[4]= {{Range[10], {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}}, x^{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, {x, x^2, x^3, x^4, x^5, x^6, x^7, x^8, x^9, x^10}}

In[5]:= g[y_] := y^2; Trace[g[1 + 1]]
Out[5]= {{1 + 1, 2}, g[2], 2^2, 4}

In[6]:= Trace[1 + 2 + 3, _Integer]
Out[6]= {6}

In[7]:= Trace[Nest[f, x, 3], _f]
Out[7]= {f[f[f[x]]]}
```

## Implementation notes

- `HoldAll`, `Protected`. The argument is held so its rewrite sequence can be

**Attributes:** `HoldAll`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/trace.c`](https://github.com/stblake/mathilda/blob/main/src/trace.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
