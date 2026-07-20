# Sequence

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Sequence[e1, e2, ...]
    represents a sequence of arguments that is automatically spliced into
    the argument list of any enclosing function. Sequence[] evaporates and
    Sequence[e] acts like the identity. Splicing is suppressed for heads
    with the attribute SequenceHold or HoldAllComplete.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= f[a, Sequence[b, c], d]
Out[1]= f[a, b, c, d]

In[2]:= {a, Sequence[b], c, Identity[d]}
Out[2]= {a, b, c, d}

In[3]:= {a, b, g[x, y], h[w], g[z, y]} /. g -> Sequence
Out[3]= {a, b, x, y, h[w], z, y}
```

## Implementation notes

- Attributes: `{Protected}`.
- Splicing happens structurally during evaluation, before `Flat`/`Listable`/

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
