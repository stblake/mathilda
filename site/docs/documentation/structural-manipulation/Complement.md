# Complement

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Complement[eall, e1, e2, ...]
    gives the sorted list of distinct elements in eall that are not in
    any of the ei (set difference). The arguments must share a head,
    which need not be List.
Complement[eall, e1, ..., SameTest -> f]
    uses f[a, b] to decide whether elements a and b are the same.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Complement[{a, b, c, d, e}, {a, c}, {d}]
Out[1]= {b, e}

In[2]:= Complement[f[a, b, c, d], f[c, a], f[b, b, a]]
Out[2]= f[d]

In[3]:= Complement[{b, e, d, a, b, c, d}, {b, c}]
Out[3]= {a, d, e}
```

## Implementation notes

- `Protected`. Unlike `Union`/`Intersection`, `Complement` is order-sensitive in

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/structural-manipulation.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/structural-manipulation.md)
