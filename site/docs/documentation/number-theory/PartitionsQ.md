# PartitionsQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
PartitionsQ[n]
    gives the number q(n) of partitions of the integer n into distinct
parts (equivalently, into odd parts). n must be an integer; q(n) = 0
for n < 0. Threads over lists. For the partitions themselves use
IntegerPartitions[n, All, Range[n]] with distinct parts.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Table[PartitionsQ[k], {k, 0, 20}]
Out[1]= {1, 1, 1, 2, 2, 3, 4, 5, 6, 8, 10, 12, 15, 18, 22, 27, 32, 38, 46, 54, 64}

In[2]:= PartitionsQ[100]
Out[2]= 444793

In[3]:= PartitionsQ[{2, 4, 6}]
Out[3]= {1, 2, 4}
```

## Implementation notes

- `Protected`, `Listable` — `PartitionsQ[{2, 4, 6}]` → `{1, 2, 4}`.
- Two engines, dispatched by the size of `n` (threshold `n = 1000`):
  - **Small `n`** — an exact GMP recurrence derived from the Euler identity

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
