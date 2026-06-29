# PartitionsP

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
PartitionsP[n]
    gives the number p(n) of unrestricted partitions of the integer n.
n must be an integer; p(n) = 0 for n < 0. Threads over lists.
For the partitions themselves use IntegerPartitions[n].
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Table[PartitionsP[k], {k, 0, 12}]
Out[1]= {1, 1, 2, 3, 5, 7, 11, 15, 22, 30, 42, 56, 77}

In[2]:= PartitionsP[100]
Out[2]= 190569292

In[3]:= PartitionsP[4096]
Out[3]= 6927233917602120527467409170319882882996950147283323368445315320451

In[4]:= Table[Times @@ PartitionsP[Last /@ FactorInteger[n]], {n, 12}]
Out[4]= {1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 2}
```

## Implementation notes

- `Protected`, `Listable` — `PartitionsP[{2, 4, 6}]` → `{2, 5, 11}`.
- Two engines, dispatched by the size of `n` (threshold `n = 1000`):
  - **Small `n`** — Euler's pentagonal-number-theorem recurrence

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
