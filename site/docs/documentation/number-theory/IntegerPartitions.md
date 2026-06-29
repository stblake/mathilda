# IntegerPartitions

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
IntegerPartitions[n]
    gives the partitions of n in reverse-lexicographic order.
IntegerPartitions[n, k] gives partitions into at most k parts;
{k} exactly k; {kmin, kmax} between; {kmin, kmax, dk} stepped.
A third argument restricts the parts (sspec; All = Range[n]); a
fourth limits the result to the first m (m>0) or last |m| (m<0).
n and the parts may be rational and negative; Length equals
PartitionsP[n] for the plain form.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= IntegerPartitions[5]
Out[1]= {{5}, {4, 1}, {3, 2}, {3, 1, 1}, {2, 2, 1}, {2, 1, 1, 1}, {1, 1, 1, 1, 1}}

In[2]:= IntegerPartitions[50, All, {6, 9, 20}]
Out[2]= {{20, 9, 9, 6, 6}, {20, 6, 6, 6, 6, 6}}

In[3]:= IntegerPartitions[5, 10, {1, -1}]
Out[3]= {{-1, -1, 1, 1, 1, 1, 1, 1, 1}, {-1, 1, 1, 1, 1, 1, 1}, {1, 1, 1, 1, 1}}
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
