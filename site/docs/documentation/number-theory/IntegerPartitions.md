# IntegerPartitions

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`IntegerPartitions[n]`**

gives the partitions of n in reverse-lexicographic order.

**`IntegerPartitions[n, k] gives partitions into at most k parts;`**

**`PartitionsP[n] for the plain form.`**

<details>
<summary>Notes</summary>

{k} exactly k; {kmin, kmax} between; {kmin, kmax, dk} stepped. A third argument restricts the parts (sspec; All = Range\[n\]); a fourth limits the result to the first m (m\>0) or last |m| (m\<0). n and the parts may be rational and negative; Length equals

</details>

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= IntegerPartitions[5]
Out[1]= {{5}, {4, 1}, {3, 2}, {3, 1, 1}, {2, 2, 1}, {2, 1, 1, 1}, {1, 1, 1, 1, 1}}

In[2]:= IntegerPartitions[50, All, {6, 9, 20}]
Out[2]= {{20, 9, 9, 6, 6}, {20, 6, 6, 6, 6, 6}}

In[3]:= IntegerPartitions[5, 10, {1, -1}]
Out[3]= {{-1, -1, 1, 1, 1, 1, 1, 1, 1}, {-1, 1, 1, 1, 1, 1, 1}, {1, 1, 1, 1, 1}}
```

## Algorithm

partitions.c — IntegerPartitions

A faithful, efficient recreation of the Wolfram-Language IntegerPartitions. The whole surface collapses onto a single count-vector enumerator over an ordered set of allowed parts, run with exact GMP rational arithmetic so that integers, big integers, rationals and negative values are all handled by the same code path.

Forms:

```text
  IntegerPartitions[n]                     all partitions of n
  IntegerPartitions[n, k]                  into at most k parts
  IntegerPartitions[n, {k}]                into exactly k parts
  IntegerPartitions[n, {kmin, kmax}]       between kmin and kmax parts
  IntegerPartitions[n, {kmin, kmax, dk}]   kmin, kmin+dk, ... parts
  IntegerPartitions[n, kspec, sspec]       parts drawn only from sspec
  IntegerPartitions[n, kspec, sspec, m]    first m (m>0) or last |m| (m<0)
```

n and the s_i may be rational and/or negative. Results are in reverse lexicographic order; within a partition the parts appear in the order of the reversed sspec (descending for the default Range[n]).

Ownership: this builtin only *reads* `res`. On every NULL return (bad arguments, ::undef, symbolic input) the evaluator keeps `res` unevaluated.

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
- Tests: [`tests/test_integer_partitions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_integer_partitions.c)
