# IntegerPartitions

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`IntegerPartitions[n] gives the partitions of the integer n -- the ways to write n as a sum of positive parts (equivalently, its Young diagrams) -- in reverse-lexicographic order.`**

**`IntegerPartitions[n, k] gives partitions into at most k parts;`**

**`PartitionsP[n] for the plain form.`**

<details>
<summary>Notes</summary>

{k} exactly k; {kmin, kmax} between; {kmin, kmax, dk} stepped. A third argument restricts the parts (sspec; All = Range\[n\]); a fourth limits the result to the first m (m\>0) or last |m| (m\<0). n and the parts may be rational and negative; Length equals

</details>

## Examples (7)

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

### Worked examples (2)

```mathematica
In[4]:= IntegerPartitions[1/2]
Out[4]= {}

In[5]:= IntegerPartitions[1/2, All, {1/6, 1/3}]
Out[5]= {{1/3, 1/6}, {1/6, 1/6, 1/6}}
```

### Applications (2)

```mathematica
In[6]:= IntegerPartitions[4]
Out[6]= {{4}, {3, 1}, {2, 2}, {2, 1, 1}, {1, 1, 1, 1}}

In[7]:= Length[IntegerPartitions[10]]
Out[7]= 42
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

**Algorithm.** `builtin_integerpartitions` collapses all of its call forms —
`IntegerPartitions[n]`, `[n, k]`, `[n, {k}]`, `[n, {kmin, kmax}]`, `[n, {kmin, kmax, dk}]`,
an optional parts set `sspec`, and a final count `m` — onto one recursive count-vector
enumerator (`ip_recurse`) over an ordered, reversed set of allowed parts. The default
`sspec` is `Range[n]` reversed (`floor(n), …, 1`); the recursion assigns a count to each
part in descending order, pruned by a length budget from a finite `kmax` and a tight numeric
bound (`ip_floor_div`) that is valid only when the remaining parts are single-signed.
Complete partitions — `remaining == 0` with a length inside `[kmin, kmax]` stepped by `dk` —
are emitted in reverse-lexicographic order, and a final `m` slices the first `m` (or, when
negative, the last `|m|`) of them.

**Data structures.** Everything runs in exact GMP rationals (`mpq_t`), so `n` and each
`s_i` may be an integer, big integer, rational or negative value handled by one code path.
The enumeration context `PartCtx` holds the reversed parts array, precomputed suffix-sign
flags (which license the numeric bound), the current count vector, and a growable `Expr**`
of emitted `List` partitions. The builtin only reads `res`; a `NULL` return (bad arguments,
symbolic input, infinite result) leaves the call unevaluated. There is no ND/packed/`Compile`
path — it is a structural list generator, and `IntegerPartitions` is `Protected` but not
`Listable`.

**Complexity / limits.** Output-sensitive: proportional to the number of partitions
enumerated, times `O(nparts)` per emission. An infinite result (a part of `0`, or mixed
signs with an unbounded part count) is detected up front and raises `IntegerPartitions::undef`;
a `take` count exceeding the number found warns `IntegerPartitions::take`; `0` or more than
`4` arguments emit `IntegerPartitions::argb`. Symbolic or real `n` / `s_i` leave the call
unevaluated.

**Attributes:** `Protected`.

## References

- G. E. Andrews, *The Theory of Partitions*, Cambridge University Press, 1998 — the standard reference on partitions and their generating functions.
- Source: [`src/partitions.c`](https://github.com/stblake/mathilda/blob/main/src/partitions.c)
- Specification: [`docs/spec/builtins/number-theory.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/number-theory.md)
- Tests: [`tests/test_integer_partitions.c`](https://github.com/stblake/mathilda/blob/main/tests/test_integer_partitions.c)

## Notes & additional examples

### Partitions of an integer

A *partition* of `n` is a way of writing it as a sum of positive integers, order disregarded
— pictured as a *Young diagram* of left-justified rows. `IntegerPartitions[n]` lists them in
reverse-lexicographic order; the second and later arguments restrict the number of parts and
the allowed parts. The count of unrestricted partitions is [`PartitionsP`](PartitionsP.md),
so `Length[IntegerPartitions[n]] == PartitionsP[n]`.
