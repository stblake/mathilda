# FixedPointList

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FixedPointList[f, expr]
    generates the list {expr, f[expr], f[f[expr]], ...} of successive
    applications of f, stopping when two consecutive results are SameQ.
    The last two elements of the result are always the same.
FixedPointList[f, expr, n]
    stops after at most n applications of f. If n is reached before
    convergence, the last two elements may not be equal.
FixedPointList[f, expr, SameTest -> s]
FixedPointList[f, expr, n, SameTest -> s]
    uses the binary predicate s instead of SameQ to test successive pairs.

FixedPointList[f, expr] is equivalent to
NestWhileList[f, expr, UnsameQ, 2].
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= FixedPointList[1 + Floor[#/2] &, 1000]
Out[1]= {1000, 501, 251, 126, 64, 33, 17, 9, 5, 3, 2, 2}

In[2]:= 1 + Floor[Last[%]/2]
Out[2]= 2

In[3]:= FixedPointList[# /. {a_, b_} /; b != 0 -> {b, Mod[a, b]} &, {28, 21}]
Out[3]= {{28, 21}, {21, 7}, {7, 0}, {7, 0}}

In[4]:= GCD[28, 21]
Out[4]= 7

In[5]:= FixedPointList[1 + Floor[#/2] &, 1000, 5]
Out[5]= {1000, 501, 251, 126, 64, 33}

In[6]:= FixedPointList[(# + 2/#)/2 &, 1.0]
Out[6]= {1.0, 1.5, 1.41667, 1.41422, 1.41421, 1.41421, 1.41421}

In[7]:= FixedPointList[(# + 2/#)/2 &, 1.0, SameTest -> (Abs[#1 - #2] < 0.01 &)]
Out[7]= {1.0, 1.5, 1.41667, 1.41422}
```

## Implementation notes

- `Protected`.
- `FixedPointList[f, expr]` is equivalent to `NestWhileList[f, expr, UnsameQ, 2]`.
- `n` (when given) must be a non-negative integer or `Infinity`. Malformed argument specs leave `FixedPointList` unevaluated.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/functional-programming.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/functional-programming.md)
