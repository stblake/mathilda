# ExactNumberQ

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`ExactNumberQ[expr]`**

gives True if expr is an exact number, False otherwise.

<details>
<summary>Notes</summary>

Exact numbers are integers, rationals, and Complex numbers whose parts are exact. Reals and MPFR numbers are inexact, so ExactNumberQ is False.

</details>

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (5)

An exact rational

```mathematica
In[1]:= ExactNumberQ[1/3]
Out[1]= True
```

A Complex with exact parts is exact

```mathematica
In[2]:= ExactNumberQ[2 + 3 I]
Out[2]= True
```

One inexact part makes the whole inexact

```mathematica
In[3]:= ExactNumberQ[2.0 + 3 I]
Out[3]= False
```

A machine real is inexact

```mathematica
In[4]:= ExactNumberQ[1.5]
Out[4]= False
```

A symbol is not a number at all

```mathematica
In[5]:= ExactNumberQ[x]
Out[5]= False
```

## Implementation notes

**Algorithm.** `builtin_exactnumberq` delegates to `core_is_exact_number`, which is
exact iff the argument is built entirely from exact pieces:

- an `EXPR_INTEGER` or `EXPR_BIGINT` — `True`;
- a `Rational[p, q]` — `True`;
- a `Complex[re, im]` — `True` only when **both** `re` and `im` are themselves
  exact (the test recurses into the two parts and ANDs them);
- anything else — `False`.

The complement among numbers is `InexactNumberQ`: a machine `Real` or MPFR number
is inexact, so `ExactNumberQ` is `False` for it.

**Data structures.** None beyond the argument tree; the recursion only descends the
two parts of a `Complex`, so it is bounded.

**Complexity / limits.** `O(1)` (a `Complex` adds two leaf checks). A non-number —
a symbol or symbolic expression such as `Pi` or `x` — is `False`, not left
symbolic, since it is a definite negative answer. A call that is not
single-argument declines and stays symbolic.

**Attributes:** `Protected`.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`ExactNumberQ[expr]` is `True` for integers, rationals, and `Complex` numbers whose
real and imaginary parts are *both* exact. Reals and arbitrary-precision (MPFR)
numbers are inexact, so they are `False` — and a single inexact part makes a whole
`Complex` inexact. It is the complement of `InexactNumberQ` among numbers; a
non-number gives `False`.
