# InexactNumberQ

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`InexactNumberQ[expr]`**

gives True if expr is an inexact number, False otherwise.

<details>
<summary>Notes</summary>

Inexact numbers are machine reals, arbitrary-precision (MPFR) reals, and Complex numbers with an inexact part. The complement of ExactNumberQ among numbers.

</details>

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

A machine real is inexact

```mathematica
In[1]:= InexactNumberQ[1.5]
Out[1]= True
```

One inexact part makes the Complex inexact

```mathematica
In[2]:= InexactNumberQ[2.0 + 3 I]
Out[2]= True
```

An exact rational is not inexact

```mathematica
In[3]:= InexactNumberQ[1/3]
Out[3]= False
```

An all-exact Complex is not inexact

```mathematica
In[4]:= InexactNumberQ[2 + 3 I]
Out[4]= False
```

## Implementation notes

**Algorithm.** `builtin_inexactnumberq` delegates to `core_is_inexact_number`,
which is the mirror of `core_is_exact_number`:

- a machine `EXPR_REAL` — `True`;
- an arbitrary-precision `EXPR_MPFR` number (when built `USE_MPFR`) — `True`;
- a `Complex[re, im]` — `True` when **either** part is inexact (the test recurses
  into the two parts and ORs them);
- anything else — `False`.

So an exact number (integer, rational, or an all-exact `Complex`) is `False`, and a
`Complex` with even one floating-point part is `True`. Among numbers it is exactly
the complement of `ExactNumberQ`.

**Data structures.** None beyond the argument; the recursion descends only the two
parts of a `Complex`.

**Complexity / limits.** `O(1)`. A non-number (symbol or symbolic expression) is
`False`. A call that is not single-argument declines and stays symbolic.

**Attributes:** `Protected`.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`InexactNumberQ[expr]` is `True` for machine reals, arbitrary-precision (MPFR)
reals, and `Complex` numbers with at least one inexact part. It is the complement
of `ExactNumberQ` among numbers, so integers, rationals and all-exact `Complex`
numbers are `False`, as is any non-number. Non-single-argument calls are left
unevaluated.
