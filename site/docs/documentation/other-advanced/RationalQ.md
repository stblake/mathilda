# RationalQ

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`RationalQ[expr]`**

gives True if expr is an exact rational number, False otherwise.

**`Rational[p, q]. Returns False on reals and on symbolic expressions.`**

<details>
<summary>Notes</summary>

True for an Integer or BigInt (an integer is rational) and for a

</details>

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

An exact Rational[3, 4]

```mathematica
In[1]:= RationalQ[3/4]
Out[1]= True
```

An integer is rational

```mathematica
In[2]:= RationalQ[7]
Out[2]= True
```

A machine real is not, even with a rational value

```mathematica
In[3]:= RationalQ[2.5]
Out[3]= False
```

An irrational constant is not rational

```mathematica
In[4]:= RationalQ[Pi]
Out[4]= False
```

## Implementation notes

**Algorithm.** `builtin_rationalq` returns `True` when its one argument is an exact
rational number. That is exactly two cases, tested directly: the argument is
integer-like (`expr_is_integer_like` — an `EXPR_INTEGER` or a GMP `EXPR_BIGINT`,
since every integer is rational) or it has head `Rational` (`core_head_is(arg,
SYM_Rational)`, a `Rational[p, q]`). Otherwise `False`.

**Data structures.** None beyond the argument. Two cheap tests — a type tag check
and an interned-pointer head comparison — then a fresh `True`/`False` symbol.

**Complexity / limits.** `O(1)`. Reals and MPFR numbers are *not* rational here
(`2.5` is `False`) even when they have an exact rational value, because the test is
on representation; an irrational constant such as `Pi`, and any non-numeric
expression, is `False`. A call that is not single-argument declines and stays
symbolic.

**Attributes:** `Protected`.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`RationalQ[expr]` is `True` for exactly the exact rationals: any integer (`Integer`
or big `BigInt`) and any `Rational[p, q]`. It is a representation test, so reals and
MPFR numbers are `False` regardless of their value, and irrational constants and
symbolic expressions are `False` as well. Non-single-argument calls are left
unevaluated.
