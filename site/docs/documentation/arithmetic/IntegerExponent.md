# IntegerExponent

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
IntegerExponent[n, b] gives the highest power of b that divides n.
IntegerExponent[n] is equivalent to IntegerExponent[n, 10] and gives the number of trailing zeros in the decimal digits of n.
IntegerExponent ignores the sign of n; IntegerExponent[0, b] is Infinity.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= IntegerExponent[1230000]
Out[1]= 4

In[2]:= IntegerExponent[2^10 + 2^7, 2]
Out[2]= 7

In[3]:= IntegerExponent[144, 2]
Out[3]= 4

In[4]:= IntegerExponent[100!, 2]
Out[4]= 97

In[5]:= IntegerExponent[0]
Out[5]= Infinity
```

```mathematica
In[1]:= IntegerExponent[]
Out[1]= IntegerExponent[]

In[2]:= IntegerExponent[1, 2, 3, 4]
Out[2]= IntegerExponent[1, 2, 3, 4]
```

```mathematica
In[1]:= IntegerExponent[1.123]
Out[1]= IntegerExponent[1.123]
```

## Implementation notes

`builtin_integerexponent` returns the largest `k` with `b^k | n` (the base-`b` valuation, default base 10 = trailing-zero count). For base 2 it uses `mpz_scan1` (position of the lowest set bit = 2-adic valuation); otherwise GMP's `mpz_remove(q, |n|, base)`, which divides out `base` repeatedly and returns the count in one library call (`intexp_count`). `IntegerExponent[0, b]` is `Infinity` (every power divides 0). Validates arity (`::argt`), numeric non-integer `n` (`::int`), and base `>= 2` (`::ibase`); symbolic `n` returns `NULL`.

- `Protected`, `Listable`. Threads element-wise over a list of integers in

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/int.c`](https://github.com/stblake/mathilda/blob/main/src/int.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= IntegerExponent[1000]
Out[1]= 3
```

```mathematica
In[1]:= IntegerExponent[1000, 2]
Out[1]= 3
```

```mathematica
In[1]:= IntegerExponent[100!, 5]
Out[1]= 24
```

```mathematica
In[1]:= IntegerExponent[20!, 2]
Out[1]= 18
```

### Notes

`IntegerExponent[n, b]` gives the highest power of `b` dividing `n`; with no
base it counts trailing decimal zeros. The factorial examples reproduce
Legendre's formula: `100!` ends in `IntegerExponent[100!, 5] = 24` zeros (the
number of factors of five), and `20!` contains `2^18`. The sign of `n` is
ignored and `IntegerExponent[0, b]` is `Infinity`.
