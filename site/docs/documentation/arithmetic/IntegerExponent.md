# IntegerExponent

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`IntegerExponent[n, b] gives the highest power of b that divides n.`**

**`IntegerExponent[n] is equivalent to IntegerExponent[n, 10] and gives the number of trailing zeros in the decimal digits of n.`**

<details>
<summary>Notes</summary>

IntegerExponent ignores the sign of n; IntegerExponent\[0, b\] is Infinity.

</details>

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

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

### Scope (2)

```mathematica
In[6]:= IntegerExponent[] IntegerExponent::argt: IntegerExponent called with 0 arguments; 1 or 2 arguments are expected.

In[7]:= IntegerExponent[1, 2, 3, 4] IntegerExponent::argt: IntegerExponent called with 4 arguments; 1 or 2 arguments are expected.
```

### Applications (4)

```mathematica
In[8]:= IntegerExponent[1000]
Out[8]= 3

In[9]:= IntegerExponent[1000, 2]
Out[9]= 3

In[10]:= IntegerExponent[100!, 5]
Out[10]= 24

In[11]:= IntegerExponent[20!, 2]
Out[11]= 18
```

## Options & behaviour

**Arity diagnostics** (`IntegerExponent::argt`). Wrong-arity calls emit the
diagnostic and echo the call back unevaluated, matching Mathematica:

**Non-integer-argument diagnostic** (`IntegerExponent::int`). A concrete
non-integer numeric (Real, Rational, Complex) at position 1 or 2 emits the
diagnostic and echoes the call back unevaluated; pure symbolic arguments
flow through silently:

## Implementation notes

`builtin_integerexponent` returns the largest `k` with `b^k | n` (the base-`b` valuation, default base 10 = trailing-zero count). For base 2 it uses `mpz_scan1` (position of the lowest set bit = 2-adic valuation); otherwise GMP's `mpz_remove(q, |n|, base)`, which divides out `base` repeatedly and returns the count in one library call (`intexp_count`). `IntegerExponent[0, b]` is `Infinity` (every power divides 0). Validates arity (`::argt`), numeric non-integer `n` (`::int`), and base `>= 2` (`::ibase`); symbolic `n` returns `NULL`.

- `Protected`, `Listable`. Threads element-wise over a list of integers in
  either argument position, e.g. `IntegerExponent[{10, 100, 1000}]` or
  `IntegerExponent[24, {2, 3, 4, 6}]`.
- Sign of `n` is discarded.
- `IntegerExponent[0, b]` is `Infinity` for any base `b` (every power of `b`
  divides 0).
- Works for both machine integers and bignums. Base 2 uses GMP's
  `mpz_scan1` (lowest-set-bit lookup, O(log n / wordsize)); other bases use
  `mpz_remove` (a single library call that strips all factors of `b`).

**Attributes:** `Listable`, `Protected`.

## References

- Source: [`src/int.c`](https://github.com/stblake/mathilda/blob/main/src/int.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_integer_exponent.c`](https://github.com/stblake/mathilda/blob/main/tests/test_integer_exponent.c)

## Notes & additional examples

### Notes

`IntegerExponent[n, b]` gives the highest power of `b` dividing `n`; with no
base it counts trailing decimal zeros. The factorial examples reproduce
Legendre's formula: `100!` ends in `IntegerExponent[100!, 5] = 24` zeros (the
number of factors of five), and `20!` contains `2^18`. The sign of `n` is
ignored and `IntegerExponent[0, b]` is `Infinity`.
