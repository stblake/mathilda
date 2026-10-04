# DigitSum

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`DigitSum[n] gives the sum of the decimal digits of the integer n.`**

**`DigitSum[n, b] gives the sum of the base-b digits of n.`**

<details>
<summary>Notes</summary>

The sign of n is discarded; DigitSum\[0\] is 0.

</details>

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= DigitSum[1234]
Out[1]= 10

In[2]:= DigitSum[255, 16]
Out[2]= 30

In[3]:= DigitSum[{1234, 0, 99}]
Out[3]= {10, 0, 18}
```

### Applications (7)

Sum of the base-10 digits

```mathematica
In[4]:= DigitSum[12345]
Out[4]= 15
```

Base 2: 255 is 11111111

```mathematica
In[5]:= DigitSum[255, 2]
Out[5]= 8
```

Base 16: FF is 15 + 15

```mathematica
In[6]:= DigitSum[255, 16]
Out[6]= 30
```

The sign is discarded

```mathematica
In[7]:= DigitSum[-987]
Out[7]= 24
```

Arbitrary-precision bignum, summed in GMP

```mathematica
In[8]:= DigitSum[2^100]
Out[8]= 115
```

Listable: threads over the list

```mathematica
In[9]:= DigitSum[{123, 4567, 89}]
Out[9]= {6, 22, 17}
```

Casting out nines: n is DigitSum[n] mod 9

```mathematica
In[10]:= Mod[n - DigitSum[n], 9] /. n -> 123456
Out[10]= 0
```

## Implementation notes

**Algorithm.** `builtin_digitsum` takes an integer `n` (a machine `int64` or a
GMP bignum) and an optional base `b >= 2` (default 10). Each argument is
validated first: a concrete-but-non-integer argument (Real, Rational, Complex)
raises `DigitSum::int` / `DigitSum::base` through the `mth_message` funnel and
returns `NULL`, while a symbolic argument is left unevaluated silently so a
downstream rewrite can still rewire the call. The sign of `n` is discarded
(`mpz_abs`), then the base-`b` digits are summed by repeated Euclidean division
with no intermediate digit list: a base that fits in `unsigned long` takes the
fast path `mpz_tdiv_q_ui`, which returns the remainder directly, and a bignum
base takes `mpz_tdiv_qr`. `DigitSum[0]` is `0`.

**Data structures.** All arithmetic is done in GMP `mpz_t`, so machine integers
and arbitrary-precision bignums are handled uniformly on one path. The running
sum accumulates in a single `mpz_t` in one pass over the digits and is demoted
back to an `EXPR_INTEGER` when it fits a signed long, otherwise returned as a
bigint via `expr_new_bigint_from_mpz`.

**Complexity / limits.** `O(d)` divisions for a `d`-digit base-`b`
representation. `DigitSum` is `Listable | NumericFunction | Protected`, so it
threads element-wise over ordinary lists (`DigitSum[{...}]`), but — unlike its
sibling `IntegerDigits` — it has **no** NDArray buffer kernel and no `Compile[]`
lowering: a visible `NDArray` integer argument is left unevaluated, and
`CompileDiagnostics` reports `Compiled -> False`. Semantically it is
`Total[IntegerDigits[n, b]]`.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

- Source: [`src/int.c`](https://github.com/stblake/mathilda/blob/main/src/int.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_digit_sum.c`](https://github.com/stblake/mathilda/blob/main/tests/test_digit_sum.c)

## Notes & additional examples

### Notes

`DigitSum[n]` adds up the base-10 digits of `n`; `DigitSum[n, b]` uses base
`b >= 2`. The sign of `n` is always discarded (`DigitSum[-987]` equals
`DigitSum[987]`), and `DigitSum[0]` is `0`. It is exactly
`Total[IntegerDigits[n, b]]`, computed without building the digit list.

Both machine integers and arbitrary-precision bignums work uniformly, since the
digit extraction runs in GMP. The casting-out-nines identity `n ≡ DigitSum[n]
(mod b-1)` (here `mod 9` for base 10) follows directly from summing digits.

`DigitSum` is `Listable`, so it maps over the elements of a list; it is **not**
accelerated on a packed buffer, and it does not lower in `Compile[]`. Its close
relatives are `IntegerDigits` (the digit list), `DigitCount` (per-value digit
tallies), and `IntegerLength` (how many digits there are).
