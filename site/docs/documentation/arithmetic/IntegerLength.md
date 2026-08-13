# IntegerLength

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`IntegerLength[n] gives the number of decimal digits in the integer n.`**

**`IntegerLength[n, b] gives the number of base b digits in n.`**

<details>
<summary>Notes</summary>

IntegerLength ignores the sign of n; IntegerLength\[0\] is 0.

</details>

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (3)

```mathematica
In[1]:= IntegerLength[123456789]
Out[1]= 9

In[2]:= IntegerLength[100!, 2]
Out[2]= 525

In[3]:= Table[IntegerLength[100!, n], {n, 2, 20}]
Out[3]= {525, 332, 263, 227, 204, 187, 175, 166, 158, 152, 147, 142, 138, 135, 132, 129, 126, 124, 122}
```

### Applications (4)

```mathematica
In[4]:= IntegerLength[12345]
Out[4]= 5

In[5]:= IntegerLength[2^1000]
Out[5]= 302

In[6]:= IntegerLength[2^1000, 2]
Out[6]= 1001

In[7]:= IntegerLength[100!]
Out[7]= 158
```

## Options & behaviour

**Arity diagnostics** (`IntegerLength::argt`). Wrong-arity calls emit the
diagnostic and echo the call back unevaluated, matching Mathematica:

**Non-integer-argument diagnostic** (`IntegerLength::int`). A concrete
non-integer numeric (Real, Rational, Complex) at position 1 or 2 emits the
diagnostic and echoes the call back unevaluated; pure symbolic arguments
flow through silently:

> **Packed arrays.** The one-argument form runs on an `int64` buffer and
> answers with a packed result.

## Implementation notes

`builtin_integerlength` returns the number of base-`b` digits of `|n|` (default base 10). For bases `<= 62` it uses GMP's `mpz_sizeinbase`, which is exact for power-of-two bases and at most one too large otherwise — corrected by comparing `|n|` against `base^(s-1)` (`intlen_count_digits`). For arbitrary-precision bases it counts via repeated `mpz_tdiv_q`. `IntegerLength[0]` is `0`. Validates arity (`IntegerLength::argt`), numeric non-integer `n` (`::int`), and base `>= 2` (`::ibase`); symbolic `n` returns `NULL`.

- `Protected`, `Listable`. Threads element-wise over a list of integers in
  either argument position, e.g. `IntegerLength[{1, 10, 100}]` or
  `IntegerLength[8, {2, 3, 4}]`.
- `IntegerLength[0]` is `0` in any base (zero has no significant digits).
- `IntegerLength[n, b]` is effectively an efficient version of
  `Floor[Log[b, |n|]] + 1` -- it never converts through floating-point and
  is exact for arbitrarily large `n`.
- Works for both machine integers and bignums. The fast path uses GMP's
  `mpz_sizeinbase` for bases `2..62` (with an exact verification step for
  non-power-of-2 bases); arbitrary-precision bases fall back to repeated
  division.

**Attributes:** `Listable`, `Protected`.

## References

- Source: [`src/int.c`](https://github.com/stblake/mathilda/blob/main/src/int.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_digit_count.c`](https://github.com/stblake/mathilda/blob/main/tests/test_digit_count.c)
- Tests: [`tests/test_integer_length.c`](https://github.com/stblake/mathilda/blob/main/tests/test_integer_length.c)
- Tests: [`tests/test_integer_string.c`](https://github.com/stblake/mathilda/blob/main/tests/test_integer_string.c)

## Notes & additional examples

### Notes

`IntegerLength[n]` returns the number of decimal digits of `n`, and
`IntegerLength[n, b]` the number of base-`b` digits. Working at arbitrary
precision, it answers questions that overflow fixed-width arithmetic: `2^1000`
has 302 decimal digits but exactly 1001 binary digits, and `100!` is a 158-digit
number. The sign of `n` is ignored and `IntegerLength[0]` is `0`.
