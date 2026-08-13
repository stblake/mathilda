# IntegerDigits

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`IntegerDigits[n] gives a list of the decimal digits in the integer n.`**

**`IntegerDigits[n, b] gives a list of the base b digits in the integer n.`**

**`IntegerDigits[n, b, len] pads the list on the left with zeros to give a list of length len; if n has more than len base-b digits, the last len least-significant digits are returned.`**

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (4)

```mathematica
In[1]:= IntegerDigits[58127]
Out[1]= {5, 8, 1, 2, 7}

In[2]:= IntegerDigits[58127, 16]
Out[2]= {14, 3, 0, 15}

In[3]:= IntegerDigits[Range[0, 7], 2, 3]
Out[3]= {{0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1}, {1, 0, 0}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1}}

In[4]:= IntegerDigits[6345354, 10, 4]
Out[4]= {5, 3, 5, 4}
```

### Applications (6)

```mathematica
In[5]:= IntegerDigits[12345]
Out[5]= {1, 2, 3, 4, 5}

In[6]:= IntegerDigits[255, 16]
Out[6]= {15, 15}

In[7]:= IntegerDigits[255, 2]
Out[7]= {1, 1, 1, 1, 1, 1, 1, 1}

In[8]:= IntegerDigits[5, 2, 8]
Out[8]= {0, 0, 0, 0, 0, 1, 0, 1}

In[9]:= Total[IntegerDigits[2^100]]
Out[9]= 115

In[10]:= IntegerDigits[100!, 10][[1 ;; 5]]
Out[10]= {9, 3, 3, 2, 6}
```

## Options & behaviour

> **Packed arrays.** The one-argument form reads an `int64` buffer directly
> rather than materialising it. The result is ragged — a list per element, of
> differing lengths — so it is an ordinary list either way.

## Implementation notes

`builtin_integerdigits` returns the base-`b` digit list of `|n|` (default base 10, optional fixed length). It coerces `n` and `base` to GMP, then peels least-significant digits with a `mpz_tdiv_qr` divmod loop into a geometrically-grown `mpz_t` buffer (`O(log_b |n|)`), and emits them most-significant-first into a `List`. A length argument left-pads with zeros (or keeps only the low-order digits if shorter). Validates arity (`IntegerDigits::argb`), non-integer numeric `n` (`::int`), base `>= 2` (`::ibase`), and a non-negative machine-sized length (`::intnn`); symbolic `n` flows through as `NULL`.

- `Protected`, `Listable`. Threading distributes element-wise over a list
  of integers in any argument position, e.g. `IntegerDigits[{6, 7, 2}, 2]`
  and `IntegerDigits[7, {2, 3, 4}]`.
- `IntegerDigits[0]` returns `{0}` (single zero digit). With an explicit
  length, `IntegerDigits[0, b, len]` returns a list of `len` zeros.
- Bases > 10 are allowed; digit values range over `{0, ..., b-1}`.
- Works seamlessly for both machine integers and arbitrary-precision
  bignums (digits are computed in GMP and demoted to machine ints
  whenever they fit).

**Attributes:** `Listable`, `Protected`.

## References

- Source: [`src/int.c`](https://github.com/stblake/mathilda/blob/main/src/int.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_digit_count.c`](https://github.com/stblake/mathilda/blob/main/tests/test_digit_count.c)
- Tests: [`tests/test_digit_sum.c`](https://github.com/stblake/mathilda/blob/main/tests/test_digit_sum.c)
- Tests: [`tests/test_from_digits.c`](https://github.com/stblake/mathilda/blob/main/tests/test_from_digits.c)
- Tests: [`tests/test_integer_digits.c`](https://github.com/stblake/mathilda/blob/main/tests/test_integer_digits.c)

## Notes & additional examples

### Notes

`IntegerDigits[n]` returns the decimal digits of `n` most-significant first;
`IntegerDigits[n, b]` works in base `b`, and a third argument left-pads (or
truncates to the least-significant) to a fixed length. Because Mathilda carries
arbitrary-precision integers, the digit list of a giant number such as `2^100`
or `100!` is exact — handy for digit-sum problems and divisibility tricks. The
sign of `n` is ignored.
