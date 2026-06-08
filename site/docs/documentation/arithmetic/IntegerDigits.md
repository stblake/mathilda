# IntegerDigits

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
IntegerDigits[n] gives a list of the decimal digits in the integer n.
IntegerDigits[n, b] gives a list of the base b digits in the integer n.
IntegerDigits[n, b, len] pads the list on the left with zeros to give a list of length len; if n has more than len base-b digits, the last len least-significant digits are returned.
```

## Examples

All examples below are verified against the current Mathilda build.

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

## Implementation notes

`builtin_integerdigits` returns the base-`b` digit list of `|n|` (default base 10, optional fixed length). It coerces `n` and `base` to GMP, then peels least-significant digits with a `mpz_tdiv_qr` divmod loop into a geometrically-grown `mpz_t` buffer (`O(log_b |n|)`), and emits them most-significant-first into a `List`. A length argument left-pads with zeros (or keeps only the low-order digits if shorter). Validates arity (`IntegerDigits::argb`), non-integer numeric `n` (`::int`), base `>= 2` (`::ibase`), and a non-negative machine-sized length (`::intnn`); symbolic `n` flows through as `NULL`.

- `Protected`, `Listable`. Threading distributes element-wise over a list

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/int.c`](https://github.com/stblake/mathilda/blob/main/src/int.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
