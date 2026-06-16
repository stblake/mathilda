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

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= IntegerDigits[12345]
Out[1]= {1, 2, 3, 4, 5}
```

```mathematica
In[1]:= IntegerDigits[255, 16]
Out[1]= {15, 15}

In[2]:= IntegerDigits[255, 2]
Out[2]= {1, 1, 1, 1, 1, 1, 1, 1}

In[3]:= IntegerDigits[5, 2, 8]
Out[3]= {0, 0, 0, 0, 0, 1, 0, 1}
```

```mathematica
In[1]:= Total[IntegerDigits[2^100]]
Out[1]= 115
```

```mathematica
In[1]:= IntegerDigits[100!, 10][[1 ;; 5]]
Out[1]= {9, 3, 3, 2, 6}
```

### Notes

`IntegerDigits[n]` returns the decimal digits of `n` most-significant first;
`IntegerDigits[n, b]` works in base `b`, and a third argument left-pads (or
truncates to the least-significant) to a fixed length. Because Mathilda carries
arbitrary-precision integers, the digit list of a giant number such as `2^100`
or `100!` is exact — handy for digit-sum problems and divisibility tricks. The
sign of `n` is ignored.
