# DigitCount

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
DigitCount[n] gives a list of the counts of digits 1, 2, ..., 9, 0 in the base-10 representation of n.
DigitCount[n, b] gives a list of the counts of digits 1, 2, ..., b-1, 0 in the base-b representation of n.
DigitCount[n, b, d] gives the number of d digits in the base-b representation of n.
The sign of n is discarded; DigitCount[0] is a list of zeros.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= DigitCount[2147, 2, 1]
Out[1]= 5

In[2]:= DigitCount[2147, 2]
Out[2]= {5, 7}

In[3]:= DigitCount[100!]
Out[3]= {15, 19, 10, 10, 14, 19, 7, 14, 20, 30}
```

## Implementation notes

`builtin_digitcount` tallies base-`b` digit occurrences of `|n|` (default base 10). The 3-argument form `DigitCount[n, b, d]` returns the scalar count of digit `d` (`dc_count_one_digit`); the 1/2-argument form returns the histogram `{c(1), c(2), ..., c(b-1), c(0)}` — digit 0 last — built by `dc_build_histogram` into a `calloc`'d `int64` array. Validates arity (`DigitCount::argb`), numeric non-integer `n` (`::int`), base `>= 2` (`::base`, also for fractional bases), digit in `[0, base)` (`::digit`), and caps the list-form base to avoid OOM (`::ovfl`); symbolic `n` returns `NULL`.

- `Protected` (intentionally not `Listable`). `DigitCount[{1,2,3}]` is left

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/int.c`](https://github.com/stblake/mathilda/blob/main/src/int.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= DigitCount[1234567890]
Out[1]= {1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
```

```mathematica
In[1]:= DigitCount[100!]
Out[1]= {15, 19, 10, 10, 14, 19, 7, 14, 20, 30}
```

```mathematica
In[1]:= DigitCount[255, 2]
Out[1]= {8, 0}
```

```mathematica
In[1]:= DigitCount[2^1000, 10, 0]
Out[1]= 28
```

### Notes

`DigitCount[n]` tallies how often each digit appears in the base-10 expansion of
`n`, in the order `1, 2, ..., 9, 0` (note that `0` is reported last). The pandigital
number `1234567890` therefore gives ten ones. The distribution of the 158 digits of
`100!` shows the histogram for a large factorial. With a base argument,
`DigitCount[n, b]` works in any radix — `DigitCount[255, 2] = {8, 0}` recovers the
binary population count (255 is `11111111`). The three-argument form
`DigitCount[n, b, d]` returns just the count of digit `d`; here `2^1000` (a 302-digit
number) contains 28 zeros. The sign of `n` is ignored and `DigitCount[0]` is all
zeros.
