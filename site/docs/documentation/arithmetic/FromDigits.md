# FromDigits

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FromDigits[list] constructs an integer from a list of decimal digits, most-significant first.`**

**`FromDigits[list, b] takes the digits to be given in base b.`**

**`FromDigits["string"] constructs an integer from a string of digits, where letters a-z and A-Z denote digit values 10-35.`**

**`FromDigits["string", b] takes the digits in the string to be given in base b.`**

<details>
<summary>Notes</summary>

Digits in list and characters in the string need not be less than the base; they are carried through Horner's method.  Symbolic digits or base expand to the polynomial sum d\[0\] b^(n-1) + ... + d\[n-1\].

</details>

## Examples (10)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= FromDigits[{5, 1, 2, 8}]
Out[1]= 5128

In[2]:= FromDigits[{1, 0, 1, 1, 0, 1, 1}, 2]
Out[2]= 91

In[3]:= FromDigits["1A3C"]
Out[3]= 2042

In[4]:= FromDigits[IntegerDigits[2^100]]
Out[4]= 1267650600228229401496703205376

In[5]:= FromDigits[{a, b, c, d, e}, x]
Out[5]= e + d x + c x^2 + b x^3 + a x^4
```

### Applications (5)

```mathematica
In[6]:= FromDigits[{1, 2, 3, 4}]
Out[6]= 1234

In[7]:= FromDigits[{1, 0, 1, 1}, 2]
Out[7]= 11

In[8]:= FromDigits["deadbeef", 16]
Out[8]= 3735928559

In[9]:= FromDigits[IntegerDigits[2^100], 10]
Out[9]= 1267650600228229401496703205376

In[10]:= FromDigits[{d2, d1, d0}, b]
Out[10]= d0 + b d1 + b^2 d2
```

## Implementation notes

`builtin_fromdigits` is the inverse of `IntegerDigits`/`IntegerString` (default base 10). For a `List` of digits it Horner-folds `value = value*base + digit` in GMP; symbolic digits or a symbolic base produce a polynomial in the base via the evaluator instead. For a `String` argument each character is mapped to a digit value in `[0, 36)` and folded with an integer base (`FromDigits["abc", b]`); a symbolic/non-integer base over a string is left unevaluated. Validates arity (`FromDigits::argb`) and integer base `>= 2` (`::ibase`).

- `Protected` (intentionally not `Listable`: the first argument *is* a list).
- Inverse of `IntegerDigits` / `IntegerString`. Since `IntegerDigits` discards
  the sign of `n`, `FromDigits[IntegerDigits[n]]` is `Abs[n]`, not `n`.
- Digits in the list (and characters in the string) need *not* be smaller
  than the base; they are carried via Horner's evaluation, matching
  Mathematica (e.g. `FromDigits[{7, 11, 0, 0, 0, 122}] == 810122`,
  `FromDigits["1A3C"] == 2042`).
- The all-integer (digits + base) case is computed exactly in GMP and
  demoted to a machine integer when it fits in `int64`; arbitrarily large
  bignums are returned otherwise.
- Symbolic, Real, or Rational digits or base trigger the polynomial Horner
  expansion `d[0] b^(n-1) + d[1] b^(n-2) + ... + d[n-1]`, which the
  evaluator simplifies normally. This single code path handles symbolic
  bases, symbolic digits, and inexact bases uniformly. For a string input,
  the base is required to be an integer; symbolic bases over strings leave
  the call unevaluated (silently).
- Edge cases: `FromDigits[{}] == 0`, `FromDigits[""] == 0`,
  `FromDigits[{0, 0, 1, 2, 3}] == 123` (leading zeros are inert).

**Attributes:** `Protected`.

## References

**See also:** [IntegerDigits](../../arithmetic/IntegerDigits/), [IntegerString](../../arithmetic/IntegerString/)

- Source: [`src/int.c`](https://github.com/stblake/mathilda/blob/main/src/int.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_from_digits.c`](https://github.com/stblake/mathilda/blob/main/tests/test_from_digits.c)
- Tests: [`tests/test_integer_string.c`](https://github.com/stblake/mathilda/blob/main/tests/test_integer_string.c)

## Notes & additional examples

### Notes

`FromDigits` evaluates a digit list most-significant-first via Horner's method.
A second argument fixes the base: `{1, 0, 1, 1}` in base 2 is `11`, and the hex
string `"deadbeef"` (letters `a`-`z` denoting 10-35) is `3735928559`. Because the
arithmetic is exact arbitrary-precision, round-tripping a 31-digit number through
`IntegerDigits`/`FromDigits` reproduces it exactly. With symbolic digits or base,
the Horner recurrence is returned as a polynomial, `d0 + b d1 + b^2 d2` — the
inverse construction to `IntegerDigits`.
