# IntegerString

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`IntegerString[n] gives a string consisting of the decimal digits in the integer n.`**

**`IntegerString[n, b] gives a string consisting of the base-b digits in n; digit values 10 to 35 use the letters a-z.`**

**`IntegerString[n, b, len] pads the string on the left with zero digits to give a string of length len; if len is less than the number of digits in n, the len least-significant digits are returned.`**

<details>
<summary>Notes</summary>

The maximum allowed base is 36; the sign of n is discarded.

</details>

## Examples (11)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (5)

```mathematica
In[1]:= IntegerString[17651, 2]
Out[1]= "100010011110011"

In[2]:= IntegerString[50!, 16]
Out[2]= "49eebc961ed279b02b1ef4f28d19a84f5973a1d2c7800000000000"

In[3]:= IntegerString[50!, 36]
Out[3]= "4q7eyp9zizmtqt0648txt4fm720cc1s00000000000"

In[4]:= IntegerString[Range[0, 7], 2, 3]
Out[4]= {"000", "001", "010", "011", "100", "101", "110", "111"}

In[5]:= IntegerString[12345, 10, 3]
Out[5]= "345"
```

### Scope (1)

```mathematica
In[6]:= IntegerString[10, 50] IntegerString::basf: 50 is not a valid base for IntegerString in IntegerString[10, 50]; the base must be an integer between 2 and 36.
Out[6]= 72.0 an and base be between integer must the
```

### Applications (5)

```mathematica
In[7]:= IntegerString[255, 16]
Out[7]= "ff"

In[8]:= IntegerString[255, 2]
Out[8]= "11111111"

In[9]:= IntegerString[42, 2, 16]
Out[9]= "0000000000101010"

In[10]:= IntegerString[123456789, 36]
Out[10]= "21i3v9"

In[11]:= IntegerString[3^50, 16]
Out[11]= "980553f0db2fd09de3c9"
```

## Options & behaviour

**Arity diagnostics** (`IntegerString::argb`). Wrong-arity calls emit the
diagnostic and echo the call back unevaluated, matching Mathematica:

**Non-integer-argument diagnostic** (`IntegerString::int`). A concrete
non-integer numeric (Real, Rational, Complex) at position 1 or 3 emits the
diagnostic and echoes the call back unevaluated; pure symbolic arguments
flow through silently:

**Invalid-base diagnostic** (`IntegerString::basf`). Bases outside the
`[2, 36]` range (integer or otherwise) trigger the diagnostic and the call
is left unevaluated:

## Implementation notes

`builtin_integerstring` renders `|n|` as a base-`b` digit string (default 10, optional left-padded length). It calls GMP's `mpz_get_str`, which emits `'0'..'9','a'..'z'` for bases up to 62; the base is capped at 36 to match the conventional surface radix. A length argument left-pads with `'0'` or keeps only the low-order digits. Sign is discarded; `IntegerString[0]` is `"0"`. Validates arity (`IntegerString::argb`), numeric non-integer `n` (`::int`), base in `[2, 36]` (`::basf`), and non-negative machine length (`::intnn`); registered `Listable` so list inputs thread.

- `Protected`, `Listable`. Threads element-wise over a list of integers in
  any argument position, e.g. `IntegerString[Range[0, 7], 2, 3]` returns
  a length-8 list of zero-padded binary strings.
- Inverse of `FromDigits` for the integer case:
  `FromDigits[IntegerString[n, b], b] == n` for any non-negative `n`.
- Defining property: `StringLength[IntegerString[n, b]] == IntegerLength[n, b]`
  for `n != 0`.
- `IntegerString[0]` returns `"0"`; `IntegerString[n, b, 0]` returns `""`.
- Works seamlessly for machine integers and arbitrary-precision bignums --
  the digit rendering goes through GMP's `mpz_get_str` in a single call,
  so even `IntegerString[100!, 36]` is essentially free.

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [FromDigits](../../arithmetic/FromDigits/)

- Source: [`src/int.c`](https://github.com/stblake/mathilda/blob/main/src/int.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_integer_string.c`](https://github.com/stblake/mathilda/blob/main/tests/test_integer_string.c)

## Notes & additional examples

### Notes

`IntegerString[n, b]` renders `n` as a base-`b` string, using the digits `0-9`
and then the letters `a-z` for values 10 through 35 (so the maximum base is 36).
A third argument left-pads with zeros to a fixed width. Because the conversion
runs at arbitrary precision, large numbers like `3^50` are formatted exactly in
hexadecimal. The sign of `n` is discarded.
