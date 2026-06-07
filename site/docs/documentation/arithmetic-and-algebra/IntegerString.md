# IntegerString

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
IntegerString[n] gives a string consisting of the decimal digits in the integer n.
IntegerString[n, b] gives a string consisting of the base-b digits in n; digit values 10 to 35 use the letters a-z.
IntegerString[n, b, len] pads the string on the left with zero digits to give a string of length len; if len is less than the number of digits in n, the len least-significant digits are returned.
The maximum allowed base is 36; the sign of n is discarded.
```

## Examples

All examples below are verified against the current Mathilda build.

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

```mathematica
In[1]:= IntegerString[]
Out[1]= IntegerString[]
```

```mathematica
In[1]:= IntegerString[11.3423]
Out[1]= IntegerString[11.3423]
```

```mathematica
In[1]:= IntegerString[10, 50]
Out[1]= IntegerString[10, 50]
```

## Implementation notes

`builtin_integerstring` renders `|n|` as a base-`b` digit string (default 10, optional left-padded length). It calls GMP's `mpz_get_str`, which emits `'0'..'9','a'..'z'` for bases up to 62; the base is capped at 36 to match Mathematica's surface convention. A length argument left-pads with `'0'` or keeps only the low-order digits. Sign is discarded; `IntegerString[0]` is `"0"`. Validates arity (`IntegerString::argb`), numeric non-integer `n` (`::int`), base in `[2, 36]` (`::basf`), and non-negative machine length (`::intnn`); registered `Listable` so list inputs thread.

- `Protected`, `Listable`. Threads element-wise over a list of integers in

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/int.c`](https://github.com/stblake/mathilda/blob/main/src/int.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)
