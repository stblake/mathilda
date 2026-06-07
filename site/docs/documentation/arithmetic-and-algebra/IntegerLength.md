# IntegerLength

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
IntegerLength[n] gives the number of decimal digits in the integer n.
IntegerLength[n, b] gives the number of base b digits in n.
IntegerLength ignores the sign of n; IntegerLength[0] is 0.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= IntegerLength[123456789]
Out[1]= 9

In[2]:= IntegerLength[100!, 2]
Out[2]= 525

In[3]:= Table[IntegerLength[100!, n], {n, 2, 20}]
Out[3]= {525, 332, 263, 227, 204, 187, 175, 166, 158, 152, 147, 142, 138, 135, 132, 129, 126, 124, 122}
```

```mathematica
In[1]:= IntegerLength[]
Out[1]= IntegerLength[]

In[2]:= IntegerLength[1, 2, 3, 4]
Out[2]= IntegerLength[1, 2, 3, 4]
```

```mathematica
In[1]:= IntegerLength[1.1234]
Out[1]= IntegerLength[1.1234]
```

## Implementation notes

`builtin_integerlength` returns the number of base-`b` digits of `|n|` (default base 10). For bases `<= 62` it uses GMP's `mpz_sizeinbase`, which is exact for power-of-two bases and at most one too large otherwise — corrected by comparing `|n|` against `base^(s-1)` (`intlen_count_digits`). For arbitrary-precision bases it counts via repeated `mpz_tdiv_q`. `IntegerLength[0]` is `0`. Validates arity (`IntegerLength::argt`), numeric non-integer `n` (`::int`), and base `>= 2` (`::ibase`); symbolic `n` returns `NULL`.

- `Protected`, `Listable`. Threads element-wise over a list of integers in

**Attributes:** `Listable`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/int.c`](https://github.com/stblake/mathilda/blob/main/src/int.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)
