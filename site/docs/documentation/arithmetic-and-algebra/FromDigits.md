# FromDigits

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
FromDigits[list] constructs an integer from a list of decimal digits, most-significant first.
FromDigits[list, b] takes the digits to be given in base b.
FromDigits["string"] constructs an integer from a string of digits, where letters a-z and A-Z denote digit values 10-35.
FromDigits["string", b] takes the digits in the string to be given in base b.
Digits in list and characters in the string need not be less than the base; they are carried through Horner's method.  Symbolic digits or base expand to the polynomial sum d[0] b^(n-1) + ... + d[n-1].
```

## Examples

All examples below are verified against the current Mathilda build.

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

## Implementation notes

`builtin_fromdigits` is the inverse of `IntegerDigits`/`IntegerString` (default base 10). For a `List` of digits it Horner-folds `value = value*base + digit` in GMP; symbolic digits or a symbolic base produce a polynomial in the base via the evaluator instead. For a `String` argument each character is mapped to a digit value in `[0, 36)` and folded with an integer base (`FromDigits["abc", b]`); a symbolic/non-integer base over a string is left unevaluated. Validates arity (`FromDigits::argb`) and integer base `>= 2` (`::ibase`).

- `Protected` (intentionally not `Listable`: the first argument *is* a list).
- Inverse of `IntegerDigits` / `IntegerString`. Since `IntegerDigits` discards

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/int.c`](https://github.com/stblake/mathilda/blob/main/src/int.c)
- Specification: [`docs/spec/builtins/arithmetic-and-algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic-and-algebra.md)
