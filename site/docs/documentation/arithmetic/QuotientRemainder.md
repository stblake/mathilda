# QuotientRemainder

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`QuotientRemainder[m, n]`**

gives the pair {Quotient\[m, n\], Mod\[m, n\]}, so the quotient is floored and the remainder carries the sign of n.

<details>
<summary>Notes</summary>

QuotientRemainder is Listable; non-numeric arguments are left unevaluated.

</details>

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

```mathematica
In[1]:= QuotientRemainder[17, 5]
Out[1]= {3, 2}

In[2]:= QuotientRemainder[-17, 5]
Out[2]= {-4, 3}

In[3]:= QuotientRemainder[17, -5]
Out[3]= {-4, -3}
```

## Implementation notes

`builtin_quotientremainder` returns the pair `{Quotient[m, n], Mod[m, n]}`, sharing the floored-division and Gaussian-integer conventions of its two components, so the remainder always carries the sign of the divisor. Registered `PROTECTED | NUMERICFUNCTION | LISTABLE`; non-numeric arguments are left unevaluated.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Mod](../../arithmetic/Mod/), [Quotient](../../arithmetic/Quotient/), [Union](../../structural-manipulation/Union/), [Tally](../../data-structures/Tally/), [DeleteDuplicates](../../data-structures/DeleteDuplicates/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_core.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core.c)
- Tests: [`tests/test_eval.c`](https://github.com/stblake/mathilda/blob/main/tests/test_eval.c)
- Tests: [`tests/test_expr_sharing.c`](https://github.com/stblake/mathilda/blob/main/tests/test_expr_sharing.c)
- Tests: [`tests/test_ndsolve_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_ndsolve_compile.c)

## Notes & additional examples

### Notes

`QuotientRemainder[m, n]` returns the quotient and remainder together as
`{Quotient[m, n], Mod[m, n]}`. Because the quotient is floored and the remainder
takes the sign of the divisor `n`, the two always reconstruct the dividend:
`n q + r == m`. So `QuotientRemainder[-17, 5] = {-4, 3}` (a non-negative
remainder) while `QuotientRemainder[17, -5] = {-4, -3}` (a non-positive one).
