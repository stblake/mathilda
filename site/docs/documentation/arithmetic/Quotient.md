# Quotient

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Quotient[m, n]`**

gives the integer quotient of m and n, rounded toward -Infinity (floored division), so that n Quotient\[m, n\] + Mod\[m, n\] == m.

**`Quotient[m, n, d]`**

uses the offset d, matching the three-argument Mod.

<details>
<summary>Notes</summary>

For complex arguments Quotient performs Gaussian-integer division, rounding the ratio to the nearest Gaussian integer. Exact inputs give exact results; Quotient is Listable.

</details>

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

```mathematica
In[1]:= Quotient[17, 5]
Out[1]= 3

In[2]:= Quotient[-17, 5]
Out[2]= -4

In[3]:= 5 Quotient[-17, 5] + Mod[-17, 5]
Out[3]= -17

In[4]:= Quotient[5 + 3 I, 2]
Out[4]= 2 + 2 I
```

## Implementation notes

`builtin_quotient` floors the ratio for real arguments (`mpz_fdiv_q` on exact integers), so that `n Quotient[m, n] + Mod[m, n] == m` holds exactly; the three-argument form applies the same `d`-offset as `Mod`. For complex `m` or `n` it switches to Gaussian-integer division — it forms the exact ratio and rounds each component to the nearest integer (ties to even), which is the quotient minimising the norm of the remainder and deliberately differs from the real (floored) branch. Registered `PROTECTED | NUMERICFUNCTION | LISTABLE`; non-numeric arguments stay symbolic.

**Attributes:** `Listable`, `NumericFunction`, `Protected`.

## References

**See also:** [Mod](../../arithmetic/Mod/), [QuotientRemainder](../../arithmetic/QuotientRemainder/), [Union](../../structural-manipulation/Union/), [Tally](../../data-structures/Tally/), [DeleteDuplicates](../../data-structures/DeleteDuplicates/)

- Knuth, "The Art of Computer Programming, Vol. 2: Seminumerical Algorithms", on the division algorithm.
- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/arithmetic.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/arithmetic.md)
- Tests: [`tests/test_compile.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile.c)
- Tests: [`tests/test_compile_coverage.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compile_coverage.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_core.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core.c)

## Notes & additional examples

### Notes

`Quotient[m, n]` floors the ratio toward `-Infinity`, so `Quotient[-17, 5] = -4`
(not `-3`) and the division identity `n Quotient[m, n] + Mod[m, n] == m` holds
exactly. For complex arguments it is Gaussian-integer division, rounding each part
of the ratio to the nearest integer, so `Quotient[5 + 3 I, 2] = 2 + 2 I` — the
Gaussian integer nearest the ratio `2.5 + 1.5 I`. The three-argument
`Quotient[m, n, d]` uses the same offset convention as the three-argument `Mod`.
