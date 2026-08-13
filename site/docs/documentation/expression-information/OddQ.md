# OddQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`OddQ[n] gives True if n is an odd integer (Integer or BigInt), False otherwise.`**

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (2)

```mathematica
In[1]:= OddQ[7]
Out[1]= True

In[2]:= OddQ[4]
Out[2]= False
```

## Implementation notes

`builtin_oddq` (`src/core.c`) returns `True` for an `EXPR_INTEGER` with `n % 2 != 0`, uses `mpz_odd_p` for an `EXPR_BIGINT`, and returns `False` for everything else.

**Attributes:** `Protected`.

## References

**See also:** [EvenQ](../../expression-information/EvenQ/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_bigint.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bigint.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_core.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core.c)

## Notes & additional examples

### Notes

`OddQ` returns `True` only for odd integers (`Integer` or `BigInt`); any non-integer argument gives `False`.
