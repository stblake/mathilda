# EvenQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`EvenQ[n] gives True if n is an even integer (Integer or BigInt), False otherwise.`**

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (2)

```mathematica
In[1]:= EvenQ[4]
Out[1]= True

In[2]:= EvenQ[7]
Out[2]= False
```

## Implementation notes

`builtin_evenq` (`src/core.c`) returns `True` for an `EXPR_INTEGER` with `n % 2 == 0`, uses `mpz_even_p` for an `EXPR_BIGINT`, and returns `False` for everything else.

**Attributes:** `Protected`.

## See also

[OddQ](../../expression-information/OddQ/)

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
- Tests: [`tests/test_backtrack.c`](https://github.com/stblake/mathilda/blob/main/tests/test_backtrack.c)
- Tests: [`tests/test_bigint.c`](https://github.com/stblake/mathilda/blob/main/tests/test_bigint.c)
- Tests: [`tests/test_compiledfunction.c`](https://github.com/stblake/mathilda/blob/main/tests/test_compiledfunction.c)
- Tests: [`tests/test_core.c`](https://github.com/stblake/mathilda/blob/main/tests/test_core.c)

## Notes & additional examples

### Notes

`EvenQ` returns `True` only for even integers (`Integer` or `BigInt`); any non-integer argument gives `False`.
