# OddQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
OddQ[n] gives True if n is an odd integer (Integer or BigInt), False otherwise.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_oddq` (`src/core.c`) returns `True` for an `EXPR_INTEGER` with `n % 2 != 0`, uses `mpz_odd_p` for an `EXPR_BIGINT`, and returns `False` for everything else.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= OddQ[7]
Out[1]= True

In[2]:= OddQ[4]
Out[2]= False
```

### Notes

`OddQ` returns `True` only for odd integers (`Integer` or `BigInt`); any non-integer argument gives `False`.
