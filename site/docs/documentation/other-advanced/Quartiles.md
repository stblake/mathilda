# Quartiles

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Quartiles[data]
    gives the {q_1/4, q_2/4, q_3/4} quantile estimates of the elements in data.
Quartiles[data,{{a,b},{c,d}}]
    uses the quantile definition specified by parameters a, b, c, d.
Quartiles[dist]
    gives the {q_1/4, q_2/4, q_3/4} quantiles of the distribution dist.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
