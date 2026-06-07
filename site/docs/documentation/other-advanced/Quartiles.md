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

**Algorithm.** `builtin_quartiles` returns `{Q1, Q2, Q3}` using a parameterised quantile estimator. For a matrix (list of equal-length lists) it `Transpose`s and recurses column-wise. For a flat list it requires all entries to be real (`is_real_numeric`, else emits `Quartiles::rectn`), sorts via the evaluator's `Sort`, and for each `q` in `{1/4, 1/2, 3/4}` computes the interpolated order statistic with the four quantile parameters `{{a, b}, {c, d}}` — defaulting to `{{1/2, 0}, {0, 1}}` (Mathematica's default). The index is `h = a + (n + b)·q`; it clamps to the endpoints when `h ≤ 1` or `h ≥ n`, takes `j = Floor[h]` and the fractional `g = h − j`, and linearly interpolates `sorted[j] + g·(c + d·g type adjustment)` between adjacent sorted elements. All arithmetic is done by building and `eval_and_free`-ing `Plus`/`Times`/`Floor` expression nodes, so exact (rational/symbolic-numeric) inputs stay exact. Non-list input returns the call unevaluated. Attribute: `ATTR_PROTECTED`.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/stats.c`](https://github.com/stblake/mathilda/blob/main/src/stats.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
