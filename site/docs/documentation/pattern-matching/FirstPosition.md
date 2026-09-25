# FirstPosition

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`FirstPosition[expr, pattern] gives the position of the first element in expr matching pattern (in depth-first order), or Missing["NotFound"] if no such element is found.`**

**`FirstPosition[expr, pattern, default] gives default if no element matching pattern is found; default is evaluated only when it is returned.`**

**`FirstPosition[expr, pattern, default, levelspec] finds only objects on the levels specified by levelspec.`**

<details>
<summary>Notes</summary>

Over an association, the position is a key, e.g. {Key\[k\], ...}. The default level specification is {0, Infinity} with Heads -\> True; a position of {} represents the whole of expr.

</details>

## Examples

_No verified examples yet for this function._

## Implementation notes

- Attributes `{HoldRest, Protected}`; `default` is evaluated only when it is returned.
- Same defaults as `Position`: levels `{0, Infinity}` with `Heads -> True`; a position of `{}` is the whole of `expr`.
- Over an association the position is a key, e.g. `{Key[k], ...}`.
- Delegates to `Position` (with the first-match cap), so traversal, level specs, the `Heads` option, and association handling match `Position` exactly.

**Attributes:** `HoldRest`, `Protected`.

## References

**See also:** [Position](../../data-structures/Position/)

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/pattern-matching.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/pattern-matching.md)
- Tests: [`tests/test_firstposition.c`](https://github.com/stblake/mathilda/blob/main/tests/test_firstposition.c)
