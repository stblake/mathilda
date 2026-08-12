# Pick

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Pick[expr, sel]`**

Picks out the elements of expr for which the corresponding element of sel is True.

**`Pick[expr, sel, patt]`**

Picks out the elements of expr for which the corresponding element of sel matches patt. Operates at all levels; sel must mirror the structure of expr, and the head of expr is preserved. Returns unevaluated if the structures disagree.

## Examples

_No verified examples yet for this function._

## Algorithm

Pick[expr, sel] / Pick[expr, sel, patt] — select elements of `expr` whose positionally-corresponding element of `sel` matches `patt` (literal `True` for the two-argument form).

The selector array mirrors the structure of `expr`, so the walk is a simultaneous recursive descent over both trees. At each element:

```text
  - selector matches `patt`            -> keep the whole `expr` element,
  - selector is compound and no match  -> recurse into the pair,
  - selector is atomic and no match    -> drop the element.
```

The head at every level comes from `expr`, never from `sel`, so Pick[f[a, b, c], {True, False, True}] is f[a, c].

Any structural disagreement between the two trees (differing lengths, or a compound selector against an atomic expression element) is not an error: the builtin returns NULL and the evaluator leaves `Pick[...]` unevaluated, matching Mathematica's Pick::incomp behaviour of returning the input. Detection is exact — a mismatch found arbitrarily deep aborts the whole call rather than yielding a partially picked result.

## Implementation notes

**Attributes:** `Protected`.

## References

- Source: [`src/list/list_init.c`](https://github.com/stblake/mathilda/blob/main/src/list/list_init.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
- Tests: [`tests/test_list.c`](https://github.com/stblake/mathilda/blob/main/tests/test_list.c)
