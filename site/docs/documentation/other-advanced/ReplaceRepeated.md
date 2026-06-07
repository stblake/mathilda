# ReplaceRepeated

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
expr //. rules or ReplaceRepeated[expr, rules]
    repeatedly applies ReplaceAll[expr, rules] until the result stops
    changing, then returns the fixed point.  Useful for chained rewrite
    systems; subject to the same recursion-limit guard as evaluator
    fixed-point iteration.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

**Algorithm.** `builtin_replace_repeated` (`//.`) is a two-argument wrapper over `apply_replace_repeated_nested`. If the rule argument is a list of *lists* of rules (no top-level `Rule`), it threads over them, producing a list of independently-repeated results. Otherwise it iterates a fixed-point loop: each pass applies `apply_replace_all_nested` (the same top-down `ReplaceAll` traversal used by `/.`), then `eval_and_free`s the result, and compares against the previous expression with `expr_eq`. It returns when the expression stops changing, with a hard cap of 65536 iterations as a divergence guard. The input is deep-copied up front so the original is not mutated.

**Attributes:** `Protected`.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/replace.c`](https://github.com/stblake/mathilda/blob/main/src/replace.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
