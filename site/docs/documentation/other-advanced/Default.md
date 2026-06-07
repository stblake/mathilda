# Default

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
Default[f]
    gives the default value supplied for a missing optional argument of
    f when the pattern _. (Optional[Blank[]]) appears in a rule.
Default[f, i]
    gives the default value for the i-th argument position of f.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`Default` is not a builtin handler — it is a symbol whose DownValues the user assigns (`Default[f] = v`) and which the pattern matcher consumes when filling an Optional pattern (`x_.`, i.e. `Optional[x_, Default[f]]`). `get_default_value(pat_head, pos, total_pats)` in `src/default_helper.c` looks up the default for a function head: if no user `Default` DownValues exist it short-circuits to the built-in fallbacks (`Plus -> 0`, `Times`/`Power -> 1`, else none); otherwise it evaluates `Default[f, pos, total_pats]`, then `Default[f, pos]`, then `Default[f]`, returning the first that the user has defined (detected by the result differing from the constructed expression), before falling back to the same Plus/Times/Power defaults. The matcher calls this once per Optional-pattern attempt, hence the no-rule fast path.

**Attributes:** none registered.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/default_helper.c`](https://github.com/stblake/mathilda/blob/main/src/default_helper.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
