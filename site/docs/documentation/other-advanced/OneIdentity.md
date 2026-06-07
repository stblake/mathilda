# OneIdentity

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
OneIdentity is an attribute that can be assigned to a symbol f to indicate that f[x], f[f[x]], etc. are all equivalent to x for the purpose of pattern matching.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`OneIdentity` is an attribute marker symbol, not a function. `attr.c` maps the name to the `ATTR_ONEIDENTITY` bitflag. The bit is consumed by the pattern matcher: in `match.c`, when matching a single argument against a sequence/optional pattern under a head carrying `ATTR_ONEIDENTITY`, `f[x]` is treated as equivalent to `x` (and an optional/defaulted argument lets `x` match `f[x, default]`). This is what lets `a + b_.` match `a` (with `b` bound to the identity 0) for `Plus`, or `a x_.` match `a` for `Times`. The matcher checks `def->attributes & ATTR_ONEIDENTITY` at the relevant decision points. `purefunc.c` maps `SYM_OneIdentity` → `ATTR_ONEIDENTITY` for pure-function attribute specs.

**Attributes:** none registered.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/match.c`](https://github.com/stblake/mathilda/blob/main/src/match.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)
