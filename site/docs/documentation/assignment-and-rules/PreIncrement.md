# PreIncrement

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
PreIncrement[x] or ++x
    increases the value of x by 1, returning the new value of x.
    ++x is equivalent to x = x + 1.

PreIncrement has attribute HoldFirst. In PreIncrement[x], x can be a
symbol or a Part expression referring to an existing value.
If x has no assigned value, PreIncrement::rvalue is emitted and the
expression is left unevaluated.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_preincrement` (`++x`) calls the shared `increment_core(lhs, 1, negate=false, pre=true, "PreIncrement")`. The `pre=true` flag makes it return the *new* value: it reads the old OwnValue via `evaluate`, computes and writes back `Plus[old, 1]` through an evaluated `Set` (which preserves lvalue shape via Set's `HoldFirst`), then returns the new value and frees the old. Requires the target to already have an OwnValue, else emits `PreIncrement::rvalue` and returns `NULL`. Carries `ATTR_HOLDFIRST | ATTR_PROTECTED`.

**Attributes:** `HoldFirst`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= j = 5
Out[1]= 5

In[2]:= ++j
Out[2]= 6

In[3]:= j
Out[3]= 6
```

### Notes

`++j` (`PreIncrement`) increases `j` by 1 and returns the *new* value. It is equivalent to `j = j + 1`. Contrast with `j++`, which returns the old value.
