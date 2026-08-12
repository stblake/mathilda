# Decrement

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Decrement[x] or x--`**

decreases the value of x by 1, returning the old value of x.

<details>
<summary>Notes</summary>

Decrement has attribute HoldFirst. In Decrement\[x\], x can be a symbol or a Part expression referring to an existing value (e.g. list\[\[2\]\]--). If x has no assigned value, Decrement::rvalue is emitted and the expression is left unevaluated.

</details>

## Examples (3)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (3)

```mathematica
In[1]:= k = 5
Out[1]= 5

In[2]:= k--
Out[2]= 5

In[3]:= k
Out[3]= 4
```

## Implementation notes

`builtin_decrement` (`src/core.c`) implements `x--` via the shared `increment_core` helper with a delta of `1`, negate=true, pre=false. `increment_core` requires the target be a symbol with an existing OwnValue (else `Decrement::rvalue`), evaluates the old value, forms and evaluates `Plus[old, Times[-1, 1]]`, and writes it back through an evaluated `Set`. Because pre=false it returns the *old* value (post-decrement). `Decrement` is `ATTR_HOLDFIRST`. The pre-form `--x` is the separate `builtin_predecrement` (pre=true).

**Attributes:** `HoldFirst`, `Protected`.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)
- Tests: [`tests/test_increment.c`](https://github.com/stblake/mathilda/blob/main/tests/test_increment.c)

## Notes & additional examples

### Notes

`k--` (`Decrement`) returns the *old* value of `k` and then decreases it by 1. Use `--k` (`PreDecrement`) to get the new value instead.
