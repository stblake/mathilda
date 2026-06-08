# Increment

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Increment[x] or x++
    increases the value of x by 1, returning the old value of x.

Increment has attribute HoldFirst. In Increment[x], x can be a symbol
or a Part expression referring to an existing value (e.g. list[[2]]++).
Increment threads over list values because Plus is Listable.
If x has no assigned value, Increment::rvalue is emitted and the
expression is left unevaluated.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`builtin_increment` (1-arg, the `x++` postfix form) builds the literal `1` and calls the shared `increment_core(lhs, 1, negate=false, pre=false, "Increment")`. `increment_core` first checks the target is an assignable symbol with an existing OwnValue (`lvalue_symbol_name` + `symtab_get_own_values`), emitting `Increment::rvalue` and returning NULL otherwise. It then `evaluate`s the current value, forms `Plus[old, ±dx]` and evaluates it to the new value, and writes the new value back by evaluating a `Set[lhs, new]` (Set has HoldFirst, so compound lvalues such as `Part[list, i]` are preserved for the assignment machinery). Because `pre` is false it returns the *old* value (post-increment semantics). The same helper backs `Decrement`, `PreIncrement`, `PreDecrement`, `AddTo`, and `SubtractFrom` via the `negate`/`pre`/`dx` flags.

**Attributes:** `HoldFirst`, `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= i = 5
Out[1]= 5

In[2]:= i++
Out[2]= 5

In[3]:= i
Out[3]= 6
```

### Notes

`i++` (`Increment`) returns the *old* value of `i` and then increases it by 1. Use `++i` (`PreIncrement`) to get the new value instead.
