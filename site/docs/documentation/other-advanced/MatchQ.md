# MatchQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
MatchQ[expr, form]
    gives True if expr matches the pattern form, False otherwise.
MatchQ[form]
    is the operator form: MatchQ[form][expr] == MatchQ[expr, form].
Pattern matching honours sequence variables (__, ___), PatternTest,
Condition, attribute-driven flattening / ordering, and the surrounding
$Assumptions / DownValues environment.
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

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= MatchQ[3, _Integer]
Out[1]= True
```

```mathematica
In[1]:= MatchQ[f[a, b], f[_, _]]
Out[1]= True
```

```mathematica
In[1]:= MatchQ[x^2, _^_]
Out[1]= True
```

```mathematica
In[1]:= MatchQ[{1, 2, 3}, {__Integer}]
Out[1]= True
```

### Notes

`MatchQ[expr, form]` tests whether `expr` matches the pattern `form`, returning
`True` or `False`. It supports the full pattern language: typed blanks
(`_Integer`), structural patterns (`f[_, _]`, `_^_`), and sequence variables
(`__Integer` for one-or-more integers, `___` for zero-or-more). Because `x^2` is
stored as `Power[x, 2]`, it matches the structural pattern `_^_`. `MatchQ` is the
predicate underlying filters like `Cases` and conditional rules.
