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

`builtin_matchq` (2-arg) is the user-facing entry into the pattern matcher. It `evaluate`s the first argument (the subject), takes the second argument as the pattern *without* evaluating it, allocates a fresh `MatchEnv` (the binding table), and calls `match(expr, pattern, env)` — the structural tree-unification engine in `match.c` that handles `Blank`/`BlankSequence`/`BlankNullSequence`, `Pattern` binding, `PatternTest`, `Condition`, `Alternatives`, `Optional`, `Repeated`, `Longest`/`Shortest`, and head-attribute-aware (`Flat`/`Orderless`/`OneIdentity`) sequence matching with backtracking. The boolean result becomes `True`/`False`; the env and the evaluated subject are freed. Bindings produced during the match are discarded — `MatchQ` reports only success/failure.

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/match.c`](https://github.com/stblake/mathilda/blob/main/src/match.c)
- Specification: [`docs/spec/builtins/pattern-matching.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/pattern-matching.md)

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
