# MatchQ

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`MatchQ[expr, form]`**

gives True if expr matches the pattern form, False otherwise.

**`MatchQ[form]`**

is the operator form: MatchQ\[form\]\[expr\] == MatchQ\[expr, form\].

<details>
<summary>Notes</summary>

Pattern matching honours sequence variables (\_\_, \_\_\_), PatternTest, Condition, attribute-driven flattening / ordering, and the surrounding $Assumptions / DownValues environment.

</details>

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (7)

```mathematica
In[1]:= MatchQ[3, _Integer]
Out[1]= True

In[2]:= MatchQ[f[a, b], f[_, _]]
Out[2]= True

In[3]:= MatchQ[x^2, _^_]
Out[3]= True

In[4]:= MatchQ[{1, 2, 3}, {__Integer}]
Out[4]= True

In[5]:= MatchQ[7, _Integer?PrimeQ]
Out[5]= True

In[6]:= MatchQ[{2, 4, 6, 8}, {p__Integer} /; And @@ (EvenQ /@ {p})]
Out[6]= True

In[7]:= MatchQ[a + b + c, x_ + y_ /; x =!= y]
Out[7]= True
```

## Implementation notes

`builtin_matchq` (2-arg) is the user-facing entry into the pattern matcher. It `evaluate`s the first argument (the subject), takes the second argument as the pattern *without* evaluating it, allocates a fresh `MatchEnv` (the binding table), and calls `match(expr, pattern, env)` — the structural tree-unification engine in `match.c` that handles `Blank`/`BlankSequence`/`BlankNullSequence`, `Pattern` binding, `PatternTest`, `Condition`, `Alternatives`, `Optional`, `Repeated`, `Longest`/`Shortest`, and head-attribute-aware (`Flat`/`Orderless`/`OneIdentity`) sequence matching with backtracking. The boolean result becomes `True`/`False`; the env and the evaluated subject are freed. Bindings produced during the match are discarded — `MatchQ` reports only success/failure.

**Attributes:** `Protected`.

## References

- Source: [`src/match.c`](https://github.com/stblake/mathilda/blob/main/src/match.c)
- Specification: [`docs/spec/builtins/pattern-matching.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/pattern-matching.md)
- Tests: [`tests/test_association.c`](https://github.com/stblake/mathilda/blob/main/tests/test_association.c)
- Tests: [`tests/test_backtrack.c`](https://github.com/stblake/mathilda/blob/main/tests/test_backtrack.c)
- Tests: [`tests/test_eval.c`](https://github.com/stblake/mathilda/blob/main/tests/test_eval.c)
- Tests: [`tests/test_match.c`](https://github.com/stblake/mathilda/blob/main/tests/test_match.c)

## Notes & additional examples

### Notes

`MatchQ[expr, form]` tests whether `expr` matches the pattern `form`, returning
`True` or `False`. It supports the full pattern language: typed blanks
(`_Integer`), structural patterns (`f[_, _]`, `_^_`), and sequence variables
(`__Integer` for one-or-more integers, `___` for zero-or-more). Because `x^2` is
stored as `Power[x, 2]`, it matches the structural pattern `_^_`. The later
examples layer on `PatternTest` (`?PrimeQ`) and `Condition` (`/;`), so the match
succeeds only when the bound pieces also satisfy an arbitrary predicate — here
"every captured element is even" and "the two `Plus` operands are structurally
distinct". `MatchQ` is the predicate underlying filters like `Cases` and
conditional rules.
