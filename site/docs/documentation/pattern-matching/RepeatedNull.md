# RepeatedNull

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`RepeatedNull[p, max] represents from 0 to max expressions matching p.`**

**`RepeatedNull[p, {min, max}] represents between min and max expressions matching p.`**

**`RepeatedNull[p, {n}] represents exactly n expressions matching p.`**

<details>
<summary>Notes</summary>

p... or RepeatedNull\[p\] is a pattern object that represents a sequence of zero or more expressions, each matching p.

</details>

## Examples (5)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (5)

```mathematica
In[1]:= MatchQ[{}, {RepeatedNull[1]}]
Out[1]= True

In[2]:= MatchQ[{1, 1}, {RepeatedNull[1]}]
Out[2]= True
```

This is exactly what makes a *trailing optional argument* expressible as a single
pattern: `f[x, RepeatedNull[y]]` matches both `f[x]` and `f[x, y, y]`:

```mathematica
In[1]:= MatchQ[f[x], f[x, RepeatedNull[y]]]
Out[1]= True

In[2]:= MatchQ[f[x, y, y], f[x, RepeatedNull[y]]]
Out[2]= True
```

A bound applies on top of the "zero allowed" base case, so `Cases` keeps the
empty and short lists alike:

```mathematica
In[1]:= Cases[{{a}, {a, a}, {a, a, a}}, {RepeatedNull[a, 2]}]
Out[1]= {{a}, {a, a}}
```

## Implementation notes

`RepeatedNull[p]` (`p...`) is the zero-or-more variant of `Repeated`, handled in the matcher. `is_repeated` in `src/match.c` recognises the head and sets the default length range to `[0, ∞)` (the only difference from `Repeated`'s `[1, ∞)`); the same optional count specs apply (`max`, `{n}`, `{min, max}` with `Infinity` permitted). The argument matcher backtracks over run lengths down to and including the empty run.

**Attributes:** none registered.

## References

- Source: [`src/match.c`](https://github.com/stblake/mathilda/blob/main/src/match.c)
- Specification: [`docs/spec/builtins/pattern-matching.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/pattern-matching.md)
- Tests: [`tests/test_match_extensive.c`](https://github.com/stblake/mathilda/blob/main/tests/test_match_extensive.c)
- Tests: [`tests/test_parse.c`](https://github.com/stblake/mathilda/blob/main/tests/test_parse.c)
- Tests: [`tests/test_replace.c`](https://github.com/stblake/mathilda/blob/main/tests/test_replace.c)

## Notes & additional examples

### Notes

`RepeatedNull[p]` (postfix `p...`) is like `Repeated` but matches *zero* or more occurrences, so an empty sequence also succeeds (Out[1]). Use it when the repeated part is optional.
