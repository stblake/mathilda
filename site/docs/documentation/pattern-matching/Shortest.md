# Shortest

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Shortest[p] is a pattern object that matches the shortest sequence consistent with the pattern p.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`Shortest[p]` is a pattern-matching modifier resolved inside the argument matcher in `src/match.c`, with no builtin. When peeling a leading pattern, the matcher unwraps wrapper heads (`Pattern`, `Optional`, `Longest`, `Shortest`) and sets an `is_shortest` flag (a following `Longest` overrides it). The flag biases the backtracking sequence matcher to prefer the *fewest* arguments consistent with the rest of the pattern (e.g. for `Optional`+`Shortest` it tries the absent/default binding first), inverting the default greedy ("longest") preference. It does not change which matches are *possible*, only the order they are tried.

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/match.c`](https://github.com/stblake/mathilda/blob/main/src/match.c)
- Specification: [`docs/spec/builtins/pattern-matching.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/pattern-matching.md)

## Notes & additional examples

### Worked examples

With `/.` only the single first match is taken; `Shortest` makes that the
fewest-element partition, while `Longest` makes it the most:

```mathematica
In[1]:= {a, b, c} /. {Shortest[x__], y__} :> {{x}, {y}}
Out[1]= {{a}, {b, c}}

In[2]:= {a, b, c} /. {Longest[x__], y__} :> {{x}, {y}}
Out[2]= {{a, b}, {c}}
```

`ReplaceList` enumerates every match, revealing the order in which partitions are
attempted — `Shortest` tries the smallest first, `Longest` the largest:

```mathematica
In[1]:= ReplaceList[{a, b, c}, {Shortest[x__], y__} :> {{x}, {y}}]
Out[1]= {{{a}, {b, c}}, {{a, b}, {c}}}

In[2]:= ReplaceList[{a, b, c}, {Longest[x__], y__} :> {{x}, {y}}]
Out[2]= {{{a, b}, {c}}, {{a}, {b, c}}}
```

### Notes

`Shortest[p]` forces a sequence pattern to consume as few elements as possible, so the shortest partition is tried first (Out[1]). It is the dual of `Longest` (Out[2]).
