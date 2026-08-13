# Optional

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

patt:def or Optional\[patt, def\] is a pattern object that matches patt if it is present; if patt is omitted from the argument sequence, def is used in its place. patt\_. (sugar for Optional\[patt\_, Default\[f\]\]) draws the default value from Default\[f\] at the call site.

## Examples (6)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (6)

```mathematica
In[1]:= f[3]
Out[1]= 4

In[2]:= f[3, 10]
Out[2]= 13

In[3]:= g[1, 2]
Out[3]= {1, 2, 0}

In[4]:= lin[2]
Out[4]= 1 + 2 x^2

In[5]:= lin[2, 3, 4]
Out[5]= 4 + 3 x + 2 x^2

In[6]:= p[a]
Out[6]= {a, 0}
```

## Implementation notes

`Optional` is a pattern wrapper, not a function. The matcher in `match.c` peels `Optional[p]` / `Optional[p, default]` off a pattern (in the same wrapper-stripping loop that handles `Pattern`, `Shortest`, `Longest`), sets `is_optional`, and records the default. When the optional argument is absent at that position, the bound symbol is filled with the explicit `default` when given (`opt_container->args[1]`), otherwise with `get_default_value(pat_head, pos, total)` — which supplies the head's identity (0 for `Plus`, 1 for `Times`, etc., the head's `Default[]` value). When the argument *is* present it matches `p` normally. This is the mechanism behind the `x_.` / `x_:def` surface syntax. `Optional` is in the set of pattern heads `eval.c` leaves unevaluated.

**Attributes:** none registered.

## References

- Source: [`src/match.c`](https://github.com/stblake/mathilda/blob/main/src/match.c)
- Specification: [`docs/spec/builtins/pattern-matching.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/pattern-matching.md)
- Tests: [`tests/test_condition_downvalue.c`](https://github.com/stblake/mathilda/blob/main/tests/test_condition_downvalue.c)

## Notes & additional examples

### Notes

`Optional[p, def]` (surface syntax `p : def`) lets a pattern argument be omitted, supplying `def` in its place — the standard way to give function definitions default-valued parameters (Out[2], Out[5]). The shorthand `patt_.` is sugar for `Optional[patt_, Default[f]]`, drawing the default from `Default[f]` at the call site so a single rule can match expressions with or without a given term.
