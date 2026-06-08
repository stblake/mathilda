# Optional

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
patt:def or Optional[patt, def]
    is a pattern object that matches patt if it is present; if patt is
    omitted from the argument sequence, def is used in its place.
patt_. (sugar for Optional[patt_, Default[f]]) draws the default value
from Default[f] at the call site.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`Optional` is a pattern wrapper, not a function. The matcher in `match.c` peels `Optional[p]` / `Optional[p, default]` off a pattern (in the same wrapper-stripping loop that handles `Pattern`, `Shortest`, `Longest`), sets `is_optional`, and records the default. When the optional argument is absent at that position, the bound symbol is filled with the explicit `default` when given (`opt_container->args[1]`), otherwise with `get_default_value(pat_head, pos, total)` — which supplies the head's identity (0 for `Plus`, 1 for `Times`, etc., the Mathematica `Default[]` value). When the argument *is* present it matches `p` normally. This is the mechanism behind the `x_.` / `x_:def` surface syntax. `Optional` is in the set of pattern heads `eval.c` leaves unevaluated.

**Attributes:** none registered.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/match.c`](https://github.com/stblake/mathilda/blob/main/src/match.c)
- Specification: [`docs/spec/builtins/pattern-matching.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/pattern-matching.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= f[x_, y_ : 1] := x + y

In[2]:= f[3]
Out[2]= 4

In[3]:= f[3, 10]
Out[3]= 13

In[4]:= g[x_, y_ : 0, z_ : 0] := {x, y, z}

In[5]:= g[1, 2]
Out[5]= {1, 2, 0}
```

### Notes

`Optional[p, def]` (surface syntax `p : def`) lets a pattern argument be omitted, supplying `def` in its place — the standard way to give function definitions default-valued parameters (Out[2], Out[5]).
