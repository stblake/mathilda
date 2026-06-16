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

A default-valued parameter: the second argument may be omitted, in which case
the default is supplied.

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

Several defaulted parameters give "optional trailing arguments" with sensible
fallbacks — here a quadratic whose linear and constant coefficients default to
`0` and `1`:

```mathematica
In[1]:= lin[a_, b_ : 0, c_ : 1] := a x^2 + b x + c

In[2]:= lin[2]
Out[2]= 1 + 2 x^2

In[3]:= lin[2, 3, 4]
Out[3]= 4 + 3 x + 2 x^2
```

The `_.` sugar draws the default from `Default[f]` at the call site, so a
pattern like `x_ + y_.` matches a bare term by treating the missing summand as
its additive identity `0` — the mechanism Mathematica uses to make rules robust
against absent structure:

```mathematica
In[1]:= p[x_ + y_.] := {x, y}

In[2]:= p[a]
Out[2]= {a, 0}
```

### Notes

`Optional[p, def]` (surface syntax `p : def`) lets a pattern argument be omitted, supplying `def` in its place — the standard way to give function definitions default-valued parameters (Out[2], Out[5]). The shorthand `patt_.` is sugar for `Optional[patt_, Default[f]]`, drawing the default from `Default[f]` at the call site so a single rule can match expressions with or without a given term.
