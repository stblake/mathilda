# Default

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`Default[f]`**

gives the default value supplied for a missing optional argument of f when the pattern \_. (Optional\[Blank\[\]\]) appears in a rule.

**`Default[f, i]`**

gives the default value for the i-th argument position of f.

## Examples (7)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (7)

```mathematica
In[1]:= Default[h] = 1
Out[1]= 1

In[2]:= h[x_, y_.] := {x, y}
Out[2]= Null

In[3]:= h[5]
Out[3]= {5, 1}

In[4]:= h[5, 9]
Out[4]= {5, 9}

In[5]:= Default[q] = 0
Out[5]= 0

In[6]:= q[a_, b_.] := {a, b}
Out[6]= Null

In[7]:= q[7]
Out[7]= {7, 0}
```

## Implementation notes

`Default` is not a builtin handler — it is a symbol whose DownValues the user assigns (`Default[f] = v`) and which the pattern matcher consumes when filling an Optional pattern (`x_.`, i.e. `Optional[x_, Default[f]]`). `get_default_value(pat_head, pos, total_pats)` in `src/default_helper.c` looks up the default for a function head: if no user `Default` DownValues exist it short-circuits to the built-in fallbacks (`Plus -> 0`, `Times`/`Power -> 1`, else none); otherwise it evaluates `Default[f, pos, total_pats]`, then `Default[f, pos]`, then `Default[f]`, returning the first that the user has defined (detected by the result differing from the constructed expression), before falling back to the same Plus/Times/Power defaults. The matcher calls this once per Optional-pattern attempt, hence the no-rule fast path.

**Attributes:** none registered.

## References

- Source: [`src/default_helper.c`](https://github.com/stblake/mathilda/blob/main/src/default_helper.c)
- Specification: [`docs/spec/builtins/pattern-matching.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/pattern-matching.md)

## Notes & additional examples

### Notes

`Default[f] = v` sets the value used for a missing optional argument matched by `_.` (`Optional[Blank[]]`). When the optional argument is supplied it is used directly, otherwise the stored default is substituted before the rule fires. This is the mechanism that lets `Plus`- and `Times`-style patterns absorb absent terms; for an ad hoc head you register the default yourself, as the two examples (additive default `0`, multiplicative-style default `1`) show. To attach a per-position default use the `Default[f, i]` form.
