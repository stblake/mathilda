# Verbatim

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

**`Verbatim[expr] is a pattern object that matches expr taken literally: the pattern constructs inside expr (Blank, Pattern, ...) are not interpreted, so Verbatim[x_] matches only the literal expression x_.`**

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

Only the literal x_ is rewritten, not the symbol x

```mathematica
In[1]:= {x_, x} /. Verbatim[x_] -> matched
Out[1]= {matched, x}
```

The literal pattern expression matches itself

```mathematica
In[2]:= MatchQ[x_, Verbatim[x_]]
Out[2]= True
```

But an ordinary value does not

```mathematica
In[3]:= MatchQ[5, Verbatim[x_]]
Out[3]= False
```

Find a literal pattern sitting inside data

```mathematica
In[4]:= Count[{a + b, x_ + y_, 1 + 2}, Verbatim[x_ + y_]]
Out[4]= 1
```

## Implementation notes

**Definition.** `Verbatim[expr]` is a pattern object that matches `expr` taken
**literally**: the pattern constructs inside `expr` (`Blank`, `Pattern`, ...) are
*not* interpreted, so `Verbatim[x_]` matches only the literal expression `x_` (i.e.
`Pattern[x, Blank[]]`), not an arbitrary expression. It has no builtin and no
rewrite rule of its own — it is a directive to the matcher — and it is `Protected`.
The docstring lives centrally in `info.c`.

**Representation.** `Verbatim` is recognised structurally inside `match_internal`
(`src/match.c`): when the pattern node is a one-argument function whose head is the
interned `SYM_Verbatim`, the matcher succeeds iff the subject is **structurally
equal** to the wrapped argument — `expr_eq(expr, Verbatim_arg)` — and then threads
the parent continuation; otherwise it fails. Because the test is `expr_eq`, no
variable is bound and no sub-pattern is interpreted. `Verbatim` sits alongside the
other matcher-transparent heads (`HoldPattern`, `Longest`/`Shortest`), and the
matcher's "is this a pattern object?" guard knows `SYM_Verbatim` so a literal
pattern expression is still reachable as subject matter.

**Usage & limits.** Use it wherever you must match or replace an expression that
*is itself* a pattern — searching a list for a literal `x_ + y_`, or rewriting the
symbol `_` as data rather than as `Blank[]`. Being an exact structural equality, it
binds nothing: `Verbatim[p]` never captures, so there is nothing to reuse on the
right-hand side of a rule beyond the literal match. It is a single-argument form.

**Attributes:** `Protected`.

## References

- Source: [`src/match.c`](https://github.com/stblake/mathilda/blob/main/src/match.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`Verbatim[expr]` matches `expr` taken literally — the pattern constructs inside it
(`Blank`, `Pattern`, ...) are not interpreted. So `Verbatim[x_]` matches only the
expression `x_` itself, which is what lets you search for, or rewrite, a pattern
*as data*. The matcher implements it as a structural-equality test (`expr_eq`), so
it binds no variables; it is transparent in the same way as `HoldPattern`, but
where `HoldPattern` keeps its argument interpretable as a pattern, `Verbatim`
freezes it to a literal.
