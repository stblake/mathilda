# SetOptions

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
SetOptions[s, name -> value, ...] sets default options for the symbol
    s and returns the new Options[s].  It can change Protected (but not
    Locked) symbols, and only changes existing options -- an unknown name
    raises SetOptions::optnf.  Use AppendTo[Options[s], ...] to add one.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= Options[LinearSolve]
Out[1]= {Method -> Automatic, Modulus -> 0, ZeroTest -> Automatic}

In[2]:= SetOptions[f, c -> 3]
Out[2]= SetOptions[f, c -> 3]
```

## Implementation notes

- `Options`, `SetOptions`, and `OptionValue` all have attribute `{Protected}`.
- Default options survive `Clear[f]` (only rules are cleared) and are removed

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/info.c`](https://github.com/stblake/mathilda/blob/main/src/info.c)
- Specification: [`docs/spec/builtins/assignment-and-rules.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/assignment-and-rules.md)
