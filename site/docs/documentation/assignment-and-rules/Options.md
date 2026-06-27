# Options

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
Options[f] gives the list of default option rules for the symbol f.
    Options[expr] gives the options explicitly set in an expression such
    as a graphics object.  Options[obj, name] gives the setting for the
    named option; Options[obj, {names}] gives a list of settings.  Assign
    to Options[f] to redefine all default options at once.
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
