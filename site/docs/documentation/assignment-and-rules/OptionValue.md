# OptionValue

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
OptionValue[name] gives the value of an option named name in the
    options matched by OptionsPattern[] in the enclosing rule.
    OptionValue[f, name] uses options associated with the head f;
    OptionValue[f, opts, name] resolves from the explicit rules opts then
    the defaults from f; a trailing Hold wraps the result in Hold.
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
