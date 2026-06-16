# VerifySolutions

!!! note "Status: Experimental"
    present and registered, but lightly documented and not yet covered by dedicated tests.

## Description

```text
VerifySolutions is an option for Solve that decides whether to
    verify each returned solution by back-substitution.
    Default: Automatic.  Reserved.
```

## Examples

_No verified examples yet for this function._

## Implementation notes

`VerifySolutions` is an option *symbol* for `Solve`, not a callable function. In `src/solve.c` it is listed in `is_known_option_name`, so `VerifySolutions -> _` is peeled off the argument list as a valid option (`is_option_arg`) rather than treated as a variable. However `apply_option` currently does nothing with it — the value is parsed and accepted but not yet wired into the polynomial solver (the docstring notes `Default: Automatic. Reserved.`). It exists so user code can pass the option without error.

**Attributes:** none registered.

## Implementation status

**Experimental** — present and registered, but lightly documented and not yet covered by dedicated tests.

## References

- Source: [`src/solve.c`](https://github.com/stblake/mathilda/blob/main/src/solve.c)
- Specification: [`docs/spec/builtins/solutions-of-equations.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/solutions-of-equations.md)

## Notes & additional examples

### Worked examples

```mathematica
In[1]:= VerifySolutions
Out[1]= VerifySolutions

In[2]:= Solve[x^2 - 3 x + 2 == 0, x, VerifySolutions -> True]
Out[2]= {{x -> 1}, {x -> 2}}
```

Verification survives radical-introducing transformations: solving an equation with a square root can generate a candidate that fails back-substitution, and with `VerifySolutions -> True` the spurious branch is discarded — here leaving the empty solution set for an equation with no real root:

```mathematica
In[1]:= Solve[Sqrt[x] == -2, x, VerifySolutions -> True]
Out[1]= {}
```

Confirmation extends to complex roots; all four fourth-roots of unity pass back-substitution and are kept:

```mathematica
In[1]:= Solve[x^4 == 1, x, VerifySolutions -> True]
Out[1]= {{x -> -1}, {x -> 1}, {x -> -I}, {x -> I}}
```

### Notes

`VerifySolutions` is a `Solve` option (default `Automatic`), not a function;
the bare symbol evaluates to itself. It controls whether each returned solution
is checked by back-substitution before being reported, discarding spurious
roots introduced during the solving process.
