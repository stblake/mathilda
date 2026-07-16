# $VersionNumber

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

```text
$VersionNumber
    gives the Mathilda version number as a real number.
```

## Examples

All examples below are verified against the current Mathilda build.

```mathematica
In[1]:= $VersionNumber
Out[1]= 0.01

In[2]:= $Version
Out[2]= "Mathilda 0.01 (16.1.0, GMP 6.3.0, MPFR 4.2.2, FLINT 3.6.0, ECM 7.0.7, Raylib 5.5, Accelerate, Readline)"
```

## Implementation notes

**Attributes:** `Protected`.

## Implementation status

**Stable** — documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## References

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
