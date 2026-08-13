# $VersionNumber

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`$VersionNumber`**

gives the Mathilda version number as a real number.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= $VersionNumber
Out[1]= 0.044

In[2]:= $Version
Out[2]= "Mathilda 0.044 (Apple LLVM 21.0.0 (clang-2100.1.1.101), GMP 6.3.0, MPFR 4.2.2, FLINT 3.6.0, ECM 7.0.7, Raylib 5.5, Accelerate, Readline)"
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [$Version](../../expression-information/$Version/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
