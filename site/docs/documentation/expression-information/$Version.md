# $Version

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`$Version`**

gives a string describing the version of Mathilda, including the versions of the libraries it was built against.

## Examples (2)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= $VersionNumber
Out[1]= 0.092

In[2]:= $Version
Out[2]= "Mathilda 0.092 (GCC 16.1.0, GMP 6.3.0, MPFR 4.2.2, FLINT 3.6.0, ECM 7.0.7, Raylib 5.5, Accelerate, Readline)"
```

## Implementation notes

**Attributes:** `Protected`.

## References

**See also:** [$VersionNumber](../../expression-information/$VersionNumber/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)
