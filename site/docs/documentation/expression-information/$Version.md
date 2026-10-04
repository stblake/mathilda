# $Version

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`$Version`**

gives a string describing the version of Mathilda, including the versions of the libraries it was built against.

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (2)

```mathematica
In[1]:= $VersionNumber
Out[1]= 0.266

In[2]:= $Version
Out[2]= "Mathilda 0.266 (GCC 16.1.0, GMP 6.3.0, MPFR 4.2.2, FLINT 3.6.0, ECM 7.0.7, Raylib 5.5, Accelerate, Readline)"
```

### Applications (2)

A descriptive string, not a number

```mathematica
In[3]:= StringQ[$Version]
Out[3]= True
```

Its head is String

```mathematica
In[4]:= Head[$Version]
Out[4]= String
```

## Implementation notes

**What it is / where defined.** `$Version` is a read-only system constant, bound as
an OwnValue in `system_constants_init` (`src/core.c`) via
`register_system_constant("$Version", expr_new_string(mathilda_version()))` and then
marked `ATTR_PROTECTED`. Its value is an `EXPR_STRING`.

**Its value.** The string is produced by `mathilda_version()` (`src/version.c`) and
assembled entirely at compile time: `"Mathilda " MATHILDA_VERSION_STRING " (<compiler>,
GMP ..., MPFR ..., FLINT ..., ...)"`. `MATHILDA_VERSION_STRING` lives in `src/version.h`
as the textual twin of `MATHILDA_VERSION_NUMBER` (the two are kept in sync by hand), and
the library versions come from each library's own preprocessor macros — so the string
names exactly what *this* binary was linked against. The pointer refers to static storage
and is never freed.

**Evaluation behaviour.** Being an OwnValue, `$Version` evaluates to its string in a
single step and, being Protected, cannot be reassigned. The exact text changes every
release, so a stable check queries it structurally (`StringQ[$Version]`) rather than
against a literal.

**Attributes:** `Protected`.

## References

**See also:** [$VersionNumber](../../expression-information/$VersionNumber/), [Real](../../other-advanced/Real/)

- The descriptive string is assembled in `src/version.c` from `MATHILDA_VERSION_STRING` (`src/version.h`) and each linked library's own version macros.
- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)

## Notes & additional examples

### Notes

`$Version` is a read-only Protected OwnValue holding a descriptive string, e.g.
`"Mathilda 0.266 (<compiler>, GMP ..., MPFR ..., FLINT ..., ...)"`, assembled at compile
time from `MATHILDA_VERSION_STRING` and each linked library's own version macros — so it
names exactly what this binary was built against. The numeric release is `$VersionNumber`.

The exact text changes every release, so the examples check it structurally rather than
against a literal.
