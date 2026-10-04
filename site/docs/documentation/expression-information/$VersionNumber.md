# $VersionNumber

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`$VersionNumber`**

gives the Mathilda version number as a real number.

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

Its head is Real

```mathematica
In[3]:= Head[$VersionNumber]
Out[3]= Real
```

A positive release number

```mathematica
In[4]:= $VersionNumber > 0
Out[4]= True
```

## Implementation notes

**What it is / where defined.** `$VersionNumber` is a read-only system constant, bound as
an OwnValue in `system_constants_init` (`src/core.c`) as
`register_system_constant("$VersionNumber", expr_new_real(MATHILDA_VERSION_NUMBER))` and
then marked `ATTR_PROTECTED`. Its value is an `EXPR_REAL`.

**Its value.** `MATHILDA_VERSION_NUMBER` (`src/version.h`) is the single source of truth
for the release — a C `double` such as `0.266` — with `MATHILDA_VERSION_STRING` its
hand-synced textual twin feeding `$Version`. Every substantive commit bumps the number
(conventionally by `+0.001`) and tags the commit `v<MATHILDA_VERSION_STRING>`.

**Evaluation behaviour.** Evaluates to the Real in one step and is Protected, so it cannot
be reassigned. The printed form drops trailing zeros (so `0.160` prints as `0.16`), and
the value changes every release — a stable check therefore queries it structurally
(`Head[$VersionNumber] -> Real`) rather than against a literal.

**Attributes:** `Protected`.

## References

**See also:** [$Version](../../expression-information/$Version/), [Real](../../other-advanced/Real/)

- Source: [`src/core.c`](https://github.com/stblake/mathilda/blob/main/src/core.c)
- Specification: [`docs/spec/builtins/expression-information.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/expression-information.md)

## Notes & additional examples

### Notes

`$VersionNumber` is a read-only Protected OwnValue giving the Mathilda release as a Real,
such as `0.266`. It comes straight from `MATHILDA_VERSION_NUMBER` in `src/version.h`, the
single source of truth for the release; the descriptive form is `$Version`.

The numeric value changes every release (and the printer drops trailing zeros, so `0.160`
would print as `0.16`), so the examples test it structurally rather than against a literal.
