---
source: src/core.c
references:
  - "The descriptive string is assembled in `src/version.c` from `MATHILDA_VERSION_STRING` (`src/version.h`) and each linked library's own version macros."
---
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
