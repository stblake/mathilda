---
source: src/core.c
---
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
