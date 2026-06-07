---
source: src/core.c
---
A read-only system constant bound as an OwnValue in `system_constants_init` (`src/core.c`) via `register_system_constant`, then marked `ATTR_PROTECTED`. Its value is `expr_new_real(DBL_EPSILON)` — the `<float.h>` machine epsilon of the local IEEE 754 `double`.
