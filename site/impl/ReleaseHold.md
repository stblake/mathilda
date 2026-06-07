---
source: src/core.c
---
`builtin_releasehold` (`src/core.c`) calls `release_hold_recursive`, which walks the tree and strips one layer off every `Hold`-family wrapper it finds (`Hold`/`HoldForm`/`HoldPattern`/`HoldComplete`, per `is_hold_head`), replacing a single-argument wrapper with its content and a multi-argument one with `Sequence[...]`. It does not recurse into the released contents.
