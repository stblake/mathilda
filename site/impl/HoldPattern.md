---
source: src/attr.c
---
`HoldPattern` has no C handler; the attribute table in `src/attr.c` gives it `ATTR_HOLDALL | ATTR_PROTECTED` so its argument is not evaluated. In pattern matching it is transparent — `HoldPattern[p]` matches exactly as `p` does, letting `p` contain otherwise-evaluating constructs on a rule's left-hand side. `ReleaseHold` strips it.
