---
source: src/attr.c
---
`HoldAllComplete` is an attribute name, not a function. It maps to the bitflag `ATTR_HOLDALLCOMPLETE` (`attr_name_to_flag` / `get_attributes` in `src/attr.c`); when set on a symbol the evaluator holds all arguments and additionally bypasses Sequence flattening, `Unevaluated` stripping, and upvalue lookup for that head.
