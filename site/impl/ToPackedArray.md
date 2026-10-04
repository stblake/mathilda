---
source: src/pack.c
---
**Algorithm.** `ToPackedArray` is Mathematica's name (`Developer`​`ToPackedArray`) for
`ToNDArray` and is **the same C builtin** (`builtin_tondarray`) registered under a second name
— not a rule that rewrites to `ToNDArray`. Registering it as a genuine alias rather than a
DownValue means it costs no extra evaluation pass, never shows up in traces, and cannot be
shadowed by a user definition on the target. Every form and option is therefore identical to
`ToNDArray`: an optional trailing `DataType -> "..."` is stripped (`pack_take_dtype`), the one
positional argument is packed by `pack_force_coerce` → `pack_build` with no size threshold and
`coerce = true`, so a mixed machine `Integer`/`Real` list widens to a `float64` buffer while
an all-integer list stays `int64`. `ToPackedArray[list] === ToNDArray[list]` by construction.

**Data structures.** Identical to `ToNDArray`: an `EXPR_NDARRAY` with `present_as =
NDA_HEAD_LIST` (the packed-`List` surface — `Head` stays `List`), backed by a dense row-major
buffer sized `ndt_elem_size(dt) * n`.

**Complexity / limits.** `O(n)` (sniff + flatten) and one allocation; ignores the
automatic-packing threshold; rejects complex, ragged, empty, and non-machine input by
returning the list unchanged. See `ToNDArray` for the dtype-coercion rules.
