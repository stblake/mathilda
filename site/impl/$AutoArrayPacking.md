---
source: src/eval.c
---
**Definition.** `$AutoArrayPacking` is a **boolean system variable** controlling whether
Mathilda silently stores large lists of machine numbers as dense packed buffers rather than
boxed `Expr` lists. It is not a function: it is one of the entries in the `EVAL_SYSFLAGS[]`
table in `src/eval.c`, wired to `pack_set_enabled` / `pack_enabled`, with the docstring in
`src/info.c`. `eval_init_sysflags` registers it as an OwnValue holding the live C-flag state
so it can be read back; it carries no attributes (not even `Protected`).

**Representation.** `$AutoArrayPacking` is a bare `EXPR_SYMBOL` whose OwnValue is `True` or
`False`, mirroring the real C flag. Reading it returns the current state; assigning goes
through `eval_sync_sysflag`, which accepts only `True` or `False` — anything else is
rejected with a `$AutoArrayPacking::flagset` message (routed through `mth_message`) and the
OwnValue is rolled back to the live state so the symbol never lies about which path is
running. Registration happens *after* the `MATHILDA_NO_PACK` environment override is read,
so it reports `False` in a session started with `MATHILDA_NO_PACK=1`.

**Usage & limits.** Setting it changes storage and speed, not answers: a packed list is an
ordinary `List` with the same head, printed form, elements, ordering and pattern matches,
and only `NDArrayQ` tells the two apart. It does **not** affect `ToNDArray`/`ToPackedArray`
(explicit requests) or the visible `NDArray[...]` head. Turn it off (`$AutoArrayPacking =
False`) to force element-at-a-time list construction, e.g. for A/B timing or to isolate a
packing-related bug.
