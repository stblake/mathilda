# $AutoArrayPacking

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`$AutoArrayPacking`**

controls whether Mathilda stores large lists of machine numbers as dense buffers (packed arrays). True by default; set it to False to build every list one element at a time.

<details>
<summary>Notes</summary>

A packed list is an ordinary List -- same Head, printed form, elements, ordering and pattern matches -- and only NDArrayQ tells the two apart. So this changes storage and speed, not answers. Does not affect ToNDArray or ToPackedArray, which are explicit requests, nor the explicit NDArray\[...\] head. Reads back False in a session started with the environment variable MATHILDA\_NO\_PACK set. Only True or False is accepted.

</details>

## Examples (4)

Every input below was run against the current Mathilda build and its output recorded.

### Applications (4)

True by default -- automatic packing is on

```mathematica
In[1]:= $AutoArrayPacking
Out[1]= True
```

A big machine-number list packs, so this is True

```mathematica
In[2]:= NDArrayQ[Range[10000]]
Out[2]= True
```

With packing off, the same list is a plain List

```mathematica
In[3]:= $AutoArrayPacking = False; NDArrayQ[Range[10000]]
Out[3]= False
```

Restore the default

```mathematica
In[4]:= $AutoArrayPacking = True
Out[4]= True
```

## Implementation notes

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

**Attributes:** none registered.

## References

- Source: [`src/eval.c`](https://github.com/stblake/mathilda/blob/main/src/eval.c)
- Specification index: [`Mathilda_spec.md`](https://github.com/stblake/mathilda/blob/main/Mathilda_spec.md)

## Notes & additional examples

### Notes

`$AutoArrayPacking` is a boolean system variable: `True` (the default) lets Mathilda store
large lists of machine numbers as dense packed buffers; `False` builds every list one
element at a time. It is a read/write OwnValue tied to the real packing flag, and only
`True` or `False` is accepted — anything else draws a `$AutoArrayPacking::flagset` message
and is rolled back.

It changes storage and speed, not answers: a packed list is an ordinary `List` — same head,
printed form, elements, ordering and pattern matches — and only `NDArrayQ` distinguishes
the two, which is why toggling it flips the `NDArrayQ` results above. It does not affect the
explicit `ToNDArray`/`ToPackedArray` or the visible `NDArray[...]` head, and it reads back
`False` in a session started with `MATHILDA_NO_PACK` set.
