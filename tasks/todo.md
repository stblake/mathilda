# Task: Implement `PackedArrayQ`

WL-faithful predicate: True only for a packed List (`is_packed_list`), False for
a visible `NDArray[...]`. Distinct from `NDArrayQ`.

## Steps
- [x] `src/ndarray.h` — declare `builtin_packedarrayq`
- [x] `src/ndarray.c` — implement `builtin_packedarrayq` + register in `ndarray_init` (ATTR_PROTECTED + docstring)
- [x] `src/pack.c` — add `"PackedArrayQ"` to `AWARE` (line 654) and `INT64_OK` (line 893)
- [x] `src/sym_names.h` — `extern const char* SYM_PackedArrayQ;`
- [x] `src/sym_names.c` — define + intern `SYM_PackedArrayQ`
- [x] `docs/spec/builtins/packed-arrays.md` — new `## PackedArrayQ` section
- [x] `docs/spec/changelog/2026-08-31.md` — changelog note
- [x] Build clean (`make -j`) — linked, no errors/warnings
- [x] REPL verification — all 11 cases correct; docstring retrieves via `Information[]`
- [x] check-c99 — clean
- [x] check-packed-aware — OK
- [x] check-array-exactness — 346 probes, 0 MIXED
- [x] check-nd-surfaces — exit 0, no disagreements

## Review

Implemented `PackedArrayQ[expr]` as a WL-faithful predicate distinct from
`NDArrayQ`: `True` only for the packed-`List` surface (`is_packed_list`),
`False` for a visible `NDArray[...]`.

Touched: `src/ndarray.{c,h}` (new `builtin_packedarrayq` + registration,
ATTR_PROTECTED + docstring), `src/pack.c` (added `"PackedArrayQ"` to `AWARE`
line 654 and `INT64_OK` line 893 — mandatory so the transparency gate doesn't
materialise a packed argument before the predicate runs), `src/sym_names.{c,h}`
(`SYM_PackedArrayQ`), `docs/spec/builtins/packed-arrays.md` (new section),
`docs/spec/changelog/2026-08-31.md` (changelog note).

Verified (all correct): plain List → False; packed float64/int64 → True;
visible `NDArray[...]` → False (vs `NDArrayQ` → True); scalars/strings → False;
wrong arity stays unevaluated; `Range[1000]` (auto-packed) → True;
`FromPackedArray[...]` round-trip → False. Docstring retrieves via
`Information[PackedArrayQ]`. Build clean; all four packed-array audits pass.

No `Compile[]`/kernel surface applies (boolean from node metadata, not a
buffer map/reduce). Single-arg form only; WL's `[expr, type]` / `[expr, type,
rank]` forms are out of scope.

## Follow-up: ToPackedArray mixed Integer/Real coercion

User reported `ToPackedArray[{1,2,3.}]` returned the list unchanged (should cast
to doubles + pack). Root cause: `pack_sniff` declines any mixed exact/inexact
list — correct for *automatic* packing (silently turning `1` into `1.` is an
observable head change, `1 === 1.` is `False`), wrong for the *explicit*
`ToNDArray`/`ToPackedArray`.

Fix (scoped to the explicit path only): threaded a `coerce` flag through
`pack_sniff`/`pack_build` (`src/pack.c`), added `pack_fold_class` (merges
Integer+Real → Real under `coerce`) and `pack_force_coerce`; `builtin_tondarray`
now calls the coercing variant. `pack_force` and all its other callers (listable
lift, part gather, linalg, dot, …) are unchanged, so `Sin[{1,2,3.}]` still stays
`{Sin[1], Sin[2], 0.1411…}`. Updated the `ToNDArray` docstring,
`docs/spec/builtins/packed-arrays.md`, and the changelog (also restored a
`## Interpolation` heading my earlier changelog edit had dropped).

Verified: `ToPackedArray[{1,2,3.}]` → `{1.,2.,3.}` float64, `PackedArrayQ` True;
all-int → int64; explicit int64 over reals → declines; nested mixed → float64;
`{True,1,2.}` → declines; `Sin[{1,2,3.}]` regression holds. Build clean;
check-c99 + check-packed-aware pass.
