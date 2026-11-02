# Automatic packed arrays

(Previous task — auto-compilation parity for the functional-programming heads —
archived to `plans/AUTOCOMPILE_FUNCPROG_TODO.md`.)

Store an ordinary `List` as a dense machine-precision buffer, invisibly: same
`Head`, printed form, elements, ordering and pattern matches, with only
`NDArrayQ` able to tell. Mathematica's packed arrays, over the storage Mathilda
already had for `NDArray[...]`.

Design: `docs/design/packed_arrays.md`. User surface:
`docs/spec/builtins/packed-arrays.md`.

The motivating measurement — on a 10^6-element list of Reals, `Length[x]` costs
30.6 ms against 2 µs packed, because the evaluator sweeps 10^6 argument nodes on
every pass. About 30 ns per element per pass, charged whenever the list is an
argument to anything.

---

## Phase 0 — prerequisites (no behaviour change) — DONE

- [x] `tests/bench_eval.c` calibration moved off `Total[Range[40000]]`, which
      packing makes ~20x faster and which would therefore have failed all six
      rows at once on the day of a large speedup. Now a list of `Rational`s,
      which a machine buffer can never hold. Baselines re-recorded; two rows had
      drifted badly *before* this change and are flagged in the file.
- [x] `Range` over an inexact iterator: deleted the per-element
      `evaluate(Plus[curr, di])` whose result that branch never reads.
      **340 ms → 93 ms (3.7x)** at n = 10^6. The exact branch is untouched.
- [x] `ndarray_part`'s scalar leaf routes through the shared element builder, so
      a fully-indexed `Part` of an integer array gives an exact `Integer`.
      Introduced `ndarray_buffer_element_to_expr` as the ONE place deciding what
      head an element materialises as.
- [x] Overflow-checked integer arithmetic hoisted from
      `src/compile/compile_internal.h` to `src/checked_int.h`, with `int64_t`
      spellings for the ndarray layer.

## Phase 1 — representation and the low-level surface — DONE

- [x] `NDPresentation present_as` on `NDArrayData`; `expr_unshare` carries it.
- [x] Constructor split: `expr_new_ndarray_raw` vs `expr_new_ndarray_like`, all
      48 sites visited (rename makes it compiler-enforced, not grep-enforced).
      Same split for the repack: `ndarray_from_nested_list_like`.
- [x] The identity trio — `expr_eq`, `expr_hash` (bit-identical to the
      materialised list's, so `Association` and `Union` work), `expr_compare`
      (List tier; buffer fast path only when shapes match exactly).
- [x] `Head` → `List`, `AtomQ` → `False`; `is_atomic`/`is_listq` deliberately
      NOT flipped (their call sites read `arg_count` immediately after).
- [x] Streaming `{...}` printing — standard, `FullForm`, `InputForm`, TeX.
- [x] `ToNDArray` / `FromNDArray`, `MATHILDA_NO_PACK`, `pack_set_enabled`,
      `pack_set_min_elements`; `DataType -> "int64"` now parses.
- [x] `tests/test_packed_list.c` — differential against the real materialised
      form, not hand-written expectations.

## Phase 2 — the transparency gate — DONE

- [x] `SymbolDef.packed_aware`, implied by the ndarray kernel setters.
- [x] The gate at `src/eval.c` step 2.7. Position is forced in both directions
      (after `flatten_sequences`, before `Listable`); see the design note.
- [x] The aware-head list in one reviewable block, `pack_mark_aware_heads()`.
- [x] Verified: `Count`, `Cases`, `Position`, `Level`, `ReplaceAll`, `Insert`,
      `Delete`, `Append`, `ListQ`, `MatchQ` and user `DownValues` with structural
      patterns are all correct on a packed list with no code of their own.
- [x] Verified packed-in/packed-out for 23 heads, and visible-in/visible-out for
      the explicit `NDArray[...]`.

### Found while doing this, not caused by it

- `MapIndexed` never repacked, on either surface — `MapIndexed[f, NDArray[…]]`
  has always returned a plain `List` while `Map` returns an `NDArray`. Pinned in
  the test file; fixing it belongs with `MapIndexed`.
- `Range[10^18, 10^18 + 3]` runs to the 10^6-element safety cap: the loop bound
  is a `double`, and one ulp at 10^18 is 128, so `val += 1` never advances.
  Phase 4's exact direct construction should fix it by looping in `int64`.
- `bench_eval`'s `replacerepeated fold` row is 2.4x slower than its recorded
  baseline and was sitting at 0.96 of the failure threshold. Only matcher-bound
  row in the set.

---

## Phase 3 — int64 exactness — DONE

Packing infers `int64` for exact lists, which put integer buffers in front of
code written when only `Compile[]` could make one. Two wrong answers were live
as soon as `ToNDArray` existed: `Total[{1,2,3}]` gave `6.` and `Sin[{1,2,3}]`
gave `{0,0,0}`.

- [x] `MATHILDA_PACK_DIAG` tripwire in the lossy accessors (`=1` warns, `=abort`
      aborts, `nd_int64_lossy_hit` is a breakpoint target). Audit by execution,
      not inspection.
- [x] `SymbolDef.packed_int64_ok` — the gate materialises an int64 packed
      argument for every head not on a short verified list, so every answer is
      the ordinary list's by construction.
- [x] The mixed-exactness write rule for `Part` assignment.
- [x] `tests/test_packed_list.c` runs 27 expressions packed-vs-plain, including
      the bigint-overflow cases.
- [x] **Exact int64 paths, which auto-packing turned from an optimisation into a
      correctness-of-performance item** — `Total[Range[10^6]]` measured 1.55x
      SLOWER packed before them, because the gate had to materialise first:
      - `nd_scalar` / `ndreduce.c`: `Total`, `Max`, `Min`, `Accumulate` exact
        Integer; `Mean`, `Median` exact reduced `Rational`.
      - `ndarray_elementwise`: exact int64 with `ci_*_i64` bailing on overflow,
        plus dtype promotion for a `Real` scalar (which had TRUNCATED:
        `Range[10] * 2.5` gave `{2, 5, 7, ...}`) and a bail for a
        Rational/BigInt/MPFR/Complex scalar.
      - `ndarray_scalar_power` / `_elementwise_power`: exact for a non-negative
        integer exponent, bail otherwise; `ndarray_base_scalar_power` bails.
      - `ndarray_dot2`: exact int64 accumulate.
      - `ndstruct.c`: int64 radix `Sort` (a double sort would reorder two
        integers past 2^53 AND round every element), `memcpy` `Transpose`,
        `Clip` degrades.
      - `Variance` / `Std` / `RMS` / `Quartiles` / moving statistics degrade on
        int64 — their exact answers are Rationals or radicals. This also fixed
        the same inexactness on the `NDArray[..., DataType -> "int64"]` surface.
      - `ndarray_int64_delist_retry` so `Plus`/`Times`/`Power` degrade to the
        List path instead of warning and leaving the call unevaluated.

## Phase 4 — producers — DONE (automatic packing is ON)

- [x] `ndbuild_open` in `Range` (both branches), `ConstantArray`, `RandomReal`,
      and `Table` gated on the compiled result type being `CT_REAL` (exposed as
      `autocompiled_result_is_real`), with `ndbuild_abandon` for the per-element
      interpreter fallback.
- [x] `pack_offer` for `Table`'s other branches, `Array`, `RandomInteger`,
      `Sort`, `Select`, `numloop`'s `reals_to_list` and `numloop_map`, and
      `ebuf_finalize` (`NestList` / `FoldList` / `NestWhileList` /
      `FixedPointList` — which also keeps the compiled and interpreted paths in
      step, so whether the body compiled does not decide the representation).
- [x] `pack_sniff` absorbs rows that are already buffers, so
      `Table[i j, {i,300}, {j,300}]` is one rank-2 array.
- [x] Threshold not observable: n = 249/250/251 agree on 13 expressions.
- [x] `Range[10^18, 10^18 + 3]` fixed (was 10^6 elements).

### Found by systematically sweeping every aware head and every kernel

Six wrong answers, each reachable from code that never mentions an array. The
sweep — packed vs `MATHILDA_NO_PACK=1` over 134 aware-head cases and 166 kernel
cases — is what found them; none was visible by reading the code.

- [x] **The evaluator's fixed-point test discarded the gate's work.** `expr_eq`
      is blind to packing by design, so a step whose only effect was
      materialising looked like no progress and the packed form was kept.
      `Table[i j, {i,300}, {j,300}]` came back as a List of 300 packed rows.
      Fixed with a gate tick counter bracketing each `evaluate_step`.
- [x] **A `Listable` head with both a plain List and a buffer broadcast instead
      of threading.** `NDArray[{1.,2.,3.}] * {10,20,30}` gave the 3x3 outer
      product. Pre-existing on the visible surface; the gate now materialises
      for a Listable head that also has a plain List argument.
- [x] **A repack coerced elements the operation introduced.**
      `Join[Range[1.,300.], {1}]` gave `1.` for the appended exact `1`.
      `pack_repack_like` re-sniffs the dtype for a packed source.
- [x] **An internal `evaluate()` bypasses the gate.** `Median` builds
      `Sort[data]`, evaluates it, then indexes the result's args — so
      `Median[Range[300]]` came back UNEVALUATED. Fixed with `pack_eval_plain`
      at the two `stats.c` sites, plus a blanket materialise in
      `internal_call_impl` (181 wrappers, ~721 call sites). The audit surface is
      tiny and greppable: only a head that can now produce a packed list matters.
- [x] **Six kernels answer with different element HEADS than the List does** —
      `Floor`, `Ceiling`, `Round`, `IntegerPart`, `Sign`, `Im` are real-closed,
      so they wrote `1.0` where the list gives the exact `1`. Plus `Clip` (clamps
      to its exact bounds) and `Precision`/`Accuracy` (Listable, so the list
      answers per element while the buffer path gave one scalar). New
      `symtab_clear_packed_aware`, and a `NOT_AWARE` list in `pack.c`.
- [x] **A resting expression could still hold a buffer.** An aware Listable head
      skips threading so its kernel can fire; when no kernel matches the arity
      (`Mod[buffer]`) the node came to rest un-threaded. Added the POST-GATE at
      `evaluate_step`'s final return.
- [x] `expr_to_latex` had no arm for a buffer, so every packed result came back
      with an empty `latex` field over the NDJSON pipe (as had every visible
      `NDArray[...]`). Now streams from the buffer.
- [x] `tests/bench_assoc.c` calibrated on `Total[Range[n]]`, which packing makes
      ~114x faster — all nine ops failed at once with nothing slower. Moved to a
      Rational list and re-recorded, exactly as Phase 0 did for `bench_eval`.
- [x] `tests/test_compile.c`'s `as_f64_buffer`: the harness fed the interpreter's
      result to `ndarray_from_nested_list`, which does not accept a value that is
      already a buffer.

### Measured (10^6 elements, packing on vs `MATHILDA_NO_PACK=1`)

| expression | packed | plain | |
|---|---|---|---|
| `Do[Length[x], {20}]` | 0.01 ms | 538 ms | 54000x |
| `Do[x[[7]], {20}]` | 0.01 ms | 558 ms | 56000x |
| `x . x` (int) | 0.95 ms | 448 ms | 471x |
| `Total[x]` (int) | 0.82 ms | 93 ms | 114x |
| `Range[10^6]` | 0.84 ms | 113 ms | 135x |
| `2 x` / `x + 1` (int) | 1.8 ms | 330 ms | ~185x |
| `Accumulate[x]` | 5.2 ms | 377 ms | 73x |
| `Total[x]` (real) | 1.4 ms | 92 ms | 67x |
| `Sin[x]` | 27.6 ms | 334 ms | 12x |
| `Sort[Reverse[x]]` (int) | 55 ms | 182 ms | 3.3x |

No row is slower packed. `bench_eval` all six within 5% of baseline;
`bench_compile` within gate.

### Recorded, NOT caused by this change

- `NDArray[{1.,2.,3.}] * {10,20,30}` still broadcasts on the **visible**
  surface. Same Listable bug; left alone because whether that form should answer
  with a `List` or an `NDArray` is a separate semantic decision.
- `MapIndexed` has never repacked, on either surface.
- `Precision` / `Accuracy` are `Listable` in Mathilda, so `Precision[{1,2,3}]`
  is a list of `Infinity` rather than a single value as in WL. Unrelated to
  packing; the packed list now matches the plain one either way.

## Flags and Mathematica-compatible names — DONE

- [x] `$AutoCompilation` (default `True`): switches off both invisible
      compilation mechanisms together — the `autocompile` adapter (Plot, Table,
      NIntegrate, NSum, FindRoot, NDSolve, the plot samplers) and `numloop`
      (Do/For/While/Map/Nest/Fold/FixedPoint bodies). `Compile[]` and any
      user-built `CompiledFunction` are deliberately unaffected.
- [x] `$AutoArrayPacking` (default `True`): switches off automatic packing at the
      producers. `ToNDArray` / `ToPackedArray` and the explicit `NDArray[...]`
      head are unaffected — those are the user asking.
- [x] Both interned in `sym_names.c`, documented in `info.c`, and registered as
      `OwnValue`s *after* the environment is read, so a session started with
      `MATHILDA_NO_PACK` / `MATHILDA_NO_AUTOCOMPILE` reports `False` rather than
      claiming `True` and being wrong.
- [x] Non-boolean assignment is refused with a `::flagset` message and the symbol
      rolled back to the live C state, following `$RecursionLimit::limset`.
      Implemented as a table (`EVAL_SYSFLAGS`) rather than a strcmp chain, behind
      a `'$'`-sigil pre-test so ordinary symbol assignment pays one byte compare.
- [x] `ToPackedArray` / `FromPackedArray` — the SAME builtin registered under
      Mathematica's names, not a rule that rewrites to the other (which would
      cost an evaluation pass, show in traces, and be shadowable). Both marked
      `packed_aware` + `packed_int64_ok` like their originals.
- [x] 5 new test functions (28 in `tests/test_packed_list.c`), including an
      `$AutoCompilation` differential over 11 expressions across both mechanisms
      and a flag-independence check.
- [x] **Regression found and fixed:** the hook first keyed on the `$` SIGIL and
      then built a probe, evaluating a delayed right-hand side. The REPL hooks
      (`$Pre`, `$PreRead`, `$Post`, `$PrePrint`, `$Epilog`) share that namespace
      with held right-hand sides, so all of them began evaluating at definition
      time. `repl_hooks_tests` caught it; the NAME is now checked before any probe
      is built, with a test asserting zero evaluations for a non-flag `$` symbol.
      A namespace prefix is not a membership test.
- [x] Leak-checked: definitely-lost is byte-for-byte identical to a `HEAD`
      baseline binary on the same corpus (13,376 B / 418 blocks, all macOS
      dyld/libobjc start-up; no definitely-lost record has a Mathilda frame).

### Note on `ToNDArray`'s return form

The request said "converts a list or nested lists to a `NDArray[...]` object",
but `ToNDArray` returns the **invisible packed form** — `Head` is `List`, it
prints as `{1., 2., 3.}`, and only `NDArrayQ` says `True`. That was settled
earlier in this work by explicit choice, and `ToPackedArray` being an alias
requires it: Mathematica's `Developer\`ToPackedArray` returns something that
prints and behaves as a `List`. `NDArray[list]` remains the way to build the
visibly different value. Flagged rather than silently reinterpreted.

## Phase 5 — the `Compile[]` boundary — DONE

Measured at n = 200000, a packed-list argument against the same value with
`MATHILDA_NO_PACK=1`, and against the visible `NDArray`:

| body | packed | plain List | visible NDArray | |
|---|---|---|---|---|
| `u^2 + 1.` x5 | 2.9 ms | 145 ms | 2.7 ms | **50x** |
| `Range[200000]` into a `_Real` slot x5 | 6.2 ms | 220 ms | — | **36x** |
| `Map[Sin[#] Exp[-#] + Sqrt[#] &, u]` x2 | 16.9 ms | 80.5 ms | 16.3 ms | **4.8x** |
| `Fold[Plus, 0., u]` x5 | 13.7 ms | 57.7 ms | 14.3 ms | **4.2x** |
| `Map[#^2 &, x]` at 10^6 | 37.9 ms | 213 ms | — | **5.6x** |
| `Map[Sin[#] Exp[-#] &, x]` at 10^6 | 32.7 ms | 180 ms | — | **5.5x** |

The packed column matches the visible-NDArray column — the boundary cost is gone,
not moved.

- [x] **The gate was the actual cost, not the marshalling.** A CompiledFunction's
      head is an `EXPR_COMPILED` with no `SymbolDef`, so it read as unaware and the
      gate materialised the very buffer the boundary exists to borrow.
      `f[Range[1., 200000.]]` was 75x slower than `f[NDArray[...]]`, two values
      differing only in `present_as`. An allowlist keyed on the symbol table cannot
      see a head that is not a symbol; `EXPR_COMPILED` needed the same explicit
      exemption the pure `Function` head already had.
- [x] `CfArgKind` splits the two questions one bool used to answer — who FREES the
      boxed value, and what the result must PRESENT as. Needed because a borrowed
      `EXPR_NDARRAY` no longer implies the visible head, and a dtype cast makes a
      temporary out of an argument whose presentation must still be honoured.
- [x] Boundary out: derived → inherit presentation **at any size** (the derived
      rule, not the producer rule); built → `pack_offer` (threshold and
      `$AutoArrayPacking`); complex → never packed, because `Complex[re, 0.]` does
      not round-trip. New `ndbuild_open_like` is the derived opener.
- [x] Boundary in: one O(n) `cf_cast_array` on a dtype mismatch instead of the
      VM's whole-call decline. Declines rather than rounds — past 2^53, narrowing
      to `_Integer`, or complex→real.
- [x] `cf_thread` (`RuntimeAttributes -> Listable`) repacks by re-SNIFFING the
      dtype via `pack_repack_like` instead of forcing `NDT_FLOAT64`, which would
      have turned an integer-valued body over a packed integer list into Reals.
      A visible source keeps `ndarray_from_nested_list_like`.
- [x] **C5, in the opposite direction from the plan.** A packed argument bypassed
      `numloop_map`'s compiled loop (it required `EXPR_FUNCTION`) and fell to the
      ndarray leading-axis walk, which runs the INTERPRETER per element. At 10^6,
      `Map[#^2 &, x]` was **424 ms** packed against 222 ms plain and
      `Map[Sin[#] Exp[-#] &, x]` **1120 ms** against 180 ms — automatic packing had
      made the most-used functional head up to 6x slower, and every value test
      passed. `numloop_map` now reads a rank-1 buffer, which also gives the O(1)
      type decision the plan wanted: the dtype answers for every element at once.
- [x] **Prerequisite fixed: a wrong compiled answer, older than packing.**
      `Compile[{{u, _Integer, 1}}, u * 2][{1, 2, 3}]` gave `{2., 4., 6.}`. An array
      opcode's scalar operand could only be Real or Complex, so the literal went
      into a real register and `ndarray_elementwise` then correctly widened the
      int64 buffer. New `AK_INT` operand kind (`ak_of`, `arr_prep`,
      `vm_scalar_pair`, `vm_box_scalar`, and the `A_PARTSET` right-hand side).
- [x] **Not done, deliberately:** the plot samplers still build ordinary Lists.
      Their output feeds renderer code that walks `data.function.args`, and
      measurement says there is nothing to win — `Plot` is at parity on and off
      (0.19 ms either way) because its cost is the compiled sampler.
- [x] 5 new test functions (33 in `tests/test_packed_list.c`); `compile_tests`
      clean; five differential corpora byte-identical on and off apart from
      `NDArrayQ`; valgrind unchanged from baseline.

## Phase 6 — differential testing — MOSTLY DONE

- [x] Table-driven sweep over every packed-aware head (134 cases) and every
      registered ndarray kernel (166 cases), packed vs `MATHILDA_NO_PACK=1`. Its
      failure list became Phase 4's fix list above; it now reports only the three
      intended observables: `NDArrayQ`, `DataType`, `ByteCount`.
- [x] Four development corpora (~110 expressions covering heads/exactness,
      mutation, aliasing, patterns, printing, identity, Association keys, nested
      producers, the repack family) agree byte for byte on and off.
- [x] `tests/test_packed_list.c`: `test_producers_pack`,
      `test_threshold_is_not_observable`,
      `test_exact_int64_arithmetic_on_a_buffer`,
      `test_regressions_automatic_packing_exposed`.
- [x] Whole suite (395 binaries) run against this tree and against a `HEAD`
      worktree baseline, and diffed.
- [ ] Table-driven sweep over the ~48 `src/list/*` builtins that are NOT aware —
      they are correct by materialisation, so this is confirmation, not a hunt.
- [ ] `tests/bench_pack.c` on the `bench_eval` normalized-baseline pattern (the
      numbers in Phase 4 were measured by hand and are not yet gated).

## Phase 7 — docs — DONE

- [x] `docs/spec/builtins/linear-algebra.md` — the "unlike Mathematica's packed
      arrays" paragraph now describes only the explicit `NDArray[...]`, and
      points at the packed form.
- [x] `docs/spec/builtins/packed-arrays.md` — "When packing happens" (the
      producer table and the threshold) and "Exact integers".
- [x] `docs/spec/builtins/lists-and-iteration.md` — `Table`, `Range`, `Array`,
      `ConstantArray`; `random-number-generation.md` — `RandomReal`,
      `RandomInteger`.
- [x] `docs/spec/changelog/2026-07-27.md` — "Automatic packing at the list
      producers".
- [ ] Complex packing, once the zero-imaginary fold is designed.
