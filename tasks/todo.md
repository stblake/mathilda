# Task: a buffer fast path for Commonest (2026-08-02)

`Commonest[Join[Table[1, {10^7}], Table[3, {10^4}], {7, 8}]]` cost 880 ms
against 21.5 ms for `Tally` of the identical (already packed) argument. The
transparency gate names the cause: `Commonest` is not on `pack.c`'s AWARE list,
so every call materialised 10^7 boxed `Expr` before the builtin could start.
Commonest IS a tally plus a top-n selection, so the whole gap is the boxing.

## Plan

- [x] Measure: confirm the gate materialises, and that Tally over the same
      buffer is ~40x cheaper.
- [x] Extract the counting substrate `ndred_tally` already has (direct-index for
      a bounded int range, inline-key hash otherwise) into one helper both heads
      call, so their first-appearance ORDER agrees by construction rather than
      by two implementations happening to match.
- [x] `ndred_commonest`: count words, sort distinct by (count desc, first
      appearance asc), take n, restore first-appearance order, return the
      selected elements PACKED at the input dtype (uniform -> one head).
- [x] Handle every shape the head can now receive: `Commonest[a, n]`,
      `Commonest[a, UpTo[n]]`, rank >= 2, complex/float32 dtypes, a non-finite
      float64, an association -- each either handled or
      `ndarray_delist_and_reeval`'d. Opting into AWARE is that commitment.
- [x] Fix the `-1` sentinel collision in the List path: `Commonest[list, -1]`
      took the "no count given" branch and answered with the most common
      elements, where every other negative count gives `{}`.
- [x] `pack.c`: AWARE + INT64_OK.
- [x] Differential tests (plain vs packed vs visible NDArray) in
      `test_packed_list.c`; Commonest cases in `test_list.c`; a `commonest`
      probe in `tools/numeric_sweep.py` so the three-surface audit covers it.
- [x] `make check-packed-aware`, `check-array-exactness`, `check-nd-surfaces`,
      `check-c99`.
- [x] Docs: `docs/spec/builtins/structural-manipulation.md`,
      `docs/spec/builtins/packed-arrays.md`, this week's changelog
      (`docs/spec/changelog/2026-07-27.md`).

## Review

`AbsoluteTiming`, minimum of five, on `v = Join[Table[1, {10^7}],
Table[3, {10^4}], {7, 8}]` (an int64 packed list):

| call | before | after | speedup |
|------|-------:|------:|--------:|
| `Commonest[v]`    | 879.7 ms | 21.8 ms | 40.4x |
| `Commonest[v, 2]` | 877.2 ms | 22.1 ms | 39.6x |
| `Tally[v]`        |  21.5 ms | 22.5 ms | 1.0x  |

Commonest now costs what tallying costs -- the first and third rows are the same
measurement. The Tally row is the control on the shared-routine refactor:
unchanged within noise.

`tools/nd_surface_audit.py --only tally,commonest` reports 40.5x with all three
surfaces (plain List, packed List, visible NDArray) agreeing on the answer and
1.09x skew.

### Two defects fixed alongside, both in the existing List path

- `n == -1` was the sentinel for "no count given", so `Commonest[list, -1]`
  answered `{1}` where `Commonest[list, -2]` answered `{}`. Now an explicit
  flag; both are `{}`. This is a real behaviour change on the plain path.
- `builtin_commonest` passed its result array to `expr_new_function`, which
  COPIES rather than adopts, and never freed it -- one leaked array per call.
  Found by valgrind on the new tests; the leak summary drops to the known macOS
  baseline noise with no stack in the new code.

### Things worth knowing

- The counting helper is shared with Tally, so a change to insertion order
  cannot desync their tie-breaking. Commonest breaks a count tie by first
  appearance, and no gate looks at tie order -- only sharing the routine makes
  the two agree by construction.
- Adding a head to AWARE is a commitment for EVERY shape it can now receive, not
  only the ones the fast path handles (the trap `DeleteDuplicates[list, test]`
  fell into on 2026-08-01). Each decline here is explicit.
- `tests/test_packed_list.c` has two pre-existing soft-assert failures
  (`NDArrayQ[Range[249]]`, `{NDArrayQ[fb[5]], fb[5]}`) that also fail at HEAD --
  verified by rebuilding the suite from the stashed tree. Not from this change.

---

Previous round's record is in git history (`tasks/todo.md` at 749ba6e) and in
`docs/design/performance.md` §15.
