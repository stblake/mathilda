# Segfaults to fix

Running list of reproducible crashes found in Mathilda, with root cause and a
suggested fix. Each entry is independent. Crashes here are **pre-existing** (present
on the committed `main` at the noted commit) unless stated otherwise.

---

## 1. `Eigenvalues` of an exact (integer/symbolic) matrix ≥ 3×3 — packed-array deref in Faddeev char-poly

- **Status:** ✅ FIXED 2026-08-13. `eigen_mat_mul` (via `dot2`) returned a packed
  NDArray for the machine matrix product; `eigen_char_poly_faddeev`'s trace/shift
  helpers index it as a nested List. Fix: `eigen_mat_mul` now `pack_unpack`s its
  result, and `eigen_char_poly_faddeev` `pack_unpack`s the initial `M = A` (a
  large integer matrix packs to int64). `eigen_tests`, `mateigen_direct_tests`,
  `lapack_builtin_tests`, and `singularvaluedecomposition_tests` now pass.
  Original report retained below.
- ~~**Status:** open. Pre-existing on `main` @ `92453b9`~~ (confirmed by reproducing
  with all local changes stashed).
- **Severity:** high — segfault (SIGSEGV, `EXC_BAD_ACCESS`) on a common operation.
- **Repro:**
  ```
  Eigenvalues[{{2, -1, 0}, {-1, 2, -1}, {0, -1, 2}}]
  ```
  (any exact matrix `n ≥ 3` that reaches the symbolic characteristic-polynomial
  path; `N[...]` inexact matrices are unaffected — they take the machine LAPACK
  path and never call this code.)
- **Crash site:** `eigen_mat_trace` (`src/linalg/eigen_common.c:154`), at
  `M->data.function.args[i]->data.function.args[i]`.
- **Backtrace:**
  ```
  eigen_mat_trace           eigen_common.c:154
  eigen_char_poly_faddeev   eigen_common.c:234   (k = 2 iteration)
  eigen_compute_eigenvalues_full  eigen.c:104
  builtin_eigenvalues       eigen.c:219 / 142
  ```
- **Root cause (confirmed under lldb):** the Faddeev–Leverrier loop computes
  `M_k = A . (M_{k-1} − p_{k-1} I)` via `eigen_mat_mul` → `dot2` → `evaluate`.
  For an exact matrix the product now comes back as a **packed `EXPR_NDARRAY`**
  (recent automatic-packing work), not the `List`-of-`List`s that
  `eigen_mat_trace` and `eigen_mat_minus_scalar_id` assume. The first trace call
  (`M_1 = expr_copy(A)`, still `EXPR_FUNCTION`) succeeds; the second call — on the
  packed `M_2` — indexes `->data.function.args` on an `EXPR_NDARRAY` union member
  and dereferences garbage (`address=0x23`).
- **Suggested fix:** make the symbolic Faddeev path packing-transparent. The
  smallest, safest change is to unpack the product inside `eigen_mat_mul` before
  returning — e.g. route the `dot2` result through the delist/unpack helper the
  other symbolic eigen paths already use (`pack_delist` / `linalg_delist_and_reeval`
  family) so `M` is always a `List` of `List`s. Alternatively, teach
  `eigen_mat_trace` / `eigen_mat_minus_scalar_id` to read an `EXPR_NDARRAY` via the
  NDArray element accessor. Prefer the former (one spot, matches the assumption the
  rest of the file already encodes).
- **Blast radius note:** this aborts every unit suite whose fixtures evaluate
  `Eigenvalues` of an exact matrix — confirmed `eigen_tests`,
  `mateigen_direct_tests` (both at `test_direct_symbolic_ignores_dispatch`),
  `lapack_builtin_tests` (`test_eigen`), and `singularvaluedecomposition_tests`
  (its eigen cross-check). All four crash at `eigen_mat_trace + 72`,
  `address=0x22/0x23`. The machine/LAPACK eigen paths and the generalized-pencil
  path are unaffected and separately verified. A full unit-suite sweep on
  2026-08-13 was otherwise clean: 407 suites pass; the only other non-passes are
  pre-existing slow-corpus timeouts (`crc_corpus`, `intrat_corpus`,
  `numeric_stress`) and one stale test expectation
  (`zero_test`: `Attributes[PossibleZeroQ]` asserts `{Listable, Protected}` but
  the value is correctly `{Protected}` — PossibleZeroQ is not Listable).
