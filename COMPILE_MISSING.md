# Heads with a numeric fast path that `Compile[]` cannot lower

**35 of 238.** Recorded 2026-08-03 by `make check-compile-coverage`
(`tools/compile_coverage.py`), which joins the NDArray kernel registry and
`src/pack.c`'s `AWARE` list against the binary's own `CompileDiagnostics`.
Closed 2026-08-03: the §1 reductions (`Tr`, `Det`, `MatrixRank`, `Norm`), the §2
single-array heads (`Inverse`, `Normalize`, `MatrixPower`, `ReverseSort`,
`ConjugateTranspose`, `PseudoInverse`), and the §3 two-array heads (`Dot`,
`LinearSolve`, `Cross`, `LeastSquares`, `ListConvolve`, `ListCorrelate`, `Join`).

## Why this is a defect list and not a wish list

A head earns a numeric fast path by being something a numeric workload runs over
machine numbers. `Compile[]` and auto-compilation exist to run exactly those
workloads. A head that is fast at the REPL and unlowerable inside the
`Compile[]` meant to speed it up is a contradiction.

It is worse than "no speedup", because **the compilable subset is a cliff, not a
gradient**. One unlowerable head anywhere in a body sends the *whole* body to
the interpreter:

```mathematica
(* Mean did not compile, so Total did not either. Same answer, 10-40x slower,
   and nothing said so. *)
Compile[{{v, _Real, 1}}, Total[v] / Mean[v]]
```

That is why the list is worth closing head by head, and why
`check-compile-coverage` ratchets against `BASELINE` rather than merely
reporting: a head that *newly* falls out fails the build.

## How to read it

Every head below is on `AWARE` or has a registered kernel, so the transparency
gate leaves its buffer alone. That is *not* the same as the interpreter having
a buffer path — §9 is the two heads that are on `AWARE` and materialise anyway,
and they are the reason to check rather than assume.

The sections are the compiler-side mechanism that would close the group, and
heads in a group close together — which is the whole reason to group them
rather than list 48 names.

Four table-driven extension points exist in `src/compile/compile.c`:

| table | shape | opcode |
|---|---|---|
| `ND_FNS` | array → array, plus trailing **integer** arguments | `A_NDFN` |
| `ND_REDS` | rank-1 (or, per entry, rank-2) array → **scalar** | `V_NDRED` |
| `ND_FN2S` | two arrays → array | `A_NDFN2` |
| `ND_FN2S` | two arrays → **scalar** (`Dot`'s inner product) | `V_NDFN2` |

Each delegates to the interpreter's own entry point, so a new row is a lowering
whose answer is bit-identical to the interpreted one by construction. Most of
what follows is a statement about why a head does not fit one of those
shapes yet.

---

## 1. Closed 2026-08-03: the four scalar reductions now lower

`Tr`, `Det`, `MatrixRank`, and `Norm` — a rank-2 (`Norm`: rank-1 *or* rank-2)
array reduced to a **scalar** — now lower through `ND_REDS` / `V_NDRED`, each
delegating to its `ndla_*` entry point (`src/linalg/ndlinalg.h`). `nd_red_result`
gained a per-entry `rank` rule (rank 2 for the linalg reductions, `0` = either
for `Norm`) and an `int_result` flag for `MatrixRank`, whose rank is an Integer.

| head | entry point | result |
|---|---|---|
| `Det` | `ndla_det` | Real |
| `Tr` | `ndla_tr` | Real |
| `MatrixRank` | `ndla_matrixrank` | Integer |
| `Norm` | `ndla_norm` | Real (rank 1 *and* rank 2) |

A **real** matrix already routes through these same `ndla_*` functions in the
interpreter (all four are on `AWARE`), so the compiled answer is identical by
construction; an int or complex operand declines to the interpreter, which is
right to answer the exact number no machine slot holds.

## 2. Closed 2026-08-03: single array → array now lowers

`Inverse`, `Normalize`, `MatrixPower`, `ReverseSort`, `ConjugateTranspose`, and
`PseudoInverse` now lower through the existing `ND_FNS` / `A_NDFN` opcode, each
delegating to its call-shaped interpreter entry point (with a one-line adapter
where the entry takes direct args — `PseudoInverse` — or prints on a structural
mismatch a cheap pre-check can pre-empt — `Inverse`, `MatrixPower`).

| head | entry point | note |
|---|---|---|
| `Inverse` | `ndla_inverse` | `rank_rule 2`; real + complex (an int matrix is exact Rationals → declines) |
| `Normalize` | `ndla_normalize` | new `rank_rule 4` (require rank 1 → rank 1) |
| `MatrixPower` | `builtin_matrixpower` | array + **integer** is the `nextra = 1` slot |
| `ReverseSort` | `builtin_reverse_sort` | `rank_rule 0`, int + real (`Sort`'s buffer path declines complex) |
| `ConjugateTranspose` | `builtin_conjugate_transpose` | `rank_rule 2`, any dtype |
| `PseudoInverse` | `nd_pseudoinverse` → `ndla_pseudoinverse_direct(a, true, 0.0)` | real only |

`NdFnSpec` gained an `elems` operand-element gate (0 = any) — a head whose fast
path converts a dtype is barred from that operand type at compile time, since
`A_NDFN` carries no result-dtype field to catch a promise/result mismatch.

## 3. Closed 2026-08-03: two ARRAY operands now lower

`Dot`, `LinearSolve`, `Cross`, `LeastSquares`, `ListConvolve`, `ListCorrelate`,
and `Join` now lower through two new opcodes sharing one `ND_FN2S` table:
`A_NDFN2` (two arrays → array) and `V_NDFN2` (two arrays → scalar, for `Dot`'s
`vector·vector` inner product). Both are `K_ARR` and rebuild the whole call, like
`A_NDFN`; the optimiser treats them opaque by kind, so no optimiser change.

| head | entry point | result rank rule |
|---|---|---|
| `Dot` | `nd_dot2` → `nd_dot_machine` (BLAS-first, non-printing) | `ra + rb − 2`; `0` → scalar (`V_NDFN2`) |
| `LinearSolve` | `nd_linearsolve` → `ndla_linearsolve` | result rank = rhs rank |
| `Cross` | `nd_cross` → `ndla_cross` | rank 1 |
| `LeastSquares` | `nd_leastsquares` → `ndla_leastsquares_direct` | result rank = rhs rank; real only |
| `ListConvolve`, `ListCorrelate` | `builtin_list_convolve` / `_correlate` | rank 1 |
| `Join` | `ndstruct_join` | require equal rank; preserve dtype |

Still open here: **`Inner`, `Outer`** — their operands include function heads
(`Inner[f, a, b, g]`), so they want §7's callback lowering, not a plain
array × array opcode. **`Riffle`** is array + a scalar/list SEPARATOR, not two
arrays, so it belongs to §6 below. **`MapThread`** is a callback (§7).

## 4. Transforms — real in, complex out

| head | note |
|---|---|
| `Fourier`, `InverseFourier`, `FourierDCT`, `FourierDST` | `machine_path_ndarray` / `dct_machine_path_ndarray` in `src/fourier.c` are **`static`** and take `(nd, a, b, sign)` rather than the call. Needs a non-static entry point of the delegating shape; the result is a complex array from a real one, which `A_NDFN`'s "preserves the element type" rule does not currently express. |

## 5. Matrix producers — rank 2 out of a rank-1 argument

| head |
|---|
| `DiagonalMatrix`, `HankelMatrix`, `ToeplitzMatrix`, `VandermondeMatrix` |

`A_NDFN` delegates and *reuses* an operand's shape; these **produce** a shape.
Wants a producing opcode, closer to `A_NEW` than to `A_NDFN`.

## 6. Extra argument is not an integer

The blocker is the table, not the head: `ND_FNS` passes trailing operands as
`CT_INT` registers only.

| head | the argument that does not fit |
|---|---|
| `Riffle` | the separator — a scalar or a cycling list, not an array (so `A_NDFN2` does not fit it) |
| `PadLeft`, `PadRight` | the padding value (any element) |
| `Partition` | offset / a list spec |
| `Append`, `Prepend` | the element being added |
| `Catenate` | a list *of arrays* |
| `ExponentialMovingAverage` | a real α |
| `Extract` | a position spec, which may be a list of positions; `Extract[v, k]` with a plain integer is `Part` and could be routed there |
| `Rescale` | rewrites to `Rescale[list, {Min[list], Max[list]}]`, i.e. it needs the two-argument form's bound list |
| `PowerMod` | the **array** form only — `PowerMod[a, b, m]` over scalars already compiles (`OP_POWMOD_I`, `K_NARY`). `ndint_powermod` is the entry point; it is ternary, so the table would have to carry two trailing integers *and* an array. |

## 7. Callback lowerings

`Map`, `Select`, `TakeWhile`, `Fold`, `FoldList`, `Nest` and `NestList` already
compile their function argument. These are the ones nobody has written yet, and
each should follow the existing pattern rather than a new mechanism.

| head |
|---|
| `Scan`, `MapAll`, `MapAt`, `MapIndexed`, `SelectFirst`, `NestWhile`, `NestWhileList` |

## 8. The SCALAR form of the integer-only kernels

Their **array** form landed on 2026-08-02; only the scalar side is open, which
is what a loop like `Do[s += MoebiusMu[k], {k, n}]` needs.

| head | kernel |
|---|---|
| `MoebiusMu`, `EulerPhi` | unary, `to_int_i` only |
| `DivisorSigma` | binary, `to_int_i` only |

The obstacle is not the table but the opcode kinds: the `K_UN` / `K_BIN` scalar
kernel opcodes write `NaN` on failure, and these must **abort** —
`EulerPhi[-5]` has no machine answer, and the program has to hand back to the
interpreter rather than invent one. Wants an int → int kernel opcode that can
reach `vm_fail` the way `OP_FAIL` does, rather than one that writes a sentinel.

## 9. Not a compiler gap: no buffer path in the interpreter either

| head | site |
|---|---|
| `RowReduce` | `src/linalg/linsolve.c:1072` |
| `NullSpace` | `src/linalg/nullspace.c:286` |

Both open with

```c
if (linalg_call_has_ndarray(res)) return linalg_delist_and_reeval(res);
```

— they are on `AWARE`, so the transparency gate leaves the buffer alone, and
then they materialise it themselves on the first line and run the **exact
fraction-free** path. That is the `packed-aware but declines anyway` class
(`docs/design/performance.md` §14): being on `AWARE` is necessary, not
sufficient.

**There is nothing to lower here yet.** A compiled `RowReduce` would delegate
to a fast path that does not exist. The order of work is a LAPACK path first,
then a row in `ND_FNS`. Listed so the next reader does not spend an afternoon
looking for the entry point to delegate to.

## 10. One-off

| head | note |
|---|---|
| `LerchPhi` | an n-ary kernel is registered, but neither the scalar nor the array form lowers. The only head on the list that is missing *both*. |

---

## What is deliberately NOT here

52 further heads are **exempt with a recorded reason** in
`tools/compile_coverage.py`'s `EXEMPT` table. The reasons that recur:

- **The result length is data-dependent** — `Union`, `Intersection`,
  `Complement`, `DeleteDuplicates`, `Commonest`, `IntegerDigits`. A compiled
  register has a static type.
- **The result is not a uniform machine array** — `Tally` (ragged
  `{value, count}` pairs), `Counts` (an `Association`), `MinMax`, `Quartiles`
  and `QuotientRemainder` (tuples), and the `LUDecomposition` /
  `QRDecomposition` / `SingularValueDecomposition` triples.
- **The result type is not static** — `Eigenvalues` may be complex for a real
  matrix; `Eigenvectors`' shape depends on whether the matrix is defective.
- **There is no boolean array dtype** (§13 gap C.1) — `Positive`, `Negative`,
  `NonNegative`, `NonPositive`, and the definite-matrix predicates. Closing
  these needs the dtype, not a lowering.
- **The head asks about the REPRESENTATION**, which the compiler erases —
  `NDArrayQ`, `DataType`, `ByteCount`, `Head`, `AtomQ`, `Normal`,
  `ToNDArray`/`FromNDArray`.
- **It is not a callable head at all** — `If`, `Which`, `Module`, `With`,
  `Print`; each is lowered as control flow, register allocation or a side
  effect.

An entry leaving `EXEMPT` shows up as a gap by itself, which is why the table
carries reasons rather than just names.

## Regenerating this list

```bash
make check-compile-coverage            # gaps only; exits non-zero on a NEW one
python3 tools/compile_coverage.py --all
python3 tools/compile_coverage.py --only Det,Tr,Norm --all
```

The tool asks the binary, not the source: `CompileDiagnostics[argspec, body]`
returns `Compiled -> False` together with the innermost subexpression that could
not be lowered, over ~19 typed argument shapes per head. If a head here starts
compiling, the run says so and names the line to delete from `BASELINE`.

Two cautions learned building it, both of which produced false entries:

- **A probe list that cannot express a head's real call shape reports the
  probe's gap as the head's.** The first run named `Map`, `Select`, `Fold`,
  `Nest` and `GCD`, all of which compile — just not as `H[v]`. Functional
  spellings (`H[f, v]`) and integer-array shapes are in the table now.
- **`array` is one family across every rank and spelling**, on purpose. The
  question is whether a head has a machine lowering *at all*; reporting "no
  rank-1 `Inverse`" would be noise.

## See also

- [`docs/design/compile.md`](docs/design/compile.md) §8a — the delegation
  tables, the element-type rule, and why the two must not drift.
- [`docs/spec/builtins/packed-arrays.md`](docs/spec/builtins/packed-arrays.md)
  — which array heads compile today.
- `tools/nd_fastpath_sweep.py`'s `OFF_BUFFER` — the companion list, for heads
  that do not reach the buffer *in the interpreter* either.
