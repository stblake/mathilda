# Heads with a numeric fast path that `Compile[]` cannot lower

**48 of 236.** Recorded 2026-08-02 by `make check-compile-coverage`
(`tools/compile_coverage.py`), which joins the NDArray kernel registry and
`src/pack.c`'s `AWARE` list against the binary's own `CompileDiagnostics`.
The §1 group (`Tr`, `Det`, `MatrixRank`, `Norm`) closed 2026-08-03.

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

Two extension points already exist in `src/compile/compile.c`:

| table | shape | opcode |
|---|---|---|
| `ND_FNS` | array → array, plus trailing **integer** arguments | `A_NDFN` |
| `ND_REDS` | rank-1 array → **scalar** | `V_NDRED` |

Both delegate to the interpreter's own entry point, so a new row is a lowering
whose answer is bit-identical to the interpreted one by construction. Most of
what follows is a statement about why a head does not fit one of those two
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

## 2. Nearly free: array → array with a single array operand

Same story for `ND_FNS`: one array in, one array out, an entry point that takes
the whole call.

| head | note |
|---|---|
| `Inverse` | `ndla_inverse`; rank 2 → rank 2, which `rank_rule 2` already expresses (`Transpose` uses it) |
| `Normalize` | `ndla_normalize`; same rank in and out |
| `MatrixPower` | `src/linalg/matpow.c` reads the buffer directly; array + **integer** is exactly `nextra = 1` |
| `ReverseSort` | `reverse_top_level` + `ndstruct_sort` (`src/sort.c`); wants a single wrapper entry point |
| `ConjugateTranspose` | delegates to `Transpose` ∘ `Conjugate`, both already lowerable; wants a wrapper or a two-opcode lowering |
| `PseudoInverse` | one array, but `ndla_pseudoinverse_direct` takes `(a, tol_automatic, tol_value)` rather than the call — deliberately, so that declining cannot re-enter the builtin that dispatched to it (`plans/HPC_IMPROVEMENT_PLAN.md` 10.2). Needs a call-shaped adapter that returns `NULL` rather than re-evaluating. |

## 3. Two ARRAY operands — needs an array × array opcode

`A_NDFN` carries **one** array register plus trailing integers, so none of these
fit however the table is extended. This is the single largest blocker on the
list: eleven heads wait on it, though `MapThread` also wants §7's callback
lowering and `LeastSquares` an adapter besides.

| head | entry point / note |
|---|---|
| `Dot` | the most valuable single entry here — matrix × vector is the commonest numeric kernel there is |
| `LinearSolve` | `ndla_linearsolve` |
| `Cross` | `ndla_cross` |
| `LeastSquares` | `ndla_leastsquares_direct` — takes `(a, b, tol…)`, not the call, so it needs an adapter as well as the opcode |
| `Inner`, `Outer` | general contraction |
| `Join`, `Riffle` | structural, n-ary over arrays |
| `ListConvolve`, `ListCorrelate` | both engines already work on flat `double` arrays |
| `MapThread` | callback *and* two arrays |

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
