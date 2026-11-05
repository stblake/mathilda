# 124 — sparse linear algebra at scale

The native sparse solvers that landed in v0.6180–v0.6189 (`src/sparse/`), on
the problems they exist for, one size class above
[`120-sparse-arrays`](../120-sparse-arrays/): direct and Krylov `LinearSolve`
on 2-D/3-D Poisson, a nonsymmetric convection-diffusion operator and a complex
Helmholtz-like system; a `LinearSolveFunction` reused over 100 right-hand
sides; the action `MatrixExp[-L, u0]` of a heat step; the 6 largest and
smallest eigenvalues of a Laplacian plus a separating potential (ARPACK); and small exact integer matrices
against Maxima and sympy.

```
HPC_PYTHON=/usr/local/bin/python3.11 python3 benchmarks/run_all.py --only 124 --check-labels --timeout 600
```

## Halves

- `sparse_linalg.m`: Mathilda, and the reference column `run_all.py` times
  with the same file.
- `sparse_linalg.py`: scipy.sparse.linalg and sympy.
  - Direct rows are `spsolve` (SuperLU, COLAMD). scipy ships neither UMFPACK
    nor CHOLMOD, so SuperLU is the direct baseline for SPD, nonsymmetric and
    complex systems alike.
  - Krylov rows are `cg` / `gmres` (restart 50) / `bicgstab` at the same
    relative tolerance. scipy has no ILU0 or ILUT proper; its only incomplete
    factorisation is SuperLU's threshold ILU (`spilu`), used in natural order
    without pivoting (drop 1e-2 / fill 1 for ILU0, drop 1e-4 / fill 5 for
    ILUT). The preconditioner build is inside the timing, as it is in the
    `.m` half.
  - Many-RHS is `splu(A)` once and `.solve` per right-hand side; the heat
    step is `expm_multiply` (the same Al-Mohy–Higham algorithm); eigenvalues
    are ARPACK `eigsh` (`sigma=0` shift-invert for the smallest).
  - The exact rows use sympy's `DomainMatrix` over ZZ/QQ, the engine under
    `Matrix`.
- `sparse_linalg.mac`: Maxima. Every machine-precision row is SKIP (Maxima
  has no sparse numeric solver); the exact rows use `determinant`, `invert`
  and `charpoly` under `ratmx : true` (fraction-free elimination in CRE form —
  with the default `ratmx : false`, and with `newdet`, Maxima expands by
  minors on matrices this sparse and does not finish 20 × 20 in minutes),
  `linsolve_by_lu` and `nullspace` from `linearalgebra`, and core `rank`.

## Inputs and checks

Operators are built from `Band` and `KroneckerProduct` (`diags_array` and
`kron` on the Python side): the 5- and 7-point Laplacians, a convection-
diffusion operator `L + 0.4 (I ⊗ S) + 0.2 (S ⊗ I)` with `S` the central
first-difference stencil, and `L + (-0.3 + 0.05 i) I`. Right-hand sides are
`(k mod 7) - 3`.

The exact matrices come from one Lehmer step
`s = (48271 (n (i - 1) + j) + seed) mod (2^31 - 1)`: entry `(s mod 19) - 9`
where `s mod 3 = 0`, else 0, plus 10 on the diagonal. A rank-r matrix is the
product of three of them, `n × n`, `n × r`, `r × n`.

- **Numeric solves** check `Round[10 Norm[x]]` or `Round[Norm[x]]`, the scale
  chosen per case so the value sits well clear of a rounding boundary.
  `Total[x]` cancels (to about -57 over 90000 entries on the 300² grid), so a
  1e-8 Krylov residual moves it; the norm it does not.
- **Heat step**: `Round[1000 Norm[y]]`. **Eigenvalues**: `Round[10^6 Total]`;
  `Eigensystem` adds the sum of |v|², 6000 for six unit vectors. The operator
  is `L + V` with `V(i, j) = 0.001 (i² + 2j²)`, not the bare Laplacian: the
  Laplacian's extreme eigenvalues come in i ↔ j pairs about 3e-4 apart
  (relative), and at that gap the reference's default Arnoldi stops early —
  on the 300² grid of experiment 120 it returns 47985493 / 4684 for the
  six-largest / six-smallest digests, against the exact 47995643 / 4357 that
  scipy and the closed form agree on. The potential separates both ends of the
  spectrum.
- **Exact**: the determinant itself, the sum of the inverse's entries and of
  the solution (a Rational), the rank, the nullity, and the characteristic
  polynomial at x = 3 (n is even, so `Det[m - x I]` and `det(x I - m)` agree).

Sizes are held where the scipy column fits the file budget. SuperLU costs
~8× Mathilda's CHOLMOD on the 500² grid and ~100× on 3-D fill (124 s at 40³,
against 1.5 s), so the 3-D direct row is 30³ and the 10⁶-unknown row is the
3-D conjugate gradient.

## Results

Measured 2026-11-05, Mathilda v0.6189, i9-9880H, on a shared and loaded host
(load average 7-25 during the run), so read ratios under ~1.5x as parity. ms,
best of the timed runs; "—" is a SKIP.

| case | Mathilda | scipy / sympy | reference | Maxima |
|---|---|---|---|---|
| poisson2d 500^2: LinearSolve direct | 697 | 10 308 | 5 241 | — |
| poisson3d 30^3: LinearSolve direct | 469 | 20 525 | 2 721 | — |
| poisson3d 100^3 (10^6): CG 1e-6 | **1 803** | 789 | 1 547 | — |
| convdiff 300^2: LinearSolve Multifrontal | 690 | 1 466 | 2 241 | — |
| convdiff 300^2: GMRES + ILUT 1e-8 | 1 350 | 20 905 | 3 202 | — |
| convdiff 300^2: BiCGSTAB + ILU0 1e-8 | 1 413 | 7 926 | 4 852 | — |
| helmholtz 300^2 complex: LinearSolve direct | 2 346 | 3 985 | 4 837 | — |
| lsf 300^2: LinearSolveFunction, factor + 100 solves | 6 178 | 5 440 | 4 837 | — |
| heat 500^2: MatrixExp[-L, u0] (action) | 253 | 230 | 1 237 | — |
| exact 40x40: Det | 2.10 | 34.8 | 6.16 | 11.7 |
| exact 30x30: Inverse | 4.89 | 19.6 | 3.88 | 12.8 |
| exact 40x40: LinearSolve | 1.86 | 27.8 | 18.8 | 27.7 |
| exact 40x40: MatrixRank (rank 30) | 2.94 | 57.5 | 15.1 | 31.1 |
| exact 30x30: NullSpace (nullity 6) | 2.37 | 12.2 | 4.22 | 116 |
| exact 20x20: CharacteristicPolynomial | **126** | 17.4 | 26.0 | 4.30 |
| exact 10x10: Inverse x100 | 28.9 | 73.3 | 136 | 57.9 |
| eigen 200^2 L + V: Eigenvalues[H, 6] (largest) | unevaluated | — | — | — |
| eigen 200^2 L + V: Eigenvalues[H, -6] (smallest) | unevaluated | — | — | — |
| eigen 200^2 L + V: Eigensystem[H, -6] | unevaluated | — | — | — |

Every check agrees across the systems that answer. `run_all.py` classifies 14
AHEAD, 2 SLOWER, 3 CHECK-FAIL:

- **SLOWER, CharacteristicPolynomial 20 x 20 exact** (29x behind Maxima,
  7x behind sympy's `DomainMatrix.charpoly`). Maxima's `charpoly` under
  `ratmx` is one fraction-free determinant of `m - x I` in CRE form. 126 ms
  for a degree-20 integer characteristic polynomial, against ~2 ms for the
  40 x 40 determinant, suggests the polynomial-entry path rather than a
  division-free O(n^4) one (Berkowitz, or Hessenberg over Z mod primes).
- **SLOWER, unpreconditioned CG at 10^6 unknowns** (2.3x behind scipy). The
  iteration count is fixed by the operator, so this is per-iteration cost: one
  SpMV plus the vector updates. scipy's loop is BLAS-1 calls on contiguous
  buffers.
- **CHECK-FAIL, the three eigen rows**: `Eigenvalues[s, k]` and
  `Eigensystem[s, k]` on a 40000 x 40000 SparseArray say `SparseArray::dense`
  and return unevaluated (the densify cap is 2^27 elements). The sparse ARPACK
  path was being built when this experiment was added; the scipy and reference
  columns agree (747907106, 2571807, {2571807, 6000}), so these rows become
  timings as soon as it lands.

Where Mathilda leads, it is the factorisation: CHOLMOD / UMFPACK against
SuperLU (15x on the 500^2 grid, 44x on 3-D fill), and natural-order ILU
against SuperLU's threshold ILU, which is both slower to build and a weaker
preconditioner here (GMRES needs ~700 iterations with it at 1e-10 on the
300^2 Poisson matrix).
