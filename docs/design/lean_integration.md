# Lean 4 integration: design note

User reference: [`docs/spec/builtins/formal-verification.md`](../spec/builtins/formal-verification.md).
Demo: [`examples/lean/`](../../examples/lean/README.md). Lean-side package:
[`lean/`](../../lean/README.md).

## Trust model: untrusted oracle, trusted checker

Mathilda is an **untrusted oracle**. The Lean 4 kernel is the **trusted checker**.

- Mathilda computes an answer (an antiderivative, a factorization, an optimum) and a
  *certificate*: cofactors, a Gram matrix, a Pratt tree, KKT multipliers, a closed form for
  an induction. Mathilda's own algorithms (Risch, Gröbner, SDP, factoring) are never trusted.
  They only have to be right often enough to be useful.
- `Certify` writes a self-contained Lean 4 file that states the result as a theorem and
  proves it from the certificate. The Lean kernel decides whether it is accepted.
- A bug in Mathilda can make `Certify` fail to find a certificate, or produce a certificate
  that Lean rejects (`"Verified" -> False`). It cannot produce a wrong `True`, except
  through the steps listed in `"TrustedSteps"` below, or through a bug in Lean itself.
- In the other direction (`lean/`), the same holds: the tactics ask Mathilda for data and
  re-check every step in Lean. A wrong oracle answer makes the tactic fail
  (`lean/test/Untrusted.lean` drives the tactics with a fake oracle that lies).

### What `"Verified" -> True` means

All of the following hold, as checked by `` Lean`Internal`VerifySource `` and by the
per-theorem attribution in `Lean.m`:

1. Lean's exit code is 0.
2. No message with severity `error` falls in the theorem's line span, or outside every
   theorem (header or prelude).
3. No message has kind `hasSorry`.
4. The theorem's `#print axioms` report (appended by the generator) names only `propext`,
   `Classical.choice` and `Quot.sound`. `sorryAx` and the auxiliary axiom that
   `native_decide` adds both fail this test.

`Missing["NotChecked"]` means `Verify -> False` (the source was generated but not run).
`Missing["NotAvailable"]` means no Lean executable was found. A missing checker is not a
refutation, so it is never `False`.

### TrustedSteps: what Lean does not check, per kind

| Kind | Tier | `"TrustedSteps"` |
|---|---|---|
| Identity (polynomial) | 0 | none |
| Identity (rational function) | 0 | the `Together` normal form of each side (used to clear denominators) |
| Factorization | 0 | irreducibility of the factors |
| Antiderivative | 0 | the elementary derivative rules and quotient rule (applied symbolically); polynomial derivatives by `D`; the `Together` forms of the integrand and of each `u_i`; that the cleared denominators are not identically zero (checked by Mathilda) |
| Antiderivative | 1 (Mathlib) | none: `HasDerivAt (fun x => F) f x` on the stated domain |
| Bezout | 0 | none |
| IdealMembership (`Implies[...]`) | 0 | none |
| IdealMembership (`PolynomialReduce`) | 0 | that the remainder is fully reduced (the division identity is proved) |
| Roots (univariate, splits over ℚ) | 0 | none: completeness `p = 0 ↔ x = r1 ∨ …` is proved |
| Roots (systems, other polynomials) | 0 | completeness (only soundness is proved) |
| MatrixInverse, Determinant, LinearSolve | 0 | none |
| NullSpace | 0 | completeness (independence and span) |
| Nonnegativity, Optimality | 0 | the statement is proved over ℚ; the identity inside it is polynomial, so the bound holds over ℝ, but that transfer is not checked |
| Primality (Pratt) | 0 | Lucas's theorem (Lean evaluates the certificate checker, it does not prove that a valid certificate implies primality) |
| Primality (`PrimeQ[2]`) | 0 | a divisor of n is at most n |
| Primality (composite, `n < 2`) | 0 | none |
| ClosedFormSum | 0 | none |

### Tiers

- **Tier 0** uses core Lean 4 only (v4.35 or newer, for `grind`). No Mathlib, no Lake. Every
  generated file inlines `src/internal/lean/Prelude.lean` (about 80 lines: `powMod`,
  `prattCheck`, `PrattCertifies`, `matMul`, `matVec`, `idMat`, `det`, `sqnn`, `mulnn`).
  Proofs are `grind` (polynomial identities, with SOS hints), `decide +kernel` (matrices,
  Pratt) or induction (sums). A warm run takes about 0.5 s.
- **Tier 1** is opt-in and needs a Lake project with Mathlib built (`$LeanProject`). It is
  used for analytic statements that core Lean cannot express: `HasDerivAt` for
  antiderivatives with `Log`/`ArcTan`. `import Mathlib` costs 10–40 s per run.

## Architecture

```
Certify[comp] ──► Lean.m generator ──► Lean source ──► lean --json ──► verdict ──► ProofObject
   (HoldFirst)      (pure .m:            (LeanForm      (system(),       (NDJSON +
                     D, Together,         strict          temp files)     #print axioms)
                     FactorInteger,       translator)
                     SOS, LinearSolve)
```

| Piece | File | Role |
|---|---|---|
| `LeanForm` | `src/print_lean.c` | display form, and the strict translator `` Lean`Internal`ToLean `` that refuses anything whose Lean text would mean something else (inexact numbers, Int division, Mathlib-only heads in Tier 0, unknown heads) |
| runner | `src/lean_run.c` | `` Lean`Internal`RunProcess ``, `` Lean`Internal`VerifySource ``: temp files, `system()`, NDJSON parsing, axiom extraction |
| package | `src/internal/lean/Lean.m` | settings and discovery, one generator per certificate kind, batching, the session cache, `ProofObject`, `` Lean`Oracle `` |
| SOS engine | `src/solve/reduce_sos.c` | `` Reduce`SOSCertificate ``: the exact Gram blocks of `Reduce`'s Putinar stage (an additive wrapper; `Reduce` is unchanged) |
| Lean package | `lean/` | `mathilda_linear_combination`, `mathilda_sos`, `mathilda_prime`, `#mathilda` |

The generators are `.m` code because they only orchestrate existing builtins and assemble
strings, which follows the repo's "C for performance, rules for mathematics" split. The
package is loaded eagerly and defines rules only: it adds about 2.5 ms to startup.

**Context hygiene.** `` Lean` `` is on `$ContextPath` after load, so any symbol it exports
captures that short name in packages loaded later. The C builtins and private helpers
therefore live in `` Lean`Internal` `` and `` Lean`Private` ``, which are not on the path.
Only the public API and the option names are in `` Lean` ``.

## Why `system()`, and why it is safe

Strict C99 has no process API other than `system()`, and the repo forbids `popen`, `fork`
and `<unistd.h>` in new code. `system()` is enough because the runner only needs an exit
code and two output files.

- **No user text reaches the shell.** The command line contains only the executable path,
  flags written by Mathilda, and temp-file paths. The Lean source, the claim and any
  `Method` tactic string go into the `.lean` file.
- **Quoting is complete on POSIX.** Every argument is single-quoted, with each embedded `'`
  written as `'\''`. No character is special inside POSIX single quotes.
- **No planted files.** The command starts with `set -C` (noclobber), so the `>` and `2>`
  redirections fail if the output file already exists, including as a symlink. The source
  file name is probed before it is written, and all three files are removed afterwards.
- **Windows.** `cmd.exe` has no complete quoting, so an argument containing any of
  `" % ! ^ & | < >` is refused (`LeanVerify::badpath`). That branch compiles under MinGW but
  has not been exercised at run time.
- **Time bounds.** Each theorem carries `set_option maxHeartbeats (40000·TimeConstraint) in`,
  which is deterministic. `system()` blocks, and `TimeConstrained` measures this process's
  CPU time, so the timer does not fire during a Lean run and no child process is orphaned.
  A wall-clock kill needs a real spawn API (F7).

## Configuration

| Setting | Meaning |
|---|---|
| `$LeanExecutable` / option `LeanExecutable` | `Automatic` or a path. `Automatic` tries the environment variable `MATHILDA_LEAN`, then the newest `~/.elan/toolchains/*/bin/lean` (called directly, so no `lean-toolchain` file is needed), then `lean` on `PATH`. |
| `$LeanProject` / option `LeanProject` | `None` or a Lake project with Mathlib, for Tier 1. Defaults to the environment variable `MATHILDA_LEAN_PROJECT`. Tier 1 runs `cd <project> && <lake> env <lean> --json FILE`, with the `lake` beside the Lean binary. |
| `MATHILDA_BIN` | Lean side: the Mathilda binary that the `lean/` tactics call (else `set_option mathilda.path`, else `./Mathilda` or `../Mathilda`, else `PATH`). |
| `MATHILDA_LAKE` | test scripts: the `lake` to use for `lean/`. |

## Limitations

- Tier 0 statements are over `Int`/`Rat`. Real-analytic content (derivatives of `Log`,
  `ArcTan`) is either reduced to a cleared polynomial identity (Tier 0, with the trusted
  steps above) or proved in Mathlib (Tier 1).
- Antiderivatives must be built from rational functions, quadratic surds, and one level of
  `Log`/`ArcTan`/`ArcTanh`. `RootSum`, `Root` and complex results have no certificate (F9).
- `Det` is limited to n ≤ 6 (Laplace expansion in the kernel; 6×6 takes about 1.5 s).
- `Domain`, `Method` and `Assumptions` apply to identity claims. Other kinds take their
  hypotheses from the claim itself (`Implies[...]`, `ForAll[vars, cond, ...]`) and pick their
  own proof script.
- Factor irreducibility, Risch non-elementarity, `NullSpace` completeness and Lucas's
  theorem are not proved (listed as trusted steps).
- The Lean-side oracle starts one Mathilda process per tactic call (about 0.05 s).
- `example`s cannot be axiom-checked by `LeanVerify`; only `theorem`s can.

## Follow-ups

| ID | Item | Why deferred |
|---|---|---|
| F1 | `LinearProgramming` with an exact rational simplex | No simplex exists today. LP optimality is already certified for every LP that `Minimize`/`Maximize` solves. |
| F2 | Gosper certificates for hypergeometric sums | The induction step needs factorial/binomial facts from Mathlib. Polynomial·cᵏ terms do not. |
| F3 | Gröbner cofactor tracking (`<G> ⊆ <F>`) | A change to the hot Gröbner kernels. |
| F4 | Irreducibility certificates; Risch non-elementarity | Research items. |
| F5 | Tier 1 primality through Mathlib's `lucas_primality` | Removes the trusted Lucas step; needs `ZMod` proof assembly. |
| F6 | `DSolve` substitution certificates (Tier 1) | Needs a second-order `HasDerivAt` builder. |
| F7 | A spawn API with a wall-clock kill, and a public `RunProcess` | Portability work on three platforms. Heartbeats already bound each run deterministically. |
| F8 | A persistent NDJSON-pipe oracle for the Lean tactics | Only pays off for proofs with many tactic calls. |
| F9 | `Root`/`RootSum`/complex antiderivatives and radical `Solve` roots in Tier 0 | Needs an algebraic-number encoding beyond quadratic surds. |
| F10 | `NullSpace` completeness (rank witness) and eigenpairs | Needs a minor-based rank certificate. |
| F11 | A core-Lean field lemma from the raw rational statement to the cleared identity | Removes the "Together normal form" trusted step; core Lean lacks `div_add_div`. |
