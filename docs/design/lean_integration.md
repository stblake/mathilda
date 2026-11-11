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

The file is header (`set_option autoImplicit false`), prelude, theorem chunks, then a
**trusted footer** with, per theorem, `example : <statement> := mathilda_cert_k` and
`#print axioms mathilda_cert_k`. All of the following hold, as checked by
`` Lean`Internal`VerifySource `` and by the per-theorem attribution in `Lean.m`:

1. No message with severity `error` or kind `hasSorry` falls in the theorem's chunk or its
   footer lines.
2. No such message falls outside every chunk and footer (header or prelude).
3. The `#print axioms` report printed on the theorem's **own footer line** names only
   `propext`, `Classical.choice` and `Quot.sound` (a report on any other line is ignored).
   `sorryAx` and the auxiliary axiom that `native_decide` adds both fail this test.
4. Lean's exit code is 0, or every error lies in some theorem's lines (then the sibling
   theorems that own those errors are `False`, and the clean ones are re-run in a file of
   their own, because a file with errors is not compiled and so cannot be replayed).
5. **Kernel replay.** The file was compiled to an `.olean`, and a fresh `lean` process
   running `src/internal/lean/Replay.lean` (none of the checked file's code) replayed every
   declaration of the module into a new kernel environment and found only standard axioms
   under the footer theorem `mathilda_cert_k_footer : S := mathilda_cert_k`. Conditions 1-4
   read the checked run's own output; condition 5 does not, so it holds even if that output
   were forged.

The heuristic engine never reports `True`: its answer is `Missing["Heuristic", b]`.

`Missing["NotChecked"]` means `Verify -> False` (the source was generated but not run).
`Missing["NotAvailable"]` means no Lean executable was found. A missing checker is not a
refutation, so it is never `False`.

### TrustedSteps: what Lean does not check, per kind

| Kind | Tier | `"TrustedSteps"` |
|---|---|---|
| Identity (polynomial) | 0 | none |
| Identity (with `Assumptions`) | 0 | none: stated over every char-0 field (or every ordered field) |
| Identity (rational function) | 0 | field algebra from the cross-multiplied identity (unreduced numerators and denominators) to the rational identity where the denominators are nonzero |
| Factorization | 0 | irreducibility of the factors |
| Antiderivative | 0 | the elementary derivative rules and quotient rule (applied symbolically); polynomial derivatives by `D`; the `Together` forms of the integrand and of each `u_i`; that the cleared denominators are not identically zero (checked by Mathilda) |
| Antiderivative | 1 (Mathlib) | none: `HasDerivAt (fun x => F) f x` on the stated domain |
| Bezout | 0 | none |
| IdealMembership (`Implies[...]`) | 0 | none (every char-0 field) |
| IdealMembership (`PolynomialReduce`) | 0 | that the remainder is fully reduced (the division identity is proved) |
| Roots (univariate, splits over ℚ) | 0 | none: completeness `p = 0 ↔ x = r1 ∨ …` is proved |
| Roots (systems, constraints, other polynomials) | 0 | completeness (only soundness is proved: every constraint, and every denominator nonzero, at every solution) |
| MatrixInverse, Determinant, LinearSolve | 0 | none |
| NullSpace | 0 | completeness (independence and span) |
| Nonnegativity, Optimality | 0 | none: stated over every ordered field, ℝ included |
| Primality (Pratt) | 0 | Lucas's theorem (Lean evaluates the certificate checker, it does not prove that a valid certificate implies primality) |
| Primality (`PrimeQ[2]`) | 0 | none: `Mathilda.IsPrime 2` |
| Primality (composite, `n < 2`) | 0 | none: `¬ Mathilda.IsPrime n` |
| ClosedFormSum | 0 | none |

### Tiers

- **Tier 0** uses core Lean 4 only (v4.35 or newer, for `grind`). No Mathlib, no Lake. Every
  generated file inlines `src/internal/lean/Prelude.lean` (about 90 lines: `powMod`,
  `prattCheck`, `PrattCertifies`, `IsPrime`, `matMul`, `matVec`, `HasShape`, `HasLength`, `idMat`,
  `det`, `sqnn`, `mulnn`).
  Proofs are `grind` (polynomial identities, with SOS hints), `decide +kernel` (matrices,
  Pratt) or induction (sums). A warm run takes about 0.5 s, plus about 1.1 s for the kernel
  replay.
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
| runner | `src/lean_run.c` | `` Lean`Internal`RunProcess ``, `` Lean`Internal`VerifySource ``: temp files, `system()`, NDJSON parsing, axiom extraction, the replay run and its report |
| replay checker | `src/internal/lean/Replay.lean` | `lean --run` script: replays the compiled module in a fresh kernel environment, prints each declaration's axioms |
| package | `src/internal/lean/Lean.m` | settings and discovery, one generator per certificate kind, batching, the session cache, `ProofObject`, `` Lean`Oracle `` |
| SOS engine | `src/solve/reduce_sos.c` | `` Reduce`SOSCertificate ``: the exact Gram blocks of `Reduce`'s Putinar stage (an additive wrapper; `Reduce` is unchanged) |
| Lean package | `lean/` | `mathilda_linear_combination`, `mathilda_sos`, `mathilda_prime`, `#mathilda`; `Mathilda/Answer.lean` parses the oracle's answer as data |

The generators are `.m` code because they only orchestrate existing builtins and assemble
strings, which follows the repo's "C for performance, rules for mathematics" split. The
package is loaded eagerly and defines rules only: it adds about 2.5 ms to startup.

**Context hygiene.** `` Lean` `` is on `$ContextPath` after load, so any symbol it exports
captures that short name in packages loaded later. The C builtins and private helpers
therefore live in `` Lean`Internal` `` and `` Lean`Private` ``, which are not on the path.
Only the public API and the option names are in `` Lean` ``.

## Why `system()`, and why it is safe

Strict C99 has no process API other than `system()`, and the repo forbids `popen` and
`fork` in new code. `system()` is enough because the runner only needs an exit code and two
output files. The only POSIX calls (`mkdtemp`, `open`, `fdopen`, `rmdir`, under
`_POSIX_C_SOURCE 200809L` plus `_DARWIN_C_SOURCE` on macOS, guarded by `#ifndef _WIN32`)
create the private scratch directory.

- **No user text reaches the shell.** The command line contains only the executable path,
  flags written by Mathilda, and temp-file paths. The Lean source, the claim and any
  `Method` tactic string go into the `.lean` file.
- **Quoting is complete on POSIX.** Every argument is single-quoted, with each embedded `'`
  written as `'\''`. No character is special inside POSIX single quotes.
- **No planted files.** Every run gets a fresh directory from `mkdtemp` (mode `0700`), and
  the source is created in it with `O_CREAT | O_EXCL | O_NOFOLLOW` (mode `0600`). No other
  user can create, replace or symlink anything there, so the `.out` that is parsed is the
  one the shell's redirection created for this run (if it is missing, the run is a
  failure); `set -C` (noclobber) stays as a second guard. The files and the directory are
  removed afterwards. (The earlier scheme -- predictable names in the shared temp
  directory, probed with `fopen` and then opened -- could be raced: a dangling symlink at
  the source name was followed, and a pre-created `.out` was parsed as Lean's output.)
- **Windows.** `cmd.exe` has no complete quoting, so an argument containing any of
  `" % ! ^ & | < >` is refused (`LeanVerify::badpath`). That branch compiles under MinGW but
  has not been exercised at run time.
- **Time bounds.** Each theorem carries `set_option maxHeartbeats (40000·TimeConstraint) in`,
  which is deterministic. `system()` blocks, and `TimeConstrained` measures this process's
  CPU time, so the timer does not fire during a Lean run and no child process is orphaned.
  A wall-clock kill needs a real spawn API (F7).

## Soundness hardening (v0.6518)

An adversarial review produced four false `"Verified" -> True` results and several weaker
holes. Each is closed at its root, and each reproducer is a regression test
(`test_review_*` in `tests/test_lean_verify.c`, `test_hardening` in
`tests/test_lean_certify.c`).

| Hole | Fix |
|---|---|
| **Method injection.** `Method` was pasted as `"by\n  " <> m`; a second line at column 0 became a top-level command, and a `macro_rules` turned a later `theorem mathilda_cert_2` of a batch into a `def`, whose empty `#print axioms` report passed `x == x + 1`. | Three layers. (1) A Method script is refused (`Certify::method`) if it contains a control character, tab, `«`, `»`, a backquote, `$`, `#`, `@[`, `--`, `/-` or `"`, or uses any command keyword or meta escape as a word (`theorem`, `def`, `example`, `macro`, `syntax`, `elab`, `set_option`, `attribute`, `open`, `namespace`, `end`, `axiom`, `run_tac`, `native_decide`, `sorry`, `Lean`, `IO`, …); every line is indented under `by`. (2) Each item with a user script runs in a Lean file of its own. (3) For every theorem the trusted footer, after all chunks, re-elaborates the statement (`example : S := mathilda_cert_k`) and prints its axioms, and only the report on that footer line counts -- so a chunk cannot pass by declaring something else under the theorem's name, and a report elsewhere (forged or not) is ignored. A deny-list was chosen over a tactic allow-list because lemma names and term arguments are open-ended; the deny-list blocks every way to leave the tactic block or extend the parser, and layers 2-3 bound the damage of anything it misses. |
| **Vacuous Int hypotheses.** `Certify[x == 0, Assumptions -> 2 x == 1]` proved `∀ x : Int, 2 * x = 1 → x = 0`. | The automatic domain is `Int` only for hypothesis-free integer polynomial identities (which then hold in every commutative ring). With `Assumptions` it is `Rat`; ideal membership (`Implies[...]`) is always `Rat` (its cofactor identity holds in every commutative ring anyway). `Domain -> Integers` is still honoured and labelled. Variables of the assumptions are bound. |
| **Solve and `x / 0 = 0`.** Solve's spurious root `x = 1` of `(x-1)/(x^2-1) == (x-1)/(2x-2)` was "verified", because both sides are `0` in Lean. Only `Equal` constraints were checked (`x^2 == 1 && x > 0` checked `1^2 = 1` only). | The soundness statement is every constraint (`==`, `!=`, `<`, `<=`, `>`, `>=`) at every solution, conjoined with `d ≠ 0` for every denominator `d` of the system. Completeness needs one polynomial equation and nothing else. Other heads (`||`, ...) get no certificate. Audit of the other kinds: antiderivative and identity statements are cleared polynomial identities; matrix, optimality and Pratt statements have no division; sums with a denominator in the variable now get no certificate. |
| **Pratt fuel.** `powMod` had 256 bits of fuel, so `PrattCertifies (2^256 + 1) [(2^256 + 1, 3, [(2, 256)])]` held for the composite Fermat number F8. | The fuel is `e + 1`, which bounds the `log2 e + 1` halving steps for every `e`; the kernel unfolds only as many steps as there are bits. F8 is rejected; `2^64 + 1`, `1000003`, `51·2^257 + 1` still verify (well under 1 s), and `NextPrime[2^256]` verifies (Lean 0.8 s; Mathilda's certificate search takes about 117 s, so it is not in the default suite). The checker already required the factor list to multiply out to `p - 1` and every factor to be an earlier step. |
| **LeanVerify.** `axiom bad : False` used from an `example` or a `def`, and `example … := by native_decide`, gave `"Success" -> True`. | Sources that mention `axiom`, `native_decide`, `implemented_by`, `extern`, `unsafe`, `debug`, `#eval`, `elab`, `macro`, `syntax`, `run_cmd`, ... are refused unrun (`LeanVerify::unsafe`). `example`s are renamed to named defs, and every theorem / lemma / def / abbrev / named instance (namespace-qualified) is axiom-checked through an appended, line-anchored `#print axioms` footer. |
| **Temp-file race.** Predictable names, probe-then-open, dangling symlinks followed, a pre-created `.out` parsed as Lean's output. | A private `mkdtemp` directory per run, `O_EXCL | O_NOFOLLOW` source creation, and a check that the shell created the output file. |
| **Shape truncation.** `matVec [[1, 1]] [1] = [1]` holds (`zipWith` truncates). | Matrix statements state `Mathilda.HasShape` / `Mathilda.HasLength`. |
| **Free identifiers.** A symbol not bound by the statement (e.g. `default` in `Solve[x^2 == 1 && default x == 0, x]`) was emitted as a free Lean identifier and resolved to `Inhabited.default`. | Generated files set `autoImplicit false`; theorem statements are translated in closed mode, which refuses unbound symbols, names that contain `«`/`»`/control characters, and two symbols with the same Lean text. |
| **Heuristic `True`.** `Engine -> "Mathilda"` used the same `"Verified" -> True`. | It now gives `Missing["Heuristic", True \| False]`; the summary box already reads `Mathilda Heuristic`. |
| **Cache key.** A 64-bit hash of the theorem text. | The text itself is the key. |

## Soundness hardening, second round (v0.6519)

A second adversarial review found two more false positives at the root of the trust model
and two weaker holes. Regression tests: `test_review2_*` in `tests/test_lean_verify.c`,
the hostile-answer cases in `lean/test/Untrusted.lean`.

| Hole | Fix |
|---|---|
| **LeanVerify deny-list bypass.** `import Lean` plus a user tactic elaborator (`@[tactic Lean.Parser.Tactic.decide]`, or `attribute [tactic …]`, which no word list catches) added `theorem one_eq_zero : (1:Nat) = 0` under `debug.skipKernelTC` with a name built from strings: `"Success" -> True`, "does not depend on any axioms". | **Kernel replay.** Every Lean run of `Certify` and `LeanVerify` compiles the source (`-R <dir> -o <dir>/MathildaCert.olean`, in the private temp directory) and then runs `lean --run src/internal/lean/Replay.lean <olean>` -- through `lake env` for a project. The script reads the module, imports its imports, replays every declaration with `Lean.Kernel.Environment.replay` (what the toolchain's `leanchecker` does per module) and computes each declaration's axioms from the replayed environment (imported constants through Lean's precomputed axiom table). It runs none of the checked file's code, so neither its verdict nor its axiom report can be forged by that file; a `skipKernelTC` declaration fails with `(kernel) declaration type mismatch`. `"Success"` needs the replay, standard axioms for every declaration of the module and for every declared name. The deny-list stays as an early refusal (extended with `addDecl`, `skipKernelTC`, `command_elab`, `@[tactic`, …). `leanchecker` itself was measured (0.6-0.7 s) but reports no axioms, which LeanVerify needs from a trusted process; the script costs about 1.1 s on core Lean (0.5 s of it elaborating the script, 0.45 s importing `Init`), below the 3 s threshold, so `Certify` replays every run too (defence in depth: its sources contain no user code except deny-listed `Method` scripts). |
| **Oracle answer elaborated as Lean.** The tactics parsed Mathilda's answer as an arbitrary `term` and elaborated it; a malicious oracle answered `by_elab do …`, which added a kernel-unchecked theorem, registered it as a `grind` lemma (so `mathilda_sos` closed `y ≤ x` and `(1:Rat) = 0` followed) and ran IO. | The answer is **data**. `lean/Mathilda/Answer.lean` is a hand-written tokenizer and Pratt parser (a pure function on `String`) for numerals, `v<i>`, `+ - * / ^`, unary `-`, lists and tuples, with size bounds; the tactics build the certificate's terms from the AST with fixed quotations, numerals and the goal's own fvars. `Lean`Oracle` now answers `[1, -v1 - v2]` instead of `fun (v1 v2 : _) => (… : List _)`. `LeanVerify` also refuses `import Mathilda` and `set_option mathilda.*`. |
| **Vacuous hypotheses over ℚ.** `Certify[x == 7, Assumptions -> x^2 + 4 == 0, Method -> …]` proved `∀ x : Rat, 4 + x^2 = 0 → x = 7`, false over Mathematica's complex semantics (`x = 2 I`); likewise `x^2 == 2`. | **Generic domains.** Core Lean's `grind` works over its own algebraic classes, so statements with hypotheses quantify over every field of characteristic 0 (`[Lean.Grind.Field K] [Lean.Grind.IsCharP K 0]`: ℚ, ℝ, ℂ) and, with an order relation, over every ordered field (`[LE K] [LT K] [Std.LawfulOrderLT K] [Std.IsLinearOrder K] [Lean.Grind.OrderedRing K]`: ℝ). The statement is translated over `Rat` and generalised (`Rat` → `K`); the theorem takes the binders as parameters, the footer states the full `∀ {K : Type} …`. Ideal membership, identities with assumptions, SOS and optimality use it; every showcase row still verifies, and the "Rat → Real transfer not checked" trusted step is gone. `Solve` soundness statements are closed rational facts and stay over `Rat`. |
| **Together cancels; weak primality statements.** `Certify[(x^2 - 1)/(x - 1) == x + 1]` proved `∀ x : Int, 1 + x = 1 + x`; `Certify[Not[PrimeQ[0]]]` proved `0 < 2`. | Rational identities are cross-multiplied with the unreduced numerators and denominators (`rawFraction`), so Lean checks `x ^ 2 - 1 = (x - 1) * (1 + x)`. The prelude defines `Mathilda.IsPrime n := 2 ≤ n ∧ ∀ d, d ∣ n → d = 1 ∨ d = n`; non-primes are `¬ Mathilda.IsPrime n` (from a factor, or from `2 ≤ n`), and `PrimeQ[2]` is `Mathilda.IsPrime 2`, all with no trusted step. |

Residual trust: the replay trusts the imported compiled modules (core, Mathlib, the
`Mathilda` package) as built, and the toolchain itself. A LeanVerify source with `import
Lean` still runs arbitrary code at elaboration time with the user's rights (the replay keeps
it from faking a verdict, not from having side effects); the deny-list refuses the common
forms.

## Configuration

| Setting | Meaning |
|---|---|
| `$LeanExecutable` / option `LeanExecutable` | `Automatic` or a path. `Automatic` tries the environment variable `MATHILDA_LEAN`, then the newest `~/.elan/toolchains/*/bin/lean` (called directly, so no `lean-toolchain` file is needed), then `lean` on `PATH`. |
| `$LeanProject` / option `LeanProject` | `None` or a Lake project with Mathlib, for Tier 1. Defaults to the environment variable `MATHILDA_LEAN_PROJECT`. Tier 1 runs `cd <project> && <lake> env <lean> --json FILE`, with the `lake` beside the Lean binary. |
| `MATHILDA_BIN` | Lean side: the Mathilda binary that the `lean/` tactics call (else `set_option mathilda.path`, else `./Mathilda` or `../Mathilda`, else `PATH`). |
| `MATHILDA_LAKE` | test scripts: the `lake` to use for `lean/`. |

## Limitations

- Tier 0 statements are over `Int`/`Rat`, or over a generic char-0 / ordered field (statements
  with hypotheses, SOS, optimality). Real-analytic content (derivatives of `Log`,
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
- `LeanVerify` renames `example`s to named defs so their axioms are checked; anonymous
  `instance`s are covered only by the error / `sorry` checks, and sources using
  metaprogramming (`elab`, `macro`, `#eval`, …) are refused rather than run.

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
| F11 | A core-Lean field lemma from the raw rational statement to the cleared identity | Removes the remaining "field algebra" trusted step of rational identities (the cross-multiplied identity now uses the unreduced denominators); core Lean lacks `div_add_div`. |
