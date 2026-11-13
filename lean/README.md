# Mathilda tactics for Lean 4

A small, core-only Lean 4 package (`import Lean`, no Mathlib) whose tactics ask
the Mathilda computer algebra system for a certificate and then check that
certificate in Lean.

**Trust model.** Mathilda is an untrusted oracle. It only ever supplies data:
cofactors, squares, or a Pratt certificate. The answer is never elaborated as
Lean syntax: `Mathilda/Answer.lean` parses it with a small hand-written parser
into numerals, variables `v<i>`, `+ - * / ^`, lists and tuples, and refuses
anything else. Each tactic builds terms from that data and turns them into
ordinary Lean proof steps (`have … := by grind`, `decide +kernel`), and the Lean
kernel checks them. A wrong, malformed or malicious answer makes the tactic
fail, and nothing in it is run. It never produces a proof of a false goal
(`test/Untrusted.lean` runs the tactics against a fake Mathilda that returns
wrong certificates and hostile answers such as Lean code).

| Tactic / command | Goal | Mathilda supplies | Lean checks with |
|---|---|---|---|
| `mathilda_linear_combination` | `a = b`, given `hᵢ : lᵢ = rᵢ` | cofactors: `a - b = Σ cᵢ (lᵢ - rᵢ)` | `grind` on the identity, then rewriting with `lᵢ - rᵢ = 0` |
| `mathilda_sos` | `0 ≤ p`, `p ≤ q`, `p ≥ q`, given `0 ≤ gᵢ` | `(d, ℓ, i)`: `q - p = Σ d ℓ² gᵢ` | the hints `0 ≤ ℓ²`, `0 ≤ ℓ² gᵢ` and `grind` |
| `mathilda_prime` | `∃ c, Mathilda.PrattCertifies p c` | a flat Pratt certificate | `decide +kernel` |
| `#mathilda "expr"` | — | `expr` evaluated, in `LeanForm` | nothing (an info message; it proves nothing) |

## Quick start

```sh
# 1. Build Mathilda (from the repository root).
make -j8

# 2. Build the package (lean-toolchain pins leanprover/lean4:v4.35.0-rc3).
cd lean
lake build

# 3. Check the examples.
MATHILDA_BIN=$PWD/../Mathilda lake env lean test/Tactics.lean
```

`tests/scripts/lean_tactics.sh` (ctest `lean_tactics_tests`) does steps 2–3 and
also runs `test/Untrusted.lean` and `test/Missing.lean`.

```lean
import Mathilda

example (x y : Rat) (h1 : x^2 + y^2 = 1) (h2 : x = y) : 2*y^2 = 1 := by
  mathilda_linear_combination

example (x y : Rat) : 2*x*y ≤ x^2 + y^2 := by
  mathilda_sos

example (x : Rat) (h : 0 ≤ x) : 0 ≤ x^3 := by
  mathilda_sos

example : ∃ c, Mathilda.PrattCertifies 274177 c := by
  mathilda_prime

#mathilda "Factor[x^10-1]"
-- (x - 1) * (1 + x) * (1 - x + x ^ 2 - x ^ 3 + x ^ 4) * (1 + x + x ^ 2 + x ^ 3 + x ^ 4)
```

To use the package from another Lake project, add it as a path dependency and
`import Mathilda`. The tactics work over any type with `Lean.Grind.CommRing`
(`mathilda_sos` also needs `Lean.Grind.OrderedRing`): `Int`, `Rat`, and `ℝ` when
the project has Mathlib.

## Finding the Mathilda binary

In order:

1. the environment variable `MATHILDA_BIN`;
2. `set_option mathilda.path "/path/to/Mathilda"`;
3. a `Mathilda` file in the working directory or its parent (so running
   `lake env lean …` inside `<repo>/lean` finds `<repo>/Mathilda`);
4. `Mathilda` on `PATH`.

If the binary cannot be run, the tactic fails with an error that names these
settings. For VS Code, export `MATHILDA_BIN` before starting the editor, or put
the `set_option` at the top of the file.

## How a call works

1. `Mathilda.Serialize` turns the goal and the usable hypotheses into Mathilda
   syntax. It accepts `+ - * /`, unary `-`, `^` with a literal exponent,
   numerals, `Nat.cast` / `Int.cast` and local variables, which are sent as
   `v1, v2, …`. A hypothesis that does not serialize is skipped. A goal that
   does not serialize is an error.
2. `Mathilda.Oracle.callMathilda` writes a temporary script
   `Print["@@MATHILDA-BEGIN@@"]; Print[Lean`Internal`Oracle[…]]; Print["@@MATHILDA-END@@"]`,
   runs `Mathilda -file <script>` with a null stdin, takes the text between the
   sentinels, and deletes the script. A call takes about 0.05 s.
3. Mathilda answers with data, for example `[1, -v1 - v2]`.
   `Mathilda.Answer.parse` reads it (bounded: 10⁶ characters, numerals of at
   most 2000 digits, exponents at most 10⁴, nesting at most 100); the tactic maps
   each `v<i>` to the goal's i-th variable and builds the terms with fixed
   quotations, elaborated at the goal's type. An answer that does not parse, or
   names a variable the goal does not have, fails with
   `<tactic>: Lean could not parse the certificate returned by Mathilda (…)`.
   An answer `MATHILDA-FAIL: <reason>` becomes the error
   `mathilda: Mathilda found no certificate: <reason>`.
4. The tactic runs the proof steps with error recovery switched off. If Lean
   rejects them, the state is restored and the tactic fails with
   `<tactic>: Lean rejected the certificate returned by Mathilda (…)`.

## Files

| File | Content |
|---|---|
| `Mathilda.lean` | re-exports the modules below |
| `Mathilda/Prelude.lean` | byte-for-byte copy of `src/internal/lean/Prelude.lean` (`prattCheck`, `PrattCertifies`, `sqnn`, `mulnn`, …); the ctest `lean_prelude_sync` enforces this |
| `Mathilda/Serialize.lean` | `toMathilda : Expr → …` |
| `Mathilda/Oracle.lean` | `callMathilda`, the `mathilda.path` option |
| `Mathilda/Answer.lean` | the parser for the oracle's answer (data only) |
| `Mathilda/Tactics.lean` | the three tactics and `#mathilda` |
| `test/Tactics.lean` | positive examples, `#print axioms` of each, false goals |
| `test/Untrusted.lean`, `test/bad_oracle.sh` | wrong certificates from a fake oracle are rejected; hostile answers (`by_elab` code, `sorry`, foreign identifiers, huge exponents / numerals, malformed tokens) are refused unread |
| `test/Missing.lean` | a missing binary is a clear error |

Each call starts a fresh Mathilda process. A persistent oracle over Mathilda's
NDJSON pipe mode is a follow-up; it only pays off for proofs with many tactic
calls.
