# Mathilda and Lean 4: showcase

Mathilda computes an answer together with a *certificate*. It writes a self-contained
Lean 4 file that states the answer as a theorem, and the Lean kernel decides whether the
theorem holds. `showcase.m` runs every kind of certificate end to end and times it.

## Run it

```sh
make -j8                                   # build ./Mathilda
./Mathilda -file examples/lean/showcase.m  # about 15-21 s with a warm Lean
sh tests/scripts/lean_tactics.sh ./Mathilda   # the Lean -> Mathilda tactics on their own
```

You need a Lean 4 toolchain, v4.35 or newer (core `grind`). Install it with
[elan](https://github.com/leanprover/elan). Mathilda finds it through `$LeanExecutable`,
then the environment variable `MATHILDA_LEAN`, then the newest `~/.elan/toolchains/*/bin/lean`,
then `lean` on `PATH`. No Mathlib and no Lake project are needed for the default run.

The script exits with status 1 if any row has an unexpected verdict. Rows 17 and 18 are
`False` by design. Row 2 is Skipped unless a Mathlib project is configured, and row 20 is
Skipped when no `lake` is found.

## What it shows

| # | Row | Lean proves |
|---|---|---|
| 1 | Risch/Rioboo antiderivative of `(x^4-3x^2+6)/(x^6-5x^4+5x^2+4)` | the integrand equals the formal derivative of `ArcTan[x] + ArcTan[x^3] + ArcTan[(x-3x^3+x^5)/2]` (cleared identity, Tier 0) |
| 2 | Mathlib Tier 1 | `HasDerivAt (fun x => ArcTan x + Log (1+x^2)/2) ((1+x)/(1+x^2)) x` over ℝ |
| 3 | `Integrate[1/(x^2+2), x]` | the same, modulo `s1^2 - 2` with `s1` standing for `Sqrt[2]` |
| 4 | `ProofObject[Factor[x^10 - 1], Engine -> "Lean"]` | `∀ x : Int, x^10 - 1 = (x-1)(1+x)(…)(…)` |
| 5 | Bezout (`PolynomialExtendedGCD`) | `s a + t b = g`, and `g` divides `a` and `b` |
| 6 | Ideal membership | `x^2+y^2 = 1 → x = y → 2y^2 = 1` |
| 7 | `Solve` completeness | `p = 0 ↔ x = 1 ∨ x = 2 ∨ x = 3` |
| 8–10 | `Inverse`, `Det`, `NullSpace` | evaluated by the kernel (`decide +kernel`) |
| 11 | Sum of squares | `x^4 - 2x^2 + 1 + y^2 ≥ 0` |
| 12 | `Minimize` (QP) | optimality, from a KKT multiplier and an SOS remainder |
| 13 | `Maximize` (LP) | optimality; the multipliers `{3, 0, 0, 1}` are the dual solution |
| 14 | `FactorInteger[2^64 + 1]` | the product, and a Pratt certificate for each prime factor |
| 15–16 | `Sum[k^3]`, `Sum[k 2^k]` | the closed form, by induction |
| 17 | **Lean catches the CAS** | `PossibleZeroQ` says the polynomial `x/10^25` is zero; `Engine -> "Mathilda"` agrees, Lean refutes it |
| 18 | A tampered certificate | one coefficient of the `Factor` certificate changed: Lean rejects it |
| 19 | Batch | rows 4–10 in one Lean run, then a cache hit |
| 20 | Lean -> Mathilda | the `lean/` tactics ask Mathilda for certificates and check them in Lean |

## Recorded output

[`showcase_output.txt`](showcase_output.txt) is a full run on an Apple M-series machine
(Lean v4.35.0-rc3, no Mathlib). Total wall time was 15.0 s; the first Lean run of a
session is about 1–3 s cold, and later runs are about 0.5 s each. Excerpt:

```
[1] Risch/Rioboo antiderivative, certified in core Lean (Tier 0)
    In:        Certify[Integrate[(x^4 - 3*x^2 + 6)/(x^6 - 5*x^4 + 5*x^2 + 4), x]]
    Result:    ArcTan[x] + ArcTan[x^3] + ArcTan[(x - 3*x^3 + x^5)/2]
    Verified:  True    Tier: 0    Kind: Antiderivative    TrustedSteps: 4
    Time:      Lean 1.15 s, whole Certify call 1.15 s

[17] Lean catches the CAS: a PossibleZeroQ false positive
    In:        PossibleZeroQ[(1 + x)^20 - Expand[(1 + x)^20] + 10^-25 x]
    Result:    True   (wrong: the expression is the polynomial x/10^25)
    Engine -> "Mathilda":  Verified: True    Tier: Heuristic
    Engine -> "Lean":      Verified: False    Tier: 0

[19] Batch: rows 4-10 in ONE Lean run (session cache cleared first)
    Verified:  {True, True, True, True, True, True, True}
    Time:      one Lean run 0.65 s, whole batch 0.66 s  vs  3.34 s for the 7 separate calls above (speedup 5.1x)
    Repeat:    Certify[Factor[x^10 - 1]] again: Cached -> True in 0.00 s
```

With a Mathlib project configured, row 2 reads:

```
[2] Mathlib Tier 1: HasDerivAt (fun x => F) f x over the reals
    In:        Certify[Integrate[(1 + x)/(1 + x^2), x], Tier -> 1]
    Result:    ArcTan[x] + Log[1 + x^2]/2
    Verified:  True    Tier: 1    Kind: Antiderivative    TrustedSteps: 0
    Time:      Lean 32.29 s, whole Certify call 32.29 s
```

## Trust model

- Mathilda is an **untrusted oracle**, and the Lean kernel is the **trusted checker**. A bug
  in Mathilda can make a certificate fail. It cannot make Lean accept a false theorem.
- `"Verified" -> True` means: Lean exited 0, reported no error and no `sorry`, and the
  theorem's `#print axioms` lists only `propext`, `Classical.choice` and `Quot.sound`
  (so `native_decide` does not count).
- `"TrustedSteps"` lists what Lean did **not** check for that kind of certificate, such as
  the irreducibility of `Factor`'s factors, or Lucas's theorem behind a Pratt certificate.
  An empty list means the theorem is the whole claim. The summary at the end of the
  showcase prints them for every row.
- Tier 0 uses core Lean only, so it needs no Mathlib, and it inlines a small prelude into every file.
  Tier 1 (Mathlib) proves analytic statements such as `HasDerivAt`, with nothing trusted.
- The Lean source, claims and tactic strings only ever go into the `.lean` file. The shell
  command holds only quoted paths and fixed flags.

More detail is in [`docs/design/lean_integration.md`](../../docs/design/lean_integration.md).
The reference for every builtin and option is in
[`docs/spec/builtins/formal-verification.md`](../../docs/spec/builtins/formal-verification.md).

## Enabling Tier 1 (Mathlib)

Point Mathilda at a Lake project that has Mathlib built **for the same toolchain** as the
Lean executable:

```sh
MATHILDA_LEAN_PROJECT=/path/to/project ./Mathilda -file examples/lean/showcase.m
```

or, inside Mathilda, set `$LeanProject = "/path/to/project"`, or pass `LeanProject -> "/path/to/project"` to a single `Certify`.
Tier 1 runs `cd <project> && <lake> env <lean> --json FILE`, using the `lake` beside the Lean binary.
Each run pays for `import Mathlib`, which takes about 10–40 s, so put Tier 1 items in one batch:
`Certify[{…}, Tier -> 1]`. Mathilda only reads the project and never writes to it.

## Lean -> Mathilda quick start

```sh
cd lean
lake build
MATHILDA_BIN=$PWD/../Mathilda lake env lean test/Tactics.lean
```

```lean
import Mathilda

example (x y : Rat) (h1 : x^2 + y^2 = 1) (h2 : x = y) : 2*y^2 = 1 := by
  mathilda_linear_combination      -- Mathilda finds the cofactors, grind checks them

example (x y : Rat) : 2*x*y ≤ x^2 + y^2 := by
  mathilda_sos                     -- Mathilda finds the squares

example : ∃ c, Mathilda.PrattCertifies 274177 c := by
  mathilda_prime                   -- decide +kernel checks the Pratt certificate

#mathilda "Factor[x^10-1]"
```

See [`lean/README.md`](../../lean/README.md) for how the tactics call Mathilda and how they
find the binary (`MATHILDA_BIN`, `set_option mathilda.path`, `./Mathilda`, `PATH`).
