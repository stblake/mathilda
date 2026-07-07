# Cancel improvement plan — closing the Goursat `Sqrt[k]` named integrals

**Status:** re-diagnosed 2026-06-30 with stack samples. The original
diagnosis (below, §0) blamed `extension_autodetect`; direct measurement
**falsifies** that. The real blocker is the **multivariate `PolynomialGCD`
blow-up** reached through `Cancel`'s *Phase E* symbolic-radicand path.
**Owner hook:** the two classic Goursat named integrals (see
`goursat/goursat_1887.tex` §Examples) and any algebraic integrand whose
splitting field is `Q(params)(Sqrt[param])` rather than a number field.

---

## 0. Superseded original diagnosis (kept for the record)

The first writeup claimed `Cancel[..., Extension -> Automatic]` hangs inside
`extension_autodetect` (`src/poly/qafactor.c`) because `Sqrt[k]` is algebraic
over the *transcendental* `Q(k)`, not over `Q`. **This is wrong.** Static
tracing + measurement show `extension_autodetect` on a `Sqrt[k]` input returns
`NULL` **immediately**: `qa_resolve_extension_tower([Sqrt[k]])` →
`qa_resolve_extension` → G8 `qa_resolve_nested_radical` →
`expr_collect_atomic_algebraics(k)` rejects the free symbol `k` and returns
`NULL` in O(1). The proposed options A/B (build `y^q - g`; prefilter to the
no-extension path) were aimed at the wrong subsystem. Option B would in fact
*hurt*: plain `Cancel` returns the input **unchanged** here (see §2), so routing
there stops the hang but produces a useless non-canonical result.

---

## 1. Symptom (unchanged)

`Cancel[..., Extension -> Automatic]` / `Together[..., Extension -> Automatic]`
(hence `canonic` in `src/calculus/integrate_goursat.c`) does not terminate
within any reasonable budget on the V4-descended cofactors of the Goursat
integrals — radicals whose radicand is a **free symbol** (`Sqrt[k]`).

### Measured (commit `2a27588`, this tree, macOS `sample`)

| Call | Input | Result | Time |
|------|-------|--------|------|
| `Cancel[e]` (plain)                 | `(a+b√k−a²k)/(b−a√k+b²k−a²k²)` | returns **input unchanged** | 0.0 s |
| `Cancel[e, Extension -> Sqrt[k]]`   | same | returns **input unchanged** (correct: coprime over Q(√k)) | 0.0 s |
| `Cancel[e, Extension -> Automatic]` | same | **hangs** | > 150 s (killed) |
| `Together[1/(b−a√k)+1/(b+a√k), Extension -> Automatic]` | small | `(2b)/(b²−a²k)` ✔ | 0.0 s |
| `Together[…, Extension -> Sqrt[k]]` (explicit) | small | **fails** — leaves `PolynomialLCM[…, Extension->Sqrt[k]]` unevaluated | 0.0 s |

Two facts to hold onto:
1. The auto path's Phase E **produces the correct rationalised answer on small
   inputs** (`(2b)/(b²−a²k)`). It is *correct*; it only *blows up by size*.
2. Explicit `Extension -> Sqrt[k]` makes **`Cancel` fast** (it routes through
   `cancel_with_extension`, which never forms a many-variable GCD) but makes
   **`Together` fail** (`PolynomialLCM` has no symbolic-radicand extension
   branch). So neither plain nor explicit is a drop-in for the auto path.

---

## 2. Root cause — multivariate `PolynomialGCD`, reached via Phase E

`builtin_cancel_compute` (`src/rat.c:995`) handles `Extension -> Automatic`:
`extension_autodetect(arg)` returns `NULL` fast (§0), `alpha` stays `NULL`, and
control reaches **Phase E** `qa_cancel_with_poly_radical` (`rat.c:1142`,
impl `src/poly/qafactor.c:4168`). Phase E implements
"substitute `α→S`, `Together`, reduce mod `S^q − radicand`, rationalise via
`PolynomialExtendedGCD`, final `PolynomialGCD`, substitute back". Every one of
those steps runs in the **full multivariate ring** `Q[a, b, k, S]`.

### Stack sample — live `Integrate` Example 1 (`sample`, abbreviated)

```
builtin_integrate → gs_guarded → TimeConstrained → gs_core → goursat_v4
  → canonic (integrate_goursat.c:143) → internal_cancel → builtin_cancel
    → qa_cancel_with_poly_radical                 (Phase E)
      → builtin_polynomialextendedgcd  poly.c:4699   (Step 5.5 rationalise)
        → poly_div_rem → internal_together → builtin_together
          → together_recursive → cancel_recursive
            → builtin_polynomialgcd  poly.c:2327
              → poly_gcd_internal → poly_content → poly_gcd_internal → …
                → pseudo_rem / is_zero_poly_depth / expr_expand   ← blow-up
```

`goursat_v4`'s `canonic` accounts for essentially the whole runtime; inside it,
`PolynomialExtendedGCD`'s internal `Together` triggers a `PolynomialGCD` over
`Q[a,b,k,S]`, and `poly_gcd_internal`'s recursive content/pseudo-remainder
descent (`poly.c:1809`) explodes — coefficient growth in ≥ 3 parameters, plus
`is_zero_poly_depth` going quadratic on the giant expanded intermediates.

### Why the existing guards don't catch it

`poly_gcd_internal` has a pseudo-remainder size budget + 50-iteration cap
(`poly.c:1865–1886`), and `builtin_polynomialgcd` has a `multivar_gcd_skipped`
guard (`v_count >= 3 && total_leaves > 200`, `poly.c:2314`). Both are **leaky
for this case**: the offending `poly_gcd_internal` call has inputs *under* 200
leaves, so it is not skipped; the explosion then happens *inside* its
**recursive `poly_content` descent** (line 1819/1827) and the
pseudo-remainder sequence, where the size budget guards the `while` loop but not
the per-content recursion that precedes it.

### The asymmetry, explained

Explicit `Extension -> Sqrt[k]` `Cancel` is fast because `cancel_with_extension`
works in `Q(α)[mainvar]` with QA arithmetic for `α = √k` — it never builds a
4-variable `Q[a,b,k,S]` GCD. Phase E throws that structure away and reduces the
whole problem to multivariate `PolynomialGCD`, which is where Mathilda is weak.

---

## 3. Fix options (re-planned)

The goal is not merely "stop the hang" — it is "**close the two integrals**",
i.e. produce the canonical rational form `canonic` needs without blowing up.
Ordered by leverage:

### A. Fast multivariate `PolynomialGCD` (root-cause, high value, large)
Replace the recursive content/pseudo-remainder `poly_gcd_internal` with a
**modular / evaluation-homomorphism GCD** (Brown's dense modular GCD, or a
GCDHEU heuristic-evaluation GCD as a fast pre-check). Evaluate the parameters
`a,b,k` at integers (and reduce mod a prime), take univariate GCDs in `S`,
interpolate/CRT back, verify by trial division. This removes the coefficient
explosion that blocks Phase E **and** the conjugate path's final GCD, and is
broadly valuable across `Together`/`Cancel`/`Apart`/`Factor`. Risk: it is core
polynomial-engine surgery; needs careful content/leading-coeff handling and a
regression sweep. This is the only option that makes the *current* Phase E
algorithm terminate as designed.

### B. q=2 conjugate rationalisation in Phase E + bounded final GCD (targeted)
For `q == 2` (both Goursat integrals), replace Step 5.5's
`PolynomialExtendedGCD` with the **explicit conjugate**: with
`den_reduced = P + Q·S` (deg ≤ 1 in S), multiply num and den by `P − Q·S`;
`den_new = P² − Q²·radicand` (S-free, **computed directly**, no ExtendedGCD,
no product-`Together`). This kills the Step-5.5 explosion. **But measurement
shows the Step-6 final `PolynomialGCD[num_new, den_new, S]` over `Q[a,b,k]`
*also* hangs** — so B must be paired with making that final GCD bounded
(tighten / correctly place the `multivar_gcd_skipped` guard so it fires here),
accepting an under-reduced result when the GCD can't be proven cheaply.
Caveat: under-reduction may leave `canonic` non-minimal; whether the Goursat
descent still closes on a non-minimal `canonic` is **unverified** and must be
tested, not assumed.

### C. Hang-only guard (cheap, does NOT close the integrals)
Add an upfront variable-count + size bail to `poly_gcd_internal` /
`builtin_polynomialgcd` (and thus Phase E) so the auto path returns the
*un-canonicalised* input instead of hanging. This restores responsiveness but,
per §1 fact, yields a non-canonical result — the integrals stay **open**. Only
worth doing as an interim safety net under a larger A/B effort.

**Recommendation:** **A** is the real fix and the only one guaranteed to close
the integrals via the existing pipeline. **B** is a smaller bet but is gated on
an unverified assumption (Goursat closing on an under-reduced `canonic`) and
still needs a bounded final GCD. **C** alone does not meet the goal. Given A is
core-engine surgery with regression surface, confirm scope before starting.

### New helpers (A)
- A modular GCD entry point alongside `poly_gcd_internal`
  (`src/poly/poly.c`), e.g. `poly_gcd_modular(A, B, vars, var_count)`, with a
  dense-interpolation driver and a univariate-mod-p kernel; dispatch from
  `builtin_polynomialgcd` when `var_count >= 2` and inputs exceed a small
  threshold, falling back to `poly_gcd_internal` for tiny/edge cases.

---

## 4. Tests

```
/* Closes (A; or B if the under-reduced form still integrates): */
Integrate[(k x^2 - 1)/((a k x + b)(b x + a) Sqrt[x (1-x)(1-k x)]), x]
Integrate[(k x^2 - 1)/((a k x + b)(b x + a) Sqrt[(1-x^2)(1-k^2 x^2)]), x]
/* Terminates fast and stays correct (regression guard): */
Cancel[(a + b Sqrt[k] - a^2 k)/(b - a Sqrt[k] + b^2 k - a^2 k^2), Extension -> Automatic]
Together[1/(b - a Sqrt[k]) + 1/(b + a Sqrt[k]), Extension -> Automatic]  (* -> (2b)/(b^2 - a^2 k) *)
/* Number-field cases unchanged: */
Cancel[(1 + Sqrt[2])/(1 - Sqrt[2]), Extension -> Automatic]   (* -> -3 - 2 Sqrt[2] *)
Cancel[x/(x - (-1)^(1/3)), Extension -> Automatic]            (* cyclotomic, unchanged *)
/* Direct GCD micro-benchmarks (A): each well under 1 s */
PolynomialGCD[a + b S - a^2 k, b - a S + b^2 k - a^2 k^2, S]
PolynomialGCD[<conjugate num_new>, b^2 - a^2 k + 2 b^3 k - 2 a^2 b k^2 + b^4 k^2 - 2 a^2 b^2 k^3 + a^4 k^4, S]
```

Wall-clock assertion on each. The number-field rows guard against A/B/C
perturbing the algebraic-over-Q path.

---

## 5. Downstream payoff — the Goursat named integrals

```
Integrate[(k x^2 - 1)/((a k x + b)(b x + a) Sqrt[x (1-x)(1-k x)]), x]
  -> -2 ArcTan[Sqrt[a+b] Sqrt[b+a k] x / (Sqrt[a] Sqrt[b] Sqrt[(-1+x) x (-1+k x)])]
       / (Sqrt[a] Sqrt[b] Sqrt[a+b] Sqrt[b+a k])

Integrate[(k x^2 - 1)/((a k x + b)(b x + a) Sqrt[(1-x^2)(1-k^2 x^2)]), x]
  -> ArcTanh[...] form  (Goursat 1887 §Examples, R = (1-t^2)(1-k^2 t^2))
```

Both are **in scope** for the existing square-root V4 descent
(`goursat_v4`/`reduce_v4_piece`) — Goursat's Theorem 2 criterion
`F + F(S1) + F(S2) + F(S3) == 0` holds (verified), and the descent reduces them
to `Integrate[ G(u)/Sqrt[quadratic(u, Sqrt[k])], u ]` with `G` rational. The
sole blocker is that `canonic` over `Q(a,b,k,Sqrt[k])` invokes the multivariate
`PolynomialGCD` blow-up above. Fixing that (A) makes `canonic` fast with **no
Goursat-side change**.

> A localised `gs_no_extension` flag in `integrate_goursat.c` (plain Cancel in
> the V4 descent) was prototyped and **reverted**: it skips the extension
> entirely, so `canonic` returns a non-canonical form and the descent does not
> close. The fix belongs in the polynomial GCD / Cancel layer, not a per-site
> flag.

---

## 6. Related notes

* `project_cyclotomic_extension_support` (memory) — made the *cyclotomic*
  (`Q(ζ_n)`) Together/Cancel fast. Distinct subsystem from this multivariate-GCD
  issue, though both surface at the same `goursat_v4` entry point.
* `SIMPLIFY_IMPROVEMENT_PLAN.md` §Phase 5 — the `reduce_v4_piece` "tower growth"
  diagnosis is the cyclotomic square-root case (ζ₆→ζ₁₂); this document is the
  **symbolic-radicand** square-root case (`Sqrt[k]`), and its true bottleneck is
  multivariate `PolynomialGCD`, not extension discovery.
* `feedback_together_no_expand` (memory) — keep Together/Cancel pre-passes from
  expanding `Power[Plus[...], n]`; any Phase E / GCD change must preserve that.
* Measurement method: macOS `sample <pid>` on the hung process (no `timeout` on
  Darwin; pipe-buffering hides `Print` output, so sample the live task).
