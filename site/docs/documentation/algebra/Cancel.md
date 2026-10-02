# Cancel

!!! success "Status: Stable"
    documented, exercised by the test suite and/or worked examples, with no known limitations recorded.

## Description

**`Cancel[expr] cancels out common factors in the numerator and denominator of expr.`**

<details>
<summary>Notes</summary>

Option Extension -\> alpha cancels factors over Q(alpha) (e.g. simplifies (x^2 - 2)/(x - Sqrt\[2\]) to x + Sqrt\[2\] when Extension -\> Sqrt\[2\]). Default Extension -\> None treats algebraic numbers as opaque.

</details>

## Examples (14)

Every input below was run against the current Mathilda build and its output recorded.

### Basic examples (7)

```mathematica
In[1]:= Cancel[(x^2 - 1) / (x - 1)]
Out[1]= 1 + x

In[2]:= Cancel[(x - y)/(x^2 - y^2) + (x^3 - 27)/(x^2 - 9)]
Out[2]= (9 + 3 x + x^2)/(3 + x) + 1/(x + y)

In[3]:= Cancel[(y - 1)/(Sqrt[y] - 1)]
Out[3]= 1 + Sqrt[y]

In[4]:= Cancel[(y - 1)/(y^(1/3) - 1)]
Out[4]= 1 + y^(1/3) + y^(2/3)

In[5]:= Cancel[1/(y^(2/3) - 1/y^(1/3))]
Out[5]= y^(1/3)/(-1 + y)

In[6]:= Cancel[Sqrt[2]/(Sqrt[2] + Sqrt[2] x^4)]
Out[6]= 1/(1 + x^4)

In[7]:= Cancel[(Sqrt[2] x + Sqrt[2] y)/(Sqrt[2] + Sqrt[2] x)]
Out[7]= (x + y)/(1 + x)
```

### Options (2)

```mathematica
In[8]:= Cancel[(x^2 - 2)/(x - Sqrt[2]), Extension -> Sqrt[2]]
Out[8]= Sqrt[2] + x

In[9]:= Cancel[(x^3 - 2)/(x - 2^(1/3)), Extension -> 2^(1/3)]
Out[9]= 2^(2/3) + 2^(1/3) x + x^2
```

### Applications (5)

```mathematica
In[10]:= Cancel[(x^2 - 1)/(x - 1)]
Out[10]= 1 + x

In[11]:= Cancel[(x^2 + 2 x + 1)/(x + 1)]
Out[11]= 1 + x

In[12]:= Cancel[(x^2 - 1)/(x^2 - 2 x + 1)]
Out[12]= (1 + x)/(-1 + x)

In[13]:= Cancel[(x^4 - 1)/(x^2 - 1)]
Out[13]= 1 + x^2

In[14]:= Cancel[(x^2 - 2)/(x - Sqrt[2]), Extension -> Sqrt[2]]
Out[14]= Sqrt[2] + x
```

## Implementation notes

**Algorithm.** `Cancel` reduces a rational expression to lowest terms by cancelling the polynomial GCD of numerator and denominator. `builtin_cancel_compute` strips an optional `Extension -> α` (with `Automatic` via `extension_autodetect`, routing single-generator and tower cases through the `Q(α)` paths `cancel_with_extension` / `qa_cancel_with_tower`) and otherwise calls `cancel_recursive`. That routine recurses through `List`/`Plus`/relational/logical heads, then for a leaf fraction: splits the expression into `num`/`den` with `extract_num_den` (which understands `Rational`, `Complex`, `Power`/`Exp` with negative or split exponents, and `Times`), strips common symbolic atoms (`rat_strip_symbolic_common`), computes `g = PolynomialGCD[num, den]`, and divides both sides by `g` using fraction-free exact polynomial division (`exact_poly_div` over the collected variable set). A soundness gate (`has_embedded_rational_subterm`) leaves the input untouched when a `Power[Plus/Times, negative]` is present, since the multivariate Euclidean GCD could blow up; the denominator's leading sign is normalised to positive.

**Data structures.** `Expr` trees throughout; the cancellation calls into the polynomial subsystem via the `PolynomialGCD` builtin (dispatching to `qaupoly_gcd` for extensions) and `exact_poly_div`, with variable sets collected by `collect_variables` and sorted by `compare_expr_ptrs`. `Numerator`/`Denominator` (same file) expose `extract_num_den` directly.

**Complexity / limits.** Dominated by the multivariate polynomial GCD; the embedded-rational gate and leaf-count fallbacks exist to avoid pathological Euclidean-blowup hangs, so some algebraically-cancellable inputs are returned uncancelled (left for later `Simplify` passes).

- `Protected`, `Listable`.
- Threads over equations, inequalities, logic functions, and sums dynamically.
- Evaluates greatest common divisors via polynomial GCD derivations avoiding extraneous expansions.
- Handles a single symbolic base appearing with rational fractional exponents (e.g. `Sqrt[y]`, `y^(1/3)`) by treating it as an algebraic generator: substitutes `y -> g^m` where `m` is the LCM of denominators, runs the polynomial cancellation in `g`, then substitutes back.
- The algebraic-generator pass runs `Together` on the substituted form (not just GCD-cancellation), so inputs whose `g`-substituted denominator is a Plus of terms with different `g`-denominators (e.g. `1/(g^2 - 1/g)` from `1/(y^(2/3) - 1/y^(1/3))`) are handled correctly.
- Extracts algebraic-constant atoms (`Sqrt[2]`, `Sqrt[3]`, `CubeRoot[5]`, `2^(2/3)`, ...) that appear in every summand of *both* numerator *and* denominator and divides them out before the polynomial GCD step. Closes a long-standing gap where `PolynomialGCD[Sqrt[2], Sqrt[2] + Sqrt[2] x^4]` returns 1 (the integer-content recursion treats `Sqrt[2]` as having content 1), so cancellations whose only shared factor was an algebraic constant survived as-is. The pass is intentionally narrow — only `Power[integer, rational/non-integer]` factors are eligible — to avoid disturbing the rational-function intermediates that the integration dispatcher pattern-matches against.
- **Option `Extension -> alpha`** (Phase 0 of the Integrate plan) cancels common factors over `Q(alpha)` instead of `Q`. Implementation: lifts numerator and denominator into `Q(alpha)[x]` via the QAUPoly machinery, runs `qaupoly_gcd`, divides both sides by `g`, and re-renders. Works for single-fraction inputs; `Plus` inputs (sums of fractions) currently fall back to the no-extension path because `PolynomialQuotient` does not yet accept `Extension` (Phase 0.5 follow-up).
- **Native `AlgebraicNumber`-coefficient path.** A single univariate fraction whose coefficients are `AlgebraicNumber[θ, {..}]` over one number field `Q(θ)` (the ParallelMixedTower assembly / incidental-`Q(i)` DSolve representation) is reduced natively over `K(x)` as a pair of FLINT `gr_poly` over the antic number-field ring (`gr_poly_gcd` + exact division; divexact reporting success certifies the result). Value-equal and fully reduced; declines to the generic path on any doubt.
- **Native MULTIVARIATE number-field GCD (`flint_field_gcd`, v0.230).** Everything above is univariate — `gr_poly` over an antic `nf_t` is `K[x]` and nothing wider — so a gcd over `K[x_1..x_n]` fell through to the classical pseudo-remainder sequence, whose content is INTEGER content; a coefficient in `K` contributes content 1, so the K-content was never stripped and the first pseudo-remainder `lc(B)*A - lc(A)*B` vanished identically, making `PolynomialGCD[f, g]` return the SECOND OPERAND (`MATHILDA_DIVERGENCES.md` A26). The replacement is the classical modular (Encarnación) algorithm — FLINT has no multivariate GCD over a number field, so the work is done by `fq_nmod_mpoly_gcd` in the residue fields of `M mod p` (which must be *split* into its irreducible factors and CRT'd back, since for a non-cyclic Galois group such as `Q(√2, √3)` no prime keeps `M` irreducible), then CRT across primes and rational reconstruction. Every result is certified before return: the candidate is monic so its leading monomial is tau-free while `M`'s is `tau^n`, which makes `{G, M}` a Gröbner basis and `fmpq_mpoly_divrem_ideal` an exact decision procedure for divisibility over `K`. Both coefficient spellings are accepted, radicals (including nested ones) being mapped into one common field — cached by atom set, since building it is a qqbar primitive-element search that dominated the cost — and rendered back through the product basis of the caller's own atoms (radicals in, radicals out); a *parametric* radical like `Sqrt[k]` is declined and stays with `flint_parametric_sqrt_gcd`. The answer is the MONIC associate — a gcd over a field is fixed only up to a constant of that field. Anything the engine declines has its classical-path answer checked for divisibility and replaced by `1` if it fails. `MATHILDA_NO_FIELD_GCD=1` disables it.
- **Stress-hardened (v0.232).** The failure mode of a modular GCD is a *decline*, not a wrong answer: the post-check above then returns `1`, which is a valid common divisor and so passes everything downstream while the real gcd is lost. Pushed on coefficient size, field degree, term count and variable count, the v0.230 engine turned out to do exactly that in two whole classes. Its coefficient ceiling was **~831 bits** — `FG_MAX_PRIMES 64` × 29-bit primes is 1856 bits of modulus, and rational reconstruction needs about twice the coefficient size — now past 13,000 bits, with 62-bit primes and the cap demoted to a backstop (safe because the certificate, not the budget, is what makes the answer correct, so extra primes only cost time). And **every degree-6 radical field** declined, from a cause upstream in the qqbar compositum: it chose a primitive element by trial membership, which does not escalate past 64 bits of working precision when the generator has degree ≤ 6, so degrees 2–5 (resolved inside 64 bits) and 7+ (allowed to escalate) worked while 6 alone failed; the primitive element is now chosen by *degree*, since `alpha + c*b` always lies in `Q(alpha, b)` and so generates it exactly when the degrees agree. The product-basis render-back also gave up whenever the operands mentioned more atoms than a basis needs — `Expand` folds `Sqrt[2] Sqrt[3]` into `Sqrt[6]`, so a `Q(√2, √3)` problem has three atoms and a product basis of 8 against `[K:Q] = 4` — and now selects atoms greedily until the product reaches `n`. Speed came from choosing primes that split *least*, one multivariate gcd running per irreducible factor of `M mod p`: 2.4–2.6× less residue-field work, and 628-term operands from 63 ms to 30 ms. `MATHILDA_FIELD_GCD_STATS=1` prints the per-stage profile; `tests/bench_field_gcd.c` is the standing gate, every case built as `d·u`, `d·v` and required to come back an associate of `d`.
- **`Root[]`-object generators (A14b).** The extension generator resolved by `qa_resolve_extension` — shared by `Cancel`, `Together`, `PolynomialGCD`/`LCM`/`Quotient`/`Remainder`/`QuotientRemainder`, `Factor`, and `IrreduciblePolynomialQ` under both `Extension -> alpha` and `Extension -> Automatic` — now recognises a `Root[Function[...], k]` object alongside `Sqrt`/`Power[c, p/q]`/`I`/roots-of-unity/nested radicals. The Root's defining polynomial (both the named-`Function[t, body]` and Slot-`Function[body]` shapes) becomes the minimal polynomial of `Q(c)`, made monic over `Q`; the result renders in terms of the Root object (index `k` preserved). E.g. `PolynomialGCD[x^5 - x + 1, tau0 - c, Extension -> Automatic]` with `c = Root[2869 #1^5 + … &, 1]` now returns the genuine degree-1 gcd over `Q(c)` instead of the trivial `1`. A **reducible** Root defining polynomial is declined (the abstract field needs an irreducible relation), so such an input falls back to the gcd over `Q` rather than risking a wrong answer.

**Attributes:** `Listable`, `Protected`.

## References

**See also:** [Together](../../algebra/Together/), [Plus](../../arithmetic/Plus/), [PolynomialQuotient](../../algebra/PolynomialQuotient/), [AlgebraicNumber](../../algebra/AlgebraicNumber/), [Expand](../../algebra/Expand/), [PolynomialGCD](../../algebra/PolynomialGCD/), [LCM](../../number-theory/LCM/), [Quotient](../../arithmetic/Quotient/)

- von zur Gathen & Gerhard, "Modern Computer Algebra", on polynomial GCD computation.
- Geddes, Czapor & Labahn, "Algorithms for Computer Algebra" (1992), on rational function simplification.
- Source: [`src/rat.c`](https://github.com/stblake/mathilda/blob/main/src/rat.c)
- Specification: [`docs/spec/builtins/algebra.md`](https://github.com/stblake/mathilda/blob/main/docs/spec/builtins/algebra.md)
- Tests: [`tests/test_crc_corpus.c`](https://github.com/stblake/mathilda/blob/main/tests/test_crc_corpus.c)
- Tests: [`tests/test_dsolve_m60_stress.c`](https://github.com/stblake/mathilda/blob/main/tests/test_dsolve_m60_stress.c)
- Tests: [`tests/test_extension_auto_builtins.c`](https://github.com/stblake/mathilda/blob/main/tests/test_extension_auto_builtins.c)
- Tests: [`tests/test_extension_options.c`](https://github.com/stblake/mathilda/blob/main/tests/test_extension_options.c)

## Notes & additional examples

### Notes

Cancel removes the common polynomial factors between numerator and denominator by
dividing out their GCD, without performing a partial-fraction split. It reduces
`(x^2-1)/(x-1)` to `1 + x` and recognises perfect-square numerators such as
`(x^2+2x+1)/(x+1)`. When a common factor remains after cancellation, the result
stays a reduced quotient: `(x^2-1)/(x^2-2x+1)` becomes `(1+x)/(-1+x)` because both
share the factor `(x-1)` but the leftover `(x+1)/(x-1)` is already in lowest
terms. Cancel does not expand the surviving factors back out.

With the `Extension -> alpha` option, cancellation is performed over the
algebraic field `Q(alpha)`: `(x^2 - 2)/(x - Sqrt[2])` factors as
`(x - Sqrt[2])(x + Sqrt[2])` over `Q(Sqrt[2])`, so the denominator divides out
and the quotient collapses to `Sqrt[2] + x`. The default `Extension -> None`
treats algebraic numbers as opaque and would leave that quotient intact.
