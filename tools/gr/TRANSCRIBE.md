# Transcribing Gradshteyn & Ryzhik integrals

This is the contract for the vision agents that read G&R (7th edition, 2007) page
images into machine records. Everything downstream depends on it: the numeric oracle
(`tools/gr/verify.m`, run in Mathilda by `tools/gr/verify.py`), the corpus builder
(`tools/gr/build_corpus.py`) and the triage (`tools/gr/triage.py`). Every check runs in
Mathilda; nothing in this pipeline uses another CAS.

The goal is to **read faithfully**. A record whose transcription is wrong never enters
the corpus, because the oracle rejects it. A *misread that happens to verify* is
impossible in practice. A *skipped* formula, however, is lost for good. So when in
doubt, transcribe it and let the oracle judge.

PDF page = printed page + 49. Read pages with the `Read` tool's `pages` parameter.

## Output

One JSON array per agent, written to `tools/gr/raw/p<first>-<last>.json`, where the
page numbers are PDF pages. Each element looks like this:

```json
{
  "id": "GR3.389.2",
  "gr": "3.389.2",
  "page": 400,
  "status": "ok",
  "lhs": "Integrate[x^(2*nu - 1)*(u^2 + x^2)^(rho - 1)*E^(-mu*x), {x, 0, Infinity}]",
  "rhs": ["u^(2*nu + 2*rho - 2)/(2*Sqrt[Pi]*Gamma[1 - rho]) * MeijerG[{{1 - nu}, {}}, {{1 - rho - nu, 0, 1/2}, {}}, mu^2*u^2/4]"],
  "conds": ["Abs[Arg[u]] < Pi/2", "Re[mu] > 0", "Re[nu] > 0"],
  "cite": "ET II 234(15)a",
  "edition_note": "7",
  "section": "3.38 Exponentials and arbitrary powers",
  "tex": "int_0^inf x^(2nu-1)(u^2+x^2)^(rho-1) e^(-mu x) dx = ... G^{31}_{13}",
  "note": ""
}
```

- **`id` / `gr`**: the G&R number plus the sub-item number (`3.389.2`). A formula
  number with no sub-items is just `3.381`. A sub-item that holds several
  equations (the 1/2/3 lines of one item, a ± pair, or a chained equality
  `∫ = A = B`, one record per form) gets the suffixes `a`, `b`, …: `GR2.148.3a`.
  `rhs` alternatives are only for ambiguous *readings* of one form, never for
  distinct forms.
- **`edition_note`**: the superscript after the item number (for example `7`, `11`
  or `*`; it marks a correction in that edition). Use `""` if there is none.
- **`page`**: the **PDF** page number (printed page + 49).
- **`lhs`**: an expression containing exactly **one** `Integrate[...]`, written in
  InputForm (the surface syntax Mathilda parses).
  - Definite: `Integrate[f, {x, lo, hi}]`.
  - Indefinite: `Integrate[f, x]`.
  - If the book writes `c ∫ … = value`, write `lhs` as `c*Integrate[...]`. The
    oracle solves for the integral.
- **`rhs`**: a list holding the value. Usually it has one element. If a reading is
  genuinely ambiguous (an unclear exponent or sign), list each candidate. The
  oracle keeps the first one that verifies.
- **`conds`**: the bracketed side conditions, one string per condition, in InputForm:
  - `Re μ > 0` → `"Re[mu] > 0"`
  - `|arg u| < π` → `"Abs[Arg[u]] < Pi"`
  - `n = 0, 1, 2, …` → `"Element[n, Integers] && n >= 0"`
  - `a > b > 0` → `"a > b > 0"`

  Conditions implied only by the result (`n != -1` when the value has
  `1/(n+1)`) may be added.  Conditions that only restate a branch convention (such as `arg ix = π/2 sign x`)
  go in `note` instead.
- **`cite`**: the source reference at the right margin (`ET I 138(14)`, `MO 213`,
  …).
- **`tex`**: a short plain-text reading of the formula, kept for the human-readable
  skip list.
- **`status`**: `"ok"`, or one of the skip statuses below. A skipped record still
  carries `id`, `gr`, `page`, `tex` and a `reason`.

| status | when |
|---|---|
| `SKIP-RECURRENCE` | Reduction or recursion formula: an integral appears on the right-hand side as well. |
| `SKIP-MULTIPLE` | Double or n-fold integral. |
| `SKIP-GENERIC` | Integrand contains an unspecified function f, φ, F, Z_ν (a generic cylinder function), etc. |
| `SKIP-PV` | Principal-value integral (marked PV or ⨍). |
| `SKIP-REFERENCE` | Entry is only "see 3.xxx", a definition, or a substitution rule with no closed value. |
| `SKIP-CONTOUR` | Contour or complex-path integral (∮, or limits that are complex paths). |
| `SKIP-UNSUPPORTED` | Uses a function with no Mathilda head (for example Lommel s_{μ,ν}, theta functions with unclear conventions, or Mathieu functions). **Still fill in lhs/rhs as well as you can.** |
| `SKIP-ILLEGIBLE` | Cannot be read. Use this sparingly. |

The general formula of a family that is followed by special cases is **not** a
recurrence. Transcribe each one.

## Notation map (G&R → InputForm)

**Elementary.**
- `log`/`ln` → `Log`; `lg` → `Log[10, x]`.
- `arc sin`/`arcsin` → `ArcSin`, and likewise for the other inverse trig functions.
- `arcsinh`/`arsinh`/`Arsh` → `ArcSinh`, and likewise for the other inverse hyperbolic functions.
- `e^x` → `E^x`; `i` → `I`; `π` → `Pi`.
- `C` (Euler's constant) → `EulerGamma`; `G` (Catalan's constant) → `Catalan`.
- Factorials `n!` → `n!`; double factorials `(2n-1)!!` → `(2*n - 1)!!`.
- `sign x` → `Sign[x]`.

**Gamma family.**
- `Γ(z)` → `Gamma[z]`; `Γ(α,x)` → `Gamma[a, x]`; `γ(α,x)` → `Gamma[a, 0, x]`.
- `B(x,y)` → `Beta[x, y]`.
- `ψ(x)` → `PolyGamma[x]`; `ψ'(x)` → `PolyGamma[1, x]`.
- `β(x)` (G&R 8.370) → `(PolyGamma[(x + 1)/2] - PolyGamma[x/2])/2`.
- `ζ(s)` → `Zeta[s]`; `ζ(s,q)` → `HurwitzZeta[s, q]`.
- `Φ(z,s,v)` (Lerch) → `LerchPhi[z, s, v]`.
- `B_n` → `BernoulliB[n]`; `E_n` → `EulerE[n]`.

**Error and Fresnel functions.**
- `Φ(x)` and `erf x` (probability integral, G&R 8.250) → `Erf[x]`; `erfc` → `Erfc`.
- `S(x)` and `C(x)` (G&R 8.251, normalised `2/√(2π) ∫₀ˣ sin t² dt`) →
  `FresnelS[Sqrt[2/Pi]*x]` and `FresnelC[Sqrt[2/Pi]*x]`.

**Exponential and trigonometric integrals.**
- `Ei` → `ExpIntegralEi`; `li` → `LogIntegral`.
- `si(x)` → `SinIntegral[x] - Pi/2`; `Si` → `SinIntegral`; `ci`/`Ci` → `CosIntegral`.
- `shi` → `SinhIntegral`; `chi` → `CoshIntegral`.

**Elliptic** (G&R uses the *modulus* k). Use these placeholder heads; the oracle
resolves the convention:
- `K(k)` → `MatharK[k]`; `E(k)` → `MatharE[k]`.
- `K'(k)` / `K(k')` → `MatharK1[k]`; `E'(k)` / `E(k')` → `MatharE1[k]`.
- `F(φ,k)` → `MatharF2[phi, k]`; `E(φ,k)` → `MatharE2[phi, k]`.
- Complete `Π(n,k)` → `MatharPi[n, k]`.
- Incomplete `Π(φ, n, k)` → `EllipticPi[n, phi, k^2]` (G&R 8.111.4 has `1 − n sin²`;
  write `-n` if the book has `1 + n sin²`).
- Jacobi `sn(u,k)` → `JacobiSN[u, k^2]`, and likewise `cn`, `dn`, `am` → `JacobiAmplitude[u, k^2]`.

**Bessel family.**
- `J_ν(z)` → `BesselJ[nu, z]`; `Y_ν`/`N_ν` → `BesselY`; `I_ν` → `BesselI`; `K_ν` → `BesselK`.
- `H^{(1)}_ν` → `HankelH1`; `H^{(2)}_ν` → `HankelH2`.
- Bold `𝐇_ν` → `StruveH`; bold `𝐋_ν` → `StruveL`.
- `𝐉_ν` (Anger) → `AngerJ`; `𝐄_ν` (Weber) → `WeberE`.
- `ber`/`bei`/`ker`/`kei` → `KelvinBer`/`KelvinBei`/`KelvinKer`/`KelvinKei`.
- `Ai`/`Bi` → `AiryAi`/`AiryBi`.
- `j_n`/`y_n` spherical → `SphericalBesselJ`/`SphericalBesselY`.

**Orthogonal polynomials.**
- `P_n(x)` → `LegendreP[n, x]`; `P^μ_ν(x)` → `LegendreP[nu, mu, x]`; `Q^μ_ν` → `LegendreQ[nu, mu, x]`.
  If the argument is off (−1, 1) (for example `√(1 + a²)`, `cosh t`, any `z > 1`), use
  the type-3 form `LegendreP[nu, mu, 3, z]` / `LegendreQ[nu, mu, 3, z]`. The 3-argument
  form is Ferrers (type 2) and is complex there.
- `C^λ_n` → `GegenbauerC[n, lam, x]`.
- `P^{(α,β)}_n` → `JacobiP[n, a, b, x]`.
- `L^α_n` → `LaguerreL[n, a, x]`; `L_n` → `LaguerreL[n, x]`.
- `H_n` → `HermiteH[n, x]`; `T_n`/`U_n` → `ChebyshevT`/`ChebyshevU`.

**Hypergeometric and related.**
- `F(α,β;γ;z)` / `₂F₁` → `Hypergeometric2F1[a, b, c, z]`.
- `Φ(α;γ;z)` / `₁F₁` / `M(a,b,z)` → `Hypergeometric1F1[a, c, z]`.
- `Ψ(a;c;z)` / `U` → `HypergeometricU[a, c, z]`.
- `pFq(a₁..;b₁..;z)` → `HypergeometricPFQ[{a1, …}, {b1, …}, z]`.
- `W_{κ,μ}(z)`/`M_{κ,μ}(z)` → `WhittakerW[k, m, z]`/`WhittakerM[k, m, z]`.
- `D_p(z)` → `ParabolicCylinderD[p, z]`.
- `G^{mn}_{pq}(z | a₁..a_p ; b₁..b_q)` →
  `MeijerG[{{a1..an}, {a(n+1)..ap}}, {{b1..bm}, {b(m+1)..bq}}, z]`.
- `Li₂(z)` / `L(z)` dilogarithm → `PolyLog[2, z]`. Lobachevsky `L(x)` → `SKIP-UNSUPPORTED`
  unless the book defines it inline.

**Derivatives and distributions.**
- A derivative of a special function with respect to its argument or order, such as
  `J'_ν(z)` or `∂J_ν/∂ν`, uses `Derivative`. For example `J'_ν(z)` →
  `Derivative[0, 1][BesselJ][nu, z]`.
- A result that is a distribution (δ, closure relations) gets `SKIP-UNSUPPORTED`, but
  fill in lhs/rhs (`DiracDelta`).
- A `±` inside a condition is split into two conditions.

**Unsupported functions.** Use these fixed placeholder heads, so that every batch
spells them the same way:
- `LommelS1[mu, nu, z]` (s_{μ,ν}) and `LommelS2[mu, nu, z]` (S_{μ,ν})
- `LommelU[nu, w, z]`, `LommelV[nu, w, z]`
- `NeumannO[n, z]`, `SchlaefliS[n, z]`
- `BatemanK[nu, z]`
- `HumbertPhi1[...]` (Φ₁, Ψ₁, Ξ₁, … accordingly)
- `ThetaN[n, v, q]` (ϑ_n)
- `MathieuCe[n, z, q]`, `MathieuSe[n, z, q]`
- `LobachevskyL[x]`
- `PolyLogLegendreChi[n, x]`
- `AppellF1[a, b1, b2, c, x, y]` (and `AppellF2`…`F4`), argument order as in DLMF 16.13

**A pole on the integration path** that the book does not mark PV still makes the
integral a principal value. Give it `SKIP-PV` with lhs and rhs filled in.

**A disjunctive condition** ("a > 0 or b < 0"): keep one branch that you can sample, and
note the other in `note`.

When the book gives a defining relation to a supported function (for example
`k_{2ν}(x) = W_{ν,1/2}(2x)/Γ(ν+1)`), **substitute it** and mark the record `ok`.

**Indexed parameter families** (β_k for k = 0..n under a Σ/Π with symbolic n) cannot be
instantiated. Give them `SKIP-GENERIC` with the reason "indexed parameter family".

**Conditions you cannot sample.** The oracle draws random parameter values. So:
- An *equality* condition (`a = b`, `b² = ac`, Δ = 0) is **substituted** into lhs and
  rhs, eliminating a parameter. Say so in `note`.
- A condition on the *integration variable* (`x > 0`, `sinh²x < a²`) goes in `note`, not
  in `conds`.
- A condition such as "p is not a negative integer" goes in `note`.

**Base integrals by reference.** When a value ends in `+ ∫dx/√R` together with "(see
2.261)", substitute the closed form of the referenced formula, choosing a branch.
Add that branch's condition (for example `c > 0`) and record the substitution in
`note`. Use SKIP-RECURRENCE only when the integral left over is a genuinely different
integral (a reduction formula).

**Page ownership.** Every numbered item that *starts* on one of your pages is yours,
even if it belongs to a section outside the one you were told about. Finish an item
that runs over onto the next page; leave alone an item that started on the previous
page.

**Book errors.** Transcribe the book *as printed*, even when you can see a misprint,
such as a parameter in the value that is absent from the integral, or a factor that
disagrees with the conditions. Add a corrected candidate as a second `rhs` entry and
explain it in `note`. The oracle decides.
A misprint in the *integral itself* (integrand, limits) cannot be offered as an `rhs`
candidate. Keep the as-printed record, and add a separate record with suffix `b` (or
the next free letter) that has the corrected lhs. Say "corrected integrand" in its
`note`.

**Symbols.**
- Greek parameters are spelled out in lower case: `nu`, `mu`, `lam`, `alpha`, `beta`, `gamma`, `rho`, `sigma`, `tau`, `phi`, `psi`, `omega`, `kappa`, `eta`, `theta`, `xi`, `zeta`, `epsilon`, `delta`.
- `p`, `q`, `a`, `b`, `c`, … stay as they are.
- **Never** use `C`, `D`, `E`, `I`, `K`, `N`, `O` as parameter names. Rename them: `C` → `CC`, `N` → `NN`, `K` → `KK`.
- Subscripted parameters: `a₁` → `a1`, `b_k` → `bk`.

**Ranges and syntax.**
- ∞ → `Infinity`.
- G&R's `Σ_{k=0}^{n}` → `Sum[..., {k, 0, n}]`.
- Products → `Product`.
- Use explicit `*` between factors and parenthesise every argument of `^`.

**Abbreviations.** G&R often abbreviates (`z = a + bx`, `Δ = …`, `R = √(a + bx + cx²)`)
at the head of a section. **Substitute the abbreviation in full** in both lhs and rhs.
Do not leave `z`, `R` or `Δ` standing. The section headings with the definitions are
on the pages around the formulas, so read back a page if needed.

## Self-check before writing

Run:

```
python3 -I tools/gr/verify.py check tools/gr/raw/<file>.json
```

It runs `tools/gr/check.m` in Mathilda, which parses every lhs/rhs/cond (held, so
nothing is evaluated) and reports syntax errors and stray abbreviation symbols.
Fix what it reports. Do **not** try to verify values yourself; the pipeline does that.
