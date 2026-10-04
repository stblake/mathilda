---
references:
  - "NIST DLMF, Chapter 9 (Airy and Related Functions), §9.2 (defining equation and values) and §9.7 (asymptotic expansions)."
source: src/special_functions/airyai.c
---
**Algorithm.** `builtin_airyaiprime` computes `Ai'(z)`, the derivative of the Airy
function. It shares the one numeric engine, `airy_ai_core`, with `AiryAi`: the core
returns `Ai(z)` and `Ai'(z)` together, and the `prime` flag selects the derivative
component. The dispatch (`airyaiprime_one_arg`) is:

1. **Exact special values.** `AiryAiPrime[0] = -1/(3^(1/3) Gamma[1/3])` (built as an
   exact expression and evaluated); `AiryAiPrime[Infinity] = 0` (the recessive
   solution's derivative decays). `-Infinity` is deliberately left unevaluated
   because `Ai'` oscillates there with growing amplitude and has no limit, unlike
   `AiryAi`; `Indeterminate` maps to `Indeterminate`.
2. **Numeric evaluation (`USE_MPFR`).** A machine `Real` evaluates at 53 bits, an
   `EXPR_MPFR` at its own precision, and a `Complex` with an inexact part at the
   working precision implied by its parts. `airy_ai_core` picks a Maclaurin series
   for small `|z|`, the asymptotic expansion for large `|z|`, and the DLMF 9.2.12
   connection relation near the negative real axis; `Ai'` is assembled as
   `-z^{1/4} exp(-zeta)/(2 sqrt(pi)) · sumAp`.
3. Otherwise it returns `NULL` and stays symbolic (`AiryAiPrime[2]`,
   `AiryAiPrime[x]`), with the Maclaurin series produced by the generic
   Taylor-via-`D` path and `D[AiryAiPrime[z], z] = z AiryAi[z]` (the Airy equation)
   in `calculus/deriv.c`.

**Data structures.** An internal `acx` arbitrary-precision complex type (a pair of
`mpfr_t`) for the core; exact special values are ordinary `Expr` trees. The machine
real path has its own `double` kernel, `sf_machine_airy_ai_prime`, registered as a
unary (`REG_U`) ND kernel, so a packed or visible real `NDArray` runs element-wise
through it.

**Complexity / limits.** `O(1)` per element (the series/asymptotic term count is
bounded by the target precision). Attributes are `Listable | NumericFunction |
Protected`, and the head lowers in `Compile[]` at scalar and rank-1 real shapes
(`CompileDiagnostics` reports `Compiled -> True`). `-Infinity` is the one real
argument left unevaluated by design.
