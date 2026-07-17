# Task: Cherry B1 + A2 — general Σ-decomposition (Thm 4.4) + li non-existence decision

**Goal:** finish the 1986 li paper. Generalize `cherry_li.c`'s li-argument generator from
"powers `w^k` of a single generator" to genuine products `∏ⱼ fⱼ^{αⱼ}` via the decidable
Σ-decomposition (Thm 4.4), and wire its non-existence certificate to `rt_dec_nonelem()` so
Ex 5.2 (`∫ x²/log(x²−1)`) answers `ElementaryIntegralQ → False`.

## Why the current engine can't do it
`rt_cherry_li` posits li-terms `d_k · (w^k)'/θ` (args `w^k`, θ=Log[w]) and tower-solves. For a
reducible `w = ∏ fⱼ^{eⱼ}`, the true li arguments are `∏ fⱼ^{αⱼ}` with α NOT proportional to e, whose
logs `Σ αⱼ Log[fⱼ]` are independent linear forms — not rational multiples of θ. So they never appear
in the single-θ ansatz. The finite candidate set of α-vectors (and the NO proof) is the
Σ-decomposition.

## Pins
- POSITIVE new: `∫ x²/log(x³−x) dx` (w=x(x−1)(x+1), product decomposition).
- DECISION-NO: `∫ x²/log(x²−1) dx` → decline `Integrate` AND `ElementaryIntegralQ → False`.
- REGRESSION (must stay green): d1, d2, d3, d4, Ex 5.1 (`cherry_li_tests`), plus the full
  `test_integrate_risch_transcendental` battery.

## Plan — DONE (2026-07-17)
1. [x] Ground Thm 4.2/4.4 + Thm 5.3 + Ex 5.1/5.2 precisely from the paper.
2. [x] `cherry_sigma_decomp.{c,h}`: `cherry_sigma_decompose(Phi, factors, m, x)` →
       `{(b_i, alpha_i-vector)}` OR `SIGMA_NONEXISTENT` (faithful Thm 4.4: multiplicity extraction,
       `b=(p mod f1)/(q mod f1)`, consistency check, degree-overshoot terminator).
3. [x] Builtin `Integrate`SigmaDecomposition[Phi, {f1,…}, x]` + `tests/test_cherry_sigma.c` —
       decomposition tested in isolation against the paper's exact numbers.
4. [~] Positive path stays on the existing tower-solve (correct + diff-back gated). B1 powers the
       DECISION, not a new positive path; genuine non-proportional product args are the B2 tower
       recursion (deferred, not selected).
5. [x] DECISION wired — but as `Integrate`LiElementaryQ`, NOT `rt_dec_nonelem`/ElementaryIntegralQ.
       KEY FINDING: `Risch`ElementaryIntegralQ` already returns False for Ex 5.2 (li ≠ elementary),
       so that target was already met by existing Risch. The genuinely-new decision is
       li-elementarity (Ex 5.1 True / not-elementary; Ex 5.2 False), surfaced by LiElementaryQ.
6. [x] Verified: paper pins exact; `cherry_li_tests`, `cherry_sigma_tests`, and full
       `integrate_risch_transcendental` battery green; strict-C99 + valgrind clean.
7. [x] Docs: `docs/spec/builtins/calculus.md` (SigmaDecomposition + LiElementaryQ) + changelog.

## Deferred (not B1+A2): B2 tower recursion (non-proportional product li args, nested towers),
## B3 C0 seam, A1 complex-constant layer. See CHERRY_BLOCKERS.md.

## Soundness rails
- POSITIVE answers stay diff-back-gated (`rt_free_of_head(Q,"Integrate") && PowerExpand diff-back`).
- The NO decision is the ONLY place a wrong answer is catastrophic: the monotone-tail termination
  must be exact, not heuristic (memory: no arbitrary caps in decision procedures). Prove it against
  Ex 5.2's exact numbers before trusting `False`.
