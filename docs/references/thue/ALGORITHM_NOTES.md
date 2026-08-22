# Tzanakis–de Weger Thue algorithm — implementation notes

Distilled from Tzanakis & de Weger, *On the practical solution of the Thue
equation*, J. Number Theory 31 (1989) 99–132 (`tzanakis_deweger_1989.pdf`).
Section/Lemma refs are to that paper. Waldschmidt's bound: Appendix II.

## Setup
`F(X,Y)=Σ f_i X^{n-i} Y^i`, deg n≥3 irreducible; `g(x)=F(x,1)`. Roots
`ξ^(1..s)` real, then `t` complex-conjugate pairs, `s+2t=n`. `K=ℚ(ξ)`.
`β=X−ξY`. Fundamental units `ε_1..ε_r`, `r=s+t−1`. Bounded-norm reps `μ`
with `f_0·N(μ)=m`; for `|f_0|=|m|=1`, `M={+1,−1}`. Then
`β = μ·∏ ε_i^{a_i}`, and we bound/enumerate `A=max_i|a_i|`.

## Constants (Lemmas 1.1, 1.2, 2.1, 2.2)
- `C1 = 2^{n-1}|m| / min_{1≤i≤s}|g'(ξ^(i))|`   (UPPER bound)
- `C2 = ½·min_{i<j}|ξ^(i)−ξ^(j)|`   (LOWER bound)
- `Y0 = ⌈(2^{n-1}|m|/(min_{i≤t}|g'(ξ^{s+i})|·min_{i≤t}|Im ξ^{s+i}|))^{1/n}⌉` (t≥1), else 1
- `Y1 = max(Y0, ⌈(4C1)^{1/(n-2)}⌉)`
- `C3 = max_{i0,j,k distinct} |(ξ^{i0}−ξ^{j})/(ξ^{i0}−ξ^{k})|`   (UPPER)
- `Y2* = max(Y1, ⌈(2 C1 C3/C2)^{1/n}⌉)`
- `U_I = (log|ε_l^{(h_i)}|)_{i,l}` for `I={h_1..h_r}⊆{1..n}`, r×r. `N[U_I^{-1}] = max row-sum |U_I^{-1}|`.
- `μ_ = min_{i,μ∈M}|μ^(i)|`, `μ_+ = max`   (for M={±1}: μ_=μ_+=1)
- `C4 = (½ + max_{i1<i2}|ξ^{i1}−ξ^{i2}|)/μ_`   (UPPER)
- `C5 = min((n−1)·min_I N[U_I^{-1}],  max_I N[U_I^{-1}])`
- `C6 = 1.39·C1·C3·C4^n / C2`
- `Y2' = max(Y2*, 2|m|^{1/n}, μ_+/C2)`
- **Upper side:** `|Λ| < K1·exp(−K2·A)` for `|Y|>Y2'`, with `K1=C6`, `K2=n/C5`.

## Linear form Λ (1.3 real / 1.4 complex)
Fix `i0∈{1..s}` (closest real root), pick `j,k` distinct from `i0` and each other.
- **Real case** (`j,k≤s`): `Λ = log|(ξ^{i0}−ξ^{j})/(ξ^{i0}−ξ^{k})| + Σ_{i=1}^r a_i log|ε_i^{(k)}/ε_i^{(j)}|`. `q=r`, `δ`=first term, `μ_i = log|ε_i^{(k)}/ε_i^{(j)}|`.
- **Complex case** (`j,k` a conjugate pair, `ξ^{(k)}=conj ξ^{(j)}`): `Λ = Arg((ξ^{i0}−ξ^{j})/(ξ^{i0}−ξ^{k})·μ^{k}/μ^{j}) + Σ a_i Arg(ε_i^{(k)}/ε_i^{(j)}) + a_0·2π`. `q=r+1`, `μ_q=2π`, `a_q=a_0`.

## Baker via Waldschmidt (Appendix II)  — N=r+1 logs
`D`=degree of splitting field (safe over-estimate `n!`). For α_j:
- unit `α_j=ε_j^{(k)}/ε_j^{(j)}`: `V_j = log H_j`, `H_j = max_h|ε_j^{(h)}| / min_h|ε_j^{(h)}|`.
- ratio `α_0=(ξ^{i0}−ξ^{j})/(ξ^{i0}−ξ^{k})`: `V_0 = log|disc(g)| + log C3`.
Sort `V_1≤…≤V_N`; `V^+=max(V,1)`; `e(N)=min(8N+51,10N+33,9N+39)`.
`C7 = 2^{e(N)}·N^{2N}·D^{N+2}·(∏V_j)·log(eD·V_{N-1}^+)`, `C8 = log(eD·V_N^+)`,
`C8' = C8` (real) or `C8+log r` (complex). `|Λ| > exp(−C7·(log A + C8'))` if Λ≠0.

**Initial bound** (Lemma 2.4): combine → `A < a + b·log A`, `a=(log C6 + C7 C8')/K2`,
`b=C7/K2`. Solve by downward fixed-point iteration from a huge seed → `K3`.
(Example: K3=3.26×10⁴⁰.)

## Reduction (Prop 3.2, inhomogeneous δ≠0)  — repeat until stable
Given `A<K3`, `|Λ|<K1 exp(−K2 A)`. Pick `c0 ≈ K3^q` (⌈·⌉ with margin;
example c0=10^140 then 10^12). Lattice Γ ⊂ ℤ^q spanned by COLUMNS of q×q `𝒜`:
col `j<q` = `e_j + ⌊c0 μ_j⌉ e_q`; col `q` = `⌊c0 μ_q⌉ e_q`.
LLL-reduce → basis `b_1..b_q` (integer q-vectors), `|b_1|`=first-vector norm.
`x=(0,…,0,−⌊c0 δ⌉)`. Solve `ℬ s = x` over ℚ (ℬ=[b_1|…|b_q] columns). `i*=max{i:s_i∉ℤ}`,
`‖s_{i*}‖`=dist to nearest int. **If** `2^{-(q-1)/2}·‖s_{i*}‖·|b_1| ≥ √(4q²+3q−¾)·K3`,
**then** `A < (1/K2)·log(c0·K1/(q·K3))`. Set K3 := that; repeat with smaller c0.
(Example → A≤72 → A≤10.)

## Completeness enumeration (my direct-exponent framing)
`B = max_{i0∈{1..s}}` reduced bound. Every solution has `|Y|≤Y2'` OR `|Y|>Y2'`:
- **|Y| ≤ Y2'**: brute-force `(X,Y)`, `|Y|≤Y2'`, `|X| ≤ ⌈max_i|ξ^(i)|·Y2' + C1⌉+1`; check `F=m`.
- **|Y| > Y2'**: `A<B` ⇒ enumerate unit exponents `|a_i|≤B` (i=1..r), reconstruct
  `β=±∏ε_i^{a_i}`, keep those in `ℤ+ℤθ`, verify `F(x,y)=m` exactly.
Union + dedup = the complete set. (Rigor: floats only propose; every kept
`(x,y)` is exact-verified, and B/Y2' are provable upper bounds.)

## Validation targets
- Constants reproduce the paper's worked quartic (C1=0.843, C3=6.645, C4=8.337,
  C5=1.211, C6=6.388e4, K3=3.26e40, A≤72, A≤10) — modulo that its order is
  non-monogenic (½ϑ² basis), so our monogenic Gate 1 declines that exact field;
  validate instead on the three monogenic targets by matching the known sets.
- Final `Solve` output must equal the reconstruction-engine sets already verified.
