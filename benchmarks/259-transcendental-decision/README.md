# 259: the transcendental decision corpus

There are 640 indefinite integrals in this benchmark, all class T: transcendental
elementary functions with no radical of x. They are the families of the C
completeness campaign's hunts (`tools/bronstein_hunts/`), rewritten as a shootout
corpus. That puts FriCAS on exactly the integrands the campaign was built against,
which benchmarks 256 and 257 do not cover: Tanh towers, tan over Q(i), Gaussian
exponent offsets, LimitedIntegrate below a logarithm, Liouvillian primitives, and
dependent logarithms and exponentials.

The harness, the verification and the verdicts are those of
[benchmark 256](../256-integration-shootout/README.md); only the corpus differs.
The comparison with FriCAS, which also counts each system's proofs of
non-elementarity, is [`../BRONSTEIN_VS_FRICAS.md`](../BRONSTEIN_VS_FRICAS.md).

```bash
make bench-integrate-decision      # Mathilda (Automatic, RischTranscendental) and FriCAS
make bench-integrate-vs-fricas     # benchmarks/BRONSTEIN_VS_FRICAS.md from 256, 257 and 259
COUNT=40 ./Mathilda -file tools/integration_shootout/gen_campaign.m > benchmarks/259-transcendental-decision/corpus.tsv
```

## The corpus

Each family has two halves of 40 cases:
- `<family>-elem`: `f = D[F]`, with F random in the family's kernels, so f has an
  elementary integral and any "not elementary" is false. The exception is `prim`,
  whose kernels are Erf, Erfi, Ei and li: there F is Liouvillian, and the half is
  labelled `expect = unknown`.
- `<family>-rand`: random integrands built from the same kernels, labelled
  `expect = unknown`; the runs themselves decide them.

| family | hunt | kernels |
|---|---|---|
| `tanh` | `tanhhunt.m` | `Tanh[v]` with exponentials and logs of it; the split special part `t^2 - 1` |
| `tan` | `f8hunt.m` | `Tan[v]` over towers (`E^Tan`, `Log[1 + Tan^2]`, `Tan[u]/u`) |
| `tani` | `tanihunt.m` | `Tan[v]` over Q(i): `Log[Tan - I]`, `1/(Tan + I)`, `E^(I x)` |
| `gauss` | `gausshunt.m` | Gaussian exponent offsets: `E^(I (1 + x))` beside `E^(I x)`, `E^I`, `Tanh[1 + x]` |
| `lim` | `limhunt.m` | `Log[v]` times functions of `E^w`, Tanh, Tan, Coth (LimitedIntegrate below a logarithm) |
| `prim` | `primhunt.m` | Erf, Erfi, Ei, li, as monomials of the tower |
| `dep` | `dephunt.m` | dependent logarithms (`Log[x^2 - 1]` beside `Log[x - 1]`, `Log[x E^x]`) |
| `f12` | `f12hunt.m` | composite exponentials (`E^(x + x^2)` beside `E^x`, `E^(x^2)`) |

The hunts' parameter `a` is `3/2` here, because the points must evaluate. The
evaluation points are as in 256: three rationals where every log argument and
radicand is positive. `f` is real at them, except in `tani` and `gauss`, where
complex values are allowed. Each system verifies its own answer, so the branch
choices of the two systems never meet.

Methods run: Mathilda's Automatic and RischTranscendental (the C implementation of
Bronstein's *Symbolic Integration I*); ParallelMixedTower is left out, and
Automatic runs it anyway. FriCAS's answers with Ei, li, erf or erfi cannot be
evaluated by `complexNumeric`, so they are verified exactly instead:
`normalize(D(r, x) - f) = 0`.

## Current run (2026-11-26)

The integrands were rewritten without their constant factors (`normalize_constants.py`;
a case that became a duplicate was left out) and every system was run
again, Mathilda at v0.6750. The easy labels were re-checked at v0.6760, whose
DerivativeDivides no longer depends on a constant factor: 4 of the run cases are
easy and left out, so 555 are compared. Closed means verified by PossibleZeroQ;
the full tables are in [`REPORT.md`](REPORT.md) and `report/report.pdf`.

| run | closed (% of the cases in the classes the run covers) | wrong |
|---|---|---|
| Mathilda Automatic | 69% (382) | 0 |
| Mathilda RischTranscendental | 62% (343) | 0 |
| FriCAS | 62% (345) | 7 |
