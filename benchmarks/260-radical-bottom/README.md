# 260: a radical at the bottom

There are 420 indefinite integrals in this benchmark, all class B: elementary
functions in `Q(x)(y_1..y_m)(t_1..t_k)`, with the `y_i` simple radicals of x and the
`t_i` logarithms and exponentials of elements of `Q(x)` and of `Q(x)(y)` (a radical
inside a kernel). Class B was the one class of benchmark 256 where FriCAS led at
v0.6646 (93% decided against Mathilda's 72%), and 256's 169 class-B cases come from
six families, and its perturbed families have no ground truth. This corpus supplies
it (Campaign D0 in `src/internal/bronstein/BRONSTEIN_PLAN.md`).

The harness, verification and verdicts are those of
[benchmark 256](../256-integration-shootout/README.md); the comparison with FriCAS
is [`../BRONSTEIN_VS_FRICAS.md`](../BRONSTEIN_VS_FRICAS.md).

```bash
make bench-integrate-bottom        # Mathilda (Automatic, ParallelMixedTower, Bronstein) and FriCAS
make bench-integrate-vs-fricas     # the comparison, now with 260
COUNT=60 ./Mathilda -file tools/integration_shootout/gen_bottom.m > benchmarks/260-radical-bottom/corpus.tsv
```

## The corpus

| family | field | halves |
|---|---|---|
| `one` | one radical, kernels over Q(x) (`E^x`, `E^(x^2)`, `Log[x]`, `Log[1 + x]`) | `-elem`, `-rand` |
| `mid` | one radical, also inside kernels (`Log[x + y]`, `Log[1 + y]`, `E^y`, `E^(x y)`), with further monomials | `-elem`, `-rand` |
| `two` | two independent radicals, a term over each and a mixed one | `-elem`, `-rand` |
| `pert` | `g + D[F]`: `g` non-elementary over `Sqrt[x]`, `F` over a second radical | `-nonelem` |

The radicals are square roots of quadratics and of x (genus 0), of cubics and
quartics (genus 1 and more), and a cube root.

- `-elem`: `f = D[F]`, F random in the family's field; the known F is recorded.
  Any "not elementary" is a false certificate.
- `-rand`: random integrands of the same field, `expect = unknown`.
- `pert-nonelem`: `g` is `E^(k x)/Sqrt[x + c]` or `E^(k x) Sqrt[x + c]` (the erfi
  type) or `Sqrt[x]/Log[x]` or `1/(Sqrt[x] Log[x])` (the li type), each classically
  non-elementary, plus the derivative of an elementary F over an independent radical.
  So f is not elementary, and an answer from any system is wrong.

Points as in 256: three rationals where every radicand and log argument is positive
and f is real.

## Current run (2026-11-26)

The integrands were rewritten without their constant factors (`normalize_constants.py`;
a case that became a duplicate was left out) and every system was run
again, Mathilda at v0.6750. The easy labels were re-checked at v0.6760, whose
DerivativeDivides no longer depends on a constant factor: 2 of the run cases are
easy and left out, so 269 are compared. Closed means verified by PossibleZeroQ;
the full tables are in [`REPORT.md`](REPORT.md) and `report/report.pdf`.

| run | closed (% of the cases in the classes the run covers) | wrong |
|---|---|---|
| Mathilda Automatic | 46% (124) | 0 |
| Mathilda ParallelMixedTower | 26% (71) | 0 |
| Mathilda Bronstein | 41% (111) | 0 |
| FriCAS | 48% (129) | 0 |
