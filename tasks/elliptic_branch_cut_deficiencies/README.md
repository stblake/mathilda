# EllipticF / EllipticPi on the branch cut

`repro.m` (`./Mathilda -file tasks/elliptic_branch_cut_deficiencies/repro.m`).

At phi = ArcSin[Sqrt[1 + I Sqrt[3]]], m = (1 - I Sqrt[3])/2, the quantity
1 - m Sin[phi]^2 is exactly -1, on the cut of the integrand 1/Sqrt[1 - m Sin[t]^2].

| call | Mathilda | correct (mpmath `ellipf`/`ellippi`, dps 30) |
|---|---|---|
| `EllipticF[N[phi], N[m]]` | 1.86639 + 0.790708 I | 1.8663879592700336655860441508 + 0.790708154959272985388940705058 I |
| `N[EllipticF[phi, m], 30]` | **0.393726134715... + 1.469405939026... I** | as above |
| `EllipticPi[N[n,30], N[phi,30], N[m,30]]`, n = -I/Sqrt[3] | 2.31182980156734... - 0.24024781417467... I | 2.31182980156734432184185434151 - 0.240247814174672231844511193491 I |
| `N[EllipticPi[n, phi, m], 20]`, phi spelled `ArcSin[Sqrt[1 + I Sqrt[3]]]` | 2.3118298015673443218 - 0.24024781417467223184 I | as above |
| the same, phi spelled `ArcSin[Sqrt[(-1 + (-1)^(2/3))/(1/2 + (-1)^(2/3))]]` (equal, by PossibleZeroQ) | **does not return** | as above |

Fix direction:
- **F.** The value on the cut is a convention. The arbitrary-precision path
  should land on the same side as the machine path, mpmath and the reference
  system. That is the limit from the side the principal branches select:
  perturb the midpoint in the direction the machine rounding takes, or
  evaluate `F` by its quasi-periodic reduction off the cut.
- **Pi.** When Arb returns a non-finite ball for an exact on-cut input, N should
  stop escalating the precision. It should decide the side the same way as `F`.

Not hit by Integrate any more: since v0.6526, ParallelMixedSpecial's CubicData
chooses the root ordering whose Legendre form has one sign on the whole real
line, and that form avoids this point.
