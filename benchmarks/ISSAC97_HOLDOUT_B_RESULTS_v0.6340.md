# ISSAC '97 holdout B (experiments 146-155): classification against the true answers

| family | cases | Mathilda correct | baseline correct | both correct | Mathilda faster | Mathilda slower | only Mathilda correct | only baseline correct | Mathilda wrong | Mathilda no answer |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 ill-conditioned linear algebra | 24 | 24 | 24 | 24 | 22 | 2 | 0 | 0 | 0 | 0 |
| 2 hard definite integrals | 21 | 16 | 16 | 12 | 12 | 0 | 4 | 4 | 2 | 3 |
| 3 slowly convergent series | 20 | 14 | 16 | 14 | 14 | 0 | 0 | 2 | 5 | 1 |
| 4 exact coefficients and residues | 30 | 24 | 27 | 23 | 21 | 2 | 1 | 4 | 1 | 5 |
| 5 zeros of special functions | 25 | 21 | 25 | 21 | 17 | 4 | 0 | 4 | 2 | 2 |
| 6 polynomial systems | 30 | 30 | 28 | 28 | 28 | 0 | 2 | 0 | 0 | 0 |
| 7 cancellation limits | 32 | 14 | 32 | 14 | 13 | 1 | 0 | 18 | 1 | 17 |
| 8 closed-form definite integrals | 30 | 29 | 29 | 29 | 29 | 0 | 0 | 0 | 0 | 1 |
| 9 operator eigenvalues | 22 | 20 | 22 | 20 | 16 | 4 | 0 | 2 | 0 | 2 |
| 10 ODE singularities | 24 | 20 | 23 | 20 | 20 | 0 | 0 | 3 | 2 | 2 |
| **all** | **258** | **212** | **242** | **205** | **192** | **13** | **7** | **37** | **13** | **33** |

## Cases

| case | Mathilda | baseline | Mathilda ms | baseline ms | ratio | speed |
|---|---|---|---|---|---|---|
| b1.1 Eigenvalues machine | CORRECT | CORRECT | 0.1 | 0.2 | 0.35 | faster |
| b1.1 Eigenvalues wp200 | CORRECT | CORRECT | 38.5 | 554.3 | 0.07 | faster |
| b1.2 SingularValueList machine | CORRECT | CORRECT | 0.1 | 0.1 | 0.48 | faster |
| b1.2 SingularValueList wp150 | CORRECT | CORRECT | 31.0 | 86.0 | 0.36 | faster |
| b1.3 LinearSolve exact | CORRECT | CORRECT | 2.5 | 31.2 | 0.08 | faster |
| b1.3 LinearSolve wp100 | CORRECT | CORRECT | 3.1 | 20.5 | 0.15 | faster |
| b1.3 LinearSolve machine | CORRECT | CORRECT | 0.2 | 0.5 | 0.33 | faster |
| b1.4 SingularValueList machine | CORRECT | CORRECT | 0.2 | 0.8 | 0.28 | faster |
| b1.4 SingularValueList wp100 | CORRECT | CORRECT | 293.0 | 1143.6 | 0.26 | faster |
| b1.5 Eigenvalues machine | CORRECT | CORRECT | 0.1 | 0.3 | 0.40 | faster |
| b1.5 Eigenvalues wp150 | CORRECT | CORRECT | 70.2 | 386.7 | 0.18 | faster |
| b1.6 LeastSquares exact | CORRECT | CORRECT | 272.9 | 48.6 | 5.62 | slower |
| b1.6 LeastSquares wp60 | CORRECT | CORRECT | 160.8 | 144.9 | 1.11 | slower |
| b1.6 LeastSquares machine | CORRECT | CORRECT | 0.4 | 1.2 | 0.36 | faster |
| b1.7 Det exact | CORRECT | CORRECT | 6.5 | 16.7 | 0.39 | faster |
| b1.7 Det wp120 | CORRECT | CORRECT | 1.6 | 6.5 | 0.25 | faster |
| b1.7 Det machine | CORRECT | CORRECT | 0.1 | 0.3 | 0.33 | faster |
| b1.8 Eigenvalues machine | CORRECT | CORRECT | 0.6 | 1.7 | 0.35 | faster |
| b1.8 Eigenvalues wp120 | CORRECT | CORRECT | 15.7 | 172.8 | 0.09 | faster |
| b1.9 Eigenvalues machine | CORRECT | CORRECT | 0.1 | 0.5 | 0.27 | faster |
| b1.9 Eigenvalues wp60 | CORRECT | CORRECT | 5.9 | 130.8 | 0.04 | faster |
| b1.10 Det exact | CORRECT | CORRECT | 3.9 | 9.7 | 0.40 | faster |
| b1.10 Det wp300 | CORRECT | CORRECT | 2.7 | 11.8 | 0.23 | faster |
| b1.10 Det machine | CORRECT | CORRECT | 0.0 | 0.1 | 0.45 | faster |
| b2.1 NIntegrate machine | CORRECT | CORRECT | 0.1 | 44.2 | 0.00 | faster |
| b2.1 NIntegrate wp30 | CORRECT | CORRECT | 4.2 | 715.8 | 0.01 | faster |
| b2.2 NIntegrate machine | CORRECT | WRONG | 7.9 |  |  |  |
| b2.2 NIntegrate wp30 | CORRECT | WRONG | 10.1 |  |  |  |
| b2.3 NIntegrate machine | CORRECT | CORRECT | 0.1 | 17.5 | 0.00 | faster |
| b2.3 NIntegrate wp30 | CORRECT | CORRECT | 2.9 | 153.9 | 0.02 | faster |
| b2.4 NIntegrate machine | CORRECT | CORRECT | 0.1 | 13.1 | 0.01 | faster |
| b2.4 NIntegrate wp30 | CORRECT | CORRECT | 4.8 | 48.7 | 0.10 | faster |
| b2.6 NIntegrate machine | CORRECT | CORRECT | 0.0 | 14.3 | 0.00 | faster |
| b2.6 NIntegrate wp30 | WRONG | CORRECT |  | 22.5 |  |  |
| b2.8 NIntegrate machine | CORRECT | CORRECT | 0.2 | 221.9 | 0.00 | faster |
| b2.8 NIntegrate wp30 | WRONG | CORRECT |  | 588.1 |  |  |
| b2.9 NIntegrate machine | CORRECT | WRONG | 10.9 |  |  |  |
| b2.9 NIntegrate wp30 | CORRECT | WRONG | 568.0 |  |  |  |
| b2.10 NIntegrate machine | CORRECT | CORRECT | 43.1 | 406.9 | 0.11 | faster |
| b2.10 NIntegrate wp30 | NO ANSWER | CORRECT |  | 10961.0 |  |  |
| b2.5 NIntegrate machine | CORRECT | CORRECT | 1.0 | 39.0 | 0.02 | faster |
| b2.5 NIntegrate machine X=25 | CORRECT | CORRECT | 0.0 | 20.7 | 0.00 | faster |
| b2.5 NIntegrate wp30 | CORRECT | CORRECT | 2.5 | 69.1 | 0.04 | faster |
| b2.7 NIntegrate machine | NO ANSWER | CORRECT |  | 1219.0 |  |  |
| b2.7 NIntegrate wp20 | NO ANSWER | NO ANSWER |  |  |  |  |
| b3.1 NSum machine | CORRECT | CORRECT | 0.1 | 91.9 | 0.00 | faster |
| b3.1 NSum wp30 | CORRECT | CORRECT | 0.7 | 123.5 | 0.01 | faster |
| b3.2 NSum machine | CORRECT | CORRECT | 0.9 | 27.7 | 0.03 | faster |
| b3.2 NSum wp30 | CORRECT | CORRECT | 38.4 | 101.0 | 0.38 | faster |
| b3.3 NSum machine | CORRECT | CORRECT | 2.6 | 58.4 | 0.04 | faster |
| b3.3 NSum wp30 | CORRECT | CORRECT | 26.0 | 99.8 | 0.26 | faster |
| b3.4 NSum machine | CORRECT | CORRECT | 0.2 | 33.6 | 0.01 | faster |
| b3.4 NSum wp30 | CORRECT | CORRECT | 6.7 | 65.8 | 0.10 | faster |
| b3.5 NSum machine | CORRECT | CORRECT | 1.5 | 246.8 | 0.01 | faster |
| b3.5 NSum wp30 | CORRECT | CORRECT | 55.7 | 632.3 | 0.09 | faster |
| b3.6 NSum machine | CORRECT | CORRECT | 0.2 | 89.8 | 0.00 | faster |
| b3.6 NSum wp30 | CORRECT | CORRECT | 0.3 | 88.1 | 0.00 | faster |
| b3.7 NSum machine | CORRECT | CORRECT | 0.2 | 99.9 | 0.00 | faster |
| b3.7 NSum wp30 | CORRECT | CORRECT | 0.9 | 105.5 | 0.01 | faster |
| b3.8 NSum machine | WRONG | WRONG |  |  |  |  |
| b3.8 NSum wp30 | WRONG | WRONG |  |  |  |  |
| b3.9 NSum machine | WRONG | WRONG |  |  |  |  |
| b3.9 NSum wp30 | WRONG | WRONG |  |  |  |  |
| b3.10 NSum machine | WRONG | CORRECT |  | 37904.7 |  |  |
| b3.10 NSum wp30 | NO ANSWER | CORRECT |  | 43368.9 |  |  |
| b4.1 SeriesCoefficient exact | CORRECT | CORRECT | 932.6 | 34800.2 | 0.03 | faster |
| b4.1 contour wp60 | CORRECT | NO ANSWER | 115757.0 |  |  |  |
| b4.1 contour machine n=100 | CORRECT | CORRECT | 1173.9 | 4654.1 | 0.25 | faster |
| b4.2 SeriesCoefficient exact | CORRECT | CORRECT | 67.3 | 3313.5 | 0.02 | faster |
| b4.2 contour wp160 | CORRECT | CORRECT | 66.0 | 649.3 | 0.10 | faster |
| b4.2 contour machine n=20 | CORRECT | CORRECT | 7.1 | 143.7 | 0.05 | faster |
| b4.3 Residue exact | CORRECT | CORRECT | 24.2 | 1168.9 | 0.02 | faster |
| b4.3 contour wp40 | CORRECT | CORRECT | 434.3 | 747.0 | 0.58 | faster |
| b4.3 contour machine | WRONG | CORRECT |  | 344.4 |  |  |
| b4.4 SeriesCoefficient exact | CORRECT | CORRECT | 45.0 | 3098.8 | 0.01 | faster |
| b4.4 contour wp80 | CORRECT | CORRECT | 84830.7 | 22417.0 | 3.78 | slower |
| b4.4 contour machine n=100 | CORRECT | CORRECT | 128.5 | 560.8 | 0.23 | faster |
| b4.5 SeriesCoefficient exact | CORRECT | CORRECT | 7.2 | 85.5 | 0.08 | faster |
| b4.5 contour wp100 | CORRECT | CORRECT | 606.7 | 767.9 | 0.79 | faster |
| b4.5 contour machine n=30 | CORRECT | CORRECT | 39.8 | 120.9 | 0.33 | faster |
| b4.6 SeriesCoefficient exact | CORRECT | CORRECT | 760.2 | 36328.8 | 0.02 | faster |
| b4.6 contour wp40 | NO ANSWER | NO ANSWER |  |  |  |  |
| b4.6 contour machine | NO ANSWER | CORRECT |  | 57545.8 |  |  |
| b4.7 Residue exact | NO ANSWER | NO ANSWER |  |  |  |  |
| b4.7 contour wp50 | CORRECT | CORRECT | 15.2 | 62.6 | 0.24 | faster |
| b4.7 contour machine z^3 | CORRECT | CORRECT | 2.7 | 44.0 | 0.06 | faster |
| b4.8 Coefficient exact | CORRECT | CORRECT | 4.5 | 129.0 | 0.03 | faster |
| b4.8 contour wp30 | CORRECT | CORRECT | 1664.7 | 318.4 | 5.23 | slower |
| b4.8 contour machine n=40 | CORRECT | CORRECT | 2.8 | 78.3 | 0.04 | faster |
| b4.9 SeriesCoefficient exact | CORRECT | CORRECT | 873.8 | 11677.9 | 0.07 | faster |
| b4.9 contour wp40 | CORRECT | CORRECT | 273.8 | 848.5 | 0.32 | faster |
| b4.9 contour machine n=30 | CORRECT | CORRECT | 7.5 | 295.4 | 0.03 | faster |
| b4.10 Residue exact | NO ANSWER | CORRECT |  | 160.1 |  |  |
| b4.10 contour wp40 | NO ANSWER | CORRECT |  | 97.0 |  |  |
| b4.10 contour machine z^10 | CORRECT | CORRECT | 7.4 | 121.9 | 0.06 | faster |
| b5.1 FindRoot wp60 | CORRECT | CORRECT | 118.7 | 213.5 | 0.56 | faster |
| b5.1 FindRoot machine n=100 | CORRECT | CORRECT | 0.1 | 3.5 | 0.02 | faster |
| b5.2 FindRoot wp40 | CORRECT | CORRECT | 2.6 | 85.5 | 0.03 | faster |
| b5.2 FindRoot machine | CORRECT | CORRECT | 0.1 | 0.3 | 0.23 | faster |
| b5.3 FindRoot wp40 | CORRECT | CORRECT | 0.2 | 1.6 | 0.14 | faster |
| b5.3 FindRoot machine | CORRECT | CORRECT | 0.0 | 0.6 | 0.05 | faster |
| b5.5 FindRoot wp40 | CORRECT | CORRECT | 0.2 | 1.1 | 0.20 | faster |
| b5.5 FindRoot machine | CORRECT | CORRECT | 0.0 | 0.3 | 0.06 | faster |
| b5.6 FindRoot wp40 | CORRECT | CORRECT | 4.6 | 2.5 | 1.85 | slower |
| b5.6 FindRoot machine | CORRECT | CORRECT | 0.2 | 0.6 | 0.36 | faster |
| b5.8 FindRoot wp40 | CORRECT | CORRECT | 8.7 | 63.0 | 0.14 | faster |
| b5.8 FindRoot machine | CORRECT | CORRECT | 0.1 | 9.7 | 0.01 | faster |
| b5.3 BesselJZero wp30 | CORRECT | CORRECT | 0.9 | 0.0 | 19.60 | slower |
| b5.4 FindRoot wp30 | CORRECT | CORRECT | 8.7 | 208.6 | 0.04 | faster |
| b5.4 FindRoot machine | CORRECT | CORRECT | 3.6 | 1.8 | 2.02 | slower |
| b5.4 ZetaZero wp30 | NO ANSWER | CORRECT |  | 0.0 |  |  |
| b5.5 AiryAiZero wp30 | NO ANSWER | CORRECT |  | 0.0 |  |  |
| b5.7 CountRoots exact | WRONG | CORRECT |  | 3.8 |  |  |
| b5.7 NSolve wp40 | CORRECT | CORRECT | 1.4 | 10.0 | 0.14 | faster |
| b5.7 FindRoot machine | CORRECT | CORRECT | 1.9 | 1.2 | 1.61 | slower |
| b5.9 CountRoots exact | WRONG | CORRECT |  | 0.9 |  |  |
| b5.9 NSolve wp80 | CORRECT | CORRECT | 5.6 | 7.6 | 0.73 | faster |
| b5.9 NSolve machine n=8 | CORRECT | CORRECT | 0.1 | 1.9 | 0.08 | faster |
| b5.10 FindRoot wp30 | CORRECT | CORRECT | 1.3 | 11.5 | 0.11 | faster |
| b5.10 FindRoot machine | CORRECT | CORRECT | 0.0 | 3.7 | 0.01 | faster |
| b6.1 GroebnerBasis exact | CORRECT | CORRECT | 18.1 | 1558.7 | 0.01 | faster |
| b6.1 NSolve wp40 | CORRECT | CORRECT | 33.5 | 3807.9 | 0.01 | faster |
| b6.1 NSolve machine | CORRECT | WRONG | 14.7 |  |  |  |
| b6.2 GroebnerBasis exact | CORRECT | CORRECT | 4.1 | 176.8 | 0.02 | faster |
| b6.2 NSolve wp40 | CORRECT | CORRECT | 16.0 | 1857.4 | 0.01 | faster |
| b6.2 NSolve machine | CORRECT | CORRECT | 6.6 | 995.8 | 0.01 | faster |
| b6.3 GroebnerBasis exact | CORRECT | CORRECT | 14.6 | 53.1 | 0.27 | faster |
| b6.3 NSolve wp40 | CORRECT | CORRECT | 12.6 | 425.6 | 0.03 | faster |
| b6.3 NSolve machine | CORRECT | CORRECT | 9.8 | 299.2 | 0.03 | faster |
| b6.4 GroebnerBasis exact | CORRECT | CORRECT | 13.3 | 73.6 | 0.18 | faster |
| b6.4 NSolve wp40 | CORRECT | CORRECT | 14.2 | 1144.3 | 0.01 | faster |
| b6.4 NSolve machine | CORRECT | CORRECT | 5.5 | 76.2 | 0.07 | faster |
| b6.7 GroebnerBasis exact | CORRECT | CORRECT | 93.1 | 719.3 | 0.13 | faster |
| b6.7 NSolve wp40 | CORRECT | CORRECT | 115.1 | 3752.4 | 0.03 | faster |
| b6.7 NSolve machine | CORRECT | CORRECT | 119.0 | 2064.9 | 0.06 | faster |
| b6.8 GroebnerBasis exact | CORRECT | CORRECT | 0.3 | 1.6 | 0.19 | faster |
| b6.8 NSolve wp40 | CORRECT | CORRECT | 6.4 | 360.8 | 0.02 | faster |
| b6.8 NSolve machine | CORRECT | CORRECT | 1.3 | 31.7 | 0.04 | faster |
| b6.9 GroebnerBasis exact | CORRECT | CORRECT | 0.2 | 1.9 | 0.11 | faster |
| b6.9 NSolve wp40 | CORRECT | CORRECT | 42.4 | 149.9 | 0.28 | faster |
| b6.9 NSolve machine | CORRECT | CORRECT | 17.1 | 83.3 | 0.21 | faster |
| b6.10 GroebnerBasis exact | CORRECT | CORRECT | 17.4 | 30.3 | 0.57 | faster |
| b6.10 NSolve wp40 | CORRECT | CORRECT | 17.1 | 162.0 | 0.11 | faster |
| b6.10 NSolve machine | CORRECT | WRONG | 11.6 |  |  |  |
| b6.5 GroebnerBasis exact | CORRECT | CORRECT | 0.3 | 1.9 | 0.17 | faster |
| b6.5 NSolve wp40 | CORRECT | CORRECT | 0.5 | 20.6 | 0.02 | faster |
| b6.5 NSolve machine | CORRECT | CORRECT | 0.5 | 10.4 | 0.05 | faster |
| b6.6 GroebnerBasis eliminate exact | CORRECT | CORRECT | 0.3 | 2.0 | 0.16 | faster |
| b6.6 NSolve wp40 | CORRECT | CORRECT | 2.8 | 125.6 | 0.02 | faster |
| b6.6 NSolve machine | CORRECT | CORRECT | 1.8 | 101.7 | 0.02 | faster |
| b7.1 Limit exact | WRONG | CORRECT |  | 242.8 |  |  |
| b7.1 approach wp300 | CORRECT | CORRECT | 2.7 | 3.8 | 0.71 | faster |
| b7.1 Neville machine | CORRECT | CORRECT | 0.1 | 0.6 | 0.21 | faster |
| b7.2 Limit exact | CORRECT | CORRECT | 114.1 | 191.7 | 0.59 | faster |
| b7.2 Neville wp600 | CORRECT | CORRECT | 2.2 | 1.7 | 1.27 | slower |
| b7.2 Neville machine | CORRECT | CORRECT | 0.0 | 0.1 | 0.60 | faster |
| b7.3 Limit exact | CORRECT | CORRECT | 2.2 | 208.4 | 0.01 | faster |
| b7.3 approach wp200 | CORRECT | CORRECT | 0.3 | 0.5 | 0.54 | faster |
| b7.3 Neville machine | CORRECT | CORRECT | 0.1 | 0.1 | 0.57 | faster |
| b7.4 Limit exact | NO ANSWER | CORRECT |  | 50.5 |  |  |
| b7.4 approach wp100 | CORRECT | CORRECT | 1.4 | 5.3 | 0.26 | faster |
| b7.4 Neville machine | CORRECT | CORRECT | 0.1 | 0.3 | 0.30 | faster |
| b7.5 Limit exact | CORRECT | CORRECT | 8.1 | 91.2 | 0.09 | faster |
| b7.5 approach wp150 | CORRECT | CORRECT | 0.5 | 0.6 | 0.80 | faster |
| b7.5 Neville machine | CORRECT | CORRECT | 0.1 | 0.2 | 0.57 | faster |
| b7.6 Limit exact | NO ANSWER | CORRECT |  | 67.7 |  |  |
| b7.6 approach wp150 | NO ANSWER | CORRECT |  | 0.9 |  |  |
| b7.6 sample machine x=1e-8 | NO ANSWER | CORRECT |  | 0.0 |  |  |
| b7.7 Limit exact | NO ANSWER | CORRECT |  | 102.7 |  |  |
| b7.7 approach wp150 | NO ANSWER | CORRECT |  | 0.8 |  |  |
| b7.7 Neville machine | NO ANSWER | CORRECT |  | 0.2 |  |  |
| b7.8 Limit exact | NO ANSWER | CORRECT |  | 77.3 |  |  |
| b7.8 approach wp250 | NO ANSWER | CORRECT |  | 1.2 |  |  |
| b7.8 Neville machine | NO ANSWER | CORRECT |  | 0.2 |  |  |
| b7.9 Limit exact | NO ANSWER | CORRECT |  | 93.8 |  |  |
| b7.9 Neville wp100 | NO ANSWER | CORRECT |  | 1.7 |  |  |
| b7.9 Neville machine | NO ANSWER | CORRECT |  | 0.2 |  |  |
| b7.10 Limit exact | NO ANSWER | CORRECT |  | 164.0 |  |  |
| b7.10 approach wp100 | NO ANSWER | CORRECT |  | 0.8 |  |  |
| b7.10 Neville machine | NO ANSWER | CORRECT |  | 0.3 |  |  |
| b7.3 Series exact | CORRECT | CORRECT | 0.7 | 17.7 | 0.04 | faster |
| b7.8 Series exact | NO ANSWER | CORRECT |  | 7.4 |  |  |
| b8.1 Integrate exact | CORRECT | CORRECT | 97.0 | 2242.2 | 0.04 | faster |
| b8.1 NIntegrate wp40 | CORRECT | CORRECT | 3.0 | 13.7 | 0.22 | faster |
| b8.1 NIntegrate machine | CORRECT | CORRECT | 0.0 | 7.4 | 0.00 | faster |
| b8.2 Integrate exact | CORRECT | CORRECT | 47.4 | 4645.6 | 0.01 | faster |
| b8.2 NIntegrate wp40 | CORRECT | CORRECT | 4.1 | 54.1 | 0.08 | faster |
| b8.2 NIntegrate machine | CORRECT | CORRECT | 0.1 | 9.8 | 0.01 | faster |
| b8.3 Integrate exact | CORRECT | CORRECT | 64.7 | 1841.2 | 0.04 | faster |
| b8.3 NIntegrate wp40 | CORRECT | CORRECT | 2.7 | 11.2 | 0.25 | faster |
| b8.3 NIntegrate machine | CORRECT | CORRECT | 0.0 | 4.9 | 0.00 | faster |
| b8.4 Integrate exact | CORRECT | CORRECT | 7.6 | 3611.3 | 0.00 | faster |
| b8.4 NIntegrate wp40 | CORRECT | CORRECT | 15.8 | 60.6 | 0.26 | faster |
| b8.4 NIntegrate machine | CORRECT | CORRECT | 0.0 | 10.0 | 0.00 | faster |
| b8.5 Integrate exact | CORRECT | CORRECT | 4.8 | 14974.5 | 0.00 | faster |
| b8.5 NIntegrate wp40 | CORRECT | CORRECT | 11.4 | 1732.6 | 0.01 | faster |
| b8.5 NIntegrate machine | CORRECT | CORRECT | 0.1 | 30.7 | 0.00 | faster |
| b8.6 Integrate exact | CORRECT | CORRECT | 23.2 | 4078.2 | 0.01 | faster |
| b8.6 NIntegrate wp40 | CORRECT | CORRECT | 5.5 | 18.8 | 0.29 | faster |
| b8.6 NIntegrate machine | CORRECT | CORRECT | 0.0 | 7.1 | 0.00 | faster |
| b8.7 Integrate exact | CORRECT | CORRECT | 53.7 | 2409.3 | 0.02 | faster |
| b8.7 NIntegrate wp40 | CORRECT | CORRECT | 4.1 | 39.2 | 0.10 | faster |
| b8.7 NIntegrate machine | CORRECT | CORRECT | 0.1 | 8.7 | 0.01 | faster |
| b8.8 Integrate exact | CORRECT | CORRECT | 192.5 | 8299.6 | 0.02 | faster |
| b8.8 NIntegrate wp40 | CORRECT | CORRECT | 5.3 | 13.0 | 0.41 | faster |
| b8.8 NIntegrate machine | CORRECT | CORRECT | 0.0 | 6.5 | 0.00 | faster |
| b8.9 Integrate exact | CORRECT | CORRECT | 0.0 | 2643.4 | 0.00 | faster |
| b8.9 NIntegrate wp40 | CORRECT | CORRECT | 14.4 | 42.9 | 0.34 | faster |
| b8.9 NIntegrate machine | CORRECT | CORRECT | 0.0 | 8.2 | 0.00 | faster |
| b8.10 Integrate exact | NO ANSWER | NO ANSWER |  |  |  |  |
| b8.10 NIntegrate wp40 | CORRECT | CORRECT | 109.9 | 427.4 | 0.26 | faster |
| b8.10 NIntegrate machine | CORRECT | CORRECT | 0.1 | 18.9 | 0.00 | faster |
| b9.1 secular FindRoot wp40 | CORRECT | CORRECT | 0.1 | 0.6 | 0.25 | faster |
| b9.1 Nystrom Richardson machine | CORRECT | CORRECT | 541.3 | 991.0 | 0.55 | faster |
| b9.2 Gauss-Hermite Nystrom wp30 | CORRECT | CORRECT | 81.1 | 372.3 | 0.22 | faster |
| b9.2 Gauss-Hermite Nystrom machine | CORRECT | CORRECT | 14.5 | 47.3 | 0.31 | faster |
| b9.3 trapezoid Nystrom wp30 N=256 | CORRECT | CORRECT | 2225.3 | 17451.6 | 0.13 | faster |
| b9.3 trapezoid Nystrom machine N=256 | CORRECT | CORRECT | 356.5 | 605.1 | 0.59 | faster |
| b9.4 Eigenvalues wp50 N=150 | CORRECT | CORRECT | 187.9 | 4667.5 | 0.04 | faster |
| b9.4 Eigenvalues machine N=500 | CORRECT | CORRECT | 80.5 | 12.6 | 6.38 | slower |
| b9.5 Eigenvalues wp60 n=60 | CORRECT | CORRECT | 50.2 | 432.0 | 0.12 | faster |
| b9.5 Eigenvalues machine n=20 | CORRECT | CORRECT | 0.4 | 0.5 | 0.81 | faster |
| b9.6 Fourier matrix wp40 | CORRECT | CORRECT | 4.6 | 78.5 | 0.06 | faster |
| b9.6 Fourier matrix machine | CORRECT | CORRECT | 1.0 | 1.4 | 0.72 | faster |
| b9.6 MathieuCharacteristicA wp30 | NO ANSWER | CORRECT |  | 0.0 |  |  |
| b9.6 MathieuCharacteristicA machine | NO ANSWER | CORRECT |  | 0.0 |  |  |
| b9.7 harmonic basis wp50 N=100 | CORRECT | CORRECT | 3154.6 | 1653.1 | 1.91 | slower |
| b9.7 harmonic basis machine N=60 | CORRECT | CORRECT | 11.0 | 7.3 | 1.50 | slower |
| b9.8 Gauss-Legendre Nystrom wp40 | CORRECT | CORRECT | 31.5 | 141.5 | 0.22 | faster |
| b9.8 Gauss-Legendre Nystrom machine | CORRECT | CORRECT | 5.9 | 7.7 | 0.76 | faster |
| b9.9 Eigenvalues wp40 N=10 | CORRECT | CORRECT | 75.0 | 638.7 | 0.12 | faster |
| b9.9 Eigenvalues machine N=40 | CORRECT | CORRECT | 313.6 | 178.2 | 1.76 | slower |
| b9.10 Eigenvalues pencil wp40 N=60 | CORRECT | CORRECT | 26.8 | 718.5 | 0.04 | faster |
| b9.10 Eigenvalues pencil machine N=400 | CORRECT | CORRECT | 370.3 | 470.6 | 0.79 | faster |
| b10.1 naive domain end machine | CORRECT | CORRECT | 0.4 | 2.3 | 0.19 | faster |
| b10.1 reciprocal switch wp40 | CORRECT | CORRECT | 2.0 | 75.0 | 0.03 | faster |
| b10.1 reciprocal switch machine | CORRECT | CORRECT | 0.2 | 1.9 | 0.13 | faster |
| b10.2 EllipticK exact | CORRECT | CORRECT | 0.0 | 0.0 | 0.72 | faster |
| b10.2 NDSolve quarter period wp40 | CORRECT | CORRECT | 2.1 | 65.9 | 0.03 | faster |
| b10.2 NDSolve quarter period machine | CORRECT | CORRECT | 0.2 | 1.1 | 0.20 | faster |
| b10.3 NDSolve x(20) wp40 | CORRECT | CORRECT | 47.0 | 1802.2 | 0.03 | faster |
| b10.3 NDSolve x(5) machine | CORRECT | CORRECT | 0.9 | 1.9 | 0.48 | faster |
| b10.4 series start + FindRoot wp30 | CORRECT | CORRECT | 2.4 | 23.3 | 0.10 | faster |
| b10.4 series start + FindRoot machine | CORRECT | CORRECT | 0.3 | 1.3 | 0.22 | faster |
| b10.5 Toepfer wp30 | CORRECT | CORRECT | 0.8 | 31.2 | 0.02 | faster |
| b10.5 Toepfer machine | CORRECT | CORRECT | 0.3 | 1.0 | 0.32 | faster |
| b10.6 NDSolve y(8) wp50 | WRONG | CORRECT |  | 1173.2 |  |  |
| b10.6 NDSolve y(3) machine | CORRECT | CORRECT | 0.1 | 0.8 | 0.15 | faster |
| b10.7 naive domain end machine | CORRECT | CORRECT | 0.8 | 2.9 | 0.26 | faster |
| b10.7 energy tail wp30 | CORRECT | CORRECT | 1.4 | 67.2 | 0.02 | faster |
| b10.7 energy tail machine | CORRECT | CORRECT | 0.3 | 1.2 | 0.21 | faster |
| b10.8 bisection shooting wp30 | WRONG | WRONG |  |  |  |  |
| b10.8 bisection shooting machine | CORRECT | CORRECT | 17.1 | 109.5 | 0.16 | faster |
| b10.9 naive domain end machine | CORRECT | CORRECT | 0.5 | 2.1 | 0.22 | faster |
| b10.9 inverse function wp30 | CORRECT | CORRECT | 3.6 | 88.5 | 0.04 | faster |
| b10.9 inverse function machine | CORRECT | CORRECT | 0.5 | 1.8 | 0.26 | faster |
| b10.10 NDSolve e=0.99 wp30 | NO ANSWER | CORRECT |  | 484.2 |  |  |
| b10.10 NDSolve e=0.9 machine | NO ANSWER | CORRECT |  | 3.1 |  |  |

## Mathilda wrong answers

| case | Mathilda check | true check |
|---|---|---|
| b2.6 NIntegrate wp30 | `{-188488201359, 2}` | `{599620322995, -1}` |
| b2.8 NIntegrate wp30 | `{360135870882, 0}` | `{360133023322, 0}` |
| b3.8 NSum machine | `{264210, 0}` | `{270581, 0}` |
| b3.8 NSum wp30 | `{264210472138, 0}` | `{270580808428, 0}` |
| b3.9 NSum machine | `{136335, 0}` | `{136843, 0}` |
| b3.9 NSum wp30 | `{136436911026, 0}` | `{136843277762, 0}` |
| b3.10 NSum machine | `{33480, -1}` | `{50000, -1}` |
| b4.3 contour machine | `{3174, -54}` | `{3189, -54}` |
| b5.7 CountRoots exact | `CountRoots[(1 - x)*(2 - x)*(3 - x)*(4 - x)*(5 - x)*(6 - x)*(7 - x)*(8 - x)*(9 - x)*(10 - x)*(11 - x)*(12 - x)*(13 - x)*(14 - x)*(15 - x)*(16 - x)*(17 - x)*(18 - x)*(19 - x)*(20 - x) - x^19/8388608, x]` | `10` |
| b5.9 CountRoots exact | `CountRoots[x^40 - 2*(-1 + 10*x)^2, {x, 0, 1/5}]` | `2` |
| b7.1 Limit exact | `{Round[10^(14 - Floor[Log[Abs[ComplexInfinity]]/Log[10]])*Abs[ComplexInfinity]]*Sign[ComplexInfinity], Floor[Log[Abs[ComplexInfinity]]/Log[10]]}` | `{270580808427785, -1}` |
| b10.6 NDSolve y(8) wp50 | `{46922036, -8}` | `{46922076, -8}` |
| b10.8 bisection shooting wp30 | `{-158807260, 0}` | `{-158807102, 0}` |

## Baseline column: cases not answered correctly

| case | baseline | baseline check | true check |
|---|---|---|---|
| b2.2 NIntegrate machine | WRONG | `{539221, 2}` | `{999760, 2}` |
| b2.2 NIntegrate wp30 | WRONG | `{999755557848, 2}` | `{999760310630, 2}` |
| b2.9 NIntegrate machine | WRONG | `{48082, -3}` | `{62513, -3}` |
| b2.9 NIntegrate wp30 | WRONG | `{480819728211, -3}` | `{625129234764, -3}` |
| b2.7 NIntegrate wp20 | NO ANSWER | `"timeout"` | `{156725, 1}` |
| b3.8 NSum machine | WRONG | `{248120, 0}` | `{270581, 0}` |
| b3.8 NSum wp30 | WRONG | `{248120186455, 0}` | `{270580808428, 0}` |
| b3.9 NSum machine | WRONG | `{134158, 0}` | `{136843, 0}` |
| b3.9 NSum wp30 | WRONG | `{134157618212, 0}` | `{136843277762, 0}` |
| b4.1 contour wp60 | NO ANSWER | `"timeout"` | `{240614678640326, 31}` |
| b4.6 contour wp40 | NO ANSWER | `"timeout"` | `{-303284129702400, 16}` |
| b4.7 Residue exact | NO ANSWER | `"unevaluated"` | `{204822323834388, -20}` |
| b6.1 NSolve machine | WRONG | `{30, {33333333, -1}}` | `{32, {33333333, -1}}` |
| b6.10 NSolve machine | WRONG | `{9, {-67290476, -2}}` | `{10, {-67290476, -2}}` |
| b8.10 Integrate exact | NO ANSWER | `"unevaluated"` | `{384373946213433, -1}` |
| b10.8 bisection shooting wp30 | WRONG | `{-158807260, 0}` | `{-158807102, 0}` |
