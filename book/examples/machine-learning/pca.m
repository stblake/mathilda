# Principal components as an eigenproblem, checked against the SVD; then a reusable reducer.
SeedRandom[4]; cloud = Map[{{2., 0.}, {1.2, 0.5}} . # + {5., 3.} &, Partition[RandomVariate[NormalDistribution[], 400], 2]];
{Mean[cloud], Covariance[cloud]}
{Eigenvalues[Covariance[cloud]], Eigenvectors[Covariance[cloud]]}
pc = PrincipalComponents[cloud]; Take[pc, 3]
{Map[Variance, Transpose[pc]], Chop[Covariance[pc]]}
{u, s, v} = SingularValueDecomposition[Map[# - Mean[cloud] &, cloud]]; Transpose[v]
Diagonal[s]^2/(Length[cloud] - 1)
Take[DimensionReduce[cloud, 1], 3]
r = DimensionReduction[cloud, 1]
{r[{7., 4.2}], r[{{5.05975, 3.03615}, {3., 1.8}}]}
Take[PrincipalComponents[cloud, Method -> "Correlation"], 3]
PrincipalComponents[{{0., 0.}, {1., 1.}, {2., 2.}, {3., 3.}, {4., 4.}}]
