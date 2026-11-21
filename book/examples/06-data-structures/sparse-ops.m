# A 10^6 x 10^6 second-difference matrix, and what runs on it natively
n = 10^6;
L = SparseArray[{Band[{1, 1}] -> 2, Band[{1, 2}] -> -1, Band[{2, 1}] -> -1}, {n, n}]
L["Density"]
Normal[L[[1 ;; 4, 1 ;; 4]]]
Normal[(L + 1)[[1, 1 ;; 4]]]
{Total[L, 2], Min[L], Max[L]}
Transpose[L] === L
Take[L . Range[n], -2]
L2 = L . L
Normal[L2[[1 ;; 4, 1 ;; 4]]]
# The transparency gate: four ways a head can see a SparseArray
v = SparseArray[{2 -> 3}, 4];
Accumulate[v]
Sort[v]
Normal[Map[f, v]]
{MemberQ[L, 2], Count[L, -1]}
Normal[KroneckerProduct[IdentityMatrix[2, SparseArray], SparseArray[{{1, 2}, {0, 1}}]]]
IdentityMatrix[n, SparseArray]
