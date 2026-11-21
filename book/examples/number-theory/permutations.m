# Permutations, cycles and permutation groups.
Cycles[{{3, 1, 2}, {5, 4}}]
PermutationList[Cycles[{{1, 3}}], 4]
PermutationCycles[{2, 3, 1, 4}]
Permute[{a, b, c, d}, Cycles[{{1, 2, 3}}]]
Signature[{3, 1, 2}]
Signature[{2, 1, 3}]
s3 = PermutationGroup[{Cycles[{{1, 2}}], Cycles[{{1, 2, 3}}]}];
GroupOrder[s3]
GroupElements[s3]
Permute[{a, b, c}, s3]
GroupOrder[PermutationGroup[{Cycles[{{1, 2, 3, 4, 5, 6, 7, 8}}], Cycles[{{1, 2}}]}]]
