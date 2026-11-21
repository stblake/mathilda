# SparseArray: the constructor, the atom, Normal and ArrayRules
s = SparseArray[{{1, 1} -> 1, {2, 3} -> x, {3, 2} -> -1}, {3, 3}]
{AtomQ[s], Head[s], Length[s], Dimensions[s]}
Normal[s]
ArrayRules[s]
InputForm[s]
s[[2, 3]]
s["NonzeroPositions"]
# Other constructor forms: a background value, pattern rules, bands, a dense array
Normal[SparseArray[{3 -> 7}, 5, 1]]
Normal[SparseArray[{{i_, i_} -> i}, {3, 3}]]
Normal[SparseArray[{Band[{1, 1}] -> 2, Band[{1, 2}] -> -1, Band[{2, 1}] -> -1}, {4, 4}]]
SparseArray[{{0, 1, 0}, {0, 0, 2}}]
