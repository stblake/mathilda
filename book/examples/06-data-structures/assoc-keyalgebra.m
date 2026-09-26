# 6.3 Key-set algebra and relational joins over associations
KeyIntersection[{<|"a" -> 1, "b" -> 2|>, <|"b" -> 9, "c" -> 4|>}]
KeyComplement[{<|"a" -> 1, "b" -> 2|>, <|"b" -> 9|>}]
KeyUnion[{<|"a" -> 1|>, <|"b" -> 2|>}]
Discard[<|"a" -> 1, "b" -> 2, "c" -> 3|>, OddQ]
JoinAcross[{<|"id" -> 1, "x" -> 10|>}, {<|"id" -> 1, "y" -> 20|>}, "id"]
PositionIndex[{"x", "y", "x", "z", "x"}]
Merge[{"x" -> 1, "x" -> 2, "y" -> 3}, Total]
