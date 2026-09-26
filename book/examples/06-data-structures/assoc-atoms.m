# 6.3 Associations are atoms: structural traversal sees the values
a = <|"x" -> 1, "y" -> 2, "z" -> 3|>
AtomQ[a]
Cases[a, _?OddQ]
Position[a, 2]
a /. ("y" -> "x")
