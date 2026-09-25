# 6.3 Associations: the key-value store
inv = <|"apples" -> 3, "pears" -> 5, "plums" -> 2|>
inv["pears"]
inv[["plums"]]
Keys[inv]
Values[inv]
Lookup[inv, "apples"]
Lookup[inv, "grapes", 0]
KeyExistsQ[inv, "plums"]
KeyExistsQ[inv, "grapes"]
inv["grapes"]
<|"a" -> 1, "b" -> 2, "a" -> 99|>
Map[#^2 &, <|"x" -> 3, "y" -> 4|>]
Total[inv]
AssociationThread[{"x", "y", "z"}, {10, 20, 30}]
