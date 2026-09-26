# 6.3 Slicing an association returns a sub-association
r = <|"a" -> 1, "b" -> 2, "c" -> 3|>
r[[2 ;; 3]]
r[[{"a", "c"}]]
r[[{1, 3}]]
