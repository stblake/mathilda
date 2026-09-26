# 6.3 Operator forms, named slots, and value threading
rows = {<|"name" -> "fig", "qty" -> 1|>, <|"name" -> "date", "qty" -> 7|>, <|"name" -> "plum", "qty" -> 3|>};
Select[rows, #qty > 2 &]
SortBy[rows, #qty &]
inv = <|"apples" -> 3, "pears" -> 5|>;
Lookup["pears"][inv]
inv + 10
Total[{<|"a" -> 1, "b" -> 2|>, <|"a" -> 10, "b" -> 20|>}]
