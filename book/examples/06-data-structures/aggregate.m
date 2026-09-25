# 6.3 Aggregating with associations: Counts, GroupBy, Merge
Counts[{"red", "blue", "red", "green", "red", "blue"}]
GroupBy[Range[10], EvenQ]
GroupBy[Range[10], EvenQ, Total]
sales = {{"fruit", 3}, {"veg", 5}, {"fruit", 7}, {"veg", 2}, {"fruit", 1}};
GroupBy[sales, First -> Last]
GroupBy[sales, First -> Last, Total]
Merge[{<|"x" -> 1, "y" -> 2|>, <|"y" -> 20, "z" -> 30|>}, Total]
