# Ch.7 -- filtering with patterns: Cases, Count, Position, Select, DeleteCases
data = {1, "two", 3.0, 4, x, 5, y^2};
Cases[data, _Integer]
Cases[data, _?NumberQ]
Count[data, _Integer]
Position[{a, b, c, b, a}, b]
Select[Range[20], PrimeQ]
DeleteCases[{1, 2, Missing[], 4, Missing[]}, _Missing]
Cases[{{1, 2}, {3, 4}, {5, 6}}, {a_, b_} :> a + b]
