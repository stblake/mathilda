# Ch.7 -- constraining a match: PatternTest (?), Condition (/;), Alternatives (|), Except
Cases[Range[12], _?EvenQ]
Cases[Range[12], n_ /; n > 7]
Cases[{1, 2, 3, 4, 5, 6}, 2 | 4 | 6]
{-3, -1, 0, 2, 5} /. n_ /; n < 0 :> -n
Cases[{{1, 2}, {3, 1}, {4, 4}}, {a_, b_} /; a < b]
Cases[{1, 2, 3, "x", 4}, Except[_String]]
