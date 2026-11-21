# Ch.7 -- messages: Check, Quiet, Message, Off/On, $MessageList
Check[1/0, "caught"]
Check[Sin[x], "caught"]
Check[{1/0, 0^0}, "caught", Power::indet]
Check[1/0, "caught", Power::indet]
Quiet[{1/0, 0^0}, {Power::infy}]
Check[Quiet[1/0], "caught"]
Quiet[Check[1/0, "caught"]]
Block[{$MessageList = {}}, 1/0; 0^0; $MessageList]
nonneg::neg = "Argument `1` should be non-negative.";
safeSqrt[x_] := (If[x < 0, Message[nonneg::neg, x]]; Sqrt[Abs[x]])
Check[safeSqrt[-4], $Failed]
Off[nonneg::neg]
Check[safeSqrt[-4], $Failed]
On[nonneg::neg]
Check[Delete[{1, 2, 3}, 5], "out of range"]
Check[Insert[{1, 2, 3}, x, 10], "out of range"]
