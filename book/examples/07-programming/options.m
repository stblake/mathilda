# Ch.7 -- options: OptionsPattern, OptionValue, OptionQ, FilterRules
Options[draw] = {Color -> Red, Width -> 1};
draw[x_, OptionsPattern[]] := {x, OptionValue[Color], OptionValue[Width]}
draw[p]
draw[p, Width -> 3]
lead[OptionsPattern[], x_] := {x, OptionValue[draw, Width]}
lead[Width -> 2, q]
{OptionQ[Width -> 3], OptionQ[{Width -> 3, Color -> Blue}], OptionQ[3]}
FilterRules[{Color -> Blue, Width -> 2, Size -> 9}, Options[draw]]
FilterRules[{Color -> Blue, Width -> 2, Size -> 9}, Except[Options[draw]]]
