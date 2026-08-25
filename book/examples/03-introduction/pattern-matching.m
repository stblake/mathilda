# Define a function only on positive integers, then try it on and off the pattern.
f[x_Integer?Positive] := x^2
f[7]
f[-3]
f[2.5]
f[3 + 2 I]
Cases[{1, 2, 3, 4, 5, 6}, _?EvenQ]
{a, b, c} /. b -> 99
area[{x_, y_}] := x y
area[{6, 7}]
