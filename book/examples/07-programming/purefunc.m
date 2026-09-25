# Ch.7 -- pure functions: Function, &, Slot, SlotSequence
(#^2 &)[5]
(#1 + #2 &)[3, 4]
Function[u, u^3][2]
Function[{x, y}, x^2 + y^2][3, 4]
Map[#^2 &, {1, 2, 3, 4}]
Select[Range[20], # > 15 &]
f[##, ##] &[a, b]
h[#2, #1] &[x, y]
