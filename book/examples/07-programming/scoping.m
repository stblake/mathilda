# Ch.7 -- scoping: Module (fresh locals), With (substitution), Block (dynamic)
x = 100;
Module[{x = 2}, x^2]
x
With[{x = 5}, x^3]
With[{a = 1 + 1}, Hold[a + z]]
Module[{c = 0}, Scan[(c = c + 1) &, Range[5]]; c]
n = 10;
g[] := n^2
Block[{n = 3}, g[]]
Module[{n = 3}, g[]]
g[]
