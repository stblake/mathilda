# Section 8.1: building a compiled function.
Compile[{x}, x^2 + 1]
sq = Compile[{{x, _Real}}, x^2 + 1];
sq[3.5]
sq[3]
Compile[{{n, _Integer}}, n^2][7]
Compile[{{z, _Complex}}, z^2][1. + 2. I]
sq[y]
