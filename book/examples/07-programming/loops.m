# Ch.7 -- loops: Do, While, For, and Break/Continue
s = 0; Do[s = s + i, {i, 1, 100}]; s
p = 1; Do[p = p i, {i, 1, 6}]; p
n = 1; While[n < 100, n = 2 n]; n
total = 0; For[k = 1, k <= 10, k++, total = total + k^2]; total
Do[If[i > 3, Break[]]; last = i, {i, 1, 10}]; last
r = 0; Do[If[EvenQ[i], Continue[]]; r += i, {i, 1, 10}]; r
