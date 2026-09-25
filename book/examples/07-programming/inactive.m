# Ch.7 -- Inactive / Activate: hold a head inert, then release it
Inactive[Plus][2, 3]
Inactive[Plus][1 + 1, 2 + 2]
Activate[Inactive[Plus][2, 3]]
Inactive[Integrate][1/Sqrt[y Log[y] + 3], y]
D[Inactive[Integrate][f[u], u], u]
Activate[Inactive[Integrate][2 y, y]]
