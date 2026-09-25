# Ch.7 -- Map, Apply, MapThread, Thread
Map[box, {a, b, c}]
box /@ {a, b, c}
Apply[Plus, {1, 2, 3, 4}]
Plus @@ {1, 2, 3, 4}
Apply[f, {a, b, c}]
f @@@ {{1, 2}, {3, 4}}
MapThread[f, {{a, b, c}, {x, y, z}}]
Thread[f[{a, b, c}, u]]
Thread[{a, b, c} == {1, 2, 3}]
