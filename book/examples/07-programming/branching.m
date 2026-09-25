# Ch.7 -- branching: If, Which, Switch
If[5 > 3, "yes", "no"]
classify[n_] := Which[n < 0, "negative", n == 0, "zero", True, "positive"]
classify /@ {-3, 0, 7}
sign[x_] := Which[x < 0, -1, x > 0, 1, True, 0]
sign /@ {-5, 0, 5}
kind[e_] := Switch[e, _Integer, "int", _Real, "real", _, "other"]
kind /@ {3, 2.5, x}
Switch[#, 1, "one", 2, "two", _, "many"] & /@ {1, 2, 3}
