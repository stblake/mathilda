# Section 8.6: a rank-1 argument compiles to one fused pass over the whole array.
sqp = Compile[{{u, _Real, 1}}, u^2 + 1.];
sqp[Range[1., 8.]]
big = sqp[Range[1., 1000000.]];
Length[big]
NDArrayQ[big]
