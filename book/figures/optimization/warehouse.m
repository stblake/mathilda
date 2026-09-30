# Featured: place a warehouse to minimize the demand-weighted straight-line distance
# to five stores (the Fermat-Weber problem).
stores = {{0, 0}, {4, 1}, {5, 5}, {1, 6}, {-2, 3}}; demand = {3, 1, 2, 1, 1};
cost = Total[demand Map[Sqrt[(x - #[[1]])^2 + (y - #[[2]])^2] &, stores]];
sol = NMinimize[cost, {x, y}]
w = {x, y} /. Last[sol];
fig = Graphics[{Gray, Line[{w, #}] & /@ stores, Blue, Disk[#, 0.12 Sqrt[#2]] & @@@ Transpose[{stores, demand}], Red, Disk[w, 0.18]}, Frame -> True, AspectRatio -> Automatic]
