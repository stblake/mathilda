Interval[{-1, 1}] - Interval[{-1, 1}]
Nest[4 # (1 - #) &, Interval[0.5], 3]
