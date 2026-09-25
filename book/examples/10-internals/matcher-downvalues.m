# DownValues are ordered by specificity, not by when they were entered.
fib[n_] := fib[n - 1] + fib[n - 2]
fib[0] = 0;
fib[1] = 1;
fib[10]
# The matcher tries the specific rules first; DownValues lists them that way:
DownValues[fib]
# The decisive test: a GENERAL rule entered AFTER a specific one still yields.
r[0] = special;
r[x_] := general[x]
DownValues[r]
r[0]
r[5]
