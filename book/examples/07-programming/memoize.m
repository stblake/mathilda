# Ch.7 -- memoization: a delayed rule whose body installs a new immediate rule
fib[0] = 0;
fib[1] = 1;
fib[n_] := fib[n] = fib[n - 1] + fib[n - 2]
fib[100]
