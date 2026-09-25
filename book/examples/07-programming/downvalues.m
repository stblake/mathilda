# Ch.7 -- defining functions as DownValues (f[x_] := ...)
square[x_] := x^2
{square[7], square[a + b], square[Pi]}
area[r_] := Pi r^2
area[b_, h_] := b h
{area[3], area[3, 4]}
fib[0] = 0;
fib[1] = 1;
fib[n_] := fib[n - 1] + fib[n - 2]
fib[10]
collatz[m_?EvenQ] := m/2
collatz[m_?OddQ] := 3 m + 1
NestList[collatz, 27, 10]
