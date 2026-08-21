Module[{s = 0}, Do[s = s + i, {i, 1, 100}]; s]
Module[{x = 1}, For[i = 1, i <= 5, i++, x = x*i]; x]
Module[{n = 27, steps = 0}, While[n > 1, If[EvenQ[n], n = n/2, n = 3 n + 1]; steps++]; steps]
