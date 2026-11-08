# Holdout B, family D10: x_0 = x_1 = 1, x_{n+1} = a x_n + b x_{n-1} with
# (a, b) drawn independently each step from a finite list with probabilities
# p.  The ratio r_n = x_{n+1}/x_n obeys r -> a + b/r; |x_n|^(1/n) -> exp(lambda),
# lambda = Int log r dnu(r), nu the stationary law of that chain (Furstenberg).
from fractions import Fraction as F
CASES = {
 "bd10.1": [(F(1, 4), 1, 0), (F(1, 4), 2, 0), (F(1, 4), 3, 0), (F(1, 4), 4, 0)],   # x_{n+1} = c x_n: 24^(1/4)
 "bd10.2": [(F(1, 2), 3, 1), (F(1, 2), 3, -1)],
 "bd10.3": [(F(1, 3), 3, 1), (F(2, 3), 3, -1)],
 "bd10.4": [(F(1, 2), 1, 1), (F(1, 2), 1, 2)],
 "bd10.5": [(F(1, 3), 3, -1), (F(1, 3), 3, 0), (F(1, 3), 3, 1)],
 "bd10.6": [(F(1, 2), 4, 2), (F(1, 2), 4, -2)],
 "bd10.7": [(F(1, 2), 1, 1), (F(1, 2), 2, 1)],
 "bd10.8": [(F(1, 2), 5, 3), (F(1, 2), 5, -3)],
 "bd10.9": [(F(1, 2), 3, 2), (F(1, 2), 1, 2)],
 "bd10.10": [(F(1, 3), 1, 1), (F(1, 3), 2, 1), (F(1, 3), 3, 1)],
}
