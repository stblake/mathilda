# 4.8.2 Spread: variance, standard deviation, RMS
# Michelson's 1879 speed-of-light measurements (first twenty), in km/s above
# 299000. Integer data, so the estimates come back exact.
michelson = {850, 740, 900, 1070, 930, 850, 950, 980, 980, 880, 1000, 980, 930, 650, 760, 810, 1000, 1000, 960, 960};
Mean[michelson]
Variance[michelson]
StandardDeviation[michelson]
N[StandardDeviation[michelson]]
RootMeanSquare[michelson]
