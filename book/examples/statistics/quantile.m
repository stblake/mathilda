# 4.8 Quantiles and robust spread
# Michelson's 1879 speed-of-light batch again (integers, so the estimates stay
# exact). n = 20, an even count -- which is exactly where the quantile
# conventions part company.
michelson = {850, 740, 900, 1070, 930, 850, 950, 980, 980, 880, 1000, 980, 930, 650, 760, 810, 1000, 1000, 960, 960};
Quantile[michelson, 1/2]
Median[michelson]
Quantile[michelson, {1/4, 1/2, 3/4}]
Quantile[michelson, 9/10]
InterquartileRange[michelson]
MeanDeviation[michelson]
MedianDeviation[michelson]
N[MeanDeviation[michelson]]
