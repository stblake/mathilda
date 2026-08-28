# 4.8.4 Shape: moments, skewness, kurtosis
# Citation counts of ten papers -- one runaway hit gives a long right tail.
cites = {1, 2, 2, 3, 3, 3, 4, 4, 5, 20};
Moment[cites, 2]
CentralMoment[cites, 2]
CentralMoment[cites, 3]
Skewness[cites]
N[Skewness[cites]]
Kurtosis[cites]
N[Kurtosis[cites]]
