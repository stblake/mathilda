# 4.8.6 Time series: moving statistics
# A short, noisy monthly series.
series = {6, 8, 5, 9, 11, 8, 13, 15, 12, 17, 16, 19};
MovingAverage[series, 3]
MovingAverage[series, {1, 2, 1}]
MovingMedian[series, 3]
ExponentialMovingAverage[series, 1/3]
N[ExponentialMovingAverage[series, 0.2]]
