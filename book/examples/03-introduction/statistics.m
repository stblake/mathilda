# Chapter 3 tour: Statistics
# Ten exam scores (out of 100) and the hours each student studied.
scores = {55, 62, 68, 71, 74, 78, 83, 88, 91, 95};
Mean[scores]
Median[scores]
StandardDeviation[scores] // N
Quartiles[scores]
hours = {2, 3, 3, 4, 5, 5, 6, 7, 8, 9};
Correlation[hours, scores] // N
MovingAverage[scores, 3]
