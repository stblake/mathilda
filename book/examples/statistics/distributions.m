# 4.8.8 From data to distributions
# Simulate a thousand IQ scores, N(100, 15), then measure the sample back.
SeedRandom[2024]
sample = RandomVariate[NormalDistribution[100, 15], 1000];
Mean[sample]
StandardDeviation[sample]
PDF[NormalDistribution[100, 15], {85., 100., 115.}]
Head[LearnDistribution[sample]]
