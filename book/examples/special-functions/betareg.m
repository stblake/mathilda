# The regularized incomplete beta function I_z(a, b)
BetaRegularized[1/2, 2, 3]
BetaRegularized[0.3, 2.5, 3.7]
N[BetaRegularized[3/10, 5/2, 37/10], 30]
N[Beta[3/10, 5/2, 37/10]/Beta[5/2, 37/10], 30]
# The reflection I_z(a, b) + I_(1-z)(b, a) = 1
BetaRegularized[0.3, 2.5, 3.7] + BetaRegularized[0.7, 3.7, 2.5]
D[BetaRegularized[z, a, b], z]
# A far tail: the binomial-tail probability of a large sample
BetaRegularized[0.3, 200., 300.]
