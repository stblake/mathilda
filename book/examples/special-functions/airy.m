# 4.7.4 The Airy functions
# The values at the origin close through the gamma function
AiryAi[0]
AiryBi[0]
AiryAiPrime[0]
# Numerically
N[AiryAi[0], 30]
N[AiryAi[1], 20]
# AiryAi differentiates to its own derivative head
D[AiryAi[x], x]
