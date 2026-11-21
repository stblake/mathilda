# Wavelets: the Haar and Daubechies scaling and wavelet functions
WaveletPsi[HaarWavelet[], x]
h = WaveletFilterCoefficients[DaubechiesWavelet[2], "PrimalLowpass"]
Total[h[[All, 2]]]
# phi at the integers: (1 + Sqrt[3])/2 and (1 - Sqrt[3])/2
{WaveletPhi[DaubechiesWavelet[2], 1.], WaveletPhi[DaubechiesWavelet[2], 2.]}
# The refinement equation phi(x) = 2 Sum h_k phi(2x - k), on the dyadic grid
{WaveletPhi[DaubechiesWavelet[2], 1.25], 2 Sum[h[[k + 1, 2]] WaveletPhi[DaubechiesWavelet[2], 2.5 - k], {k, 0, 3}]}
WaveletPhi[DaubechiesWavelet[2], 1.3]
