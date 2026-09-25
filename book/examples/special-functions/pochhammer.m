# 4.7.1 The Pochhammer symbol (rising factorial) (a)_n = a(a+1)...(a+n-1)
# A symbolic base gives a polynomial product
Pochhammer[x, 4]
# A numeric base gives an exact value; (1)_25 is 25!
Pochhammer[10, 6]
Pochhammer[1, 25]
# Half-integer arguments reduce through the gamma ratio
Pochhammer[1/2, 3]
Pochhammer[3/2, 1/2]
