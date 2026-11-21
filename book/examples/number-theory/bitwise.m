# Bitwise operations on exact integers, two's complement at unlimited width.
BitAnd[12, 10]
BitOr[12, 10]
BitXor[12, 10]
BitNot[12]
BitAnd[-8, 255]
BitSet[0, 100] == 2^100
BitGet[2^100 + 1, {0, 1, 100}]
BitAnd[255, BitNot[2^0 + 2^7]]
BitShiftLeft[1, 70]
BitShiftRight[-5, 1]
