#------------------------------------------------------------------------#
#                                                                        #
#      Copyright (C) 2013. P. Purkayastha and H. M. Kiah                 #
#      License: GPL-3 or later.                                          #
#                                                                        #
#------------------------------------------------------------------------#
from sage.all import *
from powerline_code import *
import config

#------------------------------------------------------------------------#
#                       Constant weight codes                            #
#------------------------------------------------------------------------#

#--------------------------- A(9, 4, 4) ---------------------------------#
# First - get the code corresponding to a Steiner system S(5,6,12), b=132
# This code is obtained by puncturing on one coordinate, to obtain a code
# of size 66. The codewords are the cyclic shifts of the following 6 words.
# See "Kramer Mesner -- Intersections among Steiner systems -- 1974.pdf"
a = [1,1,1,1,0,1,0,0,0,0,0]
b = [1,1,1,0,0,0,1,0,0,1,0]
c = [1,1,1,0,0,0,0,1,1,0,0]
d = [1,1,0,1,1,0,0,1,0,0,0]
e = [1,1,0,1,0,0,1,0,1,0,0]
f = [1,1,0,0,0,1,0,1,0,1,0]
C11_4_5 = []
for z in (a,b,c,d,e,f):
    for i in range(11):
        C11_4_5.append(z[:])
        z = z[1:] + [z[0]]
# Shorten this once on 0, to get to A(10,4,5) = 33.
C10_4_5 = [v[1:] for v in C11_4_5 if v[0] == 0]
# Shorten this once on 1, to get to A(9,4,4) = 18.
C9_4_4 = [v[1:] for v in C10_4_5 if v[0] == 1]
C9_4_4 = NonlinearCode(C9_4_4, q=2, d=4)


#--------------------------- A(13, 6, 5) --------------------------------#
# Getting constant weight code of size A(13, 6, 5) = 18 from
# Nordstorm-Robinson code.
# First we get the Golay code.
Gmat = ExtendedBinaryGolayCode().gen_mat()
# Make sure that the vector (1^8, 0^16) is in this code, by a permutation
# of the columns. This is obtained by a permutation of the first row.
g0 = Gmat[0]
g0indx = [i for i,x in enumerate(g0) if x == 1]
L = range(Gmat.ncols())
for i in g0indx:
    L.remove(i)
Gmat = Gmat[:, g0indx + L]
# Create the new linear code.
Golay = LinearCode(Gmat)
# Retain only the last 16 columns of those words which have either
#   1. all-zero in first eight columns, or
#   2. which have one 1 in first seven columns and 1 in eight column
# This is our nonlinear (16, 256, 8) Nordstorm-Robinson code.
NR = [v[8:] for v in Golay
            if ((v[:7]).hamming_weight() == 1 and v[7] == 1) or
               ((v[:8]).hamming_weight() == 0)
     ]
# Retain only those vectors which have weight 6.
# This is our constant weight code A(16,6,6) = 112.
C16_6_6 = [map(ZZ,v) for v in NR if v.hamming_weight() == 6]
# Shorten this once on 1, to get to A(15,6,5) = 42.
C15_6_5 = [v[1:] for v in C16_6_6 if v[0] == 1]
# Shorten this again on 0, to get to A(14,6,5) = 28.
C14_6_5 = [v[1:] for v in C15_6_5 if v[0] == 0]
# Shorten this once more on 0, to get to A(13,6,5) = 18.
C13_6_5 = [v[1:] for v in C14_6_5 if v[0] == 0]
C13_6_5 = NonlinearCode(C13_6_5, q=2, d=6)


#--------------------------- A(16, 2, 1) --------------------------------#
a = [1]+[0]*15
C16_2_1 = [a]
for i in range(15):
    a = a[1:] + a[:1]
    C16_2_1.append(a)
C16_2_1 = NonlinearCode(C16_2_1, q=2, d=2)
C16_2_1.is_cw1 = 1




#------------------------------------------------------------------------#
#                       Reed-Solomon codes                               #
#------------------------------------------------------------------------#
# Compare codes of rate 0.53
CodeRS_n15k14d2r15_q9    = PowerlineReedSolomonCode      (C9_4_4,  16, 14)
CodeRSC_n15k8d8r8_q16    = PowerlineReedSolomonCosetCode (C16_2_1, 16, 8)

# Compare codes of rate 0.3x
CodeRS_n15k14d2r15_q13   = PowerlineReedSolomonCode      (C13_6_5, 16, 14)
CodeRSC_n15k5d11r5_q16   = PowerlineReedSolomonCosetCode (C16_2_1, 16, 5)

