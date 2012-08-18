from powerline_code import *
#from code_functions import hamming_distance

#------------------------------------------------------------------------#
#                                                                        #
#             Code for DIM and binary 2^m by flips.                      #
#                                                                        #
#------------------------------------------------------------------------#
def binary_to_permutations_by_flip(cw, m):
    """
    Given the codeword cw in Z_{2^m}^n, return the permutation vector
    """
    sigma = range(1+m*len(cw))
    for i,c in enumerate(cw):
        for j,b in enumerate(get_binary(c, m)):
            if b == 1:
                indx = i*m+j
                tmp = sigma[indx+1]
                sigma[indx+1] = sigma[indx]
                sigma[indx]   = tmp
    return sigma

def bits_to_symbols(b, m):
    """
    Convert a binary list of length nm into a 2^m-ary list of length n by
    combining every m successive bits. The layout of the bits is MSB first.
    If there is an "erasure" (value 2) among the m consecutive bits, then
    simply make the symbol an erasure.
    """
    # Less computations on this one
    if m == 1: return b
    c = [0]*int(len(b)/m)
    i = 0
    j = 0
    q = 2**m # This denotes an erasure symbol
    s = 0
    twolist = [2**k for k in range(m-1,-1,-1)]
    for bj in b:
        if bj > 1:
            s = q
        else:
            s += twolist[j]*bj
        j += 1
        if j%m == 0:
            c[i] = s # If s is an erasure, this value is >= q
            i += 1
            j  = 0
            s  = 0
    return c

def decoder_for_binary_DIM(ftmatrix, tx, m, d, Clist):
    """
    Input:
        ftmatrix: the output of channel()
        tx: the q-ary Codeword()._word that is obtained from linear code
        m:  the power in 2^m
        d:  the minimum distance of the code
        Clist: the code as a list of lists in the Z_q representation
    Output:
        The number of symbols in error.

    If you want the decoded codeword then call code.<decoder function>
    directly.
    """
    # The following line could be replaced by a smarter decoder that keeps
    # all the positions which are in error (not erasure)
    y = ftmatrix_to_permutation(ftmatrix)
    b = estimate_bits(y)
    rx= bits_to_symbols(b, m)

    num_erasures = len(find(rx, '>=', 2**m))
    if 2*(hamming_distance(rx, tx) - num_erasures) <= d-1 - num_erasures:
        return 0

    # Check if it is closer to some other codeword than to tx
    for c in Clist:
        if 2*(hamming_distance(rx, c) - num_erasures) <= d-1-num_erasures:
            return hamming_distance(c, tx)

    # rx is outside any ball of radius d/2. Decoding failure.
    return hamming_distance(rx, tx)


def estimate_bits(y):
    """
    Given the received vector y, estimate the bits. This does not convert
    to the 2^m ary symbols.
    """
    N = len(y)
    b = [0]*N
    for i,yi in enumerate(y):
        if yi == i+1:
            b[i] = 1
        elif yi < i+1:
            b[i] = 0
        else:
            b[i] = 2
    return b[:N-1]

def ftmatrix_to_permutation(ftmatrix):
    """
    This function converts a frequency-time matrix to a Codeword. At time
    instance i, it looks only on the frequencies which are less than i+1.
    The code is identical to the other function ftmatrix_to_Codeword except
    for the step where we check for an erasure (the last for loop).
    The steps followed are:
    1. If there is a list of all 1's in a particular frequency, then it is
       set to an all 0 list
    2. If there is a list of all 1's in a particular time, then it is set
       to an all 0 list
    3. If there is all 0 time coordinate then it is set to q (an erasure
       symbol)
    4. If there are more than two 1's at time instance i then if one of the
       1's is at (i+1,i)-th frequency-time position, then the coordinate is
       assumed to be i+1. If there are more than two 1's in (j,i)-th
       frequency-time positions with j<i then the first frequency index of
       the 1 is assumed as the coordinate.
    5. Otherwise it returns the index of the 1 in a particular time.

    Output is a Codeword with symbols from {0,...,q} with q denoting an
    erasure symbol.
    """
    n = len(ftmatrix[0])
    q = len(ftmatrix)
    # Find all the symbols which have narrow band errors
    narrow_band_err = [e for e,l in enumerate(ftmatrix)
                        if len(find(l, '==', 1)) == n]

    # Find all the impulse noise errors
    impulse_err = [e for e,l in enumerate(zip(*ftmatrix))
                    if len(find(l, '==', 1)) == q]

    # Now to find the other errors. First delete the narrow and impulse err
    all_zeros = [0]*n
    for e in narrow_band_err:
        ftmatrix[e] = all_zeros[:]
    ftmatrix = zip(*ftmatrix) # "transpose" the matrix into a time-freq one
    all_zeros = [0]*q
    for e in impulse_err:
        ftmatrix[e] = all_zeros[:]
    # *******************************************************************#
    # Note that ftmatrix is now transposed. It is a time-freq matrix now #
    # *******************************************************************#

    # Convert ftmatrix into a Codeword (this is a time-freq matrix now)
    cw = []
    Q = [q]
    for i,f in enumerate(ftmatrix):
        E = [_ for _ in find(f, '==', 1) if _ <= i+1]
        # if the i+1 th element is in the list of received symbols then
        # assume that the i-th bit is 1
        # else if there is at least element in E (that will be at <= i)
        # assume that the i-th bit is 0
        # else E is empty, so it is an erasure.
        if i+1 in E:
            cw += [i+1]
        elif len(E) >= 1:
            cw += [E[0]]
        else:
            cw += Q
    return Codeword([ZZ(_) for _ in cw], q)

def get_binary(c, m):
    """
    get the binary representation of the symbol c as a list of length m.
    """
    b = [0]*m
    for i in range(m-1, -1, -1):
        b[i] = c%2
        c = (c - b[i])/2
    return b



#------------------------------------------------------------------------#
#                                                                        #
#             Code for DPM and q-ary obtained by additions.              #
#                                                                        #
#------------------------------------------------------------------------#

def argmax(L):
    """
    Return the index of the element which is the maximum. Assumes all
    non-negative elements.
    If there are more than one element which are equal, it will return the
    index of the first element which is the max.
    """
    return find(L, '==', max(L))[0]

def decoder_for_DPM(ftmatrix, tx, q, d, Clist):
    y = ftmatrix_to_Codeword(ftmatrix)
    if q == 2:
        s = estimate_binary_DPM(y)
    elif q == 3:
        s = estimate_ternary_DPM(y)
    else:
        raise NotImplementedError

    num_erasures = len(find(s, '>=', q))
    if 2*(hamming_distance(s, tx) - num_erasures) <= d-1 - num_erasures:
        return 0

    # Check if it is closer to some other codeword than to tx
    for c in Clist:
        if 2*(hamming_distance(s, c) - num_erasures) <= d-1-num_erasures:
            return hamming_distance(c, tx)

    # rx is outside any ball of radius d/2. Decoding failure.
    return hamming_distance(s, tx)


def estimate_binary_DPM(y):
    b = [0]*(len(y)-1)
    e = 2           # This is the erasure bit
    L = []          # Unlike the algorithm in paper, these are not the
                    # indices but the actual symbols.
    N = len(y)      # This is the symbol size of permutation space
    for i in range(len(y)-1):
        r = y[i+1]
        # Note that r points to y[i+1] and erasure => y[i] >= N
        if y[i] < N: # and len(L) < l0 % This can be used to restrict |L|
            L += [y[i]]
        if r >= N:
            b[i] = e
        else:
            t = 0
            for yl in L:
                if yl - r > 0:
                    t += 1
                else:
                    t -= 1
            if t > 0:
                b[i] = 1
            elif t < 0:
                b[i] = 0
            else:
                b[i] = e

    return b


def estimate_ternary_DPM(y):
    e = 3          # the erasure symbol for ternary
    L = []         # the list of non-erased permutation symbols
    N = len(y)     # the alphabet size of permutation space
    n = (N-1)/2    # the length of ternary vector
    s = [0]*n      # this will hold the estimated ternary symbols

    for j in range(1, n+1):
        J = 2*(j-1)
        if J-1 >= 0 and y[J-1] < N:
            L += [y[J-1]]
        if y[J] < N:
            L += [y[J]]

        J = 2*j
        if y[J] >= N or y[J-1] >= N:
            s[j-1] = e
        else:
            t = [0]*3
            for yl in L:
                p = [yl-y[J-1], yl-y[J]]
                if   p[0] < 0 and p[1] < 0: t[0] += 1
                elif p[0] > 0 and p[1] > 0: t[2] += 1
                elif p[0] < 0 and p[1] > 0: t[1] += 1

            if t[0] == 0 and t[1] == 0 and t[2] == 0:
                s[j-1] = e
            else:
                s[j-1] = argmax(t)

    return s

def qary_shift_left(L, b, m):
    """
    Shift the first m symbols cyclically by b, b in {0,...,q-1}
    """
    return [(_+ZZ(b))%m for _ in L[:m]] + L[m:]

def vector_to_permutation(v, q, shift=None):
    """
    map q-ary vector to permutation
    v = the q-ary vector
    q = q of q-ary
    shift = no. of extra coordinates that will be involved on next symbol
    """
    if shift is None:
        shift = q-1 # Thisis true for q=2,3
    n = len(v)

    def generate_permutation(L, u, m):
        if len(u) == 1:
            return qary_shift_left(L, u[0], m)
        return generate_permutation(qary_shift_left(L, u[0], m), u[1:],
                                    m+shift)

    return generate_permutation(range(q+shift*(n-1)), v, q)
