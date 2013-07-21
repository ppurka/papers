#------------------------------------------------------------------------#
#                                                                        #
#      Copyright (C) 2013. P. Purkayastha and H. M. Kiah                 #
#      License: GPL-3 or later.                                          #
#                                                                        #
#------------------------------------------------------------------------#
from sage.all import *
import config
#***********************************************************************#
#                                                                       #
# Some functions needed to work on the powerline channel.               #
#                                                                       #
#***********************************************************************#

cdef list add_ftmatrix(list f, list g):
    """
    add_ftmatrix(f, g)
    Add two ftmatrices modulo 2. The resulting matrix may have more than
    one 1 in a column or even an all zero column.

    INPUT:
        - `f`, `g` -- the ftmatrices

    EXAMPLES::

        sage: f = [[1,1]]*2; g = [[1,0]]*2; add_ftmatrix(f, g)
        [[0, 1], [0, 1]]
    """
    cdef int a, b
    cdef list f_row, g_row
    return [[(a+b)%2 for a,b in zip(f_row, g_row)]
            for f_row,g_row in zip(f, g)]

cdef list binary_symmetric_channel(list word, prob=0):
    """
    binary_symmetric_channel(word, prob=0)
    Flip the entries of word with probability prob. Hence it acts like
    a binary symmetric channel.

    INPUT:
        - ``word`` -- anything which is either one or two dimensional. It
          can be a list or tuple of list or tuple, or list, or tuple. The
          entries of it must be integers and either 0 or 1.

        - ``prob`` -- the probability with which to flip (default 0). You
          can also provide an instance of GeneralDiscreteDistribution on
          two discrete points.

    OUTPUT:
        - the resulting object as either a list or a list of list after its
          entries have been flipped.

    EXAMPLES::

        sage: binary_symmetric_channel([0]*4, prob=0.5) # random
        [1, 0, 0, 0]
        sage: D = GeneralDiscreteDistribution([0.5, 0.5])
        sage: binary_symmetric_channel([[0,0,0,0]]*2, prob=D) # random
        [[0, 1, 0, 0], [1, 1, 0, 0]]
    """
    cdef int _w, a, b
    cdef list w, errors
    if isinstance(prob, GeneralDiscreteDistribution):
        D = prob
    elif prob == 0:
        return word
    else:
        D = GeneralDiscreteDistribution([1-prob, prob])

    if isinstance(word[0], (list, tuple)):
        errors = [[D.get_random_element() for _w in w]
                  for w in word]
        return add_ftmatrix(word, errors)
    else:
        errors = [D.get_random_element() for _w in word]
        return [(a+b)%2 for a,b in zip(word, errors)]


cpdef list channel(Codeword word, impulse_noise=0, nb_noise=None, prob=0,
                   fade=0):
    """
    channel(word, impulse_noise=0, nb_noise=None, prob=0, fade=0)

    Simulate a powerline communication channel using only hard
    probabilities.

    INPUT:
        - ``word`` -- a :class:`Codeword`. This is the transmitted word.

        - ``impulse_noise -- impulse noise. If the parameter is a positive
          integer then it acts as a worst-case channel. If the parameter is
          between 0 and 1 then it acts as a random channel. The prob may be
          provided with an instance of GeneralDiscreteDistribution. This
          ensures that the discrete distribution is not reinitialized at
          every run of this function.

        - ``nb_noise`` --  narrow band noise. The parameter should be an
          instance of :class:`GeneralDiscreteDistribution`.

        - ``prob`` --  probability with which the additive white noise will
          be generated at any given time instance. This must be strictly in
          the interval [0,1]. The additive noise is generated for each
          frequency and time and so it will insert a 1 in the
          frequency-time matrix according to the frequency and time it
          corresponds to.  The prob may be provided with an instance of
          GeneralDiscreteDistribution. This ensures that the discrete
          distribution is not reinitialized at every run of this function.

        - ``fade`` -- the channel fade probability. Surprisingly, the fade
          is defined to erase a particular symbol at all time instances,
          according to the book on Powerline communications by Ferreira, et al

    OUTPUT: a frequency-time matrix containing codeword + fade + random
    noise + impulse noise + narrow-band noise, the noises being introduced
    in that order.

    In the frequency-time matrix, the alphabet is {0,1}. With 0 being
    reserved for the case when no errors occur. The position of 1 gives the
    frequency along the row and the time along the column. For instance,
    the matrix [[0,0],[1,0],[0,1]] corresponds to the codeword (1,2).

    EXAMPLES::

        sage: word = Codeword(range(2), 2); p = 0.5; word.ftmatrix
        [[1, 0], [0, 1]]
        sage: channel(word, prob=p) # random
        [[1, 0], [1, 0]]
        sage: channel(word, impulse_noise=2) # random
        [[1, 1], [1, 1]]
        sage: channel(word, nb_noise=1) # random
        [[1, 1], [0, 1]]
        sage: D = GeneralDiscreteDistribution([0.5]*2)
        sage: channel(word, prob=D) # random
        [[1, 0], [0, 0]]

    TESTS::

        sage: n = 1000; tx = Codeword(range(n), n); p = 0.003
        sage: D = GeneralDiscreteDistribution([1-p, p]); D.set_seed(0)
        sage: rx = channel(tx, prob=D); txm=tx.ftmatrix
        sage: N(sum([1 if r!=t else 0 for ri,ti in zip(rx,txm) for r,t in zip(ri,ti)])/n**2)
        0.00301000000000000

    Testing narrowband noise::

        sage: w = Codeword([7]*7, 7)
        sage: D = GeneralDiscreteDistribution([0.5]*2)
        sage: D.set_seed(0)
        sage: o = channel(w, nb_noise=D); print matrix(o)
        [0 0 1 1 1 1 1]
        [1 1 1 1 1 1 1]
        [1 1 1 1 1 1 1]
        [1 1 1 1 1 1 1]
        [0 0 0 0 0 0 0]
        [1 1 1 1 1 1 1]
        [0 0 0 0 0 0 0]
        sage: print config.nb_noise_global
        [37, 42, 37, 8, 0, 23, 0]
        sage: o = channel(w, nb_noise=D); print matrix(o)
        [1 1 1 1 1 1 1]
        [1 1 1 1 1 1 1]
        [1 1 1 1 1 1 1]
        [1 0 0 0 0 0 0]
        [0 0 0 0 0 0 0]
        [1 1 1 1 1 1 1]
        [0 0 0 0 0 0 0]
        sage: print config.nb_noise_global
        [23, 28, 23, 0, 0, 9, 0]
    """
    cdef int a, i, l, len_nb, n, q, im, ff, next_nb_indx, current_count_nb
    cdef list all_ones, output, all_zeros, partial_nb

    n = len(word)
    output = word.ftmatrix[:]
    q = len(output)

    # output after fading is added
    if isinstance(fade, GeneralDiscreteDistribution) or (0 < fade and
                                                            fade < 1):
        all_zeros = [0]*n
        for i,ff in enumerate(random_noise(q, 2, fade)):
            if ff:
                output[i] = all_zeros[:]



    # This is the channel output after random noise is added. At this step
    # the output contains the codeword+noise in the ftmatrix form.
    output = binary_symmetric_channel(output, prob)

    # Channel output after adding impulse_noise
    if isinstance(impulse_noise, GeneralDiscreteDistribution) or (
            0 < impulse_noise and impulse_noise < 1):
        all_ones = [1]*q
        output = zip(*output)
        for i,ff in enumerate(random_noise(n, 2, impulse_noise)):
            if ff:
                output[i] = all_ones[:]
        output = zip(*output)

    # Channel output after adding narrowband_noise
    if nb_noise is not None:
        # config.nb_noise_global contains all the *remaining* lengths of
        # narrowband noise at the frequencies given by its indices
        all_ones = [1]*n
        current_count_nb = len(find(config.nb_noise_global, '>', 0))
        # The max_partial_nb variable stores the maximum time index till
        # which all the *partial* narrowband noises are present. This is to
        # ensure that when we introduce the next narrowband noise, it
        # starts beyond that index. Otherwise it might lead to situation
        # where the narrowband noise is like this:
        # [1,1,1,1,1,1,0,0
        #  1,1,0,0,0,0,0,0
        #  ...
        #  At frequency 0, we might not introduce a narrowband noise, but
        #  at frequency 1, we might introduce a narrowband noise starting
        #  from the fourth position (for example). Thus the first two rows
        #  of the new matrix will look like
        # [1,1,1,1,1,1,0,0
        #  1,1,0,1,1,1,1,1
        #  ...
        # And, now at time index 3 onwards, we have got 2 narrowband
        # noises, whereas we maybe wanted only a maximum of 1 narrowband
        # noise to be present at any time instance.
        #
        partial_nb = [l for l in config.nb_noise_global if l > 0 and l < n]
        max_partial_nb = max(partial_nb) if partial_nb != [] else 0
        for i in range(q):
            len_nb = config.nb_noise_global[i]
            if len_nb >= n:
                output[i] = all_ones[:]
                config.nb_noise_global[i] = len_nb - n
                next_nb_indx = -1 # no more noise to add here
            else:
                if len_nb > 0:
                    # now we can decrease current_count_nb since there are no
                    # more nb noise at this frequency.
                    current_count_nb -= 1
                output[i] = [1]*len_nb + list(output[i][len_nb:])
                config.nb_noise_global[i] = 0
                # With probability nb_noise we decide whether to put
                # a narrowband noise at symbol i; then with probability
                # nb_noise we will again decide whether to put it at
                # len_nb+1, or len_nb+2, ..., or n.
                if nb_noise.get_random_element() == 1:
                    # len_nb is now the index from where next_nb_indx starts
                    # if len_nb was 1, then we take a random integer
                    # between 1 and n-1, and that denotes the *index* from
                    # where the next noise starts
                    if config.max_nb == q:
                        next_nb_indx = randint(len_nb, n-1)
                    else:
                        # See the comments above, before we define
                        # max_partial_nb
                        next_nb_indx = randint(max_partial_nb, n-1)
                else:
                    next_nb_indx = -1

            # Now determine a random number l between 1 and 10 that will be
            # used to generate l*n length nb noise. This is for every nb
            # noise at each symbol.
            # Add the noise only if the number of nb noise is less than the
            # max number that can be introduced.
            if current_count_nb < config.max_nb and next_nb_indx > -1:
                output[i] = output[i][:next_nb_indx] + [1]*(n-next_nb_indx)
                l = randint(1, 10)
                config.nb_noise_global[i] = l*n - (n-next_nb_indx)
                current_count_nb += 1

    # Some of the entries of output might be tuples. So, convert them back
    # to lists
    return map(list, output)


cpdef list channel_mfsk(Codeword cw, mfskchannel, additive_noise):
    tfmatrix = zip(*cw.ftmatrix) # time-freq matrix
    rx = []
    for tx in tfmatrix:          # at each time, we send a length q list
        transmitted = mfskchannel.modulate(tx)
        received    = mfskchannel.channel(transmitted, additive_noise, 0)

        rx.append(mfskchannel.demodulate(received))

    return map(list, zip(*rx))


cdef int equitable_distance(list cw, list ftmatrix):
    # This implements the "distance" in the paper
    cdef int ci, d=0, i
    for i,ci in enumerate(cw):
        if ftmatrix[ci][i] != 1:
            d += 1
    return d

cpdef int equitable_decoder(list ftmatrix, PowerlineCode code, Codeword tx,
        bool detect_nb=True, bool randomize=True):
    cdef int d, min_d, tmpd, n
    cdef Codeword c
    cdef list all_zeros, decoded_word, cww

    d = code.d
    n = code.n
    threshold = floor((n + config.r)/2)
    all_zeros = [0]*n

    # first remove the narrowband noises
    if detect_nb:
        for i,f in enumerate(ftmatrix):
            if sum(f) > threshold:
                ftmatrix[i] = all_zeros[:]

    min_d = n+1 # +1 so that decoded_word gets a vector at first run
    decoded_word = []
    for c in code._list:
        cww = c._word
        tmpd = equitable_distance(cww, ftmatrix)
        if tmpd < min_d:
            min_d = tmpd
            decoded_word = [cww]
        elif randomize and tmpd == min_d:
            decoded_word.append(cww)

    return hamming_distance(decoded_word[randint(0,len(decoded_word)-1)],
                            tx._word)

cdef int decoder_cwcode(NonlinearCode cwcode, list f):
    cdef int c, ff, min_d = cwcode.n, tmpd, k, q = cwcode.M
    cdef list cw, cwstore

    for k,cw in enumerate(cwcode._list):
        # We compare only against the non-erased coordinates and take the
        # codeword closest to the received vector.
        #
        # There is a difference. For example the output when we don't do
        # this modified rule is (for RS[15,3,13] code):
        # [(3/10, 17/150),  (1/5, 17/300),  (1/10, 17/300)]
        # with max_nb = 3, A(9,4,4) inner code, and with 100 words
        # transmitted, and with no noise other than narrowband noise.
        #
        # When we do this modified rule, by comparing only against
        # non-erased coordinates the output is a bit better:
        # [(3/10, 13/1500),  (1/5, 0),  (1/10, 0)]
        #
        # At 3/10 it should have been 0. Why are we getting 13/1500?
        # I think the answer lies in the very few cases where the number of
        # time instances where the narrowband noises are present is not
        # more than the threshold, and then the decoding fails.
        #
        # When we don't detect narrowband noise, there is no change in the
        # error rates between the unmodified and modified decoding rules.
        #
        tmpd = 0
        for c,ff in zip(cw,f):
            if ff == 2:
                continue
            elif c != ff:
                tmpd += 1
        if tmpd < min_d:
            min_d = tmpd
            q = k # this value is what we are after - the q-ary symbol
            cwstore = cw[:]

    return q

cpdef int decoder_mtfsk(list ftmatrix, code, Codeword tx,
        NonlinearCode cwcode, bool detect_nb=True, bool bdd=True,
        bool randomize=True):
    cdef int _f, fsum, min_d, tmpd, n, k, got_it = 0, count_nb = 0, q
    cdef list decoded_word, all_twos, f, cw, txw = tx._word, cww
    cdef Codeword c

    n = code.n
    q = config.Fstarlen+1
    threshold = floor((n + config.r)/2)
    all_twos = [2]*n
    cw = [0]*code.n
    orig_ftmatrix = ftmatrix[:]

    # First remove the narrowband noises, by marking the entries of
    # ftmatrix not with 0, but with a special number, in this case 2.
    if detect_nb:
        for i,f in enumerate(ftmatrix):
            if sum(f) > threshold:
                ftmatrix[i] = all_twos[:]
                count_nb += 1

    # Next, decode each column using the constant weight code.
    ftmatrix = map(list, zip(*ftmatrix))
    for k,f in enumerate(ftmatrix):
        fsum = 0
        for _f in f:
            fsum += (0 if _f == 2 else _f)
        # We need to check for impulse noise. But we could have erased some
        # rows and thereby removed some of the impulse noise. So, we need
        # to compare against the number q, less the number of narrowband
        # noise we removed.
        # Other option is if we have a constant weight code of constant
        # weight 1. Then we declare erasure if there are more than two 1's
        # in every column.
        if (fsum == q - count_nb) or (fsum != 1 and cwcode.is_cw1 == 1):
            cw[k] = q # an erasure
        else:
            cw[k] = decoder_cwcode(cwcode, f)

    # Next, decode the outer code
    if bdd:
        return code.bounded_distance_decoder(cw, txw)
        #dnew= code.bounded_distance_decoder(cw, txw)
        #if dnew > 0:
        #    print matrix(orig_ftmatrix),"\n"
        #    print "cw  =",cw,"\n","txw =",txw,"\n\n"
        #return dnew

    else:
        min_d = n+1 # +1 so that decoded_word gets a vector at first run
        for c in code._list:
            cww = c._word
            tmpd = hamming_distance(cww, cw)
            if tmpd < min_d:
                min_d = tmpd
                decoded_word = [cww]
            elif randomize and tmpd == min_d:
                decoded_word.append(cww)

        return hamming_distance(decoded_word[randint(0,len(decoded_word)-1)],
                                txw)

cdef list dft(list vec):
    cdef int i, n = len(vec), Fstarlen = config.Fstarlen
    cdef list dftvec, Fstarsorted = config.Fstarsorted
    dftvec = [sum([vec[i]*Fstarsorted[(i*j)%Fstarlen] for i in range(n)])
              for j in range(n)]
    return dftvec


cdef list idft(list vec):
    cdef int n = len(vec), Fstarlen = config.Fstarlen
    cdef list idftvec, Fstarsorted = config.Fstarsorted
    # We don't need to divide by n here if we work only over binary
    idftvec = [sum([vec[i]*Fstarsorted[(-i*j)%Fstarlen] for i in range(n)])
               for j in range(n)]
    return idftvec


cdef list find(list vec, str comp, int value):
    """
    Find the indices of vec(tor) which comp(are) to value.

    INPUT:
        - ``vec`` -- a list
        - ``comp`` -- a comparison operator provided as a string
        - ``value`` -- any number

    EXAMPLES::
        sage: vec = [1, 0, 2, 2, 1]; idx = find(vec, '==', 1); print idx
        [0, 4]
        sage: print find(vec, '>=', 1)
        [0, 2, 3, 4]
    """
    cdef int i, w
    if comp == '==':
        return [i for i,w in enumerate(vec) if w == value]
    elif comp == '>=':
        return [i for i,w in enumerate(vec) if w >= value]
    elif comp == '>':
        return [i for i,w in enumerate(vec) if w >  value]
    elif comp == '<=':
        return [i for i,w in enumerate(vec) if w <= value]
    elif comp == '<':
        return [i for i,w in enumerate(vec) if w <  value]
    else:
        raise SyntaxError("comp must be one of the comparison operators")

cdef int hamming_distance(x, y, int prevdist = -1):
    """
    Calculate the Hamming distance between two vectors.
    Both vectors must be of the same length.

    Usage: hamming_distance(list|vector, list|vector [,[prevdist =] int])

    prevdist can be used to speed up computations when the real Hamming
    distance is not desired. If it is provided, then the computation is
    stopped as soon as the value prevdist is met.

    Output: an integer.
    """
    cdef int count = 0, i, n = len(x)
    if prevdist == -1:
        prevdist = n

    for i from 0 <= i < n:
        if x[i] != y[i]:
            count += 1
            # Stop iterating if previous min is exceeded
            if count >= prevdist:
                return count
        # End of for
    return count

cdef int mylog(x, base=None):
    """
    Usage: mylog(x[, base])
    This is a very specific log function that should be used only for log
    over finite fields.  It returns 0 if the element is zero(F), and it
    returns n+1 if the element is a^n, where a is the primitive element of
    the finite field F.

    Inputs:
    - ``x`` -- The element of the field whose log is desired
    - ``base`` -- The base to which the log is computed (default: the
      multiplicative generator of the finite field)

    Output:
    - ``n+1`` if `x` is `base**n`, and `0` if `x` is the zero of the field

    EXAMPLES::
        sage: F = GF(4,'a'); mylog(F(1))
        1
        sage: mylog(F(0))
        0
        sage: mylog(F.multiplicative_generator()**2)
        3

    Note: no error checking is performed to ensure that `x` belongs to
    a finite field.
    """
    if x.is_zero():
        return 0
    if base is None:
        base = x.parent().multiplicative_generator()
    return log(x, base)+1


cdef list random_noise(int n, int q, prob=0):
    """
    random_noise(n, q, prob=0)
    Returns a vector with noise pattern from a symbol from {0,1,...,q-1}. If
    the symbol is 0, then there is no noise. This has probability 1-prob of
    occurring. The other alphabets occur with a probability of prob/(q-1).

    INPUT:
        - ``n`` -- the length of the vector
        - ``q`` -- the size of the alphabet
        - ``prob`` -- the probability with which to flip (default 0). You
          can also provide an instance of GeneralDiscreteDistribution on
          q discrete points.

    OUTPUT: a randomly generated list from {0,...,q-1}^n

    EXAMPLES::

        sage: random_noise(10, 2, 0.4)  # random
        [1, 0, 1, 0, 1, 0, 0, 0, 1, 0]
        sage: random_noise(10, 5, 0.4)  # random
        [1, 0, 0, 4, 0, 4, 0, 4, 0, 2]
    """
    cdef int i
    if isinstance(prob, GeneralDiscreteDistribution):
        D = prob
    elif prob == 0:
        return [0]*n
    else:
        D = GeneralDiscreteDistribution([1-prob] + [prob/(q-1)]*(q-1))

    return [D.get_random_element() for i in range(n)]

cdef list symbol_weight_distribution_vector(vector, q=None):
    """
    Return the symbol weight distribution of a vector.

    Usage: symbol_weight_distribution_vector(list|vector, q=None)

    Output: a list where the i-th element corresponds to the i-th element
    of the finite field (i-th element as we loop over the finite field)

    If the vector or list is not over a finite field, then the ordering is
    from 0,...,q-1, where q is the alphabet size. q must be provided as an
    argument to the function in this case.
    """
    #cdef Field F
    cdef list rvec
    cdef dict d

    if hasattr(vector[0], 'parent'):
        F = vector[0].parent()
        if sage.rings.finite_rings.constructor.is_PrimeFiniteField(F):
            # Use an even faster algo without sorting
            rvec = [0]*len(F)
            for x in vector:
                rvec[x] += 1
            return rvec

        elif sage.rings.finite_rings.constructor.is_FiniteField(F):
            rvec = [0]*len(F)
            # Initialize a dictionary
            d = {}
            for a in F:
                d[a] = 0
            for x in vector:
                d[x] += 1
            return [d[_] for _ in sorted(d)]
    else:
        if q is None:
            raise ValueError("If the vector is not over a finite field "
            "then you need to provide the alphabet size q")
        rvec = [0]*q
        for x in vector:
            rvec[x] += 1
        return rvec


cpdef int symbol_weight_matrix(code):
    cdef int n = code.n, r = 0, wt
    cdef Codeword c

    for c in code:
        wt = max(map(sum, c.ftmatrix))
        if wt == n:
            return n
        elif wt > r:
            r = wt
    return r

cpdef int symbol_weight_vector(vector, q=None):
    """
    Return the symbol weight of a vector.

    Usage: symbol_weight_vector(list|vector, q=None)

    Output: an integer.
    """
    # This is super slow, sometimes more than 6 times slower than the rest
    # Keeping it here for future recall
    #count = [ vector.count(x) for x in vector[0].parent() ]
    #return max(count)

    return max(symbol_weight_distribution_vector(vector, q))

#***********************************************************************#
#                                                                       #
# This is the Codeword class. Every element of this class is simply a   #
# python list. But it also comes with the alphabet size and most        #
# importantly, the frequency-time matrix which is the alternative       #
# representation of it.                                                 #
#                                                                       #
#***********************************************************************#
cdef class Codeword:
    r"""
    Every element of the :class:`Codeword` class is simply a python list.
    However, each element also has the attributes q, orig_word and ftmatrix
    where q denotes the size of the alphabet orig_word is the original
    input to the :class:`Codeword` class and ftmatrix is the frequency-time
    matrix representation of the orig_word.

    EXAMPLES::

        sage: c = Codeword([1,2,3,0], 4) # c/w over alphabet of size 4
        sage: 1 in c
        True
        sage: F.<a> = GF(4, 'a')

    If the coordinates are from a finite field then one need not give the
    alphabet size. The input will be automatically changed to
    a representation that contains only the powers of the primitive
    element. The zero of the finite field remains zero while the other
    elements which are a^i become i+1, for every i. ::

        sage: c = Codeword([a, F(1), a^2, F(0)]); c
        [2, 1, 3, 0]
        sage: c.ftmatrix # The frequency time representation
        [[0, 0, 0, 1], [0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0]]
        sage: c = Codeword(c, 4); c # You can again put c into Codeword
        [2, 1, 3, 0]

    If there is any coordinate >= q then it is considered an erasure
    symbol. Thus, the ftmatrix will have all zeros in that time
    coordinate. ::

        sage: c = Codeword([0,1,2], 2); c.ftmatrix
        [[1, 0, 0], [0, 1, 0]]
        sage: c+c # This adds modulo q on all but the erasure coordinates
        [0, 0, 2]

    """
    def __init__(self, word, cwcode=None, q=0):
        if q > 0:
            self.q = q
        elif hasattr(word[0], 'parent') and word[0].parent().is_finite():
            self.q = word[0].parent().order()
        else:
            raise ValueError("You must provide the alphabet size q when"
                            " the input word is not over a finite field")

        self.orig_word  = word
        self._word = self._gen_word()
        if cwcode is None: # Map the alphabet to CWC(q, 2, 1)_2, i |-> e_i
            e = [[1]+[0]*(self.q-1)]
            for _ in xrange(q-1):
                eprev = e[-1]
                e.append([eprev[-1]] + eprev[:q-1])
            cwcode = NonlinearCode(e, q=2, d=2)
            cwcode.is_cw1 = 1
        self.ftmatrix = self._ftmatrix(cwcode)

    def __add__(self, Codeword ww):
        cdef int q, ci, wi
        cdef list c,w

        c = self._word
        w = ww._word
        q = self.q
        return Codeword([(ci + wi)%q if ci != q and wi != q else q
                        for ci,wi in zip(c, w)], q)

    def __contains__(self, int b):
        return (b in self._word)

    def __getitem__(self, int i):
        return self._word[i]

    def __iter__(self):
        cdef int a
        for a in self._word:
            yield a

    def __len__(self):
        return len(self.orig_word)

    def __repr__(self):
        return str(self._word)

    def __sub__(self, Codeword w):
        cdef int q, i
        cdef list c

        c = self._word
        w = w._word
        q = self.q
        return Codeword([(c[i]-w[i])%q if c[i] != q and w[i] != q else q
                        for i in range(len(c))], q)

    cdef list _ftmatrix(self, NonlinearCode cwcode):
        cdef int w
        cdef list ftmat = [], clist = cwcode.list()
        for w in self._word:
            ftmat.append(clist[w])

        return map(list, zip(*ftmat))

    def _gen_word(self):
        word = self.orig_word
        if hasattr(word[0], 'parent') and word[0].parent().is_finite():
            base = word[0].parent().multiplicative_generator()
            return map(lambda x: mylog(x, base), word)
        else:
            return list(word)

#***********************************************************************#
#                                                                       #
# This is the class for a general NonlinearCode.                        #
# The hope is to move the min distance decoder and bounded distance     #
# decoder to this at a later stage, so that later other 'codes' such as #
# permutation codes can be defined simply by inheriting this class.     #
#                                                                       #
#***********************************************************************#
cdef class NonlinearCode:
    r"""
    This is the class for a general nonlinear code. The inputs are the list
    of codewords and optionally the alphabet size q and the minimum
    distance d of the code.

    The list of the codewords and the minimum distance are cached and not
    computed every time one needs them.

    Attributes present are the alphabet size q, the list of codewords list,
    the minimum distance d,
    a bounded distance decoder bounded_distance_decoder,
    a minimum distance decoder minimum_distance_decoder.

    Examples::

        sage: C = NonlinearCode([[1,2,3], [2,1,3], [1,4,5]], 6); C.d
        2
        sage: C.minimum_distance_decoder([1,1,1])
        [1, 2, 3]

        sage: RS = ReedSolomonCode(3,2,GF(4,'a'))
        sage: C = NonlinearCode(RS) # alphabet size not required here
        sage: C.minimum_distance()
        2
        sage: RS[0] in C
        True

    """

    def __init__(self, code, q=0, d=0):
        if isinstance(code, LinearCode):
            c = code.random_element()
        else:
            c = code[0]

        if q > 0:
            self.q  = q
        elif hasattr(c[0], 'parent') and c[0].parent().is_finite():
            self.q  = c[0].parent().order()
        else:
            raise ValueError("You must provide the alphabet size q when"
                            " the input word is not over a finite field")

        self._code  = code
        self._list  = []
        self.n      = len(c);
        self.size   = len(code);
        self.M      = self.size
        self._d     = d
        self.d      = self._get_min_dist();
        self.is_cw1 = 0     # hack for constant weight code of weight 1
                            # This is used during decoding.
                            # Read as "is (constant weight code of) weight 1"

    def __contains__(self, v):
        return (v in self._code) or (v in self.list())

    def __getitem__(self, i):
        C = self.list()
        return C[i]

    def __iter__(self):
        if isinstance(self._code, LinearCode):
            cw = self._code.random_element()
        else:
            cw = self._code[0]

        if isinstance(self._code, LinearCode):
            base = self._code.base_ring().multiplicative_generator()
            for c in self._code:
                yield map(lambda x: mylog(x, base), c)
        elif hasattr(cw[0], 'parent') and cw[0].parent().is_finite():
            base = cw[0].parent().multiplicative_generator()
            for c in self._code:
                yield map(lambda x: mylog(x, base), c)
        else:
            for c in self._code:
                yield c

    def __len__(self):
        return self.M

    def __repr__(self):
        return ("Nonlinear code of length %d and size %d"%(self.n, self.M)
                + " over an alphabet of %d elements"%(self.q))

    cpdef int _get_min_dist(self):
        cdef int d, dtmp, i, j
        if self._d > 0:
            return self._d
        if isinstance(self._code, LinearCode):
            self._d = self._code.minimum_distance()
            return self._d

        if isinstance(self._code, LinearCode):
            self._d = self._code.minimum_distance()
            return self._d

        C = self.list()
        d = self.n
        for i in xrange(self.M):
            for j in xrange(i+1, self.M):
                dtmp = hamming_distance(C[i], C[j], d)
                if dtmp < d:
                    d = dtmp
        self._d = d
        return d

    cpdef int bounded_distance_decoder(self, list rx, list tx):
        """
        Decode the received word rx using a bounded distance decoder.
        If tx is given then it can compare the decoded word with tx.

        rx, tx must be list (not a Codeword). This helps with the
        hamming_distance and find functions.

        Output is the decoded codeword if it can be decoded. Otherwise it
        is the received word.

        Warning! Exponential complexity!
        """
        cdef int num_erasures, radius

        num_erasures = len(find(rx, '>=', self.q))
        radius = self.d - 1 - num_erasures

        # 2*num_errors + num_erasures <= d-1 implies correct decoding
        if tx is not None:
            if 2*(hamming_distance(rx, tx) - num_erasures) <= radius:
                return 0

        if isinstance(self.list()[0], Codeword): # hack to make it work
            for c in self.list():
                if 2*(hamming_distance(rx, c._word) - num_erasures) <= radius:
                    return hamming_distance(c._word, tx)
        else:
            for c in self.list():
                if 2*(hamming_distance(rx, c) - num_erasures) <= radius:
                    return hamming_distance(c, tx)
        # It didn't return from above, so no codeword was found.
        return hamming_distance(rx, tx)


    def list(self):
        if self._list == []:
            self._list = list(self.__iter__())
        return self._list

    cpdef minimum_distance(self):
        return self.d

    cpdef random_element(self):
        return self.list()[randint(0, self.M-1)]

#***********************************************************************#
#                                                                       #
# This is the class for a general PowerlineCode.                        #
# The reason why this is implemented is because we need to be able to   #
# accept both linear and nonlinear codes. Also, every codeword should   #
# be an instance of the Codeword class so that we have access to the    #
# frequency-time matrix of the codeword.                                #
#                                                                       #
#***********************************************************************#
cdef class PowerlineCode(NonlinearCode):
    r"""
    This is the class for a general powerline code. The input to this can
    be either a LinearCode instance of Sage or just a list of codewords.
    One can also optionally provide the alphabet size q and the minimum
    distance d.

    The codewords are automatically converted to :class:`Codeword`
    instances so that the frequency-time matrices of each word are
    available.

    It inherits the NonlinearCode class and hence the decoders from that
    class are also available to it.

    EXAMPLES::

        sage: C = PowerlineCode(RS); C
        Powerline code of length 3 and size 16 over an alphabet of 4 elements
        sage: print C[2], RS[2]
        [3, 3, 3] (a + 1, a + 1, a + 1)
        sage: C.d # This minimum distance computation is cached
        2
        sage: C = PowerlineCode(RS, d=2) # To avoid computation, provide d
        sage: RS[2] in C # RS[2] is in C
        True

    """
    def __init__(self, code, cwcode=None, q=0, d=0):
        if cwcode is None: # Map the alphabet to CWC(q, 2, 1)_2, i |-> e_i
            if not q and hasattr(c[0], 'parent') and c[0].parent().is_finite():
                q  = c[0].parent().order()
            if not q:
                e = [[1]+[0]*(self.q-1)]
                for _ in xrange(q-1):
                    eprev = e[-1]
                    e.append([eprev[-1]] + eprev[:q-1])
                cwcode = NonlinearCode(e, q=2, d=2)
                cwcode.is_cw1 = 1
        self.cwcode = cwcode
        NonlinearCode.__init__(self, code, q, d)

    def __iter__(self):
        cw = self._code[0]
        if isinstance(self._code, LinearCode):
            for c in self._code:
                yield Codeword(c, self.cwcode, self.q)
        elif hasattr(cw[0], 'parent') and cw[0].parent().is_finite():
            for c in self._code:
                yield Codeword(c, self.cwcode, self.q)
        else:
            for c in self._code:
                yield Codeword(c, self.cwcode, self.q)

    def __repr__(self):
        return ("Powerline code of length %d and size %d"%(self.n, self.M)
                + " over an alphabet of %d elements"%(self.q)
                + " with inner code '%s'"%(self.cwcode))

    cpdef int _get_min_dist(self):
        cdef int d, dtmp, i, j
        if self._d > 0:
            return self._d
        if isinstance(self._code, LinearCode):
            self._d = self._code.minimum_distance()
            return self._d

        C = self.list()
        d = self.n
        for i in xrange(self.M):
            for j in xrange(i+1, self.M):
                dtmp = hamming_distance(C[i]._word, C[j]._word, d)
                if dtmp < d:
                    d = dtmp
        self._d = d
        return d

cdef class PowerlineLinearCode(PowerlineCode):
    """
    Special class for large linear codes.

    If you want a list of all codewords in the code C then call the
    function list(C.__iter__()).
    """
    def __init__(self, code, cwcode, q, d):
        # The distance is not an optional parameter now.
        PowerlineCode.__init__(self, code, cwcode, q, d)

    def __getitem__(self, i):
        raise NotImplementedError

    def list(self):
        # We don't want to list all the codewords.
        # So, we will keep this empty
        # self._list is set to this.
        return []

    cpdef random_element(self):
        """
        Returns a random element.
        """
        return Codeword(self._code.random_element(), self.cwcode, self.q)

cdef class PowerlineReedSolomonCode(PowerlineLinearCode):
    def __init__(self, cwcode, q, k, n=0):
        if n == 0:
            n = q-1
        d = n-k+1
        config.F = GF(q, 'a')
        config.one = config.F(1)
        config.zero = config.F(0)
        config.omega = (config.F).multiplicative_generator()
        config.Fstarlen = (config.F).order() - 1
        config.Fstarsorted = [config.omega**i for i in
                                                range(config.Fstarlen)]
        config.Rx = PolynomialRing(config.F, 'x')
        config.x = (config.Rx).gen()

        if n == q-1:
            pts = config.Fstarsorted[1:] + config.Fstarsorted[:1]
        else:
            pts = config.Fstarsorted[1:n+1]
        RScode = ReedSolomonCode(n, k, config.F, pts=pts)
        PowerlineLinearCode.__init__(self, RScode, cwcode, q, d)
        self.size = self.M = q**k

    cpdef int bounded_distance_decoder(self, list rx, list tx):
        # Berlekamp-Massey algorithm
        # Provided by Han Mao
        # Use this only with custom RS code, with j0 = 1
        # The default RS code in Sage does not have j0 = 1
        cdef int i, j, n = self.n, rho2 = 0, L, r, q = self.q
        #cdef int i, j, n = code.n, rho2 = 0, L, r, q = code.q
        cdef int num_erasures, radius
        cdef list Fstarsorted = config.Fstarsorted, LL, c
        cdef list receive2 = [config.zero]*n
        cdef object Lambx, Bx, x = config.x, one = config.one, Delta
        cdef object zero = config.zero

        # first check if it is already within the bounded distance radius
        # then we can skip the decoding
        num_erasures = len(find(rx, '>=', q))
        radius = self.d - 1 - num_erasures
        if 2*(hamming_distance(rx, tx) - num_erasures) <= radius:
            return 0

        Lambx = config.Rx(1)
        Bx = config.Rx(1)
        L = 0

        for i in range(n):
            if rx[i] >= q:
                Lambx *= (one - x*Fstarsorted[i])
                Bx = Lambx
                L += 1
                rho2 += 1
            else:
                receive2[i] = (zero if   rx[i] == 0
                                    else Fstarsorted[rx[i]-1])

        V = dft(receive2)
        S = V[:] #S = [V[(j+j0-1)%n] for j in range(n)]
        for r in range(rho2+1, n+1):
            LL = Lambx.coeffs()
            Delta = sum([LL[j] * S[(r-j)%n] for j in range(len(LL))])
            if r <= self.d-1:
                if Delta == zero:
                    Bx = x*Bx
                elif 2*L > r+rho2-1:
                    Lambx -= Delta*x*Bx
                    Bx = x*Bx
                else:
                    L = r - L + rho2
                    Lambx, Bx = (Lambx - Delta*x*Bx, (Delta**(-1))*Lambx)
            else:
                S[r%n] -= Delta

        c = idft([V[j] - S[j] for j in range(n)]) # j0 = 1
        tx = [zero if i == 0 else Fstarsorted[i-1] for i in tx]
        return hamming_distance(c, tx) #j =  hamming_distance(c, tx)
        #if j <= ceil(code.d/2): print c,"\n",tx
        #return j

cdef class PowerlineReedSolomonSubcode(PowerlineReedSolomonCode):
    """
    We remove anything that is of symbol weight n.
    """
    def __init__(self, cwcode, q, k, n=0):
        PowerlineReedSolomonCode.__init__(self, cwcode, q, k, n=n)
        v = vector(config.F, [1]*self.n)
        self.constant_words = [a*v for a in config.F]
        self.M -= q
        self.size = self.M

    def __iter__(self):
        for c in self._code:
            if c not in self.constant_words:
                yield Codeword(c, self.cwcode, self.q)

    cpdef random_element(self):
        c = self._code.random_element()
        while c in self.constant_words:
            c = self._code.random_element()
        return Codeword(c, self.cwcode, self.q)

cdef class PowerlineReedSolomonCosetCode(PowerlineReedSolomonCode):
    """
    We assume that the ReedSolomonCode got evaluated at the points
    `[a, a^2, ..., a^{(q-1)}]` where `F_q` is the finite field over
    which the ReedSolomonCode is defined. This implies that the
    ReedSolomonCode must be generated as, for example,::

        sage: F = GF(8, 'a'); a = F.multiplicative_generator()
        sage: C = ReedSolomonCode(7, 2, F, pts=[a**i for i in range(1,8)])

    """
    def __init__(self, cwcode, q, k, n=0):
        PowerlineReedSolomonCode.__init__(self, cwcode, q, k, n=n)
        # We want the codes to be code \subset outer_code, where outer_code
        # has distance one more than the distance of code. Since these are
        # Reed Solomon codes, the outer code's generator matrix will have
        # one more row: [(a^1)^k (a^2)^k (a^3)^k ... (a^(n))^k]. This row
        # will be our coset leader.
        if self.n == q-1:
            pts = config.Fstarsorted[1:] + config.Fstarsorted[:1]
        else:
            pts = config.Fstarsorted[1:self.n+1]

        # We need to generate the generator polynomial of the code of
        # distance one less. The evaluations of this polynomial will give
        # us the coset leader.
        x = config.x
        a = config.omega
        g2x = prod([x - a**i for i in range(1, self.d-1)])
        self.coset_leader = vector(config.F, [g2x.subs(x=p)**k for p in pts])

    def __iter__(self):
        for c in self._code:
            c = c + self.coset_leader
            yield Codeword(c, self.cwcode, self.q)

    cpdef int bounded_distance_decoder(self, list rx, list tx):
        # Berlekamp-Massey algorithm
        # Provided by Han Mao
        # Modified to work with coset codes.
        # Use this only with custom RS code, with j0 = 1
        # The default RS code in Sage does not have j0 = 1
        cdef int i, j, n = self.n, rho2 = 0, L, r, q = self.q
        #cdef int i, j, n = code.n, rho2 = 0, L, r, q = code.q
        cdef int num_erasures, radius
        cdef list Fstarsorted = config.Fstarsorted, LL, c
        cdef list receive2 = [config.zero]*n
        cdef object Lambx, Bx, x = config.x, one = config.one, Delta
        cdef object zero = config.zero, coset = self.coset_leader, txi, cxi

        Lambx = config.Rx(1)
        Bx = config.Rx(1)
        L = 0

        for i in range(n):
            if rx[i] >= q:
                Lambx *= (one - x*Fstarsorted[i])
                Bx = Lambx
                L += 1
                rho2 += 1
            else:
                # We need to subtract off the coset, but only from the
                # non-erased coordinates.
                receive2[i] = (zero if   rx[i] == 0
                                    else Fstarsorted[rx[i]-1]) - coset[i]

        V = dft(receive2)
        S = V[:] #S = [V[(j+j0-1)%n] for j in range(n)]
        for r in range(rho2+1, n+1):
            LL = Lambx.coeffs()
            Delta = sum([LL[j] * S[(r-j)%n] for j in range(len(LL))])
            if r <= self.d-1:
                if Delta == zero:
                    Bx = x*Bx
                elif 2*L > r+rho2-1:
                    Lambx -= Delta*x*Bx
                    Bx = x*Bx
                else:
                    L = r - L + rho2
                    Lambx, Bx = (Lambx - Delta*x*Bx, (Delta**(-1))*Lambx)
            else:
                S[r%n] -= Delta

        c = idft([V[j] - S[j] for j in range(n)]) # j0 = 1
        tx = [zero if i == 0 else Fstarsorted[i-1] for i in tx]
        tx = [txi - cxi for txi,cxi in zip(tx, coset)] # remove coset
        return hamming_distance(c, tx) #j =  hamming_distance(c, tx)
        #if j <= ceil(code.d/2): print c,"\n",tx
        #return j


    cpdef random_element(self):
        c = self._code.random_element() + self.coset_leader
        return Codeword(c, self.cwcode, self.q)

