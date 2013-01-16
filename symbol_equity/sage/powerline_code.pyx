#------------------------------------------------------------------------#
#                                                                        #
#      Copyright (C) 2013. P. Purkayastha and H. M. Kiah                 #
#      License: GPL-3 or later.                                          #
#                                                                        #
#------------------------------------------------------------------------#
from sage.all import *
import config # This file holds some global variables.

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
        - ``word`` -- anything which is either one or two dimensional. The
          entries of it must be integers 0 or 1.

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

        - ``impulse_noise -- If the parameter is between 0 and 1 then it
          acts as a random channel. The probability may also be provided
          with an instance of :class:`GeneralDiscreteDistribution`. This
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

        - ``fade`` -- the channel fade probability. This can be a real
          number in [0, 1] or a :class:`GeneralDiscreteDistribution'.
          Surprisingly, the fade is defined to erase a particular symbol at
          all time instances, according to the book on Powerline
          communications by Ferreira, et al

    OUTPUT:

        - a frequency-time matrix containing codeword + fade + random noise
          + impulse noise + narrow-band noise, the noises being introduced
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

        sage: import config
        sage: config.max_nb = 7; config.nb_noise_global = [0]*7
        sage: w = Codeword([7]*7, 7)
        sage: D = GeneralDiscreteDistribution([0.5]*2)
        sage: D.set_seed(0)
        sage: o = channel(w, nb_noise=D); print matrix(o)
        [0 0 1 1 1 1 1]
        [0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0]
        [0 0 0 0 0 1 1]
        [0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0]
        [0 0 0 0 0 1 1]
        sage: print config.nb_noise_global
        [2, 0, 0, 33, 0, 0, 12]
        sage: o = channel(w, nb_noise=D); print matrix(o)
        [1 1 0 0 1 1 1]
        [0 0 0 0 0 1 1]
        [0 0 1 1 1 1 1]
        [1 1 1 1 1 1 1]
        [0 0 0 0 1 1 1]
        [1 1 1 1 1 1 1]
        [1 1 1 1 1 1 1]
        sage: print config.nb_noise_global
        [53, 40, 16, 26, 11, 35, 5]
    """
    cdef int a, i, l, len_nb, n, q, im, ff, next_nb_indx, current_count_nb
    cdef list all_ones, output, all_zeros

    n = len(word)
    q = word.q
    output = word.ftmatrix[:]   # We need a copy of this
                                # For some reason .deepcopy() is not
                                # necessary, which is a relief since
                                # otherwise it is way slower.

    # output after fading is added
    if isinstance(fade, GeneralDiscreteDistribution) or (0 < fade and
                                                            fade < 1):
        fade = find(random_noise(q, 2, fade), '==', 1)
        if fade != []:
            all_zeros = [0]*n
            for ff in fade:
                output[ff] = all_zeros[:]


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
    # What is guaranteed is that at the end of this transmission the number
    # of narrowband noise will not exceed config.max_nb.
    if nb_noise is not None:
        # config.nb_noise_global contains all the *remaining* lengths of
        # narrowband noise at the frequencies given by its indices
        all_ones = [1]*n
        current_count_nb = len(find(config.nb_noise_global, '>', 0))
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
                    next_nb_indx = randint(len_nb, n-1)
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

cdef int equitable_distance(Codeword cw, list ftmatrix):
    """
    This implements the distance (the quasi-metric) introduced in the paper

    EXAMPLES:
    An example output of this function if it were callable from python.::

        sage: w = Codeword(range(2), 2); w.ftmatrix[0][1] = 1
        sage: print matrix(w.ftmatrix)
        [1 1]
        [0 1]
        sage: equitable_distance(Codeword(range(2),2), w.ftmatrix) #example
        0

    """
    cdef int ci, d=0, i
    for i,ci in enumerate(cw):
        if ftmatrix[ci][i] != 1:
            d += 1
    return d

cpdef int equitable_decoder(list ftmatrix, PowerlineCode code, Codeword tx,
        detect_nb=True, randomize=True):
    """
    Implements the decoder that uses the specific quasi-metric introduced
    in the paper and does a minimum distance decoding.

    INPUTS:

        - ``ftmatrix`` -- the output of the :meth:`channel`
        - ``code`` -- a :class:`PowerlineCode' instance
        - ``tx`` -- a :class:`Codeword` instance that was input to
          :meth:`channel`
        - ``detect_nb`` -- (default: True) boolean that indicates whether
          we should detect the presence of narrowband noise.
        - ``randomize`` -- (default: True) boolean that changes the
          behavior of the minimum distance decoder.
            - If True, then if we get more than one codeword at the same
              distance from the received word, then we randomly choose one
              of those codewords as the decoded word.
            - If False, then we select the first codeword out of the list.

    """
    cdef int d, min_d, tmpd, n
    cdef Codeword c
    cdef list all_zeros, decoded_word

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
        tmpd = equitable_distance(c, ftmatrix)
        if tmpd < min_d:
            min_d = tmpd
            decoded_word = [c]
        elif randomize and tmpd == min_d:
            decoded_word += [c]

    return hamming_distance(decoded_word[randint(0,len(decoded_word)-1)], tx)


cdef list find(vec, str comp, int value):
    """
    Find the indices of vec(tor) which comp(are) to value. The vec may be
    a vector or a list or a tuple.

    INPUTS:

        - ``vec`` -- a vector, tuple or list
        - ``comp`` -- a comparison operator provided as a string
        - ``value`` -- any number

    EXAMPLES:
    This is how the output would look if this were callable from python.::

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

# Usage: hamming_dist(list|vector, list|vector [,[prevdist =] int])
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

# Usage: list_codewords(linearcode, gen_matrix = False, return_list=True)
def list_codewords(linearcode, gen_matrix = False, return_list=True):
    """
    Give the list of codewords of a linear code. This function is often
    4 to 5 times faster than doing linear_code.list() in sage, and can
    sometimes be up to 100 times faster.

    Usage: list_codewords(code[,gen_matrix=True|False,return_list=True|False])

    code must be either:
        - a sage (linear) code
        - a matrix, which is a the generator matrix of the linear code.

    In the latter case one must provide the second argument to the
    function, and provide it as True. The matrix provided must be full
    ranked.

    If return_list is False, then the output of this function is
    a generator object

    Output: a list of vectors or a generator.

    See Ticket http://trac.sagemath.org/12014
    """
    codelist = []
    if gen_matrix:
        F = linearcode[0,0].parent()
        G = linearcode
        k = linearcode.nrows()
    else:
        F = linearcode.base_ring()
        G = linearcode.gen_mat()
        k = linearcode.dimension()

    # New code (after improvement by Rado). This is now almost 100 times
    # faster than Sage's implementation.
    def iterate(G, F):
        """the actual loop"""
        g = G[0]
        if G.nrows() == 1:
            for b in F:
                yield b*g
        else:
            for rest in iterate(G[1:], F):
                for b in F:
                    yield rest + b*g

    if return_list:
        return list(iterate(G, F))
    else:
        return iterate(G,F)

cpdef int mylog(x, base=None):
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

    EXAMPLES:
    If it were callable from python.::

        sage: random_noise(10, 2, 0.4)  # random
        [1, 0, 1, 0, 1, 0, 0, 0, 1, 0]
        sage: random_noise(10, 5, 0.4)  # random
        [1, 0, 0, 4, 0, 4, 0, 4, 0, 2]
    """
    if isinstance(prob, GeneralDiscreteDistribution):
        D = prob
    elif prob == 0:
        return [0]*n
    else:
        D = GeneralDiscreteDistribution([1-prob] + [prob/(q-1)]*(q-1))

    return [D.get_random_element() for _ in range(n)]


cdef list symbol_weight_distribution_vector(vector, int q=0):
    """
    Return the symbol weight distribution of the ``vector``.

    Output: a list where the i-th element corresponds to the i-th element
    of the finite field (i-th element as we loop over the finite field)

    If the vector or list is not over a finite field, then the ordering is
    from 0,...,q-1, where q is the alphabet size. q must be provided as an
    argument to the function in this case.
    """
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
        if q == 0:
            raise ValueError("If the vector is not over a finite field "
            "then you need to provide the alphabet size q")
        rvec = [0]*q
        for x in vector:
            rvec[x] += 1
        return rvec


cpdef int symbol_weight_vector(vector, int q=0):
    """
    Return the symbol weight of a vector.

    Usage: symbol_weight_vector(list|vector, q=None)

    Output: an integer.
    """
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
    cdef public int q
    cdef public list _word, ftmatrix
    cdef public object orig_word

    def __init__(self, word, q=0):
        if q > 0:
            self.q = q
        elif hasattr(word[0], 'parent') and word[0].parent().is_finite():
            self.q = word[0].parent().order()
        else:
            raise ValueError("You must provide the alphabet size q when"
                            " the input word is not over a finite field")

        self.orig_word  = word
        self._word = self._gen_word()
        self.ftmatrix = self._ftmatrix()

    def __add__(self, Codeword ww):
        cdef int q, i
        cdef list c,w

        c = self._word
        w = ww._word
        q = self.q
        # add elements mod q only if it is not q since q denotes an erasure
        return Codeword([(c[i]+w[i])%q if c[i] != q and w[i] != q else q
                        for i in range(len(c))], q)

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

    def __sub__(self, Codeword ww):
        cdef int q, i
        cdef list c, w

        c = self._word
        w = ww._word
        q = self.q
        return Codeword([(c[i]-w[i])%q if c[i] != q and w[i] != q else q
                        for i in range(len(c))], q)

    cdef list _ftmatrix(self):
        cdef int w, q
        return [[1 if w == q else 0 for w in self._word]
                for q in xrange(self.q)]

    cdef list _gen_word(self):
        word = self.orig_word
        if hasattr(word[0], 'parent') and word[0].parent().is_finite():
            return map(mylog, word)
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

    Attributes present are the alphabet size q, the list of codewords list,
    the minimum distance d,

    EXAMPLES::

        sage: C = NonlinearCode([[1,2,3], [2,1,3], [1,4,5]], 6); C.d
        2
        sage: C.minimum_distance_decoder([1,1,1])
        [1, 2, 3]

    """

    cdef public int _d, d, M, n, q, size
    cdef public list _list
    cdef public object _code

    def __init__(self, code, q=0, d=0):
        c = code[0]     # c is the first codeword
        if q > 0:
            self.q  = q
        elif hasattr(c[0], 'parent') and c[0].parent().is_finite():
            self.q  = c[0].parent().order()
        else:
            raise ValueError("You must provide the alphabet size q when"
                            " the input word is not over a finite field")

        self._code  = code
        if isinstance(code, LinearCode):
            self._list = map(lambda x: map(mylog, x),
                             list_codewords(code.gen_mat(),
                                 gen_matrix=True, return_list=False))
        elif hasattr(c[0], 'parent') and c[0].parent().is_finite():
            self._list = map(lambda x: map(mylog, x), code)
        else:
            self._list = code

        self.n      = len(c)
        self.size   = len(code)
        self.M      = self.size
        self._d     = d
        self.d      = self._get_min_dist()

    def __contains__(self, v):
        return (v in self._code) or (v in self._list)

    def __getitem__(self, i):
        C = self.list()
        return C[i]

    def __iter__(self):
        for c in self._code:
            yield c

    def __repr__(self):
        return ("Nonlinear code of length %d and size %d"%(self.n, self.M)
                + " over an alphabet of %d elements"%(self.q))

    cdef int _get_min_dist(self):
        cdef int d, dtmp, i, j
        if self._d > 0:
            return self._d
        C = self._list
        d = self.n
        for i in xrange(self.M):
            for j in xrange(i+1, self.M):
                dtmp = hamming_distance(C[i], C[j], d)
                if dtmp < d:
                    d = dtmp
        self._d = d
        return d


    cdef bounded_distance_decoder(self, rx, tx=None):
        # Removing this function since it's not used in this simulation.
        pass

    def list(self):
        return self._list

    cdef minimum_distance(self):
        return self.d

    cdef minimum_distance_decoder(self, rx, tx=None):
        # Removing this function since it's not used in this simulation.
        pass
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

    It inherits the NonlinearCode class and hence the methods from that
    class are also available to it.

    EXAMPLES::

        sage: C = PowerlineCode(RS); C
        Powerline code of length 3 and size 16 over an alphabet of 4 elements
        sage: print C[2], RS[2]
        [3, 3, 3] (a + 1, a + 1, a + 1)
        sage: C.d # This minimum distance computation is stored
        2
        sage: C = PowerlineCode(RS, d=2) # To avoid computation, provide d
        sage: RS[2] in C # RS[2] is in C
        True

    """
    def __init__(self, code, q=0, d=0):
        NonlinearCode.__init__(self, code, q, d)
        self._list = map(lambda x: Codeword(x, self.q), self._list)

    def __repr__(self):
        return ("Powerline code of length %d and size %d"%(self.n, self.M)
                + " over an alphabet of %d elements"%(self.q))
