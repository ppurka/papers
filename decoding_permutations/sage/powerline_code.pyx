from sage.all import *
#from code_functions import hamming_distance, mylog

#***********************************************************************#
#                                                                       #
# functions from code_functions.pyx                                     #
#                                                                       #
#***********************************************************************#
# Usage: hamming_dist(list|vector, list|vector [,[prevdist =] int])
cpdef int hamming_distance(x, y, int prevdist = -1):
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




#***********************************************************************#
#                                                                       #
# Some functions needed to work on the powerline channel.               #
#                                                                       #
#***********************************************************************#

cdef add_ftmatrix(f, g):
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
    return [[(a+b)%2 for a,b in zip(f_row, g_row)]
            for f_row,g_row in zip(f, g)]

cdef binary_symmetric_channel(word, prob=0):
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
    if isinstance(prob, GeneralDiscreteDistribution):
        D = prob
    elif prob == 0:
        return word
    else:
        D = GeneralDiscreteDistribution([1-prob, prob])

    if isinstance(word[0], (list, tuple)):
        errors = [[D.get_random_element() for _ in w]
                  for w in word]
        return add_ftmatrix(word, errors)
    else:
        errors = [D.get_random_element() for _w in word]
        return [(a+b)%2 for a,b in zip(word, errors)]


cpdef list channel(word, impulse_noise=0, nb_noise=0, prob=0):
    """
    channel(word, impulse_noise=0, nb_noise=0, prob=0)

    Simulate a powerline communication channel using only hard
    probabilities.

    INPUT:
        - ``word`` -- a :class:`Codeword`. This is the transmitted word.

        - ``impulse_noise -- impulse noise. If the parameter is a positive
          integer then it acts as a worst-case channel. If the parameter is
          between 0 and 1 then it acts as a random channel.

        - ``nb_noise`` --  narrow band noise. If the parameter is
          a positive integer then it acts as a worst-case channel. If the
          parameter is between 0 and 1 then it acts as a random channel.

        - ``prob`` --  probability with which the additive white noise will
          be generated at any given time instance. This must be strictly in
          the interval [0,1]. The additive noise is generated for each
          frequency and time and so it will insert a 1 in the
          frequency-time matrix according to the frequency and time it
          corresponds to.  The prob may be provided with an instance of
          GeneralDiscreteDistribution. This ensures that the discrete
          distribution is not reinitialized at every run of this function.

    OUTPUT: a frequency-time matrix containing codeword + impulse noise
    + narrow-band noise + random noise

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

    """
    cdef int n, q
    cdef list all_ones, output

    n = len(word)
    q = word.q

    # This is the channel output after random noise is added. At this step
    # the output contains the codeword+noise in the ftmatrix form.
    output = binary_symmetric_channel(word.ftmatrix, prob)

    # Channel output after adding impulse_noise
    if 0 < impulse_noise and impulse_noise < 1:
        impulse_noise = [i for i,a in enumerate(random_noise(n, 2,
                        impulse_noise)) if a == 1]
    elif 0 < impulse_noise:
        impulse_noise = random_list(n, impulse_noise)
        # This is the channel output after impulse noise is added. The
        # column corresponding to the time gets all ones.
        all_ones = [1]*q
        output = zip(*output)
        for im in impulse_noise:
            output[im] = all_ones[:]
        output = zip(*output)


    # Channel output after adding impulse_noise
    if 0 < nb_noise and nb_noise < 1:
        nb_noise = [i for i,a in enumerate(random_noise(q, 2, nb_noise))
                    if a == 1]
    elif 0 < nb_noise:
        nb_noise = random_list(q, nb_noise)
        # This is the channel output after narrowband noise is added. The
        # row corresponding to the freq gets all ones.
        all_ones = [1]*n
        for nb in nb_noise:
            output[nb] = all_ones[:]

    # Some of the entries of output might be tuples. So, convert them back
    # to lists
    return map(list, output)


cpdef int decoder(ftmatrix, code, tx, algorithm=None):
    """
    def decoder(ftmatrix, code, tx=None):

    Example decoder for powerline channel which 'decodes'.
    This will call code.bounded_distance_decoder() or
    code.minimum_distance_decoder(). Change it inside the function if you
    need.
    INPUT:
        - ``ftmatrix`` -- the output of channel()
        - ``code`` -- the code from which tx is a :class:`Codeword`
        - ``tx`` -- the original :class:`Codeword` that was input to channel()
        - ``algorithm`` -- if it is not given then do a minimum distance decoding
            Possible values are 'bounded', 'minimum'
    OUTPUT:
        The number of symbols in error.

    If you want the decoded codeword then call code.<decoder function>
    directly.
    """
    # TODO: need some tests
    # The following line could be replaced by a smarter decoder that keeps
    # all the positions which are in error (not erasure)
    cw = ftmatrix_to_Codeword(ftmatrix)
    if algorithm is None:
        algorithm = 'minimum'
    if algorithm == 'minimum':
        dec_cw = code.minimum_distance_decoder(cw)
    elif algorithm == 'bounded':
        dec_cw = code.bounded_distance_decoder(cw, tx)
    else:
        raise ValueError("algorithm can take only two values: 'bounded' "
                "or 'minimum'")

    return hamming_distance(dec_cw, tx)


cpdef list find(vec, comp, int value):
    """
    Find the indices of vec(tor) which comp(are) to value. The vec may be
    a vector or a list or a tuple.

    INPUT:
        - ``vec`` -- a vector, tuple or list
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


cpdef Codeword ftmatrix_to_Codeword(ftmatrix):
    """
    ftmatrix_to_Codeword(ftmatrix)
    Convert a frequency-time matrix to a :class:`Codeword`.
    The steps followed are:
    1. If there is a list of all 1's in a particular frequency, then it is
       set to an all 0 list
    2. If there is a list of all 1's in a particular time, then it is set
       to an all 0 list
    3. If there is all 0 all time coordinate or more than one 1 in a time
       coordinate then it is set to q (an erasure symbol)
    4. Otherwise it returns the index of the 1 in a particular time.

    INPUT:
        - ``ftmatrix`` -- an ftmatrix

    OUTPUT: a :class:`Codeword` with symbols from {0,...,q} with q denoting
    an erasure symbol.

    EXAMPLES::

        sage: ft = [[1,0,0],[0,1,0],[1,1,1]]
        sage: ftmatrix_to_Codeword(ft)
        [1, 1, 3]
    """
    cdef int n, q
    cdef list all_zeros, impulse_err, narrow_band_err

    n = len(ftmatrix[0])
    q = len(ftmatrix)
    # Find all the symbols which have narrow band errors
    narrow_band_err = [e for e,l in enumerate(ftmatrix)
                        if len(find(l, '==', 1)) == n]

    # Find all the impulse noise errors
    impulse_err = [e for e,l in enumerate(zip(*ftmatrix))
                    if len(find(l, '==', 1)) == q]

    # Now to find the other errors. First delete the narrow and impulse err
    # For "safety" let's try to copy the lists instead of copying pointers
    all_zeros = [0]*n
    for e in narrow_band_err:
        ftmatrix[e] = all_zeros[:]
    ftmatrix = zip(*ftmatrix) # "transpose" the matrix into a time-freq one
    all_zeros = [0]*q
    for e in impulse_err:
        ftmatrix[e] = all_zeros[:]
    #********************************************************************#
    # Note that ftmatrix is now transposed. It is a time-freq matrix now #
    #********************************************************************#

    # Convert ftmatrix into a Codeword (this is a time-freq matrix now)
    cw = []
    Q = [q]
    for f in ftmatrix:
        e = find(f, '==', 1)
        if len(e) == 1:
            cw += e
        else:
            cw += Q
    return Codeword([c for c in cw], q)


cdef list random_list(l, int max_num, int exact_num=0):
    """
    INPUT:
        - ``l`` -- number of elements to randomize over
        - ``max_num`` -- max number of narrow-band or impulse noise present
        - ``exact_num`` -- If it is true then return exactly max_num sized list

    OUTPUT: the frequencies/times which have narrow-band/impulse noise

    EXAMPLES::

        sage: set_random_seed(0); random_list(18, 10)
        [9]
        sage: set_random_seed(0); random_list(18, 10, True)
        [2, 9, 0, 5, 3, 12, 1, 16, 10, 15]
    """
    cdef int num, a

    if max_num == 0: return []
    if exact_num:
        num = max_num   # Use this for exactly max_num non-zero rows
    else:
        num = randint(0, max_num) # The number of non-zero rows
    llist = []          # This will hold the list of generated symbols
    if num > 0:
        while len(llist) < num:
            a = randint(0, l-1)
            if a not in llist:
                llist += [a]
    return llist

cdef list random_noise(int n, int q, double prob=0.1):
    """
    random_noise(n, q, prob=0.1)
    Returns a vector with noise pattern from a symbol from {0,1,...,q-1}. If
    the symbol is 0, then there is no noise. This has probability 1-prob of
    occurring. The other alphabets occur with a probability of prob/(q-1).

    INPUT:
        - ``n`` -- the length of the vector
        - ``q`` -- the size of the alphabet
        - ``prob`` -- the probability that a coordinate of the vector is
          non-zero

    OUTPUT: a randomly generated :class:`Codeword` from {0,...,q-1}^n

    EXAMPLES::

        sage: random_noise(10, 2, 0.4)  # random
        [1, 0, 1, 0, 1, 0, 0, 0, 1, 0]
        sage: random_noise(10, 5, 0.4)  # random
        [1, 0, 0, 4, 0, 4, 0, 4, 0, 2]
    """
    P = [1-prob] + [prob/(q-1) for _ in range(q-1)]
    D = GeneralDiscreteDistribution(P)
    return [D.get_random_element() for _ in range(n)]




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
    cdef public object __dict__, _word, orig_word

    def __add__(self, w):
        cdef int q, i

        c = self._word
        w = w._word
        q = self.q
        # add elements mod q only if it is not q since q denotes an erasure
        return Codeword([(c[i]+w[i])%q if c[i] != q and w[i] != q else q
                        for i in range(len(c))], q)

    def __contains__(self, int b):
        return (b in self._word)

    # Need this so that lazy_attribute works in cython
    def __getattr__(self, attr):
        try:
            return self.__dict__[attr]
        except KeyError:
            raise AttributeError

    def __getitem__(self, int i):
        return self._word[i]

    def __init__(self, word, q=None):
        if q is not None:
            self.q      = q
        elif hasattr(word[0], 'parent') and word[0].parent().is_finite():
            self.q      = word[0].parent().order()
        else:
            raise ValueError("You must provide the alphabet size q when"
                            " the input word is not over a finite field")

        self.orig_word  = word
        self._word      = self._gen_word()
        self.__dict__   = {}

    def __iter__(self):
        for a in self._word:
            yield a

    def __len__(self):
        return len(self.orig_word)

    def __repr__(self):
        return str(self._word)

    # Need this so that lazy_attribute works in cython
    def __setattr__(self, attr, value):
        self.__dict__[attr] = value

    def __sub__(self, w):
        cdef int q, i

        c = self._word
        w = w._word
        q = self.q
        return Codeword([(c[i]-w[i])%q if c[i] != q and w[i] != q else q
                        for i in range(len(c))], q)

    # @lazy_attribute is needed to *not* compute the ftmatrix during
    # initialization. It will be computed only when explicitly required and
    # then the computation will be cached.
    @lazy_attribute
    def ftmatrix(self):
        cdef int w, q
        return [[1 if w == q else 0 for w in self._word]
                for q in xrange(self.q)]

    cdef _gen_word(self):
        word = self.orig_word
        if hasattr(word[0], 'parent') and word[0].parent().is_finite():
            return map(mylog, word)
        else:
            return list(word)


#***********************************************************************#
# classes not included in this file: NonlinearCode, PowerlineCode       #
#***********************************************************************#
