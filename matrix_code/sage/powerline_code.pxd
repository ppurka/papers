#------------------------------------------------------------------------#
#                                                                        #
#      Copyright (C) 2013. P. Purkayastha and H. M. Kiah                 #
#      License: GPL-3 or later.                                          #
#                                                                        #
#------------------------------------------------------------------------#
from cpython cimport bool
cdef int equitable_distance(list cw, list ftmatrix)
cdef int decoder_cwcode(NonlinearCode cwcode, list f)
cdef int hamming_distance(x, y, int prevdist=*)
cdef int mylog(x, base=*)

cdef list add_ftmatrix(list f, list g)
cdef list binary_symmetric_channel(list word, prob=*)
cdef list dft(list vec)
cdef list find(list vec, str comp, int value)
cdef list idft(list vec)
cdef list random_noise(int n, int q, prob=*)
cdef list symbol_weight_distribution_vector(vector, q=*)

cpdef int equitable_decoder(list ftmatrix, PowerlineCode code, Codeword tx,
        bool detect_nb=*, bool randomize=*)
cpdef int decoder_mtfsk(list ftmatrix, code, Codeword tx,
        NonlinearCode cwcode, bool detect_nb=*, bool bdd=*,
        bool randomize=*)
cpdef int symbol_weight_matrix(code)
cpdef int symbol_weight_vector(vector, q=*)

cpdef list channel(Codeword word, impulse_noise=*, nb_noise=*, prob=*,
                   fade=*)
cpdef list channel_mfsk(Codeword cw, mfskchannel, additive_noise)



cdef class Codeword:
    cdef public int q
    cdef public list _word, ftmatrix
    cdef public object orig_word

    cdef list _ftmatrix(self, NonlinearCode cwcode)


cdef class NonlinearCode:
    cdef public int _d, d, n, q, is_cw1
    cdef public object _list
    cdef public object _code
    cdef public object M, size

    cpdef int _get_min_dist(self)
    cpdef int bounded_distance_decoder(self, list rx, list tx)
    cpdef minimum_distance(self)
    cpdef random_element(self)


cdef class PowerlineCode(NonlinearCode):
    cdef public NonlinearCode cwcode
    cpdef int _get_min_dist(self)

cdef class PowerlineLinearCode(PowerlineCode):
    cpdef random_element(self)

cdef class PowerlineReedSolomonCode(PowerlineLinearCode):
    cpdef int bounded_distance_decoder(self, list rx, list tx)

cdef class PowerlineReedSolomonSubcode(PowerlineReedSolomonCode):
    cdef public list constant_words
    cpdef random_element(self)

cdef class PowerlineReedSolomonCosetCode(PowerlineReedSolomonCode):
    cdef public object coset_leader
    cpdef int bounded_distance_decoder(self, list rx, list tx)
    cpdef random_element(self)
