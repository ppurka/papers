#------------------------------------------------------------------------#
#                                                                        #
#      Copyright (C) 2013. P. Purkayastha and H. M. Kiah                 #
#      License: GPL-3 or later.                                          #
#                                                                        #
#------------------------------------------------------------------------#
from itertools import izip
import numpy as np

class Noise:
    """docstring for Noise"""
    def __init__(self, stddev):
        self.stddev = float(stddev)

    def __add__(self, other):
        raise NotImplementedError

    def __radd__(self, other):
        raise NotImplementedError

    def __mul__(self, other):
        raise NotImplementedError

    def __rmul__(self, other):
        raise NotImplementedError

    def __call__(self, other):
        raise NotImplementedError


class AdditiveWhiteGaussianNoise(Noise):

    def __add__(self, other):
        return other + np.random.normal(0, self.stddev, len(other))

class CyclostationaryGaussianNoise(Noise):

    def __init__(self, stddev, fs, Ts):
        Noise.__init__(self, stddev)

        # The noise variance
        A = [0.23, 1.38, 7.17]
        T_AC = 1.0/60
        theta = [0, -6.0/180.0*np.pi, -35.0/180.0*np.pi]  # in radians
        nl = [0, 1.91, 1.57*10**5]
        dt = 1.0/float(fs)

        # Period of variance = T_AC/2, len=9*numt
        nt = np.arange(0, T_AC/2, dt, dtype='float')
        self.len_nt = len(nt)

        # len is integral multiple of signal len
        #print "length of noise = ", len(nt)

        # This is the standard deviation, done such that the average is
        # 1 over the time period T_AC/2
        self.sigma = np.sqrt(
                A[0] * np.ones_like(nt) # l = 0
              + A[1] * np.abs(np.sin(2*np.pi*nt/T_AC + theta[1]))**nl[1] # l = 1
              + A[2] * np.abs(np.sin(2*np.pi*nt/T_AC + theta[2]))**nl[2] # l = 2
                )

        # The noise is cyclostationary, so it must be considered like this
        self.cycle_pos = 0 # This will go up to a max of cycle_pos_max-1
        self.cycle_pos_max = np.floor(T_AC/2.0/Ts) # The max number of times
                                                   # the symbol period occurs
                                                   # within T_AC/2
        self.num_samples_Ts = np.ceil(Ts * fs)     # From 0 to floor(Ts * fs)


        # The impulse response
        a = 1.2 * 10**-5
        nht = np.append(-nt[::-1], nt[1:])         # 2*nt-1 is total length

        # This is sampled at nht points with frequency fs, and hence it is
        # truncated at frequency fs/2.
        self.ht = np.sqrt(a/2)*a/(a**2/4 + 4* np.pi**2 * nht**2)
        self.ht *= np.hanning(len(self.ht))       # Apply Hanning window
        self.ht /= np.sqrt(np.sum(self.ht**2))    # Make sure ht has energy 1

        # Remove objects from memory
        del A
        del a
        del dt
        del nl
        del nht
        del nt
        del T_AC
        del theta

    def __add__(self, other):
        if self.cycle_pos == 0:
            self.generate_noise()

        pos = self.cycle_pos
        nTs = self.num_samples_Ts
        self.cycle_pos = np.mod(self.cycle_pos + 1, self.cycle_pos_max)
        return other + self.normal[pos*nTs:(pos+1)*nTs]

    def generate_noise(self):
        self.normal = (self.stddev *
                       self.sigma * np.random.normal(0, 1, len(self.sigma))
                      )
        lnt = self.len_nt
        self.normal = np.convolve(self.ht, self.normal)[lnt:2*lnt]
