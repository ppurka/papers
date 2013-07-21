#------------------------------------------------------------------------#
#                                                                        #
#      Copyright (C) 2013. P. Purkayastha and H. M. Kiah                 #
#      License: GPL-3 or later.                                          #
#                                                                        #
#------------------------------------------------------------------------#
import numpy as np
from itertools import izip

class MFSKChannel:
    """
    - ``freqs`` -- list. The list of frequencies for MFSK.
    - ``Ts`` -- time in seconds. The time period for each symbol.
    - ``Es`` -- energy (default: 1). The energy of each symbol.
    - ``fs`` -- sampling frequency in Hz (default 500kHz).
    """

    def __init__(self, freqs, Ts, Es=1, fs=5*10**5):
        self.Es = float(Es)
        self.freqs = np.array(freqs, dtype='float')
        self.fs = float(fs)
        self.Ts = float(Ts)

        dt = 1.0/self.fs
        self._tsamples = np.arange(0, Ts, dt, dtype='float')
        self._cos = [np.sqrt(2/self.Ts) * dt *
                        np.cos(2 * np.pi * f * self._tsamples)
                     for f in self.freqs]
        self._sin = [np.sqrt(2/self.Ts) * dt *
                        np.sin(2 * np.pi * f * self._tsamples)
                     for f in self.freqs]


    def channel(self, transmitted, additive_noise, narrowband_noise):
        """
        INPUT:

        - ``transmitted`` -- The output of :meth:`modulate`.
        - ``additive_noise`` -- class. Class representing the additive
          noise.
        - ``narrowband_noise`` -- class. Class representing the narrowband
          noise.

        OUTPUT:

        - a generator for the modulated signals with added noise.

        """
        return narrowband_noise + (additive_noise + transmitted)

    def demodulate(self, received):
        """
        INPUT:

        - ``received`` -- The output of :meth:`channel`.

        OUTPUT:

        - a list taking values from `\{0,1\}`. The entry of any position in the
          list corresponds to a frequency and takes the value `0` if the energy
          in that corresponding frequency is below ``Es/2``; otherwise it takes
          the value `1`.

        """
        r = []
        for i,f in enumerate(self.freqs):
            rc = np.sum(received * self._cos[i])
            rs = np.sum(received * self._sin[i])
            r.append(1 if rc**2 + rs**2 >= self.Es/4 else 0)

        return r

    def modulate(self, transmit):
        """
        INPUT:

        - ``transmit`` -- list. The list with entries from `\{0,1\}` where
          a `1` indicates the frequency that is being transmitted.

        OUTPUT:

        - the modulated signal. Sampled with frequency ``self.fs``.

        """
        return sum(
                   np.sqrt(2*self.Es/self.Ts) *
                   np.cos(2*np.pi*f *
                       self._tsamples +
                       np.random.uniform(-np.pi, np.pi))
                   for f,tx in izip(self.freqs, transmit)
                   if tx != 0
               )

