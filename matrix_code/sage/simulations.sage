#------------------------------------------------------------------------#
#                                                                        #
#      Copyright (C) 2013. P. Purkayastha and H. M. Kiah                 #
#      License: GPL-3 or later.                                          #
#                                                                        #
#------------------------------------------------------------------------#
from sage.all import *
from powerline_code import *
import sys, gc
from mfsk import *
from noises import *


def run_simulations(Code, r, detect_nb=False, randomize=True, prob_others=0.05):
    """
    This runs the simulations with equitable symbol weight codes, and
    comparing them with cosets of Reed Solomon codes and minimum symbol
    weight codes.

    The decoding used is minimum distance decoding with the 'distance'
    defined as the one used in the paper titled "Importance of symbol
    equity in coded modulation for power line communications".

    INPUT:

    - ``Code`` -- the code whose performance is to be simulated.
    - ``r`` -- the symbol weight of the code ``Code``.
    - ``detect_nb`` -- boolean (default: False). Whether the decoder should
      detect the presence of narrowband noise.
    - ``randomize`` -- boolean (default: True). Whether the decoder should
      choose a random codeword if more than one codeword is at the same
      distance from the received word. If it is ``False``, then the decoder
      chooses the first codeword from the list of codewords that it finds.
    - ``prob_others`` -- float (default: 0.05). The probability of the
      background noise, fading and impulse noise.

    OUTPUT:

    There is no output returned by the function. It will print out the
    required information.
    """
    print Code,
    print " of symbol weight {0}, with detect_nb = {1}".format(r, detect_nb),
    print " and probability of other noise = {}".format(prob_others)
    print "Starting simulations"
    n = Code.n
    q = Code.q
    config.r = r
    config.max_nb = q

    set_random_seed(0)  # First set the random seed
    num_symbol_err  = [0]*len(prob_list)
    #prob_others = 0.05  # probability of all other noises
    Clist = Code.list()
    M = len(Clist)
    randlist = [randint(0,M-1) for _ in xrange(max(num_tx))]
    for i,p in enumerate(prob_list):
        config.nb_noise_global = [0]*q
        D = GeneralDiscreteDistribution([1-p, p])
        D.set_seed(0)
        Dothers = GeneralDiscreteDistribution([1-prob_others, prob_others])
        Dothers.set_seed(0)
        #Dothers = 0
        for j in range(num_tx[i]):
            tx = Clist[randlist[j]]        # get a random codeword
            rx = channel(tx, impulse_noise=Dothers,
                         nb_noise=D, prob=Dothers, fade=Dothers)
            num_symbol_err[i] += equitable_decoder(rx, Code, tx,
                                                   detect_nb=detect_nb,
                                                   randomize=randomize)
        num_symbol_err[i] /= num_tx[i]*n
        print (p, num_symbol_err[i]), ", ",
        sys.stdout.flush()
    print "\n----------------------- DONE ----------------------------\n\n"


def run_simulations_linear(Code, cwcode, r, detect_nb=False, bdd=True,
                    randomize=True, prev_Clist=False, max_nb=None):
    """
    This runs the simulations for the matrix codes introduced in the ISIT
    2013 paper.

    The decoding used here is bounded distance decoding for the
    outer codes and minimum distance decoding for the inner codes.


    INPUT:

    - ``Code`` -- the code whose performance is to be simulated.
    - ``cwcode`` -- the constant weight code that is the inner code of the
      concatenated code construction. This must be a ``NonlinearCode``
      instance.
    - ``r`` -- the symbol weight of the code ``Code``.
    - ``detect_nb`` -- boolean (default: False). Whether the decoder should
      detect the presence of narrowband noise.
    - ``bdd`` -- boolean (default: True). Whether the outer code should be
      decoded using bounded distance decoding.
    - ``randomize`` -- boolean (default: True). Whether the decoder should
      choose a random codeword if more than one codeword is at the same
      distance from the received word. If it is ``False``, then the decoder
      chooses the first codeword from the list of codewords that it finds.
    - ``prev_Clist`` -- boolean (default: False). Whether the random list
      of codewords should be regenerated, or the list of random codewords
      that were previously generated should be reused. This makes sense
      only if the same code is being simulated again.
    - ``max_nb`` -- integer (default: None). The maximum number of
      narrowband noise that should be introduced. If it is none, the
      maximum is set to ``m``, the number of frequencies being used.

    OUTPUT:

    There is no output returned by the function. It will print out the
    required information.
    """
    set_random_seed(0)  # First set the random seed
    output = []
    n = Code.n
    config.r = r
    config.max_nb = cwcode.n if max_nb is None else max_nb
    num_symbol_err  = [0]*len(prob_list)
    prob_others = 0.05  # probability of all other noises
    Dothers = GeneralDiscreteDistribution([1-prob_others, prob_others])
    print Code,
    print " of symbol weight {0}, with detect_nb = {1}".format(r, detect_nb),
    print " and max_nb = {0}".format(config.max_nb)
    print "Starting simulations"

    print "Generating random list of {0} codewords".format(max(num_tx)),
    sys.stdout.flush()
    global Clist
    if not prev_Clist:
        Clist = []
        gc.collect()
        Clist = [Code.random_element() for _ in xrange(max(num_tx))]
    print "Done."
    sys.stdout.flush()
    for i,p in enumerate(prob_list):
        set_random_seed(0)  # First set the random seed
        config.nb_noise_global = [0]*len(Clist[0].ftmatrix)
        D = GeneralDiscreteDistribution([1-p, p])
        D.set_seed(0)
        Dothers.set_seed(0)
        #Dothers = 0
        for j in range(num_tx[i]):
            tx = Clist[j]        # get a random codeword
            rx = channel(tx, impulse_noise=Dothers,
                         nb_noise=D, prob=Dothers, fade=Dothers)
            num_symbol_err[i] += decoder_mtfsk(rx, Code, tx, cwcode,
                                                   detect_nb=detect_nb,
                                                   bdd=bdd)
        num_symbol_err[i] /= num_tx[i]*n
        print "(%.3f, %.5f),"%(p, num_symbol_err[i]),
        sys.stdout.flush()
        output.append((p, num_symbol_err[i]))
    print "\n----------------------- DONE ----------------------------\n\n"
    return output


# The SNR in dB
EsbyNodB_min  = 0.5
EsbyNodB_max  = 4.5

# The MFSK channel parameters
df = 10800.0        # 10.8kHz
Ts = 1.0/1080       # 1/9 * T_AC/2 - this is the symbol duration
fs = 5.0*1080*10**2 # ~500kHz      - this is the sampling frequency

def run_simulations_mfsk(Code, r, detect_nb=False, randomize=True):
    """
    This runs the simulations for the MFSK channel that is simulated using
    cyclostationary noise. This is used in the revised version of the paper
    titled "Importance of symbol equity in coded modulation for power line
    communications".

    The decoding used is minimum distance decoding with the 'distance' as
    defined in the paper.


    INPUT:

    - ``Code`` -- the code whose performance is to be simulated.
    - ``r`` -- the symbol weight of the code ``Code``.
    - ``detect_nb`` -- boolean (default: False). Whether the decoder should
      detect the presence of narrowband noise.
    - ``randomize`` -- boolean (default: True). Whether the decoder should
      choose a random codeword if more than one codeword is at the same
      distance from the received word. If it is ``False``, then the decoder
      chooses the first codeword from the list of codewords that it finds.

    OUTPUT:

    There is no output returned by the function. It will print out the
    required information.
    """
    n = Code.n
    q = Code.q
    config.r = r
    dBnum = 9
    num_tx = [10**5]*dBnum

    num_symbol_err = [0]*dBnum
    Clist = Code.list()
    M = len(Clist)
    set_random_seed(0)
    randlist = [randint(0,M-1) for _ in xrange(max(num_tx))]

    mfskchannel = MFSKChannel([i*df for i in xrange(1,q+1)], Ts, fs=fs)

    for i,EsbyNodB in enumerate(xsrange(EsbyNodB_min, EsbyNodB_max, 0.5,
                                        include_endpoint=True,
                                        universe=float)):
        set_random_seed(0)
        np.random.seed(int(0))

        # Get the average noise std dev.
        EsbyNo = 10**(EsbyNodB/10)       # SNR in linear scale
        No = mfskchannel.Es/EsbyNo       # noise power spectral density (PSD)
        sigma = float(np.sqrt(No*fs/2))  # std dev of the noise.
        cyclonoise = CyclostationaryGaussianNoise(sigma, fs, Ts)

        for j in xrange(num_tx[i]):
            tx = Clist[randlist[j]]        # get a random codeword
            rx = channel_mfsk(tx, mfskchannel, cyclonoise)
            num_symbol_err[i] += equitable_decoder(rx, Code, tx,
                                                   detect_nb=detect_nb,
                                                   randomize=randomize)
            print "Computing {}/{}, errs = {}\r".format(j, num_tx[i],
                                                        num_symbol_err[i]),
            sys.stdout.flush()
        print ""
        num_symbol_err[i] /= num_tx[i]*n
        print (EsbyNodB, num_symbol_err[i])
        sys.stdout.flush()
    print "\n----------------------- DONE ----------------------------\n\n"




#------------------------------------------------------------------------#
#                                                                        #
#           Comparing ESW codes with MSW codes and RS codes              #
#                                                                        #
#------------------------------------------------------------------------#

print "Importing codes...",
sys.stdout.flush()
import codes
print "Done\n\n"
print "All the simulations are with randomized decoded selection"
print "in the case where multiple codewords are at the same"
print "distance from the received word."
sys.stdout.flush()

prob_list = ([3*10**-1, 2*10**-1, 10**-1]
            +srange(9*10**-2, 10**-3, -10**-2))
num_tx    = [10**5]*12

## Equitable code (7, 336, 5)_8 with symbol weight 1
run_simulations(codes.Code_n7d5q8r1, 1, detect_nb=True)
run_simulations(codes.Code_n7d5q8r1, 1, detect_nb=False)
## Reed Solomon [7, 2, 6]_8 with symbol weight 2
run_simulations(codes.Code_n7d6q8r2, 2, detect_nb=True)
run_simulations(codes.Code_n7d6q8r2, 2, detect_nb=False)
## Expurgated Reed Solomon [7, 3, 5]_8 with symbol weight 2, size 336
run_simulations(codes.Code_n7d5q8r2M336, 2, detect_nb=True)
run_simulations(codes.Code_n7d5q8r2M336, 2, detect_nb=False)
#
#
## Size 20160
## Equitable code of size 20160, symbol weight 1
run_simulations(codes.Code_n7d2q8r1, 1, detect_nb=True)
run_simulations(codes.Code_n7d2q8r1, 1, detect_nb=False)
## Reed Solomon coset of size 8**4, symbol weight 4
run_simulations(codes.Code_n7d4q8r4, 4, detect_nb=True)
run_simulations(codes.Code_n7d4q8r4, 4, detect_nb=False)
## Expurgated Reed Solomon code of size 20160, symbol weight 2
run_simulations(codes.Code_n7d3q8r2, 2, detect_nb=True)
run_simulations(codes.Code_n7d3q8r2, 2, detect_nb=False)
#
#
## Size 21120
## Equitable code of size 21120, symbol weight 1
run_simulations(codes.Code_n15d11q16r1, 1, detect_nb=True)
run_simulations(codes.Code_n15d11q16r1, 1, detect_nb=False)
## Expurgated Reed Solomon code of size 21120, symbol weight 3
run_simulations(codes.Code_n15d12q16r3, 3, detect_nb=True)
run_simulations(codes.Code_n15d12q16r3, 3, detect_nb=False)
## Reed Solomon coset of size 16**3, symbol weight 3
run_simulations(codes.Code_n15d13q16r3, 3, detect_nb=True)
run_simulations(codes.Code_n15d13q16r3, 3, detect_nb=False)

## Size 1000
## Equitable symbol weight code of size 1000
run_simulations(codes.Code_n11d6q10r2_esw, 2, detect_nb=True)
run_simulations(codes.Code_n11d6q10r2_esw, 2, detect_nb=False)
## Minimum symbol weight code of size 1000
run_simulations(codes.Code_n11d6q10r2, 2, detect_nb=True)
run_simulations(codes.Code_n11d6q10r2, 2, detect_nb=False)

# Size 51 -- simulation for multiple background probabilities.
num_tx    = [10**7]*12
for p in (0.1, 0.075, 0.05, 0.025, 0.01):
    # Equitable symbol weight code of size 51
    run_simulations(codes.Code_n25d24q17r2_esw, 2, detect_nb=True, prob_others=p)
    # Minimum symbol weight code of size 51
    run_simulations(codes.Code_n25d24q17r2, 2, detect_nb=True, prob_others=p)


#------------------------------------------------------------------------#
#       Size 51 -- simulation for CyclostationaryGaussianNoise           #
#------------------------------------------------------------------------#
# num_tx is set within the run_simulations_mfsk function
run_simulations_mfsk(codes.Code_n25d24q17r2_esw, 2, detect_nb=True)
run_simulations_mfsk(codes.Code_n25d24q17r2,     2, detect_nb=True)



#------------------------------------------------------------------------#
#                                                                        #
#           Comparing concatenated codes with RS codes                   #
#                                                                        #
#------------------------------------------------------------------------#

print "Importing codes...",
sys.stdout.flush()
import codes_linear as codes1
print "Done\n\n"
print "All the simulations are with randomized decoded selection"
print "in the case where multiple codewords are at the same"
print "distance from the received word."
sys.stdout.flush()

num_tx    = [10**5]*12

## Rates 0.53
CodeRS_n15k14d2r15_q9 = run_simulations_linear(codes1.CodeRS_n15k14d2r15_q9,
                                               codes1.C9_4_4, 15,
                                               detect_nb=False, bdd=True)
CodeRSC_n15k8d8r8_q16 = run_simulations_linear(codes1.CodeRSC_n15k8d8r8_q16,
                                               codes1.C16_2_1, 8,
                                               detect_nb=True, bdd=True)
# Rates 0.3x
CodeRS_n15k14d2r15_q13 = run_simulations_linear(codes1.CodeRS_n15k14d2r15_q13,
                                               codes1.C13_6_5, 15,
                                               detect_nb=False, bdd=True)
CodeRSC_n15k5d11r5_q16 = run_simulations_linear(codes1.CodeRSC_n15k5d11r5_q16,
                                               codes1.C16_2_1, 5,
                                               detect_nb=True, bdd=True)


