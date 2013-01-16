#------------------------------------------------------------------------#
#                                                                        #
#      Copyright (C) 2013. P. Purkayastha and H. M. Kiah                 #
#      License: GPL-3 or later.                                          #
#                                                                        #
#------------------------------------------------------------------------#
from sage.all import *
from powerline_code import *
import sys

print "Importing codes...",
sys.stdout.flush()
import codes
print "Done\n\n"
print "All the simulations are with randomized decoded selection"
print "in the case where multiple codewords are at the same"
print "distance from the received word."
sys.stdout.flush()

prob_list = (srange(3*10**-1, 10**-2, -10**-1)
            +srange(9*10**-2, 10**-3, -10**-2))
num_tx    = [10**5]*len(prob_list)
#prob_list       = [10**-2]
#num_tx          = [100]


def run_simulations(Code, r, detect_nb=False, randomize=True):
    print Code,
    print " of symbol weight {0}, with detect_nb = {1}".format(r, detect_nb)
    print "Starting simulations"
    n = Code.n
    q = Code.q
    config.r = r
    config.max_nb = q

    set_random_seed(0)  # First set the random seed
    num_symbol_err  = [0]*len(prob_list)
    prob_others = 0.05  # probability of all other noises
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



# Equitable code (7, 336, 5)_8 with symbol weight 1
run_simulations(codes.Code_n7d5q8r1, 1, detect_nb=True)
run_simulations(codes.Code_n7d5q8r1, 1, detect_nb=False)
# Reed Solomon [7, 2, 6]_8 with symbol weight 2
run_simulations(codes.Code_n7d6q8r2, 2, detect_nb=True)
run_simulations(codes.Code_n7d6q8r2, 2, detect_nb=False)
# Expurgated Reed Solomon [7, 3, 5]_8 with symbol weight 2, size 336
run_simulations(codes.Code_n7d5q8r2M336, 2, detect_nb=True)
run_simulations(codes.Code_n7d5q8r2M336, 2, detect_nb=False)


# Size 20160
# Equitable code of size 20160, symbol weight 1
run_simulations(codes.Code_n7d2q8r1, 1, detect_nb=True)
run_simulations(codes.Code_n7d2q8r1, 1, detect_nb=False)
# Reed Solomon coset of size 8**4, symbol weight 4
run_simulations(codes.Code_n7d4q8r4, 4, detect_nb=True)
run_simulations(codes.Code_n7d4q8r4, 4, detect_nb=False)
# Expurgated Reed Solomon code of size 20160, symbol weight 2
run_simulations(codes.Code_n7d3q8r2, 2, detect_nb=True)
run_simulations(codes.Code_n7d3q8r2, 2, detect_nb=False)


# Size 21120
# Equitable code of size 21120, symbol weight 1
run_simulations(codes.Code_n15d11q16r1, 1, detect_nb=True)
run_simulations(codes.Code_n15d11q16r1, 1, detect_nb=False)
# Expurgated Reed Solomon code of size 21120, symbol weight 3
run_simulations(codes.Code_n15d12q16r3, 3, detect_nb=True)
run_simulations(codes.Code_n15d12q16r3, 3, detect_nb=False)
# Reed Solomon coset of size 16**3, symbol weight 3
run_simulations(codes.Code_n15d13q16r3, 3, detect_nb=True)
run_simulations(codes.Code_n15d13q16r3, 3, detect_nb=False)


# Size 51
# Equitable symbol weight code of size 51
run_simulations(codes.Code_n25d24q17r2_esw, 2, detect_nb=True)
run_simulations(codes.Code_n25d24q17r2_esw, 2, detect_nb=False)
# Minimum symbol weight code of size 51
run_simulations(codes.Code_n25d24q17r2, 2, detect_nb=True)
run_simulations(codes.Code_n25d24q17r2, 2, detect_nb=False)



# Size 1000
# Equitable symbol weight code of size 1000
run_simulations(codes.Code_n11d6q10r2_esw, 2, detect_nb=True)
run_simulations(codes.Code_n11d6q10r2_esw, 2, detect_nb=False)
# Minimum symbol weight code of size 1000
run_simulations(codes.Code_n11d6q10r2, 2, detect_nb=True)
run_simulations(codes.Code_n11d6q10r2, 2, detect_nb=False)
