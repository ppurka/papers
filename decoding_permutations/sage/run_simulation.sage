
# For the plot commands this file needs sage-5.2 along with the patches in
#  1. http://trac.sagemath.org/10512
#  2. http://trac.sagemath.org/13340

from sage           import *
from powerline_code import *
from simulation     import *


# Run this block just once
num_tx          = [10**3]*9 + [10**4]*9 + [10**4]*9
prob_list       = ( srange(9*10^-1, 10^-2, -10^-1)
                  + srange(9*10^-2, 10^-3, -10^-2)
                  + srange(9*10^-3, 10^-4, -10^-3)
                  )



# Simulation for DIM over GF(2)
set_random_seed(0)  # First set the random seed
num_symbol_err  = [0]*len(prob_list)

flag = 0
m = 1
F = GF(2**m, 'a')
C = BCHCode(7, 3, F)
d = C.minimum_distance()
n = C.length()
print "Parameters of BCH: [%d, %d, %d]_%d\n"%(n, C.dimension(), d, 2^m)
Clist = [Codeword(_)._word for _ in C]
M = F.order()**C.dimension()
randlist = [randint(0,M-1) for _ in xrange(max(num_tx))]
set_random_seed(0)  # Set the random seed again. This time for channel probabilities.

for i,p in enumerate(prob_list):
    if flag == 3:
        print "Reached 0 errors. Thrice. Skipping {0}".format(p)
        continue
    D = GeneralDiscreteDistribution([1-p, p])
    D.set_seed(0)
    for j in range(num_tx[i]):
        tx  = Clist[randlist[j]]        # get a random codeword
        #tx = Codeword(tx)._word        # get the Z_q representation
        txp = binary_to_permutations_by_flip(tx, m) # get permutations
        txp = [ZZ(_) for _ in txp]      # Make Codeword() happy below
        txp = Codeword(txp, n+1)        # make this a Codeword()

        rx  = channel(txp, impulse_noise=0, nb_noise=0, prob=D)
        num_symbol_err[i] += decoder_for_binary_DIM(rx, tx, m, d, Clist)

    num_symbol_err[i] /= num_tx[i]*n
    print (p, num_symbol_err[i]), " ",
    if num_symbol_err[i] == 0:
        flag += 1
        print
    else:
        flag = 0
    sys.stdout.flush()

print
SER_2_DIM = num_symbol_err[:]
print "SER_2_DIM = ", SER_2_DIM
pSER_2_DIM = list_plot_loglog(zip(prob_list, num_symbol_err), plotjoined=True)
pSER_2_DIM.show(xmin=max(prob_list), xmax=min(prob_list),
        gridlines='minor',
        axes_labels=['Probability', 'Symbol Error Rate'], frame=True,
        axes_pad=0)






# Simulation for DIM over GF(4)
set_random_seed(0)  # First set the random seed
num_symbol_err  = [0]*len(prob_list)

flag = 0
m = 2
F = GF(2**m, 'a')
C = BCHCode(7, 3, F)
d = C.minimum_distance()
n = C.length()
print "Parameters of BCH: [%d, %d, %d]_%d\n"%(n, C.dimension(), d, 2^m)
Clist = [Codeword(_)._word for _ in C]
M = F.order()**C.dimension()
randlist = [randint(0,M-1) for _ in xrange(max(num_tx))]
set_random_seed(0)  # Set the random seed again. This time for channel probabilities.

for i,p in enumerate(prob_list):
    if flag == 3:
        print "Reached 0 errors. Thrice. Skipping {0}".format(p)
        continue
    D = GeneralDiscreteDistribution([1-p, p])
    D.set_seed(0)
    for j in range(num_tx[i]):
        tx  = Clist[randlist[j]]        # get a random codeword
        #tx = Codeword(tx)._word        # get the Z_q representation
        txp = binary_to_permutations_by_flip(tx, m) # get permutations
        txp = [ZZ(_) for _ in txp]      # Make Codeword() happy below
        txp = Codeword(txp, m*n+1)      # make this a Codeword()

        rx  = channel(txp, impulse_noise=0, nb_noise=0, prob=D)
        num_symbol_err[i] += decoder_for_binary_DIM(rx, tx, m, d, Clist)

    num_symbol_err[i] /= num_tx[i]*n
    print (p, num_symbol_err[i]), " ",
    if num_symbol_err[i] == 0:
        flag += 1
        print
    else:
        flag = 0
    sys.stdout.flush()

print
SER_4_DIM = num_symbol_err[:]
print "SER_4_DIM = ", SER_4_DIM
pSER_4_DIM = list_plot_loglog(zip(prob_list, num_symbol_err), plotjoined=True)
pSER_4_DIM.show(xmin=max(prob_list), xmax=min(prob_list),
        gridlines='minor',
        axes_labels=['Probability', 'Symbol Error Rate'], frame=True,
        axes_pad=0)














# Simulation for DPM over GF(2)
set_random_seed(0)  # First set the random seed
num_symbol_err  = [0]*len(prob_list)

flag = 0
q = 2
F = GF(q, 'a')
C = BCHCode(7, 3, F)
d = C.minimum_distance()
n = C.length()
print "Parameters of BCH: [%d, %d, %d]_%d\n"%(n, C.dimension(), d, q)
Clist = [Codeword(_)._word for _ in C]
M = F.order()**C.dimension()
randlist = [randint(0,M-1) for _ in xrange(max(num_tx))]
set_random_seed(0)  # Set the random seed again. This time for channel probabilities.

for i,p in enumerate(prob_list):
    if flag == 3:
        print "Reached 0 errors. Thrice. Skipping {0}".format(p)
        continue
    D = GeneralDiscreteDistribution([1-p, p])
    D.set_seed(0)
    for j in range(num_tx[i]):
        tx  = Clist[randlist[j]]        # get a random codeword
        #tx = Codeword(tx)._word        # get the Z_q representation
        txp = vector_to_permutation(tx, q) # get permutations
        txp = [ZZ(_) for _ in txp]      # Make Codeword() happy below
        txp = Codeword(txp, 1+(q-1)*n)  # make this a Codeword()

        rx  = channel(txp, impulse_noise=0, nb_noise=0, prob=D)
        num_symbol_err[i] += decoder_for_DPM(rx, tx, q, d, Clist)

    num_symbol_err[i] /= num_tx[i]*n
    print (p, num_symbol_err[i]), " ",
    if num_symbol_err[i] == 0:
        flag += 1
        print
    else:
        flag = 0
    sys.stdout.flush()

print
SER_2_DPM = num_symbol_err[:]
print "SER_2_DPM = ", SER_2_DPM
pSER_2_DPM = list_plot_loglog(zip(prob_list, num_symbol_err), plotjoined=True)
pSER_2_DPM.show(xmin=max(prob_list), xmax=min(prob_list),
        gridlines='minor',
        axes_labels=['Probability', 'Symbol Error Rate'], frame=True,
        axes_pad=0)











# Simulation for DPM over GF(3)
set_random_seed(0)  # First set the random seed
num_symbol_err  = [0]*len(prob_list)

flag = 0
q = 3
F = GF(q, 'a')
C = BCHCode(8, 4, F)
d = C.minimum_distance()
n = C.length()
print "Parameters of BCH: [%d, %d, %d]_%d\n"%(n, C.dimension(), d, q)
Clist = [Codeword(_)._word for _ in C]
M = F.order()**C.dimension()
randlist = [randint(0,M-1) for _ in xrange(max(num_tx))]
set_random_seed(0)  # Set the random seed again. This time for channel probabilities.

for i,p in enumerate(prob_list):
    if flag == 3:
        print "Reached 0 errors. Thrice. Skipping {0}".format(p)
        continue
    D = GeneralDiscreteDistribution([1-p, p])
    D.set_seed(0)
    for j in range(num_tx[i]):
        tx  = Clist[randlist[j]]        # get a random codeword
        #tx = Codeword(tx)._word        # get the Z_q representation
        txp = vector_to_permutation(tx, q) # get permutations
        txp = [ZZ(_) for _ in txp]      # Make Codeword() happy below
        txp = Codeword(txp, 1+(q-1)*n)  # make this a Codeword()

        rx  = channel(txp, impulse_noise=0, nb_noise=0, prob=D)
        num_symbol_err[i] += decoder_for_DPM(rx, tx, q, d, Clist)

    num_symbol_err[i] /= num_tx[i]*n
    print (p, num_symbol_err[i]), " ",
    if num_symbol_err[i] == 0:
        flag += 1
        print
    else:
        flag = 0
    sys.stdout.flush()

print
SER_3_DPM = num_symbol_err[:]
print "SER_3_DPM = ", SER_3_DPM
pSER_3_DPM = list_plot_loglog(zip(prob_list, num_symbol_err), plotjoined=True)
pSER_3_DPM.show(xmin=max(prob_list), xmax=min(prob_list),
        gridlines='minor',
        axes_labels=['Probability', 'Symbol Error Rate'], frame=True,
        axes_pad=0)











# The values obtained and generating the plot #

#SER_2_DIM =  [61669/70000, 3338/4375, 5821/8750, 10233/17500, 1427/2800, 1167/2800, 4069/14000, 17/112, 177/5000,
#              1021/35000, 2229/100000, 143/8750, 4009/350000, 5183/700000, 1549/350000, 789/350000, 297/350000, 61/350000,
#              7/50000, 81/700000, 29/350000, 3/43750, 1/17500, 3/100000, 3/175000, 3/350000, 0]
#
#SER_4_DIM =  [13849/14000, 9577/10000, 9079/10000, 59113/70000, 13333/17500, 5783/8750, 18491/35000, 24849/70000, 561/4375,
#              74643/700000, 11897/140000, 3273/50000, 33601/700000, 1143/35000, 13943/700000, 227/21875, 1439/350000, 159/175000,
#              103/140000, 61/100000, 163/350000, 7/20000, 153/700000, 107/700000, 13/140000, 13/350000, 3/175000]
#
#SER_2_DPM =  [1, 1, 69997/70000, 69977/70000, 873/875, 69357/70000, 67677/70000, 15457/17500, 38167/70000,
#              334923/700000, 282939/700000, 4067/12500,171029/700000,14611/87500,67819/700000,2219/50000, 8867/700000, 897/700000,
#              587/700000, 39/70000, 3/10000, 59/350000, 17/175000, 13/350000, 1/175000, 1/175000, 0]
#
#SER_3_DPM =  [1, 1, 1, 1, 1, 1, 1, 39989/40000, 7811/8000,
#              771579/800000, 75759/80000, 14743/16000, 176553/200000, 658877/800000, 576763/800000, 433881/800000, 219457/800000, 15993/400000,
#              22051/800000, 14543/800000, 8681/800000, 2443/400000, 2363/800000, 243/200000, 29/100000, 1/32000, 0]
#
#
#plot_2_DIM = list_plot_loglog(zip(prob_list, SER_2_DIM), plotjoined=True,
#        marker='*', legend_label="$\\Pi_0, [7, 3, 4]_2$", color='magenta')
#plot_4_DIM = list_plot_loglog(zip(prob_list, SER_4_DIM), plotjoined=True,
#        marker='o', legend_label="$\\Pi_1, [7, 3, 4]_4$", color='blue')
#plot_2_DPM = list_plot_loglog(zip(prob_list, SER_2_DPM), plotjoined=True,
#        marker='d', legend_label="$\\Pi_2, [7, 3, 4]_2$", color='black')
#plot_3_DPM = list_plot_loglog(zip(prob_list, SER_3_DPM), plotjoined=True,
#        marker='s', legend_label="$\\Pi_3, [8, 3, 5]_3$", color='red')
#
#(plot_2_DIM + plot_4_DIM + plot_2_DPM + plot_3_DPM).show(
#    axes_labels=["Probability of background noise",
#                 "Symbol error and erasure rate"],
#    axes_pad=0,
#    fig_tight=True,
#    fontsize=14,
#    frame=True,
#    gridlines='minor',
#    legend_back_color='white',
#    legend_fancybox=True,
#    legend_font_size=14,
#    legend_handlelength=2,
#    legend_markerscale=1,
#    legend_numpoints=2,
#    title="Performance of the maps using BCH code $[n,k,d]_q$",
#    xmax=min(prob_list),
#    xmin=max(prob_list),
#    )



