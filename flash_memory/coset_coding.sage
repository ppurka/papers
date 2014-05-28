#------------------------------------------------------------------------#
#                                                                        #
#      Copyright (C) 2014. P. Purkayastha and H. M. Kiah                 #
#      License: GPL-3 or later.                                          #
#                                                                        #
#------------------------------------------------------------------------#


from sage.all import *
from itertools import izip

# L: the total number of levels in the flash memory
# n: the block length of the code
# q: the alphabet size of the base field (considered as Z_q) of the code
# W: the linear subspace of F_q^n containing the all-one vector;
#    we encode into cosets F_q^n/W
# maxrewrites: make it large if a large number of rewrites are expected
global n, L, W, q, maxrewrites
maxrewrites = 80
numtimes = 1000


#======================================================================#
#                                                                      #
#                       All the encoders                               #
#                                                                      #
#======================================================================#

#----------------------------------------------------------------------#
#                       Encoder for Scheme A                           #
#----------------------------------------------------------------------#
def encode_schemeA(newword, currentstate):
    r"""
    Ideally, newword will be randomly generated like this

        ``newword = V.random_element()``

    And currentstate belongs to `ZZ^n`, where ``n`` is the dimension of
    ``V``.

    OUTPUT:

    - the new state or ValueError if it exceeds the number of levels ``L``
    """
    global L, W, n, q

    maxtempstate = max(currentstate) + q
    maxchanges = (q-1)*n+1
    newstate = currentstate # We don't need to copy [:] since we don't
                            # change any element manually.

    for w in W:
        cL = [ZZ(ci+wi) for ci,wi in izip(newword, w)]
        m = max(si-ci for si,ci in izip(currentstate, cL))
        if m == 0:
            return cL
        tempstate = [ci+m for ci in cL]
        m = max(tempstate)
        if m <= maxtempstate:        # if the maximum level < previous
            changes = sum(t - c for t,c in izip(tempstate, currentstate))
            if changes < maxchanges:
                newstate = tempstate
                maxtempstate = m
                maxchanges = changes
    
    if maxtempstate > L:
        raise ValueError('Maximum level exceeds L')
    return newstate



#----------------------------------------------------------------------#
#                       Encoder for Scheme B                           #
#----------------------------------------------------------------------#
def alpha(ci, si):
    r"""
    Return the new state of the current cell ``i`` given the current cell
    state value ``si`` and the new element ``ci`` that is to be written.
    This method is supposed to be used only with Scheme B.
    """
    global q,L
    # Note that ci <= q-1
    quo,rem = si.quo_rem(q) # si is in ZZ, so this should work
    for i in xrange(quo, ceil(L/q)+1):
        si_prime = i*q + ci
        if si_prime >= si:
            return si_prime

def encode_schemeB(newword, currentstate):
    r"""
    Return the new state given the current state vector ``currentstate``
    and the new codeword that is supposed to be written ``newword``.
    """
    global L, W, q, n

    maxtempstate = max(currentstate) + q
    maxchanges = (q-1)*n+1
    newstate = currentstate # We don't need to copy [:] since we don't
                            # change any element manually.

    for w in W:
        cL = [ZZ(ci+wi) for ci,wi in izip(newword, w)]
        tempstate = [alpha(ci,si) for ci,si in izip(cL, currentstate)]
        m = max(tempstate)
        if m <= maxtempstate:
            changes = sum(t - c for t,c in izip(tempstate, currentstate))
            if changes < maxchanges:
                newstate = tempstate
                maxtempstate = m
                maxchanges = changes
    
    if maxtempstate > L:
        raise ValueError('Maximum level exceeds L')
    return newstate



#----------------------------------------------------------------------#
#                       Encoder for Construction 18                    #
#----------------------------------------------------------------------#

def encode_perm(newword, currentstate):
    r"""
    Following the conventions in Jiang et.al, the word ``(c1,c2,...,cn)``
    encodes the state level ``(s1,s2,...,sn)`` if ``s_c1>s_c2>... >s_cn``.

    For example, a possible state (13,21,1) is encoded by (2,1,3) 

    Returns the new state, given ``newwword`` and ``currentstate``, or
    returns -1 if unable to do so.
    """
    global n,L
    
    newstate=currentstate[:]
    for i in xrange(1, n):
        clow  = newstate[newword[n-i]-1]
        chigh = newstate[newword[n-i-1]-1]
        if clow >= chigh: 
            newstate[newword[n-i-1]-1] = clow+1
            if (clow+1) > L:
                return -1
            
    return newstate


def encode_prefix(newprefix, currentstate):
    r"""
    Given a ``newprefix``, this encodes the 'best' possible state that
    yields a permutation with that ``newprefix``.
    """
    
    global n
    
    symbols = range(1, n+1)
    for i in newprefix:
        symbols.remove(i)

    maxL = max(currentstate) + n
    tempstate = currentstate[:]
    flag = 0
    
    for suffix in Permutations(symbols):
        newword  = newprefix + list(suffix)
        newstate = encode_perm(newword, currentstate)
        
        if newstate != -1:
            if max(newstate) < maxL:
                tempstate = newstate[:]
                maxL = max(newstate)
                flag = 1
    
    if flag:
        return tempstate
    else:
        return -1


def rand_perm(n, l, oldprefix):
    r"""
    Random Generator of a prefix.

    Outputs a random prefix of length ``l`` that is not equal to the
    current one ``oldprefix``.
    """
    while True:
        p = range(1, n+1)
        
        for i in xrange(l):
            t = int(random()*(n-i))
            p[i], p[i+t] = p[i+t], p[i]

        if p[:l] != oldprefix:
            return p[:l]




#----------------------------------------------------------------------#
#                       Encoder for FlipMin using Scheme B             #
#----------------------------------------------------------------------#
def encode_flipmin(newword, currentstate):
    r"""
    Return the new state vector given the new codeword ``newword`` and the
    current state vector ``currentstate``.
    """
    global L, W, q, n

    minflip = n
    zero_vec = [0]*n
    newstate = zero_vec[:]
    
    for w in W:
        cL = [(ci+wi+csi)%q for ci,wi,csi in izip(newword, w, currentstate)]
        if cL == zero_vec:
            continue

        numflips = sum(ZZ(cL[i]) for i in xrange(n))
        
        if numflips < minflip:
            newstate = [ZZ(ci)+ZZ(wi) for ci,wi in izip(currentstate, cL)]
            minflip = numflips

    if max(newstate) > L:
        raise ValueError('Maximum level exceeds L')

    return newstate




#======================================================================#
#                                                                      #
#                       All the simulation functions                   #
#                                                                      #
#======================================================================#


#----------------------------------------------------------------------#
#                       Simulation for Scheme A                        #
#----------------------------------------------------------------------#
def sim_schemeA(numtimes, V):
    r"""
    INPUT:

    - ``numtimes`` -- the number of times the simulation must be run
    - ``V`` -- the vector space `Z_q^n`

    OUTPUT:

    - a list of length ``maxrewrites`` that counts the number of times
      rewrites ``i`` happened at the index ``i``.
      
    """
    # This function shouldn't fail with an Exception
    global n, L, W, q

    Zq = V.base_ring()
    bigcount = 0
    set_random_seed(0)
    data=[0]*maxrewrites
    
    for _ in xrange(numtimes):
        currentstate = [0]*n
        count = 0
        while True:     # Keep on rewriting until it gets a ValueError
            curr_vec = vector(Zq, currentstate)
            try:
                while True:
                    v = V.random_element()
                    if v-curr_vec not in W:
                        break
                prevstate = currentstate
                currentstate = encode_schemeA(v, currentstate)
                count += 1
            except ValueError:
                break
        bigcount += count
        data[count] += 1

    print "For Scheme A:"
    print "Average number of rewrites for q={}, L={}, n={} is: {} or {}".format(
        q, L, n, bigcount/numtimes, float(bigcount/numtimes))

    return data


#----------------------------------------------------------------------#
#                       Simulation for Scheme B                        #
#----------------------------------------------------------------------#
def sim_schemeB(numtimes, V):
    r"""
    INPUT:

    - ``numtimes`` -- the number of times the simulation must be run
    - ``V`` -- the vector space `Z_q^n`

    OUTPUT:

    - a list of length ``maxrewrites`` that counts the number of times
      rewrites ``i`` happened at the index ``i``.
      
    """
    global n, L, W, q

    Zq = V.base_ring()
    bigcount = 0
    set_random_seed(0)
    data=[0]*maxrewrites
    
    for _ in xrange(numtimes):
        currentstate = [0]*n
        count = 0
        while True:     # Keep on rewriting until it gets a ValueError
            curr_vec = vector(Zq, currentstate)
            try:
                while True:
                    v = V.random_element()
                    if v-curr_vec not in W:
                        break
                prevstate = currentstate
                currentstate = encode_schemeB(v, currentstate)
                count += 1
            except ValueError:
                break
        bigcount += count
        data[count]+=1

    print "For Scheme B:"
    print "Average number of rewrites for q={}, L={}, n={} is: {} or {}".format(
        q, L, n, bigcount/numtimes, float(bigcount/numtimes))

    return data


#----------------------------------------------------------------------#
#                       Simulation for Rank Modulation                 #
#----------------------------------------------------------------------#
def sim_rm(numtimes, l):
    r"""
    INPUT:
    
    - ``numtimes`` -- number of times the simulation must be run
    - ``l`` -- the length of the prefix in Construction 18

    OUTPUT:

    - a list of length ``maxrewrites`` that counts the number of times
      rewrites ``i`` happened at the index ``i``.
      
    """
    global n, L
    
    bigcount = 0
    set_random_seed(0)
    data=[0]*maxrewrites
    
    for _ in xrange(numtimes):
        currentstate = [0]*n
        count = 0
        newprefix = 0
        
        while True:
            newprefix = rand_perm(n, l, newprefix)
            currentstate = encode_prefix(newprefix, currentstate)
            if currentstate == -1:
                break
            else:
                count += 1
                
        bigcount += count
        data[count] += 1

    print "For Construction 18:"
    print "Average number of rewrites for l={}, L={}, n={} is: {} or {}".format(
        l, L, n, bigcount/numtimes, float(bigcount/numtimes))

    # Answer is ~3.3 for n=20, and ~6.2 for n=6 for W=<j, (1..1,0..0)>
    return data


#----------------------------------------------------------------------#
#                       Simulation for FlipMin using Scheme B          #
#----------------------------------------------------------------------#
def sim_flipmin(numtimes, V):
    r"""
    INPUT:

    - ``numtimes`` -- the number of times the simulation must be run
    - ``V`` -- the vector space `Z_q^n`

    OUTPUT:

    - a list of length ``maxrewrites`` that counts the number of times
      rewrites ``i`` happened at the index ``i``.
      
    """
    # This function shouldn't fail with an Exception
    global n, L, W, q
    Zq = V.base_ring()
    bigcount = 0
    set_random_seed(0)
    data=[0]*maxrewrites
    
    for _ in xrange(numtimes):
        currentstate = [0]*n
        count = 0
        while True:     # Keep on rewriting until it gets a ValueError
            curr_vec = vector(Zq, currentstate)
            try:
                while True:
                    v = V.random_element()
                    if v-curr_vec not in W:
                        break
                prevstate = currentstate
                currentstate = encode_flipmin(v, currentstate)
                count += 1
            except ValueError:
                break
        bigcount += count
        data[count] += 1

    print "For FlipMin using Scheme B:"
    print "Average number of rewrites for q={}, L={}, n={} is: {} or {}".format(
        q, L, n, bigcount/numtimes, float(bigcount/numtimes))

    return data



#======================================================================#
#                                                                      #
#                       Simulation for Figure 2 in ISIT 2014           #
#                                                                      #
#======================================================================#

#----------------------------------------------------------------------#
#                       Some common parameters                         #
#----------------------------------------------------------------------#
l = 2
L = 16
n = 8
q = 3
V = Integers(q)**n

#----------------------------------------------------------------------#
#                       Simulate Construction 18                       #
#----------------------------------------------------------------------#
data1_rm = sim_rm(numtimes, l)

#----------------------------------------------------------------------#
#                       Simulate Scheme B with \delta = 2              #
#            W = <11111111, 11110000> = <11110000, 00001111>           #
#----------------------------------------------------------------------#
W = V.subspace(
    [
    [1]*4 + [0]*4,
    [0]*4 + [1]*4
    ]
    )
data1_qary_m2 = sim_schemeB(numtimes, V)

#----------------------------------------------------------------------#
#                       Simulate Scheme B with \delta = 4              #
#                 W = <10101010, 11001100, 11110000, 11111111>         #
#----------------------------------------------------------------------#
W = V.subspace(
    [
    [1, 0]*4,
    [1]*2 + [0]*2 + [1]*2 + [0]*2,
    [1]*4 + [0]*4,
    [1]*8
    ]
    )
data1_qary_m4 = sim_schemeB(numtimes, V)


#----------------------------------------------------------------------#
#                       Plotting Figure 2 of ISIT 2014                 #
#----------------------------------------------------------------------#
plot1 = [] # For Scheme B with \delta = 2
plot2 = [] # For Scheme B with \delta = 4
plot3 = [] # For Construction 18 with m = 2 (or l = 2, here)

for i in xrange(maxrewrites):
    y1 = data1_qary_m2[i]
    y2 = data1_qary_m4[i]
    y3 = data1_rm[i]
    
    if y1+y2+y3>0:
        plot1 += [(i, y1/numtimes)]
        plot2 += [(i, y2/numtimes)]
        plot3 += [(i, y3/numtimes)]

p1 = list_plot( plot1,
                plotjoined=True,
                marker='s',
                color='black',
                linestyle='--',
                legend_label=r'Scheme B, $\delta = 2$')
p2 = list_plot( plot2,
                plotjoined=True,
                marker='d',
                color='blue',
                linestyle='-',
                legend_label=r'Scheme B, $\delta = 4$')
p3 = list_plot( plot3,
                plotjoined=True,
                marker='o',
                color='red',
                linestyle='-.',
                legend_label=r'Construction 18, $m = 2$')

(p3 + p1 + p2).show(
        axes_pad=0,
        axes_labels=['Number of rewrites', 'Fraction of trials'],
        fontsize=16,
        frame=True,
        legend_back_color='white',
        legend_borderpad=1.0,
        legend_fancybox=True,
        legend_font_size=16,
        legend_handlelength=2.0,
        legend_loc=0,
        title=r"Comparing Scheme B with Construction 18",
        typeset='type1',
        xmin=0
        )

(p3 + p2 + p1).save(
        filename="./compare-all.eps",
        axes_pad=0,
        axes_labels=['Number of rewrites', 'Fraction of trials'],
        fontsize=16,
        frame=True,
        legend_back_color='white',
        legend_borderpad=1.0,
        legend_fancybox=True,
        legend_font_size=16,
        legend_handlelength=2.0,
        legend_loc=0,
        title=r"Comparing Scheme B with Construction 18",
        typeset='type1',
        xmin=0
        )





#======================================================================#
#                                                                      #
#                       Simulation for Figure 3 in ISIT 2014           #
#                                                                      #
#======================================================================#

#----------------------------------------------------------------------#
#                       Some common parameters                         #
#----------------------------------------------------------------------#
L = 16
n = 8
q = 3
V = Integers(q)**n

#----------------------------------------------------------------------#
#                       Simulate Scheme A & B with \delta = 1          #
#                               W = <11111111>                         #
#----------------------------------------------------------------------#
W = V.subspace([[1]*n])
data2_schemeB = sim_schemeB(numtimes, V)
data2_schemeA = sim_schemeA(numtimes, V)

#----------------------------------------------------------------------#
#                       Simulate Scheme B with \delta = 0              #
#                               W = <00000000>                         #
#----------------------------------------------------------------------#
W = V.subspace(V[0])
data2_schemeB_delta0 = sim_schemeB(numtimes, V)

#----------------------------------------------------------------------#
#        Some common parameters for second part of the plot            #
#----------------------------------------------------------------------#
L = 16
n = 8
q = 2
V = Integers(q)**n
W = V.subspace(
    [
    [1]*4 + [0]*4,
    [0]*4 + [1]*4
    ]    
    )

#----------------------------------------------------------------------#
#                Simulate Scheme B & FlipMin with \delta = 2           #
#               W = <11111111, 11110000> = <11110000, 00001111>        #
#----------------------------------------------------------------------#
data3_flipmin = sim_flipmin(numtimes, V)
data3_schemeB = sim_schemeB(numtimes, V)

#----------------------------------------------------------------------#
#                       Plotting Figure 3 of ISIT 2014                 #
#----------------------------------------------------------------------#
plot2A = []
plot2B = []
plot2C = []

for i in xrange(maxrewrites):
    y1 = data2_schemeA[i]
    y2 = data2_schemeB[i]
    y3 = data2_schemeB_delta0[i]
    
    if y1+y2>0:
        plot2A += [(i, y1/numtimes)]
        plot2B += [(i, y2/numtimes)]
        plot2C += [(i, y3/numtimes)]
        

p2A = list_plot(plot2A,
                plotjoined=True,
                marker='o',
                color='red',
                linestyle='--',
                legend_label=r'Scheme A, $\delta = 1, q=3$')
p2B = list_plot(plot2B,
                plotjoined=True,
                marker='d',
                color='magenta',
                linestyle='-',
                legend_label=r'Scheme B, $\delta = 1, q=3$')
p2C = list_plot(plot2C,
                plotjoined=True,
                marker='^',
                color='cyan',
                linestyle='-',
                legend_label=r'Scheme B, $\delta = 0, q=3$')


plot3A=[]
plot3B=[]

for i in xrange(maxrewrites):
    y1 = data3_flipmin[i]
    y2 = data3_schemeB[i]
    
    if y1+y2>0:
        plot3A += [(i, y1/numtimes)]
        plot3B += [(i, y2/numtimes)]



p3A = list_plot(plot3A,
                plotjoined=True,
                marker='*',
                color='black',
                linestyle='--',
                legend_label=r'FlipMin, $\delta = 2, q=2$')
p3B = list_plot(plot3B,
                plotjoined=True,
                marker='s',
                color='blue',
                linestyle='-',
                legend_label='Scheme B, $\delta = 2, q=2$')

(p2A + p2B + p2C + p3A + p3B).show(
        axes_pad=0,
        axes_labels=['Number of rewrites', 'Fraction of trials'],
        fontsize=16,
        frame=True,
        legend_back_color='white',
        legend_borderpad=1.0,
        legend_fancybox=True,
        legend_font_size=16,
        legend_handlelength=2.0,
        legend_loc=0,
        title=r"Comparing FlipMin, Scheme A, and Scheme B",
        #typeset='type1',
        xmin=0
        )
        
(p2A + p2B + p2C + p3A + p3B).save(
        filename="./compare-ab.eps",
        axes_pad=0,
        axes_labels=['Number of rewrites', 'Fraction of trials'],
        fontsize=16,
        frame=True,
        legend_back_color='white',
        legend_borderpad=1.0,
        legend_fancybox=True,
        legend_font_size=16,
        legend_handlelength=2.0,
        legend_loc=0,
        title=r"Comparing FlipMin, Scheme A, and Scheme B",
        typeset='type1',
        xmin=0
        )
