from sage.all import *
nb_noise_global = None
r = 1
max_nb = 3
F = None
Rx = None
x = None
one = None
zero = None
omega = None
Fstarsorted = []
Fstarlen = 0
"""
F = GF(16, 'a')
Rx = PolynomialRing(F, 'x')
x = Rx.gen()
one = F(1)
zero = F(0)
omega = F.multiplicative_generator()
Fstarlen = F.order()-1
Fstarsorted = [omega**i for i in range(Fstarlen)]
"""
