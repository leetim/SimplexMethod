import matplotlib.pyplot as plt
import scipy as sp
import numpy as np
import random as rnd
from scipy.optimize import leastsq, linprog

A = (np.random.rand(6, 8) - 0.5) * 20
delt = -np.min(A)
A_tild = A + delt
print A_tild
min_first = [np.min(i) for i in A_tild]
max_first = [np.max(i) for i in np.transpose(A_tild)]
e6 = np.ones(6)
e8 = np.ones(8)

T = np.transpose
M = np.matrix
X = linprog(e6, A_ub = -np.transpose(A_tild), b_ub = -e8)
Y = linprog(-e8, A_ub = A_tild, b_ub = e6)
print "X optimize"
print X
print "Y optimize"
print Y
print sum(X.x) - sum(Y.x)
print M(T(A_tild))*T(M(X.x))
print M(A_tild)*T(M(Y.x))


X = X.x/sum(X.x)
print X
Y = Y.x/sum(Y.x)
print Y

print np.matrix(X)*np.matrix(A)*np.transpose(np.matrix(Y))
