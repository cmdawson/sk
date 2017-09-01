from scipy import mat,linalg
from random import random
from math import sqrt

"""Bunch of useful SU(2) functions"""

ID = mat("1 0; 0 1")
X = mat("0 1; 1 0");
Y = mat("0 -1j; 1j 0");
Z = mat("1 0; 0 -1");

def random_su2():
	"""Returns a random matrix in SU(2)"""
	a = [random() for i in range(4)]
	nn = sqrt(sum(x*x for x in a))
	a = [x/nn for x in a]

	return a[0]*ID + 1j*a[1]*X + 1j*a[2]*Y + 1j*a[3]*Z

def proj_dist(A,B):
	"""Returns the 'projective distance' between of two matrices A,B"""
	return min(linalg.norm(A-B), linalg.norm(A+B))
