from math import sqrt,pi,cos,sin
from scipy import mat, linalg
from random import random
from sk import net,utils

# Set up some 'matrices'
H = [ [1j/sqrt(2), 1j/sqrt(2)], [1j/sqrt(2), -1j/sqrt(2)] ]
T = [ [cos(pi/8)-sin(pi/8)*1j, 0], [0, cos(pi/8)+sin(pi/8)*1j] ]

# Prepare the dictionary matching a label with a 'matrix'
d = { 'H' : mat(H), 'T' : mat(T) }

# Instanciate the net with a tile width of 0.15, and generate all instruction
# sequences of length up to 15
N = net(0.15)
print "... generating net ..."
N.generate(d,15);

# Optionally you could save it into a file
# N.save("ht15.net")
# and then load it later via
# N = net("ht15.net)

# Generate a random matrix in SU(2)
M = utils.random_su2()
print "M = ", M

# optain a Solovay-Kitaev approximation of depth 3 to this matrix
# This function returns a list which contains the sequence length, the 
# actual sequence, and the matrix this sequence evaluates to
A = N.solovay_kitaev(M,3)

print
print "M is approximated by the instruction sequence"
print A[1]

print "(length", A[0], ")"

# Calculate the accuracy
print "The distance between these two is"
print utils.proj_dist(M,A[2])

# You can also evaluate the instruction sequence to ensure that it is
# what it's supposed to be
# MA = N.evaluate(A[1])
# print MA
# print utils.proj_dist(MA,A[2])	# (this should be zero)
