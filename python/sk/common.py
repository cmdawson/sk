# Some common universal gate sets
from scipy import mat
from math import cos,sin,sqrt,pi

"""This defines some common universal gate sets as dictionaries that you
can pass to net.generate()

They are:
	ht_net	- Hadamard and PI/8 gate set
	lps_net	- The very efficient Lubotsky, Phillips and Sarnak set

Different gate sets require you to enumerate to different sequence lengths 
before getting a decent covering of SU(2). A rough number for the HT set is about 14. For the LPS set is it around 6.
"""

# Hadamard and PI/8 
H = [ [1j/sqrt(2), 1j/sqrt(2)], [1j/sqrt(2), -1j/sqrt(2)] ]
T = [ [cos(pi/8)-sin(pi/8)*1j, 0], [0, cos(pi/8)+sin(pi/8)*1j] ]
ht_net = { 'H' : mat(H), 'T' : mat(T) }

# Lubotsky, Phillips and Sarnak
L = [ [1/sqrt(5), 2j/sqrt(5)], [2j/sqrt(5), 1/sqrt(5)] ]
P = [ [1/sqrt(5), 2/sqrt(5)], [-2/sqrt(5), 1/sqrt(5)] ]
S = [ [(1+2j)/sqrt(5), 0], [0, (1-2j)/sqrt(5)] ]
lps_net = { 'L' : mat(L), 'P' : mat(P), 'S' : mat(S) }
