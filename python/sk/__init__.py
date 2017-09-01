from scipy import mat
from sk_base import net


def pygenerate(self,D,maxdepth):
	d = {}
	for l,M in D.iteritems():
		if M.shape != (2,2):
			raise Exception("Error: matrix is not 2x2")
			return
		d[l] = [M[0,0], M[0,1], M[1,0], M[1,1]]

	self.generate_base(d,maxdepth);


def pyapproximate(self,M):
	if M.shape != (2,2):
		raise Exception("Error: matrix is not 2x2")
		return
	
	# convert the scipy mat into a list
	m = [M[0,0], M[0,1], M[1,0], M[1,1]]

	# then pass this list to the base function 
	return self.approximate_base(m);


def pyevaluate(self,s):
	# the wrapped evaluate returns a list which we turn into a scipy mat
	M = mat(self.evaluate_base(s))
	# then reshape it and return it
	M.shape = (2,2)
	return M
	

def pysolovay_kitaev(self,M,d):
	if M.shape != (2,2):
		raise Exception("Error: matrix is not 2x2")
		return
	# Start by converting the scipy mat to a list
	m = [M[0,0], M[0,1], M[1,0], M[1,1]]
	# then call the wrapped function
	K = self.solovay_kitaev_base(m,d)

	# K is now a list containing [the string length, the string, a00, a01,
	# a10, a11]. This needs to be tidied up a bit

	#print K

	R = K[0:2]
	R.append(mat(K[2:6]))
	R[2].shape = (2,2)

	return R

net.generate = pygenerate
net.approximate = pyapproximate
net.evaluate = pyevaluate
net.solovay_kitaev = pysolovay_kitaev
