# generateGQ.py

# Computes the SRG parameters for strongly regular decompositions from 
# hemisystem in H(3,q^2), via the associated GQ(q^2,q)
# The max value of m is input at command line.

# Usage: python3.7 generateGQ.py <max m> <outfile>
# The parameters are computed up to number of vertices on the order of m^4.

#import math
from sys import argv
script, M, outfile = argv

from tabulate import tabulate

# strongly regular graph parameters: n is the order, k the valency, r and s eigenvalues.
def getSRG(n, k, r, s):  # function takes [n, k, r, s] and returns [n, k, lamda, mu]
	mu = k + r*s
	lam = mu + r + s
	
	return([ n, k, lam, mu ])

srgfile = open(outfile, 'w')
srgs = []
srds = []
chars = []

for m in range(1,int(M)+1):   # This calculates those on the order of t = 2m+1
	# formula for G_0
	t = 2*m+1
	n0 = (t**3+1)*(t+1)
	k0 = t*(t**2+1)
	lam0 = t-1
	mu0 = t**2+1
	r0 = t-1
	s0 = -t**2-1
	#srgs.append(getSRG(n, k, r, s))
	srgs.append([n0,k0,lam0,mu0,r0,s0])
	#srgs.append([n,n-k-1,n-2*k+mu-2,n-2*k+lam,-s-1,-r-1])
	
	# formula for G_1 = G_2
	n = (t**3+1)*(t+1)/2
	k = (t**2+1)*(t-1)/2
	lam = (t-3)/2
	mu = (( t-1)**2 )/2
	r = t-1
	s = -(t**2-t+2)/2
	srgs.append([n,k,lam,mu,r,s])
	# the complement
	l = n-k-1
	lamb = n-2*k+mu-2
	mub = n-2*k+lam
	rb = -s-1
	sb = -r-1
	srgs.append([n,l,lamb,mub,rb,sb])
	
	# formula for SRD parameters
	S = (t**2+1)*(t+1)/2
	a = ((t+1)**2)/2
	b = (t+1)/2
	N = t**2*(t+1)/2
	P = t*(t**2+1)/2
	srds.append([S,a,b,N,P])
	
	# character table
	alpha = S
	beta = t**2*(t**2-1)/2
	gamma = t*(t+1)/2
	z2 = (t**2+1)*(t-1)
	z3 = (t**2+1)*(t**2-t+1)/2
	chars.append([alpha,beta,gamma,z2,z3])
	
	

for row in range(0, len(srds)):
	# Gamma_0
	onesrg = srgs[3*row]
	srgfile.write('\t'.join(str(onesrg[i]) for i in range(0,6)) + '\n')
	# Gamma_1
	onesrg = srgs[3*row+1]
	srgfile.write('\t'.join(str(onesrg[i]) for i in range(0,6)) + '\n')
	# now the complement
	onesrg = srgs[3*row+2]
	srgfile.write('\t'.join(str(onesrg[i]) for i in range(0,6)) + '\n')
	onesrd = srds[row]
	srgfile.write('\t'.join(str(onesrd[i]) for i in range(0,5)) + '\n')
	onechar = chars[row]
	srgfile.write('\t'.join(str(onechar[i]) for i in range(0,5)) + '\n\n')  
	

srgfile.close()
