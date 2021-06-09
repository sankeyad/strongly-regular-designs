# generateSymplectic.py

# Computes the strongly regular graph parameters for g1 and g2 from the strongly regular
# decompositions of symplectic graphs. See 
# Haemers-Higman: Strongly regular graphs with strongly regular decomposition,
# Linear Algebra Appl., 114/115, 1989, 379--398.
# Seidel: On two-graphs, and Shult's characterization of symplectic and orthogonal geometries over GF(2)
# EUT report. WSK, Dept. of Mathematics and Computing Science, Technische Hogeschool Eindhoven, 1973.

# Number of vertices is on the order of 4^m; max m is input at command line.

# Usage: python3.7 generateSymplectic.py <max m> <outfile>


#import math
from sys import argv
script, M, outfile = argv

from tabulate import tabulate

def getSRG(n, k, r, s):  # function takes [n, k, r, s] and returns [n, k, lamda, mu]
	mu = k + r*s
	lam = mu + r + s
	
	return([ n, k, lam, mu ])

srgfile = open(outfile, 'w')
srgs = []

for m in range(2,int(M)+1):   # This calculates those on the order of n = 512
	# formula 1
	n = 2**(2*m-1) + 2**(m-1) -1
	k = 2**(2*m-2) + 2**(m-1) -2
	r = 2**(m-1) -1
	s = -2**(m-2)-1
	srgs.append(getSRG(n, k, r, s))
	
	# formula 2
	n = 2**(2*m-1) - 2**(m-1)
	k = 2**(2*m-2) -1
	r = 2**(m-2) -1
	s = -2**(m-1) -1
	srgs.append(getSRG(n, k, r, s))
	
	# formula 3
	if m > 2:
		n = 2**(2*m-1) - 2**(m-1) -1
		k = 2**(2*m-2) - 2**(m-1) -2
		r = 2**(m-2) -1
		s = -2**(m-1) -1
		srgs.append(getSRG(n, k, r, s))
	# formula 4
	n = 2**(2*m-1) + 2**(m-1)
	k = 2**(2*m-2) -1
	r = 2**(m-1) -1
	s = -2**(m-2) -1
	srgs.append(getSRG(n, k, r, s))
	

for row in range(0, len(srgs)):
	onerow = srgs[row]
	srgfile.write('\t'.join(str(onerow[i]) for i in range(0,4)) + '\n')

srgfile.close()

