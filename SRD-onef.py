# SRD-onef.py
# Reads a single text file of strongly regular graph parameters.
# Output is feasible parameter sets for strongly regular designs admitting
# strongly regular decomposition, in which two graphs from the input file are
# the SRGs on the two fibres of the SRD.

# Usage: python3.7 SRD-onef.py <infile> <outfile>
# The infile should be a tab-separated text file with one SRG parameter set in the form 'n k lambda mu' per line, e.g.:
# 204	63	22	18       
# 204	28	2	4
# 205	96	50	40
# 205	68	15	26
# 208	75	30	25
# 208	81	24	36


#import time
#begin_time = time.perf_counter()

from sys import argv

script, infile, outfile = argv

from tabulate import tabulate

import math

def is_integer_num(n):
    if isinstance(n, int):
        return True
    if isinstance(n, float):
        return n.is_integer()
    return False
    
def check_4_srd(graph1, graph2):  # checks a given pair of srgs
	n1 = graph1[0]
	k1 = graph1[1]
	lam1 = graph1[2]
	mu1 = graph1[3]
	# compute r1, s1
	c = lam1 - mu1
	delta = c**2 + 4*(k1 - mu1) 
	d = math.sqrt(delta)
	r1 = ( c + d )/2
	s1 = ( c - d )/2				
    
	n2 = graph2[0]
	k2 = graph2[1]
	lam2 = graph2[2]
	mu2 = graph2[3]
	# compute r2, s2
	c = lam2 - mu2
	delta = c**2 + 4*(k2 - mu2) 
	d = math.sqrt(delta)
	r2 = ( c + d )/2
	s2 = ( c - d )/2
	
	table = []
   
	if k2 < k1:
		S1 = 2 + k1 - k2
	else:
		S1 = 2
	n = n1 + n2
	while S1 < n1: 

		S2 = k2 + S1 - k1
		k0 = k1 + S2
	
		for a1 in range(0, S1): 
			a2 = a1 + lam2 - lam1
			lam0 = a1 + lam2
			if a2 >= 0 and lam0 < k0:  # then carry on	
			
				if a1 < S1*S2/n2:
					bot = a1+1
					top = S1
				else:
					bot = 0
					top = a1
				for b1 in range(bot, top):  # b1 < a1 in DGH but need to allow b1 < S1.
					b2 = b1 + mu2 - mu1
					mu0 = b1 + mu2
	
					c = lam0 - mu0
					delta = c**2 + 4*(k0 - mu0) 
					d = math.sqrt(delta)
					r = ( c + d )/2
					s = ( c - d )/2
					f = (( n-1 )*(-s) - k0)/(r-s)
					# check these parameters first
			
					if b2 == a2 or b2 < 0 or n-2*k0+mu0-2 < 0 or n-2*k0+lam0 < 0:
						continue #any of those would be invalid
					elif (n-k0-1)*mu0 != k0*(k0-lam0-1):
						continue
					elif r*s*s - 2*r*r*s - r*r - k0*r + k0*s*s + 2*k0*s < 0 or r*r*s - 2*r*s*s - s*s - k0*s + k0*r*r + 2*k0*r < 0:
						# above checks Krein parameters on gamma0
						continue 
					elif not is_integer_num(f):  # checks integrality of eigenvalue multiplicities
						continue 
					else: # carry on
						N1 = a2*k1/S2
						N2 = a1*k2/S1
						if is_integer_num(N1) and is_integer_num(N2) and n1 != S1 and n2 !=S2:  # then good
							P1 = ( (k1-N1)*S1 )/(n1-S1)  # check it's an integer
							P2 = ( (k2-N2)*S2 )/(n2-S2)
							if is_integer_num(P1) and is_integer_num(P2):
				
								rho1 = N1-P1
								rho2 = N2-P2
								sig1 = -(S2-b2)/(a2-b2)
								sig2 = -(S1-b1)/(a1-b1)
								if (rho1 == r1 and sig1 == s1) or (rho1 == s1 and sig1 == r1):
									# then all is good, keep going
									if (rho2 == r2 and sig2 == s2) or (rho2 == s2 and sig2 == r2):
										# print out parameters
										
										table.append([n1,k1,lam1,mu1,rho1,sig1,n,k0,lam0,mu0,r,s,S1,a1, b1, N1, P1])
										table.append([n2,k2,lam2,mu2,rho2,sig2," "," "," "," "," "," ",S2,a2, b2, N2, P2])

									else:
										continue
										# print(S1, "failed at eigenvalue stage")
								else: 
									continue
									# print(S1, "failed at eigenvalue stage")
							else:
								continue
								#print("P1 or P2 not integer")
						else:
							continue
							#print("N1 or N2 not integer")
			else:
				continue
				#print("a2 or lam0 problem")
					
		S1 += 1
	return(table)
	# end of function check_4_srd
	
	
#___________________________________	
# start of main
graphs = open(infile)
params = open(outfile, 'w')	

results = []   # will be table of feasible parameters for srds arising from a pair of srgs
g1 = [0,0,0,0]
g2 = [0,0,0,0]
lines = []
# first, extract a pair of lines from the text file of srgs
i = 1
j = 1
count = 0
lines.append(graphs.readline().split('\t'))	
	
while lines[ len(lines)-1 ] != ['']:
	lines.append(graphs.readline().split('\t') )
count = len(lines)-1   # this is because last item in list is now ''

for i in range(0,count):
	for j in range(i, count):
		line1 = lines[i]  # reads the text, separated by spaces, into list
		line2 = lines[j]
		for entry in range(0,4): # converting text to integer
			g1[entry] = int( line1[entry] )
			g2[entry] = int( line2[entry] )
			
		newones = check_4_srd(g1,g2)
		# new sets are then appended to the results table
		results += newones

# Loops through the input file to compare all pairs, including g2=g1.

# Uses tabulate to write results with a header row.
params.write(tabulate(results, headers = 
["n_i", "k_i", "lam_i", "mu_i","rho_i", "sig_i", "n", "k", "lam", "mu", "r", "s", "S_i", "a_i", "b_i", "N_i", "P_i"]))
params.close()
graphs.close()
# print( time.perf_counter() - begin_time )