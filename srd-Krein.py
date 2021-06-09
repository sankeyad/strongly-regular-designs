# SRD-Krein.py
#
# Input is a text file with strongly regular design parameters. 
# The code performs swap(s) if necessary between SRGs and their complements, so that a_i > b_i
# which is required for the Krein condition test from Hobart: 
# Krein conditions for coherent configurations, Linear Algebra Appl., 1995, 499 -- 508.

# Then the srd and its complement are both tested.
#
# Usage: Input file should have two tab-separated rows per SRD in the form
#  n_i    k_i    lam_i    mu_i    rho_i    sig_i  n    k    lam    mu    r    s        S_i    a_i    b_i    N_i    P_i
#  python3.7 srd-Krein.py <infile> <outfile>


from sys import argv
import math

script, infile, outfile = argv


#___________________________________
# start of function gettable
	
def gettable(fromfile):   # reads list of srds from text file

	srds = open(fromfile)
	lines = []
	# extract a line from the text file of srgs
	lines.append(srds.readline().split('\t'))	 # reads string, parses it to list
	
	while lines[ len(lines)-1 ] != ['']:   # keep going until you read eof
		lines.append(srds.readline().split('\t') )

	srds.close		
	return(lines)  
	# end of function gettable

	
#___________________________________	

# function computes complement of the SRD
def complement(s):	  # takes a list of parameters
	#		n		k    lam    mu  	S		a				b					N		P
	return([ s[0], s[1], s[2], s[3], s[0]-s[4], s[0]-2*s[4]+s[5], s[0]-2*s[4]+s[6], s[1]-s[8], s[1]-s[7],
	       s[9], s[10], s[11], s[12], s[9]-s[13], s[9]-2*s[13]+s[14], s[9]-2*s[13]+s[15], s[10]-s[17], s[10]-s[16] ])
# end of complement function

#___________________________________


# function swaps srg1 for its complement
def swapg1(s):   # takes a list of parameters
	#		n		k			lam					mu					S	a		b		N			P
	return([ s[0], s[0]-s[1]-1, s[0]-2*s[1]+s[3]-2, s[0]-2*s[1]+s[2], s[4], s[5], s[6], s[4]-s[7]-1, s[4]-s[8],
	         s[9], s[10], s[11], s[12], s[13], s[15], s[14], s[16], s[17] ])
# end of swapg1 function

#___________________________________


# function swaps srg2 for its complement
def swapg2(s):   # takes a list of parameters
	#		n	k		lam 	mu		S	a		b	N		P
	return([ s[0], s[1], s[2], s[3],  s[4], s[6], s[5], s[7],  s[8],
	         s[9], s[9]-s[10]-1, s[9]-2*s[10]+s[12]-2, s[9]-2*s[10]+s[11], s[13], s[14], s[15], s[13]-s[16]-1, s[13]-s[17] ])
# end of swapg2 function

#___________________________________

# function to do Krein test

def Krein(n_1, k_1, S_1, a_1, b_1, r_1, n_2, k_2, S_2, a_2, b_2, r_2): # takes srd params from n1 through a2

	X = 1 + r_1**3/(k_1*k_1) - (r_1+1)**3/(n_1-k_1-1)**2
	Y = 1 + r_2**3/(k_2*k_2) - (r_2+1)**3/(n_2-k_2-1)**2
	Z = ( S_1 + a_1*r_2 - b_1*(r_2+1))**3
	W = ( 1/(S_1*S_2) - 1/( (n_1-S_1) * (n_2-S_2) ) )**2

	return( X*Y >= Z*W -0.0000000001)  # avoiding small roundoff error
# end of Krein test
	


#___________________________________



# begin main

# will read two lines from table:
#  n_i    k_i    lam_i    mu_i    rho_i    sig_i  n    k    lam    mu    r    s        S_i    a_i    b_i    N_i    P_i


params = gettable(infile)
results = open(outfile, 'w')

# get one pair of lines
i = 0
while i<len(params)-1:
	textrow1 = params[i]    # textrow is now row i of the input file
	textrow2 = params[i+1]

	row1 = []    # will strip out only the needed values n, k, S, a, r
	row2 = []
	for entry in [0,1,2,3,12,13,14,15,16]:   # this pulls out n, k, lam, mu, S, a, b, N, P
		row1 = row1 + [int( textrow1[entry] )]
		row2 = row2 + [int( textrow2[entry] )]

	# now rearrange to  _1 followed by _2
	srd = row1 + row2

	if srd[5] < srd[6]:   	# swap if a1 < b1
		srd = swapg1(srd)
	if srd[14] < srd[15]:   # swap if a2 < b2
		srd = swapg2(srd)
		
	# Now pass r1 and r2 into the Krein test; don't need lambda, mu, b

	ans1 = Krein(srd[0], srd[1], srd[4], srd[5], srd[6], srd[7]-srd[8], srd[9], srd[10], srd[13], srd[14], srd[15], srd[16]-srd[17])
	
	srd = complement(srd)
	ans2 = Krein(srd[0], srd[1], srd[4], srd[5], srd[6], srd[7]-srd[8], srd[9], srd[10], srd[13], srd[14], srd[15], srd[16]-srd[17])

	results.write(f"Set {math.floor(i/2 + 1)} and complement Krein test: {ans1} and {ans2}\n")
	
	i += 2

results.close()





