# strongly-regular-designs
Code for generating feasible parameter sets for strongly regular designs.

# SRD-onef.py
Reads a single text file of strongly regular graph parameters. Sample input in the file 50.txt.
Output is feasible parameter sets for strongly regular designs admitting
strongly regular decomposition, in which two graphs from the input file are
the SRGs on the two fibres of the SRD.

Usage: python3.7 SRD-onef.py <infile> <outfile>
  
The infile should be a tab-separated text file with one SRG parameter set in the form 'n k lambda mu' per line, e.g.:
204	63	22	18       
204	28	2	4
205	96	50	40
205	68	15	26
208	75	30	25
208	81	24	36

  
# SRD-twof.py
Same as above, except that there are two distinct input files. Use, for example
when the first file has orders between 101 and 150, and the second between 151 and 200. Sample input
in the files 50.txt and 100.txt.

Usage: python3.7 SRD-twof.py <infile1> <infile2> <outfile>
  
  
# SRD-Krein.py
Input is a text file with strongly regular design parameters. Sample input in the file krein-tester.txt
The code performs swap(s) if necessary between SRGs and their complements, so that a_i > b_i
which is required for the Krein condition test from Hobart: 
Krein conditions for coherent configurations, Linear Algebra Appl., 1995, 499 -- 508.

Then the srd and its complement are both tested.

Usage: python3.7 srd-Krein.py <infile> <outfile>
  
Input file should have two tab-separated rows per SRD in the form
  n_i    k_i    lam_i    mu_i    rho_i    sig_i  n    k    lam    mu    r    s        S_i    a_i    b_i    N_i    P_i

  
  
# generateSymplectic.py
Computes the strongly regular graph parameters for g1 and g2 from the strongly regular
decompositions of symplectic graphs. See 
Haemers-Higman: Strongly regular graphs with strongly regular decomposition,
Linear Algebra Appl., 114/115, 1989, 379--398.
Seidel: On two-graphs, and Shult's characterization of symplectic and orthogonal geometries over GF(2)
EUT report. WSK, Dept. of Mathematics and Computing Science, Technische Hogeschool Eindhoven, 1973.

Number of vertices is on the order of 4^m; max m is input at command line.

Usage: python3.7 generateSymplectic.py <max m> <outfile>  
  

# generateGQ.py
Computes the SRG parameters for strongly regular decompositions from hemisystem in H(3,q^2), via the associated GQ(q^2,q)
The max value of m is input at command line.

Usage: python3.7 generateGQ.py <max m> <outfile>
  
The parameters are computed up to number of vertices on the order of m^4.
