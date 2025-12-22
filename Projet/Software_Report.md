# Software Report

## Program 1: Minimum Absent Word (MAW)
### Naive Approach (MAW_naif.py)
Compute all MAWs from a list of subsequences. The MAWs have a length less than `kmax`.
The first part of the code consists of several small functions required to compute the MAWs.
Then, the MAWs are computed using the algorithm described in *"Scientific_Report.pdf"*, in the section **"Naive Algorithm"**.
###Efficient Approach (MAW_dict.py)
For every read and every k lesser than kmax, we go through the read and store every couple (w,l) in a dictionary where w is a k-mer and
l is the list of letters occuring before w in the read.

Then according to the algorithm described in the scientific report, we update a set of minimal absent words we found before writing them in the output file
###Radix Tree Approach (MAW_radix.py)
For every read, we insert nodes containing every kmax-mer and the letter appearing before it. In accordance with the scientific report, after that step, 
the tree obtained is such that reading a k-mer w from the root of the tree give a node whose set of letters is the set of letter preceding w in the read.

Then, we use a depth-first search to get the preceding letters needed in the algorithm without having to calculate them for every length k
## Program 2: Minimum p-Absent Word (pMAW)
### Naive Approach (pMAW_naif.py)
Compute all pMAWs from a list of subsequences. The pMAWs have a length less than `kmax`.
The first part of the code consists of several small functions required to compute the pMAWs.
Then, the pMAWs are computed using the algorithm described in *"Scientific_Report.pdf"*, in the section **"Naive Algorithm"**.
We start by selecting a set of candidates for pMAWs among extensions of k-mers, before checking among the candidates which are actually pMAWs

We then group the pMAWs by length using groupby and write it in the output file
### Efficient Approach (pMAW_dict.py)
We have a function candidates(kmax,p,w) which does the following:
    if is of length kmax we return w if and only if it is p-absent
    else we extend w left or right, and if the resulting word x has not been looked at yet (checked by storing visited words in a set)
         we return x if it is p-absent or we try to extend it as well by calling recursively candidates(kmax,p,x)
         
In the main function, we generate candidates from each MAW, using sets as usual, and we check for each candidate if it has a subword which is a candidate

Finally, we group the pMAWs by length using groupby and write it in the output file
