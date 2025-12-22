# Readme
## Intro
This project contains programs which generate minimal absent words and minimal p-absent words for a sequence of fasta reads
## Dependencies
For the CPP version (only available for the naive MAW generation), you need to download : https://github.com/seqan/seqan
And compile with -I/path/seqan/include

For the Python version, the necessary packages are :
numpy, collections, xopen, time, itertools, os, csv, logging, sys

## Files
There are 8 different python files, here are their uses
- "radix_config.py", "radix_trees.py" and "readfa.py" are source files and don't need to be used
- "MAW_naive.py", "MAW_dict.py" and "MAW_radiw.py" implement minimal absent word generation respectively with a naive algorithm, an efficient algorithm, or an algorithm using radix trees
- "pMAW_naive.py" and "pMAW_dict.py" implement minimal p-absent word generation respectively with a naive algorithm and an efficient algorithm

## How to use
For every file, once the program is run, the user needs to run the function "read_files()" in order to load the data - If the user wants to load new data, they need to modify that function inside the file they are trying to run - Once the data is loaded, running the function "main(kmax,filename)" (resp. "main(kmax,p,filename)") will output to "filename.tsv" the MAWs (resp. pMAWs) of length lesser than kmax.

-Warning- For program "pMAW_dict.py", the algorithm relies on MAWs, which means in order to run this program, the user will need to first generate the MAWs with the program of their choice among the first three and enter the output file as an argument to "read_files()"
