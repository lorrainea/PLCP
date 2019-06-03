Longest Common Prefixes with k-errors
===

<b>Description</b>: Given an input sequence x in FATSA format and an integer value k, PLCP is a program that computes for each index i of x, the longest prefix of x[i] that exists at position j in x with k mismatches.

<b>Installation</b>: To compile plcp, please follow the instructions given in file INSTALL.

<b>INPUT</b>: A single sequence in FASTA format.

<b>OUTPUT</b>: Array PLCP[i] which holds the the longest prefix of x[i] that exists at position j in x with k mismatches and array P[i]=j.

```
 plcp <options>
 Mandatory arguments:
  -i, --input-file	<str>	Sequence input filename (FASTA format).
  -o, --output-file	<str>	Output filename.
  -k, --hamming-dist	<int>	Hamming distance between matches.
 Optional arguments:
  -t, --threads		<int>	The number of threads to be used (default: 1).
```

<b>Citation</b>
```
L. A. K. Ayad, C. Barton, P. Charalampopoulos, C. S. Iliopoulos, S. P. Pissis, "Longest Common Prefixes with k-Errors and Applications", SPIRE, Springer, 2018, pp. 27-41`
```

GNU GPLv3 License; Copyright (C) 2018 Lorraine A.K. Ayad, Carl Barton, Panagiotis 
Charalampopoulos, Costas S. Iliopoulos and Solon P. Pissis

