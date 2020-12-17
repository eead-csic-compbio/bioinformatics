# Session 1: Annotation of coding sequences

<!-- Sequence annotation (alignments, dna, proteins, domains, modelling) --> 

Analysis of sequences is one of the foundations of Genomics and Computational Biology. 
Genomic sequences are most often either nucleotide sequences and peptide sequences. 
The first can be large genome fragments such as contigs or smaller genes, transcripts, 
coding sequences or non-coding RNAs; peptides are usually proteins translated from open 
reading frames encoded in coding sequences.

In this session we will focus on coding sequences, or in other words, genes that encode proteins.

## Sequence comparison

A natural way of comparing protein or nucleic acid molecules is to align their (primary) sequences. 
This is due to the general expectation that sequence drives folding and to the fact that sequences are easy to work with in a computer or even your notebook. At least easier than structure.

![](./pics/align2.png)

## Pairwise alignment: edit distance

When two sequences are aligned residues from one sequence are matched one-by-one to residues in the other. Matches are obvious when the sequences are nearly identical, but less so when mutations accumulate. A simple way to compute how similar two sequences are is to compute their edit distance.

![](./pics/align2edit.png)

## Pairwise alignment: sequence identity

Another way to compute how similar two sequences are is to compute their % sequence identity.

![](./pics/seqidcalc.png)

## Substitution matrices

Sequence identity is a simple way of measuring conservation. However, it lacks resolution and handles all residue substitutions the same way. This is not ideal as we know that purine (A,G) and pyrimidine (C,T) nucleotides, or aromatic amino acid resides if we talk about proteins, are often not interchanged with equal probability in real genes or proteins:

![](./pics/align_substitution.png)

These preferences are captured by computing log-odds ratios of frequencies of observed mutations (a,b) with respect to estimates assuming no preference:

![equation](https://latex.codecogs.com/gif.latex?s(a,b)&space;=&space;\lambda&space;\&space;log(\frac{f_{ab}}{f_{a}&space;f_{b}})&space;\approx&space;log\frac{f_{homologues}}{f_{bychance}})

These log-odds are additive. 

## BLOSUM substitution matrices

The most frequent substitution matrices used to score protein alignments are the  [BLOSUM](https://en.wikipedia.org/wiki/BLOSUM) matrices. These matrices are described by a number X, as in BLOSUM50, derived from the analysis of alignments of protein blocks/domains with identities < X percent.  Below you can see BLOSUM50:

![](./pics/blosum50.png)

These log-odds have been scaled so that they can be accurately represented by integers.

BLOSUM matrices are scaled to 1/2-bit units [(Pearson2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3848038); a substitution score such as s(E,E) can be broken down to: 

![equation](https://latex.codecogs.com/gif.latex?s(E,E)&space;=&space;6&space;=&space;2.0&space;\&space;log_{2}(\frac{f_{EE}}{f_{E}&space;f_{E}}))

![equation](https://latex.codecogs.com/gif.latex?\frac{f_{EE}}{f_{E} f_{E}}&space;=&space;2^3)

S(E,E) is thus 8 times more likely to occur because of homology than by chance.

## Pairwise alignment: similarity

By using matrices such as BLOSUM it is possible to compute the similarity between two aligned sequences, which is added up along the alignment:

![](./pics/align2blosum50.png)

## Pairwise alignment: handling insertions and deletion (indels)

In addition to residue substitutions, insertions and deletions are usually considered while computing similarity. This can be done in many ways. The simplest is to assume a **linear cost** for insertions, proportional to their length. However, it is more accurate to compute **affine gap costs**, which charge a fix cost to openning a gap (a) and then a linear cost proportional to the gap length (bk). In this way, a gap of k residues receives a total score of -(a+bk)

![](./pics/align_affine.png)

## Multiple alignment

When more than two sequences are to be aligned we talk about multiple alignments, which can be computed in many ways. The most intuitive way, comparing them all, requires quadratically more resources as more sequences are added.

![](./pics/align3.png)


# Algorithms for sequence alignment

In this section we will visit some of the most frequent algorithms used to align sequences. The goal is to learn what you can and cannot do with each of them so that you can apply them correctly.

## Pairwise alignments

These are the simplest alignments as they involve only two sequences (of length *m* and *n*). There are several flavours, but the most important are global and local alignments, depicted in this figure taken from [(Makinen2015)](http://www.genome-scale.info):

![](./pics/global_local.jpg)

### Global alignment 

The Needleman-Wunsch (NW) algorithm is the original deterministic approach for 
the alignment of pairs of protein sequences end-to-end [(Needleman1970)](https://www.ncbi.nlm.nih.gov/pubmed/5420325). 
It was subsequently optimized by [(Gotoh1982)](https://doi.org/10.1016/0022-2836(82)90398-9). 
These algorithms belong to the family of [dynamic programming (DP) algorithms](https://en.wikipedia.org/wiki/Dynamic_programming).
The NW and related algorithms break the initial problem in subalignments which are recursively solved. In order to reach the final solution, an alignment matrix DP must be filled iteratively:

+ The top and left margins of DP are filled with gaps of increasing length, and the origin set to zero:

![equation](https://latex.codecogs.com/gif.latex?DP(0,0)&space;=&space;0)
![equation](https://latex.codecogs.com/gif.latex?DP(i,0)&space;=&space;-id&space;,&space;\text{&space;for&space;}&space;1\leq&space;i&space;\leq&space;m)
![equation](https://latex.codecogs.com/gif.latex?DP(0,j)&space;=&space;-jd&space;,&space;\text{&space;for&space;}&space;1\leq&space;j&space;\leq&space;n)

+ The alignment DP matrix can now be computed from top left to bottom right according to the following recursive function, where $s(x,y)$ is a scoring function/substitution matrix and $d$ a linear gap cost:

