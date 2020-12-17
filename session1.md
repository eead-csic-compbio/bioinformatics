# Session 1: Annotation of coding sequences

<!-- Sequence annotation (alignments, dna, proteins, domains, modelling) --> 

Analysis of sequences is one of the foundations of Genomics and Computational Biology. 
Genomic sequences are most often either nucleotide sequences and peptide sequences. 
The first can be large genome fragments such as contigs or smaller genes, transcripts, 
coding sequences or non-coding RNAs; peptides are usually proteins translated from open 
reading frames encoded in coding sequences.

In this session we will focus on coding sequences, or in other words, genes that encode proteins.

## Sequence comparison

A natural way of comparing protein or nucleic acid molecules is to align their (primary) sequences. This is due to the general expectation that sequence drives folding and to the fact that sequences are easy to work with in a computer or even your notebook. At least easier than structure.

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

$$ s(a,b) = \lambda \ log(\frac{f_{ab}}{f_{a} f_{b}}) \approx  log\frac{f_{homologues}}{f_{bychance}} $$

These log-odds are additive. 





