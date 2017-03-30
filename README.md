# CESAR 2.0

CESAR 2.0 (Coding Exon Structure Aware Realigner 2.0) is a method to realign coding exons or genes to DNA sequences using a Hidden Markov Model [1].

Compared to its predecessor [2], CESAR 2.0 is 77X times faster on average (132X times faster for large exons) and requires 30-times less memory. In addition, CESAR 2.0 improves the accuracy of the comparative gene annotation by two new features. First, CESAR 2.0 substantially improves the identification of splice sites that have shifted over a larger distance, which improves the accuracy of detecting the correct exon boundaries. 
Second, CESAR 2.0 provides a new gene mode that re-aligns entire genes at once. This mode is able to recognize complete intron deletions and will annotate larger joined exons that arose by intron deletion events. 



# Installation
Just call 

`make` 

to build CESAR2.

The code is commented in doxygen style.
To compile a doxygen documentation of this program at `doc/doxygen/index.html`, call 

`make doc` 

# Running CESAR 2.0 directly
## Minimal example

Just call

`./cesar example/example1.fa`

This will output the re-aligned exon, using the default donor/acceptor profile obtained from human. 


## Format of the input file
The input file has to be a Fasta file. It provides at least one reference and
at least one query sequence. References and queries have to be separated by a
line starting with '#'. References are the exons (together with their reading frame) that you want to align to the query sequence.

Example alignment of human exon against a mouse query sequence.
```
>human
acACGTACGTgt
####
>mouse
ACGTACGTACGTACGTACGTACGTACGTACGT
```

Example alignment of multiple human exons against multiple mouse queries.
```
>human#0
acACACGTgt
>human#1
acACGTGTgt
>human#2
acACGTACGTgt
####
>mouse-1
ACGTACGTACGTACGTACGTACGTACGTACGT
>mouse-2
ACGTACGTACGTACGTCGTCGTCGTCGTAAAAACGTACGTACGTACGTACGT
```


## Parameters

`-f/--firstexon`
The default profile for start codons is assigned to the acceptor profile of
the first given exon.

`-l/--lastexon`
The default profile for stop codons is assigned to the donor profile of the last given exon.

`-m/--matrix <matrix file>`
Set `<matrix file>` as the path to the substitution matrix.

`-p/--profiles <acceptor> <donor>`
Set acceptor and donor profiles to `<acceptor>` resp. `<donor>`.

`-c/--clade <human|mouse>` (default: `human`)
A shortcut to default sets of substitution matrix and profiles.
For example, `-c human` is synonymous to:
`-m extras/tables/human/eth_matrix.txt -p extras/tables/human/acc_profile.txt extras/tables/human/do_profile.txt`

By default, CESAR2 uses profiles obtained from human.
You can provide profiles for another species in a directory extra/tables/$species and tell CESAR 2.0 to use these profiles by
`./cesar --clade $species test/mocks/example1.fa`

If <clade> contains a slash `/` it will be interpreted as look-up directory for profiles.

**Note:** With `-l` and/or `-f`, the profiles will change accordingly.



## Special parameters

`-v/--verbosity <n>`
Print extra information to stderr.

n  | Information
------------- | -------------
1  | Input Parameters
2  | List matrices and sequences in memory
3  | Fasta parser and alignment state machine
4  | Emission table initialization and Viterbi path
5  | HMM state creation, transitions and HMM normalization
6  | Full Viterbi step
7  | Initialization and access of emission tables


`-i/--split_codon_emissions <acc split codon emissions> <do split codon emissions>`
Manually define the length of split codons for each reference at once.

**Note:** `-i` is deprecated. Use lower case letters the Fasta file to annotate
split codons and upper case letters for all other codons. Alternatively
separate split codons from full codons with the pipe character `|`.


`-s/--set <name1=value1> .. <nameN=valueN>`
Customize parameters, e.g. transition probabilities. <!--- TODO: List available parameters. -->

Use with caution!


`-V/--version`
Print the version and exit.


# References
CESAR 2.0 was implemented by Peter Schwede (MPI-CBG/MPI-PKS 2017).

[1] Sharma V, Schwede P, and Hiller M. CESAR 2.0 substantially improves speed and accuracy of comparative gene annotation. Submitted

[2] Sharma V, Elghafari A, and Hiller M. [Coding Exon-Structure Aware Realigner (CESAR) utilizes genome alignments for accurate comparative gene annotation](https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/gkw210). Nucleic Acids Res., 44(11), e103, 2016

