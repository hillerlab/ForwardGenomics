# Forward Genomics

Forward Genomics is a framework to associate phenotypic differences between species to differences in their genomes [1]. 
This framework focuses on phenotypes that were present in the common ancestor of a group of species and repeatedly lost in independent descendant species. 
Forward Genomics relies on two important principles. 
First, loss of a phenotype should over time lead to loss of the genetic information necessary for this phenotype (also called the "use it or lose it" principle). 
The reason is that the genetic information for the given phenotype is no longer under purifying selection in those descendant species that lost the phenotype. 
Consequently, neutral evolution will eventually lead to the loss of this information in the trait-loss lineages. 
In contrast, purifying selection will typically preserve the genetic information for the given phenotype in descendant species that also preserve the phenotype. 
Second, repeated loss of this phenotype should lead to loss of this genetic information in the independent trait-loss lineages, irrespective of where the first mutation that caused the loss of this trait is located. 
Over time, this gives an evolutionary signature where genomic regions are conserved in the trait-preserving lineages and diverged or completely lost in the trait-loss lineages. 
Forward Genomics uses this independent divergence signature to obtain specificity in a genome-wide search for genomic regions involved in the given independent phenotype loss. 

The R script provides three different Forward Genomics methods:
* "perfect match" searches for genomic regions where *all* trait-loss species are more diverged than *all* trait-preserving species [1]. 
* "GLS" uses a generalized least square approach [2].
* "branch method" uses per-branch divergence values [2].

In constrast to perfect-match, GLS and branch method directly control for the phylogenetic relatedness and evolutionary rate differences between species. 
On simulated data, GLS and especially the branch method, have substantially more power than perfect match [2]. 


# Installation

# Input Files
You need to have these input files:

1: The tree with branch lengths in newick format. The ancestors must be named, otherwise use tree_doctor -a to name them. 

Example: 
```
((Human:0.148845,(Rat:0.091589,Mouse:0.084509)Rat-Mouse:0.271974)Human-Rat:0,Dog:0.194257)Human-Dog
```


2: A file listing the identifiers of each element (genomic region) that should be processed. If you want to process all elements, you can get this list from the global percent-identity file with 'tail -n +2 globalPid.file | cut -f1 -d " "'

Example: 
```
ID1
ID5
ID8
```

3: A file listing which species has the phenotype (state 1) and which species have lost it (state 0). List all species (leaves in the tree) where you know the phenotypic state, omit all others. Must be a space-separated file the header 'species pheno'.

Example: 
```
species pheno
Human 1
Mouse 0
```


4: A file listing the global %id values if you want to run the GLS method and/or a file listing the local (per-branch) %id values if you want to run the branch method. 

Format globalPid.file: Space-separated table with global %id values per element (row) x species (columns). First column must be the element identifier. The first line must start with 'species'. 
Example:
```
species Dog Human Mouse Rat
ID1 0.96875 0.9375 0.91875 0.9125
ID2 0.97817 0.97101 0.87234 0.8913043
```

Format localPid.file: Space-separated file listing the local %id values for one element and one branch on a single line. Must have the header 'branch id pid'. 

Example:
```
branch id pid
Human ID1 0.99375
Rat ID1 0.97243
Mouse ID1 0.94129
Rat-Mouse ID1 0.90097
```


# HowTo


# References
[1] Hiller M, Schaar BT, Indjeian VB, Kingsley DM, Hagey LR, and Bejerano G. (2012): A "forward genomics" approach links genotype to phenotype using independent phenotypic losses among related species. Cell Reports, 2(4), 817-823

[2] Prudent X et al. in preparation
