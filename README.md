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


# HowTo


# References
[1] Hiller M, Schaar BT, Indjeian VB, Kingsley DM, Hagey LR, and Bejerano G. (2012): A "forward genomics" approach links genotype to phenotype using independent phenotypic losses among related species. Cell Reports, 2(4), 817-823
[2] Prudent X et al. in preparation
