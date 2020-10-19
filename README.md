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
You need R version 3.02 or later and the R packages [caper](https://cran.r-project.org/web/packages/caper), [xtermStyle](https://cran.r-project.org/web/packages/xtermStyle), [phangorn](https://cran.r-project.org/web/packages/phangorn), [weights](https://cran.r-project.org/web/packages/weights), and [isotone](https://cran.r-project.org/web/packages/isotone). forwardGenomics.R will try to install them. 

You also need tree_doctor from the [phast](http://compgen.cshl.edu/phast/) package. tree_doctor must be in your $PATH. Test it by running ```tree_doctor -p``` in your command line.

# Input Files
You need to have these input files:

#####1: Phylogenetic tree
This file should contain the tree with branch lengths in newick format. 
The ancestors must be named, otherwise use ```tree_doctor -a``` to name them. 

Example: 
```
((Human:0.148845,(Rat:0.091589,Mouse:0.084509)Rat-Mouse:0.271974)Human-Rat:0,Dog:0.194257)Human-Dog
```


#####2: List of element identifiers
This file should list the identifiers of each element (genomic region) that should be processed. 
If you want to process all elements, you can get this list from the global percent-identity file with 
```tail -n +2 globalPid.file | cut -f1 -d " "```

Example: 
```
ID1
ID5
ID8
```

#####3: Phenotype list
A file listing which species has the phenotype (state 1) and which species have lost it (state 0). List all species (leaves in the tree) where you know the phenotypic state, omit all others. Must be a space-separated file the header 'species pheno'.

Example: 
```
species pheno
Human 1
Mouse 0
```


#####4: %id values
A file listing the global %id values if you want to run the GLS method and/or a file listing the local (per-branch) %id values if you want to run the branch method. 

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
By default, forwardGenomics.R runs perfect-match, GLS and the branch method. This means, you need to provide both local and global %id values. Use the --methods parameter to run only perfect-match/GLS on global and only the branch method on local %id values. 

The branch method should be run separately on genomic regions that are coding (CDS) and non-coding (CNE). The reason is that the branch method needs to know the branch weights and expected %id values for every branch, and they differ for coding and non-coding regions. The directory lookUpData provides weights and expected %id values for both coding and non-coding regions. By default, the branch methods assumes you have coding regions (forwardGenomics.R will load branchWeights_CDS.txt and expPercentID_CDS.txt in lookUpData/). **If you have non-coding regions, you need to set --weights and --expectedPerIDs parameter.**

By default, GLS will normalize the %id values to control for differences in evolutionary rates and it will only consider elements where the given %id values imply at least 2 independent loss events to avoid elements with lineage-specific losses (e.g. in case 1 of 2 independent trait-loss species has missing data). 

All default values can be changed by parameters. 

## Example usage
```
forwardGenomics.R --tree=example/tree_ancestor.nh --elementIDs=example/IDlist.txt --listPheno=example/listPhenotype.txt --globalPid=example/globalPercentID_CDS.txt --localPid=example/localPercentID_CDS.txt --outFile=example/myOutput.txt
```
## All parameters
```
Mandatory arguments:
--tree          = filename
                Phylogenetic tree with branch lengths in newick format. The species names must be identical to the names in the percent ID files. Ancestors must be named (otherwise use tree_doctor -a)

--elementIDs    = filename
                File listing the ID of each element (genomic region) that should be processed. 
                If you want to process all elements, you can get this list with 'tail -n +2 globalPercentID.file | cut -f1 -d " "'

--listPheno     = filename
                File listing the phenotype for all species (all leaves in the tree). Must be a space-separated file the header 'species pheno', where 0 means trait is lost and 1 means trait is present.

--globalPid     = filename
                Only if you run GLS: Space-separated input file with the global %id values arranged in elements (row) x species (columns). First column must be the element ID. The first line must start with 'species'. 

--localPid      = filename
                Only if you run the branch method: Space-separated input file with the local %id values. Must have the header 'branch id pid'. One line per element and per branch. 

--outFile       = filename
                Output file that will contain the element ID and the P-values from the methods.


General optional arguments:
--method        = 'branch', 'GLS', 'all'
                Which method to run. GLS includes computing the margin of the perfect-match method. Default is all.

--verbose       = TRUE / FALSE
                Show much more info and create plots for each element. Default is FALSE

--outPath       = /dir/to/outputScatterPlots/
                If verbose==TRUE, this directory will be created and will contain scatter plots for each element. If verbose==FALSE, this parameter has no effect. Default directory '.'


Optional arguments for GLS:
--minLosses     = number
                In case of missing data (no %id value for some species) only consider elements where at least this many independent loss events are supported with %id values. Can be used to exclude lineage-specific losses. Default 2

--transf        = 'raw' or 'normalized'
                Whether to use raw global %id values or normalize them for the differences in evolutionary rates. Default normalized


Optional arguments for the branch method:
--weights       = filename
                File listing the weights per branch. Must have a header 'len w'. Default is the file for coding exons: lookUpData/branchWeights_CDS.txt. Use lookUpData/branchWeights_CNE.txt for non-coding genomic regions.

--expectedPerIDs = filename
                File listing the expected %id values for each branch length. Must have a header 'len mPid'. Default is the file for coding exons: lookUpData/expPercentID_CDS.txt. Use lookUpData/expPercentID_CNE.txt for non-coding genomic regions.

--thresholdConserved = floating point number
                Some element may have large indels on internal conserved branches, but descendant branches are highly conserved. We reject genomic elements if a local %id value is lower than this threshold for a conserved branch. Set to 0 to ignore this. Default: 0.5
```

# Computing %id values
This is a minimal recipe.

Setup and input: You need the binaries of the UCSC genome browser (kent) source code, in particular twoBitToFa, and prank. The path to these binaries must be included in the $PATH variable. 
Each assembly must be in a $genomePath/gbdb-HL/$assembly/$assembly.2bit directory structure. The $genomePath variable can be set as a bash environment variable and is read out by the scripts. 
The region of interest should be contained in maf file. 

To minimize the number of input and output files in big cluster jobs, we typically store and retrieve the results in a BDB file hash, which requires a installing BDB. Without BDB, a few changes in the scripts will be necessary to output into a text file, instead of a BDB file.  

```Maf2SpanningSeq_PRANK.perl maf.file ElementID -runPrank -treeFile tree.file -BDBFile output.bdb```
where maf.file is the input maf, ElementID is the name of the element, tree.file contains the phylo tree in Newick format and output.bdb determines the output BDB File that will contain the alignment with reconstructed ancestors. 

The script extracts the start and end coords of every species (the reason is that a single region is often broken into different maf blocks and there can be insertions between these maf blocks). Then it runs UCSC twoBitToFa to extract the seq of all species and runs Prank with the given species tree, which gives a full alignment including reconstructed ancestors. 


```GetGlobalAndLocalPercentID.perl output.bdb ElementID -treeFile tree.file -allowedAncestralNodes ancNodes -local -global -GlobalBDBFile global.bdb  -LocalBDBFile local.bdb```
where ancNodes is comma-separated list of node names that refer to ingroup ancestors to be used for %id calculation, and global/local.bdb refers to the BDB files containing the %id values. They can be later read out by ReadBDB.perl



# References
[1] Hiller M, Schaar BT, Indjeian VB, Kingsley DM, Hagey LR, and Bejerano G. [A "forward genomics" approach links genotype to phenotype using independent phenotypic losses among related species](http://www.cell.com/cell-reports/fulltext/S2211-1247(12)00272-0). Cell Reports, 2(4), 817-823, 2012

[2] Prudent X, Parra G, Schwede P, Roscito JG, and Hiller M. [Controlling for phylogenetic relatedness and evolutionary rates improves the discovery of associations between species’ phenotypic and genomic differences](https://academic.oup.com/mbe/article-lookup/doi/10.1093/molbev/msw098). Mol Bio Evol. 33(8):2135-50, 2016
