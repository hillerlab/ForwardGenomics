# Xavier Prudent, 2015

######################################################
## Load the libraries
######################################################
pkg=c("xtermStyle","phangorn","weights","isotone", "caper")
for( p in pkg){
  if( p %in% rownames(installed.packages()))
    {
      cat(paste("Package", p, "found\n"))
      require( p, character.only=TRUE )
    }else{
      cat(paste("Package", p, "not found\n"))
      install.packages( p, repos="http://mirrors.softliste.de/cran/" )
      require( p, character.only=TRUE )
    }
}

args = (commandArgs(TRUE))

######################################################
## Read the input arguments
######################################################
readArguments = function(){
  ## Default setting when no arguments passed
  if(length(args) < 1) {
    args <- c("--help")
  }
  
  ## Help section
  if("--help" %in% args) {
        cat(style( "\nForward Genomics:",fg="blue",font.style="bold" ) )
        cat("
 
      Mandatory arguments:
      --tree          = filename
                      Phylogenetic tree with branch lengths in newick format. The species names must be identical to the names in the percent ID files. Ancestors must be named (otherwise use tree_doctor -a)

      --elementIDs    = filename
                      File listing the ID of each element (genomic region) that should be processed. 
                      If you want to process all elements, you can get this list with 'tail -n +2 globalPercentID.file | cut -f1 -d \" \"'

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

      Example:
      forwardGenomics.R --tree=example/tree_ancestor.nh --elementIDs=example/IDlist.txt --listPheno=example/listPhenotype.txt --globalPid=example/globalPercentID_CDS.txt --localPid=example/localPercentID_CDS.txt --outFile=example/myOutput.txt\n\n")
    q(save="no")

##experimental
##      --collapseClades = 'no', 'all', 'neutralOnly',  or 'Mouse-Squirrel' or 'Mouse-Squirrel,Human-Chimp,Cow-Dog'
##                    String, should the conserved branches be collapsed? If 'all', all clades are collapsed, if only some ancestors are specificied (separated by commas without any whitespace), the corresponding clades are collapsed. If 'neutralOnly', the conserved branches are not collapsed, but the branches under neutral evolution are.

  }
  
  ## Parse arguments (we expect the form --arg=value)
  parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
  argsDF <<- as.data.frame(do.call("rbind", parseArgs(args)))

  # check for invalid arguments
  allValidParameters = c("tree","elementIDs","listPheno","globalPid","localPid","outFile","method","verbose","outPath","minLosses","transf","weights","expectedPerIDs","thresholdConserved")
  if (length(setdiff(argsDF$V1,allValidParameters)) > 0)
     stop( paste("ERROR:", setdiff(argsDF$V1,allValidParameters), " is an invalid argument.\n" ) )
  
  ## List of arguments to provide
  mandatoryArg = c( "elementIDs", "tree", "outFile", "listPheno" )
 
  ## General arguments
  for( a in mandatoryArg ){
    if( !length(which( argsDF$V1 == a )) ) stop(paste("Argument", a, "missing" ))
    
     ## Retrieve arguments
    in_elements  <<- as.character(argsDF[which( argsDF$V1 == "elementIDs" ),2])
    in_listPheno <<- as.character(argsDF[which( argsDF$V1 == "listPheno" ),2])
    in_tree      <<- as.character(argsDF[which( argsDF$V1 == "tree" ),2])
    out_data     <<- as.character(argsDF[which( argsDF$V1 == "outFile" ),2])

    ## Checks the existence of the files and directories
    for( thisFile in c( in_elements, in_listPheno, in_tree ) )
      if( ! file.exists( thisFile ) )
        stop( paste("ERROR:", thisFile, "is missing" ) )
    
  }

  ################
  ## Parse all optional arguments
  if( length(which( argsDF$V1 == "method" )) ) 
    in_analysis  <<- as.character(argsDF[which( argsDF$V1 == "method" ),2])
  else
    in_analysis  <<- "all"
  ## Checks the validity of the given inputs
  if( in_analysis != "branch" && in_analysis != "GLS" && in_analysis != "all" )
    stop("ERROR: parameter method should be either branch, GLS or all")    

  # transformation
  if( length(which( argsDF$V1 == "transf" )) ) 
    in_rawRanked  <<- as.character(argsDF[which( argsDF$V1 == "transf" ),2])
  else
    in_rawRanked  <<- "normalized"
  ## Checks the validity of the given inputs
  if( in_rawRanked != "raw" && in_rawRanked != "normalized" )
    stop( "ERROR: parameter transf should be either raw or normalized" )

  # minimum number of independent losses
  if( length(which( argsDF$V1 == "minLosses" )) ) 
    indepLoss    <<- as.numeric(as.character(argsDF[which( argsDF$V1 == "minLosses" ),2]))
  else
    indepLoss  <<- 2

  # verbose
  if( length(which( argsDF$V1 == "verbose" )) ) 
    verbose <<- as.logical(argsDF[which( argsDF$V1 == "verbose" ),2])
  else
    verbose <<- FALSE
  if( verbose != TRUE && verbose != FALSE )
     stop("ERROR: verbose should be either TRUE or FALSE")

  # conserved branch threshold
  if( length(which( argsDF$V1 == "thresholdConserved" )) ) 
    in_thresholdConserved <<- as.numeric(as.character(argsDF[which( argsDF$V1 == "thresholdConserved" ),2]))
  else
    in_thresholdConserved <<- 0.5
  if( in_thresholdConserved < 0 || in_thresholdConserved > 1 )
      stop( "ERROR: thresholdConserved should be between 0 and 1" )    

  # weights and expected %id values
  if( length(which( argsDF$V1 == "weights" )) ) 
    in_brWeights      <<- as.character(argsDF[which( argsDF$V1 == "weights" ),2])
  else
    in_brWeights      <<- paste(scriptDir,"/lookUpData/branchWeights_CDS.txt", sep="")

  if( length(which( argsDF$V1 == "expectedPerIDs" )) ) 
    in_meanPidSel     <<- as.character(argsDF[which( argsDF$V1 == "expectedPerIDs" ),2])
  else
    in_meanPidSel     <<- paste(scriptDir,"/lookUpData/expPercentID_CDS.txt", sep="")

  # output dir
  if( length(which( argsDF$V1 == "outPath" )) ) 
    out_path     <<- as.character(argsDF[which( argsDF$V1 == "outPath" ),2])
  else
    out_path     <<- "."
  ## Create the dir if missing, if not does not do anything
  dir.create( path=out_path, showWarnings=FALSE) 


  #############
  ## Display general parameters
  cat("\n")
  cat(style("Forward Genomics analysis:\n",fg="red",font.style=c("bold") ))
  cat("\n")
  cat(style(paste("Method(s).................................", in_analysis, "\n"),fg="red",font.style="bold" ))
  cat(style(paste("List of elements..........................", in_elements, "\n"),fg="red",font.style="bold" ))
  cat(style(paste("Phenotype of all species..................", in_listPheno, "\n"),fg="red",font.style="bold" ))
  cat(style(paste("Phylogenetic tree.........................", in_tree, "\n"),fg="red",font.style="bold" ))
  cat(style(paste("Output file...............................", out_data, "\n"),fg="red",font.style="bold" ))
  if (verbose == TRUE)
     cat(style(paste("Output dir................................", out_path, "\n"),fg="red",font.style="bold" ))
  cat(style(paste("Verbose...................................", verbose, "\n"),fg="red",font.style="bold" ))

  #############
  # Parse parameters for the branch method
  if( in_analysis == "branch" || in_analysis == "all" ){
    if( !length(which( argsDF$V1 == "localPid" )) ) 
       stop(paste("Argument localPid is missing if you want to run the branch method" ))
    else
       in_localPid    <<- as.character(argsDF[which( argsDF$V1 == "localPid" ),2])

    ## Checks the existence of the files and directories
    for( thisFile in c( in_localPid, in_brWeights, in_meanPidSel ) )
      if( ! file.exists( thisFile ) ) stop( paste("ERROR:", thisFile, "is missing" ) )

    ## Display
    cat(style("\nParameters for the branch method:\n",fg="red",font.style=c("bold") ))
    cat(style(paste("Local %id values..........................", in_localPid, "\n"),fg="red",font.style="bold" )) 
    cat(style(paste("Branch weights............................", in_brWeights, "\n"),fg="red",font.style="bold" ))
    cat(style(paste("Expected %id values.......................", in_meanPidSel, "\n"),fg="red",font.style="bold" ))
    cat(style(paste("%id threshold trait-preserving branches...", in_thresholdConserved, "\n"),fg="red",font.style="bold" ))

## experimental
##    in_collapseClades <<- as.character(argsDF[which( argsDF$V1 == "collapseClades" ),2] )
##    if( in_collapseClades != "neutralOnly" && in_collapseClades != "all" && in_collapseClades != "no" )
##      stop( "ERROR: collapseClades should be either neutralOnly, all, or no" )
##    cat(style(paste("Collapsing clade..........................", in_collapseClades, "\n"),fg="red",font.style="bold" ))
  }
  

  #############
  # Parse parameters for GLS
  if( in_analysis == "GLS" || in_analysis == "all" ){
    if( !length(which( argsDF$V1 == "globalPid" )) ) 
       stop(paste("Argument globalPid is missing if you want to run the GLS method" ))
    else
       in_globalPid   <<- as.character(argsDF[which( argsDF$V1 == "globalPid" ),2])

    ## Checks the existence of the files and directories
    for( thisFile in c( in_globalPid ) )
      if( ! file.exists( thisFile ) ) stop( paste("ERROR:", thisFile, "is missing" ) )

    ## Display
    cat(style("\nParameters for GLS:\n",fg="red",font.style=c("bold") ))
    cat(style(paste("Global %id values.........................", in_globalPid, "\n"),fg="red",font.style="bold" ))
    cat(style(paste("Min number of independent losses..........", indepLoss, "\n"),fg="red",font.style="bold" ))
    cat(style(paste("Raw/normalized............................", in_rawRanked, " (perfect match is computed with raw values only) \n"),fg="red",font.style="bold" ))
 }
 

  # switch off experimental feature
  in_collapseClades <<- "no"
}




########################################################
## Open input files
########################################################
inputFiles = function()
{
  cat(style("\nReading input data...\n",fg="red",font.style=c("bold") ))

  ## ---------------
  ## List of elements
  in_elements <<- read.table( in_elements )
  if( nrow( in_elements ) == 0 || ncol( in_elements ) != 1 )
    stop( "ERROR: the list of elements is not valid" )
  
  ## ---------------
  ## List of phenotypes
  
  ## Open the list of phenotypes
  listPheno <<- read.table( in_listPheno, header = TRUE )  
  ## Check the header of the list of phenotypes
  if( "species" %in% colnames(listPheno) == FALSE || "pheno" %in% colnames(listPheno) == FALSE )
    stop( "ERROR: header of the list of phenotypes must be 'species pheno'" )
  
  ## ---------------
  ## Tree
  
  ## Open the phylogeny and add nodes names using tree_doctor
  cmd = paste( "tree_doctor -n -a", in_tree )
  tree.withNodes = system( cmd, intern = TRUE ) 
  tree <<- read.tree( text = tree.withNodes )
  
  ## Branch lengths in the tree
  totLength = diag( vcv.phylo( tree ) )
  if( verbose ){
    cat(style("\nTotal branch lengths:\n",fg="red",font.style=c("bold") ))
    print(as.data.frame(totLength))
  }
  ## Propagate the phenotypes into the tree, spot the loss nodes
  in_lossNodes <<- propagatePheno( tree )

  ## Fill the branch catalog
  treePreparationOutput = prepareTree( tree, in_lossNodes )  
  vcvClade <<- treePreparationOutput[[1]]
  branchCatalog <<- treePreparationOutput[[2]]
  vcvCladeSel <<- treePreparationOutput[[3]]
  
  ## Check that each tree-leaf has a phenotype
  nbLeaf = length( tree$tip.label )
  nbLeafWithPheno = length( which( tree$tip.label %in% listPheno$species == TRUE ) )
  if( nbLeafWithPheno != nbLeaf )
    stop( "ERROR: header of the list of phenotypes must be 'species pheno'")
  
  ## Check phenotypes are 0 or 1
  nLoss <<- length(which( listPheno$pheno == 0 ))
  nCons <<- length(which( listPheno$pheno == 1 ))
  if( (nLoss + nCons) != length(listPheno$pheno) )
    stop( "ERROR: phenotypes should be 0 (loss) or 1 (conserved)'" )
  
  ## ---------------
  ## Input local percent id
  if( in_analysis == "branch" || in_analysis == "all" )
    {
      ## Open the data file
      in_localPid <<- read.table( in_localPid, header=T )

      ## Extract the data for the elements under study
     in_localPid <<- in_localPid[ in_localPid$id %in% in_elements$V1 , ]
      
      ## Check for double entries
      if( nrow(in_localPid) != nrow( in_localPid[ !duplicated(in_localPid), ]) )
        stop( "ERROR: duplicated lines were found in the input data file" )
      
      ## Check for the data type
      if( typeof(in_localPid$pid) != "double" )
        stop( "ERROR: you may have non numerical values in the pid row" )
      
      ## Check for the header
      desiredHeader = c( "branch", "id", "pid" )
      if( length( which( desiredHeader %in% colnames(in_localPid) == TRUE ) ) != length( desiredHeader ) )
        stop( "ERROR: Missing or wrong header for the data file" )

      ## Replace the name of the row 'branch'
      if( "branch"  %in% colnames(in_localPid) ) colnames(in_localPid)[ which(colnames(in_localPid) == "branch") ] <<- "br"

      ## Remove useless branch length row
      if( "len" %in% colnames(in_localPid) ) in_localPid <<- in_localPid[ -which(colnames(in_localPid) == "len") ]
    }
  
  ## ---------------
  ## Input global percentID
  if( in_analysis == "GLS" || in_analysis == "all" )
    {
      ## Open the data file
      in_globalPid <<- read.table( in_globalPid, row.names=1, header=T )

      ## Remove from the data the species that don't appear in the tree
      rm.Species()

      ## Phenotype
      tips.phenotype()

      ## Merge the raw global percentID and the phenotype dataframes
      in_rawGlobalPid_pheno <<- rbind( in_globalPid, myPheno )

      ## Ranked percentId = ( percentId - 1 )/ total branch length
      if( in_rawRanked == "normalized" ){    
        for( i in 1:length(totLength) ){
          taxa    = names(totLength)[i]
          totL    = as.numeric(totLength[i])
          taxaCol = which( colnames(in_globalPid) == taxa )
          pid     = as.numeric(as.character(in_globalPid[,taxaCol]),stringsAsFactors = FALSE)
          in_globalPid[,taxaCol] =  ( pid - 1 ) / totL
        }
      }
      ## Merge the ranked global percentID and the phenotype dataframes
      in_globalPid_pheno <<- rbind( in_globalPid, myPheno )
      
      ## Species with no phenotypes (NA)
      noPhenotype = colnames( in_globalPid_pheno )[ which( is.na( in_globalPid_pheno[ "pheno", ] ) ) ]
      cat( paste( "Removing", length(noPhenotype), "species with unknown phenotype\n" ) )
      
      ## Remove the tips with unknown phenotypes
      tree <<- drop.tip( tree, noPhenotype )      
    }
}


########################################################
## Check which species have lost their phenotypes
########################################################
tips.phenotype = function(){

  cat(style("\nPhenotype:\n",fg="red",font.style=c("bold") ))
  
  ## List of phenotypes
  myPheno <<- c()

  ## Loop over the species
  for( spc in colnames(in_globalPid) ){
    phenoValue = listPheno$pheno[ listPheno$species == spc ]
    cat( paste( spc, phenoValue, "\n", sep="\t" ) )
    myPheno <<- c( myPheno, as.numeric(phenoValue) )
  }
  myPheno = t(myPheno)
  myPheno <<- data.frame( myPheno )
  colnames(myPheno) <<- colnames(in_globalPid)
  rownames(myPheno) <<- c("pheno")
}


#############################################################
## Remove species from the data if they are not in the tree
#############################################################
rm.Species = function(){

  ## Loop over the species
  for( s in colnames(in_globalPid) ){

    ## If that species is not in the tree
    if( ! ( s %in% tree$tip.label ) ){
      in_globalPid <<- in_globalPid[ , -which( colnames(in_globalPid) == s ) ]
      cat( paste( "Species", s, " not found in the tree -> removed from data\n" ) )
    }
  }
}


########################################################
## Prepare output files
########################################################
outputFiles = function(){

  ## Header for the output per element
  nameOutput <<- out_data
  unlink( nameOutput )
  
  ## Basic infor
  if( in_analysis == "GLS" )       basicInfo = c( "elementID", "numTraitLossSpecies", "numTraitPreservingSpecies" )
  if( in_analysis == "branch"  )   basicInfo = c( "elementID", "numTraitLossBranches", "numTraitPreservingBranches" )
  if( in_analysis == "all" )       basicInfo = c( "elementID", "numTraitLossSpecies", "numTraitPreservingSpecies", "numTraitLossBranches", "numTraitPreservingBranches" )

  ## Local methods
  Res = c( "wPearson", "wPearsonPval" )
  BpB.Res = c( "weightedPearsonCorrelationCoeff", "weightedPearsonCorrelation_Pvalue" )
  if( in_collapseClades != "no" ) localRes = c( Res, BpB.Res ) else localRes = c( BpB.Res )
  ## Global method 
  globalRes = c( "PerfectMatchMargin", "GLS_Pvalue" )
  
  ## Make the header
  if( in_analysis == "GLS" ) header = c( basicInfo, globalRes, "\n" )
  if( in_analysis == "branch"  ) header = c( basicInfo, localRes, "\n" )
  if( in_analysis == "all" ) header = c( basicInfo, globalRes, localRes, "\n" )  
  cat(header, file=nameOutput, sep=" ")
} 

###################################################################
## Create the functions, interpolations of look up tables
###################################################################
prepareFn = function(){
  ## Function: Mean PID vs. branch length
  fileSimSelection = read.table( in_meanPidSel, header=T )
  fnSimSelection   <<- approxfun( fileSimSelection$len, fileSimSelection$mPid )

  ## Function: Neutral weight vs. branch length
  fileBrWeights = read.table( in_brWeights, header=T )
  fnBrWeights   <<- approxfun( fileBrWeights$len, fileBrWeights$w )
}


###################################################################
## From the leaves, propagate the phenotypes into the tree
###################################################################
propagatePheno = function( thisTree ){

  cat(style("\nPropagating the phenotypes from the leaves to the branches using Dollo parsimony...\n",fg="red",font.style=c("bold") ))

## Catalog of the status of the tips and nodes
  tips.nodes = c( thisTree$tip.label, thisTree$node.label )
  status = rep( x=-1, times=length(tips.nodes) )
  catalog = data.frame( tips.nodes, status )
  colnames( catalog )[1] = "br"
  colnames( catalog )[2] = "status"

  ## Reorder the tree so that nodes and tips follow the postorder
  thisTree = reorder.phylo( x=thisTree, order="postorder" )
  postOrderList = tips.nodes[ thisTree$edge[,2] ]  
  
  ## Loop over the tips and nodes in postorder
  for( iname in postOrderList ){

    ## Index in the catalog
    i = which( catalog$br == iname )

    ## Jump if status already set
    if( catalog$status[i] != -1 ) next
    
    ## It is a tip
    if( iname %in% thisTree$tip.label == TRUE ){

      ## Retrieve the status from the pheno file
      j = which( listPheno$species == iname )
      catalog$status[i] = listPheno$pheno[j]
    }
    
    ## It is a node
    if( iname %in% thisTree$node.label == TRUE ){

      ## Get its descendants
      j = which( tips.nodes == iname )
      desc = tips.nodes[Descendants( thisTree, j, type="children" )]
      if( length( desc ) != 2 ){
        cat("\nERROR: a node must have 2 descendants\n")
        break
      }

      ## Check the status of the two descendants
      k1 = which( catalog$br == desc[1] )
      k2 = which( catalog$br == desc[2] )
      stat1 = catalog$status[k1]
      stat2 = catalog$status[k2]
      if( stat1 == -1 || stat2 == -1 ) next

      ## Assign the status of the node according to the descendants
      ## 1,1 -> 1
      ## 0,0 -> 0
      ## 1,0 -> 1, 0,1 -> 1
      if( stat1 == 0 && stat2 == 0 )   catalog$status[i] = 0
      if( stat1 == 1 && stat2 == 1 )   catalog$status[i] = 1
      if( stat1 == 0 && stat2 == 1 )   catalog$status[i] = 1
      if( stat1 == 1 && stat2 == 0 )   catalog$status[i] = 1
    }
  }

## Identify the loss nodes and tips
vec.lossNodes = c()
for( i in 1:nrow(catalog) ){
  
  ## Name of the tip or node
  iname = as.character(catalog$br[i])
  
  ## Status
  stat = catalog$status[i]  
  if( stat == 1 ) next
  
  ## Parent node
  j = which( tips.nodes == iname )
  par = tips.nodes[ Ancestors( thisTree, j, "parent" ) ]
  
  ## Status of the parent
  k = which( catalog$br == par )
  if( length(k) == 0 ) next
  parStat = catalog$status[k]
  
  ## Conclude whether that node is a loss one
  if( parStat == 1 ) vec.lossNodes = c( vec.lossNodes, iname )
  }  

  ## Make a data frame from the loss nodes
  vec.lossTime = rep( times=length( vec.lossNodes ), -1 )
  lossNodes = data.frame( vec.lossNodes, vec.lossTime )
  colnames( lossNodes ) = c( "lossNode", "lossTime" )
  
  ## Display the loss node
  cat( "\nBranch status\n" )
  print(catalog)
  cat( "\nLoss nodes\n" )
  print(lossNodes)
  return(lossNodes)
}


###################################################################
## List the branches with their length, weight, mean percentID
###################################################################
prepareTree = function( thisTree, theseLossNodes ){

  cat("\n")
  cat(style("Preparing the tree...\n",fg="red",font.style=c("bold") ))

  ## List of tips and nodes
  tips.nodes = c( thisTree$tip.label, thisTree$node.label )
  br  = tips.nodes[thisTree$edge[,2]]

  thisvcvClade  = list()
  thisvcvCladeSel = list()

  ## -------------------------------------------------
  ## Global analysis, we only need to know the subtree
  if( in_analysis == "GLS" || in_analysis == "all" )
    {
      subTree = rep( -1, length(br) )
      thisBranchCatalog = data.frame( br, subTree )
      colnames(thisBranchCatalog)[1] = "br"
      colnames(thisBranchCatalog)[2] = "subTree"
      
      ## Loop over the loss clades
      for( i in 1:nrow( theseLossNodes ) )
        {
          n = as.character(theseLossNodes$lossNode[i])
          
          ## If the loss point is a node
          if( n %in% thisTree$node.label )
            {
              ## Extract the loss clade
              clade = extract.clade( phy=thisTree, node=n )                        
              ## Assign the subTree tag to each tip
              for( tip in clade$tip.label ) thisBranchCatalog$subTree[which( thisBranchCatalog$br == tip )] = i
            }
           ## Subtree if the loss point is a tip
          else thisBranchCatalog$subTree[which( thisBranchCatalog$br == n )] = i      
        }    
      ## Remove the useless internal nodes
      thisBranchCatalog = thisBranchCatalog[ -which( thisBranchCatalog$br %in% thisTree$node.label ), ]
    }

  ## -----------------------------------------------------------------
  ## Local analysis, a detailled description of each branch is required
  if( in_analysis == "branch" || in_analysis == "all" )
    {
      stat    = rep( -1, length(br) )
      w       = rep( -1, length(br) )
      mPid    = rep( -1, length(br) )
      thresh  = rep( in_thresholdConserved, length(br) )
      subTree = rep( -1, length(br) )
      path    = rep( -1, length(br) )
      thisBranchCatalog = data.frame( br, thisTree$edge.length, stat, w, mPid, thresh, subTree, path )
      colnames(thisBranchCatalog)[1] = "br"
      colnames(thisBranchCatalog)[2] = "len"
      colnames(thisBranchCatalog)[3] = "stat"
      colnames(thisBranchCatalog)[4] = "w"
      colnames(thisBranchCatalog)[5] = "mPid"
      colnames(thisBranchCatalog)[6] = "thresh"
      colnames(thisBranchCatalog)[7] = "subTree"
      colnames(thisBranchCatalog)[8] = "path"
      
      ## Prepare the tabulation functions
        ## Tabulation limits are [0,1]
      if( length( which( thisBranchCatalog$len > 1 ) ) != 0 ) stop("A branch in the tree has a length beyond the tabulation limits")      
      prepareFn()

      ## Loop over the loss clades
      lossClade <<- list()
      lossNode  <<- list()
      lossTime  <<- list()

      for( i in 1:nrow( theseLossNodes ) )
        {
        lossNode[[i]] <<- as.character(theseLossNodes$lossNode[i])
        lossTime[[i]] <<- as.double(theseLossNodes$lossTime[i])
        
        ## --------------------------
        ## If the loss point is a node
        if( lossNode[[i]] %in% thisTree$node.label )
          {
            ## Extract the loss clade
            clade = extract.clade( phy=thisTree, node=lossNode[[i]] )

            ## Tips of that clade
            clade.tips   = clade$tip.label
            clade.nodes  = clade$node.label
            
            ## Inverse covariance matrix for the clade
            ## Ref: "Averaging Correlated data", Michael Schmelling, CERN, CERN PPE 94-185, 1994
            thisvcvClade[[i]] = solve(vcv(clade))
              
            ## If the loss occurs before the node, consider the parent node
            if( lossTime[[i]] > 0 ){
              cat( paste( "Loss occured", lossTime[[i]], "before", lossNode[[i]], "\n" ) )
              lossNode[[i]] <<- tips.nodes[ Ancestors( x = thisTree, node = which( tips.nodes == lossNode[[i]] ), "parent" ) ]
              cat( paste( "Switching the loss node to its parent:", lossNode[[i]], "\n\n" ) )
              cladeAnc = extract.clade( phy=thisTree, node=lossNode[[i]] )             
            }
            
            ## Loop over the branches of the loss clade
            for( branch in c( clade.tips, clade.nodes ) )
              {
                ibranch = which( thisBranchCatalog$br == branch )
                ## SubTree
                thisBranchCatalog$subTree[ibranch] = i
                ## Status                
                thisBranchCatalog$stat[ibranch] = 0
                ## Branch length
                len = thisBranchCatalog$len[ibranch]
                ## Mean percentID
                thisBranchCatalog$mPid[ibranch] = fnSimSelection( len )
                ## Neutral weight
                thisBranchCatalog$w[ibranch] = fnBrWeights( len )
                ## List of the nodes along the path until the loss node
                ## Get the branch index from the original node vector - necessary for Ancestors
                if( ! branch %in% thisTree$tip.label ) next
                itips.nodes = which( tips.nodes == branch )
                ## Ancestors of that branch
                listPath = tips.nodes[ Ancestors( x = thisTree, node = itips.nodes, "all" ) ]
                ## Only the ancestors within the loss clade
                listPath = listPath[ which( listPath %in% clade$node.label ) ]
                ## If the loss starts at the loss node, don't consider the loss node
                if( lossTime[[i]] <= 0 ) listPath = listPath[ -which( listPath == lossNode[[i]] ) ]
                ## Include the branch also
                listPath = c( listPath, branch )
                thisBranchCatalog$path[[ibranch]] = list(listPath)
              }
            
          }
        else
            {        
              ## ----------------------------
              ## If the loss point is a tip
              clade = ""
              thisvcvClade[[i]] = ""
              itip = which( thisBranchCatalog$br == lossNode[[i]] )
              ## SubTree
              thisBranchCatalog$subTree[itip] = i
              ## Status        
              thisBranchCatalog$stat[itip] = 0
              ## Mean percentID
              len = thisBranchCatalog$len[itip]
              thisBranchCatalog$mPid[itip] = fnSimSelection( len )
              ## Branch length, if -1 the full branch length is returned
              if( lossTime[[i]] > 0 ) thisBranchCatalog$len[itip] = lossTime[[i]]       
              ## Neutral weight
              thisBranchCatalog$w[itip] = fnBrWeights( len )
            }
          lossClade[[i]] <<- clade
        }

      ## Check the branches under selection
      selecBranches = which( thisBranchCatalog$stat == -1 )
      thisBranchCatalog$w[ selecBranches ]    = fnBrWeights( thisBranchCatalog$len[ selecBranches ] )
      thisBranchCatalog$mPid[ selecBranches ] = fnSimSelection( thisBranchCatalog$len[ selecBranches ] )
      thisBranchCatalog$stat[ selecBranches ] = 1

      ## List of branches
      listBranches <<- unique( thisBranchCatalog$br )
      nBranches    <<- length( listBranches )
    }

  ## Gather the clades under selection
     if( in_analysis == "branch" || in_analysis == "all" )
       if( in_collapseClades != "neutralOnly"  && in_collapseClades != "no" ) thisvcvCladeSel = gatherConservedBranches( thisBranchCatalog, thisTree )
  
  ## Overview of the analysis
  options(width=300)
  print( thisBranchCatalog )

  ## Return the clades matrices and the branch catalog
  workOutput = list( thisvcvClade, thisBranchCatalog, thisvcvCladeSel)
  return( workOutput )
}


#######################################################
## Gather the conserved branches  as clades
#######################################################
gatherConservedBranches = function( thisBranchCatalog, thisTree ){

  tips.nodes = c( thisTree$tip.label, thisTree$node.label )
  
  ## =====================================
  ## Gather only the user defined clades
  if( in_collapseClades != "all" ){
    ## Ancestors are split by commas
    collapseAnc = unlist( strsplit( x=in_collapseClades, split="," ) )
    ## Tag for the conserved clade
    consClade = 1000
    ## Loop over the clades to collapse    
    for( i in 1:length(collapseAnc) ){
      anc = collapseAnc[i]
      ## Clade of that ancestor
      clade = extract.clade( phy=thisTree, node=anc )
      ## Loop over the tips of that clade
      for( tip in clade$tip.label ){
          itip0 = which( tips.nodes == tip )
          itip  = which( thisBranchCatalog$br == tip )
          ## All ancestors of that tip
          allAnc = tips.nodes[ unlist( Ancestors( node=itip0, x=thisTree, type="all" ) ) ]
          allAnc = allAnc[ which( allAnc %in% clade$node.label ) ]
          ## Path from the child node to the tip
          listPath = c( allAnc, tip )
          thisBranchCatalog$path[[itip]] <<- list(listPath)
      }
      ## Loop over all descendants
      descendantsClade = c( clade$tip.label, clade$node.label, anc )
      for( des in descendantsClade ){
        j = which( thisBranchCatalog$br == des )
        if( length(j) == 0 ) next
        if( thisBranchCatalog[ j, "subTree" ] == -1 ) thisBranchCatalog[ j, "subTree" ] <<- consClade                
      }
      consClade = consClade + 1
      ## Inverse covariance matrix for the clade
      vcvCladeSel[[ consClade - 1000 ]] <<- solve(vcv(clade))
    }       
  }

  ## =======================================================================
  ## Collapse the clades for which all ancestors have the same phenotype
    if( in_collapseClades == "all" ){
      ## All branches under selection
      consBranches = subset( thisBranchCatalog, stat == 1 )
      ## All tips under selection
      consTips = consBranches[ which( consBranches$br %in% thisTree$tip.label ), ]
      ## Tag for the conserved clade
      consClade = 1000
      ## List of the ancestors under selection
      listConsAnc = c()
      ## Loop over the tips under selection
      for( tip in consTips$br )
        {
          ## Indices of that tip in the two lists
          itip0 = which( tips.nodes == tip )
          itip  = which( thisBranchCatalog$br == tip )
          ## All ancestors of that tip
          allAnc = tips.nodes[ unlist( Ancestors( node=itip0, x=thisTree, type="all" ) ) ]
          ## Loop over these ancestors
          for( i in 1:length(allAnc) )
            {
              node = allAnc[i]
              ## Index of that node
              inode = which( tips.nodes == node )
              ## All descendants of that node
              listDescendants = tips.nodes[ unlist( Descendants( node=inode, x=thisTree, type="all" ) ) ]          
              nDescendants = length( listDescendants )
              ## All descendants under selection for that node
              consDescendants = listDescendants[ which( listDescendants %in% consBranches$br ) ]
              nConsDescendants = length( consDescendants )
              ## If some of the descendants are not neutral, that is the node we are looking for
              if( nConsDescendants != nDescendants )
                {
                  ## If that ancestor is a tip
                  if( i == 1 ) break
                  ## Consider the child node
                  lastConsNode = allAnc[i-1]
                  inode2 = which( tips.nodes == lastConsNode )
                  ## Clade of the child node
                  clade = extract.clade( phy=thisTree, node=lastConsNode )            
                  ## Path from the child node to the tip
                  listPath = allAnc[ which( allAnc %in% clade$node.label ) ]
                  listPath = c( listPath, tip )
                  thisBranchCatalog$path[[itip]] <<- list(listPath)
                  ## Set the subTree tag to all descendants of the child node including himself
                  consDescendantsClade = c( clade$tip.label, clade$node.label, lastConsNode )
                  for( consDescendants in consDescendantsClade ){
                    j = which( thisBranchCatalog$br == consDescendants )
                    if( length(j) == 0 ) next
                    if( thisBranchCatalog[ j, "subTree" ] == -1 ) thisBranchCatalog[ j, "subTree" ] <<- consClade                
                  }
                  ## Update the tag for the selection clade only if a new clade
                  if( ! lastConsNode %in% listConsAnc ){
                    consClade = consClade + 1
                    listConsAnc = c( listConsAnc, lastConsNode )
                    ## Inverse covariance matrix for the clade
                    vcvCladeSel[[ consClade - 1000 ]] = solve(vcv(clade))
                  }              
                  break
                }
            }
        }
    }
  ## ============================================================
  ## Check that the subTrees do not overlap
  
  ## List of clades
  listClades = unique(thisBranchCatalog$subTree)
  listClades = listClades[ which( listClades != -1 )]
  ## Combination of clades
  combClades = combn( listClades, m=2 )
  ## Loop over the combination
  for( comb in 1:ncol(combClades) ){
    clade1   = combClades[1,comb]
    clade2   = combClades[2,comb]
    catalog1 = subset( thisBranchCatalog, subTree == clade1 )
    catalog2 = subset( thisBranchCatalog, subTree == clade2 )
    path1 = unlist(catalog1$path)
    path1 = path1[ which( path1 != -1 )]
    path2 = unlist(catalog2$path)
    path2 = path2[ which( path2 != -1 )]
    ## Common ancestors?
    if( length( intersect( path1, path2 ) ) > 0 ){
      cat(paste( "Error in the loss node, two subTrees overlap:", clade1, "and", clade2,"\n") )
      stop()
    }
  }
return( vcvCladeSel )
}

