# Xavier Prudent, 2015

## Functions for the global analysis (perfect match, GLS)
source(paste(scriptDir, "/forwardGenomics_globalAnalysis.R", sep=""))

## Functions for the local analysis (branch method, the experimental subTree method)
source(paste(scriptDir, "/forwardGenomics_localAnalysis.R", sep=""))


#######################################################
## Loop over the elements
#######################################################
loopElements = function(){
  cat(style("\nLoop over the genomic elements...\n",fg="red",font.style=c("bold") ) )
  
  ## Loop over the elements
  for( elID in in_elements$V1 )
    {
      ## Element ID
      elID =as.character( elID )
      cat(style( paste("Element", elID, "\n"),fg="blue",font.style="bold" ) )

      ## Check if there are enough species with data
      if( in_analysis == "GLS" || in_analysis == "all" ) if( enoughSpecies( elID ) == FALSE ) next

      ## Local analysis
      if( in_analysis == "branch" || in_analysis == "all" )
        {
          ## Does every leaf have percent id data
          anyMissing = missingData( elID )

          ## If percent id is missing for some branches
          if( anyMissing == 2 ) next
          if( anyMissing == 1 ){

            if( verbose ) cat("Missing data (genome gaps) were found for that element\n")

            ## Remove branches with missing data from the tree
            prunedTree <<- removeBranches()

            ## Get phenotypes and loss nodes
            prunedLosses <<- propagatePheno( prunedTree )

            ## Fill the corresponding branch catalog
            treePreparationOutput = prepareTree( prunedTree, prunedLosses )
            prunedvcvClade <<- treePreparationOutput[[1]]
            prunedCatalog  <<- treePreparationOutput[[2]]
            prunedvcvCladeSel  <<- treePreparationOutput[[3]]

            ## Run the tree based analysis
            treeBased.res = treeBased_analysis( elID, prunedCatalog, prunedTree, prunedvcvClade, prunedvcvCladeSel )
          }

          ## If no missing data
          if( anyMissing == 0 ) treeBased.res = treeBased_analysis( elID, branchCatalog, tree, vcvClade, vcvCladeSel )

          ## Return false for no data
          #if( treeBased.res == FALSE ) next        
        }

      ## Global analysis
      if( in_analysis == "GLS" || in_analysis == "all" ){

        ## GLS analysis
        GLS_analysis( elID )

        ## Perfect match analysis    
        perfectMatch_analysis( elID )
      }

      ## Save the result in a text file
      saveResult( elID )
    }
}


######################################################
## Are there missing datapoints in the local percent ID values ?
## 0: no missing data
## 1: missing data
## 2: too much missing data
######################################################
missingData = function( elID ){

  ## Local percentID data for that element
  el.branch = subset( in_localPid, id == elID )
  
  ## Keep only the leaves
  el.branch = subset( el.branch, br %in% tree$tip.label )

  ## Which leaves don't have data
  missingTips <<- tree$tip[ which( ! tree$tip.label %in% el.branch$br ) ]
  if( length( missingTips ) == 0 ) return( 0 ) else
  {
    if( verbose ){
      cat("\nNo local percent id found for the following leaves:\n")
      print(missingTips)
    }
    
    ## Number of loss and conserved species missing  
    nLossMissing = 0
    nConsMissing = 0
    for( spc in missingTips ){
      phenoValue = listPheno$pheno[ listPheno$species == spc ]
      if( phenoValue == 0 ) nLossMissing = nLossMissing + 1
      if( phenoValue == 1 ) nConsMissing = nConsMissing + 1
    }
  
  ## At least 2 loss species and 2 conserved species should remain
    nLossRemain = nLoss - nLossMissing
    nConsRemain = nCons - nConsMissing
    if( nLossRemain >= 2 && nConsRemain >= 2 )
      return( 1 )
    else{
      if( verbose ) cat("\nNot enough species left to proceed\n")
      return( 2 )
    }
  }
}


######################################################
## Use tree_doctor from the phast package to remove branches with missing data
######################################################
removeBranches = function(){

  ## Use tree_doctor to
  ## (1) remove the leaves with missing data
  ## (2) update the nodes names accordingly
  cmd = "tree_doctor -n -a -p "
  for( spc in missingTips ){
    if( spc != tail(missingTips,n=1) ) cmd = paste( cmd, spc, ",", sep="" )
    else cmd = paste( cmd, spc, sep="" )
  }
  cmd = paste( cmd, in_tree, sep=" " )

  ## Command system
  prunedTree.newick = system( cmd, intern = TRUE )

  ## New pruned tree
  prunedTree = read.tree( text = prunedTree.newick )
  return( prunedTree )
}


######################################################
## Global analysis: are there enough species?
######################################################
enoughSpecies = function( elID ){
  
    ## Index of the element
    if( (elID %in% rownames(in_globalPid_pheno)) == FALSE ){
      cat("    > Not found\n")
      return( FALSE )
    }
    el = which( rownames(in_globalPid_pheno) == elID )
    
    ## At least 4 species with not NA pheno or phenotype
    remSpecies <<- colnames(in_globalPid_pheno)[which( !is.na(in_globalPid_pheno["pheno",]) & !is.na(in_globalPid_pheno[el,]) )]
    if( verbose ) cat(paste("Number of species:", length(remSpecies), "\n"))
    if( length(remSpecies) < 2 )
      {
        cat(paste("    > Only", length(remSpecies), "species left, element skipped\n") )
        return( FALSE )
      }else
        {        
          ## Remove from the tree the species with NA pheno or phenotype
          noPercentID = colnames(in_globalPid_pheno)[which( is.na(in_globalPid_pheno["pheno",]) || is.na(in_globalPid_pheno[el,]) )]
          finalTree   <<- drop.tip( tree, noPercentID )    
          
          ## (ranked) Dataset with not NA pheno or phenotype
          subData             <<- as.data.frame( t( in_globalPid_pheno[ c("pheno",elID), which( !is.na(in_globalPid_pheno["pheno",]) & !is.na(in_globalPid_pheno[el,]) )]))          
          subData[ ,"pheno" ] <<- as.numeric(as.character(subData[ ,"pheno" ]))
          subData[ ,elID ]    <<- as.numeric(as.character(subData[ ,elID ]))

          ## (raw) Dataset with not NA pheno or phenotype
          subData.raw             <<- as.data.frame( t( in_rawGlobalPid_pheno[ c("pheno",elID), which( !is.na(in_rawGlobalPid_pheno["pheno",]) & !is.na(in_rawGlobalPid_pheno[el,]) )]))          
          subData.raw[ ,"pheno" ] <<- as.numeric(as.character(subData.raw[ ,"pheno" ]))
          subData.raw[ ,elID ]    <<- as.numeric(as.character(subData.raw[ ,elID ]))

          ## Add the species as an additional column
          subData$species    <<- rownames(subData)
          row.names(subData) <<- NULL
          subData.raw$species    <<- rownames(subData.raw)
          row.names(subData.raw) <<- NULL
          
          if( verbose ) print(subData.raw)

          ## At least 1 species in each group (lost/conserved)
          grp1 <<- subset(subData.raw, pheno == 1 )
          n1.spec <<- nrow( grp1 )
          grp0 <<- subset(subData.raw, pheno == 0 )
          n0.spec <<- nrow( grp0 )
          if( n1.spec == 0 || n0.spec == 0 ){
            if( n1.spec == 0 ) cat("    > No species with phenotype 1 left\n", sep="" )
            if( n0.spec == 0 ) cat("    > No species with phenotype 0 left\n", sep="" )
            return( FALSE )
          }        
    
    ## Does the element have data for all required loss events
    nIndepLoss = 0
    for( i in unique(branchCatalog$subTree) )
      {
        if( i < 0 ) next
        if( length(which( subset( branchCatalog, subTree == i )$br %in% subData$species )) != 0 ) nIndepLoss = nIndepLoss + 1
        else cat(paste("    > Sub loss tree ", i, " missing\n", sep="" ))            
      }
    if( nIndepLoss < indepLoss ) return( FALSE )
        }
    
    ## If no problem
    return(TRUE)
}


######################################################
## Save all the results for that element
######################################################

saveResult = function ( elID ){

  ## Basic info on the element
 if( in_analysis == "GLS" )       basicInfo = c( elID, n0.spec, n1.spec )
 if( in_analysis == "branch"  )   basicInfo = c( elID, n0.branch, n1.branch )
 if( in_analysis == "all" )       basicInfo = c( elID, n0.spec, n1.spec, n0.branch, n1.branch )
 
  ## Tree-based method  
  if( in_analysis == "branch" || in_analysis == "all" ){
     if( in_collapseClades != "no" ) Res = c( wPearson, wPearsonPval )
     BpB.Res = c( BpB.wPearson, BpB.wPearsonPval  )
     if( in_collapseClades != "no" ) localRes = c( Res, BpB.Res ) else localRes = c( BpB.Res )
  }
  ## Global method 
  if( in_analysis == "GLS" || in_analysis == "all" ) globalRes = c( perfectMatch, slopePval )

  ## Write the result
  if( in_analysis == "GLS" )    RES = c( basicInfo, globalRes )
  if( in_analysis == "branch" ) RES = c( basicInfo, localRes )
  if( in_analysis == "all" )    RES = c( basicInfo, globalRes, localRes )
  write( RES, file=as.character(nameOutput), ncolumns=length(RES), append=TRUE )

}

