# Xavier Prudent, 2015

######################################################
## Initialization of lists
######################################################
listInitial = function(){

  ## Distance to the expected value
  totBranchesDist  <<- c()
  neutBranchesDist <<- c()
  selBranchesDist  <<- c()

  ## Weight 
  brWeights <<- c()
  neutBranchesWeights <<- c()

  ## Status
  selBranchesStatus <<- c()
  brStatus <<- c()

  ## For the breanch-per-branch method
  ## Distance to the expected value
  BpB.totBranchesDist<<- c()
  BpB.neutBranchesDist <<- c()

  ## Weight 
  BpB.brWeights <<- c()
  BpB.neutBranchesWeights <<- c()

  ## Status
  BpB.brStatus <<- c()

}


######################################################
## branch-based analysis
######################################################
treeBased_analysis = function( elID, thisCatalog, thisTree, thisvcvClade, thisvcvCladeSel ){

  listInitial()
  allOutputToNA=FALSE

   ## ===================
  ## Branch per branch analysis
  if( verbose ){
    cat(style("\n ==================================\n", font.style="bold") )
    cat(style("     Branch per branch analysis\n", font.style="bold") )
    cat(style(" ==================================\n", font.style="bold") )
  }
  branch.subTree <<- thisCatalog
  BpB.anl = na.omit(computeExpDist( elID, thisTree ))

  ## Total size of the data for that element
  n = nrow( BpB.anl )
  
  ## If no branches for that element
  if( is.null(n) || n == 0 ){
    cat("No branches data for that element\n")
    allOutputToNA = TRUE
   # return(FALSE)
  }
  else{

    ## Reject the element if a big indel occured on a conserved branch
    if( is.logical(BpB.anl) ) return(BpB.anl)

    ## Conserved branches
    BpB.anl.sel  = na.omit(subset( BpB.anl, stat == 1 ))
    n1.branch <<- nrow( BpB.anl.sel )
    if( !is.null(n1.branch) ){
      
      ## If no branches
      if( n1.branch == 0 ){
        if( verbose ) cat("No conserved branches for the analysis, all outputs set to NA\n")
        allOutputToNA = TRUE
      }

      ## Branches under selection
      if( in_collapseClades == "no" ){
        branch.subTree <<- na.omit(subset( thisCatalog, subTree == -1 ))
        selec.subTree = computeExpDist( elID, thisTree )
        selBranchesDist      <<- c( selBranchesDist, selec.subTree$Dist )
      }
      BpB.brWeights       <<- c( BpB.brWeights, rep( times = n1.branch, x = 1 ) )
      BpB.brStatus        <<- c( BpB.brStatus, rep( times = n1.branch, x = 1 ) )
      BpB.totBranchesDist  <<- c( BpB.totBranchesDist, BpB.anl.sel$Dist )

      ## Neutral
      BpB.anl.neut  = na.omit(subset( BpB.anl, stat == 0 ))
      n0.branch <<- nrow( BpB.anl.neut )
      ## If no neutral branches
      if( !is.null(n0.branch) ){      
        if( n0.branch == 0 ){
           if( verbose ) cat("No neutral branches for the analysis, all outputs set to NA\n")
           allOutputToNA = TRUE        
        }          
        BpB.brWeights       <<- c( BpB.brWeights, BpB.anl.neut$w )
        BpB.brStatus        <<- c( BpB.brStatus, rep( times = n0.branch, x = 0 ) )
        BpB.totBranchesDist  <<- c( BpB.totBranchesDist, BpB.anl.neut$Dist )
        BpB.neutBranchesWeights  <<- c( BpB.neutBranchesWeights, BpB.anl.neut$w )
        BpB.neutBranchesDist      <<- c( BpB.neutBranchesDist, BpB.anl.neut$Dist )     
      }
    }
  }
  
  ## ===================
  ## Subtree analysis
  if( in_collapseClades != "no" ){
  if( verbose ){
    cat(style("\n ==================================\n", font.style="bold") )
    cat(style("        Subtree analysis\n", font.style="bold") )
    cat(style(" ==================================\n", font.style="bold") )
  }
  
  for( clade in unique(thisCatalog$subTree) )
    {
      ## Set of branches in that subTree
      branch.subTree <<- subset( thisCatalog, subTree == clade )
      
      ## Branches under selection
      if( clade == -1 )        
        {
          ## Compute the distance to the expected value
          selec.subTree = computeExpDist( elID, thisTree )
          ## Reject the element is a big event occured on a conserved branch
          if( is.logical(selec.subTree) ){
            if( verbose ) cat("Some large deletion occured on a conserved branch, all outputs set to NA\n")
            allOutputToNA = TRUE
          }
          else{
            ## Size of the data for that element
            n = nrow( selec.subTree )
            if( n == 0 ){
              if( verbose ) cat("No branches for the analysis, all outputs set to NA\n")
              allOutputToNA = TRUE
            }              
            ## Set the quantities for these branches
            brWeights           <<- c( brWeights, rep(times=n,x=1) )
            brStatus            <<- c( brStatus, rep(times=n,x=1) )
            totBranchesDist      <<- c( totBranchesDist, selec.subTree$Dist )
            selBranchesDist      <<- c( selBranchesDist, selec.subTree$Dist )
            selBranchesStatus   <<- c( selBranchesStatus, selec.subTree$stat )                           
          }
        }

      ## SubTrees under selection
      if( clade >= 1000 )
        {
          ## Compute the distance to the expected value
          sel.subTree = computeExpDist.clade( elID, clade, thisTree )
          ## Size of the data for that element
          n = nrow( selec.subTree )
          if( !is.null(n) ){
            ## Make sure there are the same species in the data and in the covariance matrix
            sel.subTree = removeSpecies( clade, sel.subTree, thisvcvClade, thisvcvCladeSel )
            ## Average of the distance to exp value accounting for the relatedness
            ## If no species left
            if( n == 0 ){
              if( verbose ) cat("No branches for the analysis, all outputs set to NA\n")
              allOutputToNA = TRUE
              meanDist = "NA"
              meanW = "NA"
              meanL = "NA"
            }else{
              meanDist  = sum( vcvClade.el %*% sel.subTree$Dist ) / sum( vcvClade.el )
              ## Simple average of the lengths and weights of the tips
              meanW = mean( sel.subTree$w )
              meanL = mean( sel.subTree$len )
            }          
            ## Print out
            if( verbose ){
              cat(paste("    > Mean Dist to expected value :", meanDist, "\n" ) )
              cat(paste("    > Mean length                 :", meanL, "\n" ) )
              cat(paste("    > Mean weight                 :", meanW, "\n" ) )
            }          
            ## Save the results for further correlation analysis
            brWeights            <<- c( brWeights, meanW )
            brStatus             <<- c( brStatus, 1 )
            totBranchesDist       <<- c( totBranchesDist, meanDist )
          }
        }
      
      ## Neutral subTrees and branches
      if( clade > -1 && clade < 1000 )
        {        
          ## Count tips and nodes
          listBr = thisCatalog$br[ which( thisCatalog$subTree == clade ) ]
          nBr    = length( listBr )
          nTips  = length( which( listBr %in% thisTree$tip.label ) )
          
          ## Only one tip in the clade
          if( nTips == 1 && nBr == nTips ){
            oneTip.subTree = computeExpDist( elID, thisTree )
            n = nrow( oneTip.subTree )
            if( n == 0 ){
              if( verbose ) cat("No branches for the analysis, all outputs set to NA\n")
              allOutputToNA = TRUE
            }
            meanDist    = oneTip.subTree$Dist
            meanW = oneTip.subTree$w
            meanL = oneTip.subTree$len
            if( n == 0 ){
              meanDist = "NA"
              meanW = "NA"
              meanL = "NA"
            }else{
              meanDist    = oneTip.subTree$Dist
              meanW = oneTip.subTree$w
              meanL = oneTip.subTree$len
            }
          }
          
          ## Several tips
          if( nTips > 1 ){            
            manyTips.subTree = computeExpDist.clade( elID, clade, thisTree )
            n = nrow( manyTips.subTree )
            if( n == 0 ){
               if( verbose ) cat("No branches for the analysis, all outputs set to NA\n")
              allOutputToNA = TRUE
            }
            ## Remove from the covariance matrix species without data for this element
            tipsOnly.subTree = removeSpecies( clade, manyTips.subTree, thisvcvClade, thisvcvCladeSel )
            ## Average of the distance to exp value accounting for the relatedness
            ## If no species left
            if( n == 0 ){
              meanDist = "NA"
              meanW = "NA"
              meanL = "NA"
            }else{
              meanDist  = sum( vcvClade.el %*% tipsOnly.subTree$Dist ) / sum( vcvClade.el )
              ## Simple average of the lengths and weights of the tips
              meanW = mean( tipsOnly.subTree$w )
              meanL = mean( tipsOnly.subTree$len )
            }
          }
          
          ## Print out
          if( verbose ){
            cat(paste("    > Mean dist to expected value: ", meanDist, "\n" ) )
            cat(paste("    > Mean length                :", meanL, "\n" ) )
            cat(paste("    > Mean weight                :", meanW, "\n" ) )
          }
          
          ## Save the results for further correlation analysis
          brWeights            <<- c( brWeights, meanW )
          brStatus             <<- c( brStatus, 0 )
          totBranchesDist       <<- c( totBranchesDist, meanDist )
          neutBranchesWeights  <<- c( neutBranchesWeights, meanW )
          neutBranchesDist      <<- c( neutBranchesDist, meanDist )
        }
    }
}
  
  ## Perform statistical tests on the branch sample
  statTestsOnBranches( allOutputToNA )

  return(TRUE)
}


###################################################################
## Perform statistical tests on the branch sample
## p-values are recalculated to get a one-sided test
###################################################################
statTestsOnBranches = function( allOutputToNA ){

## If no data for the signal branches or species
  if( allOutputToNA ){
    n0.branch       <<- "NA"
    n1.branch       <<- "NA"
    wPearson          <<- "NA"
    wPearsonPval      <<- "NA"
    BpB.wPearson          <<- "NA"
    BpB.wPearsonPval      <<- "NA"
  } else {

     ######################################
     ## Subtree, collapsing clades
     if( in_collapseClades != "no" ){
       ## Weighted Pearson correlation
       res          = wtd.cor( x=totBranchesDist, y=brStatus, weight = brWeights )
       wPearson     <<- res[which(dimnames(res)[[2]]=="correlation")]
       ## t-value
       res.t = as.numeric(unlist(strsplit(summary(res)[16],":"))[2])
       ## Degrees of freedom
       res.n = sum(brWeights) - 2
       ## p-value
       ## 1-sided
       wPearsonPval <<- pt( q = res.t, df = res.n, lower.tail = FALSE )
       ## 2-sided
       ##wPearsonPval <<- res[which(dimnames(res)[[2]]=="p.value")]
     }

     ######################################
     ## Branch per branch, ignoring clades
     ## Weighted Pearson correlation
     res          = wtd.cor( x = BpB.totBranchesDist, y = BpB.brStatus, weight = BpB.brWeights )
     BpB.wPearson     <<- res[which(dimnames(res)[[2]]=="correlation")]
     ## t-value
     res.t = as.numeric(unlist(strsplit(summary(res)[16],":"))[2])
     ## Degrees of freedom
     res.n = sum(BpB.brWeights) - 2
     ## p-value
     ## 1-sided
     BpB.wPearsonPval <<- pt( q = res.t, df = res.n, lower.tail = FALSE )
     ## 2-sided
     ##BpB.wPearsonPval <<- res[which(dimnames(res)[[2]]=="p.value")]
   }

   ## Verbose
   if( verbose ){
      cat("\n")
      cat(style("Results of the correlation analysis:\n", font.style="bold" ) )
      cat("\n")
      if( in_collapseClades != "no" ){
         cat("    > Weighted subtree method:\n")    
         cat(paste("     > p-value =", wPearsonPval, "\n" ))
         cat("\n")
      }
      cat("    > Weighted branch method:\n")    
      cat(paste("     > p-value =", BpB.wPearsonPval, "\n" ))
      cat("\n") 
   }
}


##########################################################################
## Compute the distance to the expected value for a clade
##########################################################################
computeExpDist.clade = function( elID, clade, thisTree ){

  ## Data for all the branches of the clade for that element
  el.branch.subTree = subset( in_localPid, br %in% branch.subTree$br & id == elID )
  ## Keep the branches for which there is data for the element
  branch.subTree    = subset( branch.subTree, branch.subTree$br %in% el.branch.subTree$br )
  ## Merge the data for the element with the characteristics of the branches
  el.branch.subTree = merge( branch.subTree, el.branch.subTree, by="br" )
  ## distance to exp value: el(pid) - mean(pid(sim))
  el.branch.subTree$Dist = as.double(el.branch.subTree$pid) - as.double(el.branch.subTree$mPid)
  ## Subtrees: loop over the tips and sum the distances
  for( i in 1:nrow(el.branch.subTree) ){
    if( el.branch.subTree$br[i] %in% thisTree$node.label ) next
    
    ## Along the path to the tips
    Fi = 0
    Li = 0
    for( tip.anc in unlist(el.branch.subTree$path[i]) ){
      j  = which( el.branch.subTree$br == tip.anc )
      ## Sum the distances
      Fj = el.branch.subTree$Dist[j]
      if( length(Fj) == 0 ) next
      Fi = Fi + Fj
      ## Sum the branch lengths
      Lj = el.branch.subTree$len[j]
      Li = Li + Lj
    }
    ## If the total distance is smaller than -1, set it to -1
    el.branch.subTree$Dist[i] = ifelse( Fi < -1, -1, Fi )
    ## Update the length and the neutral weight
    el.branch.subTree$len[i] = Li
    el.branch.subTree$w[i] = fnBrWeights( Li )
  }
  if( verbose ){
    cat("\n")
    cat(style(paste("    > Clade", clade, "\n" ), font.style="bold") )
    print( el.branch.subTree )
  }
  
  return( el.branch.subTree )
}


##########################################################################
## Compute the distance to the expected value branch per branch
##########################################################################
computeExpDist = function( elID, thisTree ){

  ## Data for all the branches for that element, omit branches with NA
  el.branch = na.omit(subset( in_localPid, br %in% branch.subTree$br & id == elID ))
  
  ## Check for duplicates
  el.branch.check = data.frame( el.branch$br, el.branch$id )
  ndupl = length(which(duplicated(el.branch.check)) == TRUE)
  if( ndupl > 0 ){
    cat(style( "\nERROR: different values were found in the input data file for given branches and element\n",fg="red",font.style="bold" ) )
    break
  }  
  ## Keep the branches for which there is data for the element
  branch.subTree    = subset( branch.subTree, branch.subTree$br %in% el.branch$br )
  
  ## Merge the data for the element with the characteristics of the branches
  el.branch = merge( branch.subTree, el.branch, by="br" )

  ## distance: el(pid) - mean(pid(sim))
  el.branch$Dist = el.branch$pid - el.branch$mPid

  ## Threshold for big events (large indels) on conserved internal branches (not tips)
  bigEvent = subset( el.branch, stat == 1 & pid <= in_thresholdConserved  )
  bigDelEvent = FALSE
  internBranch = FALSE
  if( nrow(bigEvent) != 0 ){
    bigDelEvent = TRUE
  for( i in 1:nrow(bigEvent) )if( bigEvent$br[i] %in% thisTree$node.label ) internBranch = TRUE
  }
  
  if( bigDelEvent && internBranch ){
    cat(paste("Large indel occured on a conserved internal branch:\n") )
    if( verbose ){
      print( bigEvent )
      cat("\n")
    }
    return( FALSE )
  }
  
  if( verbose ){
    cat("\n")
    cat(style("    > Subset of branches\n",font.style="bold" ) )
    print( el.branch )
  }  
  return( el.branch )
}


##############################################################
## Remove from the covariance matrix the species without data
##############################################################
removeSpecies = function( clade, df, thisvcvClade, thisvcvCladeSel ){

  ## Neutral or selection
  if( clade < 1000 ) cladeMatrix = thisvcvClade else cladeMatrix = thisvcvCladeSel
  if( clade >= 1000 ) clade = clade - 999

  ## Species in the matrix
  colMatrix = colnames( cladeMatrix[[clade]] )
  
  ## Species from the matrix that are not in the data list
  l = which( ! colMatrix %in% df$br )
  
  ## Remove these species
  if( length(l) > 0 ) vcvClade.el <<- cladeMatrix[[clade]][ -l, -l ] else vcvClade.el <<- cladeMatrix[[clade]]
  
  ## Nodes in the data list
  l = which( df$br %in% tree$node.label )
  
  ## We average tips only
  if( length(l) > 0 ) df = df[ -l, ]
  
  ## If more than one species left, keep and order the species according to the matrix
  if( length(vcvClade.el) > 1 ){
    colMatrix = colnames( vcvClade.el )
    df = df[ match(colMatrix, df$br), ]
  } 
  return(df)
}


#####################################
## Complementary error function
#####################################
erfc = function(x){
  2 * pnorm(x * sqrt(2), lower = FALSE)
}

