calc_independent_likelihoods_on_each_branch_prebyte <- function(phy2, Qmat, cluster_already_open=NULL, Qmat_is_sparse=FALSE)
{
  # phy2 must have been reordered, e.g.
  phy2 = reorder(phy2, "pruningwise")
  
  if (Qmat_is_sparse==FALSE)
  {
    # Get probmats for each branch, put into a big array
    # Create empty array to store results
    #independent_likelihoods_on_each_branch = array(0, dim=c(nrow(Qmat), ncol(Qmat), length(phy2$edge.length)))
    
    independent_likelihoods_on_each_branch = vector("list", length(phy2$edge.length))
    tmpmatrix = matrix(data=0, nrow=nrow(Qmat), ncol=ncol(Qmat))
    for (m in 1:length(phy2$edge.length))
    {
      independent_likelihoods_on_each_branch[[m]] = tmpmatrix
    }
    # Calculate the conditional likelihoods for each branch
    # dgexpv NOT ALLOWED when you have a null range state
    # (maybe try very very small values here)
    
    # clusterApply and other multicore stuff (e.g. doMC) are apparently dangerous on R.app
    if (!is.null(cluster_already_open))
    {
      # 
      if (.Platform$GUI == "AQUA")
      {
        cat("In calc_loglike_sp(), cluster_already_open=", cluster_already_open, " which means you want to calculate likelihoods on branches using a multicore option.\n", sep="")
        cat("But .Platform$GUI='AQUA', which means you are running the Mac GUI R.app version of R.  Parallel multicore functions, e.g. as accessed via \n", sep="")
        cat("library(parallel), are apparently dangerous/will crash R.app (google multicore 'R.app').  So, changing to cluster_already_open=NULL.\n", sep="")
        cluster_already_open=NULL
      } # END if (.Platform$GUI == "AQUA")
    } # END if (!is.null(cluster_already_open))
    
    
    # Run on the cluster of nodes, if one is open
    # clusterApply etc. appear to NOT work on R.app
    if (!is.null(cluster_already_open))
    {
      # mcmapply
      #library(parallel)
      #independent_likelihoods_on_each_branch = mcmapply(FUN=expokit_dgpadm_Qmat, Qmat=list(Qmat), t=phy2$edge.length, transpose_needed=TRUE, SIMPLIFY="array", mc.cores=Ncores)
      independent_likelihoods_on_each_branch = clusterApply(cl=cluster_already_open, x=phy2$edge.length, fun=expokit_dgpadm_Qmat2, Qmat=Qmat, transpose_needed=TRUE)
    } else {
      # Not parallel processing
      #independent_likelihoods_on_each_branch = mapply(FUN=expokit_dgpadm_Qmat, Qmat=list(Qmat), t=phy2$edge.length, transpose_needed=TRUE, SIMPLIFY="array")
      independent_likelihoods_on_each_branch = mapply_likelihoods(Qmat, phy2, transpose_needed=TRUE)
      #independent_likelihoods_on_each_branch
    }
  } else {
    errortxt = paste("\n\nError in calc_independent_likelihoods_on_each_branch():\n Cannot be used with sparse matrix exponentiation! (Qmat cannot be sparse)\n\n", sep="")
    cat(errortxt)
    stop("\nStopping on error.\n")
  } # END if (Qmat_is_sparse==FALSE)
  
  return(independent_likelihoods_on_each_branch)
}