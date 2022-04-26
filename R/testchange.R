#' testchange
#' @export testchange
#' @description This function identifies geographic areas with significant change over time.
#' @param data a numeric matrix, each row representing a time-series
#' and each column representing a time point
#' @param time defines the time sequence
#' @param perm if perm = 'TRUE', a permutation is performed
#' @param nperm number of permuations
#' @param numclust defines the number of clusters for the parallel processing
#' @param topF number of top F values to be selected when perm = 'FALSE'
#' @return Output if perm = 'TRUE' is a list of three items:
#' \itemize{
#'  \item perm.F - F values obtained from permutation tests
#'  \item p.values - p-values obtained from permutation tests
#'  \item p.adjusted - p-values adjusted by Benjamini-Hochberg method
#'  }
#'  Output if perm = 'False' is a list of three items:
#' \itemize{
#'  \item obs.F - conventional F-statistic values
#'  \item sig.change - areas with significant change over time pattern selected by top F-statistic values
#'  \item sel.F - top F-statistic values selected
#'  }
#' @details number of permutations of >=10,000 is ideal
#' @examples
#' # This is an example not using the permutation approach
#'
#' opioid_data_noNA <- opioidData[complete.cases(opioidData), ] #remove NAs
#'
#' mydata <- as.matrix(opioid_data_noNA[,4:18])
#'
#' testchange_results <- testchange(data=mydata,perm=FALSE,time=seq(1,15,1))
#' @references 1. Song, J., Carey, M., Zhu, H., Miao, H., Ram´ırez, J. C., & Wu, H. (2018). Identifying the dynamic gene regulatory network during latent HIV-1 reactivation using high-dimensional ordinary differential equations. International Journal of Computational Biology and Drug Design, 11,135-153. doi: 10.1504/IJCBDD.2018.10011910.
#' 2. Wu, S., & Wu, H. (2013). More powerful significant testing for time course gene expression data using functional principal component analysis approaches. BMC Bioinformatics, 14:6.
#' 3. Carey, M., Wu, S., Gan, G. & Wu, H. (2016). Correlation-based iterative clustering methods for time course data: The identification of temporal gene response modules for influenza infection in humans. Infectious Disease Modeling, 1, 28-39.
#' @export testchange
#' @importFrom splines bs
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster clusterCall stopCluster
#' @importFrom stats p.adjust

testchange <- function(data,time,perm=FALSE,nperm=100,numclust=4, topF= 500){

  results <- list()

  # First, calculate an F value for each geographic area
  obs.F <- rep(NA,nrow(data))
  p.values <- rep(NA,nrow(data))

  for(i in 1:nrow(data)){
    obsData.centered <- as.numeric(scale(as.numeric(data[i,]), scale = FALSE))

    fit1 <- lm(obsData.centered~1)
    fit2 <- lm(obsData.centered~bs(time,3,Boundary.knots=range(time)))

    RSS0 <- sum(obsData.centered^2)
    RSS1 <- sum((obsData.centered-predict(fit2))^2)
    dfnum <- fit1$df.residual
    dfdenom <- fit2$df.residual
    Fval <- ((RSS0-RSS1)/dfnum)/(RSS1/dfdenom)
    obs.F[i] <- Fval
  }

  if(perm==FALSE){
    ordered_data <- data[order(-obs.F),] # order data descending
    ordered_F <- obs.F[order(-obs.F)]
    sig.change <- ordered_data[1:topF, ] # subset data with top N F values
    sel.F <- ordered_F[1:topF]

    #save results
    results$obs.F <- obs.F
    results$sig.change <- sig.change
    results$sel.F <- sel.F
  }

  if(perm==TRUE){
    print("Performing permutation test..")
    #
    # Permutation for proposed method (spline smoothing model)
    #

    cl<-makeCluster(numclust)
    registerDoParallel(cl)
    clusterCall(cl, function(x) .libPaths(x), .libPaths())

    permFval <- foreach::foreach(j=1:nperm, .combine='cbind', .packages="foreach") %dopar% {

      foreach(i=1:nrow(data), .combine='c') %do% {
        obsData.centered <- sample(scale(as.numeric(data[i,]), scale = FALSE)) # This shuffles the data at each CBSA

        # Same procedure to calculate the Fval, but now with the shuffled data
        fit1 <- lm(obsData.centered~1)
        fit2 <- lm(obsData.centered~bs(time,3,Boundary.knots=range(time)))
        RSS0 <- sum(obsData.centered^2)
        RSS1 <- sum((obsData.centered-predict(fit2))^2)
        dfnum <- fit1$df.residual
        dfdenom <- fit2$df.residual
        Fval <- ((RSS0-RSS1)/dfnum)/(RSS1/dfdenom)
        Fval

      }
    }

    stopCluster(cl)
    results$perm.F <- permFval #Saving results
    p.values <- rep(NA,nrow(data)) # Saving p-values here

    # Calculating p-values based on permutation
    for(k in 1:nrow(data)){# k <- 1
      p.values[k] <- sum(permFval >= obs.F[k])/(nperm*nrow(data))
    }

    # Obtain p-values adjusted for mutiple testing
    p.adjusted <- p.adjust(p.values, method = "BH", n = length(p.values))
    results$p.values <- p.values
    results$p.adjusted <- p.adjusted

  }
  return(results)
}
