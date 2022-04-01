#' testchange
#' @export testchange
#' @description This function identifies geographic areas with significant change over time.
#' @param data a numeric matrix, each row representing a time-series
#' and each column representing a time point
#' @param time defines the time sequence
#' @examples
#' # This is an example not using the permutation approach
#' opioid_data_noNA <- opioidData[complete.cases(opioidData), ] #remove NAs
#' mydata <- as.matrix(opioid_data_noNA[,4:18])
#' testchange_results <- testchange(data=mydata,time=seq(1,15,1))
#' # This is an example using simulated data
#' mydata <- simcurve(numareas = c(300, 300, 300), p=0.01, type="fixed", normerr = 0.01)
#' testchange_results <- testchange(data=mydata$data,time=seq(1,52,1))
#' @references 1. Song, J., Carey, M., Zhu, H., Miao, H., Ram´ırez, J. C., & Wu, H. (2018). Identifying the dynamic gene regulatory network during latent HIV-1 reactivation using high-dimensional ordinary differential equations. International Journal of Computational Biology and Drug Design, 11,135-153. doi: 10.1504/IJCBDD.2018.10011910.
#' 2. Wu, S., & Wu, H. (2013). More powerful significant testing for time course gene expression data using functional principal component analysis approaches. BMC Bioinformatics, 14:6.
#' 3. Carey, M., Wu, S., Gan, G. & Wu, H. (2016). Correlation-based iterative clustering methods for time course data: The identification of temporal gene response modules for influenza infection in humans. Infectious Disease Modeling, 1, 28-39.
#' @export testchange
#' @importFrom splines bs
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster clusterCall stopCluster

testchange <- function(data,time){

  results <- list()

  # First, calculate an F value for each geographic area
  obs.F <- rep(NA,nrow(data))

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

  # Save the observed F value
  results$obs.F <- obs.F
  return(results)
}
