#usethis::use_package(splines)
#usethis::use_package(foreach)
#usethis::use_package(doParallel)
#usethis::use_package(parallel)

#' testchange
#'
#' @param data
#' @param time
#' @param perm
#' @param nperm
#' @param numclust
#'
#' @return
#' @export
#'
#' @examples
testchange <- function(data,time,perm=FALSE,nperm=100,numclust=64){

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

  if(perm==TRUE){
    print("Performing permutation test..")
    #
    # Permutation for proposed method (spline smoothing model)
    #

    cl<-makeCluster(numclust)
    registerDoParallel(cl)
    clusterCall(cl, function(x) .libPaths(x), .libPaths())

    permFval <- foreach(j=1:nperm, .combine='cbind', .packages="foreach") %dopar% {

      foreach(i=1:nrow(data), .combine='c') %do% {

        library(splines)
        time=seq(1,15,1)
        library(foreach)
        library(doParallel)

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




##
##

#opioid_data <- read.csv("/Users/elincho/Desktop/Elin Cho/Research/FDA/data/CDC/DispenseRates_counties.csv")

# View(opioid_data)
# dim(opioid_data)
#  View(opioid_data[,4:18])

#opioid_data_noNA <- opioid_data[complete.cases(opioid_data), ]
# dim(opioid_data_noNA)


#opioid_data_NA <- opioid_data[!complete.cases(opioid_data), ]
# View(opioid_data_NA)


#mydata <- as.matrix(opioid_data_noNA[,4:18])
#  View(opioid_data_subset)
#  dim(opioid_data_subset)
# View(data)


#testchange_results <- testchange(data=mydata,perm=FALSE,nperm,time=seq(1,15,1))
# names(testchange_results)
# View(testchange_results$obs.F)

#testchange_results <- testchange(data=mydata,perm=TRUE,nperm=10,time=seq(1,15,1))
# ideally, would need about 10,000 permutations
# names(testchange_results)
# View(testchange_results$obs.F)
# View(testchange_results$p.adjusted.SP)


