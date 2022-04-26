#' simcurve
#' @description This function generates two kinds of datasets.
#' 1. Randomly generates curves with change/no change.
#' 2. Generates true curves assumed from fixed coeffecients with some random noise.
#' @param numareas number of areas to generate
#' @param p proportion of the areas that have significant change
#' @param type type of curves generated
#' @param normerr standard deviation of the Normal distribution (with mean zero) of which the coefficients are generated
#' @details If type = "random", the function generates curves with change/no change.
#' If type = "fixed", the function generates true curves assumed from fixed coefficients with some random noise.
#' If numareas is not specified, it is assumed as a vector of c(300,300,300).
#' If normerr is not specified, it is assumed as a value of 0.01. It is ignored when type= "random".
#' @examples
#' mydata_ran <- simcurve(numareas = c(300, 300, 300), p=0.01, type="random")
#'
#' mydata_fixed <- simcurve(numareas = c(300, 300, 300), p=0.01, type="fixed", normerr = 0.1)
#' @return Output from the function is a list of two items:
#' \itemize{
#'  \item data - simulated data
#'  \item parameters - parameters used to generate the data}
#' @export simcurve
#' @importFrom stats rnorm
simcurve <- function(numareas=c(300,300,300),p=0.05,type,normerr=.1){
  sumnum <- sum(numareas)
  time=seq(1,52,1)
  ntime <- length(time)
  simulated <- matrix(rep(NA,sumnum*ntime),ncol=ntime)

  ### type
  if (type == "random"){

    #
    # This part sets up random coefficients for the areas with true change
    # Coefficients are generated from standard normal
    #
    numChange <- sumnum*p ## numareas -> sumnum starting from this line
    coefMat <- matrix(rep(NA,sumnum*3),ncol=3)
    for(i in 1:numChange){
      coefMat[i,] <- c(rnorm(n=1, mean=0, sd=4),
                       rnorm(n=1, mean=0, sd=2),
                       rnorm(n=1, mean=0, sd=1))
    }

    #
    # For areas with no true change over time, we just specify coefficent of zero
    #
    for(i in (numChange+1):sumnum){coefMat[i,] <- c(0,0,0)}
  }

  else if (type == "fixed"){
    # Set up true coefficients for 3 clusters
    cluster1 <- c(4,4,4)
    cluster2 <- c(-4,-4,-4)
    cluster3 <- c(4,-4,-4)

    # Just assigning true coefficients to each area
    coefMat <- matrix(rep(NA,sumnum*3),ncol=3)
    for(i in 1:numareas[1]){coefMat[i,] <- cluster1}
    for(i in (numareas[1]+1):(numareas[1]+numareas[2])){coefMat[i,] <- cluster2}
    for(i in (numareas[1]+numareas[2]+1):sumnum){coefMat[i,] <- cluster3}
  }

  else {
    stop("'type' must be one of 'random', or 'fixed'")
  }

  #
  # This is the function to generate the true curves, based on coefficients
  #
  truth <- function(x, a, b, coef){
    value <- coef[1]*(-sqrt(2/(b-a))*cos(2*pi*(x-a)/(b-a))) +
      coef[2]*sqrt(2/(b-a))*sin(2*pi*(x-a)/(b-a)) +
      coef[3]*(-sqrt(2/(b-a))*cos(4*pi*(x-a)/(b-a)))
    value
  }

  #
  # Now generating data, based on the coefficients defined above
  #
  for(j in 1:sumnum){
    trueData <- rep(NA,length(time))
    for(i in 1:length(time)){
      trueData[i] <- truth(time[i],time[1],time[length(time)],coefMat[j,])
    }
    # We add random noise to all data (irrespective of true change or not)
    simulated[j,] <- trueData + rnorm(length(trueData),mean=0,sd=sqrt(normerr))
  }
  par <- c(sumnum,p,type)
  sim_list <- list(data = simulated, parameters = par)
}

