simcurve <- function(numareas=c(300,300,300),p,type,normerr=.01){
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

