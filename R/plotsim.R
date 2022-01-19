plotsim <- function(curve){
  if (curve$parameters[3] == "random"){
    mysimdata <- curve$data
    tiff('Fig1_simulate_CBSA_with_change.tiff', units="in", width=7, height=3, res=300, compression = 'lzw')

    par(mfrow=c(1,3))
    time=seq(1,52,1)
    plot(time, mysimdata[1,],type='l', lwd = 0.2, ylim=c(-3.5,3.5),
         ylab="Standardized number of prescriptions dispensed",
         xlab="",
         main="All CBSAs Over Time \n (n=900)",
         cex.lab=1.2, cex.axis=0.8, cex.main=1, col = "red")
    for(i in 2:45){lines(time,mysimdata[i,],ylim=c(-3.5,3.5), col = "red")}
    for(i in 46:900){lines(time,mysimdata[i,],ylim=c(-3.5,3.5), col = "deepskyblue")}

    plot(time, mysimdata[1,],type='l',ylim=c(-3.5,3.5),
         ylab="",
         xlab="Time (Weeks)",
         main="CBSAs with Change Over Time \n (n=45)",
         cex.lab=1.2, cex.axis=0.8, cex.main=1, col = "red")
    for(i in 2:45){lines(time,mysimdata[i,],ylim=c(-3.5,3.5), col = "red")}

    plot(time, mysimdata[46,],type='l',ylim=c(-3.5,3.5),
         ylab="",
         xlab="",
         main="CBSAs with No Change \n (n=855)",
         cex.lab=1.2, cex.axis=0.8, cex.main=1, col = "deepskyblue")
    for(i in 47:900){lines(time,mysimdata[i,],ylim=c(-3.5,3.5), col = "deepskyblue")}
    dev.off()
  }
  else if (curve$parameters[3] == "fixed"){
    simulated_clusters <- curve$data
    tiff('Fig2_simulation_clusters.tiff', units="in", width=6, height=2, res=300, compression = 'lzw')
    par(mfrow=c(1,3))
    time <- seq(1,52,1)
    plot(time, simulated_clusters[1,],type='l',ylim=c(-3.5,3.5),
         ylab="Standardized number of prescriptions dispensed",
         xlab = "",
         #xlab="Time",
         main="Simulated data for cluster 1",
         cex.lab=0.8, cex.axis=0.8, cex.main=1)
    for(i in 2:300){lines(time,simulated_clusters[i,],ylim=c(-3.5,3.5))}

    plot(time, simulated_clusters[301,],type='l',ylim=c(-3.5,3.5),
         #ylab="Standardized number of prescriptions dispensed",
         ylab = "",
         xlab="Time (Weeks)",
         main="Simulated data for cluster 2",
         cex.lab=0.8, cex.axis=0.8, cex.main=1)
    for(i in 301:600){lines(time,simulated_clusters[i,],ylim=c(-3.5,3.5))}

    plot(time, simulated_clusters[601,],type='l',ylim=c(-3.5,3.5),
         ylab = "",
         xlab = "",
         #ylab="Standardized number of prescriptions dispensed",
         #xlab="Time",
         main="Simulated data for cluster 3",
         cex.lab=0.8, cex.axis=0.8, cex.main=1)
    for(i in 601:900){lines(time,simulated_clusters[i,],ylim=c(-3.5,3.5))}
    dev.off()
  }}
