#' @export
#' @importFrom("stats","aggregate","as.dist","cor","cutree","hclust","lm","predict")
#' @importFrom(factoextra,fviz_nbclust,hcut)
#' @importFrom(splines,bs)
#' @importFrom(ihclust,simcurve)
#' @importFrom(ggplot2,labs)

geoclust <- function(data, smooth = TRUE, cor_criteria = 0.75, max_iteration = 100, verbose = TRUE){
  time <- seq(1,ncol(data),1)
  pred.data <- matrix(rep(NA,nrow(data)*ncol(data)),ncol=ncol(data))

  # alpha defines the correlation criteria
  alpha <- cor_criteria
  if(!((alpha>0) & (alpha<1))){
    print("Warning: Correlation criteria should be greater than 0 and less than 1. The default value alpha=0.75 is used." )
    alpha=0.75
  }

  # Whether Perform Clustering on the Fitted Curves or on the Observed Data
  if(smooth == TRUE){
    for(i in 1:nrow(data)){
      obsData.centered <- as.numeric(scale(data[i,], scale = FALSE))
      fit_sc <- lm(obsData.centered~bs(time,3,Boundary.knots=range(time)))
      pred.data[i,] <- predict(fit_sc)
    }

  }else{
    for(i in 1:nrow(data)){
      obsData.centered <- as.numeric(scale(data[i,], scale = FALSE))
      pred.data[i,] <- obsData.centered
    }

  }

  ##
  ## Clustering
  ##

  # keep a record of iteration steps
  iteration_steps <- list()
  # keep a record of the unique clusters after merging
  memb.updates <- list()

  ##
  ## Step 1: Initialization
  ##

  ## Determine the number of clusters

  # Apply Average Silhouette Method to determine the optimal number of clusters
  Silhouette_clust <- fviz_nbclust(pred.data, FUNcluster = hcut, method = "silhouette", k.max = 10,
                                   diss = as.dist(1-cor(t(pred.data), method="spearman")))+
    labs(subtitle = "Silhouette method")
  Silhouette_result <- Silhouette_clust$data
  # optimal number of cluster: the one with max Silhouette
  Silhouette_max_cluster<-as.numeric(Silhouette_result$clusters[which.max(Silhouette_result$y)])
  opt_cluster <- Silhouette_max_cluster

  ## Hierarchical clustering using optimal number of clusters

  hc <- hclust(as.dist(1-cor(t(pred.data), method="spearman")), method="complete")
  memb <- cutree(hc, k=opt_cluster)
  memb_step1 <- memb

  ## Here is when the steps 2 and 3 are iterated

  # num.not.met will be greater than 0 when there are more than 1 pair of cluster exemplars with correlation > alpha
  num.not.met <- 1
  iteration <- 0

  if(num.not.met > 0 ){

    while(num.not.met > 0){

      iteration <- iteration + 1
      if(verbose == TRUE){
        print(iteration)
      }


      ##
      ## Step 2: Merging
      ##
      exemplars <- aggregate(x = pred.data, by = list(memb), FUN = mean)

      pairwiseCor <- cor(t(exemplars[,-1]),method='spearman')
      pairwiseCor[upper.tri(pairwiseCor, diag = TRUE)] <- NA
      pairs <- which(pairwiseCor > alpha, arr.ind = TRUE)
      # View(pairs)

      memb_labels <- exemplars[,1]

      # Now merging
      num_to_merge <- dim(pairs)[1]

      if(num_to_merge>0){
        for(k in num_to_merge:1){ # k <- 1
          memb[ memb == memb_labels[pairs[k,1]] ] <- memb_labels[pairs[k,2]]
        }
      }

      rm(pairwiseCor); rm(exemplars); rm(pairs); rm(num_to_merge); rm(memb_labels)

      # keep a record of the unique clusters after merging
      memb.updates[[iteration]] <- sort(unique(memb))

      # If the cluster labels are identical to the one in the previous iteration, then stop
      if(iteration>10){
        no.change <- ifelse(any(duplicated(memb.updates)==TRUE),1,0)
        if (no.change ==1) {print(paste("Unique cluster labels after the merging step was identical to one or more prior iterations"))}
        if (no.change ==1) {print(memb.updates[[iteration-1]])}
        if (no.change ==1) {print(memb.updates[[iteration]])}
        if (no.change ==1) break
        rm(no.change)
      }

      iteration_steps[[1]] <- any(duplicated(memb.updates)==TRUE)
      iteration_steps[[2]] <- memb.updates

      ##
      ## Step 3: Pruning
      ##

      # First, recalculate exemplars
      exemplars <- aggregate(x = pred.data, by = list(memb), FUN = mean)
      memb_labels <- exemplars[,1]
      ind.cor <- rep(NA,dim(pred.data)[1])
      for(m in 1:dim(pred.data)[1]){
        ind.cor[m] <- cor(pred.data[m,], t(exemplars[which(exemplars[,1]==memb[m]),-1]),method="spearman")
      }

      # Any member with correlation <= alpha, assign to new single clusters
      num.to.prune <- sum(ind.cor <= alpha)
      memb.num.max <- max(unlist(memb.updates))

      if(num.to.prune>0){
        for(p in 1:num.to.prune){ # p <- 1
          memb[ which(ind.cor <= alpha)[p] ] <- memb.num.max + p
        }
      }

      rm(exemplars); rm(ind.cor); rm(num.to.prune); rm(memb.num.max); rm(memb_labels)

      ##
      ## Recalculate Exemplars
      ##

      # First, recalculate exemplars
      exemplars <- aggregate(x = pred.data, by = list(memb), FUN = mean)

      pairwiseCor <- cor(t(exemplars[,-1]),method='spearman')
      pairwiseCor[upper.tri(pairwiseCor, diag = TRUE)] <- NA
      pairs <- which(pairwiseCor > alpha, arr.ind = TRUE)

      rm(pairwiseCor); rm(exemplars);
      num.not.met <- dim(pairs)[1]

      if(verbose == TRUE){
        print(sort(unique(memb)))
      }

      if (iteration == max_iteration) break

    }
  }

  result = list("Cluster_Label" = memb,
             "Num_Iterations" = iteration,
             "Unique_Clusters_in_Iteration" = iteration_steps[[2]])
}

