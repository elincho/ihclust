\name{IHC_cluster}
\alias{IHC_cluster}

\title{
Iterative hierarchical clustering (IHC)
}
\description{
Identify inhomogeneous clusters using iterative hierarchical clustering (IHC) method
}
\usage{
IHC_cluster(data, smooth = TRUE, cor_criteria = 0.75, max_iteration = 100, verbose = TRUE)
}

\arguments{
  \item{data}{ numeric matrix. Each row represents a time-series and each column represents a time point.}
  \item{smooth}{ logical value. If TRUE, a smooth function is applied before clustering.}
  \item{cor_criteria}{ float, a pre-specified correlation criteria.}
  \item{max_iteration}{ integer, maximum number of iterations.}
  \item{verbose}{ logical value. If TRUE, the result of progress is printed.}
}
\details{
The IHC algorithm implements the three steps as outlined below.

First, the Initialisation step clusters the data using hierarchical clustering

Second, cluster centres are obtained as an average of all the data points in the cluster. The Merging step considers each of the cluster centres (exemplars) as ‘new data point’, and use the same procedure described in the Initialisation step to merge the exemplars into a new set of clusters.

Third, the Pruning step streamlines the clusters and removes inconsistencies by reassessing
the cluster membership by each data point.

}
\value{
Output from the function is a list of three items:
  \item{Cluster_Label}{the cluster label for each data point.}
  \item{Num_Iterations}{total number of iterations.}
  \item{Unique_Clusters_in_Iteration}{unique clusters in each iteration.}
}
%%\references{
%% ~put references to the literature/web site here ~
%%}
%%\author{
%%  ~~who you are~~
%%}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
\examples{
library(IHC1)
simulated_data <- simcurve(numareas=c(300,300,300),p=0.01,type="fixed", normerr=0.1)
data = as.matrix(simulated_data[[1]])
IHC_cluster(data, smooth = TRUE, cor_criteria = 0.75, max_iteration = 100, verbose = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
