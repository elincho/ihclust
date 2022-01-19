\name{plotsim}
\alias{plotsim}

\title{
function to plot the simulated data
}

\description{
This function can generate line plots based on types of dataset obtained from simcurve().
}

\usage{
plotsim(curve)
}

\arguments{
  \item{curve}{simulated dataset to be used}
}

\details{
For dataset containing radomly generated areas with change, line plots for 1. All areas over time, 2. areas with change, and 3. areas with no change are generated. For dataset containing curves with fixed coefficients and random noise, line plots of simulated data for 3 clusters are generated.
}

