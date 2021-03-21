\name{treePaintR}
\alias{treePaintR}
\title{Paint phylogeny based on Mk model rate variation}
\description{This function allows users to identify regions of phylogeny that support higher or lower rates than the default background rate estimated by maximum likelihood}
\usage{
treePaintR(tree = NULL, tip_states = NULL, qmat = NULL, iter = NULL, rate.class = NULL, rates = NULL, verbose = T)
}
\arguments{
  \item{tree}{phylo object (strictly bifurcating)}
  \item{tip_states}{vector of tip states in same order as tips of the phylogeny}
  \item{qmat}{matrix of rates that describe the evolution of the trait of interest. These are expected to be the maximum likelihood rates estimated in a package like Geiger, APE, or Diversitree}
  \iter{iter}{the number of iterations to perform. The default is 10,000.}
  \iter{rate.class}{Number of rate classes, should be an odd number greater than 1.}
  \iter{rates}{vector of rates that branches will be scaled by}
  \iter{verbose}{If TRUE then various messages will be printed during the fitting process}
}
\details{}

\value{Returns a list. The first element is a phylo object with an additional element providing a rate class for each branch of the phyogeny.}
\author{Terrence Sylvester and Heath Blackmon}
\references{http://coleoguy.github.io/}
\examples{}
\keyword{Mk model, trait evolution, visualization}
