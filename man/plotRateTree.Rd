\name{plotRateTree}
\alias{plotRateTree}
\title{Plot a phylogeny where branches are painted based on rates.}
\description{This function allows users to visualize regions of phylogeny that support higher or lower transition rates for an Mk model than the maximum likelihood solution based on the entire tree.}
\usage{
plotRateTree(tree = NULL,
             rates,
             scaled = F,
             cols = NULL,
             bg = "white",
             edge.width = 1)
}
\arguments{
  \item{tree}{phylo object (strictly bifurcating)}
  \item{rates}{Number of rate classes, should be 3, 5, or 7.}
  \item{scaled}{if TRUE the branches are scaled based on rate class}
  \item{cols}{a vector of colors with a length equal to rates}
  \item{bg}{background color for plot}
  \item{edge.width}{the width of the branches of the phylogeny}
}
\details{}

\value{}
\author{Terrence Sylvester and Heath Blackmon}
\references{http://coleoguy.github.io/}
\examples{}
\keyword{Mk model, trait evolution, visualization, phylogeny}
