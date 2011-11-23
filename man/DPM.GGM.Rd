\name{DPM.GGM}
\alias{DPM.GGM}
\title{Dirichlet process mixtures of Gaussian graphical models}
\description{
This function fits a Dirichlet process mixtures of Gaussian
graphical models to a set of observations X.
}
\usage{
DPM.GGM ( X,burn = 1e3,reps = 1e4, print.every = 1e2, slice=200)
}
\arguments{
\item{X}{n x r matrix. Observations}
\item{burn}{Number of burn in steps for gibbs sampling}
\item{reps}{Number of repetition steps for gibbs sampling}
\item{print.every}{After how many steps should the progress be reported?}
\item{slice}{Integer or logical. If n > slice, slice sampler is beeing used. Set to \code{TRUE}
if you want to use slice sampler anyway, set to \code{FALSE} if you don't want to use
the slice sampler.
}
}
\value{
	\item{cluster}{How many times observation i and j were in the same cluster}
	\item{edge}{How many times edge j was in the graph associated with observation i}
}
\examples{
data("exchange_rates")

X <- exchange_rates
p <- dim(X)[2]

## normalization
for(i in 1:p)X[,i] <- (X[,i] - mean(X[,i])) / sd(X[,i])

set.seed(1)
reps <- 1e2
burn <- 1e2

a <- DPM.GGM(X,burn,reps,print.every=10)
\dontrun{
	plot(a)
}
}
\seealso{
\code{\link{iHMM}}
\code{\link{plot.DPM.GGM}}
}
\references{
Abel Rodriguez, Alex Lenkoski, and Adrian Dobra: Sparse
covariance estimation in heterogeneous samples. Electronical
Journal of Statistics Volume 5 (2011), 981-1014.
\url{http://projecteuclid.org/DPubS?service=UI&version=1.0&verb=Display&handle=euclid.ejs/1316092866}
}
