\name{iHMM}
\alias{iHMM}
\title{Infinite hidden Marcov models}
\description{
This function fits an infitite hidden marcov model to a set of observations X.
}
\usage{
iHMM ( X,burn = 1e3,reps = 1e4, print.every = 1e2, predict.ahead = FALSE, all.K = FALSE)
}
\arguments{
\item{X}{n x r matrix. Observations}
\item{burn}{Number of burn in steps for gibbs sampling}
\item{reps}{Number of repetition steps for gibbs sampling}
\item{print.every}{After how many steps should the progress be reported?}
\item{predict.ahead}{Predict the next cluster the observation will land in?}
\item{all.K}{Return fitted precision matic for all observations (requires considerable RAM)}
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

\dontrun{
a <- iHMM(X,burn,reps,print.every = 10, all.K = FALSE)
plot(a)
}
}
\seealso{
\code{\link{DPM.GGM}}
\code{\link{plot.iHMM}}
}
\references{
Abel Rodriguez, Alex Lenkoski, and Adrian Dobra: Sparse
covariance estimation in heterogeneous samples. Electronical
Journal of Statistics Volume 5 (2011), 981-1014.
\url{http://projecteuclid.org/DPubS?service=UI&version=1.0&verb=Display&handle=euclid.ejs/1316092866}
}
