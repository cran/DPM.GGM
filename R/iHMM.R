iHMM <- function ( X,burn = 1e3,reps = 1e4, print.every = 1e2, predict.ahead = FALSE, all.K = FALSE)
{
	results <- run.iHMM(X=X,burn=burn,reps=reps,print.every=print.every,predict.ahead=predict.ahead,all.K=all.K)
	class(results) <- "iHMM"
	attr(results,"reps") <- reps
	return(results)
}
