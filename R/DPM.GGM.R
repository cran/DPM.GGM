DPM.GGM <- function (
	X,
	burn = 1e3,
	reps = 1e4,
	print.every = 1e2,
	slice=200
)
{
	if ( (slice != FALSE && dim(X)[1] > slice ) || slice == TRUE)
	{
		results <- run.dmp.slice (
			X=X,
			burn = burn,
			reps = reps,
			print.every = print.every
		)
	} else {
		results <- run.dmp (
			X=X,
			burn = burn,
			reps = reps,
			print.every = print.every
		)
	}
	class(results) <- "DPM.GGM"
	attr(results,"reps") <- reps
	return(results)
}
