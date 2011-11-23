plot.DPM.GGM <- function ( x, ... )
{
	## normalizing
	reps <- attr(x,"reps")
	cluster <- x$cluster / reps

	plotHeatmap(cluster)
}
