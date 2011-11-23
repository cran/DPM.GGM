plotHeatmap <- function(A,cuts=50,...)
{
	if ( all (A >= 0 ))
		levelplot(A,col.regions=grey(seq(1,0,length=100)),cuts=cuts,...)
	else {
		max <- max(A)
		min <- min(A)
		range <- max - min
		length.red <- abs(round(min/range*100))
		length.blue <- 100 - length.red
		red <- c(rep(1,length.red),seq(1,0,length=length.blue))
		blue <- c(seq(0,1,length=length.red),rep(1,length.blue))
		green <- c(seq(0,1,length=length.blue),seq(1,0,length=length.red))

		rgb <- rgb(red,green,blue)
		levelplot(A,col.regions=rgb,cuts=cuts,...)
	}	
}
