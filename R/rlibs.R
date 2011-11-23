run.iHMM <- function(X,burn = 1e3,reps = 1e4, print.every = 1e2, predict.ahead = FALSE, all.K = FALSE)
  {
    data(stirling)
    for(i in 1:(dim(stirling)[1] - 1))
      {
        stirling[i, (i + 1):dim(stirling)[1]] <- 0
      }
    
    n <- dim(X)[1]
    p <- dim(X)[2]
    ee <- choose(p,2)
    l <- .C("run_iHMM_routine", as.double(X), as.integer(n), as.integer(p), as.integer(burn), as.integer(reps),
            as.double = t(stirling), as.integer(dim(stirling)[1]), cluster = as.integer(rep(0,n * n)),
            edge = as.integer(rep(0,n * ee)), alpha = as.double(rep(0, reps)), alpha0 = as.double(rep(0, reps)), cluster.total = as.integer(rep(0,reps)),
            as.integer(predict.ahead), x.predicted = as.double(rep(0, reps * p)), mu.predicted = as.double(rep(0, reps * p)), K.predicted = as.double(rep(0,reps * p * p)),
            as.integer(print.every), as.integer(all.K), K.all = as.double(rep(0, n * p * p)))
    results <- NULL
    results$cluster <- matrix(l$cluster, n,n,byrow = TRUE)
    results$edge <- matrix(l$edge, n, ee, byrow = TRUE)
    if(predict.ahead) results$x.predicted <- matrix(l$x.predicted, reps, p, byrow = TRUE)
    if(predict.ahead) results$mu.predicted <- matrix(l$mu.predicted, reps, p, byrow = TRUE)
    if(predict.ahead) results$K.predicted <- array(l$K.predicted, dim = c(p, p, reps), data = l$K.predicted)
    if(all.K) results$K.all = array(l$K.all, dim = c(p,p,n), data = l$K.all)
    results$alpha <- l$alpha
    results$alpha0 <- l$alpha0
    results$L <- l$cluster.total
    return(results)
  }


run.dmp <- function(X,burn = 1e3,reps = 1e4, print.every = 1e2)
  {

    n <- dim(X)[1]
    p <- dim(X)[2]
    ee <- choose(p,2)
    l <- .C("run_dmp_routine", as.double(X), as.integer(n), as.integer(p), as.integer(burn), as.integer(reps), cluster = as.integer(rep(0,n * n)), edge = as.integer(rep(0,n * ee)), as.integer(print.every))
    results <- NULL
    results$cluster <- matrix(l$cluster, n,n,byrow = TRUE)
    results$edge <- matrix(l$edge, n, ee, byrow = TRUE)
    return(results)
  }

run.dmp.slice <- function(X,burn = 1e3,reps = 1e4, print.every = 1e2)
  {

    n <- dim(X)[1]
    p <- dim(X)[2]
    ee <- choose(p,2)
    l <- .C("run_dmp_slice_routine", as.double(X), as.integer(n), as.integer(p), as.integer(burn), as.integer(reps), cluster = as.integer(rep(0,n * n)), edge = as.integer(rep(0,n * ee)), as.integer(print.every))
    results <- NULL
    results$cluster <- matrix(l$cluster, n,n,byrow = TRUE)
    results$edge <- matrix(l$edge, n, ee, byrow = TRUE)
    return(results)
  }

test <- function(X)
  {
    .C("test",as.double(X),as.integer(dim(X)[2]), as.integer(dim(X)[1]))
  }

