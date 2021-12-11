Lambda.at.time.t <- function(lambda, time){
  K <- length(lambda)
  interval.len <- 1/K
  num.intervals <- floor(time/interval.len)
  if (num.intervals == 1){
    Lambda.at.t <- lambda[1]*time
  } else if (num.intervals == K) {
    Lambda.at.t <- sum(lambda)*interval.len
  } else {
    Lambda.at.t <- sum(lambda[1:num.intervals]*interval.len) + 
      lambda[num.intervals+1]*(time - interval.len*num.intervals)
  }
  return(Lambda.at.t)
}