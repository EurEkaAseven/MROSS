getsSVM <- function(x, y, w, C) {
  d <- ncol(x)
  y <- 2*y - 1
  beta <- rep(0, d)
  loop  <- 1
  Loop  <- 100
  msg <- "NA"
  {if(missing(w)){
    while (loop <= Loop) {
      sup_index <- which(1-y*(x %*% beta) > 0)
      x_sup <- x[sup_index,]
      S <- beta + 2*C*t(x_sup) %*% (x_sup %*% beta - y[sup_index])
      H <- 2*C*t(x_sup) %*% x_sup + diag(d)
      tryCatch(
        {shs <- NA
        shs <- solve(H, S) },
        error=function(e){
          cat("\n ERROR :", loop, conditionMessage(e), "\n")})
      if (is.na(shs[1])) {
        msg <- "Not converge"
        beta <- loop <- NA
        break
      }
      beta.new <- beta - shs
      tlr  <- sum((beta.new - beta)^2)
      beta  <- beta.new
      if(tlr < 0.000001) {
        msg <- "Successful convergence"
        break
      }
      if (loop == Loop)
        warning("Maximum iteration reached")
      loop  <- loop + 1
    }
  }
  else{
    while (loop <= Loop) {
      sup_index <- which(1-y*(x %*% beta) > 0)
      x_sup <- x[sup_index,]
      w_sup <- w[sup_index]
      S <- beta + 2*C*t(w_sup*x_sup) %*% (x_sup %*% beta - y[sup_index])
      H <- 2*C*t(x_sup) %*% (w_sup*x_sup) + diag(d)
      tryCatch(
        {shs <- NA
        shs <- solve(H, S) },
        error=function(e){
          cat("\n ERROR :", loop, conditionMessage(e), "\n")})
      if (is.na(shs[1])) {
        msg <- "Not converge"
        beta <- loop <- NA
        break
      }
      beta.new <- beta - shs
      tlr  <- sum((beta.new - beta)^2)
      beta  <- beta.new
      if(tlr < 0.000001) {
        msg <- "Successful convergence"
        break
      }
      if (loop == Loop)
        warning("Maximum iteration reached")
      loop  <- loop + 1
    }
  }}
  
  list(par=beta, message=msg, iter=loop)
}