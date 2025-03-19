IRLSgetMLE <- function(x, y, w = 1) {
  d <- ncol(x)
  beta <- rep(0, d)
  loop  <- 1
  Loop  <- 100
  msg <- "NA"
  while (loop <= Loop) {
    eta <- x %*% beta
    pr <- c(1 - 1 / (1 + exp(eta)))
    wdiag <- pr * (1-pr)
    z <- eta + (y-pr)/wdiag
    
    H <- t(x) %*% (wdiag * w * x)
    S <- t(x) %*% (wdiag * w * z)
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
    beta.new <- shs
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
  list(par=beta, message=msg, iter=loop, Infor = H)
}

