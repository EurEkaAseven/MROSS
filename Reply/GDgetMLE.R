GDgetMLE3 <- function(x, y, w = 1, alpha = 0.99, learning_rate0 = 1) {
  d <- ncol(x)
  beta <- rep(0, d)
  loop <- 1
  Loop <- 1000
  msg <- "NA"
  while (loop <= Loop) {
    learning_rate <- learning_rate0 * loop ^ (-alpha)
    pr <- c(1 - 1 / (1 + exp(x %*% beta)))
    S <- colSums((y - pr) * w * x)
    beta.new <- beta + learning_rate * S
    tlr <- sum((beta.new - beta)^2)
    #tlr <- sum(S^2)
    beta <- beta.new
    if (tlr < 0.000001) {
      msg <- "Successful convergence"
      break
    }
    if (loop == Loop) {
      warning("Maximum iteration reached")
    }
    loop <- loop + 1
  }
  list(par = beta, message = msg, iter = loop)
}
