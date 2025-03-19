getMLE <- function(x, y, w = 1) {
    d <- ncol(x)
    beta <- rep(0, d)
    loop  <- 1
    Loop  <- 100
    msg <- "NA"
    while (loop <= Loop) {
        pr <- c(1 - 1 / (1 + exp(x %*% beta)))
        H <- t(x) %*% (pr * (1 - pr) * w * x)
        S <- colSums((y - pr) * w * x)
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
        beta.new <- beta + shs
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

