f <- function(A, Y) {
  return(sum(diag(solve(t(Y) %*% Y) %*% t(Y) %*% A %*% Y)))
}

updateEta <- function(A, Y) {
  return(-2*(A %*% Y - Y %*% solve(t(Y) %*% Y) %*% t(Y) %*% A %*% Y))
}

updateY <- function(A, Y, alpha, beta, sigma, eta) {
  m = 0
  tay <- sigma * alpha * beta^m * sum(diag(t(eta) %*% eta))
  Y.new <- Y + alpha * beta^m * eta
  err <- f(A, Y) - f(A, Y.new)
  while(err < tay) {
    #print(m)
    m = m + 1
    Y.new <- Y + alpha * beta^m * eta
    err <- f(A, Y) - f(A, Y.new)
    tay <- sigma * alpha * beta^m * sum(diag(t(eta) %*% eta))
  }
  return(Y.new)
}

AlsRqG <- function(A, Y0, alpha = 1, beta = 0.5, sigma = 0.5, itemax = 1e3, delta = 1e-4) {
  n <- nrow(Y0)
  p <- ncol(Y0)
  if(qr(Y0)$rank!=p)
	  stop("Y0 not full rank")
  Y <- Y0
  Y.new <- Y0 + 1
  dista <- sum(diag(t(Y.new - Y) %*% (Y.new - Y)))
  ite = 1
  while(dista >= delta & ite <= itemax) {
	  eta <- updateEta(A, Y)
	  Y.new <- updateY(A, Y, alpha, beta, sigma, eta)
	  qrY.new <- qr(Y.new)
	  Y.new <- qr.Q(qrY.new)
	  dista <- sum(diag(t(Y.new - Y) %*% (Y.new - Y)))
	  cat(ite, dista, '\n', sep = '\t')
	  ite = ite + 1
	  Y <- Y.new
  }
  return(Y.new)
}

A1 <- diag(1:100)
A2 <- diag(c(1, 102:200))
A3 <- diag(c(1:5, 106:200))
Y0 <- matrix(rnorm(500), nrow = 100, ncol = 5)
Y1 <- AlsRqG(A1, Y0)
Y2 <- AlsRqG(A2, Y0)
Y3 <- AlsRqG(A3, Y0)
f(A1, Y1)
f(A2, Y2)
f(A3, Y3)