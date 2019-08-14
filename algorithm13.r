tr <- function(M) {
  return(sum(diag(M)))
}

f <- function(A, Y) {
  return(tr(solve(t(Y) %*% Y) %*% t(Y) %*% A %*% Y))
}

Grad <- function(A, Y) {
  n <- nrow(Y)
  AY <- A %*% Y
  P <- diag(1, n) - Y %*% solve(t(Y) %*% Y) %*% t(Y)
  return(2*P%*%AY)
}

UpdateAlpha <- function(A, Y, eta, alpha = 1, beta = 0.5, sigma = 0.5) {
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
  return(alpha * beta^m)
}

TransVector <- function(Y, epsilon, eta) {
  n <- nrow(Y)
  X <- Y + epsilon
  P <- diag(1, n) - X %*% solve(t(X) %*% X) %*% t(X)
  return(P%*%eta)
}

ComputeBeta <- function(A, Y, Y.new, alpha, eta) {
  top <- tr(t(Grad(A, Y.new))%*%(Grad(A, Y.new) - TransVector(Y.new, alpha*eta, Grad(A, Y))))
  bottom <- tr(t(Grad(A, Y))%*%Grad(A, Y))
  return(top/bottom)
}

GCGalg <- function(A, Y0, ite.max = 100, delta = 1e-5) {
  Y = Y0
  eta = -Grad(A, Y0)
  j = 0
  while(j <= ite.max) {
    alpha <- UpdateAlpha(A, Y, eta)
	qrY.new <- qr(Y + alpha*eta)
	Y.new <- qr.Q(qrY.new)
	if(tr(t(Y.new - Y) %*% (Y.new - Y)) <= delta) {
	  return(Y.new)
	}
	beta <- ComputeBeta(A, Y, Y.new, alpha, eta)
	eta.new <- -Grad(A, Y.new) + beta*TransVector(Y.new, alpha*eta, eta)
	j <- j + 1
	Y <- Y.new
	eta <- eta.new
  }
  return(Y.new)
}

A <- diag(c(1:5, 106:200))
Y0 <- matrix(rnorm(500), nrow = 100, ncol = 5)
Y <- GCGalg(A, Y0)
f(A, Y)