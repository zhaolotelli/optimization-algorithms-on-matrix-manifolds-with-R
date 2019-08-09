ParProj <- function(U, V) {
  n = nrow(U)
  P <- diag(1, n) - U %*% solve(t(V) %*% U) %*% t(V)
  return((t(P) + P)/2)
}

vec <- function(M) {
  return(t(t(as.vector(M))))
}

Zcompute <- function(A, B, Y) {
  n <- nrow(Y)
  d <- ncol(Y)
  P <- ParProj(B%*%Y, Y)
  A1 <- P %*% A
  A2 <- P %*% B
  A3 <- solve(t(Y) %*% B %*% Y) %*% (t(Y) %*% A %*% Y)
  A4 <- P %*% A %*% Y
  Id <- diag(1, d)
  M <- rbind((Id %x% A1) - (t(A3) %x% A2), Id %x% t(Y))
  b <- rbind(-vec(A4), matrix(0, nrow = d*d, ncol = 1))
  z <- solve(t(M) %*% M) %*% (t(M) %*% b)
  Z <- matrix(z, nrow = n, ncol = d)
  return(Z)
}

RNMalg <- function(A, B, Y0, itemax = 1e3, delta = 1e-4) {
  n <- nrow(Y0)
  p <- ncol(Y0)
  if(qr(Y0)$rank!=p)
    stop("Y0 not full rank")
  Y <- Y0
  Y.new <- Y0 + 1
  dista <- sum(diag(t(Y.new - Y) %*% (Y.new - Y)))
  ite = 1
  while(dista >= delta & ite <= itemax) {
    Z <- Zcompute(A, B, Y)
    qrY.new <- qr(Y + Z)
    Y.new <- qr.Q(qrY.new)
    dista <- sum(diag(t(Y.new - Y) %*% (Y.new - Y)))
    cat(ite, dista, '\n', sep = '\t')
    ite = ite + 1
    Y <- Y.new
  }
  return(Y)
}

A <- diag(c(1:100))
B <- diag(1, 100)
Y0 <- matrix(rnorm(500), nrow = 100, ncol = 5)
Y <- RNMalg(A, B, Y0)
sum(diag(solve(t(Y) %*% B %*% Y) %*% t(Y) %*% A %*% Y))