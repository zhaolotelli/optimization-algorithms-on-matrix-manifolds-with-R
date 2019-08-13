tr <- function(M) {
  return(sum(diag(M)))
}

Normliz <- function(Y, B) {
  BY <- t(Y) %*% B %*% Y
  sv <- svd(BY)
  U <- sv$u
  D <- sv$d
  V <- sv$v
  M <- sv$u %*% diag(sqrt(1/D)) %*% t(sv$v)
  return(Y %*% M)
}

f <- function(A, B, Y) {
  return(sum(diag(solve(t(Y) %*% B %*% Y) %*% t(Y) %*% A %*% Y)))
}

ParProj <- function(U, V) {
  n = nrow(U)
  P <- diag(1, n) - U %*% solve(t(V) %*% U) %*% t(V)
  return((t(P) + P)/2)
}

H <- function(Z, A, B, Y) {
  P <- ParProj(B%*%Y, B%*%Y)
  return(P%*%(A%*%Z - B%*%Z%*%solve(t(Y)%*%B%*%Y)%*%(t(Y)%*%A%*%Y)))
}

TCGmed <- function(A, B, Y, Delta, ite.max = 10) {
  n = nrow(Y)
  p = ncol(Y)
  eta = matrix(0, nrow = n, ncol = p)
  P = ParProj(B%*%Y, B%*%Y)
  r = P %*% A %*% Y
  delta = -r
  j = 0
  while(j <= ite.max) {
    Hy <- H(delta, A, B, Y)
    if(tr(t(delta)%*%Hy)<=0) {
	  a = tr(t(delta)%*%delta)
	  b = 2*tr(t(delta)%*%eta)
	  c = Delta - tr(t(eta)%*%eta)
	  tau = max((-b+sqrt(b^2 - 4*a*c))/(2*a),(-b-sqrt(b^2 - 4*a*c))/(2*a))
	  return(eta+tau*delta)
	}
	alpha = tr(t(r)%*%r)/tr(t(delta)%*%Hy)
	eta.new = eta + alpha*delta
	if(tr(t(eta.new)%*%eta.new)>=Delta) {
	  a = tr(t(delta)%*%delta)
	  b = 2*tr(t(delta)%*%eta)
	  c = Delta - tr(t(eta)%*%eta)
	  tau = max((-b+sqrt(b^2 - 4*a*c))/(2*a),(-b-sqrt(b^2 - 4*a*c))/(2*a))
	  return(eta+tau*delta)
	}
	r.new = r + alpha*Hy
	beta = tr(t(r.new)%*%r.new)/tr(t(r)%*%r)
	delta.new = -r.new+beta*delta
	j <- j+1
	r <- r.new
	delta <- delta.new
	eta <- eta.new
  }
  return(eta)
}

RTRalg <- function(A, B, Y0, Delta.bar, Delta0, rho.bar, ite.max = 100) {
  j = 0
  Delta = Delta0
  Y = Y0
  Y.new = Y0
  while(j <= ite.max) {
    eta <- TCGmed(A, B, Y, Delta)
	Y.new <- Normliz(Y, B)
	meta <- 2*tr(solve(t(Y)%*%B%*%Y)%*%t(eta)%*%A%*%Y) + tr(solve(t(Y)%*%B%*%Y)%*%t(eta)%*%(A%*%eta-B%*%eta%*%solve(t(Y)%*%B%*%Y)%*%t(Y)%*%A%*%Y))
	rho <- (f(A, B, Y) - f(A, B, Y.new))/(-meta)
	if(rho < 1/4) {
	  Delta.new <- 1/4*Delta
	}
	else if(rho > 3/4 & tr(t(eta)%*%eta) = Delta) {
	  Delta.new <- min(2*Delta, Delta.bar)
	}
	else {
	  Delta.new <- Delta
	}
	if(rho > rho.bar) {
	  Y <- Y.new
	}
	Delta <- Delta.new
	j <- j+1
  }
  return(Y)
}