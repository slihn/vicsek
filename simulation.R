# ----------------------------------------------------- #
# global parameters

N  = 1000 # number of particles
L  = 100  # size of the box    
v0 = 2  # speed
r0 = 1  # radius of neighbors
eta = 0.8 # intensity of noise
noise = 2 # mode of noise - 1: intrinsic, 2: extrinsic
S = 1000  # number of simulations

# ----------------------------------------------------- #
# init
theta <- runif(N, 0, 2*pi)
v <- v0 * exp(theta*1i)
x <- lapply(1:N, function(a) { complex(real=runif(1,0,L), imaginary=runif(1,0,L)) })
in0 <- list(x=x, v=v)

# ----------------------------------------------------- #
# the model
iterate_model <- function(input) {
  x <- input$x
  v <- input$v
  theta2 <- (1:N)*0
  for (i in 1:N) {
    xi <- x[[i]]
    d <- lapply(x, function(a) { abs(a-xi) }) # distance
    u <- which(d <= r0) # within radius r0
    eps <- runif(1, 0, 2*pi) # noise in angle
    if ( noise == 1 ) { 
      v_avg <- mean(v[u]) 
      theta2[i] <- Im(log(v_avg/abs(v_avg))) + eta*eps
    }
    if ( noise == 2 ) { 
      v_avg <- mean(v[u]) + eta*exp(eps*1i) 
      theta2[i] <- Im(log(v_avg/abs(v_avg)))
    }
  }
  v2 <- v0 * exp(theta2*1i)
  x2 <- mapply(function(a,b) { 
    c <- a+b
    # boundary condition
    if (Re(c) > L) { c <- c-L }
    if (Re(c) < 0) { c <- c+L }
    if (Im(c) > L) { c <- c-L*1i }
    if (Im(c) < 0) { c <- c+L*1i }
    c
  }, x, v2)
  list(x=x2, v=v2)
}
# ----------------------------------------------------- #
# run the simulation
for (k in 1:S) {
  rs <- iterate_model(in0)
  x2 <- rs$x
  v2 <- rs$v
  psi <- abs(mean(v2))/v0 # order parameter
  plot(sapply(x2,Re), sapply(x2,Im))
  in0 <- rs
  print(c(k, psi))
}
# ----------------------------------------------------- #
