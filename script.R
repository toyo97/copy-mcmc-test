library(tidyverse)
library(gtools) # for dirichlet
library(matrixStats) # for logSumExp trick
library(gridExtra) # for multiple plots on pdf

# for reproducibility
set.seed(42)

load_data <- function(path) {
  return()
}

rmc <- function(x0, length, pi0, pi) {
  S <- 0:(nrow(pi)-1)
  x0 <- sample(S, 1, prob = pi0)
  X <- rep(0, length)
  X[1] <- x0
  for (t in 2:length) {
    X[t] <- sample(S, 1, prob = pi[X[t-1] + 1,])
  }
  return(X)
}

simulate_data <- function(M, B, C) {
  # disproportion parameter
  dir.ralpha <- 10
  
  # initial probability (uniform)
  pi0 <- rep(1/M, M)
    
  # transition probability (dirichlet)
  pi <- matrix(nrow = M, ncol = M)
  for (i in 1:M) {
    dir.alpha.i <- rep(1, M)
    dir.alpha.i[i] <- dir.ralpha
    pi[i,] <- rdirichlet(1, dir.alpha.i)
  }  
  
  # read depth (gamma)
  gamma.a <- 10
  gamma.b <- 1000
  
  s <- rgamma(C, gamma.a, scale = gamma.b)
  
  # copy numbers (markov chain)
  X <- rmc(x0, B, pi0, pi)
  
  # fix per-copy expression
  mu <- rbeta(B, 1, 1000)
  
  # simulate reads (poisson)
  Y <- matrix(nrow = B, ncol = C)
  pois.mean <- t(X * t(mu %*% t(s)))
  stopifnot(all(dim(pois.mean) == c(B, C)))
  
  # sample Y from poisson and return the dataset
  for (c in 1:C) {
    Y[,c] <- rpois(B, pois.mean[,c])
  }
  return(list(pi0 = pi0, 
              pi = pi,
              s = s,
              X = X, 
              mu = mu,
              Y = Y))
}

gs_step <- function() {
  return() 
}

M <- 4
B <- 1000
C <- 30

dataset <- simulate_data(M, B, C)

n.test <- 10000 # simulations for test purposes
# TEST 1
# estimate s and fix all the rest

a0 <- .01
b0 <- 100

gamma.a <- a0 + apply(dataset$Y, MARGIN = 2, sum) # vector len C
gamma.b <- b0 / (1 + b0 * sum(dataset$mu * dataset$X)) # scalar
s.chain <- matrix(NA, nrow = n.test, ncol = C)
for (i in 1:n.test) {
  s.new <- rgamma(C, gamma.a, scale = gamma.b)
  s.chain[i,] <- s.new
}

# compare if the mean of the samples is close 
# to the true values
# TODO: use confidence intervals to see how close
# quantile(samples, (.025, .975))
# TODO: print plots of the chain
# the samples are to the true values
dataset$s
apply(s.chain, MARGIN = 2, function (x) {
  return(quantile(x, c(.025, .975)))
})

# TODO: finish this!
plts <- lapply(seq_len(C), FUN = function(i) {
  sdf <- as.data.frame(cbind(idx = 1:n.test, val = s.chain[,i]))
  ggplot(data = sdf) +
    geom_line(mapping = aes(x = idx, y = val)) +
    geom_hline(mapping = aes(yintercept = dataset$s[i], color = "red")) +
    ggtitle(cat("Cell ", i, " samples")) +
    xlab("samples") + ylab("val")
})
ggsave("s_plots.pdf", marrangeGrob(grobs = plts, nrow = 4, ncol = 2))
# TEST 2
# estimate pi and fix all the rest

dir.alpha0 <- rep(.1, M)
cn.pair.counts <- matrix(0, nrow = M, ncol = M)
for (b in 2:B) {
  i <- dataset$X[b-1] + 1
  j <- dataset$X[b] + 1
  # count the frequency of the edges/transitions
  cn.pair.counts[i, j] = cn.pair.counts[i, j] + 1
}

dir.alpha.star <- dir.alpha0 + cn.pair.counts

# pi.chain: list of M matrices, each matrix is n.test x M
pi.chain <- lapply(1:M, function(i) {
  pi.new <- rdirichlet(n.test, dir.alpha.star[i,])
  return(pi.new)
})


# compare if the mean of the samples is close 
# to the true values
do.call(rbind,
  lapply(pi.chain, function(pi.mat) {
    return(apply(pi.mat, MARGIN = 2, mean))
  })
)
dataset$pi



# TEST 3
# estimate copy numbers and fix all the rest
# Note: here we sample all the sequence together
#   with forward-backward Gibbs sampling algorithm (FBG)

pi0 <- rep(1/M, M)
# compute the forward probabilities to the end recursively
# initialize the forward prob to 1/M * b(y1)
log.fwd <- sapply(0:(M-1), function (x) {
  # compute the mean of the poisson for each copy number value
  # given s (for each cell) and the fixed mu value
  pois.mean.b <- dataset$s * x * dataset$mu[1]
  log.fwd.el <- logSumExp(dpois(dataset$Y[1,], lambda = pois.mean.b, log = T)) * pi0[x + 1]
  return(log.fwd.el)
})

log.fwd

for (b in 2:B) {
  log.fwd <- logSumExp(dpois(dataset$Y[b,], lambda = dataset$s * dataset$mu[b]))
}

# FIXME: fwd.filtered 
## fwd.filtered <- log.fwd - logSumExp(log.fwd)
