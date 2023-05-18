N <- 5000

get.rsq <- function(N){
  X<- NULL
  
  # Keep the covariates constant, and change the synthetic function 
  W1 <- rbinom(N,size=1, prob=0.5) 
  W2 <- rbinom(N, size=1, prob=plogis(0.2*W1)) 
  W3 <- runif(N, min=0, max=5) 
  W4 <- plogis(-2 +W1 + W2+ runif(N, 0, 2)) 
  W5 <- 1+ rbinom(N, size=3, p=0.3) 
  # exogenous noise
  UY <- runif(N, min=0, max=.5)
  # random noise 
  noise <- runif(N, min = 0, max = 1)
  
  A <- rbinom(N, 1,.5)
  
  #"noisy_only_predictor_2": 
  Y <- 90 + .1*A + 3*UY + noise
  g<- glm(Y~A )
  X <- c(X, with(summary(g), 1 - deviance/null.deviance) )
  #"noisy_linear_1_r_less":
  Y <- 90 + 0.07*A + .7*W1 + .3*W2 + .1*W3 + .3*W4 + .4*W5 + 0.25*A*W1 + 5*UY + noise
  g <- glm(Y~A + W1+W2+W3+W4+W5+A*W1)
  X <- c(X, with(summary(g), 1 - deviance/null.deviance) )
  
  #"noisy_multicollinear_cand1_r_less":
  Y <- 90 + .17*A + .33*(W1+W2+W3+W4+W5) - .2*W1*W3 + .5*W1*(.8-.6*W4)*W3 + 
    .25*(1-W1)*(-.2+ .15*W4) + 4.7*UY + noise
  g <- glm(Y~A + W1+W2+W3+W4+W5+A*W1+ A*W3)
  X <- c(X, with(summary(g), 1 - deviance/null.deviance) )
  
  #"noisy_polynomial_r_less": 
  Y <- 90 + .17*A + .33*(W1+W2+W3+W4+W5)-.2*W1*W3 +.5*W1*(.8-.6*W4)*W3 + .25*(1-W1)*(-.2+ .15*W4) + 
    4.7*UY + noise
  g<- glm(Y~A + W1+W2+W3+W4+W5+ W1*W3 +W1*W4*W3 + W1*W4)
  X<- c(X, with(summary(g), 1 - deviance/null.deviance) )
  X
}

set.seed(1)
Rsq <- NULL
nReps <- 5000
for(j in 1:nReps){
  Rsq <- cbind(Rsq, get.rsq(N))
}


