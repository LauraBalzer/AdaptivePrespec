rm(list=ls())


get.rsq.cont <- function(N){
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
  g <- glm(Y~A + W1+W2+W3+W4+W5+ A:W1)
  X <- c(X, with(summary(g), 1 - deviance/null.deviance) )
  
  #"noisy_multicollinear_cand1_r_less":
  Y <- 150 + .05*A + .33*W1 - .25*W2 + .5*W3 - .2*W4 + .05*W5 + .01*A*W1 + .02*A*W3 + .3*A*UY + 5.8*UY + noise
  g <- glm(Y~A + W1+W2+W3+W4+W5+ A:W1+A:W3)
  X <- c(X, with(summary(g), 1 - deviance/null.deviance) )
  
  #"noisy_polynomial_r_less": 
  Y <- 90 + .17*A + .33*(W1+W2+W3+W4+W5) - .2*W1*W3 + .5*W1*(.8-.6*W4)*W3 + .25*(1-W1)*(-.2+ .15*W4) + 4.7*UY + noise
  g<- glm(Y~A + W1+W2+W3+W4+W5+ W1:W3 +W1:W4:W3 + W1:W4)
  X<- c(X, with(summary(g), 1 - deviance/null.deviance) )
  X
}


library(rcompanion)
get.rsq.bin <- function(N){
  X<- NULL
  # Keep the covariates constant, and change the synthetic function 
  W1 <- rnorm(N, 0, 1)
  W2 <- rnorm(N, 0, 1)
  W3 <- rnorm(N, 0, 1)
  W4 <- runif(N, 0, 1)
  W5 <- runif(N, 0, 1)
  UY <- runif(N, 0, 1) # exogenous noise
  noise <- runif(N, 0, 1) # noisy variable
  A <- rbinom(N, 1,.5)
  
  #Treatment only 
  p <- plogis(0.1 * A + UY + noise) 
  Y<- as.numeric(UY< p)
  g<- glm(Y~A, family='binomial')
  X <- c(X, nagelkerke(g)$Pseudo.R.squared.for.model.vs.null[1] )
  #Noisy Linear 
  p <- plogis(A + W1 - W2 + W3 - W4 + W5 - 2*A*W1 + noise)
  Y<- as.numeric(UY< p)
  g<- glm(Y~A+W1+W2+W3+W4+W5+A:W1, family='binomial')
  X <- c(X,nagelkerke(g)$Pseudo.R.squared.for.model.vs.null[1])
  
  #Noisy multicollinear 
  p <- plogis(A + W1 + W2 + W3 + W4 + W5 + A*W1 + A*W2*W4 + A*W3 + noise*W5*A + noise)
  Y<- as.numeric(UY< p) 
  g<- glm(Y~A + W1 + W2 + W3 + W4 + W5 + A:W1 + A:W2:W4 + A:W3 + W5:A, family='binomial')
  X <- c(X, nagelkerke(g)$Pseudo.R.squared.for.model.vs.null[1] )

  # Noisy Polynomial 
  p <- plogis(A + (W1+W2+W3+W4+W5) - W1*W3 + W1*(2*W4)*W3 + (1-W1)*(-W4) + noise)
  Y<- as.numeric(UY< p ) 
  g<- glm(Y ~ A + W1+W2+W3+W4+W5 - W1:W3 + W1:W4:W3+W1:W4, family='binomial')
  X <- c(X,nagelkerke(g)$Pseudo.R.squared.for.model.vs.null[1] )
  X

}
set.seed(1)
Rsq.cont <- Rsq.bin <-  NULL
N <- 500
nReps <- 5000
for(j in 1:nReps){
  Rsq.cont<- rbind(Rsq.cont, get.rsq.cont(N))
  Rsq.bin <- rbind(Rsq.bin, get.rsq.bin(N))
}

library(xtable)
x <- rbind( colMeans(Rsq.bin), colMeans(Rsq.cont))
rownames(x) <- c('Binary', 'Continuous')
colnames(x) <- c('Treatment only', 'Linear', 'Interactive', 'Polynomial')
xtable(x, digits=c(1,2,2,2,2))
