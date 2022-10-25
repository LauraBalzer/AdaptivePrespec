# Function to call to generate synthetic data 
generate.data.wrapper <- function(n, effect, sim, stratify=F, expt_type, verbose=F){
  
  if(sim == "contY"){
    full.data <- gen.data.contY(N = n, expt_type = expt_type, effect = effect, verbose = verbose)
  } else if(sim=='binY'){
      full.data <- gen.data.binY(N=n, expt_type = expt_type, effect=effect, verbose=verbose)
  } 
  
  if(stratify){
    # stratified randomization probably a more clever way to do this
    # if binary outcome, then need to discretize into a binary W1
    strata.var <- full.data$W1 > 0
    A <- rep(NA,n)
    n.W1 <- sum(strata.var)
    A[strata.var] <- sample(c(rep(1,ceiling(n.W1/2)), rep(0,ceiling(n.W1/2))))[1:n.W1]
    n.W0 <- sum(!strata.var)
    A[!strata.var] <- sample(c(rep(1,ceiling(n.W0/2)), rep(0,ceiling(n.W0/2))))[1:n.W0]
    full.data$A <- A

    if(verbose) print( table(A))
    if(verbose) print( table(A, strata.var))
  } else{
    full.data$A <-sample(c(rep(1,n/2), rep(0,n/2)))
  }
  full.data$Y <- ifelse(full.data$A == 1, full.data$Y1, full.data$Y0)
  full.data$alpha <- 1
  full.data
}

# Function to generate data with continuous Y
gen.data.contY <- function(N, expt_type, effect, verbose=F){
  # Keep the covariates constant, and change the synthetic function 
  W1 <- rbinom(N,size=1, prob=0.5) 
  W2 <- rbinom(N, size=1, prob=plogis(0.2*W1)) 
  W3 <- runif(N, min=0, max=5) 
  W4 <- plogis(-2 +W1 + W2+ runif(N, 0, 2)) 
  W5 <- 1+ rbinom(N, size=3, p=0.3) 
  # exogenous noise
  UY <- runif(N, min=0, max=.5)
  # random noise 
  noise <- runif(n, min = 0, max = 1)
  
  
  if(expt_type == "linear_1"){
    get.contY <- function(A, W1,W2,W3,W4,W5, UY, noise){
      90 + A + 4*W1 + 3*W2 + 2*W3 + 3*W4 + 4*W5 + 2*A*W1
      # 90 + 1*A + 4*(W1+W2+W3+W4+W5+W2*W5 ) - W1*W3 + W1*(20-15*W4) + (1-W1)*(-20+ 15*W4) + A*UY/3 
    }
  } else if(expt_type == "linear_2"){
    get.contY <- function(A, W1,W2,W3,W4,W5, UY, noise){
      90 + A + 4*W1 + 3*W2 + 2*W3 + 3*W4 + 4*W5 + 2*A*W1 + A*W2*W1
    }
  } else if(expt_type == "linear_3"){
    get.contY <- function(A, W1,W2,W3,W4,W5, UY, noise){
      90 + A + 4*W1 + 3*W2 + 2*W3 + 3*W4 + 4*W5 + 2*A*W1 + A*W2*W4 + 3*A*W3
    }
  } else if(expt_type == "multicollinear"){
    get.contY <- function(A, W1,W2,W3,W4,W5,UY, noise){
      90 + A + 4*W1 + 3*W2 + 2*W3 + 3*W4 + 4*W5 + 2*A*W1 + A*W2*W4 + 3*A*W3 + noise*W5*A
    }
  } else if(expt_type == "squared"){
    get.contY <- function(A, W1,W2,W3,W4,W5, UY, noise){
      90 + 1*A + 4*(W1+W2+W3+W4+W5) - W1*W3 + W1*(20-15*W4) + (1-W1)*(-20+ 15*W4) 
    }
  } else if(expt_type == "polynomial"){
    get.contY <- function(A, W1,W2,W3,W4,W5, UY, noise){
      90 + 1*A + 4*(W1+W2+W3+W4+W5) - W1*W3 + W1*(20-15*W4)*W3 + (1-W1)*(-20+ 15*W4) 
    }
  } else if(expt_type == "noisy_linear_1"){
    get.contY <- function(A, W1,W2,W3,W4,W5, UY, noise){
      # 90 + A + 4*W1 + 3*W2 + 2*W3 + 3*W4 + 4*W5 + 2*A*W1 + UY
      90 + 0.5*A + 4*W1 + 3*W2 + 2*W3 + 3*W4 + 4*W5 + 0.75*A*W1 + UY
    } 
    }	else if(expt_type == "noisy_linear_1_cand1"){
    get.contY <- function(A, W1,W2,W3,W4,W5, UY, noise){
      180 + 20*A + W1 - 2*W2 + 3*W3 - 5*W4 + 5*A*W1 + 2*UY # Candidate 1
      # 90 + 20*A + W1 - 2*W2 + 3*W3 - 5*W4 + 5*A*W1 + 2*UY # Candidate 2
      # A + 4*W1 + W2 + 5*W3 + 10*A*W4 + 4*W5 + UY # Candidate 3
    }
  } else if(expt_type == "noisy_linear_1_cand2"){
    get.contY <- function(A, W1,W2,W3,W4,W5, UY, noise){
      #180 + 20*A + W1 - 2*W2 + 3*W3 - 5*W4 + 5*A*W1 + 2*UY # Candidate 1
       90 + 20*A + W1 - 2*W2 + 3*W3 - 5*W4 + 5*A*W1 + 2*UY # Candidate 2
      # A + 4*W1 + W2 + 5*W3 + 10*A*W4 + 4*W5 + UY # Candidate 3
    }
  } else if(expt_type == "noisy_linear_1_cand3"){
    get.contY <- function(A, W1,W2,W3,W4,W5, UY, noise){
      #180 + 20*A + W1 - 2*W2 + 3*W3 - 5*W4 + 5*A*W1 + 2*UY # Candidate 1
      # 90 + 20*A + W1 - 2*W2 + 3*W3 - 5*W4 + 5*A*W1 + 2*UY # Candidate 2
       A + 4*W1 + W2 + 5*W3 + 10*A*W4 + 4*W5 + UY # Candidate 3
    }
  } else if(expt_type == "noisy_linear_2"){
    get.contY <- function(A, W1,W2,W3,W4,W5, UY, noise){
      90 + A + 4*W1 + 3*W2 + 2*W3 + 3*W4 + 4*W5 + 2*A*W1 + A*W2*W1 + UY
    }
  } else if(expt_type == "noisy_linear_3"){
    get.contY <- function(A, W1,W2,W3,W4,W5, UY, noise){
      90 + A + 4*W1 + 3*W2 + 2*W3 + 3*W4 + 4*W5 + 2*A*W1 + A*W2*W4 + 3*A*W3 + 2*UY 
    }
  } else if(expt_type == "noisy_multicollinear_cand1"){
    get.contY <- function(A, W1,W2,W3,W4,W5, UY, noise){
      #180 + .2*A + W1 - 1.2*W2 + W3 - 1.3*W4 + W5 + .05*A*W1 + .05*A*W2 + .05*A*W3 + .05*A*W4 + noise*W5*A + .5*A*UY/3 # Candidate 1
      150 + .1*A + W1 - W2 + 4*W3 - W4 + W5 + .1*A*W1 + .1*A*W3 + A*UY/3 # Candidate 2
    }
  } else if(expt_type == "noisy_multicollinear_cand2"){
    get.contY <- function(A, W1,W2,W3,W4,W5, UY, noise){
      #180 + A + W1 - W2 + W3 - W4 + W5 + 5*A*W1 + A*W2 + 7*A*W3 + A*W4 + noise*W5*A + A*UY/3 # Candidate 1
      150 + .3*A + W1 - W2 + 4*W3 - W4 + W5 + .1*A*W1 + .2*A*W3 + A*UY/3 # Candidate 2
    }
  } else if(expt_type == "noisy_squared"){
    get.contY <- function(A, W1,W2,W3,W4,W5, UY, noise){
      90 + 1*A + 4*(W1+W2+W3+W4+W5) - W1*W3 + W1*(20-15*W4) + (1-W1)*(-20+ 15*W4) + UY 
    }
  } else if(expt_type == "noisy_polynomial"){
    get.contY <- function(A, W1,W2,W3,W4,W5, UY, noise){
      90 + 1*A + 4*(W1+W2+W3+W4+W5) - W1*W3 + W1*(20-15*W4)*W3 + (1-W1)*(-20+ 15*W4) + UY 
    }
  }
  
  # Get the function for continuous Y
  Y0 <- get.contY(A=0, W1=W1, W2=W2, W3=W3, W4=W4, W5=W5, UY=UY, noise=noise)
  if(effect){
    Y1 <- get.contY(A=1, W1=W1, W2=W2, W3=W3, W4=W4, W5=W5, UY=UY, noise=noise)
  } else{
    Y1 <- Y0 #if under the null, then set the counterfactual outcome Y1 to Y0
  }
  if(verbose) print(round(c(mean(Y1), mean(Y0)),2))
  
  data.frame(cbind(id=1:n, U=1, W1, W2, W3, W4, W5, Y1, Y0) )
}

# Function to generate data with binary Y
gen.data.binY <- function(N, expt_type, effect, verbose=F){
  # Keep the covariates constant, and change the synthetic function 
  W1 <- rnorm(N, 0, 1)
  W2 <- rnorm(N, 0, 1)
  W3 <- rnorm(N, 0, 1)
  W4 <- runif(N, 0, 1)
  W5 <- runif(N, 0, 1)
  UY <- runif(N, 0, 1) # exogenous noise
  noise <- runif(N, 0, 1) # noisy variable
  
  if(expt_type == "linear_1"){
    get.binY <- function(A, W1,W2,W3,W4,W5,UY, noise,verbose = F){
      p <- plogis(A + 4*W1 + 3*W2 + 2*W3 + 3*W4 + 4*W5 + 2*A*W1)
      if(verbose) { print(summary(p)); hist(p) }
      p <- as.numeric(UY <  p)
      p
      # 90 + 1*A + 4*(W1+W2+W3+W4+W5+W2*W5 ) - W1*W3 + W1*(20-15*W4) + (1-W1)*(-20+ 15*W4) + A*UY/3 
    }
  } else if(expt_type == "linear_2"){
    get.binY <- function(A, W1,W2,W3,W4,W5, UY, noise,verbose = F){
      p <- plogis(A + 4*W1 + 3*W2 + 2*W3 + 3*W4 + 4*W5 + 2*A*W1 + A*W2*W1)
      if(verbose) { print(summary(p)); hist(p) }
      p <- as.numeric(UY <  p)
      p
    }
  } else if(expt_type == "linear_3"){
    get.binY <- function(A, W1,W2,W3,W4,W5, UY, noise,verbose = F){
      p <- plogis(A + 4*W1 + 3*W2 + 2*W3 + 3*W4 + 4*W5 + 2*A*W1 + A*W2*W4 + 3*A*W3)
      if(verbose) { print(summary(p)); hist(p) }
      p <- as.numeric(UY <  p)
      p
    }
  } else if(expt_type == "multicollinear"){
    get.binY <- function(A, W1,W2,W3,W4,W5,UY, noise, verbose = F){
      p <- plogis(A + 4*W1 + 3*W2 + 2*W3 + 3*W4 + 4*W5 + 2*A*W1 + A*W2*W4 + 3*A*W3 + noise*W5*A)
      if(verbose) { print(summary(p)); hist(p) }
      p <- as.numeric(UY <  p)
      p
    }
  } else if(expt_type == "squared"){
    get.binY <- function(A, W1,W2,W3,W4,W5, UY, noise,verbose = F){
      p <- plogis(A + 4*(W1+W2+W3+W4+W5) - W1*W3 + W1*(20-15*W4) + (1-W1)*(-20+ 15*W4))
      if(verbose) { print(summary(p)); hist(p) }
      p <- as.numeric(UY <  p)
      p
    }
  } else if(expt_type == "polynomial"){
    get.binY <- function(A, W1,W2,W3,W4,W5, UY, noise,verbose = F){
      p <- plogis(1*A + 4*(W1+W2+W3+W4+W5) - W1*W3 + W1*(20-15*W4)*W3 + (1-W1)*(-20+ 15*W4))
      if(verbose) { print(summary(p)); hist(p) }
      p <- as.numeric(UY <  p)
      p
    }
  } else if(expt_type == "noisy_linear"){
    get.binY <- function(A, W1, W2, W3, W4, W5, UY, noise, verbose = F){
      p <- plogis(A + W1 - W2 + W3 - W4 + W5 - 2*A*W1 + UY)
      if(verbose) { print(summary(p)); hist(p) }
      p <- as.numeric(UY <  p)
      p
    }
  } else if(expt_type == "noisy_linear_2"){
    get.binY <- function(A, W1,W2,W3,W4,W5, UY, noise,verbose = F){
      p <- plogis(A + 4*W1 + 3*W2 + 2*W3 + 3*W4 + 4*W5 + 2*A*W1 + A*W2*W1 + UY)
      if(verbose) { print(summary(p)); hist(p) }
      p <- as.numeric(UY <  p)
      p
    }
  } else if(expt_type == "noisy_linear_3"){
    get.binY <- function(A, W1,W2,W3,W4,W5, UY, noise,verbose = F){
      p <- plogis(A + 4*W1 + 3*W2 + 2*W3 + 3*W4 + 4*W5 + 2*A*W1 + A*W2*W4 + 3*A*W3 + 2*UY)
      if(verbose) { print(summary(p)); hist(p) }
      p <- as.numeric(UY <  p)
      p
    }
  } else if(expt_type == "noisy_multicollinear"){
    get.binY <- function(A, W1,W2,W3,W4,W5,UY, noise, verbose = F){
      # p <- plogis(A + 4*W1 + 3*W2 + 2*W3 + 3*W4 + 4*W5 + 2*A*W1 + A*W2*W4 + 3*A*W3 + noise*W5*A + A*UY/3)
      p <- plogis(A + W1 + W2 + W3 + W4 + W5 + A*W1 + A*W2*W4 + A*W3 + noise*W5*A + A*UY/3)
      if(verbose) { print(summary(p)); hist(p) }
      p <- as.numeric(UY <  p)
      p
    }
  } else if(expt_type == "noisy_squared"){
    get.binY <- function(A, W1,W2,W3,W4,W5, UY, noise,verbose = F){
      p <- plogis(A + 4*(W1+W2+W3+W4+W5) - W1*W3 + W1*(20-15*W4) + (1-W1)*(-20+ 15*W4) + UY)
      if(verbose) { print(summary(p)); hist(p) }
      p <- as.numeric(UY <  p)
      p
    }
  } else if(expt_type == "noisy_polynomial"){
    get.binY <- function(A, W1,W2,W3,W4,W5, UY, noise,verbose = F){
      p <- plogis(A + (W1+W2+W3+W4+W5) - W1*W3 + W1*(2*W4)*W3 + (1-W1)*(-W4) + UY)
      if(verbose) { print(summary(p)); hist(p) }
      p <- as.numeric(UY <  p)
      p
    }
  }
  
  # Get the function for binary Y
  Y0 <- get.binY(A=0, W1=W1, W2=W2, W3=W3, W4=W4, W5=W5, UY=UY, noise=noise)
  if(effect){
    Y1 <- get.binY(A=1, W1=W1, W2=W2, W3=W3, W4=W4, W5=W5, UY=UY, noise=noise)
  } else{
    Y1 <- Y0 #if under the null, then set the counterfactual outcome Y1 to Y0
  }
  if(verbose){
    print(round(c(mean(Y1), mean(Y0)),2))
  }
  if(unique(is.infinite(Y1)) != FALSE){ stop("Something wrong")}
  if(unique(is.infinite(Y0)) != FALSE){ stop("Something wrong")}
  if(unique(is.nan(Y1)) != FALSE){ stop("Something wrong")}
  if(unique(is.nan(Y0)) != FALSE){ stop("Something wrong")}
  if(unique(is.na(Y1)) != FALSE){ stop("Something wrong")}
  if(unique(is.na(Y0)) != FALSE){ stop("Something wrong")}
  
  data.frame(cbind(id=1:N, U=1, W1, W2, W3, W4, W5, Y1, Y0) )
}
