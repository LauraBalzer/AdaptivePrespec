##############
# TMLE_Functions_Meta.R 
# R code to obtain point estimates with TMLE

# Modified from Stage2_Functions.R in https://github.com/LauraBalzer/TwoStageTMLE

#-----------------------------------------------------#-----------------------------------------------------
# do.TMLE: master function to run TMLE for point estimation 
#
# input: 
#		goal - aRR: arithmetic risk ratio; O/W risk difference,
#   target of inference: cluster-level ("clust") or pooled-indv effect ("indv") (target) 
#		training data (train),
#		prespecified adjustment variables for the conditional mean outcome (QAdj), 
#		prespecified adjustment variables for the propensity score (gAdj),
#		initial estimator of the conditional mean outcome (Q.out),
#		estimator of the propensity score (p.out),
#		indicator to print updates (verbose)
#
# output: list 
#		training data augmented with estimates,
#		prespecified adjustment variables for the conditional mean outcome (QAdj), 
#		prespecified adjustment variables for the propensity score (gAdj),
#		initial estimator of the conditional mean outcome (Q.out),
#		estimator of the propensity score (p.out),
#		estimated fluctuation coef (epsilon),
#-----------------------------------------------------#-----------------------------------------------------

# UPDATES
#   expand the candidate prediction functions with Qform (outcome regression) and gform (pscore)
#     (1) 'glm' for main terms; 
#     (2) 'step' for stepwise;
#     (3) 'step.interaction for stepwise interaction
#     (4) 'lasso' for LASSO; 
#     (5) 'mars' for multivariate adaptive regression splines (MARS)
#     (6) 'mars.corP' for MARS after screening 
# See do.Init.Qbar for details

do.TMLE <- function(goal, target='NOT.WORKING', sample.effect, train,
                    QAdj, Qform='glm', gAdj=NULL, gform='glm',
                    Q.out=NULL, p.out=NULL, 
                    scale_value, scale_value_min,
                    doing.CV=F, verbose=F) {	
  
  #=====================================================
  # Step1 - initial estimation of E(Y|A,W)= Qbar(A,W)
  #=====================================================
  
  # run glm on the adjustment set
  Q <- do.Init.Qbar(train=train, QAdj=QAdj, Qform=Qform, glm.out=Q.out, verbose=verbose)
  train <- Q$train
  
  #==========================================================
  # Step2: Calculate the clever covariate
  #==========================================================	
  
  G <- get.clever.cov(train=train, gAdj=gAdj, gform=gform, p.out=p.out, verbose=verbose)
  train <- G$train
  
  #==========================================================
  # Step3: Targeting
  #==========================================================
  
  eps <- get.epsilon(train=train, goal=goal, verbose=verbose)
  
  train <- do.targeting(train=train, eps=eps, goal=goal)
  
  #==========================================================
  # Step4: Parameter  estimation
  # will unscale if appropriate
  #==========================================================
  
  if(nrow(train)> length(unique(train$id)) & target=='clust')	{
    # IF DATA ARE AT THE INDV-LEVEL, BUT GOAL IS THE CLUSTER-LEVEL EFFECT 
    # get point estimates by aggregating to the cluster level 
    #   (e.g. by taking the weighted sum)
    # then take mean of cluster-level endpoints
    if(!doing.CV) print('data=indv; target=clust')
    R1<- mean( aggregate(data.frame(train$alpha*train$Qbar1W.star), by=list(train$id), sum)[,2] )
    R0<- mean( aggregate(data.frame(train$alpha*train$Qbar0W.star), by=list(train$id), sum)[,2] )
  } else{
    # OTHERWISE, JUST TAKE THE WEIGHTED MEAN ACROSS ALL ROWS
    # future work: robustify so that dont need weights that sum to J
    R1 <- mean( train$alpha*train$Qbar1W.star )
    R0 <- mean( train$alpha*train$Qbar0W.star ) 
  }
  
  # UNSCALE THE OUTCOME 
  R1 <- R1*(scale_value - scale_value_min) + scale_value_min
  R0 <- R0*(scale_value - scale_value_min) + scale_value_min
  
  #==========================================================
  # Step 5: Variance estimation
  #==========================================================
  variance.out <- get.IC.variance(goal=goal, target=target, Vdata=train, R1=R1, R0=R0,
                         scale_value = scale_value, scale_value_min = scale_value_min, 
                         doing.CV = doing.CV, sample.effect = sample.effect)
  

  RETURN<- list(train=train, 	
                QAdj=Q$QAdj, Qform=Q$Qform, Q.out=Q$glm.out,
                gAdj=G$gAdj, gform=G$gform, p.out=G$p.out, 
                eps=eps, R1=R1, R0=R0, 
                var.R1=variance.out$var.R1, 
                var.R0=variance.out$var.R0,
                var.pair=variance.out$var.pair, 
                var.break=variance.out$var.break)	
  RETURN
}

#-----------------------------------------------------#-----------------------------------------------------
# do.Init.Qbar - function to do initial estimation of E[Y|A,W] = Qbar(A,W)
# 	input: data set, adjustment variable(s), outcome regression fit, verbose
# 	output:	adjustment variable(s),	outcome regression fit 
#		  data set augmented w/ initial predictions: Qbar(A,W), Qbar(1,W) and Qbar(0,W)
#-----------------------------------------------------#-----------------------------------------------------

# UDPATE - expand the candidate functions for the outcomes regression Qform, 
do.Init.Qbar<- function(train, QAdj, Qform='glm', glm.out=NULL, verbose=F){
  
  if( is.null(QAdj) ){
    QAdj<- 'U'
  }
  
  train.temp <- train[, c(QAdj, 'A', 'Y')]
  
  X1<- X0<- train.temp
  X1$A<-1; X0$A<- 0	
  # needed for penalized regression
  Xl <- model.matrix(~-1 +. , subset(train.temp, select=-Y) )
  X1l <- model.matrix(~-1 +. ,  subset(X1, select=-Y))
  X0l <- model.matrix(~-1 +. ,  subset(X0, select=-Y))
  
  if( is.null(glm.out) ){
    # fit using the training data
    # run main terms regression
    glm.out<- suppressWarnings( glm( Y~. , family='binomial', data=train.temp, weights=train$alpha ) )	
    
    if (Qform=='step'){ # stepwise 
      glm.out <- step(glm.out, direction = "both", trace = 0, k = 2)
    } else if (Qform=='step.interaction'){
      glm.out <- step(glm.out, scope=Y~.^2, direction = "both", trace = 0, k = 2)
    } else if (Qform=='lasso'){
      familyl <- ifelse(sum(unique(train$Y))>2, 'gaussian', 'binomial')
      glm.out <- glmnet(x=Xl,  y=train$Y, weights=train$alpha,
                        family=familyl, alpha=1, nlambda = 100)
    } else if(Qform %in% c('mars', 'mars.corP')){
      # using default settings of SL.earth in SuperLearner 
      X <- subset(train.temp, select=-Y)
      if(Qform=='mars.corP'){
        X <- X[,screen.corP(Y=train$Y, X=X, family='binomial')==T]
      } 
      
      glm.out <- earth(x = X,  y=train$Y, weights=train$alpha,
                       degree = 2, penalty = 3, 
                       nk = max(21, 2 * ncol(X) + 1), pmethod = "backward", nfold = 0, 
                       ncross = 1, minspan = 0, endspan = 0,
                       glm = list(family = binomial))
    } 
    if(verbose) print(glm.out)
  }	

  # get initial predictions
  if(Qform=='lasso'){
    # from LASSO 
    QbarAW <- predict(glm.out, newx=Xl, type='response', s = min(glm.out$lambda) )
    Qbar1W <- predict(glm.out, newx=X1l, type='response', s = min(glm.out$lambda) )
    Qbar0W <- predict(glm.out, newx=X0l, type='response', s = min(glm.out$lambda) )
  }else{
    # for glms 
    QbarAW <- predict(glm.out, newdata=train.temp, type='response')
    Qbar1W <- predict(glm.out, newdata=X1, type='response')
    Qbar0W <- predict(glm.out, newdata=X0, type='response')
  }
  Qbar <- data.frame(QbarAW, Qbar1W, Qbar0W)
  colnames(Qbar) <- c('QbarAW', 'Qbar1W', 'Qbar0W')
  list(QAdj=QAdj, Qform=Qform, glm.out=glm.out, train=cbind(train, Qbar) )
}

#-----------------------------------------------------#-----------------------------------------------------
# get.clever.cov - function calculate the clever covariate
# 	input: 
#		  data set, adjustment variable(s), pscore regression, verbose
# 	output: 
#		  adjustment variable(s), pscore regression, 
#		  data set augmented with pscore & clever covariate (H.AW, H.1W, H.0W)
#-----------------------------------------------------#-----------------------------------------------------
get.clever.cov<- function(train, gAdj, gform, p.out=NULL, verbose=F){
  
  if( is.null(gAdj) ){
    gAdj <- 'U'
  }
  
  train.temp <- train[, c(gAdj, 'A')]  
  # needed for penalized regression
  Xl <- model.matrix(~-1 +. , subset(train.temp, select=-A) )
  
  if( is.null(p.out) ){ # fit pscore on training set 	
    # run main terms regression
    p.out<-   suppressWarnings( glm( A~. , family='binomial', data= train.temp, weights=train$alpha) )
    
    if (gform=='step'){ # stepwise 
      p.out <- step(p.out, direction = "both", trace = 0, k = 2)
    } else if (gform=='step.interaction'){
      p.out <- step(p.out, scope=A~.^2, direction = "both", trace = 0, k = 2)
    } else if (gform=='lasso'){
      p.out <- glmnet(x=Xl,  y=train$A, weights=train$alpha,
                        family='binomial', alpha=1, nlambda = 100)
    } else if(gform %in% c('mars', 'mars.corP')){
      # using default settings of SL.earth in SuperLearner 
      X <- subset(train.temp, select=-A)
      if(gform=='mars.corP'){
        X <- X[,screen.corP(Y=train$A, X=X, family='binomial')==T]
      } 
      p.out <- earth(x = X,  y=train$A, weights=train$alpha,
                       degree = 2, penalty = 3, 
                       nk = max(21, 2 * ncol(X) + 1), pmethod = "backward", nfold = 0, 
                       ncross = 1, minspan = 0, endspan = 0,
                       glm = list(family = binomial))
    } 
    if(verbose){ print(p.out)}
    
  }
  
  # now use p.out to get estimated pscores
  if(gform!='lasso'){
    pscore <- predict(p.out, newdata= train.temp,  type="response")
  }else{
    pscore <- predict(p.out, newx=Xl, type='response', s = min(p.out$lambda) )
  }
  
  # bound g - should not apply for a randomized trial
  pscore [pscore < 0.025] <- 0.025
  pscore [pscore > 0.975] <- 0.975
  
  A.train <- train$A
  # Clever covariate is two-dimensional; 
  H.1W <- A.train/pscore 
  H.0W <- (1-A.train)/(1-pscore )
  # via delta method
  H.AW <- H.1W - H.0W
  
  p.df <- data.frame(pscore, H.1W , H.0W , H.AW)
  colnames(p.df) <- c('pscore', 'H.1W' , 'H.0W' , 'H.AW')
  list(gAdj=gAdj, gform=gform, p.out=p.out,  train=data.frame(train, p.df) ) 
}	


#-----------------------------------------------------#-----------------------------------------------------
# get.epsilon - function calculate the fluctuation coefficient
# 	input: 
#		  data set, goal with 'aRR'=arithmetic RR, verbose
# 	output: 
#	  	estimated fluctuation coefficient (eps)
#-----------------------------------------------------#-----------------------------------------------------
get.epsilon <- function(train, goal, verbose=F){
  
  A.train<- train$A
  Y.train<- train$Y
  
  # Skip fitting if outcome=0 or outcome=1 
  # for all observations in either txt or control group
  Skip.update <-  mean(Y.train[A.train==1])==0 | mean(Y.train[A.train==0])==0 |  
    mean(Y.train[A.train==1])==1 | mean(Y.train[A.train==0])==1 
  
  if(goal=='RD'){
    
    # if going after RD, then use a 1-dim clever covariate
    if(!Skip.update){
      logitUpdate<- suppressWarnings( 
        glm(Y.train ~ -1 +offset(qlogis(train$QbarAW )) + train$H.AW, family="binomial",  weights=train$alpha))
      eps<-logitUpdate$coef
    } else{
      eps<- 0
    }
    names(eps) <- 'H.AW'
    
  }else{
    # targeting the risk or odds ratio requires a two-dimensional clever covariate
    
    if( !Skip.update  ){
      logitUpdate<- suppressWarnings(
        glm(Y.train ~ -1 +offset(qlogis(train$QbarAW )) + train$H.0W + train$H.1W, family="binomial", weights=train$alpha))
      eps<-logitUpdate$coef
    } else{
      eps <- c(0,0)
    }
    names(eps)<- c('H.0W', 'H.1W')	
  }
  if(verbose) print(eps)
  
  eps
}

#-----------------------------------------------------#-----------------------------------------------------
# do.targeting - function to update initial estimators of QbarAW
# 	input: 
#		  data set (train), fluctuation coefficient (eps), goal (aRR= arithmetic risk ratio; otherwise RD)
# 	output: data.frame w/ targeted predictions: Qbar*(A,W), Qbar*(1,W), Qbar*(0,W)
#-----------------------------------------------------#-----------------------------------------------------

do.targeting <- function(train, eps, goal){
  
  g1W<- train$pscore
  g0W<- (1 - g1W)
  
  if(goal=='RD'){
    
    # updated QbarAW estimates for training set. 
    QbarAW.star <- plogis( qlogis(train$QbarAW ) + eps*train$H.AW)	
    Qbar1W.star <- plogis( qlogis(train$Qbar1W ) + eps/g1W )
    Qbar0W.star <- plogis( qlogis(train$Qbar0W ) - eps/g0W )
    
  }else{
    # updated QbarAW estimates for training set. 
    QbarAW.star <- plogis( qlogis(train$QbarAW) + eps['H.0W']*train$H.0W + eps['H.1W']*train$H.1W)	
    Qbar0W.star <- plogis( qlogis(train$Qbar0W) + eps['H.0W']/g0W )
    Qbar1W.star <- plogis( qlogis(train$Qbar1W) + eps['H.1W']/g1W )
  }
  train <- data.frame(train, QbarAW.star, Qbar1W.star, Qbar0W.star)		
  train
}

