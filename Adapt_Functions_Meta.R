##############
# Adapt_Functions_Meta.R 
# R code to implement adaptive prespecification

# Modified from Adapt_Functions.R in https://github.com/LauraBalzer/TwoStageTMLE


#-----------------------------------------------------#-----------------------------------------------------
# do.adaptive.prespec: function to implement adaptive prespecification as described in 
#		Balzer et al. "Adaptive pre-specification in randomized trials with and without pair-matching"
#
# input: goal of analysis ('aRR' or 'RD'), 
#   target of inference: cluster-level ("clust") or pooled-indv effect ("indv") (target) 
#		indicator to break the match (break.match)
#		dataset (Ldata)
#		adjustment variables (cand.QAdj) - CHANGE these do not have to be at the cluster-level,
#		selected adjustment variable for the outcome regression (QAdj)
#		selected adjustment variable for the pscore regression (gAdj)
#   maximum value for outcome scaling (scale_value),
#   minimum value for outcome scaling (scale_value_min)
#	output: selection for candidate TMLE
#-----------------------------------------------------#-----------------------------------------------------

# UPDATE - expanded candidates for both the outcome regression & pscore estimation
#   cand.QAdj = candidate adjustment variable sets for the outcome regression (list)
#   cand.Qform = candidate adjustment algorithms for the outcome regression (vector)
#   cand.gAdj = candidate adjustment variable sets for the pscore (list)
#   cand.gform = candidate adjustment algorithms for the pscore (vector)
#
# remove.pscore: if T, remove the variable(s) selected for adjustment in the outcome regression from 
#   candidates for the pscore... should only be used if doing adaptive prespec in RCT with few indpt units

do.adaptive.prespec<- function(goal, target='indv', sample.effect=T, break.match=T, Ldata, V=5,
                               cand.QAdj, cand.Qform, cand.gAdj, cand.gform, remove.pscore=F,
                               QAdj=NULL, gAdj=NULL, scale_value, scale_value_min, verbose=F
                               ){
  
  
  
  # UPDATE: GET THE FOLDS HERE! 
  # get the indpt units (will be observation in indv RCT)
  if( !break.match ){
    Ldata$indpt.unit <- Ldata$pair
  }else{
    Ldata$indpt.unit <- Ldata$id
  }
  unique.ids <- unique(Ldata$indpt.unit)
  
  # get folds
  if(length(unique.ids)>40){
    # V-fold CV 
    folds <- get.folds(V=V, Y=Ldata$Y, ids=unique.ids)
  } else{
    # leave-one-out CV
    folds <- vector("list", length(unique.ids))
    for (v in seq(length(unique.ids))){
      folds[[v]] <- unique.ids[v]
    }
  }
  
  if( is.null(QAdj) ){
    
    # if(verbose) print(cand.Qform)
    
    # do adaptive pre-specification to select from candidate approaches for Qbar
    select.Q <- suppressWarnings( CV.selector(goal=goal, target=target,sample.effect=sample.effect,
                                              break.match=break.match, Ldata=Ldata,
                                              CAND.ADJ=cand.QAdj, CAND.FORM=cand.Qform, forQ=T, 
                                              scale_value=scale_value, scale_value_min=scale_value_min,
                                              folds=folds
                                              ) )
    
    if(verbose) print(select.Q)
    Q.index <- select.Q$adj.index
    QAdj <- select.Q$Adj
    Qform <- select.Q$form
    
    # if select unadjusted estimator for QbarAW=E(Y|A,W), then stop
    if( sum(QAdj == 'U') ){ 
      g.index<- -99; gAdj <- 'U'; gform <- 'glm'
      var.CV <- select.Q$var.CV
    } 
    
    # if did not select the unadjusted for initial estimation of Qbar(A,W), 
    if( sum(QAdj == 'U')==0 & remove.pscore){
      # then need to remove this variable from the candidate for pscore estimation
      # useful for very small trials with simple adaptive prespec
      # THIS ONLY WORKS IF cand.QAdj === cand.gAdj AND cand.gform='glm
      if(sum(unlist(cand.gAdj) != unlist(cand.QAdj))==0 &
         sum(cand.gform !='glm')==0){
         print('removing selected QAdj from candidates for pscore')
        cand.gAdj[[Q.index]] <- NULL		
        cand.gform <- cand.gform[-Q.index]
        
      }
    }
    
  }
  
  if( is.null(gAdj) ){ 		
   # if(verbose) print(cand.gform)
    
    select.G <- suppressWarnings( CV.selector(goal=goal, target=target, sample.effect=sample.effect,
                                              break.match=break.match, Ldata=Ldata,
                                              CAND.ADJ=cand.gAdj, CAND.FORM=cand.gform, forQ=F, 
                                              # input selected variables/form of the outcome regression
                                              QAdj= QAdj, Qform=Qform,
                                              scale_value=scale_value, scale_value_min=scale_value_min,
                                              folds=folds) )
    
    if(verbose) print(select.G)
    g.index <- select.G$adj.index
    gAdj <- select.G$Adj
    gform <- select.G$form
    var.CV <- select.G$var.CV		
  }		
  
  list(Q.index=Q.index, QAdj=QAdj, Qform=Qform, 
       g.index=g.index, gAdj=gAdj, gform=gform, var.CV=var.CV )
}



#-----------------------------------------------------#-----------------------------------------------------
# CV.selector: function to estimate the cross-validated risk
#		Loss function is the squared-IC; Risk is then the variance of the TMLE
#		See Balzer et al. "Adaptive pre-specification in randomized trials with and without pair-matching"
#
#	input: goal of analysis ('aRR' or 'RD), 
#   target of inference 
#		indicator to break the match (break.match)
#		dataset (Ldata)
#		candidate adjustment variables; they do not have to be at the cluster-level
#		indicator if for the conditional mean outcome (forQ)
#		selected adjustment variable for the outcome regression (QAdj)
#	output: selection for adjustment variable (corresponding to a TMLE)
#-----------------------------------------------------#-----------------------------------------------------


CV.selector <- function(goal, target, sample.effect, break.match, 
                        Ldata, CAND.ADJ, CAND.FORM, 
                        forQ, QAdj=NULL, Qform=NULL,
                        scale_value, scale_value_min, folds){
  
  if( length(CAND.FORM)==1 ){
    # if exploring only one estimation algorithm (usually GLM)
    # then need to replicate the number forms
    CAND.FORM <- rep(CAND.FORM, length(CAND.ADJ))
  }
  
  if( length(CAND.FORM) != length(CAND.ADJ)){
    print('******* PROBLEM- MISMATCH SIZE OF ADJ VAR AND QFORM')
  }
  
  
  # Number of candidate estimators is given by length Qform//gform
  num.tmles <- length(CAND.FORM)
  CV.risk <-  var.CV <-  rep(NA, num.tmles)

  for(k in 1: num.tmles){	
    
    if(forQ){
      # if selecting the adjustment approach for the outcome regression
      IC.temp<- get.IC.CV(goal=goal, target=target, break.match=break.match, Ldata=Ldata,
                          QAdj=CAND.ADJ[[k]], Qform=CAND.FORM[k], gAdj=NULL, gform='glm',
                          scale_value=scale_value, scale_value_min=scale_value_min, 
                          folds=folds, sample.effect=sample.effect)
    } else{
      # if collaboratively selecting the adjustment approach for the pscore
      IC.temp<- get.IC.CV(goal=goal, target=target, break.match=break.match, Ldata=Ldata, 
                          QAdj=QAdj, Qform=Qform, 
                          gAdj= CAND.ADJ[[k]], gform=CAND.FORM[k],
                          scale_value=scale_value, scale_value_min=scale_value_min, 
                          folds=folds, sample.effect=sample.effect)
    }
    
    # estimating the CV risk for each candidate
    CV.risk[k]<- IC.temp$CV.risk
    # estimating the CV variance for that TMLE
    var.CV[k] <- IC.temp$var.CV
    
  }
  # select the candidate estimator resulting in the smallest CV-risk
  adj.index<- which.min(CV.risk)
  list(CV.risk=CV.risk, adj.index=adj.index, 
       Adj=CAND.ADJ[[adj.index]], form=CAND.FORM[adj.index], var.CV=var.CV[adj.index])
}

#-----------------------------------------------------#-----------------------------------------------------
# getIC.CV: function to obtain a cross-validated estimate of the influence curve
#		See Balzer et al. "Adaptive pre-specification in randomized trials with and without pair-matching"
#
#	input: goal of analysis ('aRR' or 'RD), 
#   target (clust/inv)
#		indicator to break the match (break.match)
#		dataset (Ldata)
#		adjustment variable for the outcome regression (QAdj),
#   adjustment approach for outcome regression (Qform)
#		 adjustment variable for the pscore (gAdj),
#   adjustment approach for pscore regression (gform)
#	output: cross-validated estimate of the IC for pair
#-----------------------------------------------------#-----------------------------------------------------

# UPDATES
# previous version only did leave-one-out (unit or pair)
# this version generalizes to V-fold CV if V>=40
#     - can input the number of folds V (default=10)
#     - folds created stratified on binary outcomes (by default)
#     - if stratify=T and # observations in a given class is <V, 
#       then sets V=min observations in that fold
get.IC.CV<- function(goal, target, break.match, Ldata, QAdj, Qform, gAdj=NULL, gform='glm', 
                     scale_value, scale_value_min, folds, sample.effect=T){
  
  
  nFolds <- length(folds)
  DY.CV <- CV.risk <- NULL
  
  # doing a cross-validated estimate
  for(i in 1:nFolds) {
    
    these <- Ldata$indpt.unit %in% folds[[i]]  #******* IMPT
    valid <- Ldata[these, ]
    train <- Ldata[!these,]
    
    # run full TMLE algorithm on the training set
    train.out <- do.TMLE(goal=goal, target=target, sample.effect=sample.effect,
                         train=train, QAdj=QAdj, Qform=Qform, 
                         gAdj=gAdj, gform=gform,
                         scale_value=scale_value, scale_value_min=scale_value_min,
                         doing.CV=T, verbose=F)	
    
    # get the relevant components of the IC for the validation set, 
    # using fits based on the training set
    valid.out <- do.TMLE.validset(goal=goal, target=target, sample.effect=sample.effect,
                                  valid=valid, train.out=train.out,
                                  scale_value=scale_value, scale_value_min=scale_value_min)	

    # estimating the CV risk for each candidate
    # risk = Expectation of loss with loss as IC-sq
    # risk = variance of TMLE
    if(break.match){
      DY.CV <- c(DY.CV, valid.out$DY)
      CV.risk <- c(CV.risk, mean(valid.out$DY^2))
    }else{
      DY.CV <- c(DY.CV, valid.out$DY.paired)
      CV.risk <- c(CV.risk, mean(valid.out$DY.paired^2))
      
    }
  }
  
  # ave across folds
  CV.risk <- mean(CV.risk)
  
  # estimating the CV variance for that TMLE
  var.CV <- var(DY.CV)/length(DY.CV)
  
  list(CV.risk=CV.risk, var.CV=var.CV)
  
}



#-----------------------------------------------------#-----------------------------------------------------
# do.TMLE.for.valid: function to obtain a cross-validated estimate of the influence curve
#	for observations in the validation set
#		See Balzer et al. "Adaptive pre-specification in randomized trials with and without pair-matching"
#
#	input: goal of analysis ('aRR' or 'RD'),
#		validation dataset ('valid') 
#		TMLE-fits from training set (train.out)
#	output: cross-validated estimate of the IC,
#		cross-validated risk estimate (loss=IC^2)
#-----------------------------------------------------#-----------------------------------------------------

do.TMLE.validset <- function(goal, target, sample.effect, 
                             valid, train.out, scale_value, scale_value_min){
	
	# J <- length(unique(valid$id) )

	#=============================================
	# Step1 - initial estimation of E(Y|A,W)= Qbar(A,W)
	#=============================================

	valid<- do.Init.Qbar(train=valid, QAdj=train.out$QAdj, Qform=train.out$Qform, 
	                     glm.out=train.out$Q.out)$train
	
	#=============================================
	# Step2: Calculate the clever covariate
	#=============================================
	
	valid <- get.clever.cov(train=valid, gAdj=train.out$gAdj, gform=train.out$gform, 
	                        p.out=train.out$p.out)$train

	#=============================================
	# Step3: Targeting - 			
	#=============================================
	
	valid <- do.targeting(train=valid, eps=train.out$eps, goal=goal)
	
	#=============================================
	# Step5: Variance estimation using treatment-specific means from training set
	#=============================================

	 get.IC.variance(goal=goal, target=target, Vdata=valid, R1=train.out$R1, R0=train.out$R0, 
	                scale_value=scale_value, scale_value_min=scale_value_min, doing.CV=T, 
	                sample.effect=sample.effect)
	
}






#-----------------------------------------------------#-----------------------------------------------------
# get.cand.adj = function to get candidate adjustment strategies (variables + algorithms)
# for estimating the outcome regression and the propensity score

get.cand.adj <- function(all.cand, cand.Qform.fancy=NULL, cand.gform.fancy=NULL){
  
  all.cand.Ulist <- as.list(c('U', all.cand))
  
  # always consider main terms with each candidate adjustment variable
  cand.Qform <- cand.gform <-  rep('glm', length(all.cand.Ulist)) 
  
  if(is.null(cand.Qform.fancy)) {
    # simple Adaptive Prespec
    cand.QAdj <- all.cand.Ulist
  }else{
    # fancy adaptive prespec with expanded algorithms
    cand.QAdj <- c( all.cand.Ulist, rep(list(all.cand), length(cand.Qform.fancy)) )
    cand.Qform <- c(cand.Qform, cand.Qform.fancy)
  }
  
  if(is.null(cand.gform.fancy)){
    # simple adaptive prespec
    cand.gAdj <- all.cand.Ulist
  }else{
    # fancy adaptive prespec with expanded algorithms
    cand.gAdj <- c( all.cand.Ulist, rep(list(all.cand), length(cand.gform.fancy)) )
    cand.gform <- c(cand.gform, cand.gform.fancy)
  }
  
  list(cand.QAdj=cand.QAdj, cand.Qform=cand.Qform, 
       cand.gAdj=cand.gAdj, cand.gform=cand.gform)
  
  
}

#-----------------------------------------------------#-----------------------------------------------------

# adapted from .cvFolds from cvAUC package: https://CRAN.R-project.org/package=cvAUC
# by Erin LeDell 
get.folds <- function(V, Y, ids, stratify=T){
  
  if(stratify & length(unique(Y))==2 & length(Y)==length(unique(ids))){
    # stratify on the outcome
    # this option does not work for cluster randomized trials where # indpt units (IDs) < # rows in data
    classes <- tapply(1:length(Y), INDEX=Y, FUN=split, 1)
    ids.Y1 <- ids[classes$`1`]
    ids.noY1 <- ids[classes$`0`]
    if(length(ids.Y1) < V | length(ids.noY1) < V) {
      V <- min( length(ids.Y1), length(ids.noY1))
    }
    ids.Y1.split <- split(sample(length(ids.Y1)), rep(1:V, length=length(ids.Y1)))
    ids.noY1.split <- split(sample(length(ids.noY1)), rep(1:V, length=length(ids.noY1)))
    folds <- vector("list", V)
    for (v in seq(V)){
      folds[[v]] <- c(ids.Y1[ids.Y1.split[[v]]], ids.noY1[ids.noY1.split[[v]]])
    }
    
  }else{
    # dont stratify on the outcome
    ids.split <- split(sample(length(ids)), rep(1:V, length=length(ids)))
    folds <- vector("list", V)
    for (v in seq(V)){
      folds[[v]] <- ids[ids.split[[v]]]
    }
    
  }
  folds
}
