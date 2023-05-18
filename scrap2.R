rm(list=ls())
source("Stage2_Functions_Meta.R")
source("TMLE_Functions_Meta.R")
source("Adapt_Functions_Meta.R")
library("SuperLearner")
library("glmnet")
library("earth")
library('knitr')

library("speff2trial")
# help(ACTG175)
data_input <- ACTG175
# subset the data  aged 18+
data_input <- data_input[data_input$age >17,]
# create indicators
data_input$young <- as.numeric( data_input$age < 30)
data_input$cd40bin <- as.numeric(data_input$cd40 > 350)
data_input$cd80bin <- as.numeric(data_input$cd80 > 350)
data_input$recent <- as.numeric(data_input$strat==2) 

data_input$id  <- data_input$pidnum # patient id
data_input$U   <- 1
data_input$alpha <- 1
data_input$A   <- data_input$treat 


all_cand <- c("age", "young", "wtkg", "hemo",
              "karnof", "oprior", "preanti", 
              "race",  "gender", 
              "str2", "recent",  "symptom",
              "cd40", "cd40bin", "cd80", "cd80bin")

small_aps <- get.cand.adj(all.cand = all_cand, cand.Qform.fancy = NULL, cand.gform.fancy = NULL)

large_aps  <- get.cand.adj(all.cand = all_cand, 
                           cand.Qform.fancy = c("glm", "stepwise", "lasso", "mars", "mars.corP"), 
                           cand.gform.fancy = c("glm", "stepwise", "lasso", "mars", "mars.corP"))
set.seed(1)

nReps <- 5000
sim.cols <- c('est', 'bias', 'cover','reject')
UNADJ <-  data.frame(matrix(NA, nrow=nReps, ncol=length(sim.cols)))
colnames(UNADJ) <- sim.cols
FIXED <- SMALL <- BIG <- UNADJ

do.cont <- T # specify if continuous or binary outcome

if(do.cont){
  data_input$Y   <- data_input$cd420
  goal <- 'RD'
  psi <- 0 
  file.name <- 'ACTG_null_cont_5000.Rdata'
} else{
  data_input$Y <- as.numeric(data_input$cd420 > 350)
  psi <- 1
  goal <- 'aRR'
  file.name <- 'ACTG_null_bin_5000.Rdata'
}

dt <- data_input
for(r in 1:nReps){
  # randomly permute the treatment
  dt$A <- sample(dt$A)
  # implement 3 estimators 
  unadj <- suppressWarnings( Stage2(goal = goal, data.input = dt, do.data.adapt =F, psi=psi))
  fixed <- suppressWarnings( Stage2(goal = goal, data.input = dt, 
                                    do.data.adapt = F, 
                                    QAdj='age', Qform='glm', 
                                    gAdj='gender', gform='glm',  psi=psi))
  small_tmle <- suppressWarnings( Stage2(goal = goal, data.input = dt, 
                                         do.data.adapt = TRUE, V = 5, 
                                         cand.QAdj =  small_aps$cand.QAdj, cand.Qform = small_aps$cand.Qform,
                                         cand.gAdj =  small_aps$cand.gAdj, cand.gform = small_aps$cand.gform,
                                         psi=psi))
  
  large_tmle <- suppressWarnings( Stage2(goal = goal, data.input = dt, 
                                         do.data.adapt = TRUE, V = 5, 
                                         cand.QAdj =  large_aps$cand.QAdj, cand.Qform = large_aps$cand.Qform,
                                         cand.gAdj =  large_aps$cand.gAdj, cand.gform = large_aps$cand.gform,
                                         psi=psi))
  # save the output
  UNADJ[r,] <- unadj[,sim.cols]
  FIXED[r,] <- fixed[,sim.cols]
  SMALL[r,] <- small_tmle[,sim.cols]
  BIG[r, ]  <- large_tmle[,sim.cols]
  print(r)
}

save(UNADJ, FIXED, SMALL, BIG, file=file.name)
