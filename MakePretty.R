
rm(list = ls())
library(ggplot2)
library(hrbrthemes)
n <- 500
nReps <- 500
V <- 5
incl.mars <- T
verbose <- F

#sim <- "binY"
 sim <- "contY"
effect <- T

if(sim=='contY'){
  expt_type <- c('noisy_linear_1','noisy_multicollinear_cand1', 'noisy_polynomial')
  null.value=0
} else{
  expt_type <- c('noisy_linear','noisy_multicollinear', 'noisy_polynomial')
  
  null.value=1
}




get.MSE <- function(output){
  mean( (output$est - output$psi)^2 ,na.rm=TRUE )
}

get.metrics <- function(estimator){
  # cover: 95% confidence interval contained the truth?
  # reject: null hypo of no effect rejected
  # Bias: ave deviation between pt and truth
  # Variance: variance of point estimates
  # MSE
  yay <- c( colMeans(estimator[,c('cover', 'reject', 'bias')], na.rm=T),
        var(as.numeric(unlist(estimator["est"])),na.rm=TRUE),
        get.MSE(estimator)
        )
  yay <- data.frame(t(yay))
  colnames(yay) <- c('cover','power','bias', 'var', 'mse')
  yay
}

get.selection <- function(cand, this.var, this.form ){
  winner <- data.frame(matrix(0, nrow=1, ncol=length(cand)+1)) 
  colnames(winner) <- c('unadj',cand)
  winner['unadj'] <- sum(this.var==1 & this.form=='glm')
  winner['glm']<- sum(this.var!=1 & this.form=='glm')
  winner['stepwise'] <- sum(this.form=='stepwise')
  winner['step.interaction'] <- sum(this.form=='step.interaction')
  winner['lasso'] <- sum(this.form=='lasso')
  winner['mars'] <- sum(this.form=='mars')
  winner
}


YAY <- NULL
ests <- c('Unadjusted', 'Static', 'Small APS', 'Large APS')
STRATIFY <- c(F,T)
dgp <- c('Linear', 'Interactive', 'Polynomial')

WINNERQ <- WINNERG <- NULL


for(j in 1:length(expt_type)){
  for(k in 1:2){
  file.name <- paste( sim, paste0('Effect', effect),
                      paste0('N', n), paste0('V',V), paste0('mars', incl.mars),
                      paste0('nReps', nReps),paste0('stratify', STRATIFY[k]),
                      paste0('type', expt_type[j]), sep = "_")
  
  file.nameD <- paste( "OUTPUT/", file.name, paste('.RData'), sep = "_")
  print(paste0("Experiment file name is: ", file.nameD))
  load(file.nameD)
  
  SIMPLE <- OUT.AP
  FANCY <- OUT 
  yay <- data.frame(rbind( 
            get.metrics(UNADJ), get.metrics(FORCE), 
            get.metrics(SIMPLE),
            get.metrics(FANCY)))
  yay <- cbind(expt=expt_type[j], stratify=STRATIFY[k], 
               ests, yay, var.ratio=yay[1,'var']/yay[,'var'], re=yay[,'mse']/yay[1,'mse'] )
  yay <- cbind(yay, savings=(1-yay$re))
  print(round(mean(UNADJ$psi),2))
  YAY <- rbind(YAY, yay)
  
  winnerq <- cbind(expt=expt_type[j], stratify=STRATIFY[k], 
                  get.selection(cand=unique(AP.fancy$cand.Qform),
                    this.var= SELECT$QAdj,
                    this.form= SELECT$Qform)
  )
  winnerq[, 3:8] <- paste0(round(winnerq[,3:8]/500*100, 1), '%')
  WINNERQ <- rbind(WINNERQ, winnerq)
  
  
  winnerg <- cbind(expt=expt_type[j], stratify=STRATIFY[k],
                   get.selection(cand=unique(AP.fancy$cand.gform),
                                 this.var= SELECT$gAdj,
                                 this.form= SELECT$gform)
  )
  winnerg[, 3:8] <- paste0(round(winnerg[,3:8]/500*100, 1), '%')
  WINNERG <- rbind(WINNERG, winnerg)
  
  data <- data.frame(
    x=c(1:2000), 
    value1=c(UNADJ[["CI.lo"]],FORCE[["CI.lo"]],SIMPLE[["CI.lo"]],FANCY[["CI.lo"]]), 
    value2=c(UNADJ[["CI.hi"]],FORCE[["CI.hi"]],SIMPLE[["CI.hi"]],FANCY[["CI.hi"]]),
    ests=c(rep('Unadjusted',500), rep('Static',500), rep('Small APS', 500), rep('Large APS', 500))
  )
  
  psi <- mean(UNADJ$psi)
  # Plot
  ggplot(data) +
    geom_segment( aes(x=x, xend=x, y=value1, yend=value2), color="grey") +
    #geom_point( aes(x=x, y=value1), color=rgb(0.2,0.7,0.1,0.5), size=3 ) +
    #geom_point( aes(x=x, y=value2), color=rgb(0.7,0.2,0.1,0.5), size=3 ) +
    geom_point( aes(x=x, y=value1, color=factor(ests)), size=1 ) +
    geom_point( aes(x=x, y=value2, color=factor(ests)), size=1 ) +
    geom_hline(yintercept=null.value, linetype='dashed', color='black', size=0.5) +
    coord_flip()+
    theme_ipsum() +
    labs(title = dgp[j])+
    theme(
      legend.position = "none",
      #plot.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())
  ggsave(filename = paste0(file.name, '.eps'))
  rm(yay, winnerq, winnerg, data)

}
}
library(xtable)
this.order <- c('ests', 'cover', 'power', 'mse', 'bias', 'var', 're')
YAY
xtable(YAY[,this.order], digits=c(1, 1, rep(3, 6) ))


round( summary(YAY[YAY$ests=='Large APS', 're']), 3)
round( summary(YAY[YAY$ests=='Large APS', 'savings'])*100, 0)

round( summary(YAY[YAY$ests=='Small APS', 're']), 3)
round( summary(YAY[YAY$ests=='Small APS', 'savings'])*100, 0)


xtable(WINNERQ)
xtable(WINNERG)
