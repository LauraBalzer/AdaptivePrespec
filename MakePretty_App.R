# Helper function to make the output from the real data analysis look pretty

make.pretty.mini <- function(pt, lo, hi, scaler, digit, in.percent){
  this <- paste0("%.",digit,"f")
  
  paste0( sprintf(this, pt*scaler), ifelse(in.percent,'%',''),
          ' (', #' (95%CI: ',
          sprintf(this, lo*scaler), ', ',
          sprintf(this, hi*scaler), #ifelse(in.percent,'%',''), 
          ')'
  )
}   


make.pretty.app <- function(est, scaler=1, digit=1, in.percent=F, var.base){
 
  yay <- data.frame( cbind( 
    make.pretty.mini(est$Txt.est, est$Txt.CI.lo, est$Txt.CI.hi, 
                      scaler=scaler, digit=digit, in.percent=in.percent),
     make.pretty.mini(est$Con.est, est$Con.CI.lo, est$Con.CI.hi, 
                      scaler=scaler, digit=digit, in.percent=in.percent),
     make.pretty.mini(est$est, est$CI.lo, est$CI.hi, 
                      scaler=1, digit=ifelse(in.percent, (digit+1), digit), in.percent=F),
    #ifelse (est$pval<0.001, '<0.001', round(est$pval, 3))
    sprintf("%.3f", (est$se^2)/var.base ),
    paste0( round(( 1- (est$se^2)/var.base)*100, 1), '%')
  )
  )
  colnames(yay) <- c('Intervention', 'Control', 'Effect', 'Rel.Var.', 'Savings')
  yay
}

make.pretty.wrapper <- function(est, var.base, these.rows, scaler=1, digit=0, in.percent=F){
  
  yay <- NULL
  for(k in 1:nrow(est)){
    yay <- rbind(yay,
                 make.pretty.app(est[k,], scaler=scaler, digit=digit, in.percent = in.percent,
                                   var.base=var.base)
                   )
    
  }
  row.names(yay)<- these.rows
  yay
}



aps_wrapper <- function(goal, data_input, V=5, small_aps, large_aps){
  
  unadj <- Stage2(goal = goal, data.input = data_input, do.data.adapt =F)
  # fixed adjustment
  fixed <- Stage2(goal = goal, data.input = data_input, 
                  do.data.adapt = F, 
                  QAdj='age', Qform='glm', 
                  gAdj='gender', gform='glm')
 
  # TMLE with small APS
  small_tmle <- Stage2(goal = goal, data.input = data_input, 
                       do.data.adapt = TRUE, V = V, 
                       cand.QAdj =  small_aps$cand.QAdj, cand.Qform = small_aps$cand.Qform,
                       cand.gAdj =  small_aps$cand.gAdj, cand.gform = small_aps$cand.gform)
  
  # TMLE adjusting in outcome regression with small APS selection and unadjusted pscore
  small_tmle_Qonly <- Stage2(goal = goal, data.input = data_input, 
                             # do.data.adapt = F, V = 5, 
                             QAdj= unlist(small_aps$cand.QAdj[small_tmle$QAdj]), 
                             Qform=small_tmle$Qform, 
                             gAdj=NULL, gform='glm')
  
  # TMLE with large APS
  large_tmle <- Stage2(goal = goal, data.input = data_input, 
                       do.data.adapt = TRUE, V = V, 
                       cand.QAdj =  large_aps$cand.QAdj, cand.Qform = large_aps$cand.Qform,
                       cand.gAdj =  large_aps$cand.gAdj, cand.gform = large_aps$cand.gform)
  
  # TMLE adjusting in outcome regression with large APS selection and unadjusted pscore
  large_tmle_Qonly <- Stage2(goal = goal, data.input = data_input, 
                             # do.data.adapt = F, V = 5, 
                             QAdj= unlist(large_aps$cand.QAdj[large_tmle$QAdj]), 
                             Qform=large_tmle$Qform, 
                             gAdj=NULL, gform='glm')
  # data frame
  est <- data.frame(rbind(unadj, fixed, small_tmle_Qonly, small_tmle,
                          large_tmle_Qonly, large_tmle))
  rownames(est) <- c('Unadjusted','Static', 'Small TMLE', 'Small CTMLE', 
                     'Large TMLE', 'Large CTMLE')
  
  compact <- make.pretty.wrapper(est=est,
                                 # variance estimate for precision comparison
                                 var.base = (unadj$se^2), 
                                 these.rows=rownames(est),
                                 digit=ifelse(goal=='RD',1,2) )
  
  list(unadj=unadj, fixed=fixed,  small_tmle=small_tmle, small_tmle_Qonly=small_tmle_Qonly, 
       large_tmle=large_tmle, large_tmle_Qonly=large_tmle_Qonly,
       est=est, compact=compact)
}

