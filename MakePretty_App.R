# Helper function to make the output from the real data analysis look pretty

make.pretty.mini <- function(pt, lo, hi, scaler, digit, in.percent){
  paste0( round(pt*scaler, digit), ifelse(in.percent,'%',''),
          ' (', #' (95%CI: ',
          round(lo*scaler, digit), '-',
          round(hi*scaler, digit), ifelse(in.percent,'%',''), ')'
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
    round( var.base/(est$se^2), 2)
  )
  )
  colnames(yay) <- c('Intervention', 'Control', 'Effect', 'Precision')
  yay
}

make.pretty.wrapper <- function(est, var.base, these.rows, scaler=1, digit=1, in.percent=F){
  
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


