# Helper function to make the output from the real data analysis look pretty

make.pretty.mini <- function(pt, lo, hi, scaler, digit, in.percent){
  this <- paste0("%.",digit,"f")
  
  paste0( sprintf(this, pt*scaler), ifelse(in.percent,'%',''),
          ' (', #' (95%CI: ',
          sprintf(this, lo*scaler), '-',
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
    sprintf("%.2f", (est$se^2)/var.base ),
    paste0( round(( 1- (est$se^2)/var.base)*100, 0), '%')
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


