source("Stage2_Functions_Meta.R")
source("TMLE_Functions_Meta.R")
source("Adapt_Functions_Meta.R")
library(glmnet)
library(earth)
library(gt)
library(tidyverse)

library(speff2trial)
#' ACTG 175 was a randomized clinical trial to compare
#' monotherapy with zidovudine or didanosine with
#' combination therapy with zidovudine and didanosine or
#' zidovudine and zalcitabine
#' in adults infected with the human immunodeficiency virus type I
#' whose CD4 T cell counts were between
#' 200 and 500 per cubic millimeter.
help(ACTG175)
#' https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2562926/pdf/nihms46485.pdf
#' https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2600547/pdf/nihms-45052.pdf 
#' https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2574960/pdf/nihms-45063.pdf

#' The treatment effect is represented by the mean difference or
#' the log odds ratio for a quantitative or dichotomous endpoint, respectively.
#' Estimates of the treatment effect that ignore baseline covariates (naive)
#' are included in the output.

#' Using the automated model selection procedure performed by regsubsets,
#' four optimal regression models are developed for the study endpoint.
#' Initially, all baseline and postrandomization covariates
#' specified in the formula are considered for inclusion
#' by the model selection procedure
#' carried out separately in each treatment group.
#' The optimal models are used to construct predicted values of the endpoint.

#' Subsequently, in each treatment group, another regression model is fitted
#' that includes only baseline covariates
#' that were selected in the previous optimization.
#' Then predicted values of the endpoint
#' are computed based on these models.

#' If missingness occurs in the endpoint variable
#' , the model selection procedure is additionally used
#' to determine the optimal models for predicting
#' whether a subject has an observed endpoint
#' , separately in each treatment group.

#' The function regsubsets conducts optimization
#' of linear regression models only.
#' The following modification in the model selection
#' is adopted for a dichotomous variable:
#' initially, a logistic regression model is fitted
#' with all baseline and postrandomization covariates included in the formula.
#' Subsequently, an optimal model is selected
#' by using a weighted linear regression
#' with weights from the last iteration of the IWLS algorithm.
#' The optimal model is then refitted by logistic regression.

#' Besides using the built-in model selection algorithms,
#' the user has the option to explicitly enter predicted values of the endpoint
#' as well as estimated probabilities of observing the endpoint
#' if it is missing at random.

help(speff)
str(ACTG175)
mean(is.na(ACTG175$cd496))
table(is.na(ACTG175$cd420))

#========================================================================
#========================================================================
### CD4 T cell count at 96±5 weeks
#========================================================================
#========================================================================

### treatment effect estimation with a quantitative endpoint missing
### at random
fit1 <- speff(cd496 ~
    age + wtkg + hemo + homo + drugs + karnof + oprior + preanti
    + race + gender + str2 + strat + symptom + cd40 + cd80
    + cd420 + cd820 + offtrt
    , postrandom = c("cd420", "cd820", "offtrt")
    , data = ACTG175
    , trt.id = "treat")

summary(fit1)$tab

### 'fit2' adds quadratic effects of CD420 and CD820 and their
### two-way interaction
fit2 <- speff(cd496 ~
    age + wtkg + hemo + homo + drugs + karnof + oprior + preanti
    + race + gender + str2 + strat + symptom + cd40 + cd80
    + cd420 + cd820 + offtrt
    + I(cd420^2)
    + I(cd820^2) + cd420:cd820
    , postrandom = c("cd420", "cd820", "offtrt"
    , "I(cd420^2)", "I(cd820^2)", "cd420:cd820")
    , data = ACTG175,
trt.id = "treat")

### 'fit3' uses R-squared as the optimization criterion
fit3 <- speff(cd496 ~
    age + wtkg + hemo + homo + drugs + karnof + oprior + preanti
    + race + gender + str2 + strat + symptom + cd40 + cd80
    + cd420 + cd820 + offtrt
    , postrandom = c("cd420", "cd820", "offtrt")
    , data = ACTG175
    , trt.id = "treat"
    , optimal = "rsq")


###Adaptive Selection of the Optimal Strategy 

all_cand <-
    c("age", "wtkg", "hemo", "homo", "drugs"
    , "karnof", "oprior", "preanti", "race"
    , "gender", "str2", "strat", "symptom"
    , "cd40", "cd80"
    # exclude postrandomizaion variables
     #, "cd420",  "cd820", "offtrt"
    )

#' Estimator 1: Simple Adaptive Prespecification
#' This estimator will automatically consider
#' GLMs with a main term for one element of all.cand (and nothing=U)
ap_simple <- get.cand.adj(all.cand = all_cand
    , cand.Qform.fancy = NULL, cand.gform.fancy = NULL)

#' Estimator 2: Fancy Adaptive Prespecification
#' For this estimator, we need to specify candidate algorithms 
#' for adjusting for multiple covariates
#' By default, we specify "glm", "stepwise"
#' , "step.interaction", "lasso" and optionally "mars"
cand_gform_fancy <-
    cand_qform_fancy <-
        c("glm", "stepwise", "step.interaction", "lasso", "mars")
ap_fancy  <- get.cand.adj(all.cand = all_cand
    , cand.Qform.fancy = cand_qform_fancy
    , cand.gform.fancy = cand_gform_fancy)

data_input     <- ACTG175[, all_cand]
data_input$id  <- ACTG175$pidnum
data_input$U   <-
    data_input$alpha <- 1
data_input$A   <- ACTG175$treat
data_input$Y   <- ACTG175$cd496
data_input  <- na.omit(data_input)

these <- c("est", "se", "CI.lo", "CI.hi", "pval")

unadj <- Stage2(goal = "RD"
    , data.input = data_input
    , do.data.adapt = FALSE
    , do.cv.variance = FALSE
    , one.sided = FALSE
    , verbose = FALSE)

simple <- Stage2(goal = "RD"
    , data.input = data_input
    , do.data.adapt = TRUE
    , do.cv.variance = FALSE
    , V = 5
    , one.sided = FALSE
    , cand.QAdj =  ap_simple$cand.QAdj
    , cand.Qform = ap_simple$cand.Qform
    , cand.gAdj =  ap_simple$cand.gAdj
    , cand.gform = ap_simple$cand.gform
    , verbose = FALSE)


fancy <- Stage2(goal = "RD"
    , data.input  = data_input
    , do.data.adapt = TRUE
    , do.cv.variance = FALSE
    , V = 5
    , one.sided = FALSE
    , cand.QAdj =  ap_fancy$cand.QAdj
    , cand.Qform = ap_fancy$cand.Qform
    , cand.gAdj =  ap_fancy$cand.gAdj
    , cand.gform = ap_fancy$cand.gform
    , verbose = FALSE)

summary(fit1)$tab

tb <-
    dplyr::bind_rows(
        unadj[, these]
        , setNames(summary(fit1)$tab[2, ], these)
        , setNames(summary(fit2)$tab[2, ], these)
        , setNames(summary(fit3)$tab[2, ], these)
        , simple[, these]
        , fancy[, these]
        )
rownames(tb) <- c("Unadjusted"
    , "speff linear"
    , "speff quadratic"
    , "speff linear R^2"
    , "Simple"
    , "Fancy"
    )

tb %>%
    gt(rownames_to_stub = TRUE) %>%
    tab_header(title = "CD4 T cell count at 96±5 weeks") %>%
    gtsave(paste0("OUTPUT/contY96_ACTG175.html"))

#====================================================================
### a dichotomous response is created with missing values maintained
ACTG175$cd496bin <- ifelse(ACTG175$cd496 > 250, 1, 0)

### treatment effect estimation with a dichotomous endpoint missing
### at random
fit4 <- speff(cd496bin ~
    age + wtkg + hemo + homo + drugs + karnof + oprior + preanti
    + race + gender + str2 + strat + symptom + cd40 + cd80
    + cd420 + cd820 + offtrt
    , postrandom = c("cd420", "cd820", "offtrt")
    , data = ACTG175
    , trt.id = "treat"
    , endpoint = "dichotomous")


data_input$Y   <- na.omit(ACTG175$cd496bin)

unadj_b <- Stage2(goal = "aRR"
    , data.input = data_input
    , do.data.adapt = FALSE
    , do.cv.variance = FALSE
    , one.sided = FALSE
    , verbose = FALSE)
unadj_b

simple_b <- Stage2(goal = "aRR"
    , data.input = data_input
    , do.data.adapt = TRUE
    , do.cv.variance = FALSE
    , V = 5
    , one.sided = FALSE
    , cand.QAdj =  ap_simple$cand.QAdj
    , cand.Qform = ap_simple$cand.Qform
    , cand.gAdj =  ap_simple$cand.gAdj
    , cand.gform = ap_simple$cand.gform
    , verbose = FALSE)

fancy_b <- Stage2(goal = "aRR"
    , data.input  = data_input
    , do.data.adapt = TRUE
    , do.cv.variance = FALSE
    , V = 5
    , one.sided = FALSE
    , cand.QAdj =  ap_fancy$cand.QAdj
    , cand.Qform = ap_fancy$cand.Qform
    , cand.gAdj =  ap_fancy$cand.gAdj
    , cand.gform = ap_fancy$cand.gform
    , verbose = FALSE)

s4 <- summary(fit4)$tab[2, ]
s4 [1:4] <- exp(s4 [1:4])
s4

tb_b <-
    dplyr::bind_rows(
      unadj_b[, these]
    , setNames(s4, these)
    , simple_b[, these]
    , fancy_b[, these]
    )

rownames(tb_b) <- c("Unadjusted"
    , "speff"
    , "Simple"
    , "Fancy"
    )

tb_b %>%
    gt(rownames_to_stub = TRUE) %>%
    tab_header(title = "CD4 T cell count at 96±5 weeks > 250") %>%
    gtsave(paste0("OUTPUT/binY96_ACTG175.html"))

#========================================================================
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#========================================================================
### CD4 T cell count at 20±5 weeks
#========================================================================
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#========================================================================

### treatment effect estimation with a quantitative endpoint missing
### at random
f1 <- speff(cd420 ~
    age + wtkg + hemo + homo + drugs + karnof + oprior + preanti
    + race + gender + str2 + strat + symptom + cd40 + cd80
    , data = ACTG175
    , trt.id = "treat")

summary(f1)$tab

### 'f2' uses R-squared as the optimization criterion
f2 <- speff(cd420 ~
    age + wtkg + hemo + homo + drugs + karnof + oprior + preanti
    + race + gender + str2 + strat + symptom + cd40 + cd80
    , data = ACTG175
    , trt.id = "treat"
    , optimal = "rsq")


###Adaptive Selection of the Optimal Strategy
data_input     <- ACTG175[, all_cand]
data_input$id  <- ACTG175$pidnum
data_input$U   <-
    data_input$alpha <- 1
data_input$A   <- ACTG175$treat
data_input$Y   <- ACTG175$cd420
data_input  <- na.omit(data_input)

these <- c("est", "se", "CI.lo", "CI.hi", "pval")

unadj20 <- Stage2(goal = "RD"
    , data.input = data_input
    , do.data.adapt = FALSE
    , do.cv.variance = FALSE
    , one.sided = FALSE
    , verbose = FALSE)

simple20 <- Stage2(goal = "RD"
    , data.input = data_input
    , do.data.adapt = TRUE
    , do.cv.variance = FALSE
    , V = 5
    , one.sided = FALSE
    , cand.QAdj =  ap_simple$cand.QAdj
    , cand.Qform = ap_simple$cand.Qform
    , cand.gAdj =  ap_simple$cand.gAdj
    , cand.gform = ap_simple$cand.gform
    , verbose = FALSE)


fancy20 <- Stage2(goal = "RD"
    , data.input  = data_input
    , do.data.adapt = TRUE
    , do.cv.variance = FALSE
    , V = 5
    , one.sided = FALSE
    , cand.QAdj =  ap_fancy$cand.QAdj
    , cand.Qform = ap_fancy$cand.Qform
    , cand.gAdj =  ap_fancy$cand.gAdj
    , cand.gform = ap_fancy$cand.gform
    , verbose = FALSE)

tb <-
    dplyr::bind_rows(
        unadj20[, these]
        , setNames(summary(f1)$tab[2, ], these)
        , setNames(summary(f2)$tab[2, ], these)
        , simple20[, these]
        , fancy20[, these]
        )
rownames(tb) <- c("Unadjusted"
    , "speff cp"
    , "speff R^2"
    , "Simple"
    , "Fancy"
    )

tb %>%
    gt(rownames_to_stub = TRUE) %>%
    tab_header(title = "CD4 T cell count at 20±5 weeks") %>%
    gtsave(paste0("OUTPUT/contY20_ACTG175.html"))

#====================================================================
### a dichotomous response is created with missing values maintained
ACTG175$cd420bin <- ifelse(ACTG175$cd420 > 250, 1, 0)

### treatment effect estimation with a dichotomous endpoint missing
### at random
f4 <- speff(cd420bin ~
    age + wtkg + hemo + homo + drugs + karnof + oprior + preanti
    + race + gender + str2 + strat + symptom + cd40 + cd80
    , data = ACTG175
    , trt.id = "treat"
    , endpoint = "dichotomous")


data_input$Y   <- na.omit(ACTG175$cd420bin)

unadj_b20 <- Stage2(goal = "aRR"
    , data.input = data_input
    , do.data.adapt = FALSE
    , do.cv.variance = FALSE
    , one.sided = FALSE
    , verbose = FALSE)

simple_b20 <- Stage2(goal = "aRR"
    , data.input = data_input
    , do.data.adapt = TRUE
    , do.cv.variance = FALSE
    , V = 5
    , one.sided = FALSE
    , cand.QAdj =  ap_simple$cand.QAdj
    , cand.Qform = ap_simple$cand.Qform
    , cand.gAdj =  ap_simple$cand.gAdj
    , cand.gform = ap_simple$cand.gform
    , verbose = FALSE)

fancy_b20 <- Stage2(goal = "aRR"
    , data.input  = data_input
    , do.data.adapt = TRUE
    , do.cv.variance = FALSE
    , V = 5
    , one.sided = FALSE
    , cand.QAdj =  ap_fancy$cand.QAdj
    , cand.Qform = ap_fancy$cand.Qform
    , cand.gAdj =  ap_fancy$cand.gAdj
    , cand.gform = ap_fancy$cand.gform
    , verbose = FALSE)

s420 <- summary(f4)$tab[2, ]
s420 [1:4] <- exp(s420 [1:4])
s420

tb_b20 <-
    dplyr::bind_rows(
      unadj_b20[, these]
    , setNames(s420, these)
    , simple_b20[, these]
    , fancy_b20[, these]
    )

rownames(tb_b20) <- c("Unadjusted"
    , "speff"
    , "Simple"
    , "Fancy"
    )

tb_b20 %>%
    gt(rownames_to_stub = TRUE) %>%
    tab_header(title = "CD4 T cell count at 20±5 weeks > 250") %>%
    gtsave(paste0("OUTPUT/binY20_ACTG175.html"))

#=======================================================================