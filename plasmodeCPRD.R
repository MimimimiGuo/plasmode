PlasmodeBin2<- function(formulaOut=NULL, objectOut=NULL,formulaExp=NULL,objectExp=NULL,data, idVar,
                        effectOR =1, MMOut=1,MMExp=1, nsim, size, eventRate=NULL, exposedPrev=NULL)

{

  ## Code for simulating data when the outcome is binary and data set is provided to estimated the outcome and exposure.
  if(is.null(formulaOut)==FALSE & is.null(formulaExp)==FALSE & is.null(objectOut)==TRUE & is.null(objectExp)==TRUE)

  {
    outcome<- "frac_flag" ## selects the outcome variable
    exposure<- "ah_flag" ##selects the exposure variable

    x <- data[order(data[,exposure]),] # order according to exposure status, unexposed first
    n <- nrow(x)
    n1 <- sum(x[,exposure])     # number of exposed in real data
    n0 <- n - n1
    size1 <- round(ifelse(is.null(exposedPrev), n1*(size/n), size*exposedPrev))  # desired number of exposed
    size0 <- size - size1
    if(size1 > n1 | size0 > n0) stop("Number of requested exposed or unexposed exceeds observed number -- reduce size")

    # estimate logit model for probability of outcome
    modOutBin<- glm2(formulaOut, family = "binomial", data=x,control=glm.control(trace=TRUE))
    ## Design Matrix used for outcome logistic regression
    X <- model.matrix(modOutBin)

    # find event rate in base cohort
    if(is.null(eventRate)) eventRate <- mean(x[,outcome])

    # find intercept value needed to get approximate desired event rate under new parameters
    bnew <- c(coef(modOutBin)[1], MMOut*coef(modOutBin)[-1])
    bnew <- replace(bnew, names(coef(modOutBin)) == exposure, log(effectOR))
    Xbnew <- as.vector(X %*% bnew)
    fn <- function(d) mean(1 - 1/(1 + exp(d+Xbnew))) - eventRate
    delta <- uniroot(fn, interval = c(-20,20), lower = -20, upper = 20)$root

    ## Estimate logit model for probability of exposure
    modExp<- glm2(formulaExp, family = "binomial", data=x,control=glm.control(trace=TRUE))
    ## Design matrix used for exposure logistic regression
    XEXP<- model.matrix(modExp)

    # Finding the exposure prevalence in base cohort
    if(is.null(exposedPrev)) exposedPrev<- n1/n

    # Calculating intercept value needed to get approximate desired exposure prevalence under new parameters.
    bnewExp <- c(coef(modExp)[1], MMExp*coef(modExp)[-1])
    XbnewExp <- as.vector(XEXP%*%bnewExp)
    fnExp <- function(d) mean(plogis(d+XbnewExp))-exposedPrev
    deltaExp <- uniroot(fnExp, lower=-20, upper=20)$root
    Probexp <- plogis(deltaExp+XbnewExp)

    #### sample and simulate
    ids <- ynew <- expnew<-data.frame(matrix(nrow = size, ncol = nsim))
    RR<-RD<- vector('numeric', length = nsim)
    ##datasim<-list()
    for(sim in 1:nsim) {
      idxs0 <- sample(1:n0, size0, replace = TRUE) # sample unexposed (located in rows 1:n0 of x)
      ids[1:size0,sim] <- x[idxs0, idVar]
      idxs1 <- sample(n0+1:n1, size1, replace = TRUE) # sample exposed (located in rows n0 + 1:n1 of x)
      ids[size0+1:size1,sim] <- x[idxs1, idVar]
      X2 <- X[c(idxs0,idxs1),]
      expnew[,sim]<- X2[,2] <- rbinom(size,1,Probexp[c(idxs0,idxs1)])
      pnew <- plogis(delta + X2%*%bnew)
      ynew[,sim] <- rbinom(size, 1, pnew)

      datasim<-X[c(idxs0,idxs1),]
      datasim[,2]<-1
      p_1<- 1 - 1/(1+exp(as.vector(datasim %*% bnew + delta)))
      datasim[,2]<-0
      p_0<- 1 - 1/(1+exp(as.vector(datasim %*% bnew + delta)))
      RR[sim]<-mean(p_1)/mean(p_0)
      RD[sim]<-mean(p_1)-mean(p_0)
    }
    ARR<-mean(RR)
    ARD<-mean(RD)
    ## Creating simulated data for the outcome variable
    names(ids) <- paste("ID", 1:nsim, sep = "")
    names(ynew) <- paste("EVENT", 1:nsim, sep = "")
    names(expnew)<-paste("EXPOSURE",1:nsim, sep = "")
    sim_out_bin<-data.frame(ids, ynew,expnew)


    ## Finding the simulated treatment at new exposure prevalence.

    return(list(TrueOutBeta = bnew, TrueExpBeta = bnewExp, RR=ARR,RD=ARD,Sim_Data = sim_out_bin))
  }

  ## Code for simulating data when the outcome is binary and data set is provided to estimated the outcome.

  else if(is.null(formulaOut)==FALSE & is.null(formulaExp)==TRUE & is.null(objectOut)==TRUE & is.null(objectExp)==TRUE)

  {
    outcome <- all.vars(formulaOut)[1] ## selects the outcome variable
    exposure <- all.vars(formulaOut)[2] ##selects the exposure variable

    x <- data[order(data[,exposure]),] # order according to exposure status, unexposed first
    n <- nrow(x)
    n1 <- sum(x[,exposure])     # number of exposed in real data
    n0 <- n - n1
    size1 <- round(ifelse(is.null(exposedPrev), n1*(size/n), size*exposedPrev))  # desired number of exposed
    size0 <- size - size1
    if(size1 > n1 | size0 > n0) stop("Number of requested exposed or unexposed exceeds observed number -- reduce size")

    # estimate logit model for probability of outcome
    modOutBin <- glm2(formulaOut, family = "binomial", data=x,control=glm.control(trace=TRUE))
    ## Design Matrix used for outcome logistic regression
    X <- model.matrix(modOutBin)

    # find event rate in base cohort
    if(is.null(eventRate)) eventRate <- mean(x[,outcome])

    # find intercept value needed to get approximate desired event rate under new parameters
    bnew <- c(coef(modOutBin)[1], MMOut*coef(modOutBin)[-1])
    bnew <- replace(bnew, names(coef(modOutBin)) == exposure, log(effectOR))
    Xbnew <- as.vector(X %*% bnew)
    fn <- function(d) mean(1 - 1/(1 + exp(d+Xbnew))) - eventRate
    delta <- uniroot(fn, lower = -20, upper = 20)$root
    pnew <- 1 - 1/(1 + exp(delta+Xbnew))
    rm(modOutBin)

    ### sample and simulate
    ids <- ynew <-data.frame(matrix(nrow = size, ncol = nsim))
    RR <- RD <- vector('numeric', length = nsim)
    for(sim in 1:nsim) {
      idxs0 <- sample(1:n0, size0, replace = TRUE) # sample unexposed (located in rows 1:n0 of x)
      ids[1:size0,sim] <- x[idxs0, idVar]
      idxs1 <- sample(n0+1:n1, size1, replace = TRUE) # sample exposed (located in rows n0 + 1:n1 of x)
      ids[size0+1:size1,sim] <- x[idxs1, idVar]
      ynew[,sim] <- rbinom(size, 1, pnew[c(idxs0,idxs1)])
      datasim <-X[c(idxs0,idxs1),]
      datasim[,2] <- 1
      p_1 <- plogis(as.vector(datasim %*% bnew + delta))
      datasim[,2] <- 0
      p_0 <- plogis(as.vector(datasim %*% bnew + delta))
      RR[sim] <- mean(p_1)/mean(p_0)
      RD[sim] <- mean(p_1)-mean(p_0)
    }

    ARR <- mean(RR)
    ARD <- mean(RD)
    ## Creating simulated data for the outcome variable
    names(ids) <- paste("ID", 1:nsim, sep = "")
    names(ynew) <- paste("EVENT", 1:nsim, sep = "")
    sim_out_bin <- data.frame(ids, ynew)

    return(list(TrueOutBeta = bnew, RR=ARR,RD=ARD,Sim_Data = sim_out_bin))
  }

  ## Code for simulating data when the outcome is binary and data set is provided to estimated the exposure.
  else if(is.null(formulaOut)==TRUE & is.null(formulaExp)==FALSE & is.null(objectOut)==TRUE & is.null(objectExp)==TRUE)

  {
    exposure <- all.vars(formulaExp)[1] ##selects the exposure variable

    x <- data[order(data[,exposure]),] # order according to exposure status, unexposed first
    n <- nrow(x)
    n1 <- sum(x[,exposure])     # number of exposed in real data
    n0 <- n - n1
    size1 <- round(ifelse(is.null(exposedPrev), n1*(size/n), size*exposedPrev))  # desired number of exposed
    size0 <- size - size1
    if(size1 > n1 | size0 > n0) stop("Number of requested exposed or unexposed exceeds observed number -- reduce size")

    ## Estimate logit model for probability of exposure
    modExp <- glm2(formulaExp, family = "binomial", data=x,control=glm.control(trace=TRUE))
    ## Design matrix used for exposure logistic regression
    XEXP <- model.matrix(modExp)[order(data[,exposure]),]

    # Finding the exposure prevalence in base cohort
    if(is.null(exposedPrev)) exposedPrev <- n1/n

    # Calculating intercept value needed to get approximate desired exposure prevalence under new parameters.
    bnewExp <- c(coef(modExp)[1], MMExp*coef(modExp)[-1])
    XbnewExp <- as.vector(XEXP%*%bnewExp)
    fnExp <- function(d)mean(plogis(d+XbnewExp))-exposedPrev
    deltaExp <- uniroot(fnExp, lower=-20, upper=20)$root
    Probexp <- plogis(deltaExp+XbnewExp)
    rm(modExp, XEXP)

    #### sample and simulate
    ids <- expnew <- data.frame(matrix(nrow = size, ncol = nsim))
    for(sim in 1:nsim) {
      idxs0 <- sample(1:n0, size0, replace = TRUE) # sample unexposed (located in rows 1:n0 of x)
      ids[1:size0,sim] <- x[idxs0, idVar]
      idxs1 <- sample(n0+1:n1, size1, replace = TRUE) # sample exposed (located in rows n0 + 1:n1 of x)
      ids[size0+1:size1,sim] <- x[idxs1, idVar]
      expnew[,sim] <- rbinom(size,1,Probexp[c(idxs0,idxs1)])
    }

    ## Creating simulated data for the outcome variable
    names(ids) <- paste("ID", 1:nsim, sep = "")
    names(expnew) <- paste("EXPOSURE",1:nsim, sep = "")
    ##save(ids,ynew,file="simdata.RData")
    sim_out_bin <- data.frame(ids,expnew)

    return(list(TrueExpBeta = bnewExp, Sim_Data = sim_out_bin))
  }

  ## Code for simulating data when the outcome is binary and objects are provided to estimated the outcome and exposure.

  else if(is.null(formulaOut)==TRUE & is.null(formulaExp)==TRUE & is.null(objectOut)==FALSE & is.null(objectExp)==FALSE)

  {
    outcome <- all.vars(formula(objectOut))[1] ## selects the outcome variable
    exposure <- all.vars(formula(objectOut))[2] ##selects the exposure variable

    x <- data[order(data[,exposure]),] # order according to exposure status, unexposed first
    X <- model.matrix(objectOut)[order(data[,exposure]),]

    n <- nrow(x)
    n1 <- sum(x[,exposure])     # number of exposed in real data
    n0 <- n - n1
    size1 <- round(ifelse(is.null(exposedPrev), n1*(size/n), size*exposedPrev))  # desired number of exposed
    size0 <- size - size1
    if(size1 > n1 | size0 > n0) stop("Number of requested exposed or unexposed exceeds observed number -- reduce size")

    # find event rate in base cohort
    if(is.null(eventRate)) eventRate <- mean(x[,outcome])

    # find intercept value needed to get approximate desired event rate under new parameters
    bnew <- c(coef(objectOut)[1], MMOut*coef(objectOut)[-1])
    bnew <- replace(bnew, names(coef(objectOut)) == exposure, log(effectOR))
    Xbnew <- X%*%bnew
    fn <- function(d) mean(plogis(d+Xbnew)) - eventRate
    delta <- uniroot(fn, lower = -20, upper = 20)$root

    # Finding the exposure prevalence in base cohort
    if(is.null(exposedPrev)) exposedPrev<- n1/n
    XEXP <- model.matrix(objectExp)[order(data[,exposure]),]

    # Calculating intercept value needed to get approximate desired exposure prevalence under new parameters.
    ExpCoeff <- coef(objectExp)
    bnewExp <- c(ExpCoeff[1], MMExp*ExpCoeff[-1])
    XbnewExp <- as.vector(bnewExp%*%t(XEXP))
    fnExp <- function(d)mean(plogis(d+XbnewExp))-exposedPrev
    deltaExp <- uniroot(fnExp, lower=-20, upper=20)$root
    Probexp <- plogis(deltaExp+XbnewExp)

    #### sample and simulate
    ids <- ynew <- expnew <- data.frame(matrix(nrow = size, ncol = nsim))
    RR <- RD <- vector('numeric', length = nsim)
    for(sim in 1:nsim) {
      idxs0 <- sample(1:n0, size0, replace = TRUE) # sample unexposed (located in rows 1:n0 of x)
      ids[1:size0,sim] <- x[idxs0, idVar]
      idxs1 <- sample(n0+1:n1, size1, replace = TRUE) # sample exposed (located in rows n0 + 1:n1 of x)
      ids[size0+1:size1,sim] <- x[idxs1, idVar]
      X2 <- X[c(idxs0,idxs1),]
      expnew[,sim] <- X2[,2] <- rbinom(size,1,Probexp[c(idxs0,idxs1)])
      pnew <- plogis(delta + X2%*%bnew)
      ynew[,sim] <- rbinom(size, 1, pnew)

      datasim <- X[c(idxs0,idxs1),]
      datasim[,2] <- 1
      p_1 <- plogis(as.vector(datasim %*% bnew + delta))
      datasim[,2] <- 0
      p_0 <- plogis(as.vector(datasim %*% bnew + delta))
      RR[sim] <- mean(p_1)/mean(p_0)
      RD[sim] <- mean(p_1)-mean(p_0)
    }
    ARR<-mean(RR)
    ARD<-mean(RD)
    ## Creating simulated data for the outcome variable
    names(ids) <- paste("ID", 1:nsim, sep = "")
    names(ynew) <- paste("EVENT", 1:nsim, sep = "")
    names(expnew)<-paste("EXPOSURE",1:nsim, sep = "")
    sim_out_bin<-data.frame(ids, ynew,expnew)
    return(list(TrueOutBeta = bnew, TrueExpBeta = bnewExp, RR=ARR,RD=ARD,Sim_Data = sim_out_bin))
  }

  else if(is.null(formulaOut)==TRUE & is.null(formulaExp)==TRUE & is.null(objectOut)==FALSE & is.null(objectExp)==TRUE)

  {

    outcome <- all.vars(formula(objectOut))[1] ##selects the outcome variable
    exposure <- all.vars(formula(objectOut))[2] ##selects the exposure variable

    x <- data[order(data[,exposure]),]# order according to exposure status, unexposed first
    X <- model.matrix(objectOut)[order(data[,exposure]),]

    n <- nrow(x)
    n1 <- sum(x[,exposure])     # number of exposed in real data
    n0 <- n - n1
    size1 <- round(ifelse(is.null(exposedPrev), n1*(size/n), size*exposedPrev))  # desired number of exposed
    size0 <- size - size1
    if(size1 > n1 | size0 > n0) stop("Number of requested exposed or unexposed exceeds observed number -- reduce size")

    # find event rate in base cohort
    OutCoeff<-coef(objectOut)
    if(is.null(eventRate)) eventRate <- mean(x[,outcome])

    # find intercept value needed to get approximate desired event rate under new parameters
    bnew <- c(OutCoeff[1], MMOut*OutCoeff[-1])
    bnew <- replace(bnew, names(OutCoeff) == exposure, log(effectOR))
    Xbnew <- X%*%bnew
    fn <- function(d) mean(plogis(d+Xbnew)) - eventRate
    delta <- uniroot(fn, lower = -20, upper = 20)$root
    pnew <- plogis(delta+Xbnew)

    #### sample and simulate
    ids <- ynew <- data.frame(matrix(nrow = size, ncol = nsim))
    RR<-RD<- vector('numeric', length = nsim)
    for(sim in 1:nsim) {
      idxs0 <- sample(1:n0, size0, replace = TRUE) # sample unexposed (located in rows 1:n0 of x)
      ids[1:size0,sim] <- x[idxs0, idVar]
      idxs1 <- sample(n0+1:n1, size1, replace = TRUE) # sample exposed (located in rows n0 + 1:n1 of x)
      ids[size0+1:size1,sim] <- x[idxs1, idVar]
      ynew[,sim] <- rbinom(size, 1, pnew[c(idxs0,idxs1)])
      datasim <- X[c(idxs0,idxs1),]
      datasim[,2] <- 1
      p_1 <- plogis(as.vector(datasim %*% bnew + delta))
      datasim[,2] <- 0
      p_0 <- plogis(as.vector(datasim %*% bnew + delta))
      RR[sim] <- mean(p_1)/mean(p_0)
      RD[sim] <- mean(p_1)-mean(p_0)
    }
    ARR<-mean(RR)
    ARD<-mean(RD)
    ## Creating simulated data for the outcome variable
    names(ids) <- paste("ID", 1:nsim, sep = "")
    names(ynew) <- paste("EVENT", 1:nsim, sep = "")
    sim_out_bin <- data.frame(ids, ynew)
    return(list(TrueOutBeta = bnew, RR=ARR,RD=ARD,Sim_Data = sim_out_bin))
  }

  ## Code for simulating data when the outcome is binary and objects are provided to estimated the exposure.

  else if(is.null(formulaOut)==TRUE & is.null(formulaExp)==TRUE & is.null(objectOut)==TRUE & is.null(objectExp)==FALSE)

  {

    exposure <- all.vars(formula(objectExp))[1]
    x <- data[order(data[,exposure]),] # order according to exposure status, unexposed first

    n <- nrow(x)
    n1 <- sum(x[,exposure])     # number of exposed in real data
    n0 <- n - n1
    size1 <- round(ifelse(is.null(exposedPrev), n1*(size/n), size*exposedPrev))  # desired number of exposed
    size0 <- size - size1
    if(size1 > n1 | size0 > n0) stop("Number of requested exposed or unexposed exceeds observed number -- reduce size")

    # Finding the exposure prevalence in base cohort
    if(is.null(exposedPrev)) exposedPrev <- n1/n
    XEXP <- model.matrix(objectExp)[order(data[,exposure]),]

    # Calculating intercept value needed to get approximate desired exposure prevalence under new parameters.
    ExpCoeff <- coef(objectExp)
    bnewExp <- c(ExpCoeff[1], MMExp*ExpCoeff[-1])
    XbnewExp <- as.vector(bnewExp%*%t(XEXP))
    fnExp <- function(d)mean(plogis(d+XbnewExp))-exposedPrev
    deltaExp <- uniroot(fnExp, lower=-20, upper=20)$root
    Probexp <- plogis(deltaExp+XbnewExp)

    #### sample and simulate
    ids <- expnew <- data.frame(matrix(nrow = size, ncol = nsim))
    for(sim in 1:nsim) {
      idxs0 <- sample(1:n0, size0, replace = TRUE) # sample unexposed (located in rows 1:n0 of x)
      ids[1:size0,sim] <- x[idxs0, idVar]
      idxs1 <- sample(n0+1:n1, size1, replace = TRUE) # sample exposed (located in rows n0 + 1:n1 of x)
      ids[size0+1:size1,sim] <- x[idxs1, idVar]
      expnew[,sim]<- rbinom(size,1,Probexp[c(idxs0,idxs1)])
    }

    ## Creating simulated data for the outcome variable
    names(ids) <- paste("ID", 1:nsim, sep = "")
    names(expnew) <- paste("EXPOSURE",1:nsim, sep = "")
    sim_out_bin <- data.frame(ids, expnew)
    return(list(TrueExpBeta = bnewExp, Sim_Data = sim_out_bin))
  }

}



final_cohort_cov <- readRDS("D:/ah")
add_backticks = function(x) {
  paste0("`", x, "`")
}

x_lm_formula = function(x) {
  paste(add_backticks(x), collapse = "+")
}

build_lm_formula = function(x, y){
  if (length(y)>1){
    stop("y needs to be just one variable")
  }
  as.formula(
    paste0("`",y,"`", " ~ ", x_lm_formula(x))
  )
}

confColumns = colnames(final_cohort_cov)[c(5:54)] #500
instrColumns = colnames(final_cohort_cov)[c(55:74)] # 69
rfColumns = colnames(final_cohort_cov)[c(75:94)] # 69

#------------------select confounders, risk factors and instrumental variables------------------
exp_cols = "ah_flag"
out_cols = "frac_flag"
exp_x_cols = c(confColumns, instrColumns)
out_x_cols = c(confColumns, rfColumns, exp_cols)
formula_exp = build_lm_formula(exp_x_cols, exp_cols)
formula_out = build_lm_formula(out_x_cols, out_cols)
library(splitstackshape)

# final_cohort_cov_strat <- as.data.frame(stratified(final_cohort_cov, c('ah.flag'), 0.03))
final_cohort_cov$id <- c(1:nrow(final_cohort_cov))
library(glm2)

library(dplyr)
OriginalData <- final_cohort_cov %>% mutate(id = row_number())
Bin_Form1<-PlasmodeBin2(formulaOut=formula_out, objectOut=NULL,formulaExp=formula_exp,
                        objectExp= NULL,data=final_cohort_cov,idVar="id",effectOR = 1.5,
                        nsim=100,size=10000,exposedPrev=NULL,eventRate=NULL)

plasmodeData <- Bin_Form1$Sim_Data
#------------------save data------------------
for(i in 1:100)
{
  print(i)
  data <- NULL
  idname <- paste("ID", i, sep = "")
  evtname <- paste("EVENT", i, sep = "")
  expname <- paste("EXPOSURE", i, sep = "")
  data <- plasmodeData%>%select(evtname, expname, idname) %>% rename(id = idname)
  data <- OriginalData%>%left_join(data, by="id")
  form  = sprintf('data_%s.csv', i)
  write.csv(data,file = form)
}
