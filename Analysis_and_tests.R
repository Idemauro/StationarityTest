### Entomological data Analysis #####
### Idemauro Antonio Rodrigues de Lara - Cesar Taconeli ####
rm(list=ls(all=TRUE))
library(ordinal)
library(Matrix)

################################################################################
### Function for format conversion (wide to long)
widetolong <-function(basewide,varswide){
  r1<-reshape(basewide,varying=list(varswide),direction="long",
              v.names="resp")
  r1$time<-factor(r1$time)
  rownames(r1)<-1:nrow(r1)
  r1<-r1[,-which(names(r1)=="id")]
  r2<-r1[-which(r1$time==levels(r1$time)[1]),]
  r2<-cbind(r2,r1[-which(r1$time==levels(r1$time)[length(levels(r1$time))])
                  ,"resp"])
  names(r2)[ncol(r2)]<-"respp";head(r2)
  return(r2)
}


################################################################################
# Loading the dataset
setwd("C:/Users/cetac/Dropbox/Idemauro_gradiente")
dados <- read.csv2("cat3-4.csv",sep=";")
dados <- dados[,-1]


### Simulating new observations for the augmented dataset
set.seed(2020)
n <- 32

dadossomasp0 <- data.frame(cbind(dados$tratamento[1:16],
sample(1:3, n, replace = TRUE, table(dados$resp1[1:16])/16),
sample(1:3, n, replace = TRUE, table(dados$resp2[1:16])/16),
sample(1:3, n, replace = TRUE, table(dados$resp2[1:16])/16),
sample(1:3, n, replace = TRUE, table(dados$resp2[1:16])/16),
sample(1:3, n, replace = TRUE, table(dados$resp2[1:16])/16)))

colnames(dadossomasp0) <- c("tratamento","resp1","resp2","resp3",
                            "resp4","resp5")

dadossomasp1 <- data.frame(cbind(dados$tratamento[17:32],
  sample(1:3, n, replace = TRUE, table(dados$resp1[17:32])/16),
  sample(1:3, n, replace = TRUE, table(dados$resp1[17:32])/16),
  sample(1:3, n, replace = TRUE, table(dados$resp1[17:32])/16),
  sample(1:3, n, replace = TRUE, table(dados$resp1[17:32])/16),
  sample(1:3, n, replace = TRUE, table(dados$resp1[17:32])/16)))

colnames(dadossomasp1) <- c("tratamento","resp1","resp2","resp3",
                            "resp4","resp5")

dadosrep <-rbind(dados,dadossomasp0,dadossomasp1)
### dadosrep is the augmented dataset


### Long format conversion
dadoslong <- widetolong(dadosrep,c("resp1","resp2","resp3","resp4","resp5"))
dadoslong$tratamento <- as.factor(dadoslong$tratamento)
dadoslong$resp <- as.factor(dadoslong$resp)
dadoslong$respp <- as.factor(dadoslong$respp)

dadosrep$tratamento <- as.factor(dadosrep$tratamento)
dadosrep$resp1 <- as.factor(dadosrep$resp1)
dadosrep$resp2 <- as.factor(dadosrep$resp2)
dadosrep$resp3 <- as.factor(dadosrep$resp3)
dadosrep$resp4 <- as.factor(dadosrep$resp4)
dadosrep$resp5 <- as.factor(dadosrep$resp5)

################################################################################
### Fitting the stationary transition model

modE <- clm(resp ~ 1 , nominal = ~ tratamento + respp, data = dadoslong)
beta0 <- coef(modE) ### Regression parameter estimates.
vcov0 <- vcov(modE) ### Asymptotic covariance matrix.
loglik0 <- as.numeric(logLik(modE)) ### Maximized log-likelihood.


### Fitting the non-stationary transition model

### First Transition t1-t2
mod12 <- clm(resp2 ~ 1 , nominal = ~ tratamento + resp1, data=dadosrep) ### Fitted model.
beta12 <- coef(mod12) ### Regression parameter estimates.
vcov12 <- vcov(mod12) ### Asymptotic covariance matrix.
loglik12 <- as.numeric(logLik(mod12)) ### Maximized log-likelihood.

### Second Transition t2-t3
mod23 <- clm(resp3 ~ 1 , nominal = ~ tratamento + resp2, data=dadosrep)
beta23 <- coef(mod23)
vcov23 <- vcov(mod23)
loglik23 <- as.numeric(logLik(mod23))

### Third Transition t3-t4
mod34 <- clm(resp4 ~ 1 , nominal = ~ tratamento + resp3, data=dadosrep)
beta34 <- coef(mod34)
vcov34 <- vcov(mod34)
loglik34 <- as.numeric(logLik(mod34))

### Fourth Transition t4-t5
mod45 <- clm(resp5 ~ 1 , nominal = ~ tratamento + resp4, data=dadosrep)
beta45 <- coef(mod45)
vcov45 <- vcov(mod45)
loglik45 <- as.numeric(logLik(mod45))

################################################################################
### Transition probabilities for each fitted model.

dropn <- with(dadoslong,expand.grid(tratamento=levels(tratamento),
                                    respp=levels(respp)))
preds <- cbind(dropn,predict(modE,newdata=dropn)$fit)
predsA <- preds[seq(1,5,2),] ### Transition matrix for sp0;
predsP <- preds[seq(2,6,2),] ### Transition matrix for sp1.

gridpred12 <- with(dadosrep,expand.grid(tratamento=levels(tratamento),resp1=levels(resp1)))
gridpred23 <- with(dadosrep,expand.grid(tratamento=levels(tratamento),resp2=levels(resp2)))
gridpred34 <- with(dadosrep,expand.grid(tratamento=levels(tratamento),resp3=levels(resp3)))
gridpred45 <- with(dadosrep,expand.grid(tratamento=levels(tratamento),resp4=levels(resp4)))

preds12 <- cbind(gridpred23,predict(mod12,newdata=gridpred12)$fit)
preds23 <- cbind(gridpred23,predict(mod23,newdata=gridpred23)$fit)
preds34 <- cbind(gridpred34,predict(mod34,newdata=gridpred34)$fit)
preds45 <- cbind(gridpred34,predict(mod45,newdata=gridpred45)$fit)


################################################################################
### Codes for Anderson Goodman test statistic calculation

listadrogaA <- list()
listadrogaP <- list()
for(i in 1:6){
  ifelse(i%%2==1,
         listadrogaA[[(i+1)/2]]<-rbind(preds12[i,3:5],preds23[i,3:5],
                                       preds34[i,3:5],preds45[i,3:5]),
         listadrogaP[[i/2]]<-rbind(preds12[i,3:5],preds23[i,3:5],preds34[i,3:5],
                                   preds45[i,3:5]
         )
  )}
names(listadrogaA)<- c( "Fixed state 1", "Fixed state 2",
                        "Fixed state 3")
names(listadrogaP)<- c( "Fixed state 1", "Fixed state 2",
                        "Fixed state 3")

for(j in 1:3){
  row.names(listadrogaA[[j]])<- c("Time 1","Time 2","Time 3","Time 4")
  row.names(listadrogaP[[j]])<- c("Time 1","Time 2","Time 3","Time 4")
}

# Transtion probabilities matrix
dropn <- with(dadoslong,expand.grid(tratamento=levels(tratamento),
                                    respp=levels(respp)))
preds <- cbind(dropn,predict(modE,newdata=dropn)$fit)
predsA <- preds[seq(1,5,2),] ### Transition matrix for drug A;
predsP <- preds[seq(2,6,2),] ### Transition matrix for drug P.
ns <- with(dadosrep[dadosrep$tratamento=="sp0",],c(table(resp1),table(resp2),
                                                   table(resp3),table(resp4)))
ns <- matrix(ns,ncol=3,byrow=T)
colnames(ns) <- c("1","2","3")
row.names(ns) <- c("t1","t2","t3","t4")
ns2 <- with(dadosrep[dadosrep$tratamento=="sp1",],c(table(resp1),table(resp2),
                                                    table(resp3),table(resp4)))
ns2 <- matrix(ns2,ncol=3,byrow=T)
colnames(ns2) <- c("1","2","3")
row.names(ns2) <- c("t1","t2","t3","t4")
compsomaA <- list()
compsomaP <- list()
for(i in 1:3){
  difa <- matrix(unlist(apply(listadrogaA[[i]],1,"-",
                              predsA[i,3:5])),nrow=4,byrow=T)
  difp <- matrix(unlist(apply(listadrogaP[[i]],1,"-",
                              predsP[i,3:5])),nrow=4,byrow=T)
  difa2 <- difa**2
  difp2 <- difp**2
  compsomaA[[i]] <- ns[,i]*matrix(unlist(apply(difa2,1,"/",
                                               predsA[i,3:5])),nrow=4,byrow=T)
  compsomaP[[i]] <- ns2[,i]*matrix(unlist(apply(difp2,1,"/",
                                                predsP[i,3:5])),nrow=4,byrow=T)
}

################################################################################
### Tests to assess stationarity

### Anderson--Goodman test  

Xi <- sum(unlist(compsomaA))+sum(unlist(compsomaP));round(Xi,4)
pchisq(Xi, 32-8, lower.tail = FALSE)


### Wald test statistic   

beta <- c(beta12,beta23,beta34,beta45)
varcov <- as.matrix(bdiag(vcov12,vcov23,vcov34,vcov45))
beta0vec <- rep(beta0,4)
wald <- t(beta-beta0vec)%*%solve(varcov)%*%(beta-beta0vec);round(wald,4)
pchisq(wald, 32-8, lower.tail = FALSE)


### LR test 

logver <- (-2*(loglik0-sum(loglik12,loglik23,loglik34,loglik45)));round(logver,4)
pchisq(logver, 32-8, lower.tail = FALSE)


### Gradient test 

getFittedC <- ordinal:::getFittedC

clm.nll <- function(rho, par) {
  if(!missing(par)) rho$par <- par
  with(rho, {
    if(k > 0)
      sigma <- Soff * exp(drop(S %*% par[n.psi + 1:k]))
    ### NOTE: we have to divide by sigma even if k=0 since there may be an
    ### offset but no predictors in the scale model:
    eta1 <- (drop(B1 %*% par[1:n.psi]) + o1)/sigma
    eta2 <- (drop(B2 %*% par[1:n.psi]) + o2)/sigma
  })
  ### NOTE: getFitted is not found from within rho, so we have to
  ### evalueate it outside of rho
  rho$fitted <- getFittedC(rho$eta1, rho$eta2, rho$link, rho$par[length(rho$par)])
  if(all(is.finite(rho$fitted)) && all(rho$fitted > 0))
    ### NOTE: Need test here because some fitted <= 0 if thresholds are
    ### not ordered increasingly.
    -sum(rho$wts * log(rho$fitted))
  else Inf
}

clm.grad <- function(rho) {
  ### requires that clm.nll has been called prior to
  ### clm.grad.
  with(rho, {
    p1 <- if(!nlambda) dfun(eta1) else dfun(eta1, lambda)
    p2 <- if(!nlambda) dfun(eta2) else dfun(eta2, lambda)
    wtpr <- wts/fitted
    C2 <- B1*p1/sigma - B2*p2/sigma
    if(k <= 0) return(-crossprod(C2, wtpr))
    C3 <- -(eta1 * p1 - eta2 * p2) * S
    return(-crossprod(cbind(C2, C3), wtpr))
    ### NOTE: C2 and C3 are used by clm.hess
  })
}

clm.grad_direct <- function(rho, par) {
  ### does not require that clm.nll has been called prior to
  ### clm.grad.
  clm.nll(rho, par)
  clm.grad(rho)
}

# Defining the models

modEe<- update(modE, doFit=FALSE)
estimatE <- coef(modE) 

mod12a <- update(mod12, doFit=FALSE)
mod23a <- update(mod23, doFit=FALSE)
mod34a <- update(mod34, doFit=FALSE)
mod45a <- update(mod45, doFit=FALSE)

### Gradient function using liklihood estimates by clm

vetescore <- c(clm.grad_direct(mod12a, estimatE),
               clm.grad_direct(mod23a, estimatE),
               clm.grad_direct(mod34a, estimatE),
               clm.grad_direct(mod45a, estimatE))

grad <- abs(vetescore%*%(beta-beta0vec));round(grad,4)
pchisq(grad, 32-8, lower.tail = FALSE)
