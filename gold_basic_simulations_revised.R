library(MASS)
library(Iso)
library(mvtnorm)
library(optimParallel)
eta.f <- function(x){ exp(x)/(1+rowSums(exp(x))) }      
no_gold_standard_simu<- function(eta.r, theta.r, b.0.r, o.r, D, n_class=2, size=1000){
  iter = 1
  res = c()
  res[iter] = 1000 #
  # D$d_cat = rep(1, 1000)
  while((res[iter]>0.05)&(iter<100)){
    print(iter)
    par_old = c(eta.r, theta.r)
    ## now calculate the disease prevelance
    # alpha <- alpha[,1:(n_class-1),drop = FALSE]
    Zdesign <- model.matrix(formula(paste("~", colnames(D$Z))),D$Z)
    # now calculate the T = Y in this case, since the transformation is identity transformation
    x=formula(paste("~", colnames(D$X), "+", 'x1*D'))
    XDdesign0 <- model.matrix(x,cbind(D$X, data.frame(D=rep(0, size))))
    XDdesign1 <- model.matrix(x,cbind(D$X, data.frame(D=rep(1, size))))
    
    pi1.r = eta.f(Zdesign%*%alpha)
    pi0.r = 1-pi1.r
    #print(g.r)
    p0.r = dmvnorm(D$g.r-XDdesign0%*%(b.0.r), mean = c(0,0), sigma = o.r)
    p1.r = dmvnorm(D$g.r-XDdesign1%*%(b.0.r), mean = c(0,0), sigma = o.r)
    
    # p1.r = dmvnorm(g.r-as.matrix(D$X)%*%b.1.r,mean = c(0,0), sigma = o.r)
    P.r = (pi0.r*p0.r)/(pi1.r*p1.r+pi0.r*p0.r)
    # res_eta = mean(abs(D$d_prob - P.r))
    # p0.r = dmvnorm(g.r-as.matrix(D$X)%*%b.0.r,mean = c(0,0), sigma = o.r)
    # p1.r = dmvnorm(g.r-as.matrix(D$X)%*%b.1.r,mean = c(0,0), sigma = o.r)
    # P.r = (pi0.r*p0.r)/(pi1.r*p1.r+pi0.r*p0.r)
    
    # }
    
    
    # label.new = ifelse(P.r>1-P.r,1,0)
    # mean.0 = apply((g.r-as.matrix(D$X))[label.new==0,],2,mean)
    # mean.1 = apply((g.r-as.matrix(D$X))[label.new==1,],2,mean)
    #
    # P.r.cat = as.factor(ifelse(P.r>0.5,1,0))
    eta.r<- glm(P.r~., data = D$Z,family = binomial)$coef
    print(eta.r)
    # sprintf("Estimated eta: %s", eta.r)
    # sprintf("True Estimated eta: %s", eta)
    
    mytheta = function(par){
      o = matrix(0,2,2)
      o[1,1]<-par[1]
      o[2,2]<-par[2]
      o[1,2] <- par[3]*prod(sqrt(diag(o)))
      o[2,1] <- o[1,2]
      
      beta1 <- par[4:7]
      beta2 <- par[8:11]
      b0 <- matrix(c(beta1,beta2),nrow=4)
      
      # beta1 <- par[8:9]
      # beta2 <- par[10:11]
      # b1 <- matrix(c(beta1,beta2),2)
      
      
      o.det = prod(par[1],par[2])-o[1,2]^2
      o.det = det(o)
      o.inv = matrix(0,2,2)
      o.inv[1,1]<-par[2]/o.det
      o.inv[2,2]<-par[1]/o.det
      o.inv[1,2] <- -o[1,2]/o.det
      o.inv[2,1] <- -o[2,1]/o.det
      
      XDdesign0 <- model.matrix(x,cbind(D$X, data.frame(D=rep(0, size))))
      XDdesign1 <- model.matrix(x,cbind(D$X, data.frame(D=rep(1, size))))
      
      # P.e.0 = as.matrix(cbind(rep(1, 1000), rep(0, 1000), D$X, 0), ncol=4)
      # P.e.1 = as.matrix(cbind(rep(1, 1000), rep(1, 1000), D$X, D$X), ncol=4)
      sum(
        P.r*-0.5*c(apply(as.matrix(g.r)-XDdesign1%*%(b0),1,function(x) t(as.matrix(x))%*%o.inv%*%as.matrix(x)))
        ,
        (1-P.r)*-0.5*c(apply(as.matrix(g.r)-XDdesign0%*%(b0),1,function(x) t(as.matrix(x))%*%o.inv%*%as.matrix(x)))
        ,-N/2*log(o.det))
    }
    theta.opt = try(optim(theta.r,mytheta, method = 'BFGS',
                       control= list(fnscale=-1)))
    
    
    if(inherits(theta.opt, "try-error")==TRUE){
      res = res
      res = append(res, 1000, after = length(res))
    }else{
      theta.r <- theta.opt$par
      theta.conv <- theta.opt$convergence
      b.0.r = matrix(theta.r[4:11],4,2)
      # b.1.r = matrix(theta.r[8:11],2,2)
      print(b.0.r)
      # print(b.1.r)
      
      o.r = matrix(0,2,2)
      o.r[1,1]<-theta.r[1]
      o.r[2,2]<-theta.r[2]
      o.r[1,2] <- theta.r[3]*prod(sqrt(diag(o.r)))
      o.r[2,1] <- o.r[1,2]
      print(o.r)
      o.det = prod(diag(o.r))-o.r[1,2]^2
      o.inv = matrix(0,2,2)
      o.inv[1,1]<-o.r[2,2]/o.det
      o.inv[2,2]<-o.r[1,1]/o.det
      o.inv[1,2] <- -o.r[2,1]/o.det
      o.inv[2,1] <- -o.r[2,1]/o.det
      
      # g.r.1 <- matrix(0, nrow = 1000, ncol = 2)
      # g.r.1[,1] <- pava(diag(as.matrix(D$X)%*%as.matrix(apply(P.r,1, function(x) (1-x)*as.matrix(b.0.r[,1],nrow=2))+apply(P.r,1, function(x) x*as.matrix(b.1.r[,1],nrow=2))))
      #                   +0.5*(diag(as.matrix(D$X)%*%as.matrix(apply(P.r,1, function(x) (1-x)*as.matrix(b.0.r[,2],nrow=2))+apply(P.r,1, function(x) x*as.matrix(b.1.r[,2],nrow=2))))
      #                   -g.r[,2])*o.inv[1,2]/o.inv[1,1])
      # 
      # g.r.1[,2] <-  pava(diag(as.matrix(D$X)%*%as.matrix(apply(1-P.r,1, function(x) (1-x)*as.matrix(b.0.r[,2],nrow=2))+apply(P.r,1, function(x) x*as.matrix(b.1.r[,2],nrow=2))))
      # +0.5*(diag(as.matrix(D$X)%*%as.matrix(apply(P.r,1, function(x) (1-x)*as.matrix(b.0.r[,1],nrow=2))+apply(P.r,1, function(x) x*as.matrix(b.1.r[,1],nrow=2))))
      #       -g.r[,1])*o.inv[1,2]/o.inv[2,2])
      # 
      # res2 = mean(abs(g.r-g.r.1))
      # g.r <- g.r.1
      # res[iter] = mean(abs(g.r-D$Y))
      par_new = c(eta.r, theta.r)
      print(par_new)
      
      # par_new = c(res_eta, theta.r)
      res = append(res, max(abs(par_old-par_new)), after = length(res))
      # res = append(res,max(res2,abs(par_old-par_new)), after = length(res))
      # g.r <- g.r.1
      # print(res[iter])

    }
    print(o.r)
    iter = iter+1
  }
  
  
  
  # mean(as.matrix(D$X[P.r<0.5,])%*%B.0 - D$Y[P.r<0.5,])
  
  
  print(paste("Converge in %s iteration:", iter))
  print("Converge covariance matrix: ")
  print(o.r)
  
  print("Original covariance matrix: ")
  print(O)
  print(list(eta.r = eta.r, theta.r = theta.r))
  return(list(eta.r = eta.r, theta.r = theta.r))
}


eta.out = data.frame()
theta.out = data.frame()
var=c(0.1, 0.5)
rho=0.1
O = abs(diag(var,2))
O[1,2] <- rho*sqrt(var[1])*sqrt(var[2])
O[2,1] <- O[1,2]
simu_data <- function(size=1000, k_biomarkers=2, m_covariates=1, z_covariates=1, n_class=2,
                      alpha=cbind(c(0.2,0.067)), #z_covariates + 1
                      B=matrix(c(c(0.1,0.3,1,2),c(0.2,0.4,1,2)),nrow=4), eta=0.5,
                     seed=123, O=as.matrix(cbind(c(0.1, 0.02), c(0.02, 0.5)), nrow=2)){
  set.seed(seed)
  k <-k_biomarkers; #number of biomarker
  m <-m_covariates; #covariates dimension
  d <-n_class; #number of class(0-no disease 1-disease)
  N<- size #number of observation: what if we have more observations, will we have more accurate?
  ## B:  what's B, B is the covariate for each disease status. 0 as no 1 as yes
  
  ## X
  X = list()
  for (i in 1:m_covariates) {
    set.seed(seed+i)
    X = append(X, list(rnorm(size,0,1)))
    names(X) = c(names(X), paste0("x", i))
  }
  
  X_dataframe = as.data.frame(do.call(cbind, X))
  
  ## # Z: covariates for disease 
  Z = list()
  for (i in 1:m_covariates) {
    set.seed(seed*i)
    Z = append(Z, list(rnorm(size,0,1)))
    names(Z) = c(names(Z), paste0("z", i))
  }
  
  Z_dataframe = as.data.frame(do.call(cbind, Z))
  
  ## now calculate the disease prevelance
  alpha <- alpha[,1:(n_class-1),drop = FALSE]
  Zdesign <- model.matrix(formula(paste("~", colnames(Z_dataframe))),Z_dataframe)
  # exp(as.matrix(Z_dataframe%*%alpha)/(1+exp(as.matrix(Z_dataframe)%*%eta))
  # d_cat = rbinom(N,1,d_prob)
  eta.f <- function(x){ exp(x)/(1+rowSums(exp(x))) }                                       
  Pd    <- eta.f(Zdesign%*%alpha)                          # N by D
  YPd   <- cbind(1-rowSums(Pd), Pd)                        # N by D+1
  D     <- as.factor(apply(YPd,1,function(p){sum(rmultinom(1, size=1, prob=p)*(0:1))}))
  
  # now calculate the T = Y in this case, since the transformation is identity transformation
  x=formula(paste("~", colnames(X_dataframe), "+", 'x1*D'))
  XDdesign <- model.matrix(x,cbind(X_dataframe, D))
  
  ## error
  ## Control covariance matrix
  eps <- mvrnorm(n = N, c(0,0), O)
  Y = XDdesign%*%(B) + eps
  g.r = Y
  
  ## Data
  data_output <- list()
  data_output$X <- X_dataframe
  data_output$Z <- Z_dataframe
  data_output$D <-  D
  data_output$g.r <- Y
  return(data_output)
}
train = simu_data()
### Initial values:
eta.r <- cbind(c(0,0))
theta.r = runif(3,0,1)
o.r = abs(diag(theta.r[c(1,2)],2))
o.r[1,2] <- theta.r[3]*sqrt(theta.r[1])*sqrt(theta.r[2])
o.r[2,1] <- o.r[1,2]

# b.r <- matrix(rnorm(4),nrow = 2)
b.0.r <- matrix(rnorm(8),nrow = 4)
# b.1.r <- matrix(rnorm(4),nrow = 2)
theta.r = c(theta.r, c(b.0.r))
simu1_result = no_gold_standard_simu(eta.r, theta.r, b.0.r, o.r, train)
eta.out = rbind(eta.out, simu1_result$eta.r)
theta.out = rbind(theta.out, simu1_result$theta.r)
# 
# ### simulation starts here: 
# gen.latent.risk <- function(Z,alpha,latent.class.names = NULL){
#   nclass <- ncol(alpha) + 1
#   YPd   <- ctg.eta.fun(Zdesign %*% alpha)
#   try(colnames(YPd) <- c("ref",colnames(alpha)),silent = TRUE)
#   try(colnames(YPd) <- latent.class.names,silent = TRUE)
#   attr(YPd,'latent.class') <- as.factor(apply(YPd,1,function(p){sum(rmultinom(1, size=1, prob=YPd[1,])*(0:(nclass - 1)))}))
#   YPd
# }
# 
# ctg.eta.fun  <- function(x){ 
#   eta <- exp(x)/(1+rowSums(exp(x))) 
#   eta <- cbind(1-rowSums(eta),eta) 
#   ix <- which(apply(eta,1,function(x)any(is.na(x))))
#   if(length(ix)){
#     eta[ix,] <- ctg.eta.fun.careful(eta[ix,,drop = FALSE])
#   }
#   eta
# } 
# 
# for(s in seq(51,75)){
#   set.seed(s+100)
#   k <-2; #number of biomarker
#   m <-2; #covariates dimension
#   d <-2; #number of class(0-no disease 1-disease)
#   N<- 2000 #number of observation: what if we have more observations, will we have more accurate?
#   
#   ## B:  what's B, B is the covariate for each disease status. 0 as no 1 as yes
#   # beta1 <- c(0.1,0.3,1,2)
#   # beta2 <- c(0.2,0.4,1,2)
#   # B.0 <- matrix(c(beta1,beta2),nrow=4)
#   # 
#   # # beta1 <- c(-0.1,0.2,1,2)
#   # # beta2 <- c(-0.2,0.3,1,2)
#   # # B.1 <- matrix(c(beta1,beta2),2)
#   # 
#   # # B.1 = B.0
#   # 
#   # ## Delta
#   # # delta <- sample(c(0,1),N, replace = T)
#   
#   
#   ## X
#   x1 = rnorm(N,0,1)
#   # x2 = rnorm(N,0,2)
#   X = data.frame(x1 = x1)
#                  # , x2 = x2) # what's X, the values of the biomarkers. 
#   
#   ## g.y
#   eta<- c(0.5) #what's eta: the covariates for disease prevelance
#   d_prob = exp(as.matrix(X)%*%eta)/(1+exp(as.matrix(X)%*%eta))
#   d_cat = rbinom(N,1,d_prob)
#   
#   d_prob.0 = 1-d_prob[d_cat==0]
#   d_prob.1 = d_prob[d_cat==1]
#   
#   
#   X.0 = X[d_cat==0,]
#   X.1 = X[d_cat==1,]
#   
#   
#   ## D: the covarianc matrix of the residual
#   D <- list()
#   D$X <- X
#   # rbind(X.0,X.1)
#   d_cat <- c(rep(0,length(X.0)),rep(1,length(X.1)))
#   D$d_cat <-  d_cat
#   D$d_prob <- c(d_prob.0, d_prob.1)
#   
#   
#   var <- c(0.3,0.3)
#   O = abs(diag(var,2))
#   rho = 0.1
#   O[1,2] <- rho*sqrt(var[1])*sqrt(var[2])
#   O[2,1] <- O[1,2]
#   
#   ## Control covariance matrix
#   rho.r <- rho
#   var1 <- var[1]
#   var2 <- var[2]
#   
#   ## g.y
#   eps <- mvrnorm(n = N, c(0,0), O)
#   eps.0 = eps[d_cat==0]
#   eps.1 = eps[d_cat==1]
#   # 
#   Y = as.matrix(cbind(rep(1, N), X, d_cat, X*d_cat), ncol=4)%*%(B.0) + eps
#   D$Y <- as.matrix(cbind(rep(1, N), X, d_cat, X*d_cat), ncol=4)%*%(B.0) + eps
#   g.r = D$Y
#   D1 = D
#   for (name in names(D)) {
#     D1[name] = D[name][1:N/2, ]
#   }
#     
#   
#   



}
  



### now let's fit with the lag
library(latentreg) ## have to install gfortran-6.1.pkg first for Mac
# By default 'dgp()' below assumes 3 latent class statuses {0,1,2} .
loadNamespace("latentreg")
set.seed(123)
train<-latentreg:::dgp(N=1000, seed=123, nclass = 2, alpha = matrix(c(-0.5), 1), 
                       sigma = c(1,1)/2, # continuous test std.dev or sqrt(covariance matrix) 
                       sigma.x = 1, # covariate for continuous test 
                       lambda = c(1, 0.5),
                       beta = cbind(
                         c(+5.0, -1.0, +1.0, +2.5),
                         c(+5.0, -1.0, +0.5, +1.0)
                       ),
                       gamma  = list(
                         gamma1 = cbind( c(-3, -2, +5, +4),
                                         c(-4, -2, +4, +6))
                       )
                       )
train$X1 = D$X$x1
train$X2 = D$X$x2
train$Z = D$X$x1
train$Z2 = D$X$x1
train$Z3 = D$X$x2

test$Z = test$X1
test$Z2 = D$X$x1
test$Z3 = D$X$x2


test<-latentreg:::dgp(N=1000, seed=2345, , nclass = 2, alpha = matrix(c(-0.5), 1), 
                      sigma = c(1,1)/2, # continuous test std.dev or sqrt(covariance matrix) 
                      sigma.x = 1, # covariate for continuous test 
                      lambda = c(1, 0.5),
                      beta = cbind(
                        c(+5.0, -1.0, +1.0, +2.5),
                        c(+5.0, -1.0, +0.5, +1.0)
                      ),
                      gamma  = list(
                        gamma1 = cbind( c(-3, -2, +5, +4),
                                        c(-4, -2, +4, +6))
                      ))
tests.cnt <- attr(train,"continuous.model")
tests.ctg <- attr(train,"categorical.model")
prev <- attr(train,"prevalence.model")
true.param <- attr(train,"param")
## Not run:
# univariate optimization via estimateTransform in package 'car' (Box-Cox)
fit.car.bc <- lvreg(train,prev,tests.cnt,nclass = 2,verbose=2,
                    control.cnt = list(how = 'et'))
coef(fit.car.bc);coef(fit.car.bc,stratify = TRUE);summary(fit.car.bc)
test_results = cbind(apply(predict(fit.car.bc,test),1,which.max),test$'Unknown D')








# univariate optimization via estimateTransform in package 'car' (Yeo-Johnson)
fit.car.yj <- lvreg(train,prev,tests.cnt,nclass = 3,verbose=2,
                    control.cnt = list(how = 'et',family = "yeo.johnson"))
# 'univariate' minimization optim(L-BFGS-B)
umin <- lvreg(train,prev,tests.cnt,nclass = 3,verbose=2)
coef(umin);coef(umin,stratify = TRUE)

#14 lvreg
summary(umin)
# 'joint' minimization optim(L-BFGS-B)
jmin <- lvreg(train,prev,tests.cnt,nclass = 3,verbose=2,control.cnt = list(jointly = TRUE))
# 'univariate' solve by uniroot.extreme fails
# fast but requires memory ~ (nrow(data)*nclass x nsub) !
# where nsub = control.cnt$optim.control$n the number of subintervals for [lower,upper]
usolve.fail <- lvreg(train,prev,tests.cnt,nclass = 3,verbose=2)
# 'univariate' solve by uniroot.extreme, slower but low memory usage
usolve <- lvreg(train,prev,tests.cnt,nclass = 3,verbose=2,
                control.cnt = list(how = 'solve',memory.efficient = TRUE,optim.arg = list(control = list(n = 1000))))
# 'univariate' minimization with analytic gradient if lambda != 0
fit <- lvreg(train,prev,tests.cnt,control.cnt = list(analytic.gradient = TRUE),nclass = 3,verbose=2)
# joint minimization
fit <- lvreg(train,prev,tests.cnt,control.cnt = list(jointly = TRUE),nclass = 3,verbose=2)
fit <- lvreg(train,prev,tests.cnt,control.cnt = list(jointly = TRUE,optimizer = 'BBoptim'),nclass = 3,verbose=2)
fit <- lvreg(train,prev,tests.cnt,control.cnt = list(jointly = TRUE,optimizer = 'nlminb'),nclass = 3,verbose=2)
fit <- lvreg(train,prev,tests.cnt,control.cnt = list(jointly = TRUE,optimizer = 'JDEoptim'),nclass = 3,verbose=2)
# how = "solve" with jointly = TRUE,dfsane
jsolve <- lvreg(train,prev,tests.cnt,
                control.cnt = list(how = 'solve',jointly = TRUE, optimizer = 'dfsane',optim.arg = list(method=2)),nclass = 3,verbose=2)
jsolve2 <- lvreg(train,prev,tests.cnt,
                 control.cnt = list(how = 'solve',jointly = TRUE, optimizer = 'sane',optim.arg = list(method=2)),nclass = 3,verbose=2)
# accelerated convergence via squarem
fit <- lvreg(train,prev,tests.cnt,nclass = 2,verbose=2,accelerated = T)
# prediction
cbind(apply(predict(fit,test),1,which.max),test$'Unknown D')
# bootstrap a continuous model
fit$maxit <- 300
fit.boot <- boot.lvreg(fit,data = train,R = 5,sim = "o",verbose = 1);summary(fit.boot)
# fix beta of T1 and lambda of T2
fit.fixed <- lvreg(train,prev,tests.cnt,nclass = 3,verbose=2,
                   lambda = coef(fit)$lambda,beta = coef(fit)$beta,
                   fixed = list(beta = c(T,F,F),trans = c(F,T,F))
)
# bootstrap a mixed model
# mixed parametric model (Box-Cox + multinomial)
fit.mult.bc<-lvreg(train,prev,tests.cnt,tests.ctg,nclass = 3,accelerated = T,verbose = 2,control.sqem = list(tol = 1e-3))
fit.boot <- boot.lvreg(fit.mult.bc,data = train,R = 5,sim = "o",verbose = 1);summary(fit.boot)
# categorical - multinomial model
fit.mult<-lvreg(train,prev,ctg.fm = tests.ctg,nclass = 3,accelerated = T,verbose = 2)
# categorical - clm model

#lvreg 15
fit.clm<-lvreg(train,prev,ctg.fm = tests.ctg,nclass = 3,accelerated = T,verbose = 2,control.ctg = list(model = 'clm'))
# mixed parametric model (Box-Cox + clm)
fit.clm.bc<-lvreg(train,prev,tests.cnt,tests.ctg,nclass = 3,accelerated = T,verbose = 2,control.ctg = list(model = 'clm'))
fit.mix.fixed <- lvreg(train,prev,tests.cnt,tests.ctg,nclass = 3,verbose=2,
                       lambda = coef(fit.mix)$lambda,beta = coef(fit.mix)$beta,gamma = coef(fit.mix)$gamma,
                       fixed = list(beta = c(T,F,F),trans = c(F,T,F),gamma = c(T,T,F,T,F))
)
# mixed semiparametric estimation of transformation in C, of beta in R via lm
fit.sp <- lvreg(
  train,prev,tests.cnt,tests.ctg,nclass = 3,
  verbose = 2,maxit = 150,abstol = 1e-3,
  control.cnt = list(family = "sp",optim.arg = list(control = list(n = 10)),C = T)
)
# compare prediction with true 'D'
cbind(apply(predict(fit.sp,test$data),1,which.max),test$'unknown D')
##################################################################
## End(Not run) # end not run



