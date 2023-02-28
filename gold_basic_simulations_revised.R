library(MASS)
library(Iso)
library(mvtnorm)
library(optimParallel)
library(latentreg)

eta.f <- function(x){ exp(x)/(1+rowSums(exp(x))) }      
no_gold_standard_simu<- function(eta.r, theta.r, b.0.r, o.r, D, n_class, size, x_T_formula){
  iter = 1
  res = c()
  res[iter] = 1000 #
  # D$d_cat = rep(1, 1000)
  while((res[iter]>0.05)&(iter<200)){
    print(iter)
    par_old = c(eta.r, theta.r)
    ## now calculate the disease prevelance
    # alpha <- alpha[,1:(n_class-1),drop = FALSE]
    Zdesign <- model.matrix(formula(paste("~", colnames(D$Z))),D$Z)
    # now calculate the T = Y in this case, since the transformation is identity transformation
    # x=formula(paste("~", colnames(D$X), "+", 'x1*D'))
    XDdesign0 <- model.matrix(x_T_formula,cbind(D$X, data.frame(D=rep(0, size))))
    XDdesign1 <- model.matrix(x_T_formula,cbind(D$X, data.frame(D=rep(1, size))))
    
    pi1.r = eta.f(Zdesign%*%eta.r)
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
      
      beta1 <- par[4:(3+nrow(b.0.r))]
      beta2 <- par[(4+nrow(b.0.r)):(3+2*nrow(b.0.r))]
      b0 <- matrix(c(beta1,beta2),nrow=nrow(b.0.r))
      
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
      
      XDdesign0 <- model.matrix(x_T_formula,cbind(D$X, data.frame(D=rep(0, size))))
      XDdesign1 <- model.matrix(x_T_formula,cbind(D$X, data.frame(D=rep(1, size))))
      
      # P.e.0 = as.matrix(cbind(rep(1, 1000), rep(0, 1000), D$X, 0), ncol=4)
      # P.e.1 = as.matrix(cbind(rep(1, 1000), rep(1, 1000), D$X, D$X), ncol=4)
      sum(
        P.r*-0.5*c(apply(as.matrix(D$g.r)-XDdesign1%*%(b0),1,function(x) t(as.matrix(x))%*%o.inv%*%as.matrix(x)))
        ,
        (1-P.r)*-0.5*c(apply(as.matrix(D$g.r)-XDdesign0%*%(b0),1,function(x) t(as.matrix(x))%*%o.inv%*%as.matrix(x)))
        ,-size/2*log(o.det))
    }
    theta.opt = try(optim(theta.r,mytheta, method = 'BFGS',
                       control= list(fnscale=-1)))
    
    
    if(inherits(theta.opt, "try-error")==TRUE){
      res = res
      res = append(res, 1000, after = length(res))
    }else{
      theta.r <- theta.opt$par
      theta.conv <- theta.opt$convergence
      b.0.r = matrix(theta.r[4:(nrow(b.0.r)*ncol(b.0.r)+3)],nrow=nrow(b.0.r))
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
simu_data <- function(seed=123, size=1000, k_biomarkers=2, m_covariates=2, z_covariates=1, n_class=2,
                      alpha=cbind(c(0.2,0.067)), 
                      B=matrix(c(c(0.1,0.3,1,2,3),c(0.2,0.4,1,2,1)),nrow=5),
                     O=as.matrix(cbind(c(0.1, 0.02), c(0.02, 0.5)), nrow=2), x_T_formula='~x1+x1*D', 
                     transform='identity'){
  print(transform)
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
    names(X)[i] = paste0("x", i)
  }
  
  X_dataframe = as.data.frame(do.call(cbind, X))
  
  ## # Z: covariates for disease 
  Z = list()
  for (i in 1:z_covariates) {
    set.seed(seed*i)
    Z = append(Z, list(rnorm(size,0,1)))
    names(Z)[i] = paste0("z", i)
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
  # x=formula(paste("~", paste(colnames(X_dataframe), collapse = '+'), "+", 'x1*D'))
  XDdesign <- model.matrix(x_T_formula,cbind(X_dataframe, D))
  
  ## error
  ## Control covariance matrix
  eps <- mvrnorm(n = N, c(0,0), O)
  g.r = XDdesign%*%(B) + eps
  if (transform=='identity'){
    Y = g.r
    }
  
  ## Data
  data_output <- list()
  data_output$X <- X_dataframe
  data_output$Z <- Z_dataframe
  data_output$D <-  D
  data_output$g.r <- g.r
  data_output$Y <- Y
  return(data_output)
}

my_pred <- function(test, response, covariates, alpha, B, O, x_T_formula, prev){
  Zdesign <- model.matrix(prev, test)
  # now calculate the T = Y in this case, since the transformation is identity transformation
  # x=formula(paste("~", colnames(D$X), "+", 'x1*D'))
  XDdesign0 <- model.matrix(x_T_formula,cbind(test[covariates], data.frame(D=rep(0, nrow(test)))))
  XDdesign1 <- model.matrix(x_T_formula,cbind(test[covariates], data.frame(D=rep(1, nrow(test)))))
  
  pi1.r = eta.f(Zdesign%*%alpha)
  pi0.r = 1-pi1.r
  #print(g.r)
  p0.r = dmvnorm(test[c(response)]-XDdesign0%*%(B), mean = c(0,0), sigma = O)
  p1.r = dmvnorm(test[response]-XDdesign1%*%(B), mean = c(0,0), sigma = O)
  
  # p1.r = dmvnorm(g.r-as.matrix(D$X)%*%b.1.r,mean = c(0,0), sigma = o.r)
  P.r = (pi0.r*p0.r)/(pi1.r*p1.r+pi0.r*p0.r)
  
  return(P.r)
}



x_T_formula=formula(paste("~", paste(c('x1', 'x2'), collapse = '+'), "+", 'x1*D'))
var=c(0.1, 0.5)
rho=0.1
O = abs(diag(var,2))
O[1,2] <- rho*sqrt(var[1])*sqrt(var[2])
O[2,1] <- O[1,2]

size=2000; k_biomarkers=2; m_covariates=2; z_covariates=1; n_class=2;
alpha=cbind(c(0.2,0.067)); #z_covariates + 1
B=matrix(c(c(0.1,0.3,1,2,3),c(0.2,0.4,1,2,1)),nrow=5); eta=0.5;
seed=123; O=as.matrix(cbind(c(0.1, 0.02), c(0.02, 0.5)), nrow=2)
train = simu_data(seed=123, size=size, k_biomarkers=k_biomarkers, m_covariates=m_covariates, z_covariates=z_covariates, n_class=2,
                  alpha=alpha, #z_covariates + 1
                  B=B, 
                  O=O,
                  x_T_formula, transform='identity')

test = simu_data(seed+1, 2000, k_biomarkers, m_covariates, z_covariates, n_class,
                  alpha,B, O,x_T_formula)

### Initial values:
eta.r <- cbind(c(0,0))
theta.r = runif(3,0,1)
o.r = abs(diag(theta.r[c(1,2)],2))
o.r[1,2] <- theta.r[3]*sqrt(theta.r[1])*sqrt(theta.r[2])
o.r[2,1] <- o.r[1,2]

# b.r <- matrix(rnorm(4),nrow = 2)
b.0.r <- matrix(rnorm(nrow(B)*ncol(B)),nrow = nrow(B))
# b.1.r <- matrix(rnorm(4),nrow = 2)
theta.r = c(theta.r, c(b.0.r))
simu1_result = {}
for (boots in c(1,100)){
  subsample = {}
  subsample_index = sample(size,size=500)
  for (item in names(train)){
    print(item)
    if(is.null(dim(train[[item]]))){
      print('yes')
      subsample[[item]] = train[[item]][subsample_index]
    }else{
      subsample[[item]] = as.matrix(train[[item]][subsample_index,])
      print(head(subsample[[item]]))
      if(!is.null(names(train[[item]]))){
        subsample[[item]] = data.frame(subsample[[item]])
        names(subsample[[item]]) = names(train[[item]])
      }
      print(head(subsample[[item]]))
    }
  }
  
  simu1_result = no_gold_standard_simu(eta.r, theta.r, b.0.r, o.r, subsample,
                                              n_class=n_class,
                                              size=500, x_T_formula=x_T_formula)
  lapply(simu1_result, write, paste0("test_", boots, ".txt"), append=TRUE, ncolumns=1000)
  fpath = paste0('Documents/no_gold_standard/', "test_", boots, ".txt")
  sink(fpath)
  #print my_list to file
  print(simu1_result)
  #close external connection to file 
  sink()
} 


lapply(simu1_result, write, "test.txt", append=TRUE, ncolumns=1000)
sink('Documents/no_gold_standard/simu1_result_revised.txt')
#print my_list to file
print(simu1_result)
#close external connection to file 
sink()


eta.out = rbind(eta.out, simu1_result$eta.r)
theta.out = rbind(theta.out, simu1_result$theta.r)
o.out = abs(diag(simu1_result$theta.r[c(1,2)],2))
o.out[1,2] <- simu1_result$theta.r[3]*sqrt(simu1_result$theta.r[1])*sqrt(simu1_result$theta.r[2])
o.out[2,1] <- o.out[1,2]

o.out[1,1]<-simu1_result$theta.r[1]
o.out[2,2]<-simu1_result$theta.r[2]
o.out[1,2] <- simu1_result$theta.r[3]*prod(sqrt(diag(o.out)))
o.out[2,1] <- o.out[1,2]


colnames(train$g.r) <- c('t1', 't2')
colnames(train$Y) <- c('y1', 'y2')
colnames(test$g.r) <- c('t1', 't2')
colnames(test$Y) <- c('y1', 'y2')

train_lvreg = cbind(train$X, train$Z, train$g.r, train$Y)
test_lvreg = cbind(test$X, test$Z, test$g.r, test$Y)
cnt.manifest <- formula(paste(paste(colnames(train$g.r), collapse = '+'), "~", paste(c('x1', 'x2'), collapse = '+'), 
                                  "+", 'x1*D'))
prev <- formula(paste("D ~", colnames(head(train$Z))))
fit.rl <- lvreg(train_lvreg,prev,cnt.manifest,nclass = 2)
y = list()
y <- data.frame(cbind(pred_lvreg=c(predict(fit.rl,test_lvreg,type = "class")), D=test$D))
y_pred_my <-my_pred(test_lvreg, response=c('t1', 't2'),covariates=c('x1', 'x2'), alpha, B, O, 
                    x_T_formula=formula(paste("~", paste(c('x1', 'x2'), collapse = '+'), "+", 'x1*D')), 
                    prev=formula(paste('~', paste(colnames(head(train$Z)), collapse = '+'))))
y_fit_my <- my_pred(train_lvreg, response=c('t1', 't2'),covariates=c('x1', 'x2'), alpha, B, O, 
                     x_T_formula=formula(paste("~", paste(c('x1', 'x2'), collapse = '+'), "+", 'x1*D')), 
                     prev=formula(paste('~', paste(colnames(head(train$Z)), collapse = '+'))))
y_fit_my = apply(y_fit_my,1, function(x) if(x>(1-x)){0}else{1})
mean(y_fit_my==train$D)

y$D_pred_my =  apply(y_pred_my,1, function(x) if(x>(1-x)){1}else{2})
mean(y$D_pred_my == y$D)
mean(y$pred_lvreg== y$D, na.rm = TRUE)
sum(y$pred_lvreg== y$D, na.rm = TRUE)


