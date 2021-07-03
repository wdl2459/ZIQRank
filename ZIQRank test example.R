
## install the following packages if needed
# install.packages("quantreg")
# install.packages("MASS")
# install.packages("fastDummies")


##### zero-inflated quantile rank-score based test (ZIQRank) #####

library(quantreg)
library(MASS)
library(fastDummies)



### marginal tests ###

ZIQRank <- function(formula.logistic, formula.quantile, C, y_CorD="C", data, taus=c(0.1, 0.25, 0.5, 0.75, 0.9), seed=2020){
  
  ## formulas
  
  # logistic model
  
  # arrange logistic model
  mf.logistic = model.frame(formula.logistic, data=data)
  y = model.response(mf.logistic, "numeric")
  b = 1*(y > 0)
  formula.logistic = update(formula.logistic, b ~ .)
  data.logistic = cbind(data, b)
  mf.logistic = model.frame(formula.logistic, data=data.logistic)
  
  # locate C in logistic model
  namesx = all.vars(formula.logistic)[-1]
  condition.loc = which(namesx %in% C) 
  
  # arrange logistic null model if C is not a single covariate
  if (length(condition.loc) > 1){
    mul.logistic = T
    namesx.null = setdiff(namesx, C)
    if (length(namesx.null) == 0){
      formula.logistic.null = as.formula( "b ~ 1" )
    } else formula.logistic.null = as.formula( paste( "b ~", paste(namesx.null, collapse = "+") ) )
    mf.logistic.null = model.frame(formula.logistic.null, data=data.logistic)
  } else mul.logistic = F
  
  
  # quantile model
  
  # determine elements in quantile model, create the "positive subset"
  namey = all.vars(formula.quantile)[1]
  namesx = all.vars(formula.quantile)[-1]
  namesx.score = setdiff(namesx, C)
  data.quantile = data[b==1, ]
  if (y_CorD == "D"){ # perturbation if response is count
    set.seed(seed)
    data.quantile[, namey] = dither(data.quantile[, namey], type = "right", value = 1)
  } 
  
  # extract C
  formula.quantile = as.formula( paste( namey, "~", paste(C, collapse = "+") ) )
  mf.quantile = model.frame(formula.quantile, data=data.quantile)
  c = model.matrix(attr(mf.quantile, "terms"), data=mf.quantile)[, -1]
  if (is.null(dim(c))){ # determine whether C is a single covariate, also continuous or binary
    single_CorB=T
  } else single_CorB=F
  
  # arrange quantile null model, and extract Z
  if (length(namesx.score) == 0){
    formula.quantile = as.formula( paste( namey, "~ 1" ) )
  } else formula.quantile = as.formula( paste( namey, "~", paste(namesx.score, collapse = "+") ) )
  mf.quantile = model.frame(formula.quantile, data=data.quantile)
  z = model.matrix(attr(mf.quantile, "terms"), data=mf.quantile)
  
  
  
  ## set up parameters
  m = length(y) # total sample size
  width = length(taus) # size of tau
  zerorate = length(which(b == 0)) / m # rate of 0's
  
  
  
  ## compute p-values from the marginal tests
  
  if (single_CorB == T){ # when C is a single covariate, either continuous or binary
    
    # logistic, wald test
    mod.logistic = glm(mf.logistic, family=binomial(link = 'logit'))
    pvalue.logistic = summary(mod.logistic)$coef[condition.loc+1, 4]
    
    # estimate quantiles of y|y>0 | H0
    rq0 = rq(mf.quantile, tau=taus)
    qpred0 = predict(rq0)
    
    # project C on the space of intercept and Z
    C.star = c - z %*% solve( (t(z) %*% z) ) %*% t(z) %*% c
    
    # compute the rank-score test stats, and its covariance matrix
    RS = unlist( lapply(1:width, function(kk){ sum( (taus[kk] - (data.quantile[, namey] < as.matrix(qpred0, ncol=width)[, kk]))*C.star ) / sqrt(m) }) )
    
    if (width == 1){
      cov.RS = taus*(1 - taus)
    } else {
      cov.RS = matrix(0, ncol=width, nrow=width)
      for (kk in 1:(width-1)){
        for (ll in (kk+1):width){
          cov.RS[kk, ll] = min(taus[kk], taus[ll]) - taus[kk]*taus[ll]
        }
      }
      cov.RS = cov.RS + t(cov.RS) + diag(taus*(1 - taus))
    }
    
    Sigma.hat = cov.RS * sum( C.star^2 ) / m
    if (width == 1){
      sigma.hat = sqrt( Sigma.hat )
    } else {
      sigma.hat = sqrt( diag(Sigma.hat) )
    }
    
    # marginal p-value in quantile regression 
    pvalue.quantile = 2*( 1 - pnorm( abs( RS / sigma.hat ) ) ) 
    
  } else { # when C is a set of covariates, or a single covariate with multiple categories
    
    # logistic, score test
    if (mul.logistic != T){
      mod.logistic = glm(mf.logistic, family=binomial(link = 'logit'))
      pvalue.logistic = anova(mod.logistic, test="Rao")$`Pr(>Chi)`[condition.loc+1] 
    } else {
      mod.logistic = glm(mf.logistic, family=binomial(link = 'logit'))
      mod.logistic.null = glm(mf.logistic.null, family=binomial(link = 'logit'))
      pvalue.logistic = anova(mod.logistic.null, mod.logistic, test="Rao")$`Pr(>Chi)`[2]
    }
    
    # estimate quantiles of y|y>0 | H0
    rq0 = rq(mf.quantile, tau=taus)
    qpred0 = predict(rq0)
    
    # project C on the space of intercept and Z
    C.star = c - z %*% solve( (t(z) %*% z) ) %*% t(z) %*% c
    
    # compute the rank-score test stats, and its covariance matrix
    RS = lapply(1:width, function(kk){ apply( (taus[kk] - (data.quantile[, namey] < as.matrix(qpred0, ncol=width)[, kk]))*C.star, 2, sum ) / sqrt(m) })
    
    df = ncol(C.star)
    tmp = t(C.star) %*% C.star / m
    
    var.RS = NULL
    for (kk in 1:width){
      var.RS[[kk]] = taus[kk] * (1 - taus[kk]) * tmp
    }
    
    Sigma.hat  = NULL
    for (kk in 1:width){
      temp = NULL
      for (ll in 1:width){
        temp = cbind( temp, ( min(taus[kk], taus[ll]) - taus[kk]*taus[ll] )*tmp )
      }
      Sigma.hat = rbind(Sigma.hat, temp)
    }
    
    # marginal p-value in quantile regression 
    pvalue.quantile = NULL
    for (kk in 1:width){
      stat = t( RS[[kk]] ) %*% solve(var.RS[[kk]]) %*% RS[[kk]]
      pvalue.quantile[kk] = 1 - pchisq(stat, df=df)
    }
    
  }
  
  return(list(pvalue.logistic=pvalue.logistic, pvalue.quantile=pvalue.quantile, Sigma.hat=Sigma.hat, zerorate=zerorate, taus=taus))
  
}



### tests combination ###

Combination <- function(input, method="MinP", taus=c(0.1, 0.25, 0.5, 0.75, 0.9), M=10000){
  
  ## check
  
  # choose either MinP or Cauchy
  if (method != "MinP" & method != "Cauchy"){
    stop("Please choose 'MinP' or 'Cauchy', no other options.")
  }
  
  # taus and ind should match
  if (!all(taus %in% input$taus)){
    stop("taus should be a subset of that taus used to produce input.")
  }
  
  
  ## whether from single_CorB=T or not
  if (!is.null( ncol(input$Sigma.hat) )){
    if (length(input$pvalue.quantile) != ncol(input$Sigma.hat)){
      single_CorB = F
      df = ncol(input$Sigma.hat) / length(input$pvalue.quantile)
    } else single_CorB = T
  } else single_CorB = T
  
  
  ## compute the aggregate p-value
  width = length(taus)
  
  ind = match(taus, input$taus)
  pvalue.quantile = input$pvalue.quantile[ind]
  
  
  if (method == "MinP"){
    
    # t.obs in the minp test
    t.obs = min(input$pvalue.logistic, pvalue.quantile)
    
    if (single_CorB != F){
      
      if (is.null( ncol(input$Sigma.hat) )){
        Sigma.hat = input$Sigma.hat[ind]
      } else Sigma.hat = input$Sigma.hat[ind, ind]
      
      # the (1 - t.obs/2)th percentile of the statistics in quantile regression, normal
      if (is.null( ncol(input$Sigma.hat) )){
        sigma.hat = sqrt( Sigma.hat )
      } else sigma.hat = sqrt( diag(Sigma.hat) )
      qmin.quantile = qnorm((1-t.obs/2), mean=0, sd=sigma.hat)
      
      # MC, to estimate the probability that the absolute of joint statistics in quantile regression < each threshold
      beta.sim = mvrnorm(n=M, mu=rep(0, width), Sigma=Sigma.hat)
      prob.quantile = mean( apply(beta.sim, 1, function(z){ all( abs(z) < qmin.quantile ) }) )
      
    } else {
      
      index = unlist( lapply(ind, function(kk){ ((kk-1)*df+1):(kk*df) }) ) 
      Sigma.hat = input$Sigma.hat[index, index]
      
      # the (1 - t.obs)th percentile of the statistics in quantile regression, chisq
      qmin.quantile = rep(qchisq(1-t.obs, df=df), width)
      
      # MC, to estimate the probability that the joint statistics in quantile regression < each threshold
      beta.sim = mvrnorm(n=M, mu=rep(0, width*df), Sigma=Sigma.hat)
      
      prob.quantile = mean( apply(beta.sim, 1, function(z){ 
        obs.sim = NULL
        for (kk in 1:width){
          obs.sim[kk] = t( z[((kk-1)*df+1):(kk*df)] ) %*% solve(Sigma.hat[((kk-1)*df+1):(kk*df), ((kk-1)*df+1):(kk*df)]) %*% z[((kk-1)*df+1):(kk*df)]
        }
        all( obs.sim < qmin.quantile )
      }) )
      
    }
    
    # final p-value
    pvalue = 1 - (1 - t.obs) * prob.quantile
    
  } else {
    
    # weights for quantile tests
    w = taus*(taus <= 0.5) + (1-taus)*(taus > 0.5)
    w = w / sum(w) * (1 - input$zerorate)
    
    stats.cauchy = input$zerorate * tan( (0.5-input$pvalue.logistic)*pi ) + sum ( w * tan( (0.5-pvalue.quantile)*pi ) )
    
    # final p-value
    pvalue = 1 - pcauchy(stats.cauchy)
    
  }
  
  return(pvalue)
  
}





### example in real application


# load the 100th gene in the pre-processed GSE84465 and cell-level covariates
load("example_scRNA_data.Rdata")

# analyze the 100th gene
result = ZIQRank(formula.logistic = expression ~ condition + patient,
                 formula.quantile = expression ~ condition + patient,
                 C = "condition", 
                 y_CorD = "C", 
                 data = Dat,
                 taus=seq(0.05, 0.95, by=0.05))

Combination(result, taus=seq(0.05, 0.95, by=0.05))
Combination(result, method="Cauchy", taus=seq(0.05, 0.95, by=0.05))  
  
  