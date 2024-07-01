library(glmnet)
library(plyr)
source("utils.R")
source("data_fun.R")

# oracle lasso c(alpha)
LassoQuantile = function(X, Y, beta, beta_fun, lam1,  alpha, ite, sig)
{
  n = dim(X)[1]
  p = dim(X)[2]
  c = rep(0, ite)
  s = sum(beta != 0)
  for (i in 1:ite)
  {
    beta_prime = beta_fun()
    Y_prime = X %*% beta_prime + matrix(rnorm(n, sd=sig), ncol = 1)
    lassoModel = glmnet(X, Y_prime, family = 'gaussian', intercept = F, lambda = lam1)
    beta_hat = as.matrix(lassoModel$beta)
    c[i] = L2NormSquare(X %*% (beta_prime - beta_hat)) / s / log(p) / (sig^2)
  }
  cAlpha = quantile(c, probs=1 - alpha)
  
  lassoModel = glmnet(X, Y, family = 'gaussian', intercept = F, lambda = lam1)
  r = sig * sqrt(cAlpha * s * log(p) / n)
  rv = L2NormSquare(X %*% (beta - lassoModel$beta)) / n / (r^2)
  coverage = rv <= 1
  logVol = n * log(r)
  result = list(r = r, coverage=coverage, logVol=logVol)
  return(result)
}

oracle_test = function(n, p, sig, alpha, ite, arg)
{
  # add ite_quantile as argument eventually
  ite_quantile = 100
  r = rep(0, ite)
  coverage = rep(0, ite)
  logVol = rep(0, ite)
  
  if (arg[['lam1Type']] == 'val') {
      lam1 = rep(2 * sqrt(2) * sqrt(log(p) / n), ite)
  }
  for (i in 1:ite)
  {
    X = generateX(n, p, arg[['Xtype']])
    beta = arg[['beta_fun']]()
    Y = X %*% beta + matrix(rnorm(n, sd = sig), ncol = 1)
    # calculate necessary statistics
    tempRes = LassoQuantile(X, Y, beta, arg[['beta_fun']], lam1[i],  alpha, ite_quantile, sig)
    r[i] = tempRes$r
    coverage[i] = tempRes$coverage
    logVol[i] = tempRes$logVol
  }
  return(data.frame(r, coverage, logVol))
}

oracle_single = function(n, p, svalue, bvalue, sig, lam1Types, alpha, ite, Xtype)
{
  arg = list()
  res_list = list()
  count = 1
  for (i in 1:length(lam1Types))
  {
    for (j in 1:length((svalue)))
    {
      for (k in 1:length(bvalue))
      {
        arg[['lam1Type']] = lam1Types[i]
        arg[['Xtype']] = Xtype
        arg[['beta_fun']] = function() {
            return(generateBeta(type = 3, p, svalue[j], bvalue[k]))
        }
        # calculate necessary statistics
        tempRes = oracle_test(n, p, sig, alpha, ite, arg)
        tempRes['lam1Type'] = lam1Types[i]
        tempRes['s'] = svalue[j]
        tempRes['b'] = bvalue[k]
        tempRes['method'] = 'oracle'
        res_list[[count]] = tempRes
        count = count + 1
      }
    }
  }
  return(ldply(res_list))
}
 
oracle_multiple = function(n, p, s_s, b_s, s_w, b_w, sig, lam1Types, alpha, ite, Xtype)
{
  arg = list()
  res_list = list()
  count = 1
  for (i in 1:length(lam1Types))
  {
    for (j in 1:length((b_s)))
    {
      for (k in 1:length(b_w))
      {
        arg[['lam1Type']] = lam1Types[i]
        arg[['Xtype']] = Xtype
        arg[['beta_fun']] = function() {
            return(generateBeta(type = 4, p, s = s_s, b = b_s[j], s2 = s_w, b2 = b_w[k]))
        }
        # calculate necessary statistics
        tempRes = oracle_test(n, p, sig, alpha, ite, arg)
        tempRes['lam1Type'] = lam1Types[i]
        tempRes['bs'] = b_s[j]
        tempRes['bw'] = b_w[k]
        tempRes['method'] = 'oracle'
        res_list[[count]] = tempRes
        count = count + 1
      }
    }
  }
  return(ldply(res_list))
}

