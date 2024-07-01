library(glmnet)
library(plyr)
library(foreach)

# Robin's adaptive method
adaptive = function(X1, Y1, X2, Y2, beta, lam, sig, alpha)
{
  # calculate betaHat for first part
  lasso = glmnet(X1, Y1, family="gaussian", lambda=lam, intercept=F)
  betaHat = lasso$beta
  # calculate a, b, c for second part
  n2 = dim(X2)[1]
  a = max(0, L2NormSquare(Y2 - X2 %*% betaHat) / n2 - sig^2)
  b = 2 * sig^4 / n2
  c = 4 * sig^2 / n2
  zAlpha = qnorm(1 - alpha)
  
  # calculate the radius/bound
  delta = (2 * a + zAlpha^2 * c)^2 - 4 * (a^2 - zAlpha^2 * b)
  rSquare = ((2 * a + zAlpha^2 * c) + sqrt(delta)) / 2
  r = sqrt(rSquare)
  
  rv = L2NormSquare(X2 %*% (beta - betaHat)) / n2 / (r^2)
  coverage = rv <= 1
  logVol = n2 * log(r)
  
  result = list(r = r, coverage=coverage, logVol=logVol)
  return(result)
}

# warp a group of test
# arg: a named list contains all paths
adaptive_test = function(n, p, sig, alpha, ite, arg)
{
  r = rep(0, ite)
  coverage = rep(0, ite)
  logVol = rep(0, ite)

  beta_pool = t(as.matrix(read.table(arg[['beta']])))
  Y1_pool = t(as.matrix(read.table(arg[['Y1']])))
  Y2_pool = t(as.matrix(read.table(arg[['Y2']]))) 

  if (arg[['lam1Type']] == '1se') {
      lam1 = read.table(arg[['1se']])[,1]
  }
  else if (arg[['lam1Type']] == 'cv') {
      lam1 = read.table(arg[['cv']])[,1]
  }
  else if (arg[['lam1Type']] == 'val') {
      lam1 = rep(2 * sqrt(2) * sqrt(log(p) / n), ite)
  }

  #foreach (i = 1:ite, .export = c("adaptive", "L2NormSquare"), .packages = c("glmnet")) %dopar%
  for(i in 1:ite)
  {
    fileX1 = paste(arg[['root']], '/X1/', i, '.txt', sep = "") 
    fileX2 = paste(arg[['root']], '/X2/', i, '.txt', sep = "") 

    X1 = as.matrix(read.table(fileX1))
    X2 = as.matrix(read.table(fileX2))    
    beta = beta_pool[,i]
    Y1 = Y1_pool[,i]
    Y2 = Y2_pool[,i]   
    # necessary statistics
    tempRes = adaptive(X1, Y1, X2, Y2, beta, lam1[i], sig, alpha) 
    r[i] = tempRes$r
    coverage[i] = tempRes$coverage
    logVol[i] = tempRes$logVol
  }
  return(data.frame(r, coverage, logVol))
}

# warp all cases given a design matrix
adaptive_single = function(n, p, svalue, bvalue, sig, lam1Types, alpha, ite, rawDir)
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
        arg[['root']] = rawDir
        arg[['beta']] = paste(rawDir, '/beta/', 'beta', '_s', svalue[j], '_b', bvalue[k], '.txt', sep = "")
        arg[['Y1']] = paste(rawDir, '/Y1/', 'Y1', '_s', svalue[j], '_b', bvalue[k], '.txt', sep = "")
        arg[['Y2']] = paste(rawDir, '/Y2/', 'Y2', '_s', svalue[j], '_b', bvalue[k], '.txt', sep = "")  
        arg[['lam1Type']] = lam1Types[i]
        if (arg[['lam1Type']] == '1se') {
            arg[['1se']] = paste(rawDir, '/1se/', '1se', '_s', svalue[j], '_b', bvalue[k], '.txt', 
                               sep = "")
        }
        else if (arg[['lam1Type']] == 'cv') {
            arg[['cv']] = paste(rawDir, '/cv/', 'cv', '_s', svalue[j], '_b', bvalue[k], '.txt', 
                               sep = "")
        }
        # calculate necessary statistics
        tempRes = adaptive_test(n, p, sig, alpha, ite, arg)
        tempRes['lam1Type'] = lam1Types[i]
        tempRes['s'] = svalue[j]
        tempRes['b'] = bvalue[k]
        tempRes['method'] = 'adaptive'
        res_list[[count]] = tempRes
        count = count + 1
      }
    }
  }
  return(ldply(res_list))
}

# mulitple signals
adaptive_multiple = function(n, p, s_s, b_s, s_w, b_w, sig, lam1Types, alpha, ite, rawDir)
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
        arg[['root']] = rawDir
        arg[['beta']] = paste(rawDir, '/beta', '/beta', '_ss', s_s, '_bs', b_s[j], '_sw', s_w, 
                              '_bw', b_w[k], '.txt', sep="")
        arg[['Y1']] = paste(rawDir, '/Y1', '/Y1', '_ss', s_s, '_bs', b_s[j], '_sw', s_w, 
                            '_bw', b_w[k], '.txt', sep="")
        arg[['Y2']] = paste(rawDir, '/Y2', '/Y2', '_ss', s_s, '_bs', b_s[j], '_sw', s_w, 
                            '_bw', b_w[k], '.txt', sep="")
        arg[['lam1Type']] = lam1Types[i]
        if (arg[['lam1Type']] == '1se') {
            arg[['1se']] = paste(rawDir, '/1se', '/1se', '_ss', s_s, '_bs', b_s[j], '_sw', s_w, 
                                 '_bw', b_w[k], '.txt', sep="")
        }
        else if (arg[['lam1Type']] == 'cv') {
            arg[['cv']] = paste(rawDir, '/cv', '/cv', '_ss', s_s, '_bs', b_s[j], '_sw', s_w, 
                                '_bw', b_w[k], '.txt', sep="")
        }
        # calculate necessary statistics
        tempRes = adaptive_test(n, p, sig, alpha, ite, arg)
        tempRes['lam1Type'] = lam1Types[i]
        tempRes['bs'] = b_s[j]
        tempRes['bw'] = b_w[k]
        tempRes['method'] = 'adaptive'
        res_list[[count]] = tempRes
        count = count + 1
      }
    }
  }
  return(ldply(res_list))
}

# library(plyr)

