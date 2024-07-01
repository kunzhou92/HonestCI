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
  k = sum(lassoModel$beta != 0)
  r = sig * sqrt(cAlpha * s * log(p) / n)
  normalizedSquaredLoss = L2NormSquare(X %*% (beta - lassoModel$beta)) / n
  rv = normalizedSquaredLoss / (r^2)
  coverage = rv <= 1
  logVol = n * log(r)
  result = list(r = r, coverage=coverage, logVol=logVol, cAlpha=cAlpha, 
                normalizedSquaredLoss=normalizedSquaredLoss, k=k)
  return(result)
}

oracle_test = function(n, p, sig, alpha, ite, arg)
{
  # add ite_quantile as argument eventually
  ite_quantile = 100
  r = rep(0, ite)
  coverage = rep(0, ite)
  logVol = rep(0, ite)
  cAlpha = rep(0, ite)
  normalizedSquaredLoss = rep(0, ite)
  k = rep(0, ite)

  beta_pool = t(as.matrix(read.table(arg[['beta']])))
  Y1_pool = t(as.matrix(read.table(arg[['Y1']])))
  Y2_pool = t(as.matrix(read.table(arg[['Y2']]))) 
  
  if (arg[['lam1Type']] == 'val') {
      lam1 = rep(2 * sqrt(2) * sqrt(log(p) / n), ite)
  }
  for (i in 1:ite)
  {
    fileX1 = paste(arg[['root']], '/X1/', i, '.txt', sep = "") 
    fileX2 = paste(arg[['root']], '/X2/', i, '.txt', sep = "") 

    X1 = as.matrix(read.table(fileX1))
    X2 = as.matrix(read.table(fileX2))    
    beta = beta_pool[, i, drop = F]
    Y1 = Y1_pool[, i, drop = F]
    Y2 = Y2_pool[, i, drop = F]    
    X = rbind(X1, X2)
    Y = rbind(Y1, Y2)
    # calculate necessary statistics
    tempRes = LassoQuantile(X, Y, beta, arg[['beta_fun']], lam1[i],  alpha, ite_quantile, sig)
    r[i] = tempRes$r
    coverage[i] = tempRes$coverage
    logVol[i] = tempRes$logVol
    cAlpha[i] = tempRes$cAlpha
    normalizedSquaredLoss[i] = tempRes$normalizedSquaredLoss
    k[i] = tempRes$k
  }
  return(data.frame(r, coverage, logVol, cAlpha, normalizedSquaredLoss, k))
}

oracle_single = function(n, p, svalue, bvalue, sig, lam1Types, alpha, ite, rawDir)
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
 
oracle_multiple = function(n, p, s_s, b_s, s_w, b_w, sig, lam1Types, alpha, ite, rawDir)
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


# a sequence of beta is stored in matrix by column
# arg: nfolds
# beta_pool[, i, drop = F]
# betas, Ys, squaredR: same order

adjust_lasso = function(X, betas, Ys, val_df, arg)
{
  n = dim(X)[1]
  total = dim(betas)[2]
  normalizedSquaredLoss_cv = rep(0, total)
  normalizedSquaredLoss_1se = rep(0, total)
  k_cv = rep(0, total)
  k_1se = rep(0, total)
  squaredR = (val_df$r) ^ 2

  for(i in 1:total)
  {
    beta = betas[, i, drop = F]
    Y = Ys[, i, drop = F]
    # compute beta_hat
    cv_glmnet = cv.glmnet(X, Y, nfolds = arg[['nfolds']], intercept = F)
    lambda_cv = cv_glmnet$lambda.min
    lambda_1se = cv_glmnet$lambda.1se
    model_cv = glmnet(X, Y, family = 'gaussian', intercept = F, lambda = lambda_cv)
    model_1se = glmnet(X, Y, family = 'gaussian', intercept = F, lambda = lambda_1se)
    # compute k
    k_cv[i] = sum(model_cv$beta != 0)
    k_1se[i] = sum(model_1se$beta != 0)
    # compute loss
    normalizedSquaredLoss_cv[i] = L2NormSquare(X %*% (beta - model_cv$beta)) / n
    normalizedSquaredLoss_1se[i] = L2NormSquare(X %*% (beta - model_1se$beta)) / n
  }
  count = 0
  ratio_cv = normalizedSquaredLoss_cv / squaredR
  coef_cv = adjust_coef(ratio_cv, 1, count, arg)
  ratio_1se = normalizedSquaredLoss_1se / squaredR
  coef_1se = adjust_coef(ratio_1se, 1, count, arg)
  # summarize

  r_adjusted_cv = sqrt(coef_cv * squaredR)
  coverage_cv = normalizedSquaredLoss_cv / (coef_cv * squaredR) <= 1
  logVol_cv = n * log(r_adjusted_cv)
  cAlpha_cv = coef_cv * val_df$cAlpha
  df_cv = data.frame(r = r_adjusted_cv,
                     coverage = coverage_cv,
                     logVol = logVol_cv,
                     cAlpha = cAlpha_cv,
                     normalizedSquaredLoss = normalizedSquaredLoss_cv,
                     k = k_cv)

  r_adjusted_1se = sqrt(coef_1se * squaredR)
  coverage_1se = normalizedSquaredLoss_1se / (coef_1se * squaredR) <= 1
  logVol_1se = n * log(r_adjusted_1se)
  cAlpha_1se = coef_1se * val_df$cAlpha
  df_1se = data.frame(r = r_adjusted_1se,
                      coverage = coverage_1se,
                      logVol = logVol_1se,
                      cAlpha = cAlpha_1se,
                      normalizedSquaredLoss = normalizedSquaredLoss_1se,
                      k = k_1se)
  return(list('df_cv' = df_cv, 'df_1se' = df_1se))
}

# arg: facotr, lower, upper, stop
# count for stop
# trick: let r be divided by max_r
if (F)
{
  adjust_coef = function(ratio, coef, count, arg)
  {
    if (count > arg[['stop']])
      return(coef)
    coverage = mean(ratio <= 1)
    if (coverage < arg[['lower']])
    {
      return(adjust_coef(ratio / arg[['factor']], coef * arg[['factor']], count + 1, arg))
    }
    else if (coverage > arg[['upper']])
    {
      return(adjust_coef(ratio * arg[['factor']], coef / arg[['factor']], count + 1, arg))
    }
    else
    {
      return(coef)
    }
  }
}

adjust_coef = function(ratio, coef, count, arg)
{
  coverage = mean(ratio <= 1)
  while((coverage < arg[['lower']] || coverage > arg[['upper']]) && (count <= arg[['stop']]))
  {
    if (coverage < arg[['lower']])
    {
      ratio = ratio / arg[['factor']]
      coef = coef * arg[['factor']]
    }
    else if (coverage > arg[['upper']])
    {
      ratio = ratio * arg[['factor']]
      coef = coef / arg[['factor']]
    }
    count = count + 1
    coverage = mean(ratio <= 1)
  }
  return(coef)
}




# n should be doubled
adjust_lasso_single = function(n, p, svalue, bvalue, sig, ite, rawDir,
                               summary_val, factor, lower, upper, stop, nfolds, design)
{
  arg = list('factor' = factor, 'lower' = lower, 'upper' = upper,
             'stop' = stop, 'nfolds' = nfolds)
  res_list_cv = list()
  res_list_1se = list()
  
  for (i in 1:ite)
  {
    # get X
    fileX1 = paste(rawDir, '/X1/', i, '.txt', sep = "") 
    fileX2 = paste(rawDir, '/X2/', i, '.txt', sep = "") 
    X1 = as.matrix(read.table(fileX1))
    X2 = as.matrix(read.table(fileX2)) 
    X = rbind(X1, X2)
    # get Ys, betas
    total = length(svalue) * length(bvalue)
    Ys = matrix(0, nrow = n, ncol = total)
    betas = matrix(0, nrow = p, ncol = total)
    count = 1
    for (j in 1:length((svalue)))
    {
      for (k in 1:length(bvalue))
      {
        fileBeta = paste(rawDir, '/beta/', 'beta', '_s', svalue[j], '_b', bvalue[k], '.txt', sep = "")
        fileY1 = paste(rawDir, '/Y1/', 'Y1', '_s', svalue[j], '_b', bvalue[k], '.txt', sep = "")
        fileY2 = paste(rawDir, '/Y2/', 'Y2', '_s', svalue[j], '_b', bvalue[k], '.txt', sep = "") 
        betas[, count] = as.numeric(read.table(fileBeta, skip = i - 1, nrows = 1))
        Y1 = as.numeric(read.table(fileY1, skip = i - 1, nrows = 1))
        Y2 = as.numeric(read.table(fileY2, skip = i - 1, nrows = 1))
        Ys[, count] = c(Y1, Y2)
        count = count + 1
      }
    }
    val_df = summary_val[summary_val$design == design,][((1:total) - 1) * ite + i,]
    res = adjust_lasso(X, betas, Ys, val_df, arg)
    res_list_cv[[i]] = res[['df_cv']]
    res_list_1se[[i]] = res[['df_1se']]
  }
  df_s_b = expand.grid(bvalue, svalue)
  colnames(df_s_b) = c('b', 's')

  df_cv = ldply(res_list_cv)
  df_cv['lam1Type'] = 'cv'
  df_cv['method'] = 'oracle'
  df_cv = cbind(df_cv, df_s_b)

  df_1se = ldply(res_list_1se)
  df_1se['lam1Type'] = '1se'
  df_1se['method'] = 'oracle'
  df_1se = cbind(df_1se, df_s_b)
  return(list(df_cv = df_cv, df_1se = df_1se))
}

########## adjust cAlpha so that coverage is around 0.95 across b > 0.3 ##########################


adjust_lasso_version2 = function(X, betas, Ys, val_df, arg)
{
  n = dim(X)[1]
  total = dim(betas)[2]
  normalizedSquaredLoss_cv = rep(0, total)
  normalizedSquaredLoss_1se = rep(0, total)
  k_cv = rep(0, total)
  k_1se = rep(0, total)
  squaredR = (val_df$r) ^ 2

  for(i in 1:total)
  {
    beta = betas[, i, drop = F]
    Y = Ys[, i, drop = F]
    # compute beta_hat
    cv_glmnet = cv.glmnet(X, Y, nfolds = arg[['nfolds']], intercept = F)
    lambda_cv = cv_glmnet$lambda.min
    lambda_1se = cv_glmnet$lambda.1se
    model_cv = glmnet(X, Y, family = 'gaussian', intercept = F, lambda = lambda_cv)
    model_1se = glmnet(X, Y, family = 'gaussian', intercept = F, lambda = lambda_1se)
    # compute k
    k_cv[i] = sum(model_cv$beta != 0)
    k_1se[i] = sum(model_1se$beta != 0)
    # compute loss
    normalizedSquaredLoss_cv[i] = L2NormSquare(X %*% (beta - model_cv$beta)) / n
    normalizedSquaredLoss_1se[i] = L2NormSquare(X %*% (beta - model_1se$beta)) / n
  }
  count = 0
  ratio_cv = normalizedSquaredLoss_cv / squaredR
  coef_cv = adjust_coef(ratio_cv[arg[['index']]], 1, count, arg)
  ratio_1se = normalizedSquaredLoss_1se / squaredR
  coef_1se = adjust_coef(ratio_1se[arg[['index']]], 1, count, arg)

  # summarize

  r_adjusted_cv = sqrt(coef_cv * squaredR)
  coverage_cv = normalizedSquaredLoss_cv / (coef_cv * squaredR) <= 1
  logVol_cv = n * log(r_adjusted_cv)
  cAlpha_cv = coef_cv * val_df$cAlpha
  df_cv = data.frame(r = r_adjusted_cv,
                     coverage = coverage_cv,
                     logVol = logVol_cv,
                     cAlpha = cAlpha_cv,
                     normalizedSquaredLoss = normalizedSquaredLoss_cv,
                     k = k_cv)
  
  r_adjusted_1se = sqrt(coef_1se * squaredR)
  coverage_1se = normalizedSquaredLoss_1se / (coef_1se * squaredR) <= 1
  logVol_1se = n * log(r_adjusted_1se)
  cAlpha_1se = coef_1se * val_df$cAlpha
  df_1se = data.frame(r = r_adjusted_1se,
                      coverage = coverage_1se,
                      logVol = logVol_1se,
                      cAlpha = cAlpha_1se,
                      normalizedSquaredLoss = normalizedSquaredLoss_1se,
                      k = k_1se)
  return(list('df_cv' = df_cv, 'df_1se' = df_1se))
}


# adjust so that coverage > 0.95 across all b > 0.3
adjust_lasso_single_version2 = function(n, p, svalue, bvalue, sig, ite, rawDir,
                               summary_val, factor, lower, upper, stop, nfolds, design)
{
  arg = list('factor' = factor, 'lower' = lower, 'upper' = upper,
             'stop' = stop, 'nfolds' = nfolds)
  res_list_cv = list()
  res_list_1se = list()
  
  for (i in 1:ite)
  {
    # get X
    fileX1 = paste(rawDir, '/X1/', i, '.txt', sep = "") 
    fileX2 = paste(rawDir, '/X2/', i, '.txt', sep = "") 
    X1 = as.matrix(read.table(fileX1))
    X2 = as.matrix(read.table(fileX2)) 
    X = rbind(X1, X2)
    # get Ys, betas
    total = length(svalue) * length(bvalue)
    Ys = matrix(0, nrow = n, ncol = total)
    betas = matrix(0, nrow = p, ncol = total)
    count = 1
    index = c()
    for (j in 1:length((svalue)))
    {
      for (k in 1:length(bvalue))
      {
        if (bvalue[k] > 0.3)
        {
          index = c(index, count)
        }
        fileBeta = paste(rawDir, '/beta/', 'beta', '_s', svalue[j], '_b', bvalue[k], '.txt', sep = "")
        fileY1 = paste(rawDir, '/Y1/', 'Y1', '_s', svalue[j], '_b', bvalue[k], '.txt', sep = "")
        fileY2 = paste(rawDir, '/Y2/', 'Y2', '_s', svalue[j], '_b', bvalue[k], '.txt', sep = "") 
        betas[, count] = as.numeric(read.table(fileBeta, skip = i - 1, nrows = 1))
        Y1 = as.numeric(read.table(fileY1, skip = i - 1, nrows = 1))
        Y2 = as.numeric(read.table(fileY2, skip = i - 1, nrows = 1))
        Ys[, count] = c(Y1, Y2)
        count = count + 1
      }
    }
    arg[['index']] = index
    val_df = summary_val[summary_val$design == design,][((1:total) - 1) * ite + i,]
    res = adjust_lasso_version2(X, betas, Ys, val_df, arg)
    res_list_cv[[i]] = res[['df_cv']]
    res_list_1se[[i]] = res[['df_1se']]
  }
  df_s_b = expand.grid(bvalue, svalue)
  colnames(df_s_b) = c('b', 's')

  df_cv = ldply(res_list_cv)
  df_cv['lam1Type'] = 'cv'
  df_cv['method'] = 'oracle'
  df_cv = cbind(df_cv, df_s_b)

  df_1se = ldply(res_list_1se)
  df_1se['lam1Type'] = '1se'
  df_1se['method'] = 'oracle'
  df_1se = cbind(df_1se, df_s_b)
  return(list(df_cv = df_cv, df_1se = df_1se))
}


adjust_lasso_multiple_version2 = function(n, p, s_s, b_s, s_w, b_w, sig, ite, rawDir,
                               summary_val, factor, lower, upper, stop, nfolds, design)
{
  arg = list('factor' = factor, 'lower' = lower, 'upper' = upper,
             'stop' = stop, 'nfolds' = nfolds)
  res_list_cv = list()
  res_list_1se = list()
  
  for (i in 1:ite)
  {
    # get X
    fileX1 = paste(rawDir, '/X1/', i, '.txt', sep = "") 
    fileX2 = paste(rawDir, '/X2/', i, '.txt', sep = "") 
    X1 = as.matrix(read.table(fileX1))
    X2 = as.matrix(read.table(fileX2)) 
    X = rbind(X1, X2)
    # get Ys, betas
    total = length(b_s) * length(b_w)
    Ys = matrix(0, nrow = n, ncol = total)
    betas = matrix(0, nrow = p, ncol = total)
    count = 1
    index = c()
    for (j in 1:length((b_s)))
    {
      for (k in 1:length(b_w))
      {
        if (b_s[j] > 0.3)
        {
          index = c(index, count)
        }
        fileBeta = paste(rawDir, '/beta', '/beta', '_ss', s_s, '_bs', b_s[j], '_sw', s_w, 
                         '_bw', b_w[k], '.txt', sep="")
        fileY1 = paste(rawDir, '/Y1', '/Y1', '_ss', s_s, '_bs', b_s[j], '_sw', s_w, 
                       '_bw', b_w[k], '.txt', sep="")
        fileY2 = paste(rawDir, '/Y2', '/Y2', '_ss', s_s, '_bs', b_s[j], '_sw', s_w, 
                       '_bw', b_w[k], '.txt', sep="")
        betas[, count] = as.numeric(read.table(fileBeta, skip = i - 1, nrows = 1))
        Y1 = as.numeric(read.table(fileY1, skip = i - 1, nrows = 1))
        Y2 = as.numeric(read.table(fileY2, skip = i - 1, nrows = 1))
        Ys[, count] = c(Y1, Y2)
        count = count + 1
      }
    }
    arg[['index']] = index
    val_df = summary_val[summary_val$design == design,][((1:total) - 1) * ite + i,]
    res = adjust_lasso_version2(X, betas, Ys, val_df, arg)
    res_list_cv[[i]] = res[['df_cv']]
    res_list_1se[[i]] = res[['df_1se']]
  }
  df_s_b = expand.grid(b_w, b_s)
  colnames(df_s_b) = c('bw', 'bs')

  df_cv = ldply(res_list_cv)
  df_cv['lam1Type'] = 'cv'
  df_cv['method'] = 'oracle'
  df_cv = cbind(df_cv, df_s_b)

  df_1se = ldply(res_list_1se)
  df_1se['lam1Type'] = '1se'
  df_1se['method'] = 'oracle'
  df_1se = cbind(df_1se, df_s_b)
  return(list(df_cv = df_cv, df_1se = df_1se))
}
