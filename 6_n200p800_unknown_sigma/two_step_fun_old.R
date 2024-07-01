# the main change is function "threshold" now decide the number of candidate sets we want and then
# evenly space from empty to supp(\beta).
# remove tau_chosen
# qr$qr is the combination of Q and R. use qr.Q and qr.R to recover

library(glmnet)
library(plyr)
library(LDRTools)
source("data_fun.R")

# Shrink Y towards zero by Stein estimate
SteinToZero = function(Y, df, sig)
{
  B = sig^2 * df / L2NormSquare(Y)
  mu = (1 - B) * Y
  sure = sig^2 * max(1 - B, 0)
  return(list(mu = mu, sure = sure))
}

# calculate c(alpha) for Stein. (for squared radius)
quantile_stein = function(mu, df, alpha, sig, iteNum, projection=NULL)
{
  n = dim(mu)[1]
  if (is.null(projection))
    projection = diag(rep(1, n))
  z = rep(0, iteNum)
  for (i in 1:iteNum)
  {
    Y_sample = mu + tcrossprod(projection, matrix(rnorm(n, sd = sig), nrow = 1))
    res = SteinToZero(Y_sample, df, sig)
    z[i] = abs(res$sure - L2NormSquare(res$mu - mu) / df) * sqrt(df) / (sig^2)
  }
  c = quantile(z, probs = 1 - alpha)
  return (c)
}

# truncate beta
trunate_beta = function(beta, s)
{
  idx = order(abs(beta), decreasing = TRUE)[1:s] # 
  beta[-idx,] = 0
  return(beta)
}

# a practical method of calculating c(alpha) for Lasso. (for squared radius)
quantile_lasso = function(X1, Y1, X2, s, sig, lam2, alpha, iteNum)
{
  n2 = dim(X2)[1]
  p = dim(X1)[2]
  # compute a reasonable bound of b
  b_vec = rep(0, p)
  for (i in 1:p)
  {
    b_vec[i] = crossprod(X1[, i], Y1)[1,1] / L2NormSquare(X1[, i])
  }
  b = max(abs(b_vec))
  # simluation q
  z = rep(0, iteNum)
  for (i in 1:iteNum)
  {
    beta_sample = generateBeta(type = 3, p, s, b)
    Y_sample = tcrossprod(X2, t(beta_sample)) + matrix(rnorm(n2, sd = sig), ncol = 1)
    lassoModel2 = glmnet(X2, Y_sample, family='gaussian', intercept = F, lambda = lam2)
    beta_tilde = trunate_beta(as.matrix(lassoModel2$beta), s)
    z[i] = L2NormSquare(tcrossprod(X2, t(beta_tilde - beta_sample))) / s / log(p) / (sig^2)
  }
  c = quantile(z, probs = 1 - alpha)
  return(c)
}

weak_stein = function(X2, Y2, A, sig, alpha_s, c_alpha_w, upper)
{
  n2 = dim(X2)[1]
  p = dim(X2)[2]
  # get k
  qr_decom = qr(X2[,A])
  k = qr_decom$rank
  if (k == 0)
  {
    P_A = matrix(0, ncol = n2, nrow = n2)
    P_perp = diag(rep(1, n2))
  } else
  {
    Q = qr.Q(qr_decom)
    R = qr.R(qr_decom)[,1:k]
    P_A = B2P(Q %*% R)
    P_perp = diag(rep(1, n2)) - P_A
  }
  
  # strong signals
  if (k == 0)
  {
    mu_s = 0
    r_s = 0
  } else
  {
    mu_s = P_A %*% Y2
    r_s = sqrt(qchisq(1-alpha_s, df = k) / n2) * sig
  }
  # weak signals
  if (k == n2)
  {
    r_w = 0
    mu_w = 0
  } else {
    stein = SteinToZero(tcrossprod(P_perp, t(Y2)), n2 - k, sig)                                  
    r_w = sqrt((n2 - k) / n2 * (c_alpha_w * sig^2 / sqrt(n2 - k) + stein$sure))
    mu_w = stein$mu
  }
  
  # update r_s_ratio r_w_ratio if A is not empty nor rank(X_A) == n2
  if (k == 0)
  {
    r_s_volume = 0
    r_w_volume = r_w
    logVol_volume = n2 * log(r_w_volume)
    r_radius = r_w
    logVol_radius = n2 * log(r_radius)
  } else if (k == n2)
  {
    r_s_volume = r_s
    r_w_volume = 0
    logVol_volume = n2 * log(r_s_volume)
    r_radius = r_s
    logVol_radius = n2 * log(r_radius)
  } else
  {
    # by volume
    lower = upper / (upper - 1) 
    c_s = max(lower, min(n2 / k, upper))
    c_w = c_s / (c_s - 1)
    r_s_volume = r_s * sqrt(c_s)
    r_w_volume = r_w * sqrt(c_w)
    logVol_volume = k * log(r_s_volume) + (n2 - k) * log(r_w_volume)
    # by radius
    r_radius = sqrt(r_s^2 + r_w^2)
    logVol_radius = n2 * log(r_radius)
  }
  return(list(r_s_volume = r_s_volume,
              r_w_volume = r_w_volume,
              logVol_volume = logVol_volume,
              r_radius = r_radius,
              logVol_radius = logVol_radius,
              mu_s = mu_s,
              mu_w = mu_w))
}

# ensure |A| < s before implementing
weak_lasso = function(X2, Y2, A, sig, alpha_s, c_alpha_w, lam2Coef, s, upper)
{
  n2 = dim(X2)[1]
  p = dim(X2)[2]
  # get k
  qr_decom = qr(X2[,A])
  k = qr_decom$rank
  if (k == 0)
  {
    P_A = matrix(0, ncol = n2, nrow = n2)
    P_perp = diag(rep(1, n2))
  } else
  {
    Q = qr.Q(qr_decom)
    R = qr.R(qr_decom)[,1:k]
    P_A = B2P(Q %*% R)
    P_perp = diag(rep(1, n2)) - P_A
  }
  
  # strong signals
  if (k == 0)
  {
    mu_s = 0
    r_s = 0
  } else
  {
    mu_s = P_A %*% Y2
    r_s = sqrt(qchisq(1 - alpha_s, df = k) / n2) * sig
  }
  # weak signals
  # sparsity = max(1, s - k) # min_sparsity = 1 
  if (k == n2)
  {
    r_w = 0
    mu_w = 0
  } else
  {
    sparsity = s - k # min_sparsity = 1 
    lam2 = lam2Coef * 2 * sqrt(2) * sig * sqrt(log(p - k) / (n2 - k))
    lasso_weak = glmnet(P_perp %*% X2, P_perp %*% Y2, intercept = F, lambda = lam2)
    mu_w = P_perp %*% X2 %*% trunate_beta(as.matrix(lasso_weak$beta), sparsity)
    r_w = sig * sqrt(c_alpha_w * sparsity * log(p - k) / n2)  # then r_w only depends on the size of A 
  }
  
  # update r_s_ratio r_w_ratio if A is not empty
  if (k == 0)
  {
    r_s_volume = 0
    r_w_volume = r_w
    logVol_volume = n2 * log(r_w_volume)
    r_radius = r_w
    logVol_radius = n2 * log(r_radius)
  } else if (k == n2)
  {
    r_s_volume = r_s
    r_w_volume = 0
    logVol_volume = n2 * log(r_s_volume)
    r_radius = r_s
    logVol_radius = n2 * log(r_radius)
  }
  else
  {
    # by volume
    lower = upper / (upper - 1) 
    c_s = max(lower, min(n2 / k, upper))
    c_w = c_s / (c_s - 1)
    r_s_volume = r_s * sqrt(c_s)
    r_w_volume = r_w * sqrt(c_w)
    logVol_volume = k * log(r_s_volume) + (n2 - k) * log(r_w_volume)
    # by radius
    r_radius = sqrt(r_s^2 + r_w^2)
    logVol_radius = n2 * log(r_radius)
  }
  return(list(r_s_volume = r_s_volume,
              r_w_volume = r_w_volume,
              logVol_volume = logVol_volume,
              r_radius = r_radius,
              logVol_radius = logVol_radius,
              mu_s = mu_s,
              mu_w = mu_w))
}

threshold = function(beta, numCandidate, X2, max_size = NULL) 
{
  beta = as.numeric(beta)
  nonzero_beta_indices = (1:length(beta))[beta != 0]
  beta = beta[beta != 0]
  if (is.null(max_size))
  {
    size_of_candidates = unique(round(seq(from = 0, to = length(beta), length.out = numCandidate)))
  } else
  {
    size_of_candidates = unique(round(seq(from = 0, to = min(length(beta), max_size), length.out = numCandidate)))
  }
  
  order_nonzero_beta = order(abs(beta), decreasing = TRUE)
  candidates = list()
  
  for (i in 1:length(size_of_candidates))
  {
    A = nonzero_beta_indices[order_nonzero_beta[0:size_of_candidates[i]]]
    qr_decom = qr(X2[,A])
    k = qr_decom$rank
    if (i != 1 && k == pre_k)
    {
      next
    }
    pre_k = k
    candidates[[i]] = A
  }
  return(list(candidates = candidates,
              count = length(size_of_candidates)))
}

two_step = function(X1, X2, Y1, Y2, lam1, lam2Coef, s, sig, numCandidate, alpha, upper, ite)
{
  n1 = dim(X1)[1] 
  n2 = dim(X2)[1]
  p = dim(X1)[2]
  # stein c(alpha)
  mu = matrix(rep(0, n2), ncol = 1)
  c_alpha_stein = quantile_stein(mu, n2, alpha / 2, sig, ite)
  # lasso c(alpha)
  lam2 = lam2Coef * 2 * sqrt(2) * sig * sqrt(log(p) / n2)
  c_alpha_lasso = quantile_lasso(X1, Y1, X2, s, sig, lam2, alpha / 2, ite)
  
  # candidates
  lasso_split = glmnet(X1, Y1, family='gaussian', intercept = F, lambda = lam1)
  res_candidate = threshold(as.matrix(lasso_split$beta), numCandidate, X2)
  # trunate beta to s-1 nonzero coefs, so |A| <= s-1
  res_candidate_lasso = threshold(as.matrix(lasso_split$beta), numCandidate, X2, s-1)
  

  # statistics
  res_stein_pool = list()
  logVol_volume_stein = rep(0, res_candidate[['count']])
  logVol_radius_stein = rep(0, res_candidate[['count']])

  res_lasso_pool = list()
  logVol_volume_lasso = rep(0, res_candidate_lasso[['count']])
  logVol_radius_lasso = rep(0, res_candidate_lasso[['count']])

  # record
  for (i in 1:res_candidate[['count']])
  {
    res_stein = weak_stein(X2, Y2, res_candidate[['candidates']][[i]], sig, alpha / 2, c_alpha_stein, upper)
    res_stein_pool[[i]] = res_stein
    logVol_volume_stein[i] = res_stein[['logVol_volume']]
    logVol_radius_stein[i] = res_stein[['logVol_radius']]
  }

  for (i in 1:res_candidate_lasso[['count']])
  {
    res_lasso = weak_lasso(X2, Y2, res_candidate_lasso[['candidates']][[i]], sig, alpha / 2, c_alpha_lasso, lam2Coef, s, upper)
    res_lasso_pool[[i]] = res_lasso
    logVol_volume_lasso[i] = res_lasso[['logVol_volume']]
    logVol_radius_lasso[i] = res_lasso[['logVol_radius']]
  }

  i_volume_stein = which.min(logVol_volume_stein)
  i_radius_stein = which.min(logVol_radius_stein)
  i_volume_lasso = which.min(logVol_volume_lasso)
  i_radius_lasso = which.min(logVol_radius_lasso)
  
  volume_stein = list(r_s = res_stein_pool[[i_volume_stein]][['r_s_volume']], 
                      r_w = res_stein_pool[[i_volume_stein]][['r_w_volume']],
                      mu_s = res_stein_pool[[i_volume_stein]][['mu_s']],
                      mu_w = res_stein_pool[[i_volume_stein]][['mu_w']],
                      logVol = res_stein_pool[[i_volume_stein]][['logVol_volume']],
                      A = res_candidate[['candidates']][[i_volume_stein]],
                      k = length(res_candidate[['candidates']][[i_volume_stein]]))

  radius_stein = list(r_s = res_stein_pool[[i_radius_stein]][['r_radius']],
                      r_w = res_stein_pool[[i_radius_stein]][['r_radius']],
                      mu_s = res_stein_pool[[i_radius_stein]][['mu_s']],
                      mu_w = res_stein_pool[[i_radius_stein]][['mu_w']],
                      logVol = res_stein_pool[[i_radius_stein]][['logVol_radius']],
                      A = res_candidate[['candidates']][[i_radius_stein]],
                      k = length(res_candidate[['candidates']][[i_radius_stein]]))

  volume_lasso = list(r_s = res_lasso_pool[[i_volume_lasso]][['r_s_volume']],
                      r_w = res_lasso_pool[[i_volume_lasso]][['r_w_volume']],
                      mu_s = res_lasso_pool[[i_volume_lasso]][['mu_s']],
                      mu_w = res_lasso_pool[[i_volume_lasso]][['mu_w']],
                      logVol = res_lasso_pool[[i_volume_lasso]][['logVol_volume']],
                      A = res_candidate_lasso[['candidates']][[i_volume_lasso]],
                      k = length(res_candidate_lasso[['candidates']][[i_volume_lasso]]))
                      
  radius_lasso = list(r_s = res_lasso_pool[[i_radius_lasso]][['r_radius']],
                      r_w = res_lasso_pool[[i_radius_lasso]][['r_radius']],
                      mu_s = res_lasso_pool[[i_radius_lasso]][['mu_s']],
                      mu_w = res_lasso_pool[[i_radius_lasso]][['mu_w']],
                      logVol = res_lasso_pool[[i_radius_lasso]][['logVol_radius']],
                      A = res_candidate_lasso[['candidates']][[i_radius_lasso]],
                      k = length(res_candidate_lasso[['candidates']][[i_radius_lasso]]))

  return(list(volume_stein = volume_stein, radius_stein = radius_stein, 
              volume_lasso = volume_lasso, radius_lasso = radius_lasso))
}

convert = function(m_list, X2, beta) {
  n2 = dim(X2)[1]
  p = dim(X2)[2]
  if (m_list[['k']] == 0) 
  {
    part1 = m_list[["mu_w"]] - X2 %*% beta
    part2 = diag(rep(1, n2)) / m_list[["r_w"]]^2
  } else if (m_list[['k']] == n2)
  {
    part1 = m_list[["mu_s"]] - X2 %*% beta
    part2 = diag(rep(1, n2)) / m_list[["r_s"]]^2
  }
  else
  {
    P_A = B2P(X2[,m_list[['A']]])
    P_perp = diag(rep(1, n2)) - P_A  
    part1 = m_list[["mu_s"]] + m_list[["mu_w"]] - X2 %*% beta
    part2 = P_A / m_list[["r_s"]]^2  + P_perp / m_list[["r_w"]]^2
  }
  temp = (t(part1) %*% part2 %*% part1)/ n2 
  return(list(coverage = as.numeric(temp) <= 1,
              r_s = m_list[['r_s']],
              r_w = m_list[['r_w']],
              k = m_list[['k']],
              logVol = m_list[['logVol']]))
}

# ite for computing quantile is set here
two_step_test = function(n, p, sig, alpha, ite, arg)
{
  ite_quantile = 100
  volume_stein_pool = list()
  radius_stein_pool = list()
  volume_lasso_pool = list()
  radius_lasso_pool = list()

  beta_pool = t(as.matrix(read.table(arg[['beta']])))
  Y1_pool = t(as.matrix(read.table(arg[['Y1']])))
  Y2_pool = t(as.matrix(read.table(arg[['Y2']]))) 

  if (arg[['lam1Type']] == '1se') {
      lam1 = read.table(arg[['1se']])[,1]
  } else if (arg[['lam1Type']] == 'cv') {
      lam1 = read.table(arg[['cv']])[,1]
  } else if (arg[['lam1Type']] == 'val') {
      lam1 = rep(2 * sqrt(2) * sqrt(log(p) / n), ite)
  }
  lam2Coef = arg[['lam2Coef']]
  s = arg[['s']]
  numCandidate = arg[['numCandidate']]
  upper = arg[['upper']]
  #foreach (i = 1:ite, .export = c("SteinToZero", "quantile_stein", "trunate_beta", "quantile_lasso", "weak_stein", "weak_lasso", "threshold", "convert", "two_step", "L2NormSquare", "generateBeta"), .packages = c("glmnet")) %dopar%
  for( i in 1:ite)
  {
    fileX1 = paste(arg[['root']], '/X1/', i, '.txt', sep = "") 
    fileX2 = paste(arg[['root']], '/X2/', i, '.txt', sep = "") 

    X1 = as.matrix(read.table(fileX1))
    X2 = as.matrix(read.table(fileX2))
    beta = beta_pool[, i, drop = F]
    Y1 = Y1_pool[, i, drop = F]
    Y2 = Y2_pool[, i, drop = F]   
    # necessary statistics
    res = two_step(X1, X2, Y1, Y2, lam1[i], lam2Coef, s, sig, numCandidate, alpha, upper, ite_quantile)
    volume_stein = convert(res[['volume_stein']], X2, beta)
    volume_stein_pool[[i]] = volume_stein
    radius_stein = convert(res[['radius_stein']], X2, beta)
    radius_stein_pool[[i]] = radius_stein
    volume_lasso = convert(res[['volume_lasso']], X2, beta)
    volume_lasso_pool[[i]] = volume_lasso
    radius_lasso = convert(res[['radius_lasso']], X2, beta)
    radius_lasso_pool[[i]] = radius_lasso
  }

  volume_stein_df = ldply(volume_stein_pool, data.frame)
  radius_stein_df = ldply(radius_stein_pool, data.frame)
  volume_lasso_df = ldply(volume_lasso_pool, data.frame)
  radius_lasso_df = ldply(radius_lasso_pool, data.frame)
  return(list(volume_stein_df = volume_stein_df, 
              radius_stein_df = radius_stein_df,
              volume_lasso_df = volume_lasso_df,
              radius_lasso_df = radius_lasso_df))
}

# set the oracle sparsity as the true sparsity in arg[['s']]
two_step_single = function(n, p, svalue, bvalue, sig, lam1Types, alpha, ite, rawDir, 
                           lam2Coef, numCandidate, upper)
{
  arg = list()
  volume_stein_list = list()
  radius_stein_list = list()
  volume_lasso_list = list()
  radius_lasso_list = list()
  count = 1
  for (i in 1:length(lam1Types))
  {
    for (j in 1:length((svalue)))
    {
      for (k in 1:length(bvalue))
      {
        print(paste(lam1Types[i], svalue[j], bvalue[k]))
        arg[['root']] = rawDir
        arg[['beta']] = paste(rawDir, '/beta/', 'beta', '_s', svalue[j], '_b', bvalue[k], '.txt', sep = "")
        arg[['Y1']] = paste(rawDir, '/Y1/', 'Y1', '_s', svalue[j], '_b', bvalue[k], '.txt', sep = "")
        arg[['Y2']] = paste(rawDir, '/Y2/', 'Y2', '_s', svalue[j], '_b', bvalue[k], '.txt', sep = "")  
        arg[['lam1Type']] = lam1Types[i]
        if (arg[['lam1Type']] == '1se') {
            arg[['1se']] = paste(rawDir, '/1se/', '1se', '_s', svalue[j], '_b', bvalue[k], '.txt', 
                               sep = "")
        } else if (arg[['lam1Type']] == 'cv') {
            arg[['cv']] = paste(rawDir, '/cv/', 'cv', '_s', svalue[j], '_b', bvalue[k], '.txt', 
                               sep = "")
        }
        arg[['lam2Coef']] = lam2Coef
        arg[['s']] = svalue[j]
        arg[['numCandidate']] = numCandidate
        arg[['upper']] = upper

        # calculate necessary statistics
        tempRes = two_step_test(n, p, sig, alpha, ite, arg)
        volume_stein_df = tempRes$volume_stein_df
        radius_stein_df = tempRes$radius_stein_df
        volume_lasso_df = tempRes$volume_lasso_df
        radius_lasso_df = tempRes$radius_lasso_df
        
        volume_stein_df['lam1Type'] = lam1Types[i]
        volume_stein_df['s'] = svalue[j]
        volume_stein_df['b'] = bvalue[k]
        volume_stein_df['method'] = 'twoStepSteinVolume'

        radius_stein_df['lam1Type'] = lam1Types[i]
        radius_stein_df['s'] = svalue[j]
        radius_stein_df['b'] = bvalue[k]
        radius_stein_df['method'] = 'twoStepSteinRadius'

        volume_lasso_df['lam1Type'] = lam1Types[i]
        volume_lasso_df['s'] = svalue[j]
        volume_lasso_df['b'] = bvalue[k]
        volume_lasso_df['method'] = 'twoStepLassoVolume'

        radius_lasso_df['lam1Type'] = lam1Types[i]
        radius_lasso_df['s'] = svalue[j]
        radius_lasso_df['b'] = bvalue[k]
        radius_lasso_df['method'] = 'twoStepLassoRadius'
        
        volume_stein_list[[count]] = volume_stein_df
        radius_stein_list[[count]] = radius_stein_df
        volume_lasso_list[[count]] = volume_lasso_df
        radius_lasso_list[[count]] = radius_lasso_df
        count = count + 1
      }
    }
  }
  return(list(volume_stein_df = ldply(volume_stein_list),
              radius_stein_df = ldply(radius_stein_list),
              volume_lasso_df = ldply(volume_lasso_list),
              radius_lasso_df = ldply(radius_lasso_list)))
}


# set the oracle sparsity as the true sparsity in arg[['s']]
two_step_multiple = function(n, p, s_s, b_s, s_w, b_w, sig, lam1Types, alpha, ite, rawDir, 
                           lam2Coef, numCandidate, upper)
{
  arg = list()
  volume_stein_list = list()
  radius_stein_list = list()
  volume_lasso_list = list()
  radius_lasso_list = list()
  count = 1
  for (i in 1:length(lam1Types))
  {
    for (j in 1:length((b_s)))
    {
      for (k in 1:length(b_w))
      {
        print(paste(lam1Types[i], b_s[j], b_w[k]))
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
        arg[['lam2Coef']] = lam2Coef
        arg[['s']] = s_s + s_w
        arg[['numCandidate']] = numCandidate
        arg[['upper']] = upper
        # calculate necessary statistics
        tempRes = two_step_test(n, p, sig, alpha, ite, arg)
        volume_stein_df = tempRes$volume_stein_df
        radius_stein_df = tempRes$radius_stein_df
        volume_lasso_df = tempRes$volume_lasso_df
        radius_lasso_df = tempRes$radius_lasso_df

        volume_stein_df['lam1Type'] = lam1Types[i]
        volume_stein_df['bs'] = b_s[j]
        volume_stein_df['bw'] = b_w[k]
        volume_stein_df['method'] = 'twoStepSteinVolume'

        radius_stein_df['lam1Type'] = lam1Types[i]
        radius_stein_df['bs'] = b_s[j]
        radius_stein_df['bw'] = b_w[k]
        radius_stein_df['method'] = 'twoStepSteinRadius'

        volume_lasso_df['lam1Type'] = lam1Types[i]
        volume_lasso_df['bs'] = b_s[j]
        volume_lasso_df['bw'] = b_w[k]
        volume_lasso_df['method'] = 'twoStepLassoVolume'

        radius_lasso_df['lam1Type'] = lam1Types[i]
        radius_lasso_df['bs'] = b_s[j]
        radius_lasso_df['bw'] = b_w[k]
        radius_lasso_df['method'] = 'twoStepLassoRadius'

        volume_stein_list[[count]] = volume_stein_df
        radius_stein_list[[count]] = radius_stein_df
        volume_lasso_list[[count]] = volume_lasso_df
        radius_lasso_list[[count]] = radius_lasso_df
        count = count + 1
      }
    }
  }
  return(list(volume_stein_df = ldply(volume_stein_list),
              radius_stein_df = ldply(radius_stein_list),
              volume_lasso_df = ldply(volume_lasso_list),
              radius_lasso_df = ldply(radius_lasso_list)))
}

