# the main change is function "threshold" now decide the number of candidate sets we want and then
# evenly space from empty to supp(\beta).
# remove tau_chosen
# qr$qr is the combination of Q and R. use qr.Q and qr.R to recover

library(glmnet)
library(plyr)
library(LDRTools)
source("data_fun.R")
################################# New part #############################################
library(scalreg)
########################################################################################


# Shrink Y towards zero by Stein estimate
SteinToZero = function(Y, df, sig)
{
  B = sig^2 * df / L2NormSquare(Y)
  mu = max(1 - B, 0) * Y
  sure = sig^2 * max(1 - B, 0)

  return(list(mu = mu, sure = sure))
}

# calculate c(alpha) for Stein. (for squared radius)
quantile_stein = function(n, alpha, sig, iteNum)
{
  z = rep(0, iteNum)
  for (i in 1:iteNum)
  {
    Y_sample =  matrix(rnorm(n, sd = sig), ncol = 1)
    res = SteinToZero(Y_sample, n, sig)
    z[i] = abs(res$sure - L2NormSquare(res$mu) / n) * sqrt(n) / (sig^2)
  }
  c = quantile(z, probs = 1 - alpha)
  return (c)
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
    P_A = B2P(Q[,1:k])
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
    r_w = sqrt((n2 - k) / n2 * (c_alpha_w * sig^2 / sqrt(n2 - k) + abs(stein$sure)))
    mu_w = stein$mu
  }
  

  # update r_s_ratio r_w_ratio if A is not empty nor rank(X_A) == n2
  if (k == 0)
  {
    r_s_volume = 0
    r_w_volume = r_w
    logVol_volume = n2 * log(r_w_volume)
  } else if (k == n2)
  {
    r_s_volume = r_s
    r_w_volume = 0
    logVol_volume = n2 * log(r_s_volume)
  } else
  {
    # by volume
    lower = upper / (upper - 1) 
    c_s = max(lower, min(n2 / k, upper))
    c_w = c_s / (c_s - 1)
    r_s_volume = r_s * sqrt(c_s)
    r_w_volume = r_w * sqrt(c_w)
    logVol_volume = k * log(r_s_volume) + (n2 - k) * log(r_w_volume)
  }
  return(list(r_s_volume = r_s_volume,
              r_w_volume = r_w_volume,
              logVol_volume = logVol_volume,
              mu_s = mu_s,
              mu_w = mu_w))
}

# beta: numeric vector
threshold = function(beta, numCandidate) {
  nonzero_beta_indices = (1:length(beta))[beta != 0]
  beta = beta[beta != 0]
  size_of_candidates = unique(round(seq(from = 0, to = length(beta), length.out = numCandidate)))
  
  order_nonzero_beta = order(abs(beta), decreasing = TRUE)
  candidates = list()
  for (i in 1:length(size_of_candidates))
  {
    A = nonzero_beta_indices[order_nonzero_beta[0:size_of_candidates[i]]]
    candidates[[i]] = A
  }
  
  return(list(candidates = candidates,
              count = length(size_of_candidates)))
}


two_step = function(X1, X2, Y1, Y2, lam1, s, sig, numCandidate, alpha, upper, ite)
{
  n1 = dim(X1)[1] 
  n2 = dim(X2)[1]
  p = dim(X1)[2]
  ################################# sigma #############################################
  # scaled lasso
  # sig = scalreg(X1, Y1, LSE = T)$lse$hsigma
  # theoretical value
  ########################################################################################
  
  
  ############################## lambda ################################################
  lasso_split = glmnet(X1, Y1, family='gaussian', intercept = F, lambda = lam1)
  #########################################################################################
  
  ############################## candidate ################################################
  # multiple candidates
  res_candidate = threshold(as.vector(lasso_split$beta), numCandidate)
  # single candidate
  candidates = list()
  beta = as.numeric(lasso_split$beta)
  nonzero_beta_indices = (1:length(beta))[beta != 0]
  beta = beta[beta != 0]
  order_nonzero_beta = order(abs(beta), decreasing = TRUE)
  A = nonzero_beta_indices[order_nonzero_beta[0:min(8, length(nonzero_beta_indices))]]
  # A = which(lasso_split$beta != 0)
  candidates[[1]] = A
  res_candidate = list(candidates = candidates, count = 1)
  ########################################################################################
  
  
  # stein c(alpha)
  ######################################  new part ##################
  c_alpha_stein = quantile_stein(n2, alpha / 2, sig, ite)
  # fixed alpha
  # c_alpha_stein = 3.070031
  ##################################################################
  
  
  # statistics
  res_stein_pool = list()
  logVol_volume_stein = rep(0, res_candidate[['count']])

  # record
  for (i in 1:res_candidate[['count']])
  {
    res_stein = weak_stein(X2, Y2, res_candidate[['candidates']][[i]], sig, alpha / 2, c_alpha_stein, upper)
    res_stein_pool[[i]] = res_stein
    logVol_volume_stein[i] = res_stein[['logVol_volume']]
  }
  # i is chosen
  i_volume_stein = which.min(logVol_volume_stein)
  
  volume_stein = list(r_s = res_stein_pool[[i_volume_stein]][['r_s_volume']], 
                      r_w = res_stein_pool[[i_volume_stein]][['r_w_volume']],
                      mu_s = res_stein_pool[[i_volume_stein]][['mu_s']],
                      mu_w = res_stein_pool[[i_volume_stein]][['mu_w']],
                      logVol = res_stein_pool[[i_volume_stein]][['logVol_volume']],
                      A = res_candidate[['candidates']][[i_volume_stein]],
                      k = length(res_candidate[['candidates']][[i_volume_stein]]),
                      hsigma = sig)
  volume_stein
}

convert = function(m_list, X2, beta) {
  n2 = dim(X2)[1]
  p = dim(X2)[2]
  if (m_list[['k']] == 0) 
  {
    part1 = m_list[["mu_w"]] - X2 %*% beta
    part2 = diag(rep(1, n2)) / m_list[["r_w"]]^2
    
    mu_s_norm = 0
    mu_w_norm = L2NormSquare(X2 %*% beta) / n2
    diff_s_norm = 0
    diff_w_norm = L2NormSquare(m_list[["mu_w"]] - X2 %*% beta) / n2
    cover_s = NA
    cover_w = L2NormSquare(m_list[["mu_w"]] - X2 %*% beta) / n2 / m_list[["r_w"]]^2
    r_s_raw = NA
    r_w_raw = m_list[['r_w']]
    c_s = NA
    c_w = 1
    
  } else if (m_list[['k']] == n2)
  {
    part1 = m_list[["mu_s"]] - X2 %*% beta
    part2 = diag(rep(1, n2)) / m_list[["r_s"]]^2
    
    mu_s_norm = L2NormSquare(X2 %*% beta) / n2
    mu_w_norm = 0
    diff_s_norm = L2NormSquare(m_list[["mu_s"]] - X2 %*% beta) / n2
    diff_w_norm = 0
    cover_s = L2NormSquare(m_list[["mu_s"]] - X2 %*% beta) / n2 / m_list[["r_s"]]^2
    cover_w = NA
    r_s_raw = m_list[["r_s"]]
    r_w_raw = NA
    c_s = 1
    c_w = NA
  }
  else
  {
    
    ######################## new part #################
    qr_decom = qr(X2[,m_list[['A']]])
    k = qr_decom$rank
    Q = qr.Q(qr_decom)
    P_A = B2P(Q[,1:k])
    # k = length(m_list[['A']])
    # P_A = B2P(X2[,m_list[['A']]])
    #################################################
    
    P_perp = diag(rep(1, n2)) - P_A  
    part1 = m_list[["mu_s"]] + m_list[["mu_w"]] - X2 %*% beta
    part2 = P_A / m_list[["r_s"]]^2  + P_perp / m_list[["r_w"]]^2
    
    mu_s_norm = L2NormSquare(P_A %*% X2 %*% beta) / n2
    mu_w_norm = L2NormSquare(P_perp %*% X2 %*% beta) / n2
    diff_s_norm = L2NormSquare(m_list[["mu_s"]] - P_A %*% X2 %*% beta) / n2
    diff_w_norm = L2NormSquare(m_list[["mu_w"]] - P_perp %*% X2 %*% beta) / n2
    c_s = NA
    c_w = 1
    
    
    ###################### tmp ##################################
    upper = 3
    lower = upper / (upper - 1) 
    c_s = max(lower, min(n2 / k, upper))
    c_w = c_s / (c_s - 1)
    # print(c_s)
    # print(c_w)
    cover_s = L2NormSquare(m_list[["mu_s"]] - P_A %*% X2 %*% beta) / n2 / m_list[["r_s"]]^2 * c_s
    cover_w = L2NormSquare(m_list[["mu_w"]] - P_perp %*% X2 %*% beta) / n2 / m_list[["r_w"]]^2 * c_w
    r_s_raw = m_list[['r_s']] / sqrt(c_s)
    r_w_raw = m_list[['r_w']] / sqrt(c_w)
  }
  temp = (t(part1) %*% part2 %*% part1)/ n2 
  list(coverage = as.numeric(temp) <= 1,
       r_s = m_list[['r_s']],
       r_w = m_list[['r_w']],
       k = m_list[['k']],
       logVol = m_list[['logVol']],
       hsigma = m_list[['hsigma']],
       mu_s_norm = mu_s_norm,
       mu_w_norm = mu_w_norm,
       cover_s = cover_s,
       cover_w = cover_w,
       diff_s_norm = diff_s_norm,
       diff_w_norm = diff_w_norm,
       r_s_raw = r_s_raw,
       r_w_raw = r_w_raw,
       c_s = c_s,
       c_w = c_w,
       A = do.call(paste, as.list(as.character(m_list[['A']])))
       )
         
}

# ite for computing quantile is set here
two_step_test = function(n, p, sig, alpha, ite, arg)
{
  ite_quantile = 1000
  volume_stein_pool = list()

  beta_pool = t(as.matrix(read.table(arg[['beta']])))
  Y1_pool = t(as.matrix(read.table(arg[['Y1']])))
  Y2_pool = t(as.matrix(read.table(arg[['Y2']]))) 
  ################################## new #########################################
  # hsigma = as.matrix(read.table(arg[['hsigma']]))
  ################################################################################

  if (arg[['lam1Type']] == '1se') {
      lam1 = read.table(arg[['1se']])[,1]
  } else if (arg[['lam1Type']] == 'cv') {
      lam1 = read.table(arg[['cv']])[,1]
  } else if (arg[['lam1Type']] == 'val') {
      lam1 = rep(2 * sqrt(2) * sqrt(log(p) / n), ite)
  }
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
    ######################### new part #####################################
    # X2 = X1
    # Y2 = Y1
    # sig = hsigma[i,]
    ###########################################################################
    
    # necessary statistics
    res = two_step(X1, X2, Y1, Y2, lam1[i], s, sig, numCandidate, alpha, upper, ite_quantile)
    volume_stein = convert(res, X2, beta)
    volume_stein_pool[[i]] = volume_stein
  }

  volume_stein_df = ldply(volume_stein_pool, data.frame)
  volume_stein_df
}

# set the oracle sparsity as the true sparsity in arg[['s']]
two_step_single = function(n, p, svalue, bvalue, sig, lam1Types, alpha, ite, rawDir, 
                           numCandidate, upper)
{
  arg = list()
  volume_stein_list = list()
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
        arg[['s']] = svalue[j]
        arg[['numCandidate']] = numCandidate
        arg[['upper']] = upper
        ############################# new ##############################
        # arg[['hsigma']] = paste(rawDir, '/hsigma_scaled_lasso/', 'hsigma', '_s', svalue[j], '_b', bvalue[k], '.txt', sep = "")
        ####################################################################
        # calculate necessary statistics
        volume_stein_df = two_step_test(n, p, sig, alpha, ite, arg)

        volume_stein_df['lam1Type'] = lam1Types[i]
        volume_stein_df['s'] = svalue[j]
        volume_stein_df['b'] = bvalue[k]
        volume_stein_df['method'] = 'twoStepSteinVolume'
        
        volume_stein_list[[count]] = volume_stein_df
        count = count + 1
      }
    }
  }
  ldply(volume_stein_list)
}


quantile_stein2 = function(n, alpha, sig, iteNum, coef)
{
  z = rep(0, iteNum)
  for (i in 1:iteNum)
  {
    mu = rnorm(n)
    mu = coef * matrix(mu / sqrt(L2NormSquare(mu)), ncol = 1)
    
    Y_sample =  mu + matrix(rnorm(n, sd = sig), ncol = 1)
    res = SteinToZero(Y_sample, n, sig)
    z[i] = abs(res$sure - L2NormSquare(res$mu - mu) / n) * sqrt(n) / (sig^2)
  }
  c = quantile(z, probs = 1 - alpha)
  return (c)
}
