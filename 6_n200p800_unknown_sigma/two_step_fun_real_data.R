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

two_step = function(X1, X2, Y1, Y2, lam1, numCandidate, alpha, upper, ite, sig = NULL)
{
  n1 = dim(X1)[1] 
  n2 = dim(X2)[1]
  p = dim(X1)[2]
  
  ################################# New part #############################################
  # estimate unknwon sigma
  sig = scalreg(X1, Y1, LSE = T)$lse$hsigma
  ########################################################################################

  # candidates
  lasso_split = glmnet(X1, Y1, family='gaussian', intercept = F, lambda = lam1)
  ################################# new part ###########################
  # res_candidate = threshold(as.matrix(lasso_split$beta), numCandidate, X2)
  # res_candidate$count = length(res_candidate$candidates)
  
  candidates = list()
  candidates[[1]] = which(lasso_split$beta != 0)
  # candidates[[1]] = integer()
  res_candidate = list(candidates = candidates, count = 1)
  ###########################################################################
  # stein c(alpha)

  A = which(lasso_split$beta != 0)
  qr_decom = qr(X2[,A])
  k = qr_decom$rank
  if (k == 0) {
    P_perp = diag(rep(1, n2))
  } else {
    Q = qr.Q(qr_decom)
    P_A = B2P(Q[,1:k])
    P_perp = diag(rep(1, n2)) - P_A
  }
  coef = sqrt(L2NormSquare(P_perp %*% Y2))
  
  c_alpha_stein = quantile_stein2(n2, alpha / 2, sig, ite, coef)

  # statistics
  res_stein_pool = list()
  logVol_volume_stein = rep(0, res_candidate[['count']])
  logVol_radius_stein = rep(0, res_candidate[['count']])
 
  
  # record
  for (i in 1:res_candidate[['count']])
  {
    res_stein = weak_stein(X2, Y2, res_candidate[['candidates']][[i]], sig, alpha / 2, c_alpha_stein, upper)
    res_stein_pool[[i]] = res_stein
    logVol_volume_stein[i] = res_stein[['logVol_volume']]
    logVol_radius_stein[i] = res_stein[['logVol_radius']]
  }
  
  i_volume_stein = which.min(logVol_volume_stein)
  i_radius_stein = which.min(logVol_radius_stein)
  
  volume_stein = list(r_s = res_stein_pool[[i_volume_stein]][['r_s_volume']], 
                      r_w = res_stein_pool[[i_volume_stein]][['r_w_volume']],
                      mu_s = res_stein_pool[[i_volume_stein]][['mu_s']],
                      mu_w = res_stein_pool[[i_volume_stein]][['mu_w']],
                      logVol = res_stein_pool[[i_volume_stein]][['logVol_volume']],
                      A = res_candidate[['candidates']][[i_volume_stein]],
                      k = length(res_candidate[['candidates']][[i_volume_stein]]),
                      hsigma = sig)
  
  radius_stein = list(r_s = res_stein_pool[[i_radius_stein]][['r_radius']],
                      r_w = res_stein_pool[[i_radius_stein]][['r_radius']],
                      mu_s = res_stein_pool[[i_radius_stein]][['mu_s']],
                      mu_w = res_stein_pool[[i_radius_stein]][['mu_w']],
                      logVol = res_stein_pool[[i_radius_stein]][['logVol_radius']],
                      A = res_candidate[['candidates']][[i_radius_stein]],
                      k = length(res_candidate[['candidates']][[i_radius_stein]]),
                      hsigma = sig)
  
 
  return(list(volume_stein = volume_stein, radius_stein = radius_stein))
}

adaptive = function(X1, Y1, X2, Y2, mu2, lam, sig, alpha)
{
  ####################################### new part ###################################
  sig = scalreg(X1, Y1, LSE = T)$lse$hsigma
  # print(sig)
  ##################################################################################  
  # calculate betaHat for first part
  lasso = glmnet(X1, Y1, family="gaussian", lambda=lam, intercept=F)
  betaHat = as.matrix(lasso$beta)
  k = sum(betaHat != 0)
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

  rv = L2NormSquare(mu2 - X2 %*% betaHat) / n2 / (r^2)
  coverage = rv <= 1
  logVol = n2 * log(r)
  
  result = list(r = r, coverage=coverage, logVol=logVol, k=k, hsigma = sig)
  return(result)
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
