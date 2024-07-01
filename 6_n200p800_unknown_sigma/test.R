source("two_step_fun_test.R")
source("utils.R")
source("adaptive_fun.R")
N = 200
P = 800
SVALUE = c(10)
BVALUE = c(3)


ITE = 100
NFOLDS = 5
NUM_CANDIDATES = 40
SIG = 1
toe_dir = "data/toeplitz"
exp_dir = 'data/exp_decay'
equi_dir = 'data/equi_corr'
alpha = 0.05
UPPER = 10

# lam1Types = c('val', 'cv', '1se')
# dir_pool = c(toe_dir, exp_dir, equi_dir)
# design = c("toeplitz", "exp. decay", "equal cor.")

lam1Types = c('1se')
dir_pool = c(toe_dir)
design = c("toeplitz")

volume_stein_single = list()
for (i in 1:length(dir_pool)) {
  volume_stein_single[[i]] = two_step_single(N, P, SVALUE, BVALUE, sig = 1, lam1Types, alpha, 
                                             ITE, dir_pool[i], NUM_CANDIDATES, UPPER)
}
summary_volume_stein_single = ldply(volume_stein_single)

source("two_step_fun_test.R")
result = list()
sig = c(1, 1.2)
for (i in 1:length(sig)) {
  
  result[[i]] = two_step_single(N, P, SVALUE, BVALUE, sig[i], '1se', alpha, 
                           ITE, toe_dir, NUM_CANDIDATES, 3)
  col = c("coverage", "r_s", "r_w", "k", "logVol", "hsigma", 
          "mu_s_norm", "mu_w_norm", "diff_s_norm", "diff_w_norm",
          "r_s_raw", "r_w_raw", "c_s", "c_w")
  print('================================================')
  print(colMeans(result[[i]][col]))
  print(mean(result[[i]]['cover_s'] < 1))
  print(mean(result[[i]]['cover_w'] < 1))
  
}

# theoretical analysis of sig and coverage \hmu and r_\perp
# nwec_alpha


result2 = list()
# sig = c(0.8, 0.9, 1, 1.1, 1.2, 100)
sig = c(2, 3, 4)
for (i in 1:length(sig)) {
  
  result2[[i]] = adaptive_single(N, P, SVALUE, BVALUE, sig[i], '1se', alpha, 
                                ITE, toe_dir)
  col = c("coverage", "r", "k", "logVol", "hsigma")
  print('================================================')
  print(colMeans(result2[[i]][col]))
  
}




quantile_stein = function(n, alpha, sig, iteNum, coef)
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

coef = 10 ^ seq(from = 0, to = 1, length.out = 20)
q = rep(0, 20)
for (i in 1:20){
  q[i] = quantile_stein(200, 0.05, 1, 10000, coef = coef[i])
}
plot(coef, q)


quantile_stein(200, 0.05, 1, 10000, coef = 1)
