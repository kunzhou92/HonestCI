# ite
library(plyr)
source("data_fun.R")
source("utils.R")
source("adaptive_fun.R")
source("oracle_fun.R")
source("two_step_fun.R")


N = 200
P = 800
SVALUE = c(5, 10)
BVALUE = c(0.1, 0.3, 0.5, 0.7, 0.9, 1, 2, 3, 4, 5)
S_S = 5
B_S = c(0.1, 0.3, 0.5, 0.7, 0.9, 1, 2, 3, 4, 5)
S_W = 5
B_W = c(0.1, 0.2)
ITE = 100
NFOLDS = 5
RATIOS = seq(from = 0, to = 4, by = 0.05)


######################## data ###########################

dir.create("data", showWarnings = F)
# toeplitz
toe_dir = "data/toeplitz"
dir.create(toe_dir)
X_TYPE = 1
gen_data(N, P, toe_dir, ITE, X_TYPE)
covariate_single(N, P, SVALUE, BVALUE, toe_dir, ITE, NFOLDS)
covariate_multiple(N, P, S_S, B_S, S_W, B_W, toe_dir, ITE, NFOLDS)

# exponential decay
exp_dir = 'data/exp_decay'
dir.create(exp_dir)
X_TYPE = 2
gen_data(N, P, exp_dir, ITE, X_TYPE)
covariate_single(N, P, SVALUE, BVALUE, exp_dir, ITE, NFOLDS)
covariate_multiple(N, P, S_S, B_S, S_W, B_W, exp_dir, ITE, NFOLDS)

# equi. corr.
equi_dir = 'data/equi_corr'
dir.create(equi_dir)
X_TYPE = 3
gen_data(N, P, equi_dir, ITE, X_TYPE)
covariate_single(N, P, SVALUE, BVALUE, equi_dir, ITE, NFOLDS)
covariate_multiple(N, P, S_S, B_S, S_W, B_W, equi_dir, ITE, NFOLDS)



############################ summary #####################
toe_dir = "data/toeplitz"
exp_dir = 'data/exp_decay'
equi_dir = 'data/equi_corr'
lam1Types = c('val', 'cv', '1se')
alpha = 0.05
dir_pool = c(toe_dir, exp_dir, equi_dir)
design = c("toeplitz", "exp. decay", "equal cor.")
RATIOS = seq(from = 0, to = 4, by = 0.05)
UPPER = 10
lam2Coef = 0.5


# adaptive
res_single = list()
res_multiple = list()
for (i in 1:length(dir_pool)) {
  temp = adaptive_single(N, P, SVALUE, BVALUE, sig = 1, lam1Types, alpha, ITE, dir_pool[i])
  temp['design'] = design[i]
  res_single[[i]] = temp

  temp = adaptive_multiple(N, P, S_S, B_S, S_W, B_W, sig = 1, lam1Types, alpha, ITE, dir_pool[i])
  temp['design'] = design[i]
  res_multiple[[i]] = temp
}
summary_adaptive_single = ldply(res_single)
summary_adaptive_multiple = ldply(res_multiple)


# two step
volume_stein_single = list()
radius_stein_single = list()
volume_lasso_single = list()
radius_lasso_single = list()

volume_stein_multiple = list()
radius_stein_multiple = list()
volume_lasso_multiple = list()
radius_lasso_multiple = list()
for (i in 1:length(dir_pool)) {
  temp = two_step_single(N, P, SVALUE, BVALUE, sig = 1, lam1Types, alpha, ITE, dir_pool[i], 
                         lam2Coef, RATIOS, UPPER)
  temp[['volume_stein_df']]['design'] = design[i]
  temp[['radius_stein_df']]['design'] = design[i]
  temp[['volume_lasso_df']]['design'] = design[i]
  temp[['radius_lasso_df']]['design'] = design[i]
  volume_stein_single[[i]] = temp[['volume_stein_df']]
  radius_stein_single[[i]] = temp[['radius_stein_df']]
  volume_lasso_single[[i]] = temp[['volume_lasso_df']]
  radius_lasso_single[[i]] = temp[['radius_lasso_df']]

  temp = two_step_multiple(N, P, S_S, B_S, S_W, B_W, sig = 1, lam1Types, alpha, ITE, toe_dir, 
                           lam2Coef, RATIOS, UPPER)
  temp[['volume_stein_df']]['design'] = design[i]
  temp[['radius_stein_df']]['design'] = design[i]
  temp[['volume_lasso_df']]['design'] = design[i]
  temp[['radius_lasso_df']]['design'] = design[i]
  volume_stein_multiple[[i]] = temp[['volume_stein_df']]
  radius_stein_multiple[[i]] = temp[['radius_stein_df']]
  volume_lasso_multiple[[i]] = temp[['volume_lasso_df']]
  radius_lasso_multiple[[i]] = temp[['radius_lasso_df']]
}
summary_volume_stein_single = ldply(volume_stein_single)
summary_radius_stein_single = ldply(radius_stein_single)
summary_volume_lasso_single = ldply(volume_lasso_single)
summary_radius_lasso_single = ldply(radius_lasso_single)

summary_volume_stein_multiple = ldply(volume_stein_multiple)
summary_radius_stein_multiple = ldply(radius_stein_multiple)
summary_volume_lasso_multiple = ldply(volume_lasso_multiple)
summary_radius_lasso_multiple = ldply(radius_lasso_multiple)

# oracle lasso
# N is twice as much as adaptive and two-step for fariness
lam1Types_oracle = c('val')
res_single = list()
res_multiple = list()
for (i in 1:length(dir_pool)) {
  temp = oracle_single(2 * N, P, SVALUE, BVALUE, sig = 1, lam1Types_oracle, alpha, ITE, i)
  temp['design'] = design[i]
  res_single[[i]] = temp

  temp = oracle_multiple(2 * N, P, S_S, B_S, S_W, B_W, sig = 1, lam1Types_oracle, alpha, ITE, i)
  temp['design'] = design[i]
  res_multiple[[i]] = temp
}
summary_oracle_single = ldply(res_single)
summary_oracle_multiple = ldply(res_multiple)

save_list = c("summary_adaptive_single", "summary_adaptive_multiple",
              "summary_oracle_single", "summary_oracle_multiple",
              "summary_volume_stein_single", "summary_radius_stein_single",
              "summary_volume_lasso_single", "summary_radius_lasso_single",
              "summary_volume_stein_multiple", "summary_radius_stein_multiple",
              "summary_volume_lasso_multiple", "summary_radius_lasso_multiple")

save(list = save_list, file = 'summary.Rdata')