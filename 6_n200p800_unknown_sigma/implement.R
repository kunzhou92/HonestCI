# ite
library(plyr)
source("data_fun.R")
source("utils.R")
source("adaptive_fun.R")
source("oracle_fun2.R")
source("two_step_fun.R")

N = 200
P = 800
SVALUE = c(10)
BVALUE = c(0.1, 0.3, 0.5, 0.7, 0.9, 1, 2, 3, 4, 5)

# S_S = 100
# B_S = c(0.1, 0.3, 0.5, 0.7, 0.9, 1, 2, 3, 4, 5)
# S_W = 100
# B_W = c(0.2)
S_S = 0
B_S = c()
S_W = 0
B_W = c()



ITE = 100
NFOLDS = 5
NUM_CANDIDATES = 40
SIG = 1
dir.create("data", showWarnings = F)
toe_dir = "data/toeplitz"
exp_dir = 'data/exp_decay'
equi_dir = 'data/equi_corr'
######################## data ###########################

# toeplitz
dir.create(toe_dir)
X_TYPE = 1
#gen_data(N, P, toe_dir, ITE, X_TYPE)
covariate_single(N, P, SVALUE, BVALUE, toe_dir, ITE, NFOLDS, sig = SIG)
covariate_multiple(N, P, S_S, B_S, S_W, B_W, toe_dir, ITE, NFOLDS, sig = SIG)

# exponential decay
dir.create(exp_dir)
X_TYPE = 2
#gen_data(N, P, exp_dir, ITE, X_TYPE)
covariate_single(N, P, SVALUE, BVALUE, exp_dir, ITE, NFOLDS, sig = SIG)
covariate_multiple(N, P, S_S, B_S, S_W, B_W, exp_dir, ITE, NFOLDS, sig = SIG)

# equi. corr.
dir.create(equi_dir)
X_TYPE = 3
#gen_data(N, P, equi_dir, ITE, X_TYPE)
covariate_single(N, P, SVALUE, BVALUE, equi_dir, ITE, NFOLDS, sig = SIG)
covariate_multiple(N, P, S_S, B_S, S_W, B_W, equi_dir, ITE, NFOLDS, sig = SIG)

hsigma_single(SVALUE, BVALUE, ITE, toe_dir)
hsigma_single(SVALUE, BVALUE, ITE, exp_dir)
hsigma_single(SVALUE, BVALUE, ITE, equi_dir)

hsigma_scaled_lasso(SVALUE, BVALUE, ITE, toe_dir)
hsigma_scaled_lasso(SVALUE, BVALUE, ITE, exp_dir)
hsigma_scaled_lasso(SVALUE, BVALUE, ITE, equi_dir)

############################ summary #####################
toe_dir = "data/toeplitz"
exp_dir = 'data/exp_decay'
equi_dir = 'data/equi_corr'
alpha = 0.05
NUM_CANDIDATES = 40
UPPER = 10
lam2Coef = 0.5

lam1Types = c('val', 'cv', '1se')
dir_pool = c(toe_dir, exp_dir, equi_dir)
design = c("toeplitz", "exp. decay", "equal cor.")

# lam1Types = c('cv')
# dir_pool = c(toe_dir)
# design = c("toeplitz")


# adaptive
res_single = list()
# res_multiple = list()
for (i in 1:length(dir_pool)) {
  temp = adaptive_single(N, P, SVALUE, BVALUE, sig = SIG, lam1Types, alpha, ITE, dir_pool[i])
  temp['design'] = design[i]
  res_single[[i]] = temp

  # temp = adaptive_multiple(N, P, S_S, B_S, S_W, B_W, sig = SIG, lam1Types, alpha, ITE, dir_pool[i])
  # temp['design'] = design[i]
  # res_multiple[[i]] = temp
}
summary_adaptive_single = ldply(res_single)
# summary_adaptive_multiple = ldply(res_multiple)


# two step
volume_stein_single = list()
radius_stein_single = list()
volume_lasso_single = list()
radius_lasso_single = list()

# volume_stein_multiple = list()
# radius_stein_multiple = list()
# volume_lasso_multiple = list()
# radius_lasso_multiple = list()
for (i in 1:length(dir_pool)) {
  temp = two_step_single(N, P, SVALUE, BVALUE, sig = SIG, lam1Types, alpha, ITE, dir_pool[i], 
                         lam2Coef, NUM_CANDIDATES, UPPER)
  temp[['volume_stein_df']]['design'] = design[i]
  temp[['radius_stein_df']]['design'] = design[i]
  temp[['volume_lasso_df']]['design'] = design[i]
  temp[['radius_lasso_df']]['design'] = design[i]
  volume_stein_single[[i]] = temp[['volume_stein_df']]
  radius_stein_single[[i]] = temp[['radius_stein_df']]
  volume_lasso_single[[i]] = temp[['volume_lasso_df']]
  radius_lasso_single[[i]] = temp[['radius_lasso_df']]

  # temp = two_step_multiple(N, P, S_S, B_S, S_W, B_W, sig = SIG, lam1Types, alpha, ITE, dir_pool[i], 
  #                          lam2Coef, NUM_CANDIDATES, UPPER)
  # temp[['volume_stein_df']]['design'] = design[i]
  # temp[['radius_stein_df']]['design'] = design[i]
  # temp[['volume_lasso_df']]['design'] = design[i]
  # temp[['radius_lasso_df']]['design'] = design[i]
  # volume_stein_multiple[[i]] = temp[['volume_stein_df']]
  # radius_stein_multiple[[i]] = temp[['radius_stein_df']]
  # volume_lasso_multiple[[i]] = temp[['volume_lasso_df']]
  # radius_lasso_multiple[[i]] = temp[['radius_lasso_df']]
}
summary_volume_stein_single = ldply(volume_stein_single)
summary_radius_stein_single = ldply(radius_stein_single)
summary_volume_lasso_single = ldply(volume_lasso_single)
summary_radius_lasso_single = ldply(radius_lasso_single)

# summary_volume_stein_multiple = ldply(volume_stein_multiple)
# summary_radius_stein_multiple = ldply(radius_stein_multiple)
# summary_volume_lasso_multiple = ldply(volume_lasso_multiple)
# summary_radius_lasso_multiple = ldply(radius_lasso_multiple)

# oracle lasso
# N is twice as much a  adaptive and two-step for fariness
lam1Types_oracle = c('val')
res_single = list()
# res_multiple = list()
for (i in 1:length(dir_pool)) {
  temp = oracle_single(N, P, SVALUE, BVALUE, sig = SIG, lam1Types_oracle, alpha, ITE, dir_pool[i])
  temp['design'] = design[i]
  res_single[[i]] = temp

  # temp = oracle_multiple(2 * N, P, S_S, B_S, S_W, B_W, sig = SIG, lam1Types_oracle, alpha, ITE, dir_pool[i])
  # temp['design'] = design[i]
  # res_multiple[[i]] = temp
}
summary_oracle_single = ldply(res_single)
# summary_oracle_multiple = ldply(res_multiple)


res_single_cv = list()
res_single_1se = list()
# res_multiple_cv = list()
# res_multiple_1se = list()

for (i in 1:length(dir_pool)) {
  temp = adjust_lasso_single_version2(N, P, SVALUE, BVALUE, SIG, ITE, dir_pool[i],
                                      summary_oracle_single, factor = 1.02, lower = 0.9, 
                                      upper = 0.95, stop = 1000, nfolds = 5, design[i])
  temp$df_cv['design'] = design[i]
  temp$df_1se['design'] = design[i]
  res_single_cv[[i]] = temp$df_cv
  res_single_1se[[i]] = temp$df_1se  

  # temp = adjust_lasso_multiple_version2(2 * N, P, S_S, B_S, S_W, B_W, SIG, ITE, dir_pool[i],
  #                                       summary_oracle_multiple, factor = 1.02, lower = 0.9, 
  #                                       upper = 0.95, stop = 1000, nfolds = 5, design[i])
  # temp$df_cv['design'] = design[i]
  # temp$df_1se['design'] = design[i]
  # res_multiple_cv[[i]] = temp$df_cv
  # res_multiple_1se[[i]] = temp$df_1se 
  
}
summary_oracle_single_cv = ldply(res_single_cv)
summary_oracle_single_1se = ldply(res_single_1se)
# summary_oracle_multiple_cv = ldply(res_multiple_cv)
# summary_oracle_multiple_1se = ldply(res_multiple_1se)

summary_oracle_single_version2 = rbind(summary_oracle_single, summary_oracle_single_cv, summary_oracle_single_1se)
# summary_oracle_multiple_version2 = rbind(summary_oracle_multiple, summary_oracle_multiple_cv, summary_oracle_multiple_1se)

save_list = c("summary_adaptive_single", "summary_adaptive_multiple",
               "summary_oracle_single", "summary_oracle_multiple",
               "summary_volume_stein_single", "summary_radius_stein_single",
               "summary_volume_lasso_single", "summary_radius_lasso_single",
               "summary_volume_stein_multiple", "summary_radius_stein_multiple",
               "summary_volume_lasso_multiple", "summary_radius_lasso_multiple",
               "summary_oracle_single_version2", "summary_oracle_multiple_version2")
save(list = save_list, file = 'summary_non_sparse.Rdata')





save_list = c( "summary_volume_stein_single", "summary_radius_stein_single",
               "summary_volume_lasso_single", "summary_radius_lasso_single",
               "summary_volume_stein_multiple", "summary_radius_stein_multiple",
               "summary_volume_lasso_multiple", "summary_radius_lasso_multiple")
save(list = save_list, file = 'summary_two_step.Rdata')



save_list = c( "summary_volume_stein_single", "summary_radius_stein_single",
               "summary_volume_lasso_single", "summary_radius_lasso_single")
save(list = save_list, file = 'unknown_sigma_fixed_A_least_square.Rdata')


