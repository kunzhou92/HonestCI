source("utils.R")
source("data_fun.R")

N = 200
P = 800
SVALUE = c(5, 10)
BVALUE = c(0.1, 0.3, 0.5, 0.7, 0.9, 1, 2, 3, 4, 5)
S_S = 5
B_S = c(0.1, 0.3, 0.5, 0.7, 0.9, 1, 2, 3, 4, 5)
S_W = 5
B_W = c(0.1, 0.2)
ITE = 2
NFOLDS = 5
RATIOS = seq(from = 0, to = 4, by = 0.05)
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


# B_W = c(0.1, 0.2)
# S_S S_W

# ite
