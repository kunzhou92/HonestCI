library(plyr)
library(hdi)
library(scalreg)
source("data_fun_real_simulation.R")
source("utils.R")
source("adaptive_fun.R")
source("two_step_fun.R")

# read data

data("riboflavin")
x = riboflavin$x
y = riboflavin$y

X = colScale(x)
Y = y - mean(y)

############################################## clustering ####################################
dist = as.dist(1 - abs(cor(t(x))))

model = hclust(dist)
den = as.dendrogram(model)
plot(den, horiz = T, leaflab  = 'none', xlab = expression(paste("1-|", rho, "|")))


sub_grp = cutree(model, k = 2)
table(sub_grp)

X1 = x[sub_grp == 1, ]
Y1 = matrix(y[sub_grp == 1], ncol = 1)
N1 = dim(X1)[1]
X2 = x[sub_grp == 2, ]
N2 = dim(X2)[1]
Y2 = matrix(y[sub_grp == 2], ncol = 1)

# choose lambda
# cv_model = cv.glmnet(x, y, intercept = F)
# lam = cv_model$lambda.1se

model = scalreg(X, Y, LSE = T)
sigma = model$lse$hsigma
s = sum(model$lse$coefficients != 0)
b = max(abs(model$lse$coefficients))
print(s)


P = dim(x)[2]
# SVALUE = c(9, 14)
# BVALUE = c(0.1, 0.3, 0.5, 0.7, 0.9, 1, 2, 3, 4, 5) * sigma
# ITE = 100
# lam1Types = c('val', 'cv', '1se')

BVALUE = 1 * sigma
ITE = 100
SVALUE = 14
lam1Types = c('cv')

NFOLDS = 5
NUM_CANDIDATES = 80
SIG = sigma

alpha = 0.05
lam2Coef = 0.5
UPPER = 10


dir.create("data", showWarnings = F)
cluster_dir = "data/cluster"

# X1:44; X2:27, should swap later
# cluster
dir.create(cluster_dir)
X_TYPE = 1
gen_data_clustering(NaN, P, cluster_dir, ITE, X1, X2)
covariate_single_clustering(N1, N2, P, SVALUE, BVALUE, cluster_dir, ITE, NFOLDS, sig = SIG)
hsigma_single(SVALUE, BVALUE, ITE, cluster_dir)
# hsigma_scaled_lasso(SVALUE, BVALUE, ITE, cluster_dir)


# adaptive method
summary_adaptive_single = adaptive_single(N1, N2, P, SVALUE, BVALUE, sig = SIG, lam1Types, alpha, ITE, cluster_dir)
summary_adaptive_single['design'] = 'clustering'

# two-step stein
temp = two_step_single(N1, N2, P, SVALUE, BVALUE, sig = SIG, lam1Types, alpha, ITE, cluster_dir, 
                       lam2Coef, NUM_CANDIDATES, UPPER)
temp[['volume_stein_df']]['design'] = 'clustering'
temp[['radius_stein_df']]['design'] = 'clustering'
temp[['volume_lasso_df']]['design'] = 'clustering'
temp[['radius_lasso_df']]['design'] = 'clustering'
summary_volume_stein_single = temp[['volume_stein_df']]
summary_radius_stein_single = temp[['radius_stein_df']]
summary_volume_lasso_single = temp[['volume_lasso_df']]
summary_radius_lasso_single = temp[['radius_lasso_df']]


save_list = c("summary_adaptive_single",
              "summary_volume_stein_single", "summary_radius_stein_single",
              "summary_volume_lasso_single", "summary_radius_lasso_single")
save(list = save_list, file = 'clustering_X44_X27.Rdata')



lasso_split = glmnet(X1, Y1, family='gaussian', intercept = F, lambda = 0.01)

candidates = list()
if (sum(lasso_split$beta != 0) <= n2) {
  candidates[[1]] = which(lasso_split$beta != 0)
} else {
  hbeta_tmp = as.numeric(lasso_split$beta)
  candidates[[1]] = order(hbeta_tmp, decreasing = T)[1:n2]
}
