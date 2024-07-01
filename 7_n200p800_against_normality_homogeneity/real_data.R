library(plyr)
library(hdi)
library(scalreg)
source("data_fun.R")
source("utils.R")
source("adaptive_fun.R")
source("oracle_fun2.R")
source("two_step_fun_real_data.R")

# read data

data("riboflavin")
x = riboflavin$x
y = riboflavin$y


# clustering
dist = as.dist(1 - abs(cor(t(x))))

model = hclust(dist)
den = as.dendrogram(model)
plot(den, horiz = T, leaflab  = 'none', xlab = expression(paste("1-|", rho, "|")))


sub_grp = cutree(model, k = 2)
table(sub_grp)


library(glmnet)
# set.seed(42)
# fit.cv <- cv.glmnet(x = x, y = y)
# b <- as.matrix(coef(fit.cv))
# 
# rownames(b)[b != 0]
# sum(b != 0)

source("two_step_fun_real_data.R")
x = colScale(x)
y = y - mean(y)

X1 = x[sub_grp == 1, ]
Y1 = matrix(y[sub_grp == 1], ncol = 1)
X2 = x[sub_grp == 2, ]
Y2 = matrix(y[sub_grp == 2], ncol = 1)

# single cnadidate set
# numCandidate has no effect now
cv_model = cv.glmnet(X1, Y1, intercept = F)
lam1 = cv_model$lambda.1se
tmp = two_step(X1, X2, Y1, Y2, lam1, numCandidate = 0, alpha = 0.05, upper = 3, ite = 10000)
a = tmp$volume_stein

cv_model = cv.glmnet(X2, Y2, intercept = F)
lam2 = cv_model$lambda.1se
tmp = two_step(X2, X1, Y2, Y1, lam2, numCandidate = 0, alpha = 0.05, upper = 3, ite = 10000)
b = tmp$volume_stein
  
n1 = dim(X1)[1]
n2 = dim(X2)[1]
a$r_s
a$r_w
a$A
a$hsigma
L2NormSquare(a$mu_s + a$mu_w - Y2) / n2
yhat = as.vector(a$mu_s + a$mu_w)

b$r_s
b$r_w

res = scalreg(X1, Y1, LSE = T)

########################## band ############################################


CI = function(mu_i, mat, i) {
  diff = matrix(0, nrow = dim(mat)[1], ncol = 1)
  diff[i,1] = mu_i
  as.numeric(t(diff) %*% mat %*% diff) - 1
}

theme2 = theme_bw() + 
  theme(panel.grid = element_blank(), 
        text=element_text(size=14), 
        plot.title = element_text(hjust = 0.5),
        legend.background=element_blank(),
        legend.key=element_blank(), 
        legend.key.height=unit(1, 'cm'),
        strip.background = element_blank(),
        legend.position = 'bottom',
        strip.placement = "outside")

#################################### confidence band 1 ############################
hmu = a$mu_s + a$mu_w

qr_decom = qr(X2[,a$A])
k = qr_decom$rank
Q = qr.Q(qr_decom)
P_A = B2P(Q[,1:k])
P_perp = diag(rep(1, dim(X2)[1])) - P_A  
mat = P_A / a$r_s^2  + P_perp / a$r_w^2

n2 = dim(X2)[1]
left = rep(0, n2)
right = rep(0, n2)
center = as.numeric(a$mu_s + a$mu_w)
for (i in 1:n2) {
  left[i] = uniroot(CI, interval = c(-1, 0), mat = mat, i = i)$root
  right[i] = uniroot(CI, interval = c(0, 1), mat = mat, i = i)$root
}


df = data.frame(sample = factor(which(sub_grp == 2)), lower = center + left, upper = center + right, center = center,
                color = factor((center + left) <= 0 & 0 <= (center + right)))
library(ggplot2)

f1 = ggplot(df, aes(x = sample, color = color)) + geom_errorbar(aes(x = sample, ymin = lower, ymax = upper)) +
  geom_point(aes(y = center)) + 
  geom_hline(yintercept = 0, linetype = 3) +
  guides(color = FALSE) + ylab('confidence band') + xlab('individual index') +
  ggtitle("(a) group 1") + theme2

###################################### confidence band 2 ##############################

hmu = b$mu_s + b$mu_w

qr_decom = qr(X1[,b$A])
k = qr_decom$rank
Q = qr.Q(qr_decom)
P_A = B2P(Q[,1:k])
P_perp = diag(rep(1, dim(X1)[1])) - P_A  
mat = P_A / b$r_s^2  + P_perp / b$r_w^2

n1 = dim(X1)[1]
left = rep(0, n1)
right = rep(0, n1)
center = as.numeric(b$mu_s + b$mu_w)
for (i in 1:n1) {
  left[i] = uniroot(CI, interval = c(-1, 0), mat = mat, i = i)$root
  right[i] = uniroot(CI, interval = c(0, 1), mat = mat, i = i)$root
}


df2 = data.frame(sample = factor(which(sub_grp == 1)), lower = center + left, upper = center + right, center = center,
                color = factor((center + left) <= 0 & 0 <= (center + right)))

f2 = ggplot(df2[1:22,], aes(x = sample, color = color)) + geom_errorbar(aes(x = sample, ymin = lower, ymax = upper)) +
  geom_point(aes(y = center)) + 
  geom_hline(yintercept = 0, linetype = 3) +
  guides(color = FALSE) + ylab('confidence band') + xlab('individual index') +
  ggtitle("(a) group 2, part 1") + 
  theme2

f3 = ggplot(df2[23:44,], aes(x = sample, color = color)) + geom_errorbar(aes(x = sample, ymin = lower, ymax = upper)) +
  geom_point(aes(y = center)) + 
  geom_hline(yintercept = 0, linetype = 3) +
  guides(color = FALSE) + ylab('confidence band') + xlab('individual index') + 
  ggtitle("(b) group 2, part 2") + 
  theme2

library(gridBase)
library(gridExtra)
fig_final = grid.arrange(f1, f2, f3, nrow = 3, heights = c(1, 1, 1))

fig_final = grid.arrange(f2, f3, nrow = 2)


ggsave('confidence_band.pdf', fig_final, width = 6.67, height = 4.44, units = 'in')

#############################################################################

# normalize X and y
x = colScale(x)
y = y - mean(y)
sig = scalreg(x, y, LSE = T)$lse$hsigma
n = length(y)

source("two_step_fun_real_data.R")

# generate data
rep_int = 10
row = rep(1:n, rep_int)
X_large = x[row, ]
mu_large = matrix(y[row], ncol = 1)
# output
tot = 1
stein1 = list()
stein2 = list()
adapt1 = list()
adapt2 = list()
for (i in 1:tot){
  # index = which(sub_grp == 1)
  index = sample(1:n, ceiling(n / 2))
  index = rowSums(expand.grid(index, (0:(rep_int - 1)) * n))
  # generate data
  Y_large = mu_large + matrix(sig * rnorm(n * rep_int), ncol = 1)
  
  X1 = X_large[index, ]
  X2 = X_large[-index, ]
  Y1 = Y_large[index, ,drop=F]
  Y2 = Y_large[-index, ,drop=F]
  mu1 = mu_large[index, ,drop=F]
  mu2 = mu_large[-index, ,drop=F]
  
  adapt1[[i]] = adaptive_real(X1, X2, Y1, Y2, mu1, mu2, sig)
  adapt2[[i]] = adaptive_real(X2, X1, Y2, Y1, mu2, mu1, sig)
  stein1[[i]] = two_step_stein_real_data(X1, X2, Y1, Y2, mu1, mu2, sig)
  stein2[[i]] = two_step_stein_real_data(X2, X1, Y2, Y1, mu2, mu1, sig)
}

adapt_df1 = ldply(adapt1, data.frame)
adapt_df2 = ldply(adapt2, data.frame)
stein_df1 = ldply(stein1, data.frame)
stein_df2 = ldply(stein2, data.frame)

# colMeans(df1, na.rm = T)
# colMeans(df2, na.rm = T)


# twostep_unknownSigma_realSingle_A_groupSplit1 = df1
# twostep_unknownSigma_realSingle_A_groupSplit2 = df2
# twostep_unknownSigma_realSingle_A_evenSplit = rbind(df1, df2)
# adaptive_unknownSigma_groupSplit1 = df1
# adaptive_unknownSigma_groupSplit2 = df2
# adaptive_unknownSigma_evenSplit = rbind(df1, df2)



mean(df1$coverage)
mean(df1$r_s)
mean(df1$r_w)
# mean(df1$r)
mean(df1$k)
mean(df1$hsigma)
mean(df1$logVol)


mean(df2$coverage)
# mean(df2$r_s)
# mean(df2$r_w)
mean(df2$r)
mean(df2$k)
mean(df2$hsigma)
mean(df2$logVol)




# a wrap function
two_step_stein_real_data = function(X1, X2, Y1, Y2, mu1, mu2, sig) {
  
  # constants
  alpha = 0.05
  upper = 10
  ite = 400
  numCandidate = 20

  # lam1 setup
  lam1 = cv.glmnet(X1, Y1, intercept = F)$lambda.1se

  model = two_step(X1, X2, Y1, Y2, lam1, numCandidate, alpha, upper, ite, sig)$volume_stein
  
  n2 = dim(Y2)[1]
  if (model[['k']] == 0) 
  {
    part1 = model[["mu_w"]] - mu2
    part2 = diag(rep(1, n2)) / model[["r_w"]]^2
    
    mu_s_norm = 0
    mu_w_norm = L2NormSquare(model[['mu_w']]) / n2
    diff_s_norm = 0
    diff_w_norm = L2NormSquare(model[["mu_w"]] - mu2) / n2
    w_norm = L2NormSquare(mu2) / n2
    cover_s = NA
    cover_w = L2NormSquare(model[["mu_w"]] - mu2) / n2 / model[["r_w"]]^2
    r_s_raw = NA
    r_w_raw = model[['r_w']]
    
  } else if (model[['k']] >= n2)
  {
    part1 = model[["mu_s"]] - mu2
    part2 = diag(rep(1, n2)) / model[["r_s"]]^2
    
    mu_s_norm = L2NormSquare(model[['mu_s']]) / n2
    mu_w_norm = 0
    diff_s_norm = L2NormSquare(model[["mu_s"]] - mu2) / n2
    diff_w_norm = 0
    w_norm = 0
    cover_s = L2NormSquare(model[["mu_s"]] - mu2) / n2 / model[["r_s"]]^2
    cover_w = NA
    r_s_raw = model[["r_s"]]
    r_w_raw = NA
    
  } else
  {
    ########################## new #########################
    qr_decom = qr(X2[,model[['A']]])
    k = qr_decom$rank
    Q = qr.Q(qr_decom)
    P_A = B2P(Q[,1:k])
    #########################################################
   
    P_perp = diag(rep(1, n2)) - P_A  
    part1 = model[["mu_s"]] + model[["mu_w"]] - mu2
    part2 = P_A / model[["r_s"]]^2  + P_perp / model[["r_w"]]^2
    
    mu_s_norm = L2NormSquare(model[['mu_s']]) / n2
    mu_w_norm = L2NormSquare(model[['mu_w']]) / n2
    diff_s_norm = L2NormSquare(model[["mu_s"]] - P_A %*% mu2) / n2
    diff_w_norm = L2NormSquare(model[["mu_w"]] - P_perp %*% mu2) / n2
    w_norm = L2NormSquare(P_perp %*% mu2) / n2
    lower = upper / (upper - 1) 
    c_s = max(lower, min(n2 / k, upper))
    c_w = c_s / (c_s - 1)
    # print(c_s)
    # print(c_w)
    cover_s = L2NormSquare(model[["mu_s"]] - P_A %*% mu2) / n2 / model[["r_s"]]^2 * c_s
    cover_w = L2NormSquare(model[["mu_w"]] - P_perp %*% mu2) / n2 / model[["r_w"]]^2 * c_w
    r_s_raw = model[['r_s']] / sqrt(c_s)
    r_w_raw = model[['r_w']] / sqrt(c_w)
    
  }
  covered = (t(part1) %*% part2 %*% part1)/ n2 <= 1
  
  list(r_s = model[['r_s']], r_w = model[['r_w']], k = model[['k']],
       hsigma = model[['hsigma']], logVol = model[['logVol']], coverage = covered,
       cover_s = cover_s, cover_w = cover_w,
       diff_s_norm = diff_s_norm, diff_w_norm = diff_w_norm,
       r_s_raw = r_s_raw, r_w_raw = r_w_raw,
       w_norm = w_norm)
  
}

adaptive_real = function(X1, X2, Y1, Y2, mu1, mu2, sig) {

  # constants
  alpha = 0.05
  
  # lam1 setup
  lam1 = cv.glmnet(X1, Y1, intercept = F)$lambda.1se
  
  adaptive(X1, Y1, X2, Y2, mu2, lam1, sig, alpha)
}















# save_list = c( "twostep_unknownSigma_realSingle_A_even",
#                "twostep_unknownSigma_realSingle_A_groupSplit1",
#                "twostep_unknownSigma_realSingle_A_groupSplit2",
#                "adaptive_unknownSigma_evenSplit",
#                "adaptive_unknownSigma_groupSplit1",
#                "adaptive_unknownSigma_groupSplit2")
# save(list = save_list, file = 'real_data_simulation3.Rdata')



colMeans(twostep_knownSigma_realSingle_A_groupSplit1, na.rm = T)
colMeans(twostep_knownSigma_realSingle_A_groupSplit2, na.rm = T)

