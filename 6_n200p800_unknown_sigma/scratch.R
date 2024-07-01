X = matrix(rnorm(1000), ncol = 5)
beta = matrix(1:5, ncol = 1)

y = X %*% beta + rnorm(200, sd = 2)


output = glmnet(X, y, family='gaussian', intercept = F, lambda = 1e-3)


library(scalreg)

output2 = scalreg(X, y)

output2$hsigma


# c(alpha): center at 0, sig: 0.7 coverage=0.87 r_w=0.5782494
y
# c(alpha): center at 0, sig: 0.8 coverage=1 r_w=0.5769122
z 
# c(alpha): center at 0, sig: 0.9 coverage=1 r_w=0.5198486
a 
# c(alpha): center at 0, sig: 1 coverage=1 r_w=0.4360677
b 
# c(alpha): center at 0, sig: 1.1 coverage=0.73, r_w=0.4337906
c 
# c(alpha): center at 0, sig: 1.2 coverage=0.29 r_w=0.4702923
d 


# c(alpha): center at mu, sig: 1.0 coverage=0.96, r_w=0.4324172
e  
# c(alpha): center at mu, sig: 1.1 coverage=0.76, r_w=0.4349259
f 
# c(alpha): center at mu, sig: 1.2 coverage=0.34, r_w=0.5140029
g 
# c(alpha): center at mu, sig: 0.9 coverage=1 r_w=0.5227806
h 
# c(alpha): center at mu, sig: 0.8 coverage=1 r_w=0.6068316
i 
# c(alpha): center at mu, sig: 0.7 coverage=0.9 r_w=0.6129564
j  


# c(alpha): center at 0, pure STEIN, sig :1 coverage=1 r_w=0.6285648
t1 
# c(alpha): center at 0, pure STEIN, sig :1.1 coverage=0.97 r_w=0.7156935
t2 
# c(alpha): center at 0, pure STEIN, sig :1.2 coverage=0.81 r_w=0.8194488
t3 
# c(alpha): center at 0, pure STEIN, sig :1.3 coverage=0.6 r_w=0.880606
t4 
# c(alpha): center at 0, pure STEIN, sig :0.9 coverage=1 r_w=0.6393519
t5
# c(alpha): center at 0, pure STEIN, sig :0.8 coverage=1 r_w=0.6691045
t6 
# c(alpha): center at 0, pure STEIN, sig :0.7 coverage=0.96 r_w=0.6451372
t7 



############ Stein method #################################

result = rep(0, 5000)
radius = rep(0, 5000)
for (i in 1:5000) {
  df = 200
  n = 200
  alpha = 0.05 / 2
  sig = 1
  mu = matrix(0, nrow = 200, ncol = 1)
  Y_sample = mu + matrix(rnorm(n, sd = sig), ncol = 1)
  sig2 = 1.3
  res = SteinToZero(Y_sample, df, sig2)
  c_alpha = quantile_stein(mu, df, alpha, sig2, iteNum = 100)
  r_w = sqrt((c_alpha * sig2^2 / sqrt(n) + res$sure))
  
  result[i] = mean((mu - res$mu)^2) / (r_w * r_w) 
  radius[i] = r_w
  
}

mean(result < 1)
mean(radius)

# center at hmu, sig=0.9 -> coverage=1, r=0.5002384
# center at hmu, sig=1.0 -> coverage=1, 
# center at hmu, sig=1.1 -> coverage=0.9888, r=0.4020781
# center at hmu, sig=1.2 -> coverage=0.7916, r=0.4833993
# center at hmu, sig=1.3 -> coverage=0.2406, r=0.569993
# center at 0, sig=0.9 -> coverage=1, r=0.4981471
# center at 0, sig=1.0 -> coverage=0.9998, r=0.4020002
# center at 0, sig=1.1 -> coverage=0.9418, r=0.3987881
# center at 0, sig=1.2 -> coverage=0.4704, r=0.4341105
# center at 0, sig=1.3 -> coverage=0.044, r=0.4704693


#######################################################################
n = 200
p = 800
Xtype = 1
sig = 1
X1 = generateX(n, p, Xtype)
beta = generateBeta(p, 10, 3, type = 3)
Y1 = tcrossprod(X1, t(beta)) + matrix(rnorm(n, sd = sig), ncol = 1)

library(scalreg)

tmp = scalreg(X1, Y1)
tmp$hsigma

scalreg(X1, Y1, lam0 = 'univ')$hsigma
scalreg(X1, Y1, lam0 = sqrt(2 * log(p) / n))$hsigma

tmp = scalreg(X1, Y1, lam0 = 'quantile', LSE = T)


n = 200
p = 800
s = 10
Xtype = 1
sig = 1



source("two_step_fun.R")
volume_stein_list = list()

for (i in 1:50) {
  X1 = generateX(n, p, Xtype)
  X2 = generateX(n, p, Xtype)
  beta = generateBeta(type = 3, p, s = 10, b = 3)
  Y1 = tcrossprod(X1, t(beta)) + matrix(rnorm(n, sd = sig), ncol = 1)
  Y2 = tcrossprod(X2, t(beta)) + matrix(rnorm(n, sd = sig), ncol = 1)
  
  # lam = 2 * sqrt(2 * log(p) / n)
  lam = cv.glmnet(X1, Y1, family='gaussian', intercept = F)$lambda.min
  # lam = cv.glmnet(X1, Y1, family='gaussian', intercept = F)$lambda.1se
  
  res = two_step(X1, X2, Y1, Y2, lam, 1, s, sig, 30, 0.05, 10, 100)
  volume_stein = convert(res[['volume_stein']], X2, beta)
  volume_stein_list[[i]] = volume_stein
}

volume_stein_df = ldply(volume_stein_list, data.frame)
mean(volume_stein_df$coverage)
mean(volume_stein_df$r_s)
mean(volume_stein_df$r_w)
mean(volume_stein_df$k)

# sig: 1; coverage: 1; r_s: 1.061; r_w: 0.596; k: 12.14
# sig: 1.1; coverage: 0.94; r_s: 1.116; r_w: 0.689; k: 10.62
# sig: 1.2; coverage: 0.86; r_s: 1.156; r_w: 0.841; k: 9.14
# sig: 0.8; coverage: 1; r_s: 0.855219; r_w: 0.6427057; k: 12.54



# adaptive sig:1.0 coverage=0.96 r=0.9217437
ada1
# adaptive sig:1.1 coverage=0.79 r=0.8125239
ada2
# adaptive sig:1.2 coverage=0.34 r=0.7044679
ada3
# adaptive sig:0.9 coverage=1 r=1.009321
ada4
# adaptive sig:0.8 coverage=1 r=1.078587
ada5


# adaptive
res_single = list()
for (i in 1:length(dir_pool)) {
  temp = adaptive_single(N, P, SVALUE, BVALUE, sig = 0, lam1Types, alpha, ITE, dir_pool[i])
  temp['design'] = design[i]
  res_single[[i]] = temp

}
ada6 = ldply(res_single)


# lam=val, sig=1.1, A=fixed, cov=0.73
# lam=1se, sig=1.1, A=fixed, cov=0.94
# lam=1se, sig=1.2, A=fixed, cov=0.47
# lam=val, sig=1.2, A=fixed, cov=0.26
# lam=1se, sig=1.1, A=data, cov=0.68
# lam=1se, sig=1.2, A=data, cov=0.25

# lam=1se, sig=0.8, A=fixed, cov=1
# lam=1se, sig=0.7, A=fixed, cov=0.79
# lam=1se, sig=1.2, A=cv, cov=0.61
