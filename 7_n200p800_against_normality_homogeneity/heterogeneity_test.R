# X and beta
X1 = generateX(n = 200, p = 800, type = 1)
X2 = generateX(n = 200, p = 800, type = 1)
beta = generateBeta(type = 3, p = 800, s = 10, b = 5)
mu1 = as.vector(tcrossprod(X1, t(beta)))
mu2 = as.vector(tcrossprod(X2, t(beta)))

# mu and sigma
lo = quantile(c(mu1, mu2), 0.025)
hi = quantile(c(mu1, mu2), 0.975)
incre = 15
coef = incre / (hi - lo)^2

sig = 1
sigma_var1 = sig + sig * coef * (mu1 - lo)^2
sigma_var1[mu1 < lo] = sig
# sigma_var1[hi < mu1] = sig + incre * sig

sigma_var2 = sig + sig * coef * (mu2 - lo)^2
sigma_var2[mu2 < lo] = sig
# sigma_var2[hi < mu2] = sig + incre * sig

coef = 2 / max_size # 50% perturbation
sig1 = 1 + coef * mu1

# Y
Y1 = tcrossprod(X1, t(beta)) + matrix(sigma_var1 * rnorm(200), ncol = 1)
Y2 = tcrossprod(X2, t(beta)) + matrix(sigma_var2 * rnorm(200), ncol = 1)

# plot 
plot(mu1, sigma_var1)
plot(mu2, sigma_var2)
plot(mu1, sigma_var1 * rnorm(200))
plot(mu2, sigma_var2 * rnorm(200))

# two step
lam1 = cv.glmnet(X1, Y1, intercept = F)$lambda.1se
hsig = scalreg(X1, Y1, LSE = T)$lse$hsigma
res = two_step(X1, X2, Y1, Y2, lam1, lam2Coef = 0, s = 10, sig = hsig,
               numCandidate = 40, alpha = 0.05, upper = 10, ite = 100)
  
res2 = convert(res$volume_stein, X2, beta)
res2$hsigma
res2$coverage

# adaptive
res3 = adaptive(X1, Y1, X2, Y2, beta, lam1, sig = hsig, alpha = 0.05)
res3$coverage


# test its residual under linear models

A2 = which(beta != 0 )
res3 = lm(Y1~X1[,A2] - 1)
plot(res3$fitted.values, res3$residuals)




mu = mu1
max_size = max(abs(mu1))
# max_size = max(abs(quantile(mu, c(0.025, 0.975)))) 
coef = 0.5 / max_size # 50% perturbation
sig = 1
Y1 = tcrossprod(X1, t(beta)) + matrix(sig * (1 + coef * mu1) * rnorm(200), ncol = 1)
res = scalreg(X1, Y1, LSE = T)
res$lse$hsigma

veps = rnorm(200)
scaled_veps = (1 + coef * mu1) * veps
plot(mu1, veps)
points(mu1, scaled_veps, col = 'red')


plot(mu1, scaled_veps)


hbeta = matrix(res$coefficients, ncol = 1)

x_axis = as.numeric(X1 %*% hbeta)
y_axis = as.numeric(Y1 - X1 %*% hbeta)
plot(x_axis, y_axis)

hist((1 + coef * mu1), breaks = 20)


res2 = cv.glmnet(X1, Y1, intercept = F)
res3 = glmnet(X1, Y1, lambda = res2$lambda.1se, intercept = F)
A = which(res3$beta != 0)
X1_sub = X1[,A]
res4 = lm(Y1~X1_sub - 1)
plot(res4$fitted.values, res4$residuals)

x_axis2 = as.numeric(X1 %*% res3$beta)
y_axis2 = as.numeric(Y1 - X1 %*% res3$beta)
plot(x_axis2, y_axis2)


# test generated data
X1 = as.matrix(read.table("data/toeplitz_coef_0.5/X1/10.txt"))
Y1 = t(as.matrix(read.table("data/toeplitz_coef_0.5/Y1/Y1_s10_b3.txt")[10,]))
beta = t(as.matrix(read.table("data/toeplitz_coef_0.5/beta/beta_s10_b3.txt")[10,]))
# true A
A2 = which(beta != 0 )
res5 = lm(Y1~X1[,A2] - 1)
plot(res5$fitted.values, res5$residuals)

quantile(mu1, 0.975)
