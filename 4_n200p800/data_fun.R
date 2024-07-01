library(MASS)
library(glmnet)

# generate design matrix
# n: value
# p: value
# type: 1 toeplitz; 2 exponential decay; 3 equi. corr 
# return: matrix: a n*p scaled matrix. 
generateX = function(n, p, type)
{
  mu = rep(0, p)
  if (type == 1)
  {
    firstRow = 0.5^(0:(p-1))
    Sigma = toeplitz(firstRow)
  }
  if (type == 2)
  {
    firstRow = 0.4^(0:(p-1))
    Sigma = solve(toeplitz(firstRow))
  }
  if (type == 3)
  {
    Sigma = matrix(0.8, ncol = p, nrow = p)
    diag(Sigma) = rep(1, p)
  }
  X = colScale(mvrnorm(n, mu = mu, Sigma = Sigma))
  return(X)
}


# generate beta
# b: strong signal; b2: weak signal
# type: 1 first s, 2 spaced evenly, 3 random position , 4 combined signals
# return: matrix: a p*1 matrix
generateBeta = function(type, p, s, b, s2=NULL, b2=NULL)
{
  beta = rep(0, p)
  if (type == 1)
    beta[1:s] = runif(s, -b, b)
  else if (type == 2)
    beta[round(seq(from = 1, to = p, length.out = s))] = runif(s, -b, b)
  else if (type == 3)
    beta[sample(1:p, size = s, replace = F)] = runif(s, -b, b)
  else if (type == 4)
  {
  	pos_tot = sample(1:p, size = s + s2, replace = F)
  	pos1 = sample(pos_tot, size = s, replace = F)
  	pos2 = pos_tot[!(pos_tot %in% pos1)]
  	beta[pos1] = runif(s, -b, b)
  	beta[pos2] = runif(s2, -b2, b2)
  }
  return(matrix(beta, ncol = 1))
}


# store "ite" number of matrices
gen_data = function(n, p, dir, ite, Xtype)
{
  pathX1 = file.path(dir, "X1")
  pathX2 = file.path(dir, "X2")
  dir.create(pathX1, showWarnings = F)
  dir.create(pathX2, showWarnings = F)
  for (i in 1:ite)
  {
    name = paste(i, ".txt", sep="")
    fileX1 = file.path(pathX1, name)
    fileX2 = file.path(pathX2, name)
    X1 = generateX(n, p, Xtype)
    X2 = generateX(n, p, Xtype)
    write.table(X1, fileX1, row.names = F, col.names = F)
    write.table(X2, fileX2, row.names = F, col.names = F)
  }
}

# store corresponding beta, Y, lambda
covariate_single = function(n, p, svalue, bvalue, dir, ite, n_folds)
{
  pathBeta = file.path(dir, "beta")
  pathY1 = file.path(dir, "Y1")
  pathY2 = file.path(dir, "Y2")
  pathCv = file.path(dir, "cv")
  path1se = file.path(dir, "1se")
  dir.create(pathBeta, showWarnings = F)
  dir.create(pathY1, showWarnings = F)
  dir.create(pathY2, showWarnings = F)
  dir.create(pathCv, showWarnings = F)
  dir.create(path1se, showWarnings = F)
  for (j in 1:length(svalue))
  {
    for (k in 1:length(bvalue))
    {
      fileBeta = paste(pathBeta, '/beta', '_s', svalue[j], '_b', bvalue[k], '.txt', sep="")
      fileY1 = paste(pathY1, '/Y1', '_s', svalue[j], '_b', bvalue[k], '.txt', sep="")
      fileY2 = paste(pathY2, '/Y2', '_s', svalue[j], '_b', bvalue[k], '.txt', sep="")
      fileCv = paste(pathCv, '/cv', '_s', svalue[j], '_b', bvalue[k], '.txt', sep="")
      file1se = paste(path1se, '/1se', '_s', svalue[j], '_b', bvalue[k], '.txt', sep="")
      
      betaMat = matrix(0, nrow = ite, ncol = p)
      Y1Mat = matrix(0, nrow = ite, ncol = n)
      Y2Mat = matrix(0, nrow = ite, ncol = n)
      lamCvMat = matrix(0, nrow = ite, ncol = 1)
      lam1seMat = matrix(0, nrow = ite, ncol = 1)
      
      for (i in 1:ite)
      {
      	fileX1 = paste(dir, '/X1/', i, '.txt', sep = "") 
        fileX2 = paste(dir, '/X2/', i, '.txt', sep = "") 
        X1 = as.matrix(read.table(fileX1))
        X2 = as.matrix(read.table(fileX2))
        beta = generateBeta(p, svalue[j], bvalue[k], type = 3)
        betaMat[i,] = t(beta)
        Y1 = tcrossprod(X1, t(beta)) + matrix(rnorm(n), ncol = 1)
        Y1Mat[i,] = Y1
        Y2Mat[i,] = tcrossprod(X2, t(beta)) + matrix(rnorm(n), ncol = 1)
        cv_glmnet = cv.glmnet(X1, Y1, nfolds = n_folds, intercept = F)
        lamCvMat[i,] = cv_glmnet$lambda.min
        lam1seMat[i,] = cv_glmnet$lambda.1se
      }
     
      write.table(betaMat, fileBeta, row.names = F, col.names = F)
      write.table(Y1Mat, fileY1, row.names = F, col.names = F)
      write.table(Y2Mat, fileY2, row.names = F, col.names = F)
      write.table(lamCvMat, fileCv, row.names = F, col.names = F)
      write.table(lam1seMat, file1se, row.names = F, col.names = F)
    }
  }
}

covariate_multiple = function(n, p, s_s, b_s, s_w, b_w, dir, ite, n_folds)
{
  pathBeta = file.path(dir, "beta")
  pathY1 = file.path(dir, "Y1")
  pathY2 = file.path(dir, "Y2")
  pathCv = file.path(dir, "cv")
  path1se = file.path(dir, "1se")
  dir.create(pathBeta,showWarnings = F)
  dir.create(pathY1, showWarnings = F)
  dir.create(pathY2, showWarnings = F)
  dir.create(pathCv, showWarnings = F)
  dir.create(path1se, showWarnings = F)
  for (j in 1:length(b_s))
  {
    for (k in 1:length(b_w))
    {
      fileBeta = paste(pathBeta, '/beta', '_ss', s_s, '_bs', b_s[j], '_sw', s_w, '_bw', b_w[k], '.txt', sep="")
      fileY1 = paste(pathY1, '/Y1', '_ss', s_s, '_bs', b_s[j], '_sw', s_w, '_bw', b_w[k], '.txt', sep="")
      fileY2 = paste(pathY2, '/Y2', '_ss', s_s, '_bs', b_s[j], '_sw', s_w, '_bw', b_w[k], '.txt', sep="")
      fileCv = paste(pathCv, '/cv', '_ss', s_s, '_bs', b_s[j], '_sw', s_w, '_bw', b_w[k], '.txt', sep="")
      file1se = paste(path1se, '/1se', '_ss', s_s, '_bs', b_s[j], '_sw', s_w, '_bw', b_w[k], '.txt', sep="")
      
      betaMat = matrix(0, nrow = ite, ncol = p)
      Y1Mat = matrix(0, nrow = ite, ncol = n)
      Y2Mat = matrix(0, nrow = ite, ncol = n)
      lamCvMat = matrix(0, nrow = ite, ncol = 1)
      lam1seMat = matrix(0, nrow = ite, ncol = 1)
      
      for (i in 1:ite)
      {
      	fileX1 = paste(dir, '/X1/', i, '.txt', sep = "") 
        fileX2 = paste(dir, '/X2/', i, '.txt', sep = "") 
        X1 = as.matrix(read.table(fileX1))
        X2 = as.matrix(read.table(fileX2))
        beta = generateBeta(p, s_s, b_s[j], type = 4, s_w, b_w[k])
        betaMat[i,] = t(beta)
        Y1 = tcrossprod(X1, t(beta)) + matrix(rnorm(n), ncol = 1)
        Y1Mat[i,] = Y1
        Y2Mat[i,] = tcrossprod(X2, t(beta)) + matrix(rnorm(n), ncol = 1)
        cv_glmnet = cv.glmnet(X1, Y1, nfolds = n_folds, intercept = F)
        lamCvMat[i,] = cv_glmnet$lambda.min
        lam1seMat[i,] = cv_glmnet$lambda.1se
      }
      write.table(betaMat, fileBeta, row.names = F, col.names = F)
      write.table(Y1Mat, fileY1, row.names = F, col.names = F)
      write.table(Y2Mat, fileY2, row.names = F, col.names = F)
      write.table(lamCvMat, fileCv, row.names = F, col.names = F)
      write.table(lam1seMat, file1se, row.names = F, col.names = F)
    }
  }
}

# the name of each file _s_ -> _s
# utils readlambda 
# gen_data, covariate_single, covariate_multiple