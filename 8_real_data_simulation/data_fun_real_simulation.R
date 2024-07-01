library(MASS)
library(glmnet)

# to sav coding work, it simply copies X1, X2 for ite times
gen_data_random_splitting = function(n, p, dir, ite, X, n1)
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

    index = sample(1:n, n1)
    X1 = X[index, ]
    X2 = X[-index, ]
    
    write.table(X1, fileX1, row.names = F, col.names = F)
    write.table(X2, fileX2, row.names = F, col.names = F)
  }
}


# to sav coding work, it simply copies X1, X2 for ite times
gen_data_clustering = function(n, p, dir, ite, X1, X2)
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
    write.table(X1, fileX1, row.names = F, col.names = F)
    write.table(X2, fileX2, row.names = F, col.names = F)
  }
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



# store corresponding beta, Y, lambda
covariate_single_clustering = function(n1, n2, p, svalue, bvalue, dir, ite, n_folds, sig = 1)
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
      Y1Mat = matrix(0, nrow = ite, ncol = n1)
      Y2Mat = matrix(0, nrow = ite, ncol = n2)
      lamCvMat = matrix(0, nrow = ite, ncol = 1)
      lam1seMat = matrix(0, nrow = ite, ncol = 1)
      
      for (i in 1:ite)
      {
        fileX1 = paste(dir, '/X1/', i, '.txt', sep = "") 
        fileX2 = paste(dir, '/X2/', i, '.txt', sep = "") 
        X1 = as.matrix(read.table(fileX1))
        X2 = as.matrix(read.table(fileX2))
        beta = generateBeta(type = 3, p = p, s = svalue[j], b = bvalue[k])
        betaMat[i,] = t(beta)
        Y1 = tcrossprod(X1, t(beta)) + matrix(rnorm(n1, sd = sig), ncol = 1)
        Y1Mat[i,] = Y1
        Y2Mat[i,] = tcrossprod(X2, t(beta)) + matrix(rnorm(n2, sd = sig), ncol = 1)
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



# set the oracle sparsity as the true sparsity in arg[['s']]
hsigma_single = function(svalue, bvalue, ite, rawDir)
{
  
  dir.create(paste(rawDir, '/hsigma', sep=''), showWarnings = T)
  count = 1
  for (j in 1:length((svalue)))
  {
    for (k in 1:length(bvalue))
    {
      beta_file = paste(rawDir, '/beta/', 'beta', '_s', svalue[j], '_b', bvalue[k], '.txt', sep = "")
      Y1_file = paste(rawDir, '/Y1/', 'Y1', '_s', svalue[j], '_b', bvalue[k], '.txt', sep = "")
      Y2_file = paste(rawDir, '/Y2/', 'Y2', '_s', svalue[j], '_b', bvalue[k], '.txt', sep = "")  
      hsigma_file = paste(rawDir, '/hsigma/', 'hsigma', '_s', svalue[j], '_b', bvalue[k], '.txt', sep = "")
      
      beta_pool = t(as.matrix(read.table(beta_file)))
      Y1_pool = t(as.matrix(read.table(Y1_file)))
      Y2_pool = t(as.matrix(read.table(Y2_file)))
      
      hsigma = matrix(0, nrow = ite, ncol = 1)
      for( i in 1:ite)
      {
        fileX1 = paste(rawDir, '/X1/', i, '.txt', sep = "") 
        X1 = as.matrix(read.table(fileX1))
        Y1 = Y1_pool[, i, drop = F]
        hsigma[i,1] = scalreg(X1, Y1, LSE = T)$lse$hsigma
      }
      write.table(hsigma, hsigma_file, row.names = F, col.names = F)
    }
  }
}




# set the oracle sparsity as the true sparsity in arg[['s']]
hsigma_scaled_lasso = function(svalue, bvalue, ite, rawDir)
{
  
  dir.create(paste(rawDir, '/hsigma_scaled_lasso', sep=''), showWarnings = T)
  count = 1
  for (j in 1:length((svalue)))
  {
    for (k in 1:length(bvalue))
    {
      beta_file = paste(rawDir, '/beta/', 'beta', '_s', svalue[j], '_b', bvalue[k], '.txt', sep = "")
      Y1_file = paste(rawDir, '/Y1/', 'Y1', '_s', svalue[j], '_b', bvalue[k], '.txt', sep = "")
      Y2_file = paste(rawDir, '/Y2/', 'Y2', '_s', svalue[j], '_b', bvalue[k], '.txt', sep = "")  
      hsigma_file = paste(rawDir, '/hsigma_scaled_lasso/', 'hsigma', '_s', svalue[j], '_b', bvalue[k], '.txt', sep = "")
      
      beta_pool = t(as.matrix(read.table(beta_file)))
      Y1_pool = t(as.matrix(read.table(Y1_file)))
      Y2_pool = t(as.matrix(read.table(Y2_file)))
      
      hsigma = matrix(0, nrow = ite, ncol = 1)
      for( i in 1:ite)
      {
        fileX1 = paste(rawDir, '/X1/', i, '.txt', sep = "") 
        X1 = as.matrix(read.table(fileX1))
        Y1 = Y1_pool[, i, drop = F]
        hsigma[i,1] = scalreg(X1, Y1)$hsigma
      }
      write.table(hsigma, hsigma_file, row.names = F, col.names = F)
    }
  }
}



# the name of each file _s_ -> _s
# utils readlambda 
# gen_data, covariate_single, covariate_multiple



