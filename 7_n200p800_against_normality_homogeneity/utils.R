library(matrixStats)

# normalize matrix columnwise
# x: matrix
# return: scaled matrix
colScale = function(x) {
  cm = colMeans(x, na.rm = TRUE)
  csd = colSds(x, center = cm)
  x = t( (t(x) - cm) / csd )
  return(x)
}

# squared l2-norm
# x: numeric vector
# return: the square of l2-norm of x
L2NormSquare = function(x)
{
  return(crossprod(x)[1,1])
}

# Calculate baseline
# n: degree of freedom
# sig: standard deviation
# return: baseline
baseline = function(n, sig, percent=0.95)
{
  return (sqrt(qchisq(percent, df = n) / n) * sig)
}

###################### Later #################################


# read particular lambda
# n: number of rows of X
# p: number of columns of X
# it



readLambda = function(n, p, ite, s, b, lamType, rawDir = NULL)
{
  # cv
  if (lamType == "cv")
  {
    
    path_cv = paste(rawDir, "/cv/","cv", "_s_", s, "_b_", b, ".txt", sep = "")
    lam = as.numeric(read.table(path_cv, skip = (ite - 1), nrows = 1))
  }
  # 1se
  else if (lamType == "1se")
  {
    path_1se = paste(rawDir, "/1se/","1se", "_s_", s, "_b_", b, ".txt", sep = "")
    lam = as.numeric(read.table(path_1se, skip = (ite - 1), nrows = 1))
  }
  # theoretical value
  else
  {
    lam = lamType * 2 * sqrt(2) * sqrt(log(p) / n)
  }
  return(lam)
}

####################################################
# return the  label of lam
# lamType: numeric
# return: string, description
####################################################
lamTypeLabel = function(lamType)
{
  if (lamType == 1)
  {
    label = 'const'
  } else if (lamType == 2)
  {
    label = '0.8_const' 
  } else if (lamType == 3) 
  {
    label = '0.2_const'
  } else if (lamType == 4)
  {
    label = 'cv' 
  } else if (lamType == 5)
  {
    label = '1se'
  }
  return(label)
}
