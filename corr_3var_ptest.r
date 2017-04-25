
# Check whether the correlation difference of two variables to a third is
# significant
# Adapted from: https://stats.stackexchange.com/questions/151219/is-one-variable-more-correlated-than-another-with-a-third-in-r

cortest <- function(my.data, a, b, c, n.iter=500) {
  stat <- function(a, b, c) abs(cor(a, b) - cor(a, c)) # Test statistic function
  d <- scale(my.data[c(a, b, c)])             # Standardize all variables
  n <- dim(d)[1]                  # Find the number of observations
  x <- d[, 1]                     # Reference the first standardized variable
  s <- stat(x, d[,2], d[,3])      # Compute the test statistic
  sim <- replicate(n.iter, {      # Estimate the permutation distribution:
    i <- runif(n) < 1/2           #   Choose which values to swap
    y <- d[, 2]; y[i] <- d[i, 3]  #   Swap them
    z <- d[, 3]; z[i] <- d[i, 2]
    stat(x, y, z)                 #   Compute the test statistic
  })
  p.value <- sum(sim >= s) / s # Find the estimated p-value
  return(list(stat=s, p.value=p.value, sim=sim))
}