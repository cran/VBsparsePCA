#' Function for the PX-CAVI algorithm using the multivariate normal slab
#'
#'
#' This function employs the PX-CAVI algorithm proposed in Ning (2020).
#' The \eqn{g} in the slab density of the spike and slab prior is chosen to be the multivariate normal distribution, i.e.,
#' \eqn{N(0, \sigma^2/\lambda_1 I_r)}.
#' Details of the model and the prior can be found in the Details section in the description of the `VBsparsePCA()` function.
#'
#'
#'@param x Data an \eqn{n*p} matrix.
#'@param r Rank.
#'@param lambda Tuning parameter for the density \eqn{g}.
#'@param max.iter The maximum number of iterations for running the algorithm.
#'@param eps The convergence threshold; the default is \eqn{10^{-4}}.
#'@param jointly.row.sparse The default is true, which means that the jointly row sparsity assumption is used; one could not use this assumptio by changing it to false.
#'@param sig2.true The default is false, \eqn{\sigma^2} will be estimated; if sig2 is known and its value is given, then \eqn{\sigma^2} will not be estimated.
#'@param threshold The threshold to determine whether \eqn{\gamma_j} is 0 or 1; the default value is 0.5.
#'@param theta.int The initial value of theta mean; if not provided, the algorithm will estimate it using PCA.
#'@param theta.var.int The initial value of theta.var; if not provided, the algorithm will set it to be 1e-3*diag(r).
#'@param kappa.para1 The value of \eqn{\alpha_1} of \eqn{\pi(\kappa)}; default is 1.
#'@param kappa.para2 The value of \eqn{\alpha_2} of \eqn{\pi(\kappa)}; default is \eqn{p+1}.
#'@param sigma.a The value of \eqn{\sigma_a} of \eqn{\pi(\sigma^2)}; default is 1.
#'@param sigma.b The value of \eqn{\sigma_b} of \eqn{\pi(\sigma^2)}; default is 2.
#'
#'
#'@return \item{iter}{The number of iterations to reach convergence.}
#' \item{selection}{A vector (if \eqn{r = 1} or with the jointly row-sparsity assumption) or a matrix (if otherwise) containing the estimated value for \eqn{\boldsymbol \gamma}.}
#' \item{theta.mean}{The loadings matrix.}
#' \item{theta.var}{The covariance of each non-zero rows in the loadings matrix.}
#' \item{sig2}{Variance of the noise.}
#' \item{obj.fn}{A vector contains the value of the objective function of each iteration. It can be used to check whether the algorithm converges}
#'
#'@examples
#' #In this example, the first 20 rows in the loadings matrix are nonzero, the rank is 1
#' set.seed(2021)
#' library(MASS)
#' library(pracma)
#' n <- 200
#' p <- 1000
#' s <- 20
#' r <- 1
#' sig2 <- 0.1
#' # generate eigenvectors
#' U.s <- randortho(s, type = c("orthonormal"))
#' U <- rep(0, p)
#' U[1:s] <- as.vector(U.s[, 1:r])
#' s.star <- rep(0, p)
#' s.star[1:s] <- 1
#' eigenvalue <- seq(20, 10, length.out = r)
#' # generate Sigma
#' theta.true <- U * sqrt(eigenvalue)
#' Sigma <- tcrossprod(theta.true) + sig2*diag(p)
#' # generate n*p dataset
#' X <- t(mvrnorm(n, mu = rep(0, p), Sigma = Sigma))
#' result <- spca.cavi.mvn(x = X, r = 1)
#' loadings <- result$theta.mean
#' @export

spca.cavi.mvn <- function(
  x, r, lambda = 1, max.iter = 100, eps = 1e-4,
  jointly.row.sparse = TRUE, sig2.true = NA,
  threshold = 0.5, theta.int = NA, theta.var.int = NA,
  kappa.para1 = NA, kappa.para2 = NA, sigma.a = NA,
  sigma.b = NA
  ) {

  #------ Dimension of input matrices ------#
  n <- dim(x)[1]
  p <- dim(x)[2]

  # initialize theta.hat and theta.var
  svd.res <- svd(x)
  if (is.na(theta.int) == TRUE) {
    if (r == 1) {
      theta.mean <- as.matrix((t(svd.res$u) %*% x)[1, ]) / sqrt(n-1)
    } else {
      theta.mean <- t((t(svd.res$u) %*% x)[1:r, ]) / sqrt(n-1)
    }
  } else {
    theta.mean <- theta.int
  }

  if (is.na(theta.var.int) == TRUE) {
    theta.var <- 1e-3*diag(r)
  } else {
    theta.var <- theta.var.int
  }

  # initialize z
  if (jointly.row.sparse == F) {
    if (r == 1) {
      z <- rep(1,p)
    } else {
      z <- matrix(1, p, r)
    }
  } else {
    z <- rep(1, p)
  }

  # initialize sigma - choose sig2 to be the second smallest eigenvalue of the covariance matrix
  if (is.na(sig2.true) == TRUE) {
    sig2 <- (svd.res$d[length(svd.res$d)-1])^2 / (n-1)
  } else {
    sig2 <- sig2.true
  }

  # choose hyperparameter for priors
  # para1 for the prior of theta
  if (is.na(kappa.para1) == TRUE) {
    kappa.para1 <- 1
  } else {
    kappa.para1 <- kappa.para1
  }

  # para2 for the prior of theta
  if (is.na(kappa.para2) == TRUE) {
    kappa.para2 <- p+1
  } else {
    kappa.para2 <- kappa.para2
  }

  # para1 for the prior of sigma
  if (is.na(sigma.a) == TRUE) {
    sigma.a <- 1
  } else {
    sigma.a <- sigma.a
  }

  # para2 for the prior of sigma
  if (is.na(sigma.b) == TRUE) {
    sigma.b <- 2
  } else {
    sigma.b <- sigma.b
  }

  # create empty matrix to collect results
  theta.res <- array(NA, dim = c(p, r, max.iter))
  theta.var.res <- rep(NA, max.iter)
  selection.res <- array(NA, dim = c(p,r, max.iter))
  obj.fn.res <- matrix(NA, max.iter)
  sig2.res <- rep(NA, max.iter)

  # starting the loop
  for (iter in 1:max.iter) {

    ########################## E-step: update q(w) #############################
    # sum.theta.var <- Reduce("+", theta.var)
    if (iter == 1) {
      theta.mean.old <- theta.mean
      z.hat.old <- z.theta <- z
    }

    z.theta.mean <- theta.mean * z.theta # obtain theta.mean * z

    if (jointly.row.sparse == TRUE) {
      var.w <- solve(crossprod(z.theta.mean)/sig2 + sum(z.theta)*(theta.var) + diag(r)) # var of w
    } else {
      z.theta.var <- diag(theta.var) * z.theta # obtain z * theta.var
      if (r == 1) {
        var.w <- solve(crossprod(z.theta.mean)/sig2 + sum(z.theta.var) + diag(r))
      } else {
        var.w <- solve(crossprod(z.theta.mean)/sig2 + diag(colSums(z.theta.var)) + diag(r)) # var of w
      }
    }

    # mean of w
    mean.w <- sapply(1:n, FUN = function(i) {
      x[i,] %*% ((z.theta.mean) %*% var.w) / sig2
    })

    # MSE of w
    MSE.w.i <- lapply(1:n, FUN = function(i) {
      if (r == 1) {
        (mean.w[i]^2 + var.w)
      } else {
        (tcrossprod(mean.w[, i]) + var.w)
      }
    })
    MSE.w <- Reduce("+", MSE.w.i) # sum of n second moments of w_i

    ########################## update q(beta) #############################
    ## Mean of beta
    if (r == 1) {
      beta.term1 <- (MSE.w + lambda)
      beta.term2 <- t(mean.w %*% x)
      beta.mean <- as.matrix(c(beta.term2) / c(beta.term1))
      beta.var <- 1/beta.term1
    } else {
      beta.term1 <- (MSE.w + lambda * diag(r))
      beta.term2 <- t(mean.w %*% x)
      beta.var <- solve(beta.term1)
      beta.mean <- beta.term2 %*% beta.var
    }

    ########################### update z ################################
    h.hat <- sapply(1:p, FUN = function(j) {
      if (r == 1) {

        term1 <- - 1 /(sig2) * sum((mean.w %*% x[, j]) * beta.mean[j, ])
        term2 <- 1 / (2*sig2) * beta.mean[j, ]^2 * c(MSE.w)
        term3.1 <- sig2 * beta.var * c(MSE.w)
        term3 <- term3.1  / (2*sig2)
        term4 <- lambda / 2 * (sig2 * beta.var + (beta.mean[j])^2)
        term5 <-  - r * log (lambda) / 2 - log(beta.var)/2 -
          1/2 + log(kappa.para2/kappa.para1)

      } else {

        if (jointly.row.sparse == TRUE) {
          term1 <- - 1/(sig2) * sum((mean.w %*% x[, j]) * beta.mean[j, ])
          term2 <- c(1/(2*sig2) * beta.mean[j, ] %*% MSE.w %*% beta.mean[j, ] )
          term3.1 <- sig2 * beta.var %*% MSE.w
          term3 <- sum(diag(term3.1)) / (2*sig2)
          term4 <- lambda/2 * (sig2 * sum(diag(beta.var)) + crossprod(beta.mean[j, ]))
          term5 <- - r * log (lambda) / 2 - log(det(beta.var))/2 -
            1/2  + log(kappa.para2/kappa.para1)

        } else {
          # if jointly.row.sparse = False, then each term is a r*1 vector
          term1 <- - 1/(sig2) * (mean.w %*% x[, j]) * beta.mean[j, ]
          term2.1 <- tcrossprod(beta.mean[j, ]) %*% MSE.w
          term2 <- diag( term2.1 )/(2*sig2)
          term3.1 <- sig2 * beta.var %*% MSE.w
          term3 <- diag(term3.1) / (2*sig2)
          term4 <- lambda/2 * (sig2 * diag(beta.var) + beta.mean[j, ]^2)
          term5 <- - log(lambda)/2 - log(diag(beta.var))/2 - 1/2 +
             log(kappa.para2/kappa.para1)
        }
      }
      sum <- -(c(term1) + c(term2) + c(term3) + c(term4) + c(term5))
    })

    if (jointly.row.sparse == F) {
      if (r == 1) {
        z.hat <- pracma::sigmoid(h.hat)
      } else {
        z.hat <- pracma::sigmoid(t(h.hat))
      }
    } else {
      z.hat <- pracma::sigmoid(h.hat)
    }
    z <- ifelse(z.hat < threshold, 0, 1) # the estimated z for beta
    z.beta.mean <- beta.mean * z

    ######################## update sig2 ###########################
    if (is.na(sig2.true) == TRUE) {

      # evaluate the sum for each j
      sigma.foreach.j <- sapply(1:p, FUN = function(j) {
        sum(x[, j]^2) - 2*sum(mean.w %*% x[, j] * z.beta.mean[j, ]) +
          c(z.beta.mean[j, ] %*% MSE.w %*% z.beta.mean[j, ]) +
          lambda*sum(z.beta.mean[j, ]^2)
      })

      sig2 <- (sum(sigma.foreach.j) + 2*sigma.b) / (p*n + 2*sigma.a+2)
    }

    #################### Rotate beta back to theta ###################
    # first rotation obtain D
    D <- MSE.w/n
    chol.D <- t(chol(D))
    beta.mean.tilde <- beta.mean %*% chol.D
    theta.var <- chol.D %*% beta.var %*% t(chol.D) # theta.var = beta.var

    # second rotation obtain A
    beta.mean.tilde.svd <- svd(beta.mean.tilde)
    A <- beta.mean.tilde.svd$u # obtain A
    if (r == 1) {
      theta.mean <- A %*% beta.mean.tilde.svd$d
    } else {
      theta.mean <- A %*% diag(beta.mean.tilde.svd$d)
    }

    # when r > 1 and jointly.row.sparsity == F, then the support of theta
    # can be different from it of beta; an extra step is needed to obtain the
    # support of theta
    ####################### obtain sparsity for theta #######################
    if (jointly.row.sparse == F && r > 1) {
      h.theta.hat <- sapply(1:p, FUN = function(j) {
        # if not jointly.row.sparse, then each term is an r*1 vector
        # one need to check if the sign of theta.mean is the same as it of beta.mean
        if (sign(theta.mean[j, 1]) == sign(beta.mean[j, 1])) {
          term1 <- - 1/(sig2) * (mean.w %*% x[, j]) * theta.mean[j, ]
        } else {
          term1 <- 1/(sig2) * (mean.w %*% x[, j]) * theta.mean[j, ]
        }
        term2.1 <- tcrossprod(theta.mean[j, ]) %*% MSE.w
        term2 <- diag(term2.1)/(2*sig2)
        term3.1 <- sig2 * theta.var %*% MSE.w
        term3 <- diag(term3.1) / (2*sig2)
        term4 <- lambda/2 * diag(sig2 * theta.var + tcrossprod(theta.mean[j, ]))
        term5 <- - log(lambda)/2 - log(diag(theta.var))/2 - 1/2 +
          log(kappa.para2/kappa.para1)
        sum <- -(c(term1) + c(term2) + c(term3) + c(term4) + c(term5))
      })
      z.theta.hat <- pracma::sigmoid(t(h.theta.hat))
      z.theta <- ifelse(z.theta.hat < threshold, 0, 1)
    } else {
      z.theta <- z
    }

    ################# obtain the objective function #####################
    obj.fn1 <- sum((tcrossprod(theta.mean) - tcrossprod(theta.mean.old))^2)
    obj.fn2 <- sum(abs(z.hat - z.hat.old))
    obj.fn <- max(obj.fn1, obj.fn2)

    ######################### collect results ############################
    if (jointly.row.sparse == TRUE) {
      selection.res[, 1,iter] <- z.theta
    } else {
      selection.res[, , iter] <- z.theta
    }
    theta.res[, , iter] <- theta.mean * z.theta
    obj.fn.res[iter] <- obj.fn
    sig2.res[iter] <- sig2

    if (obj.fn < eps) {
      break
    } else {
      theta.mean.old <- theta.mean
      z.hat.old <- z.hat
    }
  }

  if (jointly.row.sparse == TRUE) {
    return(list(iter = iter, selection = selection.res[,1,iter],
                theta.mean = theta.res[,,iter], theta.var  = theta.var,
                sig2 = sig2, obj.fn = obj.fn.res[1:iter]))
  } else {
    return(list(iter = iter, selection = selection.res[,,iter],
                theta.mean = theta.res[,,iter], theta.var  = theta.var,
                sig2 = sig2, obj.fn = obj.fn.res[1:iter]))
  }
}
