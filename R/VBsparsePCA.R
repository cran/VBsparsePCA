#' The main function for the variational Bayesian method for sparse PCA
#'
#'
#' This function employs the PX-CAVI algorithm proposed in Ning (2021).
#' The method uses the sparse spiked-covariance model and the spike and slab prior (see below).
#' Two different slab densities can be used: independent Laplace densities and a multivariate normal density.
#' In Ning (2021), it recommends choosing the multivariate normal distribution.
#' The algorithm allows the user to decide whether she/he wants to center and scale their data.
#' The user is also allowed to change the default values of the parameters of each prior.
#'
#'
#' The model is
#' \deqn{X_i = \theta w_i + \sigma \epsilon_i}
#' where \eqn{w_i \sim N(0, I_r), \epsilon \sim N(0,I_p)}.
#'
#' The spike and slab prior is given by
#'
#' \deqn{\pi(\theta, \boldsymbol \gamma|\lambda_1, r) \propto \prod_{j=1}^p \left(\gamma_j \int_{A \in V_{r,r}} g(\theta_j|\lambda_1, A, r) \pi(A) d A+ (1-\gamma_j) \delta_0(\theta_j)\right)}
#' \deqn{g(\theta_j|\lambda_1, A, r) = C(\lambda_1)^r \exp(-\lambda_1 \|\beta_j\|_q^m)}
#' \deqn{\gamma_j| \kappa \sim Bernoulli(\kappa)}
#' \deqn{\kappa \sim Beta(\alpha_1, \alpha_2)}
#' \deqn{\sigma^2 \sim InvGamma(\sigma_a, \sigma_b)}
#' where \eqn{V_{r,r} = \{A \in R^{r \times r}: A'A = I_r\}} and \eqn{\delta_0} is the Dirac measure at zero.
#' The density \eqn{g} can be chosen to be the product of independent Laplace distribution (i.e., \eqn{q = 1, m =1}) or the multivariate normal distribution (i.e., \eqn{q = 2, m = 2}).
#’
#'
#' @references Ning, B. (2021). Spike and slab Bayesian sparse principal component analysis. arXiv:2102.00305.
#'
#'
#'@param dat Data an \eqn{n*p} matrix.
#'@param r Rank.
#'@param lambda Tuning parameter for the density \eqn{g}.
#'@param slab.prior The density \eqn{g}, the default is "MVN", the multivariate normal distribution. Another choice is "Laplace".
#'@param max.iter The maximum number of iterations for running the algorithm.
#'@param eps The convergence threshold; the default is \eqn{10^{-4}}.
#'@param jointly.row.sparse The default is true, which means that the jointly row sparsity assumption is used; one could not use this assumptio by changing it to false.
#'@param sig2.true The default is false, \eqn{\sigma^2} will be estimated; if sig2 is known and its value is given, then \eqn{\sigma^2} will not be estimated.
#'@param threshold The threshold to determine whether \eqn{\gamma_j} is 0 or 1; the default value is 0.5.
#'@param center.scale The default if false. If true, then the input date will be centered and scaled.
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
#' \item{loadings}{The loadings matrix.}
#' \item{uncertainty}{The covariance of each non-zero rows in the loadings matrix.}
#' \item{scores}{Score functions for the \eqn{r} principal components.}
#' \item{sig2}{Variance of the noise.}
#' \item{obj.fn}{A vector contains the value of the objective function of each iteration. It can be used to check whether the algorithm converges}
#'
#'
#' @examples
#' #In this example, the first 20 rows in the loadings matrix are nonzero, the rank is 2
#' set.seed(2021)
#' library(MASS)
#' library(pracma)
#' n <- 200
#' p <- 1000
#' s <- 20
#' r <- 2
#' sig2 <- 0.1
#' # generate eigenvectors
#' U.s <- randortho(s, type = c("orthonormal"))
#' if (r == 1) {
#'   U <- rep(0, p)
#'   U[1:s] <- as.vector(U.s[, 1:r])
#' } else {
#'   U <- matrix(0, p, r)
#'   U[1:s, ] <- U.s[, 1:r]
#' }
#' s.star <- rep(0, p)
#' s.star[1:s] <- 1
#' eigenvalue <- seq(20, 10, length.out = r)
#' # generate Sigma
#' if (r == 1) {
#'   theta.true <- U * sqrt(eigenvalue)
#'   Sigma <- tcrossprod(theta.true) + sig2*diag(p)
#' } else {
#'   theta.true <- U %*% sqrt(diag(eigenvalue))
#'   Sigma <- tcrossprod(theta.true) + sig2 * diag(p)
#' }
#' # generate n*p dataset
#' X <- t(mvrnorm(n, mu = rep(0, p), Sigma = Sigma))
#' result <- VBsparsePCA(dat = t(X), r = 2, jointly.row.sparse = TRUE, center.scale = FALSE)
#' loadings <- result$loadings
#' scores <- result$scores
#' @export

VBsparsePCA <- function(
  dat, r, lambda = 1, slab.prior = "MVN", max.iter = 100, eps = 1e-3,
  jointly.row.sparse = TRUE, center.scale = FALSE,
  sig2.true = NA, threshold = 0.5, theta.int = NA,
  theta.var.int = NA, kappa.para1 = NA, kappa.para2 = NA,
  sigma.a = NA, sigma.b = NA
) {

  # if center and scale the data is needed
  if (center.scale == TRUE) {
    x <- t(scale(t(dat), center = TRUE, scale = TRUE))
  } else {
    x <- dat
  }

  if (r == 1) {
    if (slab.prior == "MVN") {

      result <- spca.cavi.mvn(x, r, lambda, max.iter, eps, jointly.row.sparse,
                              sig2.true, threshold, theta.int, theta.var.int,
                              kappa.para1, kappa.para2, sigma.a, sigma.b)

    } else if (slab.prior == "Laplace") {

      result <- spca.cavi.Laplace(x, r, lambda, max.iter, eps,
                                  sig2.true, threshold, theta.int, theta.var.int,
                                  kappa.para1, kappa.para2, sigma.a, sigma.b)
    } else {
      warning("The input slab prior is invalid")
    }


    # collect result
    selection <- result$selection
    loadings <- as.matrix(result$theta.mean)
    scores <- as.matrix(dat %*% svd(loadings)$u)
    uncertainty <- result$theta.var
    sig2 <- result$sig2
    obj.fn <- result$obj.fn

    # rename the loadings and scores
    # rename the loadings and scores
    if (is.null(colnames(dat)) == FALSE) {
      rownames(loadings) <- colnames(dat)
    }
    coln <- paste("PC", rep(1:r), sep = "")
    colnames(loadings) <- coln

    if (is.null(rownames(dat)) == FALSE) {
      rownames(scores) <- rownames(dat)
    }
    colnames(scores) <- coln

    return(list(iter = result$iter, selection = selection,
                loadings = loadings, uncertainty  = uncertainty,
                scores = scores, sig2 = sig2, obj.fn = obj.fn))

  } else if (r > 1) {

    if (slab.prior == "MVN") {

      result <- spca.cavi.mvn(x, r, lambda, max.iter, eps, jointly.row.sparse,
                              sig2.true, threshold, theta.int, theta.var.int,
                              kappa.para1, kappa.para2, sigma.a, sigma.b)

    } else if (slab.prior == "Laplace") {
      warning("The solution using the Laplace prior when r > 2 is unstable,
              the program will automatically switch to use the multivariate normal slab.")

      result <- spca.cavi.mvn(x, r, lambda, max.iter, eps, jointly.row.sparse,
                              sig2.true, threshold, theta.int, theta.var.int,
                              kappa.para1, kappa.para2, sigma.a, sigma.b)
    } else {
      warning("The input slab prior is invalid")
    }
  }

  # collect result
  selection <- result$selection
  loadings <- result$theta.mean
  scores <- dat %*% svd(loadings)$u
  uncertainty <- result$theta.var
  sig2 <- result$sig2
  obj.fn <- result$obj.fn

  # rename the loadings and scores
  if (is.null(colnames(dat)) == FALSE) {
    rownames(loadings) <- colnames(dat)
  }
  coln <- paste("PC", rep(1:r), sep = "")
  colnames(loadings) <- coln

  if (is.null(rownames(dat)) == FALSE) {
    rownames(scores) <- rownames(dat)
  }
  colnames(scores) <- coln

  return(list(iter = result$iter, selection = selection,
              loadings = loadings, uncertainty  = uncertainty,
              scores = scores, sig2 = sig2, obj.fn = obj.fn))
}
