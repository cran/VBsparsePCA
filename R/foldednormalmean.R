#' The function for obtaining the mean of a folded normal distribution
#'
#'
#' This function calculates the mean of the folded normal distribution given its location and scale parameters.
#'
#'
#' The mean of the folded normal distribution with location \eqn{\mu} and scale \eqn{\sigma^2} is
#' \deqn{\sigma \sqrt{2/\pi} \exp(-\mu^2/(2\sigma^2)) + \mu (1-2\Phi(-\mu/\sigma))}.
#'
#'
#'@param mean Location parameter of the folded normal distribution.
#'@param var Scale parameter of the folded normal distribution.
#'
#'
#'@return \item{foldednorm.mean}{The mean of the folded normal distribution of iterations to reach convergence.}
#'
#' @examples
#' #Calculates the mean of the folded normal distribution with mean 0 and var 1
#' mean <- foldednorm.mean(0, 1)
#' print(mean)
#' @export

foldednorm.mean <- function(mean, var) {
  sqrt(var*2/pi) * exp(-mean^2/(2*var)) + mean*(1-2*pnorm(-mean/sqrt(var)))
}
