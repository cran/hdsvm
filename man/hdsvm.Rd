% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hdsvm.R
\name{hdsvm}
\alias{hdsvm}
\title{Solve Penalized SVM}
\usage{
hdsvm(
  x,
  y,
  nlambda = 100,
  lambda.factor = ifelse(nobs < nvars, 0.01, 1e-04),
  lambda = NULL,
  lam2 = 0,
  hval = 1,
  pf = rep(1, nvars),
  pf2 = rep(1, nvars),
  exclude,
  dfmax = nvars + 1,
  pmax = min(dfmax * 1.2, nvars),
  standardize = TRUE,
  eps = 1e-08,
  maxit = 1e+06,
  sigma = 0.9,
  is_exact = FALSE
)
}
\arguments{
\item{x}{Matrix of predictors, with dimensions (\eqn{n} observations by \eqn{p} variables).}

\item{y}{Response variable vector of length \eqn{n}.}

\item{nlambda}{Number of \code{lambda} values to consider (default is 100).}

\item{lambda.factor}{The factor for getting the minimal value
in the \code{lambda} sequence, where
\code{min(lambda)} = \code{lambda.factor} * \code{max(lambda)}
and \code{max(lambda)} is the smallest value of \code{lambda}
for which all coefficients (except the intercept when it is present)
are penalized to zero. The default depends on the relationship
between \eqn{n} (the number of rows in the design matrix) and
\eqn{p} (the number of predictors). If \eqn{n < p}, it defaults to
\code{0.05}. If \eqn{n > p}, the default is \code{0.001},
closer to zero.  A very small value of \code{lambda.factor} will
lead to a saturated fit. The argument takes no effect if there is a
user-supplied \code{lambda} sequence.}

\item{lambda}{A user-supplied \code{lambda} sequence. Typically,
by leaving this option unspecified, users can have the program
compute its own \code{lambda} sequence based on \code{nlambda}
and \code{lambda.factor}. It is better to supply, if necessary,
a decreasing sequence of \code{lambda} values than a single
(small) value. The program will ensure that the user-supplied
\code{lambda} sequence is sorted in decreasing order before}

\item{lam2}{Regularization parameter \code{lambda2} for the
quadratic penalty of the coefficients. Unlike \code{lambda},
only one value of \code{lambda2} is used for each fitting process.}

\item{hval}{Smoothing parameter for the smoothed hinge loss, default is 1.}

\item{pf}{L1 penalty factor of length \eqn{p} used for the adaptive
LASSO or adaptive elastic net. Separate L1 penalty weights can be
applied to each coefficient to allow different L1 shrinkage.
Can be 0 for some variables (but not all), which imposes no
shrinkage, and results in that variable always being included
in the model. Default is 1 for all variables (and implicitly
infinity for variables in the \code{exclude} list).}

\item{pf2}{L2 penalty factor of length \eqn{p} used for adaptive
elastic net. Separate L2 penalty weights can be applied to
each coefficient to allow different L2 shrinkage.
Can be 0 for some variables, which imposes no shrinkage.
Default is 1 for all variables.}

\item{exclude}{Indices of variables to be excluded from the model.
Default is none. Equivalent to an infinite penalty factor.}

\item{dfmax}{The maximum number of variables allowed in the model.
Useful for very large \eqn{p} when a partial path is desired.
Default is \eqn{p+1}.}

\item{pmax}{Maximum count of non-zero coefficients across the solution path.}

\item{standardize}{Logical flag for variable standardization,
prior to fitting the model sequence. The coefficients are
always returned to the original scale. Default is \code{TRUE}.}

\item{eps}{Convergence criterion for stopping the algorithm.}

\item{maxit}{Maximum number of iterations permitted.}

\item{sigma}{Penalty parameter in the quadratic term of the augmented Lagrangian.}

\item{is_exact}{If \code{TRUE}, solutions are computed exactly; otherwise, approximations are used.}
}
\value{
An object with S3 class \code{hdsvm} consisting of
  \item{call}{the call that produced this object}
  \item{b0}{intercept sequence of length \code{length(lambda)}}
  \item{beta}{a \code{p*length(lambda)} matrix of coefficients,
              stored as a sparse matrix (\code{dgCMatrix} class,
              the standard class for sparse numeric matrices in
              the \code{Matrix} package.). To convert it into
              normal type matrix, use \code{as.matrix()}.}
  \item{lambda}{the actual sequence of \code{lambda} values used}
  \item{df}{the number of nonzero coefficients for each value
            of \code{lambda}.}
  \item{npasses}{the number of iterations for every lambda value}
  \item{jerr}{error flag, for warnings and errors, 0 if no error.}
}
\description{
Fits a penalized support vector machine (SVM) model using a range of \code{lambda} values,
allowing for detailed control over regularization parameters and model complexity.
}
\details{
The function utilizes the hinge loss function combined with elastic net penalization:
\deqn{1'[\max\{1 - y_i (\beta_0 + X_i^\top \beta), 0\}]/N + \lambda_1 \cdot |pf_1 \circ \beta|_1 +
0.5 \cdot \lambda_2 \cdot (\sqrt{pf_2} \circ \beta)^2,}
where \eqn{\circ} denotes the Hadamard product. 

For faster computation, if the algorithm is not converging or
running slow, consider increasing \code{eps}, increasing
\code{sigma}, decreasing \code{nlambda}, or increasing
\code{lambda.factor} before increasing \code{maxit}.
}
\examples{
set.seed(315)
n <- 100
p <- 400
x1 <- matrix(rnorm(n / 2 * p, -0.25, 0.1), n / 2)
x2 <- matrix(rnorm(n / 2 * p, 0.25, 0.1), n / 2)
x <- rbind(x1, x2)
beta <- 0.1 * rnorm(p)
prob <- plogis(c(x \%*\% beta))
y <- 2 * rbinom(n, 1, prob) - 1
lam2 <- 0.01
fit <- hdsvm(x, y, lam2=lam2)
}
\keyword{classification}
\keyword{svm}
