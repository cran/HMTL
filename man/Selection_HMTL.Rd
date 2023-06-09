% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HMTL.R
\name{Selection_HMTL}
\alias{Selection_HMTL}
\title{Model Selection for Multi-task Feature Learning based on Bayesian Information Criterion (BIC)}
\usage{
Selection_HMTL(
  y,
  x,
  lambda,
  Kn,
  p,
  n,
  beta = 0.1,
  tau = 1.45,
  Cont_Model = "adaptive Huber",
  type = "Heterogeneity",
  import_w = 1,
  tol = 0.05,
  max_iter = 100,
  Complete = "True",
  diagnostics = FALSE,
  gamma = 1,
  alpha = 1
)
}
\arguments{
\item{y}{List. A list of responses vectors for all tasks. The order of the list put the continuous responses before the binary responses.}

\item{x}{List. Listing matrices of the predictors for all tasks align with the same order as in y.}

\item{lambda}{Numeric. The penalty parameter used for block-wise regularization (\eqn{\ell_{2,1}} norm).}

\item{Kn}{Vector of two elements. First element is the number of tasks with continuous responses, and the second element is the number of tasks with binary responses.}

\item{p}{Numeric. The number of features.}

\item{n}{Numeric or vector. If only one numeric value is provided, equal sample size will be assumed for each task. If a vector is provided, then the elements are the sample sizes for all tasks.}

\item{beta}{(\strong{optional}). Numeric or matrix. An initial value or matrix of values \eqn{p} by \eqn{K} for the estimation. The default value is 0.1.}

\item{tau}{Numeric or vector. The robustification parameter used for methods "Huber regression" or "Adaptive Huber". The default value is 1.45.}

\item{Cont_Model}{Character("regression", "Huber regression", or "adaptive Huber"). The models used for tasks with continuous responses.}

\item{type}{Character("Heterogeneity", "Continuous" or "Binary"). \code{type = "Heterogeneity"} used for \code{ MTL_hetero}; \code{type = "Continuous"} used for \code{ MTL_reg}; and \code{type = "Binary"} used for \code{ MTL_class}.}

\item{import_w}{Numeric or vector. The weights assigned to different tasks. An equal weight is set as the default.}

\item{tol}{(\strong{optional}). Numeric. The tolerance level of optimation.}

\item{max_iter}{(\strong{optional}). Numeric. The maximum number of iteration steps.}

\item{Complete}{Logic input. If the predictors in each task are all measured, set `Complete == TRUE`; If some predictors in some but not all task are all measured, set`Complete == FALSE`, and the missing values are imputed by column mean. The adjustment weights will be assigned based on the completeness of the predictors.}

\item{diagnostics}{Logic input. If `diagnostics == TRUE`, the function provides Bayesian information criterion, and the selected model performance is evalued by the MSE and MAE for tasks with continuous response and the AUC and deviance for tasks with binary responses.}

\item{gamma}{(\strong{optional}). Numeric. Step size for each inner iteration. The default is equal to 1.}

\item{alpha}{(\strong{optional}). Numeric. A tuning parameter for BIC penalty. The default is equal to 1.}
}
\value{
A table of Bayesian Information Criterion (BIC)
\describe{
\item{\code{lambda}}{A list of penalty parameters.}
\item{\code{Compo_likelihood}}{Sum of empirical loss functions estimated based on the selected parameters .}
\item{\code{Degree of freedom}}{Penalty component based on the selected parameters.}
\item{\code{Info criterion}}{Bayesian Information Criterion (BIC): robust BIC or pseudo BIC.}
}
}
\description{
\code{Selection_HMTL} can be used to search the optimal candidate model based on Bayesian Information Criterion (BIC).
}
\details{
The Bayesian information criterion is given by
\deqn{  BIC(s) = 2\mathcal{L}_s(\hat{\theta}) + d_s^{*} \gamma_n, }
where \eqn{\hat{\theta}} is the estimated coefficients and \eqn{s} is denoted the selected support set of \eqn{\hat{\theta}}.
In addition, \eqn{\mathcal{L}_s(\hat{\theta})} denoted the estimated composite quasi-likelihood or adaptive Huber loss evaluated as the \eqn{\hat{\theta}}.
The degree of freedom \eqn{d_s^{*}} can measure the model complexity, which is estimated by \eqn{tr(H_s^{-1}(\hat{\theta}) J_s(\hat{\theta}) )}. The sensitivity matrix and specificity matrix can be given by \eqn{H(\theta) = E(\nabla^2 \mathcal{L}( {\theta}))} and \eqn{J(\theta) = -Cov(\nabla \mathcal{L}( {\theta}))}.
The penalty term \eqn{\gamma_n} can be defined by users.
}
\examples{
lambda <- c(1.2, 2.0, 2.3, 2.8)
cv.mtl <- Selection_HMTL(mockdata2,mockdata1, lambda =lambda, Kn = c(2,2), p=500,
               n = c(500,250,250,250),gamma = 2, Complete = FALSE, diagnostics = TRUE, alpha = 2)
plot_HMTL(cv.mtl)
cv.mtl$Selection_Results
}
\references{
Y. Zhong, W. Xu, and X. Gao (2023) Robust Multi-task Feature Learning. Submitted

Gao, X and Carroll, R. J. (2017) Data integration with high dimensionality. Biometrika, 104, 2, pp. 251-272
}
