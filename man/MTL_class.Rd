% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HMTL.R
\name{MTL_class}
\alias{MTL_class}
\title{Multiple Classification Task Feature Learning}
\usage{
MTL_class(
  y,
  x,
  lambda,
  Kn,
  p,
  n,
  beta = 0.1,
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
\item{y}{List. A list of binary responses vectors for all tasks.}

\item{x}{List. Listing matrices of the predictors for all tasks align with the same order as in y.}

\item{lambda}{Numeric. The penalty parameter used for block-wise regularization (\eqn{\ell_{2,1}} norm).}

\item{Kn}{Numeric. The number of tasks with binary responses.}

\item{p}{Numeric. The number of features.}

\item{n}{Numeric or vector. If only one numeric value is provided, equal sample size will be assumed for each task. If a vector is provided, then the elements are the sample sizes for all tasks.}

\item{beta}{(\strong{optional}). Numeric or matrix. An initial value or matrix of values \eqn{p} by \eqn{K} for the estimation. The default value is 0.1.}

\item{import_w}{Numeric or vector. The weights assigned to different tasks. An equal weight is set as the default.}

\item{tol}{(\strong{optional}). Numeric. The tolerance level of optimation.}

\item{max_iter}{(\strong{optional}). Numeric. The maximum number of iteration steps.}

\item{Complete}{Logic input. If the predictors in each task are all measured, set `Complete == TRUE`; If some predictors in some but not all task are all measured, set`Complete == FALSE`, and the missing values are imputed by column mean. The adjustment weights will be assigned based on the completeness of the predictors.}

\item{diagnostics}{Logic input. If `diagnostics == TRUE`, the function provides Bayesian information criterion, and the selected model performance is evalued by the MSE and MAE for tasks with continuous response and the AUC and deviance for tasks with binary responses.}

\item{gamma}{(\strong{optional}). Numeric. Step size for each inner iteration. The default is equal to 1.}

\item{alpha}{(\strong{optional}). Numeric. A tuning parameter for BIC penalty. The default is equal to 1.}
}
\value{
A list including the following terms will be returned:
\describe{
\item{\code{beta}}{A \eqn{p} by \eqn{K} matrix of estimated sparse parameters.}
\item{\code{Task type}}{The models used in each task.}
\item{\code{Task weights}}{The weights assigned to each task.}
\item{\code{Selected_List}}{The index of non-zero parameters.}
}

If `diagnostics = TRUE`, the following terms will be returned:
\describe{
\item{\code{Bayesian_Information}}{Table of the information criterion: Composite likelihood, Degree of freedom, and (peudo or robust) Bayesian informtion criterion.}
\item{\code{Class_Perform}}{Table of the model performance for classification tasks: the area under ROC curve (AUC), and the deviance (DEV) estimated by `glm`.}
\item{\code{Residuals}}{The residuals for all tasks.}
}
}
\description{
\code{MTL_class} conducts multi-tasks feature learning to the learning tasks with binary response variables, namely logistic regression. The penalty function applies a mixed \eqn{\ell_{2,1}} norm to combine regression coefficients of predictor shared across all tasks.
}
\examples{

x_class <- list(mockdata1[[3]],mockdata1[[4]])
y_class <- list(mockdata2[[3]],mockdata2[[4]])
model <- MTL_class(y_class,x_class, lambda = 2/11 , Kn = 2, p=500,
                  n = 250 ,gamma = 1, Complete = FALSE, diagnostics = TRUE, alpha = 2)
# Selected non-zero coefficients
model$beta[model$Selected_List,]
# Estimated Pseudo-BIC
model$Bayesian_Information
# Classification accuracy
model$Class_Perform
}
\references{
Zhong, Y., Xu, W., and Gao X., (2023) Heterogeneous multi-task feature learning with mixed \eqn{\ell_{2,1}} regularization. Submitted
}
