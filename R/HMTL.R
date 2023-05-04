#' @title Heterogeneous Multi-task Feature Learning
#' @description \code{MTL_hetero} conducts multi-tasks feature learning to different types of learning tasks, including linear regression, Huber regression, adaptive Huber, and logistic regression. The penalty function applies a mixed \eqn{\ell_{2,1}} norm to combine regression coefficients of predictor shared across all tasks.
#' @param y List. A list of responses vectors for all tasks. The order of the list put the continuous responses before the binary responses.
#' @param x List. Listing matrices of the predictors for all tasks align with the same order as in y.
#' @param lambda Numeric. The penalty parameter used for block-wise regularization (\eqn{\ell_{2,1}} norm).
#' @param Kn Vector of two elements. First element is the number of tasks with continuous responses, and the second element is the number of tasks with binary responses.
#' @param p Numeric. The number of features.
#' @param n Numeric or vector. If only one numeric value is provided, equal sample size will be assumed for each task. If a vector is provided, then the elements are the sample sizes for all tasks.
#' @param beta (\strong{optional}). Numeric or matrix. An initial value or matrix of values \eqn{p} by \eqn{K} for the estimation. The default value is 0.1.
#' @param tau Numeric or vector. The robustification parameter used for methods "Huber regression" or "Adaptive Huber". The default value is 1.45.
#' @param Cont_Model Character("regression", "Huber regression", or "adaptive Huber"). The models used for tasks with continuous responses.
#' @param import_w Numeric or vector. The weights assigned to different tasks. An equal weight is set as the default.
#' @param tol (\strong{optional}). Numeric. The tolerance level of optimation.
#' @param max_iter (\strong{optional}). Numeric. The maximum number of iteration steps.
#' @param Complete Logic input. If the predictors in each task are all measured, set `Complete == TRUE`; If some predictors in some but not all task are all measured, set`Complete == FALSE`, and the missing values are imputed by column mean. The adjustment weights will be assigned based on the completeness of the predictors.
#' @param diagnostics Logic input. If `diagnostics == TRUE`, the function provides Bayesian information criterion, and the selected model performance is evalued by the MSE and MAE for tasks with continuous response and the AUC and deviance for tasks with binary responses.
#' @param gamma (\strong{optional}). Numeric. Step size for each inner iteration. The default is equal to 1.
#' @param alpha (\strong{optional}). Numeric. A tuning parameter for BIC penalty. The default is equal to 1.
#' @return A list including the following terms will be returned:
#' \describe{
#' \item{\code{beta}}{A \eqn{p} by \eqn{K} matrix of estimated sparse parameters.}
#' \item{\code{Task type}}{The models used in each task.}
#' \item{\code{Task weights}}{The weights assigned to each task.}
#' \item{\code{Selected_List}}{The index of non-zero parameters.}
#' }
#' @return If `diagnostics = TRUE`, the following terms will be returned:
#' \describe{
#' \item{\code{Bayesian_Information}}{Table of the information criterion: Composite likelihood, Degree of freedom, and (peudo or robust) Bayesian informtion criterion.}
#' \item{\code{Reg_Error}}{Table of the model performance for (Huber) regressions: the mean square error (MSE), and the mean absolute error (MAE).}
#' \item{\code{Class_Perform}}{Table of the model performance for classification tasks: the area under ROC curve (AUC), and the deviance (DEV) estimated by `glm`.}
#' \item{\code{Residuals}}{The residuals for all tasks.}
#' }
#' @references Zhong, Y., Xu, W., and Gao X., (2023) Heterogeneous multi-task feature learning with mixed \eqn{\ell_{2,1}} regularization. Submitted
#' @note When penalty parameter is too small, the estimated coefficients may have \eqn{p \ge n}. The algorithm can provide the estimated values, but the \code{diagnostics}'s results will not be given.
#' @importFrom stats na.omit binomial cov glm var
#' @importFrom pROC roc auc
#' @importFrom Matrix bdiag
#' @examples
#' model <- MTL_hetero(mockdata2,mockdata1, lambda =2.5, Kn = c(2,2), p=500,n = c(500,250,250,250),
#'                     gamma = 2, Complete = FALSE, diagnostics = TRUE, alpha = 2)
#' # Selected non-zero coefficients
#' model$beta[model$Selected_List,]
#' # Estimated Pseudo-BIC
#' model$Bayesian_Information
#' # Regression error
#' model$Reg_Error
#' # Classification accuracy
#' model$Class_Perform
#' @export
MTL_hetero <- function(y, x, lambda, Kn, p, n, beta = 0.1, tau = 1.45,
                       Cont_Model = "adaptive Huber" , import_w = 1,
                       tol = 0.05, max_iter = 100, Complete = "True",
                       diagnostics = FALSE, gamma = 1, alpha = 1){

  K1 <- Kn[1]
  K2 <- Kn[2]

  K <- K1+K2

  methods <- Cont_Model

  if (length(methods ) == 1){
    methods <- rep(methods, K1)
  }
  link = rep("logit",K2)
  if( length(tau ) == 1 ) {
    tau <-  rep(tau, K1)
  }

  if (length(n) == 1) {
    n <- rep(n, K)
  }
  N <- max(n)

  sample <- dimension <-  numeric(K)
  for( k in 1:K){
    sample[k] <- nrow(x[[k]])
    dimension[k] <- ncol(x[[k]])
    if(sample[k]!=n[k])
      warning(paste("The samples of the task", k, "not match with input values" ))
    if (dimension[k]!=p)
      warning(paste("The dimension of predictors in the task", k, "not match with input values" ))
    if(sample[k]!=length(y[[k]]))
      warning(paste("The number of observations in the predictors not aligned with the responses in task", k))
    if (length(y)!= K)
      warning("The number of tasks not match with input values")
    if (length(x)!= K)
      warning("The number of tasks not match with input values")
  }

  Output <- .Main_Algo(x, y, Kn = Kn, p = p, n = n, lambda = lambda, beta = beta, tau = tau,gamma =gamma, methods ,
                      import_w = import_w, tol = tol, max_iter = max_iter, Complete = Complete)

  beta <- Output$beta
  beta[is.na(beta)] <- 0


  info_fitting<-NULL
  if (diagnostics == TRUE){
    info_fitting <-  .info_fitting(x, y, beta, Kn =Kn, n=n, p=p, tau=tau, gamma= gamma, methods,  import_w = import_w,
                                  alpha = alpha,tol=tol,Complete = Complete)
  }


  if( Output[5] == max_iter)
    warning(paste("Algorithm did not converge after", max_iter, "iterations"))


  Results <- c(Output, info_fitting)
  return(Results)

}
#'
#' @title Robust Multi-Task Feature Learning
#' @description \code{MTL_reg} conducts multi-tasks feature learning to the learning tasks with continous response variables, such as the linear regression, Huber regression, adaptive Huber. The adaptive Huber method is based on Sun, Q., Zhou, W.-X. and Fan, J. (2020) and Wang, L., Zheng, C., Zhou, W. and Zhou, W.-X. (2021). The penalty function applies a mixed \eqn{\ell_{2,1}} norm to combine regression coefficients of predictor shared across all tasks.
#' The Huber regression and adaptive Huber regression need the robustification parameter \eqn{\tau_k} to strike a balance between the unbiasedness and robustness, and the adaptive method can determine this parameter by a tuning-free principle.
#' @param y List. A list of continuous responses vectors for all tasks.
#' @param x List. Listing matrices of the predictors for all tasks align with the same order as in y.
#' @param lambda Numeric. The penalty parameter used for block-wise regularization (\eqn{\ell_{2,1}} norm).
#' @param Kn Numeric. The number of tasks with continuous responses.
#' @param p Numeric. The number of features.
#' @param n Numeric or vector. If only one numeric value is provided, equal sample size will be assumed for each task. If a vector is provided, then the elements are the sample sizes for all tasks.
#' @param beta (\strong{optional}). Numeric or matrix. An initial value or matrix of values \eqn{p} by \eqn{K} for the estimation. The default value is 0.1.
#' @param tau Numeric or vector. The robustification parameter used for methods "Huber regression" or "Adaptive Huber". The default value is 1.45.
#' @param Cont_Model Character("regression", "Huber regression", or "adaptive Huber"). The models used for tasks with continuous responses.
#' @param import_w Numeric or vector. The weights assigned to different tasks. An equal weight is set as the default.
#' @param tol (\strong{optional}). Numeric. The tolerance level of optimation.
#' @param max_iter (\strong{optional}). Numeric. The maximum number of iteration steps.
#' @param Complete Logic input. If the predictors in each task are all measured, set `Complete == TRUE`; If some predictors in some but not all task are all measured, set`Complete == FALSE`, and the missing values are imputed by column mean. The adjustment weights will be assigned based on the completeness of the predictors.
#' @param diagnostics Logic input. If `diagnostics == TRUE`, the function provides Bayesian information criterion, and the selected model performance is evalued by the MSE and MAE for tasks with continuous response and the AUC and deviance for tasks with binary responses.
#' @param gamma (\strong{optional}). Numeric. Step size for each inner iteration. The default is equal to 1.
#' @param alpha (\strong{optional}). Numeric. A tuning parameter for BIC penalty. The default is equal to 1.
#' @return A list including the following terms will be returned:
#' \describe{
#' \item{\code{beta}}{A \eqn{p} by \eqn{K} matrix of estimated sparse parameters.}
#' \item{\code{Task type}}{The models used in each task.}
#' \item{\code{Task weights}}{The weights assigned to each task.}
#' \item{\code{Selected_List}}{The index of non-zero parameters.}
#' }
#' @return If `diagnostics = TRUE`, the following terms will be returned:
#' \describe{
#' \item{\code{Bayesian_Information}}{Table of the information criterion: Composite likelihood, Degree of freedom, and (peudo or robust) Bayesian informtion criterion.}
#' \item{\code{Reg_Error}}{Table of the model performance for (Huber) regressions: the mean square error (MSE), and the mean absolute error (MAE).}
#' \item{\code{Residuals}}{The residuals for all tasks.}
#' }
#' @references Huber, P. J. (1964). Robust estimation of a location parameter. Ann. Math. Statist., 35, 73–101.
#' @references Sun, Q., Zhou, W.-X. and Fan, J. (2020). Adaptive Huber regression. J. Amer. Statist. Assoc., 115, 254-265.
#' @references Wang, L., Zheng, C., Zhou, W. and Zhou, W.-X. (2021). A new principle for tuning-free Huber regression. Stat. Sinica, 31, 2153-2177.
#' @references Zhong, Y., Xu, W., and Gao X., (2023) Robust Multi-task Feature Learning. Submitted
#' @importFrom stats na.omit binomial cov glm var
#' @importFrom Matrix bdiag
#' @examples
#' x_reg <- list(mockdata1[[1]],mockdata1[[2]])
#' y_reg <- list(mockdata2[[1]],mockdata2[[2]])
#' model <- MTL_reg(y_reg,x_reg, lambda = 2.5  , Kn = 2, p=500,
#'                 n = c(500,250 ),gamma = 2, Complete = FALSE, diagnostics = TRUE, alpha = 2)
#'
#' # Selected non-zero coefficients
#' model$beta[model$Selected_List,]
#' # Estimated Pseudo-BIC
#' model$Bayesian_Information
#' # Regression error
#' model$Reg_Error
#' @export
MTL_reg <- function(y, x, lambda, Kn, p, n, beta = 0.1, tau = 1.45,
                    Cont_Model = "adaptive Huber" , import_w = 1,
                       tol = 0.05, max_iter = 100, Complete = "True", diagnostics = FALSE,gamma = 1, alpha = 1){
  Results  <- MTL_hetero(y=y, x=x, lambda=lambda, Kn = c(Kn,0), p = p, n = n, beta = beta, tau = tau, import_w = import_w,Cont_Model = Cont_Model,
                         tol = tol, max_iter = max_iter, Complete = Complete, diagnostics = diagnostics, gamma = gamma, alpha = alpha)
  return(Results)
}
#'
#' @title Multiple Classification Task Feature Learning
#' @description \code{MTL_class} conducts multi-tasks feature learning to the learning tasks with binary response variables, namely logistic regression. The penalty function applies a mixed \eqn{\ell_{2,1}} norm to combine regression coefficients of predictor shared across all tasks.
#' @param y List. A list of binary responses vectors for all tasks.
#' @param x List. Listing matrices of the predictors for all tasks align with the same order as in y.
#' @param lambda Numeric. The penalty parameter used for block-wise regularization (\eqn{\ell_{2,1}} norm).
#' @param Kn Numeric. The number of tasks with binary responses.
#' @param p Numeric. The number of features.
#' @param n Numeric or vector. If only one numeric value is provided, equal sample size will be assumed for each task. If a vector is provided, then the elements are the sample sizes for all tasks.
#' @param beta (\strong{optional}). Numeric or matrix. An initial value or matrix of values \eqn{p} by \eqn{K} for the estimation. The default value is 0.1.
#' @param import_w Numeric or vector. The weights assigned to different tasks. An equal weight is set as the default.
#' @param tol (\strong{optional}). Numeric. The tolerance level of optimation.
#' @param max_iter (\strong{optional}). Numeric. The maximum number of iteration steps.
#' @param Complete Logic input. If the predictors in each task are all measured, set `Complete == TRUE`; If some predictors in some but not all task are all measured, set`Complete == FALSE`, and the missing values are imputed by column mean. The adjustment weights will be assigned based on the completeness of the predictors.
#' @param diagnostics Logic input. If `diagnostics == TRUE`, the function provides Bayesian information criterion, and the selected model performance is evalued by the MSE and MAE for tasks with continuous response and the AUC and deviance for tasks with binary responses.
#' @param gamma (\strong{optional}). Numeric. Step size for each inner iteration. The default is equal to 1.
#' @param alpha (\strong{optional}). Numeric. A tuning parameter for BIC penalty. The default is equal to 1.
#' @return A list including the following terms will be returned:
#' \describe{
#' \item{\code{beta}}{A \eqn{p} by \eqn{K} matrix of estimated sparse parameters.}
#' \item{\code{Task type}}{The models used in each task.}
#' \item{\code{Task weights}}{The weights assigned to each task.}
#' \item{\code{Selected_List}}{The index of non-zero parameters.}
#' }
#' @return If `diagnostics = TRUE`, the following terms will be returned:
#' \describe{
#' \item{\code{Bayesian_Information}}{Table of the information criterion: Composite likelihood, Degree of freedom, and (peudo or robust) Bayesian informtion criterion.}
#' \item{\code{Class_Perform}}{Table of the model performance for classification tasks: the area under ROC curve (AUC), and the deviance (DEV) estimated by `glm`.}
#' \item{\code{Residuals}}{The residuals for all tasks.}
#' }
#' @references Zhong, Y., Xu, W., and Gao X., (2023) Heterogeneous multi-task feature learning with mixed \eqn{\ell_{2,1}} regularization. Submitted
#' @importFrom stats na.omit binomial cov glm var
#' @importFrom Matrix bdiag
#' @importFrom pROC auc roc
#' @examples
#'
#' x_class <- list(mockdata1[[3]],mockdata1[[4]])
#' y_class <- list(mockdata2[[3]],mockdata2[[4]])
#' model <- MTL_class(y_class,x_class, lambda = 2/11 , Kn = 2, p=500,
#'                   n = 250 ,gamma = 1, Complete = FALSE, diagnostics = TRUE, alpha = 2)
#' # Selected non-zero coefficients
#' model$beta[model$Selected_List,]
#' # Estimated Pseudo-BIC
#' model$Bayesian_Information
#' # Classification accuracy
#' model$Class_Perform
#' @export
MTL_class <- function(y, x, lambda, Kn, p, n, beta = 0.1,  import_w = 1,
                       tol = 0.05, max_iter = 100, Complete = "True", diagnostics = FALSE, gamma = 1, alpha = 1){
  Results  <- MTL_hetero(y=y, x=x, lambda=lambda, Kn = c(0,Kn), p=p, n=n, beta = beta, tau = NULL, import_w = import_w,
                         tol = tol, max_iter = max_iter, Complete = Complete, diagnostics = diagnostics, gamma = gamma, alpha = alpha)
  return(Results)

}
#'
#' @title Model Selection for Multi-task Feature Learning based on Bayesian Information Criterion (BIC)
#' @description \code{Selection_HMTL} can be used to search the optimal candidate model based on Bayesian Information Criterion (BIC).
#' @param y List. A list of responses vectors for all tasks. The order of the list put the continuous responses before the binary responses.
#' @param x List. Listing matrices of the predictors for all tasks align with the same order as in y.
#' @param lambda Numeric. The penalty parameter used for block-wise regularization (\eqn{\ell_{2,1}} norm).
#' @param Kn Vector of two elements. First element is the number of tasks with continuous responses, and the second element is the number of tasks with binary responses.
#' @param p Numeric. The number of features.
#' @param n Numeric or vector. If only one numeric value is provided, equal sample size will be assumed for each task. If a vector is provided, then the elements are the sample sizes for all tasks.
#' @param beta (\strong{optional}). Numeric or matrix. An initial value or matrix of values \eqn{p} by \eqn{K} for the estimation. The default value is 0.1.
#' @param tau Numeric or vector. The robustification parameter used for methods "Huber regression" or "Adaptive Huber". The default value is 1.45.
#' @param Cont_Model Character("regression", "Huber regression", or "adaptive Huber"). The models used for tasks with continuous responses.
#' @param type Character("Heterogeneity", "Continuous" or "Binary"). \code{type = "Heterogeneity"} used for \code{ MTL_hetero}; \code{type = "Continuous"} used for \code{ MTL_reg}; and \code{type = "Binary"} used for \code{ MTL_class}.
#' @param import_w Numeric or vector. The weights assigned to different tasks. An equal weight is set as the default.
#' @param tol (\strong{optional}). Numeric. The tolerance level of optimation.
#' @param max_iter (\strong{optional}). Numeric. The maximum number of iteration steps.
#' @param Complete Logic input. If the predictors in each task are all measured, set `Complete == TRUE`; If some predictors in some but not all task are all measured, set`Complete == FALSE`, and the missing values are imputed by column mean. The adjustment weights will be assigned based on the completeness of the predictors.
#' @param diagnostics Logic input. If `diagnostics == TRUE`, the function provides Bayesian information criterion, and the selected model performance is evalued by the MSE and MAE for tasks with continuous response and the AUC and deviance for tasks with binary responses.
#' @param gamma (\strong{optional}). Numeric. Step size for each inner iteration. The default is equal to 1.
#' @param alpha (\strong{optional}). Numeric. A tuning parameter for BIC penalty. The default is equal to 1.
#' @return A table of Bayesian Information Criterion (BIC)
#' \describe{
#' \item{\code{lambda}}{A list of penalty parameters.}
#' \item{\code{Compo_likelihood}}{Sum of empirical loss functions estimated based on the selected parameters .}
#' \item{\code{Degree of freedom}}{Penalty component based on the selected parameters.}
#' \item{\code{Info criterion}}{Bayesian Information Criterion (BIC): robust BIC or pseudo BIC.}
#' }
#' @details The Bayesian information criterion is given by
#' \deqn{  BIC(s) = 2\mathcal{L}_s(\hat{\theta}) + d_s^{*} \gamma_n, }
#' where \eqn{\hat{\theta}} is the estimated coefficients and \eqn{s} is denoted the selected support set of \eqn{\hat{\theta}}.
#' In addition, \eqn{\mathcal{L}_s(\hat{\theta})} denoted the estimated composite quasi-likelihood or adaptive Huber loss evaluated as the \eqn{\hat{\theta}}.
#' The degree of freedom \eqn{d_s^{*}} can measure the model complexity, which is estimated by \eqn{tr(H_s^{-1}(\hat{\theta}) J_s(\hat{\theta}) )}. The sensitivity matrix and specificity matrix can be given by \eqn{H(\theta) = E(\nabla^2 \mathcal{L}( {\theta}))} and \eqn{J(\theta) = -Cov(\nabla \mathcal{L}( {\theta}))}.
#' The penalty term \eqn{\gamma_n} can be defined by users.
#' @references Y. Zhong, W. Xu, and X. Gao (2023) Robust Multi-task Feature Learning. Submitted
#' @references Gao, X and Carroll, R. J. (2017) Data integration with high dimensionality. Biometrika, 104, 2, pp. 251-272
#' @importFrom  stats na.omit binomial cov glm var
#' @importFrom pROC roc auc
#' @importFrom Matrix bdiag
#' @examples
#' lambda <- c(1.2, 2.0, 2.3, 2.8)
#' cv.mtl <- Selection_HMTL(mockdata2,mockdata1, lambda =lambda, Kn = c(2,2), p=500,
#'                n = c(500,250,250,250),gamma = 2, Complete = FALSE, diagnostics = TRUE, alpha = 2)
#' plot_HMTL(cv.mtl)
#' cv.mtl$Selection_Results
#' @export
Selection_HMTL <-  function(y, x, lambda, Kn, p, n, beta = 0.1, tau = 1.45,Cont_Model = "adaptive Huber",
                            type = "Heterogeneity", import_w = 1,  tol = 0.05,
                            max_iter = 100, Complete = "True", diagnostics = FALSE,  gamma = 1,  alpha = 1){

  if (length(lambda) != 1){
    M <- length(lambda)
  }
  else {
    M <- 100
    lambda <- sqrt(log(p)/n*seq(1,10, length.out = M))
  }
  methods = Cont_Model
  output <- matrix(NA, ncol = 4, nrow = M)
  output[,1] <- lambda
  for(i in 1:M){

    if("Heterogeneity" %in% type ){

      fit <-  MTL_hetero(y, x, lambda[i], Kn=Kn, p=p, n=n, beta = beta, tau = tau, Cont_Model = methods ,
                        import_w = import_w, tol = tol, max_iter = max_iter,
                         Complete = Complete, diagnostics = TRUE, gamma = gamma, alpha = alpha)

      if ( is.na(fit$Bayesian_Information[3]) == FALSE){
        output[i,-1] <- fit$Bayesian_Information
      }

    }

    if("Continuous" %in% type ){

      fit <-  MTL_reg(y, x, lambda[i], Kn=Kn, p=p, n=n,  beta = beta, tau = tau, Cont_Model = methods ,
                      import_w = import_w,  tol = tol, max_iter = max_iter, Complete = Complete,
                      diagnostics = TRUE, gamma = gamma, alpha = alpha)

      if ( is.na(fit$Bayesian_Information[3]) == FALSE){
        output[i,-1] <- fit$Bayesian_Information
      }

    }

    if("Binary" %in% type ){

      fit <-  MTL_class(y, x, lambda[i], Kn=Kn, p=p, n=n,  beta = beta,
                      import_w = import_w,  tol = tol, max_iter = max_iter, Complete = Complete,
                      diagnostics = TRUE, gamma = gamma, alpha = alpha)

      if ( is.na(fit$Bayesian_Information[3]) == FALSE){
        output[i,-1] <- fit$Bayesian_Information
      }

    }
}

  lambda_min <-  output[which.min(output[,4]),]

  colnames(output) <- c("Lambda",  "Compo_likelihood", "Degree of freedom ",   "Info criterion")
  results <- list( "Selection_Results" = output,"Lambda_min" = lambda_min)
  class(results) <- "Model Selection"
  return(results)

}
#' @title Plot diagram of the information criterion vs. penalty parameters
#' @description Plot a diagram to illustrate the change of Bayesian information criterion with different penalty paameters for model selection
#' @param x Object of class "Model Selection" created by \code{Selection_HMTL}
#' @return X axis represents the value of penalty parameters, and y axis represents the estimated values of the composite likelihood, degree of freedom for model complexity, and the (robust) Bayesian information criterion.
#' @importFrom graphics plot legend points
#' @export
plot_HMTL <- function(x){
  if ( inherits(x, "Model Selection",TRUE) != 1)
    stop("argument must be of class Model Selection")

   x$Selection_Results <- na.omit(x$Selection_Results)
  if ( dim(x$Selection_Results )[1] ==1  )
    stop("number of penalty parameter must be more than one")

    plot(x$Selection_Results[,1], x$Selection_Results[,4],xlab="Lambda",ylab="BIC", type = "l",
       main = "Model Selection", ylim = c(min(x$Selection_Results[,-1]),max(x$Selection_Results[,-1])))
    points(x$Selection_Results[,1],x$Selection_Results[,2],col="blue",lty = 2,type = "l")
    points(x$Selection_Results[,1],x$Selection_Results[,3],col="red",lty = 2,type = "l")
    legend("bottomleft", c("robust-BIC","Comp Likelihood","Degree of Freedom"),
         cex = 0.5,lty =c(1,2,2),col = c("black","blue","red"),
         pch=c(NA,NA,NA))

  }

