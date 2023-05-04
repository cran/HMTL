\name{HMTL-package}
\alias{HMTL}
\docType{package}
\title{
Heterogeneous Multi-task Feature Learning
}
\description{
\code{HMTL} package implements the block-wise sparse estimation by grouping the coefficients of related predictors across multiple tasks. The tasks can be either regression, Huber regression, adaptive Huber regression, and logistic regression, which provide a wide variety of data types for the integration. The robust methods, such as the Huber regression and adaptive Huber regression, can deal with outlier contamination based on Sun, Q., Zhou, W.-X. and Fan, J. (2020), and Wang, L., Zheng, C., Zhou, W. and Zhou, W.-X. (2021). The model selection applies the modified form of Bayesian information criterion to measure the mdoel performance, which has similar formulation as previous work developed by Gao, X.,  and Carroll, R. J., (2017).
}
\details{
In the context of multi-task learning, there are \eqn{K} different data sets obtained from \eqn{K} related sources. The data sets can be modeled by different types of learning tasks based on the data distributions.  Let the candidate features be denoted as \eqn{\{M_1,M_2,...,M_j,...,M_p \}}. When the integrated data sets have different measurements, we assume the predictors to share some similarities. For example, the \eqn{j}th predictors collected as \eqn{M_j = (X_{1j}, X_{2j}, \cdots, X_{Kj})} in the table below represent the same type of feature in all related studies. In some cases, the tasks can share same set of predictor, then \eqn{X_{1j} = X_{2j} = \cdots = X_{Kj}}.
\tabular{llcccccc}{
Tasks \tab  Formula \tab \eqn{M_1} \tab \eqn{M_2} \tab \eqn{\dots} \tab \eqn{M_j} \tab \eqn{\dots} \tab \eqn{M_p} \cr
1 \tab \eqn{y_1: g_1(\mu_1) \sim} \tab \eqn{x_{11}\theta_{11}+} \tab \eqn{x_{12}\theta_{12}+} \tab \eqn{\dots}  \tab \eqn{x_{1j}\theta_{1j}+} \tab \eqn{\dots} \tab \eqn{x_{1p}\theta_{1p}} \cr
2 \tab \eqn{y_2: g_2(\mu_2) \sim} \tab \eqn{x_{21}\theta_{21}+} \tab \eqn{x_{22}\theta_{22}+} \tab \eqn{\dots} \tab \eqn{x_{2j}\theta_{2j}+} \tab \eqn{\dots}  \tab \eqn{x_{2p}\theta_{2p}}  \cr
\tab  ... \tab  \cr
K \tab \eqn{y_K: g_K(\mu_K) \sim} \tab \eqn{x_{K1}\theta_{K1}+} \tab \eqn{x_{K2}\theta_{K2}+} \tab \eqn{\dots} \tab \eqn{x_{Kj}\theta_{Kj}+} \tab \eqn{\dots} \tab \eqn{x_{Kp}\theta_{Kp}}  \cr
}
The coefficients can be grouped as the vector \eqn{\theta^{(j)}} for the feature \eqn{M_j}.
\tabular{llll}{
Platforms \tab \eqn{\bold{M_j}} \tab \tab \eqn{\bold{\theta^{(j)}}} \cr
1 \tab \eqn{x_{1j}} \tab \tab \eqn{\theta_{1j}} \cr
2 \tab \eqn{x_{2j}} \tab \tab \eqn{\theta_{2j}} \cr
 \tab ... \tab \tab ... \cr
k \tab \eqn{x_{Kj}} \tab \tab \eqn{\theta_{Kj}}
}

%Consider the following examples.

%\emph{Example 1.} Suppose \eqn{k} different types of experiments are conducted to study the genetic mechanism of a disease. The predictors in this research are different facets of individual genes, such as mRNA expression, protein expression, RNAseq expression and so on. The goal is to select the genes which affect the disease, while the genes are assessed in a number of ways through different measurement processes across \eqn{k} experimental platforms.

%\emph{Example 2.} The predictive models for three different financial indices are simultaneously built from a panel of stock index predictors. In this case, the predictor values across different models are the same, but the regression coefficients are different.

%In the conventional approach, the model for each of the \eqn{k} platforms is analyzed separately.
The heterogeneous multi-task feature learning \code{HMTL} can select significant features through the overall objective function:
\deqn{Q(\theta)=  \mathcal{L}(\theta) + \mathcal{R}(\theta).}
The loss function is defined as \eqn{\mathcal{L}(\theta) = \sum_{k=1}^K w_k \ell_k(\theta_k)}, which can be the composite quasi-likelihood or the composite form of (adaptive) Huber loss with additional robustification parameter \eqn{\tau_k}. The penalty function is the mixed \eqn{\ell_{2,1}} regularization, such that \eqn{\mathcal{R}(\theta) =  \lambda \sum_{j=1}^p (\sum_{k=1}^K \theta_{kj}^2)^{1/2}}.


This package also contains functions to provide the Bayesian information criterion:

\deqn{  BIC(s) = 2\mathcal{L}_s(\hat{\theta}) + d_s^{*} \gamma_n }

with \eqn{\mathcal{L}_s(\hat{\theta})} denoting the composite quasi-likelihood or adaptive Huber loss, \eqn{d_s^{*}} measuring the model complexity and \eqn{\gamma_n} being the penalty on the model complexity.

In this package, the function \code{\link{MTL_reg}} deals with regression tasks, which can be outlier contaminated. The function \code{\link{MTL_class}} is applied to model multiple classification tasks, and the function \code{\link{MTL_hetero}} can integrate different types of tasks together.
}
\author{
Yuan Zhong, Wei Xu, and Xin Gao

Maintainer: Yuan Zhong <aqua.zhong@gmail.com>
}
\references{
Zhong, Y., Xu, W., and Gao X., (2023) Heterogeneous multi-task feature learning with mixed \eqn{\ell_{2,1}} regularization. Submitted

Zhong, Y., Xu, W., and Gao X., (2023) Robust Multi-task Feature Learning. Submitted

Gao, X.,  and Carroll, R. J., (2017) Data integration with high dimensionality. Biometrika, 104, 2, pp. 251-272

Huber, P. J. (1964). Robust estimation of a location parameter. Ann. Math. Statist., 35, 73–101.

Sun, Q., Zhou, W.-X. and Fan, J. (2020). Adaptive Huber regression. J. Amer. Statist. Assoc., 115, 254-265.

Wang, L., Zheng, C., Zhou, W. and Zhou, W.-X. (2021). A new principle for tuning-free Huber regression. Stat. Sinica, 31, 2153-2177.
}




