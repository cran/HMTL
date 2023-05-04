\name{mockdata}
\alias{mockdata1}
\alias{mockdata2}
\docType{data}
\title{
Mock Gene Data
}
\description{
This data set is a mock version of two related research outcomes.
}
\usage{data("mockdata1")
data("mockdata2")}
\format{
The \code{mockdata1} contains all predictor variables for four data sets, and the \code{mockdata2} is the list of four columns of response variables. The sample sizes are \eqn{n_1 = 500, n_2 = 250, n_3 = 250}, and \eqn{n_4 = 250}.

The response varibles are heterogeneous, such that the first two columns are continuous data and last two are binary data. The predictors of each data matrix in \code{mockdata1} have 500 columns of variables, and the columns in different data matrices are matched. For example, the first column in the first data matrix represents the same type of feature as the first columns in other three data matrices. Among all canadidate predictors, the response variables are set to be associated with the first 22 columns, and the remaining columns are not important predictors.

}
\details{
This data set is used as an example to implement functions in the package.
}


