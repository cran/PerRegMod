\name{LSE_Reg_per}
\alias{LSE_Reg_per}
\title{Least squares estimator for periodic coefficients regression model}
\usage{
LSE_Reg_per(x,y,s)
}
\arguments{
  \item{x}{A list of independent variables with dimension \eqn{p}.}
  \item{y}{A response variable.}
  \item{s}{A period of the regression model.}
}

\description{
LSE_Reg_per() function gives the least squares estimation of parameters of a periodic coefficients regression model.
}
\value{
  \item{beta}{ Parameters to be estimated.}
  \item{X}{Matrix of predictors.}
  \item{Y}{The response vector.}
}

\examples{
set.seed(6)
n=400
s=4
x1=rnorm(n,0,1.5)
x2=rnorm(n,0,0.9)
x3=rnorm(n,0,2)
x4=rnorm(n,0,1.9)
y=rnorm(n,0,2.5)
x=list(x1,x2,x3,x4)
LSE_Reg_per(x,y,s)
}
