\name{estimate_para_adaptive_method}
\alias{estimate_para_adaptive_method}
\title{Adaptive estimator for periodic coefficients regression model}
\usage{
estimate_para_adaptive_method(n,s,y,x)
}
\arguments{
   \item{n}{The length of vector  \eqn{y}.}
   \item{s}{A period of the regression model.}
  \item{y}{A response variable.}
  \item{x}{A list of independent variables with dimension \eqn{p}.}
}

\description{
estimate_para_adaptive_method() function gives the adaptive estimation of parameters of a periodic coefficients regression model.
}
\value{
  \item{beta_ad}{ Parameters to be estimated.}
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
model=lm(y~x1+x2+x3+x4)
z=model$residuals
estimate_para_adaptive_method(n,s,y,x)
}
