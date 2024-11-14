\name{sd_estimation_for_each_s}
\alias{sd_estimation_for_each_s}
\title{Estimating periodic variances in a periodic coefficients regression model}
\usage{
sd_estimation_for_each_s(x,y,s,beta_hat)
}
\arguments{
  \item{x}{A list of independent variables with dimension \eqn{p}.}
  \item{y}{A response variable.}
  \item{s}{A period of the regression model.}
  \item{beta_hat}{The least squares estimation using \link[PerRegMod]{LSE_Reg_per}.}
}

\description{
sd_estimation_for_each_s() function gives the estimation of variances, \eqn{\widehat{\sigma}_s^2=\frac{1}{m-p-1}\sum\limits_{\underset{ }{r=0}}^{m-1}\widehat{\varepsilon}^2_{s+Sr}} for all \eqn{s=1,...,S},in a periodic coefficients regression model.
}
\value{
returns the value of \eqn{\widehat{\sigma}_s^2}.
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
beta_hat=LSE_Reg_per(x,y,s)$beta
sd_estimation_for_each_s(x,y,s,beta_hat)
}
