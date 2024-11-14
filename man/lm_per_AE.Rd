\name{lm_per_AE}
\alias{lm_per_AE}
\title{Fitting periodic coefficients regression model by using Adaptive estimation method}
\description{
lm_per_AE() function gives the adaptive estimation of parameters, intercept \eqn{\mu_s}, slope \eqn{\boldsymbol{\beta}_s}, and standard deviation \eqn{\sigma_s}, of a periodic coefficients regression model. \eqn{\widehat{\boldsymbol{\theta}}_{AE} ={\widehat{\boldsymbol{\vartheta} }_{LSE} }+\frac{1}{\sqrt{n}}{\mathbf{\Gamma}}^{-1}\mathbf{\Delta}}.
}

\usage{
lm_per_AE(x,y,s)
}
\arguments{
  \item{x}{A list of independent variables with dimension \eqn{p}.}
  \item{y}{A response variable.}
  \item{s}{A period of the regression model.}
}
\value{
\item{Residuals}{
the residuals, that is response minus fitted values}
\item{Coefficients}{a named vector of coefficients}
\item{Root mean square error}{The root mean square error}

}


\examples{
set.seed(6)
n=200
s=2
x1=rnorm(n,0,1.5)
x2=rnorm(n,0,0.9)
x3=rnorm(n,0,2)
x4=rnorm(n,0,1.9)
y=rnorm(n,0,2.5)
x=list(x1,x2,x3,x4)
lm_per_AE(x,y,s)
}
