\name{pseudo_gaussian_test}
\alias{pseudo_gaussian_test}
\title{Detecting periodicity of parameters in the regression model}
\usage{
pseudo_gaussian_test(x,z,s)
}
\arguments{
  \item{x}{A list of independent variables with dimension \eqn{p}.}
  \item{z}{The residuals vector.}
  \item{s}{A period of the regression model.}
}

\description{
pseudo_gaussian_test() function gives the value of the statistic test, \eqn{T^{(n)}}, for detecting periodicity of parameters in the regression model. See \link[PerRegMod]{check_periodicity} function.
}
\value{
returns the value of the statistic test, \eqn{T^{(n)}}.
}

