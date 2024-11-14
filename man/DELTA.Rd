\name{DELTA}
\alias{DELTA}
\title{Calculating the component of vector DELTA}
\usage{
DELTA(x,phi,s,e,sigma)
}
\arguments{
  \item{x}{A list of independent variables with dimension \eqn{p}.}
  \item{phi}{\link[PerRegMod]{phi_n}.}
  \item{s}{A period of the regression model.}
  \item{e}{The residuals vector.}
  \item{sigma}{\link[PerRegMod]{sd_estimation_for_each_s}.}
}

\description{
DELTA() function gives the value of the component of vector DELTA \eqn{\boldsymbol{\Delta}}. See \emph{Regui et al. (2024)} for periodic simple regression model.
\eqn{\mathbf{\Delta}=
\left[\begin{array}{c}
\mathbf{\Delta}_1 \\
 \mathbf{\Delta}_2\\
 \mathbf{\Delta}_3
\end{array}\right]\ }, where \eqn{ \mathbf{\Delta}_1} is a vector of dimension \eqn{S} with component \eqn{
   \frac{n^{\frac{-1}{2} } }{\widehat{ \sigma}_s}\sum\limits_{\underset{ }{r=0}}^{m-1}\widehat{\phi}(Z_{s+Sr,t})},  \eqn{\mathbf{\Delta}_2} is a vector of dimension \eqn{pS} with component \eqn{\frac{ n^{\frac{-1}{2} } }{\widehat{\sigma}_{s}}\sum\limits_{\underset{ }{r=0}}^{m-1} \widehat{\phi}(Z_{s+Sr})K_{s}^{(n)} \mathbf{X}_{s+Sr} },
   \eqn{\mathbf{\Delta}_3} is a vector of dimension \eqn{S} with component \eqn{ \frac{n^{\frac{-1}{2} } }{2\widehat{\sigma}_{s}^{2}}\sum\limits_{\underset{ }{r=0}}^{m-1}{Z_{s+Sr} \widehat{\phi}(Z_{s+Sr})-1 }}.

}
\value{
\item{DELTA()}{ returns the values of \eqn{\mathbf{\Delta}}. See \emph{Regui et al. (2024)} for simple periodic coefficients regression model.}
}

\references{
Regui, S., Akharif, A., & Mellouk, A. (2024). "Locally optimal tests against periodic linear regression in short panels." Communications in Statistics-Simulation and Computation, 1--15. \doi{https://doi.org/10.1080/03610918.2024.2314662}
}



