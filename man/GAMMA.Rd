\name{GAMMA}
\alias{GAMMA}
\title{Calculating the component of matrix GAMMA}
\usage{
GAMMA(x,phi,s,z,sigma)
}
\arguments{
  \item{x}{A list of independent variables with dimension \eqn{p}.}
  \item{phi}{\link[PerRegMod]{phi_n}.}
  \item{s}{A period of the regression model.}
  \item{z}{The residuals vector.}
  \item{sigma}{\link[PerRegMod]{sd_estimation_for_each_s}.}
}

\description{
GAMMA() function gives the value of the component of matrix GAMMA \eqn{\boldsymbol{\Gamma}}. See \emph{Regui et al. (2024)} for periodic simple regression model.
\eqn{\mathbf{\Gamma}=\frac{1}{S}
\left[\begin{array}{ccc}
\left(\mathbf{\Gamma}_{11}\right)_{S \times S }&\mathbf{0} & \mathbf{\Gamma}_{13}  \\
\mathbf{0} &\left(\mathbf{\Gamma}_{22} \right)_{pS\times pS } &\mathbf{0} \\
\mathbf{\Gamma}_{13} & \mathbf{0}& \left(\mathbf{\Gamma}_{33} \right)_{S\times S}
\end{array}\right]\ }, where \eqn{\mathbf{\Gamma}_{11}=\widehat{I}_{n}\text{diag}(\frac{1}{\widehat{\sigma}_{1}^{2}},...,\frac{1}{\widehat{\sigma}_{S}^{2}} )}, \eqn{\mathbf{\Gamma}_{13}=\frac{\widehat{N}_{n}}{2}\text{diag}(\frac{1}{\widehat{\sigma}_{1}^{3}},...,\frac{1}{\widehat{\sigma}_{S}^{3}} )},
 \eqn{
\mathbf{\Gamma}_{22}=\widehat{I}_{n}\text{diag}(\frac{1}{\widehat{\sigma}_{1}^{2}},...,\frac{1}{\widehat{\sigma}_{S}^{2}} ) \otimes \mathbf{I}_{p}},
\eqn{\mathbf{\Gamma}_{33}=\frac{\widehat{J}_{n}}{4}\text{diag}(\frac{1}{\widehat{\sigma}_{1}^{4}},...,\frac{1}{\widehat{\sigma}_{S}^{4}} )}, \eqn{\widehat{I}_n=\frac{1}{nT}\sum\limits_{\underset{ }{s=1}}^{S}\sum\limits_{\underset{}{r=0}}^{m-1}{\widehat{\phi}^{2}\left(\frac{\widehat{Z}_{s+Sr}}{\widehat{
 \sigma}_s} \right)}}, \eqn{\widehat{N}_n=\frac{1}{nT}\sum\limits_{\underset{ }{s=1}}^{S}\sum\limits_{\underset{ }{r=0}}^{m-1}{\widehat{\phi}}^{2}\left( \frac{\widehat{Z}_{s+Sr}}{\widehat{
 \sigma}_s}\right)\frac{\widehat{Z}_{s+Sr}}{\widehat{
 \sigma}_s}}, \eqn{\widehat{J}_n=\frac{1}{nT}\sum\limits_{\underset{ }{s=1}}^{S}\sum\limits_{\underset{}{r=0}}^{m-1}\widehat{\phi}^{2}\left( \frac{\widehat{Z}_{s+Sr}}{\widehat{
 \sigma}_s}\right)\left( \frac{\widehat{Z}_{s+Sr}}{\widehat{
 \sigma}_s}\right)^{2}-1}, and

    \eqn{ \widehat{\phi}(x)=\frac{1}{b^2_n}\frac{\sum\limits_{\underset{ }{s=1}}^{S}\sum\limits_{\underset{}{r=0}}^{m-1}\left(x-Z_{s+Sr}\right)\exp\left(-\frac{\left(x-Z_{s+Sr} \right)^2}{2b_n^2}\right)   }{\sum\limits_{\underset{}{s=1}}^{S}\sum\limits_{\underset{}{r=0}}^{m-1}\exp\left(-\frac{\left(x-Z_{s+Sr} \right)^2}{2b_n^2}\right)  } \text{ with }b_n\rightarrow 0}.

}

\value{
\item{GAMMA()}{returns the matrix \eqn{\mathbf{\Gamma}}. See \emph{Regui et al. (2024)} for simple periodic coefficients regression model.}
}
\references{
Regui, S., Akharif, A., & Mellouk, A. (2024). "Locally optimal tests against periodic linear regression in short panels." Communications in Statistics-Simulation and Computation, 1--15. \doi{https://doi.org/10.1080/03610918.2024.2314662}
}



