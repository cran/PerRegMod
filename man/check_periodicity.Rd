\name{check_periodicity}
\alias{check_periodicity}
\title{Checking the periodicity of parameters in the regression model}
\usage{
check_periodicity(x,y,s)
}
\arguments{
  \item{x}{A list of independent variables with dimension \eqn{p}.}
  \item{y}{A response variable.}
  \item{s}{A period of the regression model.}
}

\description{
check_periodicity() function allows to detect the periodicity of parameters in the regression model using \link[PerRegMod]{pseudo_gaussian_test}. See \emph{Regui et al. (2024)} for periodic simple regression model.
\eqn{ T^{(n)}=\left(\mathbf{\Delta}_{1}^{\circ(n)'},\mathbf{\Delta}_{2}^{\circ(n)'},\mathbf{\Delta}_{3}^{\circ(n)'} \right) \left(\begin{array}{ccc}
    \mathbf{\Gamma}^{\circ} _{1} & \mathbf{\Gamma}^{\circ}_{12} & \mathbf{0} \\
    \mathbf{\Gamma}^{\circ}_{12} &\mathbf{\Gamma}^{\circ}_{22}  & \mathbf{0} \\
    \mathbf{0} &\mathbf{0} & \mathbf{\Gamma}^{\circ}_{33}
\end{array}  \right)^{-1} \left(\begin{array}{c}
     \mathbf{\Delta}_{1}^{\circ(n)} \\
     \mathbf{\Delta}_{2}^{\circ(n)}\\
     \mathbf{\Delta}_{3}^{\circ(n)}
\end{array} \right)},
where
\eqn{\boldsymbol{\Delta}_{1}^{\circ(n)}= n^{\frac{-1}{2}} \sum\limits_{\underset{ }{r=0}}^{m-1} \left(\begin{array}{c}
 \widehat{\phi}(Z_{1+Sr})-\widehat{\phi}(Z_{S+Sr})
 \\
 \vdots\\
  \widehat{\phi}(Z_{S-1+Sr})-\widehat{\phi}(Z_{S+Sr})
 \end{array} \right)},

  \eqn{ \mathbf{\Delta}_{2}^{\circ(n)}= \frac{n^{\frac{-1}{2}}}{2\widehat{\sigma} }\sum\limits_{\underset{ }{r=0}}^{m-1} \left(\begin{array}{c}
   \widehat{\psi}(Z_{1+Sr})-  \widehat{\psi}(Z_{S+Sr}) \\
  \vdots\\
\widehat{\psi}(Z_{S-1+Sr})- \widehat{\psi}(Z_{S+Sr})  \\
  \end{array}\right)},

\eqn{\mathbf{\Delta}_{3}^{\circ(n)}=n^{\frac{-1}{2}} \sum\limits_{\underset{ }{r=0}}^{m-1} \left(
  \begin{array}{c}
 \widehat{\phi}(Z_{1+Sr}) \mathbf{K}_1^{(n)}\mathbf{X}_{1+Sr}- \widehat{\phi}(Z_{S+Sr}) \mathbf{K}_S^{(n)}\mathbf{X}_{S+Sr}\\ \vdots\\
 \widehat{\phi}(Z_{S-1+Sr})\mathbf{K}_{S-1}^{(n)}\mathbf{X}_{S-1+Sr}- \widehat{\phi}(Z_{S+Sr})\mathbf{K}_S^{(n)}\mathbf{X}_{S+Sr}
 \end{array} \right)},
\eqn{\mathbf{\Gamma}^{\circ} _{11}=\frac{\widehat{I}_n }{S} \Sigma },  \eqn{\mathbf{\Gamma}^{\circ} _{22}=\dfrac{\widehat{I}_n}{4S\widehat{\sigma}^2}
\Sigma},  \eqn{\mathbf{\Gamma}^{\circ} _{12}=\frac{ \widehat{N}_n }{2S\widehat{\sigma}} \Sigma}, and
\eqn{\mathbf{\Gamma}^{\circ} _{33}=\frac{\widehat{I}_n }{S} \Sigma \otimes \mathbf{I}_{p\times p}}
with
\eqn{\widehat{I}_n=\frac{1}{nT}\sum\limits_{\underset{ }{s=1}}^{S}\sum\limits_{\underset{}{r=0}}^{m-1}{\widehat{\phi}^{2}\left(\frac{\widehat{Z}_{s+Sr}}{\widehat{
 \sigma}_s} \right)}}, \eqn{\widehat{N}_n=\frac{1}{nT}\sum\limits_{\underset{ }{s=1}}^{S}\sum\limits_{\underset{ }{r=0}}^{m-1}{\widehat{\phi}}^{2}\left( \frac{\widehat{Z}_{s+Sr}}{\widehat{
 \sigma}_s}\right)\frac{\widehat{Z}_{s+Sr}}{\widehat{
 \sigma}_s}},


  \eqn{ \Sigma=\left[\begin{array}{cccc}
2 & 1& \ldots&1 \\
1&\ddots & \ddots& \vdots\\
  \vdots& \ddots &\ddots & 1 \\
 1&\ldots &1 & 2
\end{array}\right]\ },
 \eqn{Z_{s+Sr}=\frac{y_{s+Sr}-\widehat{\mu}_s-\sum\limits_{\underset{}{j=1}}^{p}\widehat{\beta}^j_{s}x^j_{s+Sr}}{\widehat{\sigma}_s}}, \eqn{\mathbf{ X}_{s+Sr}=\left(x^1_{s+Sr},...,x^p_{s+Sr}  \right)^{'}}, \eqn{  \mathbf{K}^{(n)}_{s}=\left[\begin{array}{ccc}
\overline{(x^1_{s})^2 } & &\overline{x^i_{s}x^j_{s} }\\
  &\ddots & \\
\overline{x^j_{s}x^i_{s} } & &\overline{(x^p_{s})^2 }
\end{array}\right]^{\frac{-1}{2} } },

\eqn{\overline{x^i_{s}x^j_{s} } =\frac{1}{m}\sum\limits_{\underset{ }{r=0}}^{m-1}{x^i_{s+Sr}x^j_{s+Sr}}},   \eqn{\overline{(x^i_{s})^2 } =\frac{1}{m}\sum\limits_{\underset{ }{r=0}}^{m-1}{(x^i_{s+Sr})^2 }}, \eqn{\widehat{\psi}(x)=x\widehat{\phi}(x)-1}, and

  \eqn{   \widehat{\phi}(x)=\frac{1}{b^2_n}\frac{\sum\limits_{\underset{ }{s=1}}^{S}\sum\limits_{\underset{}{r=0}}^{m-1}\left(x-Z_{s+Sr}\right)\exp\left(-\frac{\left(x-Z_{s+Sr} \right)^2}{2b_n^2}\right)   }{\sum\limits_{\underset{}{s=1}}^{S}\sum\limits_{\underset{}{r=0}}^{m-1}\exp\left(-\frac{\left(x-Z_{s+Sr} \right)^2}{2b_n^2}\right)  }} with \eqn{b_n\rightarrow 0}.

}

\references{
Regui, S., Akharif, A., & Mellouk, A. (2024). "Locally optimal tests against periodic linear regression in short panels." Communications in Statistics-Simulation and Computation, 1--15. \doi{https://doi.org/10.1080/03610918.2024.2314662}
}
\value{
\item{check_periodicity()}{returns the value of observed statistic,  \eqn{T^{(n)}}, degrees of freedom, \eqn{(S-1)\times(p+2)},  and p-value}
}

\examples{
library(expm)
set.seed(6)
n=400
s=4
x1=rnorm(n,0,1.5)
x2=rnorm(n,0,0.9)
x3=rnorm(n,0,2)
x4=rnorm(n,0,1.9)
y=rnorm(n,0,2.5)
x=list(x1,x2,x3,x4)
check_periodicity(x,y,s)
}
