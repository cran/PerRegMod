\name{lm_per}
\alias{lm_per}
\title{Fitting periodic coefficients regression model by using LSE}
\usage{
lm_per(x,y,s)
}
\arguments{
  \item{x}{A list of independent variables with dimension \eqn{p}.}
  \item{y}{A response variable.}
  \item{s}{A period of the regression model.}
}

\description{
lm_per() function gives the least squares estimation of parameters, intercept \eqn{\mu_s}, slope \eqn{\boldsymbol{\beta}_s}, and standard deviation \eqn{\sigma_s}, of a periodic coefficients regression model using \link[PerRegMod]{LSE_Reg_per} and \link[PerRegMod]{sd_estimation_for_each_s} functions.
\eqn{\widehat{\boldsymbol{\vartheta}}=\left(X^{'}X\right)^{-1}X^{'} Y} where \eqn{X=
\left[\begin{array}{ccccccccccc}
&\mathbf{X}^1_{1}&0&\ldots & 0& &\mathbf{X}^p_{1}&0&\ldots & 0 \\
 & 0&\mathbf{X}^1_{2} &\ldots &0 &  &0&\mathbf{X}^p_{2} &\ldots &0\\

 \textbf{I}_{S}\otimes \mathbf{1}_{m} &0&0&  \ddots&\vdots&\ldots&0&  0&\ddots&\vdots \\

& 0 &0&0 &\mathbf{X}^1_{S}& &0 &0&0 &\mathbf{X}^p_{S}
\end{array}\right]\ },

\eqn{ \mathbf{X}^j_{s}=\left(x^j_{s},...,x^j_{s+(m-1)S}\right)^{'}},
\eqn{Y=(\mathbf{Y}_1^{'},...,\mathbf{Y}_S^{'})^{'}}, \eqn{\mathbf{Y}_{s} =(y_{s},...,y_{(m-1)S+s})^{'}},
\eqn{\mathbf{\epsilon}=(\mathbf{\epsilon}_{1}^{'},...,\mathbf{\epsilon}_{S}^{'})^{'}},
\eqn{\mathbf{\epsilon}_{s} =(\varepsilon_{s},...,\varepsilon_{(m-1)S+s})^{'}}, \eqn{\mathbf{1}_{m}} is a vector of ones of dimension \eqn{m}, \eqn{\textbf{I}_{S}} is the identity matrix of dimension \eqn{S}, \eqn{\otimes} denotes the Kronecker product, and \eqn{\boldsymbol{\vartheta} =\left(\boldsymbol{\mu}^{'} ,{\boldsymbol{\beta}}^{'}\right)^{'}} with \eqn{\boldsymbol{\mu}=(\mu_1,...,\mu_S)^{'}} and \eqn{\boldsymbol{\beta}=(\beta^1_{1},...,\beta^1_{S};...;\beta^p_{1},...,\beta^p_{S})^{'}}.

}
\value{
\item{Residuals}{
the residuals, that is response minus fitted values}
\item{Coefficients}{a named vector of coefficients}
\item{Root mean square error}{The root mean square error}
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
lm_per(x,y,s)
}
