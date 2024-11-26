\name{phi_n}
\alias{phi_n}
\title{Calculating the value of \eqn{\phi} function }
\usage{
phi_n(x)
}
\arguments{
  \item{x}{A numeric value.}
}

\description{
phi_n() function gives the value of \eqn{\widehat{\phi}(x)=\frac{1}{b^2_n}\frac{\sum\limits_{\underset{ }{s=1}}^{S}\sum\limits_{\underset{}{r=0}}^{m-1}\left(x-Z_{s+Sr}\right)\exp\left(-\frac{\left(x-Z_{s+Sr} \right)^2}{2b_n^2}\right)   }{\sum\limits_{\underset{}{s=1}}^{S}\sum\limits_{\underset{}{r=0}}^{m-1}\exp\left(-\frac{\left(x-Z_{s+Sr} \right)^2}{2b_n^2}\right)  } } with \eqn{b_n=0.002}.
}
\value{
 returns the value of \eqn{\widehat{\phi}(x)}
}



