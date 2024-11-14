\name{A_x_B}
\alias{A_x_B}
\title{A Kronecker product B}
\usage{
A_x_B(A,B)
}
\arguments{
  \item{A}{A matrix.}
  \item{B}{A matrix.}
}

\description{
A_x_B() function gives A Kronecker product B
}
\value{\item{A_x_B(A,B)}{returns the matrix A Kronecker product B, \eqn{A\otimes B } } }
\examples{
A=matrix(rep(1,6),3,2)
B=matrix(seq(1,8),2,4 )
A_x_B(A,B)

}
