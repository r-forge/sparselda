\name{sda}
\alias{sda}
\alias{sda.default}
\title{Sparse discriminant analysis}
\description{Performs sparse linear discriminant analysis. Using an alternating minimization algorithm to minimize the SDA
criterion.}
\usage{
sda (x, ...)

\method{sda}{default}(x,y,lambda=1e-6,stop,maxIte=100,trace=FALSE,tol=1e-6, ...)
}
\arguments{
\item{x}{A matrix of the training data with observations down the rows and variables in the columns.}
\item{y}{A matrix initializing the dummy variables
representing the groups.}

\item{lambda}{The weight on the L2-norm for elastic net
regression. Default: 1e-6.}

\item{stop}{If STOP is negative, its absolute value
corresponds to the desired number of variables. If STOP is
positive, it corresponds to an upper bound on the L1-norm of the b
coefficients. There is a one to one correspondence between stop
and t.}

\item{maxIte}{Maximum number of iterations. Default: 100.}

\item{trace}{If TRUE, prints out its progress. Default:
FALSE.}

\item{tol}{Tolerance for the stopping criterion (change in RSS). Default is 1e-6.}

\item{\ldots}{additional arguments}
}
\value{Returns a list with the following attributes:
\item{beta}{The loadings of the sparse discriminative
directions.}

\item{theta}{The optimal scores.}

\item{rss}{A vector of the Residual Sum of Squares at
each iteration.}
}
\details{
The function finds sparse directions for linear classification.
}
\references{
Clemmensen, L., Hastie, T. and Ersboell, K. (2007) "Sparse discriminant
analysis", Technical report, IMM, Technical University of Denmark
}
\seealso{
\code{\link{normalize}, \link{normalizetest}, \link{smda}}
}
\author{Line Clemmensen}
\examples{
## load data
data(penicilliumYES)

X <- penicilliumYES$X
Y <- penicilliumYES$Y
colnames(Y) <- c("P. Melanoconidium",
                 "P. Polonicum",
                 "P. Venetum")


## test samples
Iout<-c(3,6,9,12)
Iout<-c(Iout,Iout+12,Iout+24)

## training data
Xtr<-X[-Iout,]
k<-3
n<-dim(Xtr)[1]

## Normalize data
Xc<-normalize(Xtr)
Xn<-Xc$Xc
p<-dim(Xn)[2]

## Perform SDA with one non-zero loading for each discriminative
## direction with Y as matrix input
out <- sda(Xn, Y,
           lambda = 1e-6,
           stop = -1,
           maxIte = 25,
           trace = TRUE)

## predict training samples
train <- predict(out, Xn)

## testing
Xtst<-X[Iout,]
Xtst<-normalizetest(Xtst,Xc)

test <- predict(out, Xtst)
print(test$class)

## Factor Y as input
Yvec <- factor(rep(colnames(Y), each = 8))
out2 <- sda(Xn, Yvec,
            lambda = 1e-6,
            stop = -1,
            maxIte = 25,
            trace = TRUE)
}
\keyword{classif}
\keyword{multivariate}
