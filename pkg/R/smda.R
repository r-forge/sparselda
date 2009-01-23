.packageName <- "sparseLDA"

smda <- function(X,Z,Rj,lambda=1e-6,stop,maxIte=50,trace=FALSE,tol=1e-4){
#
# smda performs Sparse Mixture Disciminant Analysis
# Solving: argmin{|(Y*theta-X*b)|_2^2 + t*|beta|_1 + lambda*|beta|_2^2}
#
# INPUT:
# X      : matrix of n observations down the rows and p variable columns. The
#          columns are assumed normalized
# Z      : matrix initializing the probabilities representing the groups
# Rj     : K length vector containing the number of subclasses in each of
#          the K classes
# lambda : the weight on the L2-norm for elastic net regression. Default: 1e-6
# stop   : nonzero STOP will perform
#          elastic net regression with early stopping. If STOP is negative, its 
#          absolute value corresponds to the desired number of variables. If STOP
#          is positive, it corresponds to an upper bound on the L1-norm of the
#          b coefficients. There is a one to one correspondence between stop
#          and t.
# maxIte : Maximum number of iterations. Default: 50.
# trace  : trace = FALSE turns printing of RSS off and trace = TRUE turns it on.
# tol    : Tolerance for the stopping criterion (change in RSS). Default: 1e-4
#
# OUTPUT:
# $beta   : The regression parameters
# $theta  : Optimal scores
# $Z      : Updated subclass probabilities
# $rss    : Residual Sum of Squares at each itearation
#
# Author: Line H. Clemmensen, IMM, DTU, lhc@imm.dtu.dk
# Based on the elastic net algorithm by Hui Zou and Trevor Hastie
#

N <- dim(X)[1]
p <- dim(X)[2]
K <- length(Rj) # number of classes
R <- dim(Z)[2] # number of subclasses
RSSold <- 1e8 
RSS <- 1e6
ite <- 0
Zhat <- matrix(0,N,R-1)
Dp <- apply(Z,2,sum)
Dp_inv <- diag(1/sqrt(Dp/N)) # R x R
theta <- 1/sum(diag(Dp/N))*diag(rep(1,R))[,1:(R-1)]/R
Ztheta <- Z%*%theta  # N x R-1
rss <- rep(0,maxIte)
b <- matrix(0,p,R-1)
if (length(stop)< (R-1)){
    stop <- rep(stop[1],1,R-1)
}
if (stop[1]<0) sparse <- "varnum" else sparse <- "penalty" 


while (abs(RSSold-RSS)/RSS > tol & ite < maxIte){ 
    RSSold <- RSS
    ite <- ite + 1
    # 1. Estimate beta:    
    for (j in 1:(R-1)){
        Zc <- Ztheta[,j]
        beta<- solvebeta(X, Zc, paras=c(lambda, abs(stop[j])),sparse=sparse)
        b[,j] <- t(beta)
        Zhat[,j] <- X%*%b[,j]
    }    

    # 2. Optimal scores: (balanced Procrustes problem)
    B <- t(Z)%*%Zhat
    sb <- svd(B,nu=R-1,nv=R-1)
    theta.old <- theta
    theta <- Dp_inv%*%sb$u%*%t(sb$v)
    Ztheta <- Z%*%theta
    RSS <- sum((Ztheta-Zhat)*(Ztheta-Zhat))
    rss[ite] <- RSS
    if (trace){
    cat('ite: ', ite, ' RSS: ', RSS,'\n')
    }

    # 3. update parameter estimates:
    Sigma <- matrix(0,R-1,R-1)
    mu <- matrix(0,(R-1)*R,K)
    dim(mu) <- c(R-1,R,K)
    for (i in 1:K){
        IK <- (sum(Rj[1:i-1])+1):(sum(Rj[1:i-1])+Rj[i])
        Ik <- apply(Z[,IK]>0,1,any)
        Ik.length <- sum(Ik)
        for (j in 1:Rj[i]){
            mu[,IK[j],i] = apply(matrix(1,Ik.length,1)%*%t(Z[Ik,IK[j]])%*%Zhat[Ik,],2,sum)/Dp[IK[j]]
            Sigma = Sigma + t(Zhat[Ik,]-matrix(1,Ik.length,1)%*%t(matrix(mu[,IK[j],i])))%*%(Z[Ik,IK[j]]%*%matrix(1,1,Ik.length))%*%(Zhat[Ik,]-
                    matrix(1,Ik.length,1)%*%t(matrix(mu[,IK[j],i])))/(Ik.length-Rj[i])
        }
    }
    Sigma_inv <- solve(Sigma + 1e-2*diag(rep(1,R-1)))

    for (i in 1:K){
        IK <- (sum(Rj[1:i-1])+1):(sum(Rj[1:i-1])+Rj[i])
        Ik <- apply(Z[,IK]>0,1,any)
        Ik.length <- sum(Ik)
        Dmahal_K <- matrix(0,Ik.length,Rj[i])
        for (j in 1:Rj[i]){
        Dmahal_K[,j] <- diag((Zhat[Ik,]-matrix(1,Ik.length,1)%*%t(matrix(mu[,IK[j],i])))%*%Sigma_inv%*%t(Zhat[Ik,]-
                        matrix(1,Ik.length,1)%*%t(matrix(mu[,IK[j],i]))))
        }
        sum_K <- apply(Z[Ik,IK]*exp(-Dmahal_K/2),1,sum)
        for (j in 1:Rj[i]){
            Z[Ik,IK[j]] <- Z[Ik,IK[j]]*exp(-Dmahal_K[,j]/2)/(sum_K+1e-6)
        }
        Z[Ik,IK] <- Z[Ik,IK]/(apply(Z[Ik,IK],1,sum)*rep(1,1,Rj[i]))
    }
    Ztheta <- Z%*%theta
    Dp <- apply(Z,2,sum)
    Dp_inv <- diag(1/sqrt(Dp/N)) # R x R
}

# Remove trivial directions
Ik <- sb$d > 1e-6
M <- sum(Ik)
theta <- theta[,1:M]
Ztheta <- Z%*%theta
b <- b[,1:M]
Zhat <- Zhat[,1:M]
for (j in 1:M){
     Zc <- Ztheta[,j]
     beta<- solvebeta(X, Zc, paras=c(lambda, abs(stop[j])),sparse=sparse)
     b[,j] <- t(beta)
     Zhat[,j] <- X%*%b[,j]
}
if (trace){
  RSS <- sum((Ztheta-Zhat)*(Ztheta-Zhat))
  cat('final update, RSS: ', RSS,'\n')
}
list(beta=b,theta=theta,rss=rss[1:ite])
}

