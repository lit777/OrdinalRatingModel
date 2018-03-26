# required libraries
library(truncnorm)
library(data.table)

# Data
beamdata <- fread("/users/chanminkim/dropbox/kerrie_R01/data_code/beam_updated_aam_v2_lt3.csv")



# ORDINAL CLASSIFICATIONS W_IJs - THE MAX OF THE LEFT (RECL) AND RIGHT BREAST (RECR) RATINGS
beamdata$ordcat = ifelse(beamdata$RECR > beamdata$RECL,beamdata$RECR,beamdata$RECL)

K <- length(unique(beamdata$ordcat)) # <- Num. of categoris
m <- length(unique(beamdata$CASE)) # <- Num. of cases
n <- length(unique(beamdata$rater)) # <- Num. of raters

D <- beamdata[,(uniquetruth=unique(truth)),by="CASE"]$V1


w <- matrix(ncol=n,nrow=m)
ii <- 0
for(i in unique(beamdata$CASE)){
    jj <- 0
    ii <- ii + 1
    for(j in unique(beamdata$rater)){
        jj <- jj + 1
        w[ii,jj] <- beamdata$ordcat[which(beamdata$CASE==i & beamdata$rater==j)]
    }
}

# prior settings
gamma <- 1
tau <- 0.01

# Initial Values
A <- seq(0, 5, length=K)
A[K] <- Inf # set A[K] <- Inf

z0 <- rep(0, m)
z <- matrix(1, nrow=m, ncol=n)
u <- rnorm(m, 0, 1)

a <- rnorm(n, 0, 1)
a[1] <- 0 # reference rater
b <- rtruncnorm(n, a=0, b=Inf, 0, 1)

mu_a <- 0
tau_a <- 1
mu_b <- 0
tau_b <- 1


##################################
## Main MCMC Function ############
##################################


MCMC <- function(m, n, z0, z, u, w, D, a, b, mu_a, tau_a, mu_b, tau_b, A, gamma = 1, tau = 0.01){

    # sample z_i0

    z0.test <- (D)*rtruncnorm(m, a=0, b=Inf, mean=u, sd=1)+(1-D)*rtruncnorm(m, a=-Inf, b=0, mean=u, sd=1)
 

    # sample z_ij (a_1 <- 0)

    W1 <- which(w==1, arr.ind=T)
    z[W1] <- rtruncnorm(dim(W1)[1], a=-Inf, b= A[1], mean=a[W1[,2]]+b[W1[,2]]*u[W1[,1]], sd=1)
    W2 <- which(w==2, arr.ind=T)
    z[W2] <- rtruncnorm(dim(W2)[1], a=A[1], b= A[2], mean=a[W2[,2]]+b[W2[,2]]*u[W2[,1]], sd=1)
    W3 <- which(w==3, arr.ind=T)
    z[W3] <- rtruncnorm(dim(W3)[1], a=A[2], b= A[3], mean=a[W3[,2]]+b[W3[,2]]*u[W3[,1]], sd=1)
    W4 <- which(w==4, arr.ind=T)
    z[W4] <- rtruncnorm(dim(W4)[1], a=A[3], b= A[4], mean=a[W4[,2]]+b[W4[,2]]*u[W4[,1]], sd=1)
    W5 <- which(w==5, arr.ind=T)
    z[W5] <- rtruncnorm(dim(W5)[1], a=A[4], b= Inf, mean=a[W5[,2]]+b[W5[,2]]*u[W5[,1]], sd=1)


    # sample u_i (!! include z_0 !!)

    u <- rnorm(m,((z-matrix(a, nrow=m, ncol=n, byrow=T))%*%b+z0)/(sum(b^2)+1), sqrt(1/(sum(b^2)+1)))
    u <- (u-mean(u))/sd(u)

    # sample a_j (a_1 : reference group)

    for(j in 2:n){
        a[j] <- rnorm(1,(sum(z[,j]-b[j]*u)+mu_a*tau_a)/(m+tau_a), sqrt(1/(m+tau_a)))
    }

    # if Num. of raters (n) is large, then use the following lines
    # a <- rnorm(n,(colSums(z-u%*%t(b))+mu_a*tau_a)/(m+tau_a), sqrt(1/(m+tau_a)))
    # a[1] <- 0


    # sample b_j
    b <- rtruncnorm(n, a=0, b=Inf, (u%*%(z-matrix(a, nrow=m, ncol=n, byrow=T))+mu_b*tau_b)/(sum(u^2)+tau_b), sqrt(1/(sum(u^2)+tau_b)))
    
    # sample mu_a

    mu_a <- rnorm(1, tau_a*sum(a)/(n*tau_a+tau) , sqrt(1/(n*tau_a+tau)))


    # sample tau_a

    tau_a <- rgamma(1, gamma+n/2, gamma+0.5*sum((a-mu_a)^2))

    # sample mu_b

    mu_b_prop <- rnorm(1, mu_b, 0.1)
    rat1 <- sum(log(dtruncnorm(b, a=0, b=Inf, (u%*%(z-matrix(a, nrow=m, ncol=n, byrow=T))+mu_b_prop*tau_b)/(sum(u^2)+tau_b), sqrt(1/(sum(u^2)+tau_b)))))+dnorm(mu_b_prop, 0, sqrt(1/tau), log=TRUE)+dnorm(mu_b, mu_b_prop, 0.1, log=TRUE)
    rat2 <- sum(log(dtruncnorm(b, a=0, b=Inf, (u%*%(z-matrix(a, nrow=m, ncol=n, byrow=T))+mu_b*tau_b)/(sum(u^2)+tau_b), sqrt(1/(sum(u^2)+tau_b)))))+dnorm(mu_b, 0, sqrt(1/tau), log=TRUE)+dnorm(mu_b_prop, mu_b, 0.1, log=TRUE)
    rat <- rat1 - rat2
    if(is.na(rat) | log(runif(1))>rat){
        mu_b_prop <- mu_b
    }else{
        mu_b <- mu_b_prop
    }


    # sample tau_b

    tau_b_prop <- rgamma(1, 1*tau_b, 1)
    rat1 <- sum(log(dtruncnorm(b, a=0, b=Inf, (u%*%(z-matrix(a, nrow=m, ncol=n, byrow=T))+mu_b*tau_b_prop)/(sum(u^2)+tau_b_prop), sqrt(1/(sum(u^2)+tau_b_prop)))))+dgamma(tau_b_prop, gamma, gamma, log=TRUE)+dgamma(tau_b, 1*tau_b_prop, 1, log=TRUE)
    rat2 <- sum(log(dtruncnorm(b, a=0, b=Inf, (u%*%(z-matrix(a, nrow=m, ncol=n, byrow=T))+mu_b*tau_b)/(sum(u^2)+tau_b), sqrt(1/(sum(u^2)+tau_b)))))+dgamma(tau_b, gamma, gamma, log=TRUE)+dgamma(tau_b_prop, 1*tau_b, 1, log=TRUE)
    rat <- rat1 - rat2
    if(is.na(rat) | log(runif(1))>rat){
        tau_b_prop <- tau_b
    }else{
        tau_b <- tau_b_prop
    }

    # sample A
    U <- NULL
    L <- NULL
    for(k in 1:(K-1)){
        U[k] <- min(z[which(w==k+1, arr.ind=T)])
        L[k] <- max(z[which(w==k, arr.ind=T)])
        A[k] <- rtruncnorm(1, a=L[k], b=U[k], mean=0, sd=100)
    }


    return(list(z0=z0,
                z=z,
                u=u,
                w=w,
                D=D,
                a=a,
                b=b,
                mu_a=mu_a,
                tau_a=tau_a,
                mu_b=mu_b,
                tau_b=tau_b,
                A=A))

}


#################################
update <- list()
update[[1]] <- list(z0=z0,
                     z=z,
                     u=u,
                     w=w,
                     D=D,
                     a=a,
                     b=b,
                     mu_a=mu_a,
                     tau_a=tau_a,
                     mu_b=mu_b,
                     tau_b=tau_b,
                     A=A)

for(t in 2:15000){
  update[[t]] <- MCMC(m=m,
                      n=n,
                      z0=update[[t-1]]$z0,
                      z=update[[t-1]]$z,
                      u=update[[t-1]]$u,
                      w=update[[t-1]]$w,
                      D=update[[t-1]]$D,
                      a=update[[t-1]]$a,
                      b=update[[t-1]]$b,
                      mu_a=update[[t-1]]$mu_a,
                      tau_a=update[[t-1]]$tau_a,
                      mu_b=update[[t-1]]$mu_b,
                      tau_b=update[[t-1]]$tau_b,
                      A=update[[t-1]]$A)
  print(t);print(update[[t]]$tau_b);print(update[[t]]$mu_b);print(update[[t]]$tau_a);print(update[[t]]$mu_a)
}
