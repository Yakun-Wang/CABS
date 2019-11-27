
library(pracma)
library(trust)
library(mvtnorm)
library(MASS)
library(NPflow)
library(plotly)


#########################
#   1. simulate data    #
#########################

# abrupt

library(forecast)
set.seed(1001)
L=8
Nobs=1000
w_l=rep(0,8)                                 ## u
ts.sim1=matrix(rep(0,L*Nobs/2),L,Nobs/2)
ts.sim2=matrix(rep(0,L*Nobs/2),L,Nobs/2)


for (l in 1: L)
{
  w_l[l]=(l-1)/(L-1)
}
t=seq(1,1000,1)

for(l in 1:L/2)
{
  ts.sim1[l,]=arima.sim(list(order=c(1,0,0),ar=-0.5),n=Nobs/2)
  ts.sim2[l,]=arima.sim(list(order=c(1,0,0),ar=0.5),n=Nobs/2)
}
for (l in ((L/2)+1):L)
{
  ts.sim1[l,]=arima.sim(list(order=c(1,0,0),ar=-0.9),n=Nobs/2)
  ts.sim2[l,]=arima.sim(list(order=c(1,0,0),ar=0.9),n=Nobs/2)
}
ts.sim=cbind(ts.sim1,ts.sim2)              ## x
#ts.plot(ts.sim)

x=t(ts.sim)
u=w_l

# standardized

for (i in 1:L)
{
  xmat=cbind(matrix(1,dim(x)[1],1), matrix(seq(1,dim(x)[1],1),dim(x)[1],1))
  linfit=solve(t(xmat)%*%xmat)%*%t(xmat)%*%x[,i]
  x[,i]=x[,i]-xmat%*%linfit
}
ts.plot(x[,1])


# slowly varying data
nrep=8
nobs=1000
sd_true=1
e=matrix(rnorm(nrep*nobs,0,sd_true),nrep,nobs)
x=matrix(0,nrep,nobs)


w_l=rep(0,8)
for (l in 1: nrep)
{
  w_l[l]=(l-1)/(nrep-1)
}
u=w_l


for (i in 1:nrep)
{
  if (u[i] <= 0.5)
  { # slowly varying -0.5 to 0.5
    phi_true=-0.5+(1:nobs)/nobs
    for (j in 1:nobs)
      if (j==1)
      {
        x[i,j]=e[i,j]
      }else
      {
        x[i,j]=phi_true[j]*x[i,j-1]+e[i,j]
      }
  }else
  {
    # slowly varying from -0.9 to 0.9
    phi_true=-0.9+((1:nobs)/nobs)*1.8
    for (j in 1:nobs)
    {
      if (j==1)
      {
        x[i,j]=e[i,j]
      }else
      {
        x[i,j]=phi_true[j]*x[i,j-1]+e[i,j]
      }
    }
  }
}
x=t(x)

# standardized

for (i in 1:nrep)
{
  xmat=cbind(matrix(1,dim(x)[1],1), matrix(seq(1,dim(x)[1],1),dim(x)[1],1))
  linfit=solve(t(xmat)%*%xmat)%*%t(xmat)%*%x[,i]
  x[,i]=x[,i]-xmat%*%linfit
}
ts.plot(x)

#########################
# 2. call CABS function #
#########################
source('~/Desktop/R_code/CABS.R')
cabs(x,u,nloop=5000,seed=6666,nwarmup=2000,nexp_tmax=10,nexp_umax=8,tmin=50,umin=1,nbasis=7,
     Rev_Jump_t=1,Rev_Jump_u=1,prob_mm1=0.8,prob_mm1_u=0.8,sigmasqalpha=100,tau_prior_a=-1,
     tau_prior_b=0,var_inflate=1.0,plotting=TRUE)


## Calculation of spectrum
spec_est <- array(0, c(51, dim(x)[1], dim(x)[2]))

for (k in 1:dim(x)[2])
{
  for (p in (nwarmup + 1):nloop)
  {
    xi_t_curr <- xi[[nexp_tcurr[p + 1]]][, (p + 1)]
    nseg_t_curr <- nseg_tcurr[[nexp_tcurr[p + 1]]][, (p + 1)]
    log_spec_hat <- log_spec_hat_curr[[nexp_tcurr[p + 1], nexp_ucurr[p +
                                                                       1]]][, , which(u[k] <= ui[[nexp_tcurr[p + 1]]][, (p + 1)])[1], (p +
                                                                                                                                              1)]
    for (g in 1:nexp_tcurr[p + 1])
    {
      if (g == 1)
      {
        spec_est[, 1:xi_t_curr[g], k] <- spec_est[, 1:xi_t_curr[g], k] + repmat(as.matrix(log_spec_hat[, g]), 1, nseg_t_curr[g]) /
          (nloop - nwarmup)
        
      } else
      {
        spec_est[, (xi_t_curr[g - 1] + 1):xi_t_curr[g], k] <- spec_est[, (xi_t_curr[g -
                                                                                      1] + 1):xi_t_curr[g], k] +
          repmat(as.matrix(log_spec_hat[, g]), 1, nseg_t_curr[g]) / (nloop -
                                                                       nwarmup)
      }
      
    }
  }
}

# plot estimation of spectrum of the first observation
  spec_mean <-  apply(spec_est[, , 1], c(1, 2), mean)
  plot_ly(z = ~ spec_mean) %>% add_surface() %>% layout(title = paste(
    "Estimated Log Spectrum for data with the
    first observation"))
  
# Time partition plots for a partition of 2
  kk=which(nexp_tcurr[1:nloop]==2)
  if ( length(kk) != 0 )
  {
    for (g in 1:1)
    {
      
      plot(xi[[2]][g, kk], ylim=c(490,510),type = "l",
           main = "Time partition plots for a partition of 2")
      
    }
  }

