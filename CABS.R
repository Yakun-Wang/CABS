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

########################
#   2. main function   #
########################
#   seed: random number seed
#   nloop: number of total MCMC iterations 
#   nwarmup: number of iterations to use as burn-in
#	nexp_tmax: maximum number of time-based partition segments(if
#                Rev_Jump_t=0, this will be the exact number of time-
#                based partition segments since no birth and death 
#                moves are allowed)
# nexp_umax: maximum number of covariate-based partition segments(if
#            Rev_Jump_u=0, this will be the exact number of covariate-
#            based partition segments since no birth and death 
#            moves are allowed)
# Rev_Jump_t: 1 if want to have reversible jumps in time partition
#             and 0 otherwise (i.e. 1 enables birth and death moves
#             while 0 limits you only to within model moves) 
# Rev_Jump_u: 1 if want to have reversible jumps in covariate
#             partition and 0 otherwise (i.e. 1 enables birth and 
#             death moves while 0 limits you only to within model moves) 
# tmin: minimum number of time observations in each time segment
# umin: minimum number of unique covariate values in each covariate segment
# prob_mm1: probability of a small move (+/- 1 value) in
#          time-based within move step
# prob_mm1_u: probability of a small move (+/- 1 value) in
#          covariate-based within move step
# nbasis: number of basis functions to use for spectral estimation  
# sigmasqalpha: prior known variance for constant term in spectral estimation
# tau_prior_a: tau prior distribution parameter 
# tau_prior_b: tau prior distribution parameter
# var_inflate: inflation factor for variance of beta distribution


library(pracma)
library(trust)
library(mvtnorm)
library(MASS)
library(NPflow)
library(plotly)

cabs <-
  function(x,
           u,
           nloop,
           seed,
           nwarmup,
           nexp_tmax,
           nexp_umax,
           tmin,
           umin,
           nbasis,
           Rev_Jump_t,
           Rev_Jump_u,
           prob_mm1,
           prob_mm1_u,
           sigmasqalpha,
           tau_prior_a,
           tau_prior_b,
           var_inflate,
           plotting)
  {
    lin_basis_func <-
      function(freq, nbeta)
        # Function calculates basis function values for specified frequencies
      {
        n <- length(freq)
        omega <- matrix(1, n, nbeta)
        for (j in 2:nbeta)
        {
          omega[, j] <- sqrt(2) * cos((j - 1) * pi * freq) / (pi * (j - 1))
        }
        
        return(omega)
        
      }
    
    whittle_like <- function(y, fhat, nseg) {
      nfreq <- floor(nseg / 2)
      log_prop_spec_dens <- 0
      
      
      if (nseg %% 2 == 1)
      {
        for (i in 1:dim(y)[2])
        {
          log_prop_spec_dens <- log_prop_spec_dens - sum(fhat[2:(nfreq + 1)] + exp(y[2:(nfreq +
                                                                                          1), i] - fhat[2:(nfreq + 1)])) -
            0.5 * (fhat[1] + exp(y[1, i] - fhat[1])) - 0.5 * nfreq * log(2 *
                                                                           pi)
          
        }
        
      } else
      {
        for (i in 1:dim(y)[2])
        {
          log_prop_spec_dens <- log_prop_spec_dens - sum(fhat[2:nfreq] + exp(y[2:nfreq, i] -
                                                                               fhat[2:nfreq])) -
            0.5 * (fhat[1] + exp(y[1, i] - fhat[1])) - 0.5 * (fhat[nfreq +
                                                                     1] + exp(y[nfreq + 1, i] - fhat[nfreq + 1])) -
            0.5 * nfreq * log(2 * pi)
        }
      }
      
      return(log_prop_spec_dens)
    }
    
    whittle_derivs2 <-
      function(param,
               n,
               nu_mat,
               ytemp,
               tau_temp,
               nbeta,
               nbasis,
               sigmasqalpha)
      {
        ydev <- ytemp - repmat(nu_mat %*% param, 1, dim(ytemp)[2]) # Apply element-wise operation to two arrays with implicit expansion enabled
        # ydev=logY-logf
        # ytemp=yy: log periodograms
        # n is the segment length in the time domain
        n1 <- floor(n / 2)
        if (n %% 2 == 1)
          # odd n
        {
          f <- sum(sum(repmat(nu_mat[2:(n1 + 1), ] %*% param, 1, dim(ytemp)[2]) + exp(ydev[2:(n1 +
                                                                                                1), ]))) +
            0.5 * (sum(as.vector(nu_mat[1, ] %*% param) + exp(ydev[1, ]))) +
            0.5 * (t(param[2:nbeta]) %*% param[2:nbeta] / tau_temp + param[1] ^
                     2 / sigmasqalpha)
          
        } else
        {
          f <- sum(sum(repmat(nu_mat[2:n1, ] %*% param, 1, dim(ytemp)[2]) + exp(ydev[2:n1, ]))) +
            0.5 * (sum(as.vector(nu_mat[1, ] %*% param) + exp(ydev[1, ]))) +
            0.5 * (sum(as.vector(nu_mat[n1 + 1, ] %*% param) + exp(ydev[n1 +
                                                                          1, ]))) +
            0.5 * (t(param[2:nbeta]) %*% param[2:nbeta] / tau_temp + param[1] ^
                     2 / sigmasqalpha)
        }
        
        g <- matrix(0, nbeta, 1)
        g[1, ] <- param[1] / sigmasqalpha
        g[2:nbeta, ] <- param[2:nbeta] / tau_temp
        
        h <- matrix(0, nbeta, nbeta)
        h[1, 1] <- 1 / sigmasqalpha
        h[2:nbeta, 2:nbeta] <- 1 / tau_temp * diag(nbasis)
        
        
        if (n %% 2 == 1)
        {
          for (i in 1:dim(ydev)[2])
          {
            temp_mat <- nu_mat[2:(n1 + 1), ] * repmat((1 - exp(as.matrix(ydev[2:(n1 +
                                                                                   1), i]))), 1, dim(nu_mat)[2])
            g <- g + as.matrix(apply(temp_mat, 2, sum)) + t(0.5 * (1 - exp(ydev[1, i])) %*%
                                                              nu_mat[1, ])
          }
          
          jj <- seq(1:nbeta)
          
          
          big_mat <- repmat(t(nu_mat[2:(n1 + 1), ]), nbeta, 1) * t(nu_mat[2:(n1 +
                                                                               1), repmat(jj, nbeta, 1)])
          for (i in 1:dim(ydev)[2])
          {
            coef_mat <- repmat(exp(t(ydev[2:(n1 + 1), i])), nbeta ^ 2, 1)
            h <- h + apply(array(big_mat * coef_mat, dim = c(nbeta, nbeta, n1)), c(1, 2), sum) +
              t(0.5 * exp(ydev[1, i]) %*% nu_mat[1, ]) %*% nu_mat[1, ]
          }
          
        } else
        {
          for (i in 1:dim(ydev)[2])
          {
            temp_mat <- nu_mat[2:n1, ] * repmat((1 - exp(as.matrix(ydev[2:n1, i]))), 1, dim(nu_mat)[2])
            g <- g + as.matrix(apply(temp_mat, 2, sum)) + t(0.5 * (1 - exp(ydev[1, i])) %*%
                                                              nu_mat[1, ]) +
              t(0.5 * (1 - exp(ydev[n1 + 1, i])) %*% nu_mat[n1 + 1, ])
          }
          
          jj <- seq(1:nbeta)
          big_mat <- repmat(t(nu_mat[2:n1, ]), nbeta, 1) * t(nu_mat[2:n1, repmat(jj, nbeta, 1)])
          for (i in 1:dim(ydev)[2])
          {
            coef_mat <- repmat(exp(t(ydev[2:n1, i])), nbeta ^ 2, 1)
            h <- h + apply(array(big_mat * coef_mat, dim = c(nbeta, nbeta, n1 -
                                                               1)), c(1, 2), sum) +
              t(0.5 * exp(ydev[1, i]) %*% nu_mat[1, ]) %*% nu_mat[1, ] +
              t(0.5 * exp(ydev[n1 + 1, i]) %*% nu_mat[n1 + 1, ]) %*% nu_mat[n1 +
                                                                              1, ]
          }
          
        }
        list(value = f,
             gradient = g,
             hessian = h)
      }
    
    
    postbeta <-
      function(j,
               i,
               nseg_time_temp,
               x,
               u,
               xi_temp,
               ui_temp,
               tau_temp,
               nbasis,
               sigmasqalpha)
      {
        
        
        # i is the index for time segment
        # j is the index for covariate segment
        
        nfreq <- floor(nseg_time_temp / 2)
        freq <- (0:nfreq) / (2 * nfreq)
        
        # create set of indices for replications in jth covariate segment
        if (i == 1) {
          uu <- which(u <= ui_temp[i])   
        } else {
          uu <- which(u <= ui_temp[i] & u > ui_temp[i - 1])
        }
        
        # create log periodograms for replications in ith time seg and jth cov seg
        y <- matrix(0, nseg_time_temp, length(uu))
        yy <- matrix(0, nfreq + 1, length(uu))
        jj <- 0
        for (k in uu[1]:uu[length(uu)])
        {
          jj <- jj + 1
          if (j > 1) {
            y[, jj] <- log(abs(fft(x[(xi_temp[j - 1] + 1):xi_temp[j], k])) ^ 2 / nseg_time_temp)
            yy[, jj] <- y[1:(nfreq + 1), jj]
          } else
            # j==1
          {
            y[, jj] <- log(abs(fft(x[1:xi_temp[j], k])) ^ 2 / nseg_time_temp)
            yy[, jj] <- y[1:(nfreq + 1), jj]
          }
        }
        
        
        # pass log periodograms into optimizer to obtain beta_mean and beta_var for normal approximation to beta posterior
        nbeta <- nbasis + 1
        nu_mat <- lin_basis_func(freq, nbeta) # basis function
        nn <- nseg_time_temp
        ytemp <- yy
        param <- rep(0, nbeta)
        
        
        post <- trust(
          whittle_derivs2,
          param,
          rinit = 1,
          rmax = 100,
          parscale = rep(1, nbeta),
          iterlim = 100,
          fterm = sqrt(.Machine$double.eps),
          mterm = sqrt(.Machine$double.eps),
          minimize = TRUE,
          blather = FALSE,
          nn,
          nu_mat,
          yy,
          tau_temp,
          nbeta,
          nbasis,
          sigmasqalpha
        )
        beta_mean <- post$argument
        beta_var <-  solve(post$hessian)
        list(
          beta_mean = beta_mean,
          beta_var = beta_var,
          nu_mat = nu_mat,
          y = y
        )
        
      }
    
    move <- function(kk, nexp_curr, nexp_max)
    {
      if (kk == 0)
        # Stay where you (if nexp_curr=1) or join segments if there are no available segments to cut
      {
        if (nexp_curr == 1)
          #  nexp_curr=number of segments
        {
          nexp_prop <- nexp_curr # Stay where you are
          log_move_prop <- 0
          log_move_curr <- 0
        } else
        {
          nexp_prop <- nexp_curr - 1  # join segments,so number of segments is reduced
          log_move_prop <- 0
          if (nexp_prop == 1)
          {
            log_move_curr <- 0
          } else
          {
            log_move_curr <- log(0.5)
          }
        }
      } else
      {
        if (nexp_curr == 1)
        {
          nexp_prop <- nexp_curr + 1
          log_move_prop <- 0
          if (nexp_prop == nexp_max)
          {
            log_move_curr <- 0
          } else
          {
            log_move_curr <- log(0.5)
          }
          
        } else if (nexp_curr == nexp_max)
        {
          nexp_prop <- nexp_curr - 1
          log_move_prop <- 0
          if (nexp_prop == 1)
          {
            log_move_curr <- 0
          } else
          {
            log_move_curr <- log(0.5)
          }
          
        } else
        {
          r <- runif(1, 0, 1)
          if (r < 0.5)
          {
            nexp_prop <- nexp_curr + 1
            if (nexp_prop == nexp_max)
            {
              log_move_curr <- 0
              
              log_move_prop = log(0.5)
            } else
            {
              log_move_curr <- log(0.5)
              log_move_prop <- log(0.5)
            }
            
          } else
          {
            nexp_prop <- nexp_curr - 1
            if (nexp_prop == 1)
            {
              log_move_curr <- 0
              log_move_prop <- log(0.5)
            } else
            {
              log_move_curr <- log(0.5)
              log_move_prop <- log(0.5)
            }
          }
          
        }
      }
      
      list(
        nexp_prop = nexp_prop,
        log_move_curr = log_move_curr,
        log_move_prop = log_move_prop
      )
      
    }
    
    death_t <-
      function(x,
               u,
               nexp_tcurr_temp,
               nexp_ucurr_temp,
               nexp_prop,
               tau_curr_temp,
               xi_curr_temp,
               ui_curr_temp,
               nseg_curr_temp,
               beta_curr_temp,
               log_move_curr,
               log_move_prop,
               nbasis,
               sigmasqalpha,
               tmin)
      {
        nobs <- dim(x)[1]
        nbeta <- nbasis + 1
        beta_prop <- array(0, c(nbeta, nexp_prop, nexp_ucurr_temp))
        tau_prop <- matrix(1, nexp_prop, nexp_ucurr_temp)
        nseg_prop <- matrix(0, nexp_prop, 1)
        xi_prop <- matrix(0, nexp_prop, 1)
        
        # Drawing cut point to delete
        cut_del <- sample(1:(nexp_tcurr_temp - 1), 1, replace = T)
        j <- 0
        for (k in 1:nexp_prop)
        {
          j <- j + 1
          if (k == cut_del)
          {
            # PROPOSED VALUES
            xi_prop[k] <- xi_curr_temp[j + 1]
            tau_prop[k, ] <- sqrt(tau_curr_temp[j, ] * tau_curr_temp[j + 1, ]) #  Combining 2 taus into 1
            nseg_prop[k] <- nseg_curr_temp[j] + nseg_curr_temp[j + 1] # Combining two segments into 1
            # Evaluating the Likelihood, Proposal and Prior Densities at the Proposed values
            
            loglike_prop <- 0
            log_beta_prop <- 0
            log_beta_prior_prop <- 0
            log_tau_prior_prop <- 0
            
            # Computing mean and variances for beta proposals
            for (i in 1:nexp_ucurr_temp)
            {
              
              postbeta_1 <- postbeta(
                k,
                i,
                nseg_prop[k],
                x,
                u,
                xi_prop,
                ui_curr_temp,
                tau_prop[k, i],
                nbasis,
                sigmasqalpha
              )
              beta_prop[, k, i] <- rmvnorm(1,
                                           postbeta_1$beta_mean,
                                           0.5 * (postbeta_1$beta_var + t(postbeta_1$beta_var))) # Drawing a new value of beta
              # Loglikelihood  at proposed values
              fhat <- postbeta_1$nu_mat %*% beta_prop[, k, i]
              log_prop_spec_dens <- whittle_like(postbeta_1$y, fhat, nseg_prop[k])
              loglike_prop <- loglike_prop + log_prop_spec_dens
              # Evaluating the Prior Densities at the Proposed values for tau and beta
              # Beta
              log_beta_prior_prop <- log_beta_prior_prop + dmvnorm(beta_prop[, k, i], matrix(0, nbeta, 1), diag(c(
                sigmasqalpha, tau_prop[k, i] * matrix(1, nbasis, 1)
              )), log = T)
              # Tau
              log_tau_prior_prop <- log_tau_prior_prop - log(tau_prop[k, i])
              # Evaluating the Proposal Densities at the Proposed values of beta
              # the cut points
              # Beta
              log_beta_prop <- log_beta_prop + dmvnorm(beta_prop[, k, i],
                                                       postbeta_1$beta_mean,
                                                       0.5 * (postbeta_1$beta_var + t(postbeta_1$beta_var)),
                                                       log = T)
            }
            # Segment
            log_seg_prop <- -log(nexp_tcurr_temp - 1)
            # Calculating Jacobian(sum of log Jacobian for each covariate seg)
            log_jacobian <- -sum(log(2 * (
              sqrt(tau_curr_temp[j, ]) + sqrt(tau_curr_temp[j + 1, ])
            ) ^ 2))
            # Calculating log proposal density at proposed values
            log_proposal_prop <- log_beta_prop + log_seg_prop + log_move_prop
            
            
            ## CURRENT VALUES	 (current values have two points j and j+1)
            # Evaluating the Likelihood, Proposal and Prior Densities at the current values
            # Beta proposal and prior
            
            loglike_curr <- 0
            log_beta_curr <- 0
            log_beta_prior_curr <- 0
            log_tau_prior_curr <- 0
            
            for (i in 1:nexp_ucurr_temp)
            {
              for (jj in j:(j + 1))
              {
                
                postbeta_2 <- postbeta(
                  jj,
                  i,
                  nseg_curr_temp[jj],
                  x,
                  u,
                  xi_curr_temp,
                  ui_curr_temp,
                  tau_curr_temp[jj, i],
                  nbasis,
                  sigmasqalpha
                )
                log_beta_curr <- log_beta_curr + dmvnorm(
                  beta_curr_temp[, jj, i],
                  postbeta_2$beta_mean,
                  0.5 * (postbeta_2$beta_var + t(postbeta_2$beta_var)),
                  log = T
                )
                log_beta_prior_curr <- log_beta_prior_curr + dmvnorm(beta_curr_temp[, jj, i],
                                                                     matrix(0, nbeta, 1),
                                                                     diag(c(
                                                                       sigmasqalpha, tau_curr_temp[jj, i] * matrix(1, nbasis, 1)
                                                                     )),
                                                                     log = T)
                log_tau_prior_curr <- log_tau_prior_curr - log(tau_curr_temp[jj, i])
                # Log likelihood at current values
                fhat <- postbeta_2$nu_mat %*% beta_curr_temp[, jj, i]
                log_curr_spec_dens <- whittle_like(postbeta_2$y, fhat, nseg_curr_temp[jj])
                loglike_curr <- loglike_curr + log_curr_spec_dens
                
              }
            }
            # Calculating log proposal density at current values
            log_proposal_curr <- log_move_curr + log_beta_curr
            # Calculating priors at current values
            log_prior_curr <- log_beta_prior_curr + log_tau_prior_curr
            
            j <- j + 1
            
          } else
          {
            xi_prop[k] <- xi_curr_temp[j]
            tau_prop[k, ] <- tau_curr_temp[j, ]
            nseg_prop[k] <- nseg_curr_temp[j]
            beta_prop[, k, ] <- beta_curr_temp[, j, ]
          }
          
        }
        
        # Evaluating target density at proposed values
        log_prior_cut_prop <- 0
        
        if (nexp_prop == 1)
        {
          log_prior_cut_prop <- 0
        } else{
          for (k in 1:(nexp_prop - 1))
          {
            if (k == 1)
            {
              log_prior_cut_prop <- -log(nobs - (nexp_prop - k + 1) * tmin + 1)
            } else
            {
              log_prior_cut_prop <- log_prior_cut_prop - log(nobs - xi_prop[k - 1] - (nexp_prop -
                                                                                        k + 1) * tmin + 1)
            }
          }
        }
        log_target_prop <- loglike_prop + log_tau_prior_prop + log_beta_prior_prop +
          log_prior_cut_prop
        # Evaluating target density at current values
        log_prior_cut_curr <- 0
        for (k in 1:(nexp_tcurr_temp - 1))
        {
          if (k == 1)
          {
            log_prior_cut_curr <- -log(nobs - (nexp_tcurr_temp - k + 1) * tmin + 1)
          } else
          {
            log_prior_cut_curr <- log_prior_cut_curr - log(nobs - xi_curr_temp[k - 1] -
                                                             (nexp_tcurr_temp - k + 1) * tmin + 1)
          }
          
        }
        
        log_target_curr <- loglike_curr + log_prior_curr + log_prior_cut_curr
        
        met_rat <- min(
          1,
          exp(
            log_target_prop - log_target_curr + log_proposal_curr - log_proposal_prop +
              log_jacobian
          )
        )
        
        list(
          met_rat = met_rat,
          nseg_prop = nseg_prop,
          xi_prop = xi_prop,
          tau_prop = tau_prop,
          beta_prop = beta_prop
        )
        
      }
    
    death_u <-
      function(x,
               u,
               nexp_tcurr_temp,
               nexp_ucurr_temp,
               nexp_prop,
               tau_curr_temp,
               xi_curr_temp,
               ui_curr_temp,
               nseg_curr_temp,
               nseg_cov_curr_temp,
               beta_curr_temp,
               log_move_curr,
               log_move_prop,
               nbasis,
               sigmasqalpha,
               umin)
      {
        nbeta <- nbasis + 1
        nu_unique <- length(unique(u))
        beta_prop <- array(0, c(nbeta, nexp_tcurr_temp, nexp_prop))
        tau_prop <- matrix(1, nexp_tcurr_temp, nexp_prop)
        nseg_prop <- matrix(0, nexp_prop, 1)
        ui_prop <- matrix(0, nexp_prop, 1)
        
        # Drawing cut point to delete
        cut_del <- sample(1:(nexp_ucurr_temp - 1), 1, replace = T)
        j <- 0
        for (k in 1:nexp_prop)
        {
          j <- j + 1
          if (k == cut_del)
          {
            # PROPOSED VALUES
            
            ui_prop[k] <- ui_curr_temp[j + 1]
            tau_prop[, k] <- sqrt(tau_curr_temp[, j] * tau_curr_temp[, j + 1]) # Combining 2 taus into 1
            nseg_prop[k] <- nseg_cov_curr_temp[j] + nseg_cov_curr_temp[j + 1] # Combining two segments into 1
            
            # Evaluating the Likelihood, Proposal and Prior Densities at the Proposed values
            
            loglike_prop <- 0
            log_beta_prop <- 0
            log_beta_prior_prop <- 0
            log_tau_prior_prop <- 0
            
            # Computing mean and variances for beta proposals
            for (i in 1:nexp_tcurr_temp)
            {
              #source('~/Desktop/R_code/postbeta.R')
              postbeta_3 <- postbeta(
                i,
                k,
                nseg_curr_temp[i],
                x,
                u,
                xi_curr_temp,
                ui_prop,
                tau_prop[i, k],
                nbasis,
                sigmasqalpha
              )
              beta_prop[, i, k] <- rmvnorm(1,
                                           postbeta_3$beta_mean,
                                           0.5 * (postbeta_3$beta_var + t(postbeta_3$beta_var))) # Drawing a new value of beta
              # Loglikelihood  at proposed values
              fhat <- postbeta_3$nu_mat %*% beta_prop[, i, k]
              log_prop_spec_dens <- whittle_like(postbeta_3$y, fhat, nseg_curr_temp[i])
              loglike_prop <- loglike_prop + log_prop_spec_dens
              
              # Evaluating the Prior Densities at the Proposed values for tau and
              
              # Beta
              log_beta_prior_prop <- log_beta_prior_prop + dmvnorm(beta_prop[, i, k], matrix(0, nbeta, 1), diag(c(
                sigmasqalpha, tau_prop[i, k] * matrix(1, nbasis, 1)
              )), log = T)
              #  Tau
              log_tau_prior_prop <- log_tau_prior_prop - log(tau_prop[i, k])
              
              # Evaluating the Proposal Densities at the Proposed values of beta
              # the cut points
              
              # Beta
              log_beta_prop <- log_beta_prop + dmvnorm(beta_prop[, i, k],
                                                       postbeta_3$beta_mean,
                                                       0.5 * (postbeta_3$beta_var + t(postbeta_3$beta_var)),
                                                       log = T)
              
            }
            
            # Segment
            log_seg_prop <- -log(nexp_ucurr_temp - 1)
            # Calculating Jacobian(sum of log Jacobian for each covariate seg)
            log_jacobian <- -sum(log(2 * (
              sqrt(tau_curr_temp[, j]) + sqrt(tau_curr_temp[, j + 1])
            ) ^ 2))
            # Calculating log proposal density at proposed values
            log_proposal_prop <- log_beta_prop + log_seg_prop + log_move_prop
            
            
            # CURRENT VALUES
            # Evaluating the Likelihood, Proposal and Prior Densities at the Current values
            # Beta proposal and prior
            loglike_curr <- 0
            log_beta_curr <- 0
            log_beta_prior_curr <- 0
            log_tau_prior_curr <- 0
            for (i in 1:nexp_tcurr_temp)
            {
              for (jj in j:(j + 1))
              {
                #source('~/Desktop/R_code/postbeta.R')
                postbeta_4 <- postbeta(
                  i,
                  jj,
                  nseg_curr_temp[i],
                  x,
                  u,
                  xi_curr_temp,
                  ui_curr_temp,
                  tau_curr_temp[i, jj],
                  nbasis,
                  sigmasqalpha
                )
                log_beta_curr <- log_beta_curr + dmvnorm(
                  beta_curr_temp[, i, jj],
                  postbeta_4$beta_mean,
                  0.5 * (postbeta_4$beta_var + t(postbeta_4$beta_var)),
                  log = T
                )
                log_beta_prior_curr <- log_beta_prior_curr + dmvnorm(beta_curr_temp[, i, jj],
                                                                     matrix(0, nbeta, 1),
                                                                     diag(c(
                                                                       sigmasqalpha, tau_curr_temp[i, jj] * matrix(1, nbasis, 1)
                                                                     )),
                                                                     log = T)
                log_tau_prior_curr <- log_tau_prior_curr - log(tau_curr_temp[i, jj])
                # Loglikelihood  at current values
                fhat <- postbeta_4$nu_mat %*% beta_curr_temp[, i, jj]
                log_curr_spec_dens <- whittle_like(postbeta_4$y, fhat, nseg_curr_temp[i])
                loglike_curr <- loglike_curr + log_curr_spec_dens
              }
              
            }
            
            # Calculating log proposal density at current values
            log_proposal_curr <- log_move_curr + log_beta_curr
            # Calculating priors at current values
            log_prior_curr <- log_beta_prior_curr + log_tau_prior_curr
            j = j + 1
          } else
          {
            ui_prop[k] <- ui_curr_temp[j]
            tau_prop[, k] <- tau_curr_temp[, j]
            nseg_prop[k] <- nseg_cov_curr_temp[j]
            beta_prop[, , k] <- beta_curr_temp[, , j]
          }
          
        }
        
        # Evaluating Target density at proposed values(taking duplicates into consideration)
        u_sort <- unique(sort(u))
        log_prior_cut_prop <- 0
        
        if (nexp_prop==1)
        {
          log_prior_cut_prop=0
        }else{
          
          for (k in 1:(nexp_prop - 1))
          {
            if (k == 1)
            {
              log_prior_cut_prop <- -log(nu_unique - (nexp_prop - k + 1) * umin + 1)
            } else
            {
              log_prior_cut_prop <- log_prior_cut_prop - log(nu_unique - sum(u_sort <= ui_prop[k -
                                                                                                 1]) - (nexp_prop - k + 1) * umin + 1)
            }
          }
        }
        
        log_target_prop <- loglike_prop + log_tau_prior_prop + log_beta_prior_prop +
          log_prior_cut_prop
        # Evaluating Target density at current values(taking duplicates into consideration)
        log_prior_cut_curr <- 0
        for (k in 1:(nexp_ucurr_temp - 1))
        {
          if (k == 1)
          {
            log_prior_cut_curr <- -log(nu_unique - (nexp_ucurr_temp - k + 1) * umin +
                                         1)
          } else
          {
            log_prior_cut_curr <- log_prior_cut_curr - log(nu_unique - sum(u_sort <= ui_curr_temp[k -
                                                                                                    1]) - (nexp_ucurr_temp - k + 1) * umin + 1)
          }
        }
        
        log_target_curr <- loglike_curr + log_prior_curr + log_prior_cut_curr
        
        met_rat <- min(
          1,
          exp(
            log_target_prop - log_target_curr + log_proposal_curr - log_proposal_prop +
              log_jacobian
          )
        )
        
        list(
          met_rat = met_rat,
          nseg_prop = nseg_prop,
          ui_prop = ui_prop,
          tau_prop = tau_prop,
          beta_prop = beta_prop
        )
        
        
      }
    
    birth_t <-
      function(x,
               u,
               nexp_tcurr_temp,
               nexp_ucurr_temp,
               nexp_prop,
               tau_curr_temp,
               xi_curr_temp,
               ui_curr_temp,
               nseg_curr_temp,
               beta_curr_temp,
               log_move_curr,
               log_move_prop,
               nbasis,
               sigmasqalpha,
               tmin)
      {
        nobs <- dim(x)[1]
        nbeta <- nbasis + 1
        beta_prop <- array(0, c(nbeta, nexp_prop, nexp_ucurr_temp))
        tau_prop <- matrix(1, nexp_prop, nexp_ucurr_temp)
        nseg_prop <- matrix(0, nexp_prop, 1)
        xi_prop <- matrix(0, nexp_prop, 1)
        
        # Drawing  segment to split
        kk <- which(nseg_curr_temp >= 2 * tmin)  #  Number of segments available for splitting
        nposs_seg <- length(kk)
        seg_cut <- kk[sample(1:nposs_seg, 1, replace = T)] # Drawing segment to split
        nposs_cut <- nseg_curr_temp[seg_cut] - 2 * tmin + 1 # Drawing new cutpoint
        
        for (jj in 1:nexp_tcurr_temp)
        {
          if (jj < seg_cut)
          {
            xi_prop[jj] <- xi_curr_temp[jj]
            tau_prop[jj, ] <- tau_curr_temp[jj, ]
            nseg_prop[jj] <- nseg_curr_temp[jj]
            beta_prop[, jj, ] <- beta_curr_temp[, jj, ]
          } else if (jj == seg_cut)
          {
            index <- sample(1:nposs_cut, 1, replace = T)
            if (seg_cut == 1)
            {
              xi_prop[seg_cut] <- index + tmin - 1
            } else
            {
              xi_prop[seg_cut] <- xi_curr_temp[jj - 1] - 1 + tmin + index
            }
            xi_prop[seg_cut + 1] <- xi_curr_temp[jj]
            zz <- matrix(runif(nexp_ucurr_temp), 1, nexp_ucurr_temp) # Drawing new tausq
            
            tau_prop[seg_cut, ] <- tau_curr_temp[seg_cut, ] * (zz / (1 - zz))
            tau_prop[(seg_cut + 1), ] <- tau_curr_temp[seg_cut, ] * ((1 - zz) /
                                                                       zz)
            nseg_prop[seg_cut] <- index + tmin - 1
            nseg_prop[seg_cut + 1] <- nseg_curr_temp[jj] - nseg_prop[seg_cut]
            
            for (k in jj:(jj + 1))
            {
              for (i in 1:nexp_ucurr_temp)
              {
                
                postbeta_5 <- postbeta(
                  k,
                  i,
                  nseg_prop[k],
                  x,
                  u,
                  xi_prop,
                  ui_curr_temp,
                  tau_prop[k, i],
                  nbasis,
                  sigmasqalpha
                )
                beta_prop[, k, i] <- rmvnorm(1,
                                             postbeta_5$beta_mean,
                                             0.5 * (postbeta_5$beta_var + t(postbeta_5$beta_var))) # Drawing a new value of beta
                
              }
            }
            
          } else
          {
            xi_prop[jj + 1] <- xi_curr_temp[jj]
            tau_prop[jj + 1, ] <- tau_curr_temp[jj, ]
            nseg_prop[jj + 1] <- nseg_curr_temp[jj]
            beta_prop[, jj + 1, ] <- beta_curr_temp[, jj, ]
          }
          
          
        }
        
        # Calculating Jacobian(sum of log Jacobian for each covariate seg)
        log_jacobian <- sum(log(2 * tau_curr_temp[seg_cut, ] / (zz * (1 - zz))))
        
        
        # Evaluating the Likelihood, Proposal and Prior Densities at the Proposed values
        
        loglike_prop <- 0
        log_beta_prop <- 0
        log_beta_prior_prop <- 0
        log_tau_prior_prop <- 0
        
        for (jj in seg_cut:(seg_cut + 1))
        {
          for (i in 1:nexp_ucurr_temp)
          {
            
            postbeta_6 <- postbeta(jj,
                                   i,
                                   nseg_prop[jj],
                                   x,
                                   u,
                                   xi_prop,
                                   ui_curr_temp,
                                   tau_prop[jj, i],
                                   nbasis,
                                   sigmasqalpha)
            log_beta_prop <- log_beta_prop + dmvnorm(beta_prop[, jj, i],
                                                     postbeta_6$beta_mean,
                                                     0.5 * (postbeta_6$beta_var + t(postbeta_6$beta_var)),
                                                     log = T)
            log_beta_prior_prop <- log_beta_prior_prop + dmvnorm(beta_prop[, jj, i], matrix(0, nbeta, 1), diag(c(
              sigmasqalpha, tau_prop[jj, i] * matrix(1, nbasis, 1)
            )), log = T) # Prior Density of beta
            log_tau_prior_prop <- log_tau_prior_prop - log(tau_prop[jj, i]) # Prior density of tausq
            fhat <- postbeta_6$nu_mat %*% beta_prop[, jj, i]
            log_prop_spec_dens <- whittle_like(postbeta_6$y, fhat, nseg_prop[jj]) # Log likelihood  at proposed values
            loglike_prop <- loglike_prop + log_prop_spec_dens
            
          }
          
        }
        
        log_seg_prop <- -log(nposs_seg) # Proposal for segment choice
        log_cut_prop <- -log(nposs_cut) # Proposal for cut point choice
        # Evaluating prior density for cut points at proposed values
        log_prior_cut_prop <- 0
        for (k in 1:(nexp_prop - 1))
        {
          if (k == 1)
          {
            log_prior_cut_prop <- -log(nobs - (nexp_prop - k + 1) * tmin + 1)
          } else
          {
            log_prior_cut_prop <- log_prior_cut_prop - log(nobs - xi_prop[k - 1] - (nexp_prop -
                                                                                      k + 1) * tmin + 1)
          }
        }
        
        # Calculating Log Proposal density at proposed values
        log_proposal_prop <- log_beta_prop + log_seg_prop + log_move_prop + log_cut_prop
        # Calculating Log Prior density at proposed values
        log_prior_prop <- log_beta_prior_prop + log_tau_prior_prop + log_prior_cut_prop
        # Calculating target density at proposed values
        log_target_prop <- loglike_prop + log_prior_prop
        
        # CURRENT VALUES
        # Evaluating the Likelihood, Proposal and Prior Densities at the Current values
        loglike_curr <- 0
        log_beta_curr <- 0
        log_beta_prior_curr <- 0
        log_tau_prior_curr <- 0
        
        for (i in 1:nexp_ucurr_temp)
        {
          # Beta proposal and prior
          #source('~/Desktop/R_code/postbeta.R')
          postbeta_7 <- postbeta(
            seg_cut,
            i,
            nseg_curr_temp[seg_cut],
            x,
            u,
            xi_curr_temp,
            ui_curr_temp,
            tau_curr_temp[seg_cut, i],
            nbasis,
            sigmasqalpha
          )
          log_beta_curr <- log_beta_curr + dmvnorm((beta_curr_temp[, seg_cut, i]),
                                                   (postbeta_7$beta_mean),
                                                   0.5 * (postbeta_7$beta_var + t(postbeta_7$beta_var)),
                                                   log = T
          )
          log_beta_prior_curr <- log_beta_prior_curr + dmvnorm(beta_curr_temp[, seg_cut, i], matrix(0, nbeta, 1), diag(c(
            sigmasqalpha, tau_curr_temp[seg_cut, i] * matrix(1, nbasis, 1)
          )), log = T)
          log_tau_prior_curr <- log_tau_prior_curr - log(tau_curr_temp[seg_cut, i])
          # Log likelihood  at current values
          fhat <- postbeta_7$nu_mat %*% beta_curr_temp[, seg_cut, i]
          log_curr_spec_dens <- whittle_like(postbeta_7$y, fhat, nseg_curr_temp[seg_cut])
          loglike_curr <- loglike_curr + log_curr_spec_dens
          
        }
        
        # Calculating log proposal density at current values
        log_proposal_curr <- log_beta_curr + log_move_curr
        # Evaluating  prior density for cut points at current values
        log_prior_cut_curr <- 0
        
        for (k in 1:(nexp_tcurr_temp - 1))
        {
          if (k == 1)
          {
            log_prior_cut_curr <- -log(nobs - (nexp_tcurr_temp - k + 1) * tmin + 1)
          } else
          {
            log_prior_cut_curr <- log_prior_cut_curr - log(nobs - xi_curr_temp[k - 1] -
                                                             (nexp_tcurr_temp - k + 1) * tmin + 1)
          }
        }
        
        # Calculating priors at current values
        log_prior_curr <- log_beta_prior_curr + log_tau_prior_curr + log_prior_cut_curr
        # Evaluating target densities at current values
        log_target_curr <- loglike_curr + log_prior_curr
        
        met_rat <- min(
          1,
          exp(
            log_target_prop - log_target_curr + log_proposal_curr - log_proposal_prop +
              log_jacobian
          )
        )
        
        list(
          met_rat = met_rat,
          nseg_prop = nseg_prop,
          xi_prop = xi_prop,
          tau_prop = tau_prop,
          beta_prop = beta_prop
        )
        
      }
    
    birth_u <-
      function(x,
               u,
               nexp_tcurr_temp,
               nexp_ucurr_temp,
               nexp_prop,
               tau_curr_temp,
               xi_curr_temp,
               ui_curr_temp,
               nseg_curr_temp,
               nseg_cov_curr_temp,
               beta_curr_temp,
               log_move_curr,
               log_move_prop,
               nbasis,
               sigmasqalpha,
               umin)
      {
        nbeta <- nbasis + 1
        nu_unique <- length(unique(u))
        beta_prop <- array(0, c(nbeta, nexp_tcurr_temp, nexp_prop))
        tau_prop <- matrix(1, nexp_tcurr_temp, nexp_prop)
        nseg_prop <- matrix(0, nexp_prop, 1)
        ui_prop <- matrix(0, nexp_prop, 1)
        
        # Drawing  segment to split
        # Available segments (taking duplicate covariate values into consideration)
        u_sort <- unique(sort(u))
        kk <- matrix(0, nexp_ucurr_temp, 1)
        for (i in 1:nexp_ucurr_temp)
        {
          if ((i == 1) && (sum(u_sort <= ui_curr_temp[i]) >= 2 * umin))
          {
            kk[i] <- i
          } else if ((i > 1) &&
                     (sum(u_sort <= ui_curr_temp[i] &
                          u_sort > ui_curr_temp[i - 1]) >= 2 * umin))
          {
            kk[i] <- i
          }
        }
        
        kk <- as.matrix(kk[kk != 0])
        
        nposs_seg <- length(kk)
        seg_cut <- kk[sample(1:nposs_seg, 1, replace = T)] # Drawing segment to split
        
        # Drawing new cutpoint
        if (seg_cut == 1)
        {
          nposs_cut <- sum(u_sort <= ui_curr_temp[seg_cut]) - 2 * umin + 1
        } else
        {
          nposs_cut <- sum(u_sort <= ui_curr_temp[seg_cut] &
                             u_sort > ui_curr_temp[seg_cut - 1]) - 2 * umin + 1
        }
        
        
        for (jj in 1:nexp_ucurr_temp)
        {
          if (jj < seg_cut)
          {
            ui_prop[jj] <- ui_curr_temp[jj]
            tau_prop[, jj] <- tau_curr_temp[, jj]
            nseg_prop[jj] <- nseg_cov_curr_temp[jj]
            beta_prop[, , jj] <- beta_curr_temp[, , jj]
          } else if (jj == seg_cut)
          {
            index <- sample(1:nposs_cut, 1, replace = T)
            
            if (seg_cut == 1)
            {
              ui_prop[seg_cut] <- u_sort[index + umin - 1]
            } else
            {
              ui_prop[seg_cut] <- u_sort[sum(u_sort <= ui_curr_temp[seg_cut - 1]) - 1 +
                                           umin + index]
            }
            
            ui_prop[seg_cut + 1] <- ui_curr_temp[jj]
            zz <- matrix(runif(nexp_tcurr_temp), nexp_tcurr_temp, 1) # Drawing new tausq
            tau_prop[, seg_cut] <- tau_curr_temp[, seg_cut] * (zz / (1 - zz))
            tau_prop[, seg_cut + 1] <- tau_curr_temp[, seg_cut] * ((1 - zz) /
                                                                     zz)
            
            if (seg_cut == 1)
            {
              nseg_prop[seg_cut] <- sum(u_sort <= ui_prop[seg_cut])
            } else
            {
              nseg_prop[seg_cut] <- sum(u_sort <= ui_prop[seg_cut] &
                                          u_sort > ui_prop[seg_cut - 1])
            }
            
            nseg_prop[seg_cut + 1] <- sum(u_sort <= ui_prop[seg_cut + 1] &
                                            u_sort > ui_prop[seg_cut])
            
            for (k in jj:(jj + 1))
            {
              for (i in 1:nexp_tcurr_temp)
              {
                postbeta_8 <- postbeta(
                  i,
                  k,
                  nseg_curr_temp[i],
                  x,
                  u,
                  xi_curr_temp,
                  ui_prop,
                  tau_prop[i, k],
                  nbasis,
                  sigmasqalpha
                )
                beta_prop[, i, k] <- rmvnorm(1,
                                             postbeta_8$beta_mean,
                                             0.5 * (postbeta_8$beta_var + t(postbeta_8$beta_var))) # Drawing a new value of beta
              }
            }
            
          } else
          {
            ui_prop[jj + 1] <- ui_curr_temp[jj]
            tau_prop[, jj + 1] <- tau_curr_temp[, jj]
            nseg_prop[jj + 1] <- nseg_cov_curr_temp[jj]
            beta_prop[, , jj + 1] <- beta_curr_temp[, , jj]
          }
          
        }
        
        # Calculating Jacobian(sum of log Jacobian for each covariate seg)
        log_jacobian <- sum(log(2 * tau_curr_temp[, seg_cut] / (zz * (1 - zz))))
        # Evaluating the Likelihood, Proposal and Prior Densities at the Proposed values
        loglike_prop <- 0
        log_beta_prop <- 0
        log_beta_prior_prop <- 0
        log_tau_prior_prop <- 0
        
        for (jj in seg_cut:(seg_cut + 1))
        {
          for (i in 1:nexp_tcurr_temp)
          {
            postbeta_9 <- postbeta(
              i,
              jj,
              nseg_curr_temp[i],
              x,
              u,
              xi_curr_temp,
              ui_prop,
              tau_prop[i, jj],
              nbasis,
              sigmasqalpha
            )
            log_beta_prop <- log_beta_prop + dmvnorm(beta_prop[, i, jj],
                                                     postbeta_9$beta_mean,
                                                     0.5 * (postbeta_9$beta_var + t(postbeta_9$beta_var)),
                                                     log = T)
            log_beta_prior_prop <- log_beta_prior_prop + dmvnorm(beta_prop[, i, jj], matrix(0, nbeta, 1), diag(c(
              sigmasqalpha, tau_prop[i, jj] * matrix(1, nbasis, 1)
            )), log = T) # Prior Density of beta
            log_tau_prior_prop <- log_tau_prior_prop - log(tau_prop[i, jj]) # Prior Density of Tausq
            fhat <- postbeta_9$nu_mat %*% beta_prop[, i, jj]
            log_prop_spec_dens <- whittle_like(postbeta_9$y, fhat, nseg_curr_temp[i]) # Loglikelihood  at proposed values
            loglike_prop <- loglike_prop + log_prop_spec_dens
          }
          
        }
        
        
        log_seg_prop <- -log(nposs_seg) # Proposal for segment choice
        log_cut_prop <- -log(nposs_cut) # Proposal for cut point choice
        # Evaluating prior density for cut points at proposed values
        log_prior_cut_prop = 0
        
        for (k in 1:(nexp_prop - 1))
        {
          if (k == 1)
          {
            log_prior_cut_prop <- -log(nu_unique - (nexp_prop - k + 1) * umin + 1)
          } else
          {
            log_prior_cut_prop <- log_prior_cut_prop - log(nu_unique - sum(u_sort <= ui_prop[k -
                                                                                               1]) - (nexp_prop - k + 1) * umin + 1)
          }
        }
        
        
        # Calculating log proposal density at proposed values
        log_proposal_prop <- log_beta_prop + log_seg_prop + log_move_prop + log_cut_prop
        # Calculating log prior density at proposed values
        log_prior_prop <- log_beta_prior_prop + log_tau_prior_prop + log_prior_cut_prop
        # Calculating target density at proposed values
        log_target_prop <- loglike_prop + log_prior_prop
        
        # CURRENT VALUES
        # Evaluating the Likelihood, Proposal and Prior Densities at the Current values
        
        loglike_curr <- 0
        log_beta_curr <- 0
        log_beta_prior_curr <- 0
        log_tau_prior_curr <- 0
        
        for (i in 1:nexp_tcurr_temp)
        {
          # Beta proposal and prior
          postbeta_10 <- postbeta(
            i,
            seg_cut,
            nseg_curr_temp[i],
            x,
            u,
            xi_curr_temp,
            ui_curr_temp,
            tau_curr_temp[i, seg_cut],
            nbasis,
            sigmasqalpha
          )
          log_beta_curr <- log_beta_curr + dmvnorm(beta_curr_temp[, i, seg_cut],
                                                   postbeta_10$beta_mean,
                                                   0.5 * (postbeta_10$beta_var + t(postbeta_10$beta_var)),
                                                   log = T)
          log_beta_prior_curr <- log_beta_prior_curr + dmvnorm(beta_curr_temp[, i, seg_cut], matrix(0, nbeta, 1), diag(c(
            sigmasqalpha, tau_curr_temp[i, seg_cut] * matrix(1, nbasis, 1)
          )), log = T)
          log_tau_prior_curr <- log_tau_prior_curr - log(tau_curr_temp[i, seg_cut])
          # Log likelihood  at current values
          fhat <- postbeta_10$nu_mat %*% beta_curr_temp[, i, seg_cut]
          log_curr_spec_dens <- whittle_like(postbeta_10$y, fhat, nseg_curr_temp[i])
          loglike_curr <- loglike_curr + log_curr_spec_dens
          
        }
        
        # Calculating log proposal density at current values
        log_proposal_curr <- log_beta_curr + log_move_curr
        # Evaluating  prior density for cut points at current values
        log_prior_cut_curr <- 0
        for (k in 1:(nexp_ucurr_temp - 1))
        {
          if (k == 1)
          {
            log_prior_cut_curr <- -log(nu_unique - (nexp_ucurr_temp - k + 1) * umin +
                                         1)
          } else
          {
            log_prior_cut_curr <- log_prior_cut_curr - log(nu_unique - sum(u_sort <= ui_curr_temp[k -
                                                                                                    1]) - (nexp_ucurr_temp - k + 1) * umin + 1)
          }
        }
        
        
        # Calculating priors at current values
        log_prior_curr <- log_beta_prior_curr + log_tau_prior_curr + log_prior_cut_curr
        # Evaluating target densities at current values
        log_target_curr <- loglike_curr + log_prior_curr
        
        met_rat <- min(
          1,
          exp(
            log_target_prop - log_target_curr + log_proposal_curr - log_proposal_prop +
              log_jacobian
          )
        )
        
        list(
          met_rat = met_rat,
          nseg_prop = nseg_prop,
          ui_prop = ui_prop,
          tau_prop = tau_prop,
          beta_prop = beta_prop
        )
        
      }
    
    within_t <- function(x,
                         u,
                         nexp_tcurr_temp,
                         nexp_ucurr_temp,
                         xi_curr_temp,
                         ui_curr_temp,
                         beta_curr_temp,
                         nseg_curr_temp,
                         tau_temp,
                         nbasis,
                         sigmasqalpha,
                         tmin,
                         prob_mm1)
    {
      nobs <- dim(x)[1]
      nbeta <- nbasis + 1
      xi_prop <- xi_curr_temp
      beta_prop <- beta_curr_temp
      nseg_new <- nseg_curr_temp
      
      if (nexp_tcurr_temp > 1)
      {
        seg_temp <- sample(1:(nexp_tcurr_temp - 1), 1, replace = T)  # Drawing Segment to cut
        r <- runif(1)
        cut_poss_curr <- xi_curr_temp[seg_temp]
        nposs_prior <- nseg_curr_temp[seg_temp] + nseg_curr_temp[seg_temp +
                                                                   1] - 2 * tmin + 1
        
        if (r < prob_mm1)
        {
          # small move
          size_w <- 1
          if (nseg_curr_temp[seg_temp] == tmin &&
              nseg_curr_temp[seg_temp + 1] == tmin)
          {
            nposs <- 1  # Number of possible locations for new cutpoint
            new_index <- sample(1:nposs, 1, replace = T)  # Drawing index of new cutpoint
            cut_poss_new <- xi_curr_temp[seg_temp] - 1 + new_index
          } else if (nseg_curr_temp[seg_temp] == tmin)
          {
            nposs <- 2  # Number of possible locations for new cutpoint
            new_index <- sample(1:nposs, 1, replace = T)  # Drawing index of new cutpoint
            cut_poss_new <- xi_curr_temp[seg_temp] - 1 + new_index
          } else if (nseg_curr_temp[seg_temp + 1] == tmin)
          {
            nposs <- 2 # Number of possible locations for new cutpoint
            new_index <- sample(1:nposs, 1, replace = T) # Drawing index of new cutpoint
            cut_poss_new <- xi_curr_temp[seg_temp] + 1 - new_index
          } else
          {
            nposs <- 3 #  Number of possible locations for new cutpoint
            new_index <- sample(1:nposs, 1, replace = T)  # Drawing index of new cutpoint
            cut_poss_new <- xi_curr_temp[seg_temp] - 2 + new_index
          }
          
        } else
          # big move
        {
          size_w <- 2
          new_index <- sample(1:nposs_prior, 1, replace = T)
          # cut_poss_new=sum(nseg_curr_temp[1:(seg_temp-1)])-1+opt$tmin+new_index
          if (seg_temp > 1) {
            cut_poss_new <-
              sum(nseg_curr_temp[1:(seg_temp - 1)]) - 1 + tmin + new_index
          }
          else {
            cut_poss_new <- -1 + tmin + new_index
          }
          
        }
        
        xi_prop[seg_temp] <- cut_poss_new
        if (seg_temp > 1)
        {
          nseg_new[seg_temp] <- xi_prop[seg_temp] - xi_curr_temp[seg_temp - 1]  # Number of observations in lower part of new cutpoint
          
        } else
        {
          nseg_new[seg_temp] <- xi_prop[seg_temp]
        }
        
        nseg_new[seg_temp + 1] <- nseg_curr_temp[seg_temp] + nseg_curr_temp[seg_temp +
                                                                              1] - nseg_new[seg_temp]  #Number of observations in upper part of new cutpoint
        
        # Evaluating the Proposal density for the cut-point at the current and proposed values
        if (abs(cut_poss_new - cut_poss_curr) > 1)
        {
          log_prop_cut_prop <- log(1 - prob_mm1) - log(nposs_prior)
          log_prop_cut_curr <- log(1 - prob_mm1) - log(nposs_prior)
        } else if (nseg_curr_temp[seg_temp] == tmin &&
                   nseg_curr_temp[seg_temp + 1] == tmin)
        {
          log_prop_cut_prop <- 0
          log_prop_cut_curr <- 0
        } else
        {
          if (nseg_curr_temp[seg_temp] == tmin ||
              nseg_curr_temp[seg_temp + 1] == tmin)
          {
            log_prop_cut_prop <- log(1 - prob_mm1) - log(nposs_prior) + log(1 / 2) +
              log(prob_mm1)
          } else
          {
            log_prop_cut_prop <- log(1 - prob_mm1) - log(nposs_prior) + log(1 / 3) +
              log(prob_mm1)
          }
          if (nseg_new[seg_temp] == tmin ||
              nseg_new[seg_temp + 1] == tmin)
          {
            log_prop_cut_curr <- log(1 - prob_mm1) - log(nposs_prior) + log(1 / 2) +
              log(prob_mm1)
          } else
          {
            log_prop_cut_curr <- log(1 - prob_mm1) - log(nposs_prior) + log(1 / 3) +
              log(prob_mm1)
          }
        }
        
        # Evaluating the Loglikelihood, Priors and Proposals at the current values
        loglike_curr <- 0
        log_beta_curr_temp <- 0
        log_prior_curr <- 0
        
        for (i in 1:nexp_ucurr_temp)
        {
          for (j in seg_temp:(seg_temp + 1))
          {
            
            postbeta_11 <- postbeta(
              j,
              i,
              nseg_curr_temp[j],
              x,
              u,
              xi_curr_temp,
              ui_curr_temp,
              tau_temp[j, i],
              nbasis,
              sigmasqalpha
            )
            # Compute log proposal density of beta at current  values
            log_beta_curr_temp <- log_beta_curr_temp + dmvnorm(beta_curr_temp[, j, i],
                                                               postbeta_11$beta_mean,
                                                               0.5 * (postbeta_11$beta_var + t(postbeta_11$beta_var)),
                                                               log = TRUE)
            fhat = postbeta_11$nu_mat %*% beta_curr_temp[, j, i]
            # Compute Loglike at current values
            log_curr_spec_dens <- whittle_like(postbeta_11$y, fhat, nseg_curr_temp[j])
            loglike_curr <- loglike_curr + log_curr_spec_dens
            # Compute priors at current values
            log_prior_curr <- log_prior_curr + dmvnorm(beta_curr_temp[, j, i], matrix(0, nbeta, 1), diag(c(
              sigmasqalpha, tau_temp[j, i] * matrix(1, nbasis, 1)
            )), log = T)
            
          }
        }
        
        # Evaluating the Loglikelihood, Priors and Proposals at the proposed values
        # Likelihood
        loglike_prop <- 0
        log_beta_prop <- 0
        log_prior_prop <- 0
        
        for (i in 1:nexp_ucurr_temp)
        {
          for (j in seg_temp:(seg_temp + 1))
          {
            
            postbeta_12 <- postbeta(j,
                                    i,
                                    nseg_new[j],
                                    x,
                                    u,
                                    xi_prop,
                                    ui_curr_temp,
                                    tau_temp[j, i],
                                    nbasis,
                                    sigmasqalpha)
            beta_prop[, j, i] <- rmvnorm(1, postbeta_12$beta_mean, 0.5 * (postbeta_12$beta_var + t(postbeta_12$beta_var)))
            # Compute log proposal density of beta at proposed  values
            log_beta_prop <- log_beta_prop + dmvnorm(beta_prop[, j, i],
                                                     postbeta_12$beta_mean,
                                                     0.5 * (postbeta_12$beta_var + t(postbeta_12$beta_var)),
                                                     log = T)
            fhat <- postbeta_12$nu_mat %*% beta_prop[, j, i]
            # Compute Loglike at proposed values
            log_prop_spec_dens <- whittle_like(postbeta_12$y, fhat, nseg_new[j])
            loglike_prop <- loglike_prop + log_prop_spec_dens
            # Compute priors at proposed values
            log_prior_prop <- log_prior_prop + dmvnorm(beta_prop[, j, i], matrix(0, nbeta, 1), diag(c(
              sigmasqalpha, tau_temp[j, i] * matrix(1, nbasis, 1)
            )), log = T)
            
          }
        }
        
        # Proposal for beta
        log_proposal_curr <- log_beta_curr_temp + log_prop_cut_curr
        log_proposal_prop <- log_beta_prop + log_prop_cut_prop
        log_prior_cut_prop <- 0
        log_prior_cut_curr <- 0
        for (k in 1:(nexp_tcurr_temp - 1))
        {
          if (k == 1)
          {
            log_prior_cut_prop <- -log(nobs - (nexp_tcurr_temp - k + 1) * tmin + 1)
            log_prior_cut_curr <- -log(nobs - (nexp_tcurr_temp - k + 1) * tmin +
                                         1)
          } else
          {
            log_prior_cut_prop <- log_prior_cut_prop - log(nobs - xi_prop[k - 1] - (nexp_tcurr_temp -
                                                                                      k + 1) * tmin + 1)
            log_prior_cut_curr <- log_prior_cut_curr - log(nobs - xi_curr_temp[k -
                                                                                 1] - (nexp_tcurr_temp - k + 1) * tmin + 1)
          }
        }
        
        log_target_prop <- loglike_prop + log_prior_prop + log_prior_cut_prop
        log_target_curr <- loglike_curr + log_prior_curr + log_prior_cut_curr
        
      } else
      {
        # no move since only one segment
        size_w <- 0
        nseg_new <- nobs
        log_beta_prop <- 0
        log_beta_curr_temp <- 0
        loglike_prop <- 0
        loglike_curr <- 0
        log_prior_prop <- 0
        log_prior_curr <- 0
        
        for (i in 1:nexp_ucurr_temp)
        {
          
          postbeta_13 <- postbeta(1,
                                  i,
                                  nobs,
                                  x,
                                  u,
                                  xi_prop,
                                  ui_curr_temp,
                                  tau_temp[1, i],
                                  nbasis,
                                  sigmasqalpha)
          beta_prop[, 1, i] <- rmvnorm(1, postbeta_13$beta_mean, 0.5 * (postbeta_13$beta_var + t(postbeta_13$beta_var)))
          # Compute log proposal density of beta at proposed  values
          log_beta_prop <- log_beta_prop + dmvnorm(beta_prop[, 1, i],
                                                   postbeta_13$beta_mean,
                                                   0.5 * (postbeta_13$beta_var + t(postbeta_13$beta_var)),
                                                   log = T)
          # Compute log proposal density of beta at current  values
          log_beta_curr_temp <- log_beta_curr_temp + dmvnorm(
            beta_curr_temp[, 1, i],
            postbeta_13$beta_mean,
            0.5 * (postbeta_13$beta_var + t(postbeta_13$beta_var)),
            log = T
          )
          # Compute Loglike at proposed values
          fhat <- postbeta_13$nu_mat %*% beta_prop[, 1, i]
          log_prop_spec_dens <- whittle_like(postbeta_13$y, fhat, nobs)
          loglike_prop <- loglike_prop + log_prop_spec_dens
          # Compute Loglike at proposed values
          fhat <- postbeta_13$nu_mat %*% beta_curr_temp[, 1, i]
          log_curr_spec_dens <- whittle_like(postbeta_13$y, fhat, nobs)
          loglike_curr <- loglike_curr + log_curr_spec_dens
          # Compute Priors at proposed values
          log_prior_prop <- log_prior_prop + dmvnorm(beta_prop[, 1, i], matrix(0, nbeta, 1), diag(c(
            sigmasqalpha, tau_temp[1, i] * matrix(1, nbasis, 1)
          )), log = T)
          # Compute Priors at current values
          log_prior_curr <- log_prior_curr + dmvnorm(beta_curr_temp[, 1, i], matrix(0, nbeta, 1), diag(c(
            sigmasqalpha, tau_temp[1, i] * matrix(1, nbasis, 1)
          )), log = T)
          
        }
        
        log_proposal_curr <- log_beta_curr_temp
        log_proposal_prop <- log_beta_prop
        log_target_prop <- loglike_prop + log_prior_prop
        log_target_curr <- loglike_curr + log_prior_curr
        
      }
      
      epsilon <- min(1,
                     exp(
                       log_target_prop - log_target_curr + log_proposal_curr - log_proposal_prop
                     ))
      
      list(
        epsilon = epsilon,
        xi_prop = xi_prop,
        beta_prop = beta_prop,
        nseg_new = nseg_new,
        size_w = size_w
      )
    }
    
    within_u <-
      function(x,
               u,
               nexp_tcurr_temp,
               nexp_ucurr_temp,
               xi_curr_temp,
               ui_curr_temp,
               beta_curr_temp,
               nseg_curr_temp,
               nseg_cov_curr_temp,
               tau_temp,
               nbasis,
               sigmasqalpha,
               umin,
               prob_mm1_u)
      {
        nbeta <- nbasis + 1
        nrep <- dim(x)[2]
        nu_unique <- length(unique(u))
        
        ui_prop <- ui_curr_temp
        beta_prop <- beta_curr_temp
        nseg_new <- nseg_cov_curr_temp
        
        u_sort <- unique(sort(u))
        
        if (nexp_ucurr_temp > 1)
        {
          seg_temp <- sample(1:(nexp_ucurr_temp - 1), 1, replace = T)  # Drawing Segment to cut
          r <- runif(1)
          cut_poss_curr <- sum(u_sort <= ui_curr_temp[seg_temp]) # index of current cut point
          if (seg_temp == 1)
          {
            nposs_prior <- sum(u_sort <= ui_curr_temp[seg_temp + 1]) - 2 * umin + 1
          } else
          {
            nposs_prior <- sum(u_sort > ui_curr_temp[seg_temp - 1] &
                                 u_sort <= ui_curr_temp[seg_temp + 1]) - 2 * umin + 1
          }
          if (seg_temp == 1)
          {
            nseg_low <- sum(u_sort <= ui_curr_temp[seg_temp])
            
          } else
          {
            nseg_low <- sum(u_sort > ui_curr_temp[seg_temp - 1] &
                              u_sort <= ui_curr_temp[seg_temp])
          }
          
          nseg_high <- sum(u_sort > ui_curr_temp[seg_temp] &
                             u_sort <= ui_curr_temp[seg_temp + 1])
          
          if (r < prob_mm1_u)
          {
            # small move
            size_w <- 1
            
            if (nseg_low == umin && nseg_high == umin)
            {
              nposs <- 1 # Number of possible locations for new cutpoint
              new_index <- sample(1:nposs, 1, replace = T) # Drawing index of new cutpoint
              cut_poss_new <- cut_poss_curr - 1 + new_index
            } else if (nseg_low == umin)
            {
              nposs <- 2 # Number of possible locations for new cutpoint
              new_index <- sample(1:nposs, 1, replace = T) # Drawing index of new cutpoint
              cut_poss_new <- cut_poss_curr - 1 + new_index
            } else if (nseg_high == umin)
            {
              nposs <- 2 # Number of possible locations for new cutpoint
              new_index <- sample(1:nposs, 1, replace = T) # Drawing index of new cutpoint
              cut_poss_new <- cut_poss_curr + 1 - new_index
            } else
            {
              nposs <- 3 # Number of possible locations for new cutpoint
              new_index <- sample(1:nposs, 1, replace = T) # Drawing index of new cutpoint
              cut_poss_new <- cut_poss_curr - 2 + new_index
            }
          } else
          {
            # big move
            size_w <- 2
            new_index <- sample(1:nposs_prior, 1, replace = T)
            if (seg_temp == 1)
            {
              cut_poss_new <- 0 - 1 + umin + new_index
            } else
            {
              cut_poss_new <- sum(u_sort <= ui_curr_temp[seg_temp - 1]) - 1 + umin + new_index
            }
            
          }
          
          ui_prop[seg_temp] <- u_sort[cut_poss_new]
          
          if (seg_temp > 1)
          {
            nseg_new[seg_temp] <- sum(u <= ui_prop[seg_temp] &
                                        u > ui_prop[seg_temp - 1]) # Number of observations in lower part of new cutpoint
            
          } else
          {
            nseg_new[seg_temp] <- sum(u <= ui_prop[seg_temp])
          }
          
          nseg_new[seg_temp + 1] <- sum(u <= ui_prop[seg_temp + 1] &
                                          u > ui_prop[seg_temp]) # Number of observations in upper part of new cutpoint
          
          # Evaluating the Proposal density for the cut-point at the current and proposed values
          if (seg_temp == 1)
          {
            nseg_new_low <- sum(u_sort <= ui_prop[seg_temp])
          } else
          {
            nseg_new_low <- sum(u_sort > ui_prop[seg_temp - 1] &
                                  u_sort <= ui_prop[seg_temp])
          }
          
          nseg_new_high <- sum(u_sort > ui_prop[seg_temp] &
                                 u_sort <= ui_prop[seg_temp + 1])
          
          if (abs(cut_poss_new - cut_poss_curr) > 1)
          {
            log_prop_cut_prop <- log(1 - prob_mm1_u) - log(nposs_prior)
            log_prop_cut_curr <- log(1 - prob_mm1_u) - log(nposs_prior)
          } else if (nseg_low == umin && nseg_high == umin)
          {
            log_prop_cut_prop <- 0
            log_prop_cut_curr <- 0
          } else
          {
            if (nseg_low == umin || nseg_high == umin)
            {
              log_prop_cut_prop <- log(1 - prob_mm1_u) - log(nposs_prior) + log(1 / 2) +
                log(prob_mm1_u)
            } else
            {
              log_prop_cut_prop <- log(1 - prob_mm1_u) - log(nposs_prior) + log(1 / 3) +
                log(prob_mm1_u)
            }
            if (nseg_new_low == umin || nseg_new_high == umin)
            {
              log_prop_cut_curr <- log(1 - prob_mm1_u) - log(nposs_prior) + log(1 / 2) +
                log(prob_mm1_u)
            } else
            {
              log_prop_cut_curr <- log(1 - prob_mm1_u) - log(nposs_prior) + log(1 / 3) +
                log(prob_mm1_u)
            }
          }
          
          # Evaluating the Loglikelihood, Priors and Proposals at the current values
          loglike_curr <- 0
          log_beta_curr_temp <- 0
          log_prior_curr <- 0
          
          for (i in seg_temp:(seg_temp + 1))
          {
            for (j in 1:nexp_tcurr_temp)
            {
              
              postbeta_14 <- postbeta(
                j,
                i,
                nseg_curr_temp[j],
                x,
                u,
                xi_curr_temp,
                ui_curr_temp,
                tau_temp[j, i],
                nbasis,
                sigmasqalpha
              )
              # Compute log proposal density of beta at current  values
              log_beta_curr_temp <- log_beta_curr_temp + dmvnorm(
                beta_curr_temp[, j, i],
                postbeta_14$beta_mean,
                0.5 * (postbeta_14$beta_var + t(postbeta_14$beta_var)),
                log = T
              )
              fhat <- postbeta_14$nu_mat %*% beta_curr_temp[, j, i]
              # Compute Loglike at current values
              log_curr_spec_dens <- whittle_like(postbeta_14$y, fhat, nseg_curr_temp[j])
              loglike_curr <- loglike_curr + log_curr_spec_dens
              # Compute priors at current values
              log_prior_curr <- log_prior_curr + dmvnorm(beta_curr_temp[, j, i], matrix(0, nbeta, 1), diag(c(
                sigmasqalpha, tau_temp[j, i] * matrix(1, nbasis, 1)
              )), log = T)
              
            }
          }
          
          # Evaluating the Loglikelihood, Priors and Proposals at the proposed values
          # Likelihood
          loglike_prop <- 0
          log_beta_prop <- 0
          log_prior_prop <- 0
          
          for (i in seg_temp:(seg_temp + 1))
          {
            for (j in 1:nexp_tcurr_temp)
            {
              
              postbeta_15 <- postbeta(
                j,
                i,
                nseg_curr_temp[j],
                x,
                u,
                xi_curr_temp,
                ui_prop,
                tau_temp[j, i],
                nbasis,
                sigmasqalpha
              )
              beta_prop[, j, i] <- rmvnorm(1,
                                           postbeta_15$beta_mean,
                                           0.5 * (postbeta_15$beta_var + t(postbeta_15$beta_var)))
              # Compute log proposal density of beta at proposed  values
              log_beta_prop <- log_beta_prop + dmvnorm(beta_prop[, j, i],
                                                       postbeta_15$beta_mean,
                                                       0.5 * (postbeta_15$beta_var + t(postbeta_15$beta_var)),
                                                       log = T)
              fhat <- postbeta_15$nu_mat %*% beta_prop[, j, i]
              # Compute Loglike at proposed values
              log_prop_spec_dens <- whittle_like(postbeta_15$y, fhat, nseg_curr_temp[j])
              loglike_prop <- loglike_prop + log_prop_spec_dens
              # Compute priors at proposed values
              log_prior_prop <- log_prior_prop + dmvnorm(beta_prop[, j, i], matrix(0, nbeta, 1), diag(c(
                sigmasqalpha, tau_temp[j, i] * matrix(1, nbasis, 1)
              )), log = T)
            }
          }
          
          # Proposal for beta
          log_proposal_curr <- log_beta_curr_temp + log_prop_cut_curr
          log_proposal_prop <- log_beta_prop + log_prop_cut_prop
          log_prior_cut_prop <- 0
          log_prior_cut_curr <- 0
          
          for (k in 1:(nexp_ucurr_temp - 1))
          {
            if (k == 1)
            {
              log_prior_cut_prop <- -log(nu_unique - (nexp_ucurr_temp - k + 1) * umin +
                                           1)
              log_prior_cut_curr <- -log(nu_unique - (nexp_ucurr_temp - k + 1) *
                                           umin + 1)
            } else
            {
              log_prior_cut_prop <- log_prior_cut_prop - log(nu_unique - sum(u_sort <= ui_prop[k -
                                                                                                 1]) - (nexp_ucurr_temp - k + 1) * umin + 1)
              log_prior_cut_curr <- log_prior_cut_curr - log(nu_unique - sum(u_sort <=
                                                                               ui_curr_temp[k - 1]) - (nexp_ucurr_temp - k + 1) * umin + 1)
            }
            
          }
          
          log_target_prop <- loglike_prop + log_prior_prop + log_prior_cut_prop
          log_target_curr <- loglike_curr + log_prior_curr + log_prior_cut_curr
          
        } else
        {
          # no move since only one segment
          size_w <- 0
          nseg_new <- nrep
          log_beta_prop <- 0
          log_beta_curr_temp <- 0
          loglike_prop <- 0
          loglike_curr <- 0
          log_prior_prop <- 0
          log_prior_curr <- 0
          
          for (i in 1:nexp_tcurr_temp)
          {
            
            postbeta_16 <- postbeta(
              i,
              1,
              nseg_curr_temp[i],
              x,
              u,
              xi_curr_temp,
              ui_prop,
              tau_temp[i, 1],
              nbasis,
              sigmasqalpha
            )
            beta_prop[, i, 1] <- rmvnorm(1,
                                         postbeta_16$beta_mean,
                                         0.5 * (postbeta_16$beta_var + t(postbeta_16$beta_var)))
            # Compute log proposal density of beta at proposed  values
            log_beta_prop <- log_beta_prop + dmvnorm(beta_prop[, i, 1],
                                                     postbeta_16$beta_mean,
                                                     0.5 * (postbeta_16$beta_var + t(postbeta_16$beta_var)),
                                                     log = T)
            # Compute Loglike at proposed values
            fhat <- postbeta_16$nu_mat %*% beta_prop[, i, 1]
            log_prop_spec_dens <- whittle_like(postbeta_16$y, fhat, nseg_curr_temp[i])
            loglike_prop <- loglike_prop + log_prop_spec_dens
            # Compute Loglike at proposed values
            fhat <- postbeta_16$nu_mat %*% beta_curr_temp[, i, 1]
            log_curr_spec_dens <- whittle_like(postbeta_16$y, fhat, nseg_curr_temp[i])
            loglike_curr <- loglike_curr + log_curr_spec_dens
            # Compute Priors at proposed values
            log_prior_prop <- log_prior_prop + dmvnorm(beta_prop[, i, 1], matrix(0, nbeta, 1), diag(c(
              sigmasqalpha, tau_temp[i, 1] * matrix(1, nbasis, 1)
            )), log = T)
            # Compute Priors at current values
            log_prior_curr <- log_prior_curr + dmvnorm(beta_curr_temp[, i, 1], matrix(0, nbeta, 1), diag(c(
              sigmasqalpha, tau_temp[i, 1] * matrix(1, nbasis, 1)
            )), log = T)
            
          }
          
          log_proposal_curr <- log_beta_curr_temp
          log_proposal_prop <- log_beta_prop
          log_target_prop <- loglike_prop + log_prior_prop
          log_target_curr <- loglike_curr + log_prior_curr
          
          
        }
        
        epsilon <- min(1,
                       exp(
                         log_target_prop - log_target_curr + log_proposal_curr - log_proposal_prop
                       ))
        
        list(
          epsilon = epsilon,
          ui_prop = ui_prop,
          beta_prop = beta_prop,
          nseg_new = nseg_new,
          size_w = size_w
        )
        
      }
    
    
    
    ## set seed
    set.seed(seed)
    
    ## set up initial values  create matrix
    nobs <- dim(x)[1]
    nrep <- dim(x)[2]
    nu_unique <- length(unique(u))
    nbeta <- nbasis + 1
    u_sort <- sort(unique(u))
    
    
    nexp_tcurr <- rep(1, nloop + 1)
    nexp_ucurr <- rep(1, nloop + 1)
    
    tausq_curr <- matrix(list(), nexp_tmax, nexp_umax)
    beta_curr <- matrix(list(), nexp_tmax, nexp_umax)
    xi_curr <- matrix(list(), nexp_tmax, 1)
    ui_curr <- matrix(list(), nexp_umax, 1)
    nseg_tcurr <- matrix(list(), nexp_tmax, 1)
    nseg_ucurr <- matrix(list(), nexp_umax, 1)
    log_spec_hat_curr <- matrix(list(), nexp_tmax, nexp_umax)
    
    
    
    # nu mat hat (current)
    nfreq_hat <- 50
    freq_hat <- (0:nfreq_hat) / (2 * nfreq_hat)
    
    nu_mat_hat <- lin_basis_func(freq_hat, nbeta)
    
    for (j in 1:nexp_tmax) {
      for (i in 1:nexp_umax) {
        tausq_curr[[j, i]] <-
          array(rep(1, i * j * (nloop + 1)), dim = c(j, i, nloop + 1))
        beta_curr[[j, i]] <-
          array(rep(1, nbeta * i * j * (nloop + 1)), dim = c(nbeta, j,
                                                             i, nloop + 1))
        log_spec_hat_curr[[j, i]] <-
          array(rep(0, (nfreq_hat + 1) * j * i *
                      (nloop + 1)), dim = c(nfreq_hat + 1, j, i, nloop + 1))
        
      }
    }
    
    for (j in 1:nexp_tmax) {
      xi_curr[[j]] <- matrix(1, j, nloop + 1)
      nseg_tcurr[[j]] <- matrix(1, j, nloop + 1)
    }
    
    for (i in 1:nexp_umax) {
      ui_curr[[i]] <- matrix(1, i, nloop + 1)
      
      nseg_ucurr[[i]] <- matrix(1, i, nloop + 1)
    }
    
    
    
    ## number of segments for first iteration
    if (Rev_Jump_t == 1)
    {
      nexp_tcurr[1] <- sample(1:nexp_tmax, 1, replace = T)
    } else {
      nexp_tcurr[1] <- nexp_tmax
    }
    
    if (Rev_Jump_u == 1)
    {
      nexp_ucurr[1] <- sample(1:nexp_umax, 1, replace = T)
    } else {
      nexp_ucurr[1] <- nexp_umax
    }
    
    # tau sq  is uniform U(0,1000) distribution for first iteration
    tau_up_limit <- 1000
    
    for (j in 1:nexp_tcurr[1]) {
      for (i in 1:nexp_ucurr[1]) {
        tausq_curr[[nexp_tcurr[1], nexp_ucurr[1]]][j, i, 1] <-
          runif(1) * tau_up_limit
      }
    }
    
    
    
    ## partitions
    
    # time t partition
    
    for (j in 1:nexp_tcurr[1])
    {
      if (nexp_tcurr[1] == 1)
      {
        xi_curr[[nexp_tcurr[1]]][j, 1] <- nobs
        nseg_tcurr[[nexp_tcurr[1]]][j, 1] <- nobs
      } else
      {
        if (j == 1)
        {
          nposs <- nobs - nexp_tcurr[1] * tmin + 1
          xi_curr[[nexp_tcurr[1]]][j, 1] <- tmin + sample(1:nposs, 1, replace = T) -
            1
          nseg_tcurr[[nexp_tcurr[1]]][j, 1] <- xi_curr[[nexp_tcurr[1]]][j, 1]
        }
        else if (j > 1 && j < nexp_tcurr[1])
        {
          nposs <- nobs - xi_curr[[nexp_tcurr[1]]][j - 1, 1] - tmin * (nexp_tcurr[1] -
                                                                         j + 1) + 1
          xi_curr[[nexp_tcurr[1]]][j, 1] <- tmin + sample(1:nposs, 1, replace = T) +
            xi_curr[[nexp_tcurr[1]]][j - 1, 1] - 1
          nseg_tcurr[[nexp_tcurr[1]]][j, 1] <- xi_curr[[nexp_tcurr[1]]][j, 1] -
            xi_curr[[nexp_tcurr[1]]][j - 1, 1]
        } else
          # i==nexp_tcurr
        {
          xi_curr[[nexp_tcurr[1]]][j, 1] <- nobs
          nseg_tcurr[[nexp_tcurr[1]]][j, 1] <- xi_curr[[nexp_tcurr[1]]][j, 1] -
            xi_curr[[nexp_tcurr[1]]][j - 1, 1]
        }
      }
    }
    
    # covariate u partition
    
    for (i in 1:nexp_ucurr[1])
    {
      if (nexp_ucurr[1] == 1)
      {
        ui_curr[[nexp_ucurr[1]]][i, 1] <- max(u)
        nseg_ucurr[[nexp_ucurr[1]]][i, 1] <- nrep
      } else if (i == 1)
      {
        nposs <- nu_unique - nexp_ucurr[1] * umin + 1
        ui_curr[[nexp_ucurr[1]]][i, 1] <- u_sort[umin + sample(1:nposs, 1, replace = T) -
                                                   1]
        nseg_ucurr[[nexp_ucurr[1]]][i, 1] <- sum(u <= ui_curr[[nexp_ucurr[1]]][1, 1])
      } else if (i > 1 && i < nexp_ucurr[1])
      {
        nposs <- nu_unique - which(u_sort == ui_curr[[nexp_ucurr[1]]][i - 1, 1]) -
          umin * (nexp_ucurr[1] - i + 1) + 1
        ui_curr[[nexp_ucurr[1]]][i, 1] <- u_sort[umin + sample(1:nposs, 1, replace = T) +
                                                   which(u_sort == ui_curr[[nexp_ucurr[1]]][i - 1, 1]) - 1]
        nseg_ucurr[[nexp_ucurr[1]]][i, 1] <- sum(u <= ui_curr[[nexp_ucurr[1]]][i, 1]) -
          sum(u <= ui_curr[[nexp_ucurr[1]]][i - 1, 1])
      } else
      {
        ui_curr[[nexp_ucurr[1]]][i, 1] <- max(u)
        nseg_ucurr[[nexp_ucurr[1]]][i, 1] <- nrep - sum(u <= ui_curr[[nexp_ucurr[1]]][i -
                                                                                        1, 1])
      }
    }
    
    
    # beta (current)
    
    for (j in 1:nexp_tcurr[1])
    {
      for (i in 1:nexp_ucurr[1])
      {
        postbeta_curr <- postbeta(
          j,
          i,
          nseg_tcurr[[nexp_tcurr[1]]][j, 1],
          x,
          u,
          xi_curr[[nexp_tcurr[1]]][, 1],
          ui_curr[[nexp_ucurr[1]]][, 1],
          tausq_curr[[nexp_tcurr[1], nexp_ucurr[1]]][j, i, 1],
          nbasis,
          sigmasqalpha
        )
        beta_curr[[nexp_tcurr[1], nexp_ucurr[1]]][, j, i, 1] <- rmvnorm(1, postbeta_curr$beta_mean, 0.5 * (postbeta_curr$beta_var + t(postbeta_curr$beta_var)))
      }
    }
    
    
    met_rat_t <- rep(0, nloop)
    met_rat_u <- rep(0, nloop)
    epsilon_t <- rep(0, nloop)
    epsilon_u <- rep(0, nloop)
    
    
    ### movement
    for (p in 1:nloop)
      # opt.nloop=5000
    {
      if (Rev_Jump_t == 1)
      {
        # BETWEEN MODEL MOVE
        # Number of available segments
        kk <- length(which(nseg_tcurr[[nexp_tcurr[p]]] >= 2 * tmin))  #  how many segments contain time points larger than 2*opt.tmin
        # Deciding on birth or death
        
        move_t <- move(kk, nexp_tcurr[p], nexp_tmax)
        
        if (move_t$nexp_prop < nexp_tcurr[p])
        {
          # btwn=-1;
          # Death
          
          death_time <- death_t(
            x,
            u,
            nexp_tcurr[p],
            nexp_ucurr[p],
            move_t$nexp_prop,
            matrix(tausq_curr[[nexp_tcurr[p], nexp_ucurr[p]]][, , p], nexp_tcurr[p], nexp_ucurr[p]),
            xi_curr[[nexp_tcurr[p]]][, p],
            ui_curr[[nexp_ucurr[p]]][, p],
            nseg_tcurr[[nexp_tcurr[p]]][, p],
            array(beta_curr[[nexp_tcurr[p], nexp_ucurr[p]]][, , , p], dim = c(nbeta, nexp_tcurr[p], nexp_ucurr[p])),
            move_t$log_move_curr,
            move_t$log_move_prop,
            nbasis,
            sigmasqalpha,
            tmin
          )
          
          met_rat_t[p] <- death_time$met_rat
          nseg_tprop <- death_time$nseg_prop
          xi_prop <- death_time$xi_prop
          tau_prop <- death_time$tau_prop
          beta_prop <- death_time$beta_prop
          
        }
        else if (move_t$nexp_prop > nexp_tcurr[p])
        {
          # btwn=1;
          # Birth
          
          birth_time <-
            birth_t(
              x,
              u,
              nexp_tcurr[p],
              nexp_ucurr[p],
              move_t$nexp_prop,
              matrix(tausq_curr[[nexp_tcurr[p], nexp_ucurr[p]]][, , p], nexp_tcurr[p], nexp_ucurr[p]),
              xi_curr[[nexp_tcurr[p]]][, p],
              ui_curr[[nexp_ucurr[p]]][, p],
              nseg_tcurr[[nexp_tcurr[p]]][, p],
              array(beta_curr[[nexp_tcurr[p], nexp_ucurr[p]]][, , , p], dim = c(nbeta, nexp_tcurr[p], nexp_ucurr[p])),
              move_t$log_move_curr,
              move_t$log_move_prop,
              nbasis,
              sigmasqalpha,
              tmin
            )
          
          met_rat_t[p] <- birth_time$met_rat
          nseg_tprop <- birth_time$nseg_prop
          xi_prop <- birth_time$xi_prop
          tau_prop <- birth_time$tau_prop
          beta_prop <- birth_time$beta_prop
          
        }
        else
        {
          # btwn=0;
          xi_prop <- xi_curr[[nexp_tcurr[p]]][, p]
          nseg_tprop <- nseg_tcurr[[nexp_tcurr[p]]][, p]
          tau_prop <- as.matrix(tausq_curr[[nexp_tcurr[p], nexp_ucurr[p]]][, , p])
          beta_prop <- array(beta_curr[[nexp_tcurr[p], nexp_ucurr[p]]][, , , p], dim = c(nbeta, nexp_tcurr[p], nexp_ucurr[p]))
          met_rat_t[p] <- 1
        }
        
        # Update current values
        r <- runif(1)
        if (r < met_rat_t[p])
        {
          nexp_tcurr[p + 1] <- move_t$nexp_prop
          xi_curr[[nexp_tcurr[p + 1]]][, (p + 1)] <- xi_prop
          nseg_tcurr[[nexp_tcurr[p + 1]]][, (p + 1)] <- nseg_tprop
          beta_curr[[nexp_tcurr[p + 1], nexp_ucurr[p]]][, , , (p + 1)] <- beta_prop
          tausq_curr[[nexp_tcurr[p + 1], nexp_ucurr[p]]][, , (p + 1)] <- tau_prop
        } else {
          nexp_tcurr[p + 1] <- nexp_tcurr[p]
          xi_curr[[nexp_tcurr[p + 1]]][, (p + 1)] <- xi_curr[[nexp_tcurr[p +
                                                                           1]]][, p]
          nseg_tcurr[[nexp_tcurr[p + 1]]][, (p + 1)] <- nseg_tcurr[[nexp_tcurr[p +
                                                                                 1]]][, p]
          beta_curr[[nexp_tcurr[p + 1], nexp_ucurr[p]]][, , , (p + 1)] <- beta_curr[[nexp_tcurr[p +
                                                                                                  1], nexp_ucurr[p]]][, , , p]
          tausq_curr[[nexp_tcurr[p + 1], nexp_ucurr[p]]][, , (p + 1)] <- tausq_curr[[nexp_tcurr[p +
                                                                                                  1], nexp_ucurr[p]]][, , p]
          
        }
        
      } else {
        nexp_tcurr[p + 1] <- nexp_tcurr[p]
        xi_curr[[nexp_tcurr[p + 1]]][, (p + 1)] <- xi_curr[[nexp_tcurr[p +
                                                                         1]]][, p]
        nseg_tcurr[[nexp_tcurr[p + 1]]][, (p + 1)] <- nseg_tcurr[[nexp_tcurr[p +
                                                                               1]]][, p]
        beta_curr[[nexp_tcurr[p + 1], nexp_ucurr[p]]][, , , (p + 1)] <- beta_curr[[nexp_tcurr[p +
                                                                                                1], nexp_ucurr[p]]][, , , p]
        tausq_curr[[nexp_tcurr[p + 1], nexp_ucurr[p]]][, , (p + 1)] <- tausq_curr[[nexp_tcurr[p +
                                                                                                1], nexp_ucurr[p]]][, , p]
        
      }
      
      # WITHIN MODEL MOVE
      # Drawing a new cut point and betas simultaneously
      # First draw the size of the move
      
      
      within_time <- within_t(
        x,
        u,
        nexp_tcurr[p + 1],
        nexp_ucurr[p],
        xi_curr[[nexp_tcurr[p + 1]]][, (p + 1)],
        ui_curr[[nexp_ucurr[p]]][, p],
        array(beta_curr[[nexp_tcurr[p + 1], nexp_ucurr[p]]][, , , (p + 1)], dim = c(nbeta, nexp_tcurr[p +
                                                                                                        1], nexp_ucurr[p])),
        nseg_tcurr[[nexp_tcurr[p + 1]]][, (p + 1)],
        matrix(tausq_curr[[nexp_tcurr[p + 1], nexp_ucurr[p]]][, , (p + 1)], nexp_tcurr[p +
                                                                                         1], nexp_ucurr[p]),
        nbasis,
        sigmasqalpha,
        tmin,
        prob_mm1
      )
      
      r <- runif(1)
      if (r < within_time$epsilon || p == 1)
      {
        epsilon_t[p] <- within_time$epsilon
        xi_curr[[nexp_tcurr[p + 1]]][, (p + 1)] <- within_time$xi_prop
        nseg_tcurr[[nexp_tcurr[p + 1]]][, (p + 1)] <- within_time$nseg_new
        beta_curr[[nexp_tcurr[p + 1], nexp_ucurr[p]]][, , , (p + 1)] <- within_time$beta_prop
      } else
      {
        for (j in 1:nexp_tcurr[p + 1])
        {
          for (i in 1:nexp_ucurr[p])
          {
            
            postbeta_18 <- postbeta(
              j,
              i,
              nseg_tcurr[[nexp_tcurr[p + 1]]][j, (p + 1)],
              x,
              u,
              xi_curr[[nexp_tcurr[p + 1]]][, (p + 1)],
              ui_curr[[nexp_ucurr[p]]][, p],
              tausq_curr[[nexp_tcurr[p + 1], nexp_ucurr[p]]][j, i, (p + 1)],
              nbasis,
              sigmasqalpha
            )
            beta_curr[[nexp_tcurr[p + 1], nexp_ucurr[p]]][, j, i, (p + 1)] <- 
              rmvnorm(1, postbeta_18$beta_mean, 0.5 * (postbeta_18$beta_var + t(postbeta_18$beta_var)))
          }
        }
        xi_curr[[nexp_tcurr[p + 1]]][, (p + 1)] <-
          xi_curr[[nexp_tcurr[p + 1]]][, (p + 1)]
        nseg_tcurr[[nexp_tcurr[p + 1]]][, (p + 1)] <-
          nseg_tcurr[[nexp_tcurr[p + 1]]][, (p + 1)]
        
      }
      
      
      # Drawing tau
      for (j in 1:nexp_tcurr[p + 1])
      {
        for (i in 1:nexp_ucurr[p])
        {
          tau_a <- nbasis / 2 + tau_prior_a
          tau_b <- sum(beta_curr[[nexp_tcurr[p + 1], nexp_ucurr[p]]][2:nbeta, j, i, (p +
                                                                                       1)] ^ 2) / 2 + tau_prior_b
          tausq_curr[[nexp_tcurr[p + 1], nexp_ucurr[p]]][j, i, (p + 1)] <- 
            1 / rgamma(1, tau_a, tau_b)
        }
        
      }
      
      
      if (Rev_Jump_u == 1)
      {
        # BETWEEN MODEL MOVE
        # Number of available segments (taking duplicate covariate values into consideration)
        kk <- 0
        for (i in 1:nexp_ucurr[p])
        {
          if (i == 1)
          {
            kk <- kk + (sum(u_sort <= ui_curr[[nexp_ucurr[p]]][i, p]) >= 2 * umin)
          } else
          {
            kk <- kk + (sum(u_sort <= ui_curr[[nexp_ucurr[p]]][i, p] &
                              u_sort > ui_curr[[nexp_ucurr[p]]][i - 1, p]) >= 2 * umin)
          }
        }
        
        # Deciding on birth or death
        
        move_u <- move(kk, nexp_ucurr[p], nexp_umax)
        if (move_u$nexp_prop < nexp_ucurr[p])
        {
          # btwn=-1;
          # Death
          
          death_uvar <- death_u(
            x,
            u,
            nexp_tcurr[p + 1],
            nexp_ucurr[p],
            move_u$nexp_prop,
            matrix(tausq_curr[[nexp_tcurr[p + 1], nexp_ucurr[p]]][, , p + 1], nexp_tcurr[p +
                                                                                           1], nexp_ucurr[p]),
            xi_curr[[nexp_tcurr[p + 1]]][, p + 1],
            ui_curr[[nexp_ucurr[p]]][, p],
            nseg_tcurr[[nexp_tcurr[p + 1]]][, p + 1],
            nseg_ucurr[[nexp_ucurr[p]]][, p],
            array(beta_curr[[nexp_tcurr[p + 1], nexp_ucurr[p]]][, , , p +
                                                                  1], dim = c(nbeta, nexp_tcurr[p + 1], nexp_ucurr[p])),
            move_u$log_move_curr,
            move_u$log_move_prop,
            nbasis,
            sigmasqalpha,
            umin
          )
          
          ui_prop <- death_uvar$ui_prop
          nseg_uprop <- death_uvar$nseg_prop
          tau_prop <- death_uvar$tau_prop
          beta_prop <- death_uvar$beta_prop
          met_rat_u[p] <- death_uvar$met_rat
          
        } else if (move_u$nexp_prop > nexp_ucurr[p])
        {
          # btwn=1;
          # Birth,
          
          birth_uvar <- birth_u(
            x,
            u,
            nexp_tcurr[p + 1],
            nexp_ucurr[p],
            move_u$nexp_prop,
            matrix(tausq_curr[[nexp_tcurr[p + 1], nexp_ucurr[p]]][, , p + 1], nexp_tcurr[p +
                                                                                           1], nexp_ucurr[p]),
            xi_curr[[nexp_tcurr[p + 1]]][, p + 1],
            ui_curr[[nexp_ucurr[p]]][, p],
            nseg_tcurr[[nexp_tcurr[p + 1]]][, p + 1],
            nseg_ucurr[[nexp_ucurr[p]]][, p],
            array(beta_curr[[nexp_tcurr[p + 1], nexp_ucurr[p]]][, , , p +
                                                                  1], dim = c(nbeta, nexp_tcurr[p + 1], nexp_ucurr[p])),
            move_u$log_move_curr,
            move_u$log_move_prop,
            nbasis,
            sigmasqalpha,
            umin
          )
          
          ui_prop <- birth_uvar$ui_prop
          nseg_uprop <- birth_uvar$nseg_prop
          tau_prop <- birth_uvar$tau_prop
          beta_prop <- birth_uvar$beta_prop
          met_rat_u[p] <- birth_uvar$met_rat
          
        } else
        {
          # btwn=0;
          ui_prop <- ui_curr[[nexp_ucurr[p]]][, p]
          nseg_uprop <- nseg_ucurr[[nexp_ucurr[p]]][, p]
          tau_prop <- as.matrix(tausq_curr[[nexp_tcurr[p + 1], nexp_ucurr[p]]][, , p +
                                                                                 1])
          beta_prop <- array(beta_curr[[nexp_tcurr[p + 1], nexp_ucurr[p]]][, , , p +
                                                                             1], dim = c(nbeta, nexp_tcurr[p + 1], nexp_ucurr[p]))
          met_rat_u[p] <- 1
          
        }
        
        r <- runif(1)
        if (r < met_rat_u[p])
        {
          nexp_ucurr[p + 1] <- move_u$nexp_prop
          ui_curr[[nexp_ucurr[p + 1]]][, (p + 1)] <- ui_prop
          nseg_ucurr[[nexp_ucurr[p + 1]]][, (p + 1)] <- nseg_uprop
          tausq_curr[[nexp_tcurr[p + 1], nexp_ucurr[p + 1]]][, , (p + 1)] <- 
            tau_prop
          beta_curr[[nexp_tcurr[p + 1], nexp_ucurr[p + 1]]][, , , (p + 1)] <- 
            beta_prop
          
        } else{
          nexp_ucurr[p + 1] <- nexp_ucurr[p]
          ui_curr[[nexp_ucurr[p + 1]]][, (p + 1)] <- ui_curr[[nexp_ucurr[p +
                                                                           1]]][, p]
          nseg_ucurr[[nexp_ucurr[p + 1]]][, (p + 1)] <- nseg_ucurr[[nexp_ucurr[p +
                                                                                 1]]][, p]
          tausq_curr[[nexp_tcurr[p + 1], nexp_ucurr[p + 1]]][, , (p + 1)] <- 
            tausq_curr[[nexp_tcurr[p + 1], nexp_ucurr[p + 1]]][, , p]
          beta_curr[[nexp_tcurr[p + 1], nexp_ucurr[p + 1]]][, , , (p + 1)] <- 
            beta_curr[[nexp_tcurr[p + 1], nexp_ucurr[p + 1]]][, , , p]
        }
        
      } else {
        nexp_ucurr[p + 1] <- nexp_ucurr[p]
        ui_curr[[nexp_ucurr[p + 1]]][, (p + 1)] <- ui_curr[[nexp_ucurr[p +
                                                                         1]]][, p]
        nseg_ucurr[[nexp_ucurr[p + 1]]][, (p + 1)] <- nseg_ucurr[[nexp_ucurr[p +
                                                                               1]]][, p]
        beta_curr[[nexp_tcurr[p + 1], nexp_ucurr[p + 1]]][, , , (p + 1)] <- 
          beta_curr[[nexp_tcurr[p + 1], nexp_ucurr[p + 1]]][, , , p]
        tausq_curr[[nexp_tcurr[p + 1], nexp_ucurr[p + 1]]][, , (p + 1)] <- 
          tausq_curr[[nexp_tcurr[p + 1], nexp_ucurr[p + 1]]][, , p]
        
      }
      
      
      # WITHIN MODEL MOVE
      # Drawing a new cut point and betas simultaneously
      # First draw the size of the move
      
      within_uvar <- within_u(
        x,
        u,
        nexp_tcurr[p + 1],
        nexp_ucurr[p + 1],
        xi_curr[[nexp_tcurr[p + 1]]][, (p + 1)],
        ui_curr[[nexp_ucurr[p + 1]]][, (p + 1)],
        array(beta_curr[[nexp_tcurr[p + 1], nexp_ucurr[p + 1]]][, , , (p + 1)], dim = c(nbeta, nexp_tcurr[p +
                                                                                                            1], nexp_ucurr[p + 1])),
        nseg_tcurr[[nexp_tcurr[p + 1]]][, (p + 1)],
        nseg_ucurr[[nexp_ucurr[p + 1]]][, (p + 1)],
        matrix(tausq_curr[[nexp_tcurr[p + 1], nexp_ucurr[p +
                                                           1]]][, , (p + 1)], nexp_tcurr[p + 1], nexp_ucurr[p + 1]),
        nbasis,
        sigmasqalpha,
        umin,
        prob_mm1_u
      )
      
      r <- runif(1)
      if (r < within_uvar$epsilon || p == 1)
      {
        epsilon_u[p] <- within_uvar$epsilon
        ui_curr[[nexp_ucurr[p + 1]]][, (p + 1)] <- within_uvar$ui_prop
        nseg_ucurr[[nexp_ucurr[p + 1]]][, (p + 1)] <- within_uvar$nseg_new
        beta_curr[[nexp_tcurr[p + 1], nexp_ucurr[p + 1]]][, , , (p + 1)] <- 
          within_uvar$beta_prop
      } else
      {
        for (j in 1:nexp_tcurr[p + 1])
        {
          for (i in 1:nexp_ucurr[p + 1])
          {
            postbeta_19 <- postbeta(
              j,
              i,
              nseg_tcurr[[nexp_tcurr[p + 1]]][j, (p + 1)],
              x,
              u,
              xi_curr[[nexp_tcurr[p + 1]]][, (p + 1)],
              ui_curr[[nexp_ucurr[p + 1]]][, (p + 1)],
              tausq_curr[[nexp_tcurr[p + 1], nexp_ucurr[p + 1]]][j, i, (p + 1)],
              nbasis,
              sigmasqalpha
            )
            beta_curr[[nexp_tcurr[p + 1], nexp_ucurr[p + 1]]][, j, i, (p +
                                                                         1)] <- rmvnorm(1, postbeta_19$beta_mean, 0.5 * (postbeta_19$beta_var + t(postbeta_19$beta_var)))
          }
          
        }
        ui_curr[[nexp_ucurr[p + 1]]][, (p + 1)] <- ui_curr[[nexp_ucurr[p +
                                                                         1]]][, (p + 1)]
        nseg_ucurr[[nexp_ucurr[p + 1]]][, (p + 1)] <- nseg_ucurr[[nexp_ucurr[p +
                                                                               1]]][, (p + 1)]
        
      }
      
      # Drawing tau
      for (j in 1:nexp_tcurr[p + 1])
      {
        for (i in 1:nexp_ucurr[p + 1])
        {
          tau_a <- nbasis / 2 + tau_prior_a
          tau_b <- sum(beta_curr[[nexp_tcurr[p + 1], nexp_ucurr[p + 1]]][2:nbeta, j, i, (p +
                                                                                           1)] ^ 2) / 2 + tau_prior_b
          tausq_curr[[nexp_tcurr[p + 1], nexp_ucurr[p + 1]]][j, i, (p + 1)] <- 
            1 / rgamma(1, tau_a, tau_b)
        }
        
      }
      
      # Estimating Spectral Density
      for (j in 1:nexp_tcurr[p + 1])
      {
        for (i in 1:nexp_ucurr[p + 1])
        {
          log_spec_hat_curr[[nexp_tcurr[p + 1], nexp_ucurr[p + 1]]][, j, i, (p + 1)] <- 
            nu_mat_hat %*% as.matrix(beta_curr[[nexp_tcurr[p + 1], nexp_ucurr[p + 1]]][, j, i, p +
                                                                                         1])      ####
        }
        
      }
      
      
      cat("p:", p, "\n")
      cat("t:", nexp_tcurr[p + 1], "\n")
      cat("u:", nexp_ucurr[p + 1], "\n")
      cat("xi:", xi_curr[[nexp_tcurr[p + 1]]][, (p + 1)], "\n")
      cat("ui:", ui_curr[[nexp_ucurr[p + 1]]][, (p + 1)], "\n")
      
    }
    
    
    
    ### plot
    
    if (plotting == TRUE)
    {
      post_prob_nexp_t=matrix(0,nexp_tmax,1)
      for (k in 1:nexp_tmax)
      {
        kk_t <- which(nexp_tcurr[(nwarmup + 1):(nloop+1)] == k)
        post_prob_nexp_t[k]=length(kk_t)/(nloop-nwarmup)
        
        if (length(kk_t) != 0)
        {
          xi_mat=matrix(0,length(kk_t),k)
          for (g in 1:length(kk_t))
          {
            xi_mat[g,]=xi_curr[[k]][,kk_t[g]+nwarmup]
          }
          for (g in 1:(k-1))
          {
            plot(xi_mat[,g],type = "l",main = paste("Time Partition Points for a Partition of",k))
          }
          for (g in 1:(k-1))
          {
            hist(xi_mat[,g],main = paste("Distribution of",g,"Time Partition Point for a mixture of",k))
          }
          
        }
        
      }
      barplot(t(post_prob_nexp_t),main = "Posterior Probability of Number of Time Segments")
      
      post_prob_nexp_u <- matrix(0,nexp_umax,1)
      for (k in 1:nexp_umax)
      {
        kk_u <- which(nexp_ucurr[(nwarmup + 1):(nloop+1)] == k)
        post_prob_nexp_u[k]=length(kk_u)/(nloop-nwarmup)
        
        if (length(kk_u) != 0)
        {
          ui_mat=matrix(0,length(kk_u),k)
          for (g in 1:length(kk_u))
          {
            ui_mat[g,]=ui_curr[[k]][,kk_u[g]+nwarmup]
          }
          for (g in 1:(k-1))
          {
            plot(ui_mat[,g],type = 'l', main = paste("Covariate Partition Points for a Partition of",k))
          }
          for (g in 1:(k-1))
          {
            hist(ui_mat[,g],main = paste("Distribution of",g,"Covariate Partition Point for a mixture of",k))
          }
          
        }
        
      }
      barplot(t(post_prob_nexp_u), main = "Posterior Probability of Number of Covariate Segments")
      
    }
    
    
    
    z <- list(xi = xi_curr, ui=ui_curr, log_spec_hat_curr = log_spec_hat_curr, nloop = nloop,
              nwarmup = nwarmup, nexp_tcurr = nexp_tcurr, nexp_ucurr = nexp_ucurr, x = x, u = u,
              nseg_tcurr = nseg_tcurr, nseg_ucurr = nseg_ucurr,
              nexp_tmax = nexp_tmax, nexp_umax = nexp_umax, tmin = tmin, umin = umin, sigmasqalpha = sigmasqalpha,
              tau_prior_a = tau_prior_a, tau_prior_b = tau_prior_b, prob_mm1 = prob_mm1,
              var_inflate = var_inflate, nbasis = nbasis, nfreq_hat = nfreq_hat)
    return(z)
    
  }



#########################
# 3. call CABS function #
#########################
cabs(x,u,nloop=100,seed=6666,nwarmup=20,nexp_tmax=10,nexp_umax=8,tmin=50,umin=1,nbasis=7,
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

