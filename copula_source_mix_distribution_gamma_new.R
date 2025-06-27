library(rugarch)
library(xts)
library(ADGofTest)
library(qqtest)
library(copula)
library(qrmdata)
library(qrmtools)
library(LaplacesDemon)
library(EnvStats)
############################### t- Copula##################################

MCMC_tcopula <- function(N = 10000,
                         kappa = 0.1,
                         component = 10,
                         alpha = 0.2,
                         sda = 0.05,
                         sdb = 0.05,
                         data_source) {
  
  #### 1) Initialization ####
  data <- data_source  # Make sure data_source > 0 if using Gamma
  
  n  <- nrow(data)
  d  <- ncol(data)
  
  v_acceptance_rate <- 0
  sigma_acceptance_rate <- 0
  vv0 <- numeric(N)
  K   <- numeric(N)
  
  # Initialize cluster labels:
  k_init    <- rep(seq_len(5), length.out = n)
  k_matrix  <- matrix(0, nrow = N, ncol = n)
  k_matrix[1, ] <- k_init
  
  # --- Gamma parameters: shape and rate ---
  dgamma_shape <- array(1, dim = c(N, d, component))
  dgamma_rate  <- array(1, dim = c(N, d, component))
  
  # --- Correlation / Covariance arrays ---
  rho <- array(NA, dim = c(d, d, component, N))
  rho[,,, 1] <- diag(d)
  
  cov <- rho
  Lij <- cov
  
  # --- Degrees of freedom for t-copula ---
  v_init    <- sapply(rep(0, component), function(x) 1 + rgamma(1, shape = 5, rate = 1))
  v_matrix  <- matrix(0, nrow = N, ncol = component)
  v_matrix[1, ] <- v_init
  
  # --- DP concentration parameter ---
  alpha_matrix <- rep(alpha, N)
  
  # --- rate_matrix for slice threshold ---
  #   = (1-kappa) * kappa^(j-1)
  rate_matrix <- (1 - kappa) * kappa^(0:(component - 1))
  
 
  
  # --- Stick-breaking weights ---
  V_matrix <- rep(1 - kappa, component)
  V_matrix[component] <- 1  # The last stick
  w_matrix_init <-  rep(c(1,0),c(5,component-5))
  w_matrix      <- matrix(0, nrow = N, ncol = component)
  w_matrix[1, ] <- w_matrix_init/sum(w_matrix_init)
  
  # --- Slice variable init ---
  z_init    <- runif(n,min=0,max = 1/5)
  z_matrix  <- matrix(0, nrow = N, ncol = n)
  z_matrix[1, ] <- z_init
  
  
  
  # --- Bookkeeping for output ---
  y_pred    <- matrix(0, nrow = N, ncol = d)
  K_number  <- numeric(N)
  
  #### 2) MCMC Sampling ####
  for (t in seq_len(N - 1)) {
    
    ########## (a) Update correlation in each used cluster ##########
    for (ii in unique(k_matrix[t, ])) {
      # For correlation, we do random-walk proposals in the Cholesky factor:
      for (iii in 2:d) {
        for (jjj in 1:(iii - 1)) {
          
          updating_index <- which(k_matrix[t, ] == ii)
          
          propose <- cov[,, ii, t]
          # Propose a small random step in the lower-triangular entry:
          propose[iii, jjj] <- rnorm(1, mean = cov[iii, jjj, ii, t], sd = 0.05)
          
          # Old inverse covariance:
          old_L   <- cov[,, ii, t]
          mpo     <- solve(old_L %*% t(old_L))
          
          # Proposed inverse covariance:
          mp      <- solve(propose %*% t(propose))
          
          # Convert these to correlation matrices:
          Sigma_new <- (diag(mp)^(-1) * diag(d))^(1/2) %*% mp %*% (diag(mp)^(-1) * diag(d))^(1/2)
          Sigma_old <- (diag(mpo)^(-1) * diag(d))^(1/2) %*% mpo %*% (diag(mpo)^(-1) * diag(d))^(1/2)
          
          new_correlation_param <- P2p(Sigma_new)
          correlation_param     <- P2p(Sigma_old)
          
          # Build t-copulas for old vs new
          tcop_new <- ellipCopula("t", new_correlation_param,
                                  dim = d, df = v_matrix[t, ii], dispstr = "un")
          tcop_old <- ellipCopula("t", correlation_param,
                                  dim = d, df = v_matrix[t, ii], dispstr = "un")
          
          # Evaluate copula likelihood of the subset of data in cluster ii:
          dup <- sapply(seq_len(d), function(x) 
            pgamma(data[updating_index, x],
                   shape = dgamma_shape[t, x, ii],
                   rate  = dgamma_rate[t, x, ii]))
          
          cop_new <- dCopula(dup, tcop_new, log = TRUE)
          cop_old <- dCopula(dup, tcop_old, log = TRUE)
          
          dcopula_new <- sum(cop_new[!is.na(cop_new)])
          dcopula_old <- sum(cop_old[!is.na(cop_old)])
          
          # Prior on correlation entry ~ Normal(0, 0.5^2):
          dcorr_old <- dnorm(cov[iii, jjj, ii, t], mean = 0, sd = 0.5, log = TRUE)
          dcorr_new <- dnorm(propose[iii, jjj],      mean = 0, sd = 0.5, log = TRUE)
          
          acceptance_sigma <- (dcopula_new - dcopula_old) + (dcorr_new - dcorr_old)
          
          if (is.finite(acceptance_sigma) &&
              runif(1) < exp(acceptance_sigma)) {
            cov[,, ii, t] <- propose
            # Optionally store Sigma_new or not
         #   cat("update corr\n")
          } else {
         #   cat("no update corr\n")
          }
        }
      }
      
      # After iterating over all off-diagonal proposals, copy forward to t+1
      Lij[,, ii, t + 1] <- cov[,, ii, t + 1] <- cov[,, ii, t]
      newL     <- Lij[,, ii, t + 1]
      news_inv <- solve(newL %*% t(newL))
      rho[,, ii, t + 1] <- (diag(news_inv)^(-1) * diag(d))^(1/2) %*% news_inv %*%
        (diag(news_inv)^(-1) * diag(d))^(1/2)
    }
    
    # For clusters not used, carry over old values:
    no_up <- setdiff(seq_len(component), unique(k_matrix[t, ]))
    Lij[,, no_up, t + 1] <- cov[,, no_up, t + 1] <- cov[,, no_up, t]
    rho[,, no_up, t + 1] <- rho[,, no_up, t]
    
    ########## (b) Update degrees of freedom v for used clusters ##########
    for (ii in unique(k_matrix[t, ])) {
      updating_index <- which(k_matrix[t, ] == ii)
      
      # Propose log(df-1):
      proposal_r <- runif(1,
                          min = log(v_matrix[t, ii] - 1) - 0.1,
                          max = log(v_matrix[t, ii] - 1) + 0.1)
      proposal_r <- exp(proposal_r) + 1
      
      correlation_param <- P2p(rho[,, ii, t + 1])
      
      tcop_new <- ellipCopula("t", correlation_param,
                              dim = d, df = proposal_r, dispstr = "un")
      tcop_old <- ellipCopula("t", correlation_param,
                              dim = d, df = v_matrix[t, ii], dispstr = "un")
      
      # Evaluate copula likelihood
      dup <- sapply(seq_len(d), function(x)
        pgamma(data[updating_index, x],
               shape = dgamma_shape[t, x, ii],
               rate  = dgamma_rate[t, x, ii]))
      
      cop_new <- dCopula(dup, tcop_new, log = TRUE)
      cop_old <- dCopula(dup, tcop_old, log = TRUE)
      
      dcopula_new <- sum(cop_new[!is.na(cop_new)])
      dcopula_old <- sum(cop_old[!is.na(cop_old)])
      
      # Log proposal densities in the df-1 space:
      dprop_top    <- log(proposal_r - 1)
      dprop_bottom <- log(v_matrix[t, ii] - 1)
      
      # Gamma prior on (df - 1):
      d_v_prior_new <- dgamma(proposal_r - 1, shape = 5, rate = 1, log = TRUE)
      d_v_prior_old <- dgamma(v_matrix[t, ii] - 1, shape = 5, rate = 1, log = TRUE)
      
      acceptance_r <- (dcopula_new - dcopula_old) +
        (dprop_top - dprop_bottom) +
        (d_v_prior_new - d_v_prior_old)
      
      if (is.finite(acceptance_r) &&
          runif(1) < exp(acceptance_r)) {
        v_matrix[t + 1, ii] <- proposal_r
        # cat("accept df!\n")
      } else {
        v_matrix[t + 1, ii] <- v_matrix[t, ii]
        # cat("refuse df!\n")
      }
    }
    # For clusters not used, carry forward v:
    v_matrix[t + 1, no_up] <- v_matrix[t, no_up]
    
    ########## (c) Slice variable & cluster re-assignment ##########
    z_matrix[t + 1, ] <- sapply(k_matrix[t, ],
                                function(x) runif(1, 0, w_matrix[t, x]))
    
    for (i in seq_len(n)) {
      probs <- numeric(component)
      
      # Among clusters j s.t. z_i <= rate_matrix[j], compute unnormalized probs
      feasible_js <- which(z_matrix[t + 1, i] < w_matrix[t, ])
     # print(feasible_js)
      for (j in feasible_js) {
        # transform data to CDF under Gamma params
        dup <- sapply(seq_len(d), function(x)
          pgamma(data[i, x],
                 shape = dgamma_shape[t, x, j],
                 rate  = dgamma_rate[t, x, j]))
        
        # Copula part
        cop_val <- dCopula(dup,
                           copula = tCopula(param = P2p(rho[,, j, t + 1]),
                                            dim = d,
                                            dispstr = "un",
                                            df  = v_matrix[t + 1, j]),
                           log = FALSE)  # We'll take log manually if we want
        # Marginal part
        marg_val <- prod(sapply(seq_len(d), function(x)
          dgamma(data[i, x],
                 shape = dgamma_shape[t, x, j],
                 rate  = dgamma_rate[t, x, j])))
        
        # Multiply by w_matrix[t, j] * 1/rate_matrix[j], as you had:
        probs[j] <- (1) * cop_val * marg_val
        
      }
      
      probs[is.na(probs)] <- 1e-300
     # print(probs)
      # sample a new cluster label among feasible_js
      if(length(feasible_js)==1){k_matrix[t+1,i] <- feasible_js}else if(length(feasible_js)==0){
        k_matrix[t + 1, i] <- 1}else{
          if(sum(probs[feasible_js])<=0){probs[feasible_js] <- probs[feasible_js] + 1e-300}
      k_matrix[t + 1, i] <- sample(feasible_js, size = 1,
                                   replace = TRUE,
                                   prob = probs[feasible_js])}
    }
    
    ########## (d) Stick-Breaking update of w ##########
    for (i in seq_len(component)) {
      V_matrix[i] <- rbeta(
        1,
        1 + sum(k_matrix[t + 1, ] == i),
        alpha_matrix[t + 1] + sum(k_matrix[t + 1, ] > i)
      )
    }
    w_matrix[t + 1, 1] <- V_matrix[1]
    for (i in 2:component) {
      w_matrix[t + 1, i] <- V_matrix[i] * prod(1 - V_matrix[1:(i - 1)])
    }
    
    ########## (e) Predict y_pred from a randomly chosen cluster ##########
    cs <- sample(seq_len(component),
                 size = 1,
                 replace = FALSE,
                 prob = w_matrix[t + 1, ])
    y_pred[t + 1, ] <- rCopula(
      1,
      copula = tCopula(param = P2p(rho[,, cs, t + 1]),
                       dim = d, dispstr = "un",
                       df  = v_matrix[t + 1, cs])
    )
    y_pred[t+1,] <- sapply(seq_len(d), function(x) qgamma(y_pred[t+1,x], 
                                                          shape = dgamma_shape[t, x, cs],
                                                          rate = dgamma_rate[t,x,cs]))
    
    ########## (f) Update alpha (concentration) ##########
    a <- 1
    b <- 1
    k <- length(unique(k_matrix[t + 1, ]))
    
    eta <- rbeta(1,
                 shape1 = alpha_matrix[t] + 1,
                 shape2 = n)
    mix <- (a + k - 1) / (a + k - 1 + n * (b - log(eta)))
    
    if (is.finite(mix) && runif(1) < mix) {
      alpha <- rgamma(1, shape = a + k,     rate = b - log(eta))
    } else {
      alpha <- rgamma(1, shape = a + k - 1, rate = b - log(eta))
    }
    alpha_matrix[t + 1] <- alpha
    
    ########## (g) Gamma-parameter update (was Beta update) ##########
    dgamma_shape[t + 1, , ] <- dgamma_shape[t, , ]
    dgamma_rate[t + 1, , ]  <- dgamma_rate[t, , ]
    
    for (i_clust in unique(k_matrix[t + 1, ])) {
      updating_index <- which(k_matrix[t + 1, ] == i_clust)
      if (length(updating_index) == 0) next
      
      for (j_dim in seq_len(d)) {
        
        ################# Part 1: Update SHAPE alone #################
        # Current shape, rate (already in [t+1] arrays from previous iteration)
        a_old <- dgamma_shape[t + 1, j_dim, i_clust]
        b_old <- dgamma_rate[t + 1, j_dim, i_clust]
        
        # Proposal scale
        if (t >= 200) {
          sda <- sda  # shape proposal sd
        } else {
          sda <- sda
        }
        
        # Propose new shape in log‐space, keep rate fixed
        log_shape_curr <- log(a_old)
        shape_prop     <- rnorm(1, mean = log_shape_curr, sd = sda)
        a_prop         <- exp(shape_prop)  # Proposed shape
        
        # Construct "dup" = old CDFs for all dims in this cluster
        dup <- sapply(seq_len(d), function(x)
          pgamma(data[updating_index, x],
                 shape = dgamma_shape[t + 1, x, i_clust],
                 rate  = dgamma_rate[t + 1, x, i_clust])
        )
        # We'll replace only the j_dim column with the new shape
        data_prop <- dup
        if (!is.matrix(data_prop)) {
          data_prop <- matrix(dup,
                              nrow = length(updating_index),
                              ncol = d, byrow = FALSE)
        } else {
          data_prop <- matrix(data_prop,
                              nrow = length(updating_index),
                              ncol = d)
        }
        
        # Use the new shape, but keep old rate b_old
        data_prop[, j_dim] <- pgamma(data[updating_index, j_dim],
                                     shape = a_prop, rate = b_old)
        
        # Evaluate new vs old copula densities
        tcop_new <- tCopula(
          param   = P2p(rho[,, i_clust, t + 1]),
          dim     = d,
          dispstr = "un",
          df      = v_matrix[t + 1, i_clust]
        )
        denprop     <- dCopula(data_prop, copula = tcop_new, log = TRUE)
        denprop_sum <- sum(denprop[!is.na(denprop)])
        
        denprop_cur <- dCopula(dup, copula = tcop_new, log = TRUE)
        denprop_cur_sum <- sum(denprop_cur[!is.na(denprop_cur)])
        
        # Evaluate new vs old marginal log-likelihood for dimension j_dim
        marg_new <- dgamma(data[updating_index, j_dim],
                           shape = a_prop, rate = b_old,
                           log   = TRUE)
        marg_cur <- dgamma(data[updating_index, j_dim],
                           shape = a_old,    rate = b_old,
                           log   = TRUE)
        marg_new_sum <- sum(marg_new[!is.na(marg_new)])
        marg_cur_sum <- sum(marg_cur[!is.na(marg_cur)])
        
        # Log-prior on log(shape) (if you keep the same normal prior as before)
        prop_prior <- dnorm(shape_prop,       mean = 0, sd = 2, log = TRUE)
        cur_prior  <- dnorm(log_shape_curr,   mean = 0, sd = 2, log = TRUE)
        
        # MH acceptance ratio (rate is unchanged => no prior/likelihood change in rate)
        accp <- (denprop_sum  - denprop_cur_sum) +
          (marg_new_sum - marg_cur_sum)   +
          (prop_prior   - cur_prior)
        
        if (is.finite(accp) && runif(1) < exp(accp)) {
          # Accept the shape update
          dgamma_shape[t + 1, j_dim, i_clust] <- a_prop
          # cat("Accept shape!\n")
        } else {
          # Keep old shape
          # cat("Reject shape!\n")
        }
        
        ################# Part 2: Update RATE alone #################
        a_final <- dgamma_shape[t + 1, j_dim, i_clust]  # use the (possibly updated) shape
        b_old   <- dgamma_rate[t + 1, j_dim, i_clust]
        
        if (t >= 200) {
          sdb <- sdb
        } else {
          sdb <- sdb
        }
        
        log_rate_curr <- log(b_old)
        rate_prop     <- rnorm(1, mean = log_rate_curr, sd = sdb)
        b_prop        <- exp(rate_prop)  # Proposed rate
        
        # Rebuild 'dup' with the updated shape 'a_final' but old rate (b_old)
        # because we want a baseline. Then we only replace j_dim with new rate b_prop:
        dup <- sapply(seq_len(d), function(x)
          pgamma(data[updating_index, x],
                 shape = dgamma_shape[t + 1, x, i_clust],
                 rate  = dgamma_rate[t + 1, x, i_clust])
        )
        if (!is.matrix(dup)) {
          dup <- matrix(dup,
                        nrow = length(updating_index),
                        ncol = d, byrow = FALSE)
        } else {
          dup <- matrix(dup, nrow = length(updating_index), ncol = d)
        }
        
        data_prop <- dup
        data_prop[, j_dim] <- pgamma(data[updating_index, j_dim],
                                     shape = a_final, rate = b_prop)
        
        # Evaluate new vs old copula density
        denprop     <- dCopula(data_prop, copula = tcop_new, log = TRUE)
        denprop_sum <- sum(denprop[!is.na(denprop)])
        
        denprop_cur <- dCopula(dup, tcop_new, log = TRUE)
        denprop_cur_sum <- sum(denprop_cur[!is.na(denprop_cur)])
        
        # Evaluate new vs old marginal log-likelihood for dimension j_dim
        marg_new <- dgamma(data[updating_index, j_dim],
                           shape = a_final, rate = b_prop,
                           log   = TRUE)
        marg_cur <- dgamma(data[updating_index, j_dim],
                           shape = a_final, rate = b_old,
                           log   = TRUE)
        marg_new_sum <- sum(marg_new[!is.na(marg_new)])
        marg_cur_sum <- sum(marg_cur[!is.na(marg_cur)])
        
        # Log-prior on log(rate) only
        prop_prior <- dnorm(rate_prop,     mean = 0, sd = 2, log = TRUE)
        cur_prior  <- dnorm(log_rate_curr, mean = 0, sd = 2, log = TRUE)
        
        # MH acceptance ratio (shape is now fixed => no shape prior difference)
        accp <- (denprop_sum  - denprop_cur_sum) +
          (marg_new_sum - marg_cur_sum)   +
          (prop_prior   - cur_prior)
        
        if (is.finite(accp) && runif(1) < exp(accp)) {
          dgamma_rate[t + 1, j_dim, i_clust] <- b_prop
          # cat("Accept rate!\n")
        } else {
          # cat("Reject rate!\n")
        }
        
      } # end for j_dim
    } # end for i_clust
    
    ########## (h) Post-processing for label-switching (optional) ##########
    # Re-order by some criterion:
    ord <- order( -( w_matrix[t + 1, ] + (1 - 0.7) * 0.7^((1:component) - 1) ) )
    
    k_matrix[t + 1, ] <- sapply(k_matrix[t + 1, ], function(x) which(ord == x))
    
    v_matrix[t + 1, ]                 <- v_matrix[t + 1, ord]
    rho[,, seq_len(component), t + 1] <- rho[,, ord, t + 1]
    cov[,, seq_len(component), t + 1] <- cov[,, ord, t + 1]
    Lij[,, seq_len(component), t + 1] <- Lij[,, ord, t + 1]
    
    dgamma_shape[t + 1, , seq_len(component)] <-
      dgamma_shape[t + 1, , ord]
    dgamma_rate[t + 1, , seq_len(component)]  <-
      dgamma_rate[t + 1, , ord]
    w_matrix[t + 1, ] <- w_matrix[t + 1, ord]
    
    ########## (i) Diagnostics / simple plots ##########
    K_number[t + 1] <- length(unique(k_matrix[t + 1, ]))
    # cat("k_matrix", k_matrix, "\n")
    cat("k = 1 ", sum(k_matrix[t + 1, ] == 1),
        " ; k = 2 ", sum(k_matrix[t + 1, ] == 2), "\n")
    cat("shape1 =  ", dgamma_shape[t+1,,1],
        " ;shape2 =  ",  dgamma_shape[t+1,,2], "\n")
    cat("v3 =  ", v_matrix[(t + 1), 3],
        " ;w3 =  ",  w_matrix[(t + 1), 3], "\n")
    cat("shape3 =  ", dgamma_shape[t+1,,3],
        " ;rate3 =  ",  dgamma_rate[t+1,,3], "\n")
    par(mfrow = c(3, 3))
    plot(K_number[1:(t + 1)], type = "l", main = "Number of clusters")
    plot(v_matrix[1:(t + 1), 1], type = "l", main = "v_matrix[,1]")
    plot(v_matrix[1:(t + 1), 2], type = "l", main = "v_matrix[,2]")
    plot(dgamma_shape[1:(t + 1), 1, 1], type = "l", main = "Gamma shape [1,1]")
    plot(dgamma_rate[1:(t + 1), 1, 1],  type = "l", main = "Gamma rate [1,1]")
    plot(rho[1, 2, 1, 1:(t + 1)], type = "l", main = "rho(1,2) cluster=1")
    plot(rho[1, 2, 2, 1:(t + 1)], type = "l", main = "rho(1,2) cluster=2")
    plot(w_matrix[1:(t + 1), 1], type = "l", main = "w_matrix[,1]")
    plot(w_matrix[1:(t + 1), 2], type = "l", main = "w_matrix[,2]")
    
    cat("loop ", t, "\n")
    
  } # end of main MCMC loop
  
  #### 3) Return MCMC output ####
  return(list(
    v_matrix      = v_matrix,
    rho           = rho,
    w_matrix      = w_matrix,
    y_pred        = y_pred,
    K_number      = K_number,
    k_matrix      = k_matrix,
    dgamma_shape  = dgamma_shape,
    dgamma_rate   = dgamma_rate
  ))
}
