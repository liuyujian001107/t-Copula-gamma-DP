#######################################################################
##                 0.  PACKAGES + YOUR SAMPLER                       ##
#######################################################################
pkgs <- c("mclust", "mvtnorm", "sn",      # simulation & GMM fitting
          "clue", "mclust",              # metrics  (solve_LSAP & adjustedRand)
          "copula", "LaplacesDemon")     # required by your sampler
need <- setdiff(pkgs, rownames(installed.packages()))
if(length(need)) install.packages(need, dependencies = TRUE)
invisible(lapply(pkgs, require, character.only = TRUE))

## ---- source your long routine (name exactly as you asked) ----------
source("copula_source_mix_distribution_gamma_new.R") # brings MCMC_tcopula()

#######################################################################
##                 1.  HELPERS                                        ##
#######################################################################
shift_positive <- function(X){
  sweep(X, 2, apply(X, 2, min), FUN = "-") + 0.01
}

## ---------- (fixed) performance from confusion matrix --------------
#######################################################################
##  1.  safest Hungarian-based macro metrics                          ##
#######################################################################
best_map_accuracy_f1 <- function(tab) {
  if (nrow(tab) == 0 || ncol(tab) == 0)
    return(c(Accuracy = NA, F1 = NA))
  
  ## --- shortcut: same labels, same order ----------------------------
  if (nrow(tab) == ncol(tab) &&
      all(sort(rownames(tab)) == sort(colnames(tab)))) {
    acc <- sum(diag(tab)) / sum(tab)
    
    prec <- diag(tab) / pmax(colSums(tab), 1)
    rec  <- diag(tab) / pmax(rowSums(tab), 1)
    f1   <- ifelse(prec + rec == 0, 0, 2 * prec * rec / (prec + rec))
    
    return(c(Accuracy = acc, F1 = mean(f1)))
  }
  
  ## --- otherwise: pad to square then Hungarian ----------------------
  r <- nrow(tab); c <- ncol(tab)
  if (r > c) {
    tab <- cbind(tab, matrix(0, r, r - c))   # add empty columns
  } else if (c > r) {
    tab <- rbind(tab, matrix(0, c - r, c))   # add empty rows
  }
  
  perm <- solve_LSAP(tab, maximum = TRUE)
  tab2 <- tab[, perm]
  
  acc <- sum(diag(tab2)) / sum(tab2)
  
  prec <- diag(tab2) / pmax(colSums(tab2), 1)
  rec  <- diag(tab2) / pmax(rowSums(tab2), 1)
  f1   <- ifelse(prec + rec == 0, 0, 2 * prec * rec / (prec + rec))
  
  c(Accuracy = acc, F1 = mean(f1))
}

#######################################################################
##  2.  one-call wrapper  (unchanged signature)                       ##
#######################################################################
metrics <- function(truth, pred) {
  tab <- table(truth, pred)
  c(ARI = adjustedRandIndex(truth, pred),
    best_map_accuracy_f1(tab))
}


#######################################################################
##               2.  2-D & 3-D *true* parameter sets                  ##
#######################################################################
true_params <- function(d = 2, K = 3){
  stopifnot(d %in% c(2,3), K %in% 2:4)
  
  ## --- weights ------------------------------------------------------
  if(K==2) w <- c(2/3, 1/3)
  if(K==3) w <- c(1/2, 1/3, 1/6)
  if(K==4) w <- c(.40, .25, .20, .15)  
  
  ## --- correlation matrices  ---------------------------------------
  if(d==2){
    S1 <- matrix(c( 1,  0.7,
                    0.7, 1), 2)
    S2 <- matrix(c( 1, -0.4,
                    -0.4, 1), 2)
    S3 <- matrix(c( 1, -0.8,
                    -0.8,  1), 2)
    S4 <- diag(2)                       # extra comp. for K = 4
    Sig <- list(S1,S2,S3,S4)[1:K]
  } else { # d == 3
    S1 <- matrix(c( 1,  0.7,  0.5,
                    0.7, 1,    0.6,
                    0.5, 0.6,  1), 3, byrow = TRUE)
    S2 <- matrix(c( 1, -0.4, -0.2,
                    -0.4,  1,  -0.1,
                    -0.2, -0.1,  1), 3, byrow = TRUE)
    S3 <- matrix(c( 1, -0.8, -0.6,
                    -0.8,  1,   0.3,
                    -0.6,  0.3, 1), 3, byrow = TRUE)
    S4 <- diag(3)
    Sig <- list(S1,S2,S3,S4)[1:K]
  }
  
  ## --- degrees of freedom & Gamma marginals -------------------------
  nu_vec <- c(8,6,4,10)[1:K]       # df for t / skew-t
  shape  <- c(5,2,3,4)[1:K]        # Γ shapes (same for every dim per comp)
  rate   <- c(2,5,3,4)[1:K]        # Γ rates
  
  list(w = w, Sig = Sig, nu = nu_vec,
       shape = shape, rate = rate)
}

#######################################################################
##               3.  DATA GENERATORS                                  ##
#######################################################################
## -- helper: sample from t-copula + Gamma marginals ------------------
sample_tcop_gamma <- function(n, d, Sigma, df, shape, rate){
  tc <- tCopula(param = P2p(Sigma), dim = d,
                dispstr = "un", df = df)
  U  <- rCopula(n, tc)
  X  <- sapply(seq_len(d),
               function(j) qgamma(U[,j], shape = shape, rate = rate))
  X
}

sim_tc_gamma_mix <- function(n, d, K){
  par <- true_params(d,K)
  labs <- sample(seq_len(K), n, TRUE, par$w)
  X <- matrix(NA, n, d)
  for(k in seq_len(K)){
    idx <- which(labs==k)
    if(length(idx)){
      X[idx,] <- sample_tcop_gamma(length(idx), d,
                                   par$Sig[[k]], par$nu[k],
                                   par$shape[k], par$rate[k])
    }
  }
  list(x = X, y = labs)
}

## GMM, t-mix and skew-t mix  (unchanged except parameters come from screenshot)
sim_GMM <- function(n,d,K){
  par <- true_params(d,K)
  labs <- sample(seq_len(K), n, TRUE, par$w)
  mus  <- lapply(1:K, function(k) rep( k*3, d))
  X <- matrix(NA,n,d)
  for(k in seq_len(K)){
    nk <- sum(labs==k)
    if(nk) X[labs==k,] <- mvtnorm::rmvnorm(nk, mus[[k]], par$Sig[[k]])
  }
  list(x=X,y=labs)
}

sim_tmix <- function(n,d,K){
  par <- true_params(d,K)
  labs <- sample(seq_len(K), n, TRUE, par$w)
  mus  <- lapply(1:K, function(k) rep( k*3, d))
  X <- matrix(NA,n,d)
  for(k in seq_len(K)){
    nk <- sum(labs==k)
    if(nk) X[labs==k,] <- mvtnorm::rmvt(nk, sigma = par$Sig[[k]],
                                        df = par$nu[k], delta = mus[[k]])
  }
  list(x=X,y=labs)
}

sim_skewtmix <- function(n, d, K){
  par <- true_params(d, K)                       # weights, Sigmas, ν, …
  labs <- sample(seq_len(K), n, TRUE, par$w)
  
  mus    <- lapply(1:K, function(k) rep(3*k, d))           # centroids
  alphas <- lapply(1:K, function(k) rep((-1)^(k+1) * 3, d)) # ± skewness
  
  X <- matrix(NA, n, d)
  for(k in seq_len(K)){
    nk <- sum(labs == k)
    if(nk){
      X[labs == k, ] <- sn::rmst(nk,
                                 xi    = mus[[k]],
                                 Omega = par$Sig[[k]],
                                 alpha = alphas[[k]],
                                 nu    = par$nu[k])        # *** skew-t ***
    }
  }
  list(x = X, y = labs)
}

########################################################################
##               4.  FITTING ROUTINES                                  ##
########################################################################
fit_gmm <- function(X) {
  ## run unconstrained Gaussian mixture (VVV = full Σ_k)
  g <- Mclust(X, modelNames = "VVV")
  
  ## mixture weights  (π_k, k = 1…G)  returned by mclust
  w  <- g$parameters$pro              # numeric vector length G
  G  <- length(w)                     # number of components actually fitted
  
  ## order component indices by decreasing weight
  ## e.g. if weights = c(0.60, 0.25, 0.15)  →  ord = c(1,2,3)
  ##      if weights = c(0.10, 0.70, 0.20)  →  ord = c(2,3,1)
  ord <- order(w, decreasing = TRUE, na.last = NA)
  
  ## build a mapping old_label -> rank(weight)
  ##   old 1,2,3  (example above)  →  new 2,3,1
  map <- match(seq_len(G), ord)      # match() returns the rank positions
  
  ## re-label every observation
  cl_new <- map[g$classification]
  
  list(cl = cl_new)
}

################################################################################
##  REPLACEMENT:  use the *last 5 000* post-burn-in draws                      ##
################################################################################
fit_tcopula_gamma <- function(X,
                              component = 10,
                              n_iter    = 20000,   # <- long chain
                              burn      = 18000) { # <- leaves 5 000 rows
  #stopifnot(n_iter - burn >= 5000)
  
  Xs  <- shift_positive(X)
  out <- MCMC_tcopula(N = n_iter,
                      component   = component,
                      data_source = Xs)
  
  ## ---------- majority vote over the last 5 000 iterations ------------
  keep <- (burn + 1):n_iter                      # rows to use
  Kmat <- out$k_matrix[keep, , drop = FALSE]     # (5000 × n) integer matrix
  
  ## mode for each observation (ties broken at random)
  mode_of_vec <- function(v){
    tbl <- table(v); as.integer(sample(names(tbl[tbl == max(tbl)]), 1))
  }
  cl <- apply(Kmat, 2, mode_of_vec)              # length-n vector
  
  list(cl = cl)
}

########################################################################
##               5.  ONE REPLICATE                                     ##
########################################################################
run_once <- function(datagen = c("gmm","tmix","skewt","tcgamma"),
                     n = 1000, d = 2, K = 3){
  datagen <- match.arg(datagen)
  dd <- switch(datagen,
               gmm     = sim_GMM     (n,d,K),
               tmix    = sim_tmix    (n,d,K),
               skewt   = sim_skewtmix(n,d,K),
               tcgamma = sim_tc_gamma_mix(n,d,K))
  X <- dd$x; y <- dd$y
  
  g_res  <- fit_gmm(X)
  print(g_res)
  tc_res <- fit_tcopula_gamma(X,
                              component = 10,
                              n_iter    = 20000,   # <-- 20 000 total
                              burn      = 18000)   # <-- 15 000 burn-in
  
  rbind(GMM         = metrics(y, g_res$cl),
        TC_Gamma_DP = metrics(y, tc_res$cl))
}

########################################################################
##               6.  FULL EXPERIMENT                                   ##
########################################################################
experiment <- function(datagen, d, K, n=1000, reps=10){
  arr <- replicate(reps,
                   run_once(datagen,n,d,K),
                   simplify = FALSE)
  A <- array(unlist(arr), dim=c(2,3,reps),
             dimnames=list(Method=c("GMM","TC_Gamma_DP"),
                           Metric=c("ARI","Accuracy","F1"),
                           Rep=NULL))
  round(apply(A,1:2,mean),3)
}

########################################################################
##               7.  RUN EVERYTHING                                    ##
########################################################################
#set.seed(2025)
settings <- expand.grid(Data  = c("gmm","skewt","tcgamma"),
                       Dim   = c(2,3),
                       K     = c(2,3,4),
                       N     = c(500,1000,2000),         # << NEW COLUMN
                       stringsAsFactors = FALSE)


run_setting <- function(set, reps = 5){
  cat("\n>>>  Data:",set$Data,
      " d=",set$Dim,
      " K=",set$K,
      " n=",set$N,"\n")
  experiment(datagen = set$Data,
             d       = set$Dim,
             K       = set$K,
             n       = set$N,
             reps    = reps)             # ← keep reps small for a quick test
}
## ---- iterate over the grid -----------------------------------------------
###############################################################################
##  A.  EXTRA PACKAGES + 100-core parallel plan                               ##
###############################################################################
pkgs_extra <- c("future.apply", "data.table")
need <- setdiff(pkgs_extra, rownames(installed.packages()))
if(length(need)) install.packages(need, dependencies = TRUE)
library(future.apply)
library(data.table)

plan(multisession, workers = min(120, future::availableCores()-2))
message("Using ", nbrOfWorkers(), " parallel workers")

dir.create("data_dumps",  showWarnings = FALSE)
dir.create("metrics_raw", showWarnings = FALSE)

###############################################################################
##  B.  SINGLE REPLICATION  (also dumps the data set)                         ##
###############################################################################
run_one_and_save <- function(model, n, d, K, rep_id, tag){
  dat <- switch(model,
                gmm   = sim_GMM     (n,d,K),
                tmix  = sim_tmix    (n,d,K),
                skewt = sim_skewtmix(n,d,K))
  ## -- record the synthetic sample -----------------------------------------
  saveRDS(dat,
          file = sprintf("data_dumps/%s_rep%03d.rds", tag, rep_id),
          compress = "xz")
  
  ## -- fit the two models ---------------------------------------------------
  g_res <- fit_gmm(dat$x)
  t_res <- fit_tcopula_gamma(dat$x,
                             component = 20,
                             n_iter    = 20000,
                             burn      = 18000)
  
  ## -- three metrics --------------------------------------------------------
  rbind(GMM         = metrics(dat$y, g_res$cl),
        TC_Gamma_DP = metrics(dat$y, t_res$cl))
}

###############################################################################
##  C.  100 REPETITIONS IN PARALLEL FOR ONE SCENARIO                          ##
###############################################################################
run_100 <- function(model, d, K, n, reps = 100){
  tag <- sprintf("%s_d%d_K%d_n%d", model, d, K, n)
  
  mats <- future_lapply(seq_len(reps), function(r)
    run_one_and_save(model, n, d, K, r, tag))
  
  A <- array(unlist(mats), dim = c(2, 3, reps),
             dimnames = list(Method = c("GMM", "TC_Gamma_DP"),
                             Metric = c("ARI", "Accuracy", "F1"),
                             Rep    = NULL))
  
  saveRDS(A, file = sprintf("metrics_raw/%s_metrics.rds", tag))
  
  list(tag  = tag,
       mean = apply(A, 1:2, mean),
       se   = apply(A, 1:2, sd) / sqrt(reps))
}

###############################################################################
##  D.  MONTE-CARLO GRID (only gmm / tmix / skewt)                            ##
###############################################################################


###############################################################################
## 1.  build (scenario, rep) job table                                        ##
###############################################################################
grid <- expand.grid(Data = c("gmm","skewt","tcgamma"),   # no "tcgamma"
                    Dim  = c(2,3),
                    K    = c(2,3),
                    N    = c(1000,2000,3000),
                    Rep  = 1:50,                     # 50 repetitions
                    stringsAsFactors = FALSE)

###############################################################################
## 2.  single job: generate data, fit, dump, return metrics                  ##
###############################################################################
job_fun <- function(s){
  
  tag <- sprintf("%s_d%d_K%d_n%d_rep%03d",
                 s$Data, s$Dim, s$K, s$N, s$Rep)
  
  ## ------ simulate ---------------------------------------------------------
  dat <- switch(s$Data,
                gmm   = sim_GMM     (s$N, s$Dim, s$K),
                tmix  = sim_tmix    (s$N, s$Dim, s$K),
                skewt = sim_skewtmix(s$N, s$Dim, s$K),
                tcgamma = sim_tc_gamma_mix(s$N, s$Dim, s$K))
  
  saveRDS(dat, file = file.path("data_dumps", paste0(tag, ".rds")),
          compress = "xz")
  
  ## ------ fit models -------------------------------------------------------
  g_res <- fit_gmm(dat$x)
  t_res <- fit_tcopula_gamma(dat$x,
                             component = 20,
                             n_iter    = 20000,
                             burn      = 18000)
  
  ## ------ metrics ----------------------------------------------------------
  cbind(
    Scenario = sprintf("%s_d%d_K%d_n%d", s$Data, s$Dim, s$K, s$N),
    Rep      = s$Rep,
    Method   = c('GMM','TC_Gamma_DP'),
    ARI  = c(metrics(dat$y, g_res$cl)[1],metrics(dat$y, t_res$cl)[1]),
    Accuracy = c(metrics(dat$y, g_res$cl)[2],metrics(dat$y, t_res$cl)[2]),
    F1 =  c(metrics(dat$y, g_res$cl)[3],metrics(dat$y, t_res$cl)[3])
  )
}

###############################################################################
## 3.  run all jobs in parallel (100 cores fully utilised)                    ##
###############################################################################
res_list <- future_lapply(split(grid, seq_len(nrow(grid))), job_fun)
metrics_raw <- do.call(rbind, res_list)         #  (5400 × 5) matrix

saveRDS(metrics_raw, "metrics_raw/all_metrics2.rds")
fwrite(metrics_raw, "metrics_raw2.csv")

###############################################################################
## 4.  summary: mean ± SE per scenario & method                               ##
###############################################################################
library(data.table)
DT <- as.data.table(metrics_raw)


DT[, `:=`(ARI       = as.numeric(ARI),
          Accuracy  = as.numeric(Accuracy),
          F1        = as.numeric(F1))]

summary_dt <- DT[, .(
  Mean_ARI      = mean(ARI),
  SE_ARI        = sd(ARI)/sqrt(.N),
  Mean_Accuracy = mean(Accuracy),
  SE_Accuracy   = sd(Accuracy)/sqrt(.N),
  Mean_F1       = mean(F1),
  SE_F1         = sd(F1)/sqrt(.N)
), by = .(Scenario, Method)]

fwrite(summary_dt, "results_summary2.csv")
print(summary_dt[])



