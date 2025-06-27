
source("copula_source_mix_distribution_gamma_new.R")


#######################################Sports application###########################

###################################################################################

library(readxl)
library(dplyr)

# Define the file paths
file1 <- "advplayer2324.xlsx"
file2 <- "advplayer2223.xlsx"
file3 <-"advplayer2122.xlsx"
# Read each file
df1 <- read_xlsx(file1)
df2 <- read_xlsx(file2)
df3 <- read_xlsx(file3)

# Combine them by stacking rows
merged_df <- bind_rows(df1, df2, df3)
#merged_df_23 <- df1
# Take a look at the first few rows
head(merged_df)


merged_df <- as.data.frame(merged_df)


#s_tcop2 <- MCMC_tcopula(N=2,data = df_pobs)
merged_df_233 <- merged_df[,c('Age','PER','TS%','ORB%','DRB%','AST%','BLK%','TOV%','USG%',
                                'OWS','DWS','OBPM','DBPM')]

merged_df_233 <- scale(na.omit(merged_df_233))

merged_df_233_z <- merged_df_233

merged_df_233[,1:ncol(merged_df_233)] <- sapply(1:ncol(merged_df_233),function(x) merged_df_233[,x] -min(merged_df_233[,x]) + 0.01)



library(ggplot2)
library(tidyverse)

# Sample data (uncomment and modify if you want to test quickly):
# df_normalized <- data.frame(
#   Points = rnorm(100, mean = 20, sd = 5),
#   Rebounds = rnorm(100, mean = 5, sd = 2),
#   Assists = rnorm(100, mean = 6, sd = 3)
# )
df_center <- as.data.frame(merged_df_233)
# 1. Reshape data from wide to long format
df_long <- df_center %>%
  pivot_longer(cols = everything(), 
               names_to = "Stat", 
               values_to = "Value")

# 2. Plot histograms + density curves in facets
ggplot(df_long, aes(x = Value)) +
  # Histogram (bins = 30 can be adjusted)
  geom_histogram(aes(y = ..density..), 
                 fill = "steelblue", 
                 color = "white", 
                 bins = 20, 
                 alpha = 0.6) +
  # Overlay density curve
  geom_density(color = "firebrick", 
               alpha = 0.5, 
               size = 1) +
  # Separate plot for each stat
  facet_wrap(~ Stat, scales = "free_x") +
  # Nice theme
  theme_minimal(base_size = 14) +
  labs(
    title = "Distribution of NBA Player Stats",
    x = "Value",
    y = "Density"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )


pkgs <- c("mclust", "mvtnorm", "sn",      # simulation & GMM fitting
          "clue", "mclust",              # metrics  (solve_LSAP & adjustedRand)
          "copula", "LaplacesDemon")     # required by your sampler
need <- setdiff(pkgs, rownames(installed.packages()))
if(length(need)) install.packages(need, dependencies = TRUE)
invisible(lapply(pkgs, require, character.only = TRUE))

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
  
  list(cl = cl_new,G=g)
}


gmmfit <- fit_gmm(merged_df_233_z)


################Copuola model runs#####################
t_cop_sports <- MCMC_tcopula(N=30000, data = merged_df_233)
#######################################################

G <- gmmfit$G

## ------------------------------------------------------------
## 1.  BIC value
## ------------------------------------------------------------
cat("Model BIC =", G$bic, "\n")   # scalar (overall) BIC

## ------------------------------------------------------------
## 2.  mixture weights, ordered largest → smallest
## ------------------------------------------------------------
w <- G$parameters$pro            # numeric vector
ord <- order(w, decreasing = TRUE)
w_sorted <- data.frame(
  Component = paste0("Cluster ", ord),
  Weight    = w[ord]
)
print(w_sorted, row.names = FALSE)

## ------------------------------------------------------------
## 3.  correlation matrix heat-map
##    (for VVV model, Σ_k differ → show the global covariance
##     OR one matrix per cluster; here we plot the pooled Σ̂)
## ------------------------------------------------------------

## ------------------------------------------------------------
## 0.  libraries
## ------------------------------------------------------------
library(corrplot)       # for quick correlation heat-maps
library(scales)         # for percent_format()

## ------------------------------------------------------------
## 1.  get weights & covariance array p×p×G
## ------------------------------------------------------------
w        <- G$parameters$pro
SigmaArr <- G$parameters$variance$sigma   # array (p × p × G)
G_hat    <- length(w)                     # number of clusters

## sort clusters by weight
ord <- order(w, decreasing = TRUE)
w   <- w[ord]
SigmaArr <- SigmaArr[,, ord, drop = FALSE]

## ------------------------------------------------------------
## 2.  loop over clusters, draw a plot each
## ------------------------------------------------------------
pal <- colorRampPalette(c("#2166ac","white","#b2182b"))(200)

## ------------------------------------------------------------
## 0.  libraries
## ------------------------------------------------------------
library(corrplot)       # for quick correlation heat-maps
library(scales)         # for percent_format()

## ------------------------------------------------------------
## 1.  get weights & covariance array p×p×G
## ------------------------------------------------------------
w        <- G$parameters$pro
SigmaArr <- G$parameters$variance$sigma   # array (p × p × G)
G_hat    <- length(w)                     # number of clusters

## sort clusters by weight
ord <- order(w, decreasing = TRUE)
w   <- w[ord]
SigmaArr <- SigmaArr[,, ord, drop = FALSE]

## ------------------------------------------------------------
## 2.  loop over clusters, draw a plot each
## ------------------------------------------------------------
## ------------------------------------------------------------
## 0.  libraries
## ------------------------------------------------------------
library(corrplot)       # for quick correlation heat-maps
library(scales)         # for percent_format()


## ------------------------------------------------------------
##  extract & order by weight
## ------------------------------------------------------------
w        <- G$parameters$pro
SigmaArr <- G$parameters$variance$sigma   # p × p × G array
ord      <- order(w, decreasing = TRUE)
w        <- w[ord]
SigmaArr <- SigmaArr[,, ord, drop = FALSE]
G_hat    <- length(w)

pal <- colorRampPalette(c("#2166ac","white","#b2182b"))(200)

## ------------------------------------------------------------
##  loop and write files
## ------------------------------------------------------------
out_dir <- "cluster_corrplots"
dir.create(out_dir, showWarnings = FALSE)

for (k in seq_len(G_hat)) {
  
  rho <- cov2cor(SigmaArr[,,k])
  weight_lbl <- percent_format(accuracy = .1)(w[k])
  
  ## ----------- PDF ------------
  pdf(file = file.path(out_dir,
                       sprintf("corr_cluster%02d.pdf", k)),
      width = 5, height = 5, family = "Times")
  corrplot(rho,
           method      = "color",
           col         = pal,
           addCoef.col = "black",
           tl.col      = "black",
           tl.srt      = 45,
           mar         = c(0,0,2,0),
           title       = sprintf("Cluster %d  (weight = %s)",
                                 k, weight_lbl))
  dev.off()
  
  ## ----------- PNG (optional) ------------
  png(file = file.path(out_dir,
                       sprintf("corr_cluster%02d.png", k)),
      width = 1800, height = 1800, res = 300)
  corrplot(rho,
           method      = "color",
           col         = pal,
           addCoef.col = "black",
           tl.col      = "black",
           tl.srt      = 45,
           mar         = c(0,0,2,0),
           title       = sprintf("Cluster %d  (weight = %s)",
                                 k, weight_lbl))
  dev.off()
}
