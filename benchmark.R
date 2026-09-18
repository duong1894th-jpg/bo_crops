library(dplyr)
library(parallel)

dir.create("./results", showWarnings = FALSE)

# ---------------------------------------------------------
# CORE MATH & UTILITIES
# ---------------------------------------------------------
minMax <- function(x) {
  if (max(x, na.rm = TRUE) == min(x, na.rm = TRUE)) return(rep(0, length(x)))
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

rbf_kernel <- function(train, test, l = 1) {
  train <- as.matrix(train); test <- as.matrix(test)
  sum_sq_train <- rowSums(train^2); sum_sq_test <- rowSums(test^2)
  dist_sq <- outer(sum_sq_train, sum_sq_test, "+") - 2 * (train %*% t(test))
  exp(-0.5 * pmax(dist_sq, 0) / (l^2))
}

compute_joint_kernel <- function(X1, X2, C1, C2, l_x, l_c, sigma) {
  (sigma^2) * (rbf_kernel(X1, X2, l_x) * rbf_kernel(C1, C2, l_c))
}

get_ei_max <- function(mu, sd, y_best, alpha) {
  imp <- mu - (1 + alpha) * y_best
  sd_safe <- pmax(sd, 1e-9)
  Z <- imp / sd_safe
  ei  <- imp * pnorm(Z) + sd * dnorm(Z)
  ei[sd == 0] <- 0
  return(ei)
}

get_expected_loss <- function(mu, sd, xi_n, remaining_budget) {
  sd_safe <- pmax(sd, 1e-9)
  Z <- (xi_n - mu) / sd_safe
  loss <- ((xi_n - mu) * pnorm(Z) + sd * dnorm(Z)) / remaining_budget
  loss[sd == 0] <- 0
  return(loss)
}

# ---------------------------------------------------------
# MODEL 1: IVAN'S HETEROSKEDASTIC BO
# ---------------------------------------------------------
neg_log_likelihood_ivan <- function(params, X_train, C_train, y) {
  l_x <- exp(params[1]); l_c <- exp(params[2]); sigma <- exp(params[3])
  alpha_noise <- params[4]
  nx <- ncol(X_train); nc <- ncol(C_train)
  beta_X <- params[5:(4+nx)]; beta_C <- params[(5+nx):(4+nx+nc)]
  
  K <- compute_joint_kernel(X_train, X_train, C_train, C_train, l_x, l_c, sigma)
  noise_log <- alpha_noise + as.matrix(X_train) %*% beta_X + as.matrix(C_train) %*% beta_C
  Ky <- K + diag(as.numeric(exp(noise_log)) + 1e-6, nrow(X_train))
  
  L <- tryCatch(chol(Ky), error = function(e) NULL)
  if (is.null(L)) return(Inf)
  alpha_vec <- backsolve(L, forwardsolve(t(L), y))
  return(as.numeric(0.5 * t(y) %*% alpha_vec + sum(log(diag(L))) + (nrow(X_train) / 2) * log(2 * pi)))
}

GP_ivan <- function(X_train, X_test, C_train, C_test, params, Y_train) {
  l_x <- exp(params[1]); l_c <- exp(params[2]); sigma <- exp(params[3])
  alpha_noise <- params[4]
  nx <- ncol(X_train); nc <- ncol(C_train)
  beta_X <- params[5:(4+nx)]; beta_C <- params[(5+nx):(4+nx+nc)]
  
  noise_log_tr <- alpha_noise + as.matrix(X_train) %*% beta_X + as.matrix(C_train) %*% beta_C
  K_tr <- compute_joint_kernel(X_train, X_train, C_train, C_train, l_x, l_c, sigma) + 
          diag(as.numeric(exp(noise_log_tr)) + 1e-6, nrow(X_train))
  Ks <- compute_joint_kernel(X_train, X_test, C_train, C_test, l_x, l_c, sigma)
  Kss <- compute_joint_kernel(X_test, X_test, C_test, C_test, l_x, l_c, sigma)
  
  L <- tryCatch(chol(K_tr), error=function(e) NULL)
  if (is.null(L)) return(list(mu=rep(0,nrow(X_test)), sd=rep(Inf,nrow(X_test))))
  
  alpha_vec <- backsolve(L, forwardsolve(t(L), Y_train))
  mu <- as.numeric(t(Ks) %*% alpha_vec)
  v <- forwardsolve(t(L), Ks)
  Cov <- Kss - t(v) %*% v
  return(list(mu = mu, sd = sqrt(pmax(diag(Cov), 0))))
}

# ---------------------------------------------------------
# MODELS 2 & 3: STANDARD BO (CBO & SBO)
# ---------------------------------------------------------
neg_log_likelihood_cgp <- function(params, X_train, C_train, y) {
  l_x <- exp(params[1]); l_c <- exp(params[2]); sigma <- exp(params[3]); noise_var <- exp(params[4])
  Ky <- compute_joint_kernel(X_train, X_train, C_train, C_train, l_x, l_c, sigma) + 
        diag(noise_var + 1e-6, nrow(X_train))
  L <- tryCatch(chol(Ky), error = function(e) NULL)
  if (is.null(L)) return(Inf)
  alpha <- backsolve(L, forwardsolve(t(L), y))
  return(as.numeric(0.5 * t(y) %*% alpha + sum(log(diag(L))) + (nrow(X_train) / 2) * log(2 * pi)))
}

GP <- function(X1, X2, C1, C2, l_x, l_c, Y, sigma, noise_var) {
  n_train <- nrow(as.matrix(X1))
  K_noise <- compute_joint_kernel(X1, X1, C1, C1, l_x, l_c, sigma) + diag(noise_var + 1e-6, n_train)
  Ks <- compute_joint_kernel(X1, X2, C1, C2, l_x, l_c, sigma)
  Kss <- compute_joint_kernel(X2, X2, C2, C2, l_x, l_c, sigma)
  
  L <- tryCatch(chol(K_noise), error = function(e) NULL)
  if (is.null(L)) return(list(mu = rep(0, nrow(as.matrix(X2))), sd = rep(Inf, nrow(as.matrix(X2)))))
  
  alpha <- backsolve(L, forwardsolve(t(L), Y))
  mu <- as.numeric(t(Ks) %*% alpha)
  v <- forwardsolve(t(L), Ks)
  Cov <- Kss - t(v) %*% v
  return(list(mu = mu, sd = sqrt(pmax(diag(Cov), 0))))
}

# ---------------------------------------------------------
# SIMULATION LOGIC
# ---------------------------------------------------------
cat("Loading Data...\n")
nutrient_full <- read.csv('Soil Nutrients.csv')

N_budget <- 20
batch_size <- 5
seeds <- 1:100

run_seed_for_crop <- function(seed, X_bo, C_bo, Y_bo, actual_max) {
  set.seed(seed)
  init_indices <- sample(1:nrow(X_bo), 5)
  X_init <- X_bo[init_indices, ]; C_init <- C_bo[init_indices, ]; Y_init <- as.numeric(Y_bo[init_indices])
  X_pool <- X_bo[-init_indices, ]; C_pool <- C_bo[-init_indices, ]; Y_pool <- as.numeric(Y_bo[-init_indices])
  C_dummy <- matrix(0, nrow=nrow(C_init), ncol=ncol(C_init))
  C_pool_dummy <- matrix(0, nrow=nrow(C_pool), ncol=ncol(C_pool))
  
  run_model <- function(model_type) {
    X_tr <- X_init; C_tr <- if(model_type=="SBO") C_dummy else C_init; Y_tr <- Y_init
    X_te <- X_pool; C_te <- if(model_type=="SBO") C_pool_dummy else C_pool; Y_te <- Y_pool
    history_batch_max <- c()
    epsilon <- 0.02
    
    while(nrow(X_tr) < N_budget && max(Y_tr) < actual_max && length(Y_te) > 0) {
      if (model_type == "RANDOM") {
        sel_idx <- sample(1:nrow(X_te), min(batch_size, nrow(X_te)))
        batch_y <- Y_te[sel_idx]
        X_tr <- rbind(X_tr, X_te[sel_idx,]); Y_tr <- c(Y_tr, batch_y)
        X_te <- X_te[-sel_idx,]; Y_te <- Y_te[-sel_idx]
        history_batch_max <- c(history_batch_max, max(batch_y))
        next
      }
      
      best_val <- Inf; best_par <- NULL
      for (r in 1:3) {
        if (model_type == "IVAN") {
          init <- c(log(runif(1,0.05,2)), log(runif(1,0.001,1)), log(0.01), rep(0, ncol(X_tr)), rep(0, ncol(C_tr)))
          res <- optim(init, neg_log_likelihood_ivan, X_train=X_tr, C_train=C_tr, y=Y_tr, method="BFGS")
        } else {
          init <- c(log(runif(1,0.05,2)), log(runif(1,0.05,2)), log(runif(1,0.001,1)), log(0.01))
          res <- optim(init, neg_log_likelihood_cgp, X_train=X_tr, C_train=C_tr, y=Y_tr, method="BFGS")
        }
        if (res$value < best_val) { best_val <- res$value; best_par <- res$par }
      }
      
      if (model_type == "IVAN") {
        pred <- GP_ivan(X_tr, X_te, C_tr, C_te, best_par, Y_tr)
        pred_tr <- GP_ivan(X_tr, X_tr, C_tr, C_tr, best_par, Y_tr)
      } else {
        pred <- GP(X_tr, X_te, C_tr, C_te, exp(best_par[1]), exp(best_par[2]), Y_tr, exp(best_par[3]), exp(best_par[4]))
        pred_tr <- GP(X_tr, X_tr, C_tr, C_tr, exp(best_par[1]), exp(best_par[2]), Y_tr, exp(best_par[3]), exp(best_par[4]))
      }
      
      xi_n <- max(pred_tr$mu)
      ei_pool <- get_ei_max(pred$mu, pred$sd, xi_n, 0)
      
      if (model_type != "IVAN") {
        rem_budget <- max(1, N_budget - nrow(X_tr))
        loss_pool <- get_expected_loss(pred$mu, pred$sd, xi_n, rem_budget)
        valid_mask <- (ei_pool >= loss_pool)
        if (any(valid_mask)) ei_pool[!valid_mask] <- -Inf
      }
      
      batch_X <- NULL; batch_C <- NULL; batch_y <- c()
      X_p <- X_te; C_p <- C_te; y_p <- Y_te; ei_p <- ei_pool
      curr_eps <- epsilon
      
      l_x <- exp(best_par[1]); l_c <- if(model_type=="IVAN") exp(best_par[2]) else exp(best_par[2]); sig <- if(model_type=="IVAN") exp(best_par[3]) else exp(best_par[3])
      
      for (b in 1:min(batch_size, nrow(X_p))) {
        ord <- order(ei_p, decreasing = TRUE)
        sel_idx <- NULL
        for (idx in ord) {
          cx <- X_p[idx,,drop=F]; cc <- C_p[idx,,drop=F]
          if (is.null(batch_X)) { sel_idx <- idx; break }
          dists <- sapply(1:nrow(batch_X), function(j) compute_joint_kernel(cx, batch_X[j,,drop=F], cc, batch_C[j,,drop=F], l_x, l_c, sig))
          if (sum(dists) < curr_eps) { sel_idx <- idx; break }
        }
        if (is.null(sel_idx)) { curr_eps <- curr_eps * 2; next }
        batch_X <- rbind(batch_X, X_p[sel_idx,]); batch_C <- rbind(batch_C, C_p[sel_idx,]); batch_y <- c(batch_y, y_p[sel_idx])
        X_p <- X_p[-sel_idx,]; C_p <- C_p[-sel_idx,]; y_p <- y_p[-sel_idx]; ei_p <- ei_p[-sel_idx]
      }
      
      X_tr <- rbind(X_tr, batch_X); C_tr <- rbind(C_tr, batch_C); Y_tr <- c(Y_tr, batch_y)
      X_te <- X_p; C_te <- C_p; Y_te <- y_p
      history_batch_max <- c(history_batch_max, max(batch_y))
    }
    
    # Pad history if stopped early
    expected_rounds <- ceiling((N_budget - 5) / batch_size)
    if (length(history_batch_max) < expected_rounds) {
      pad_val <- if(length(history_batch_max) > 0) history_batch_max[length(history_batch_max)] else actual_max
      history_batch_max <- c(history_batch_max, rep(pad_val, expected_rounds - length(history_batch_max)))
    }
    
    return(actual_max - cummax(history_batch_max))
  }
  
  ivan_res <- run_model("IVAN")
  cbo_res <- run_model("CBO")
  sbo_res <- run_model("SBO")
  rand_res <- run_model("RANDOM")
  
  return(list(IVAN = ivan_res, CBO = cbo_res, SBO = sbo_res, RANDOM = rand_res))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  crops <- args
} else {
  crops <- c("Lettuce", "Strawberry", "Spinach", "Asparagus", "Cabbage")
}


for (crop in crops) {
  cat(sprintf("\n=========================================\n"))
  cat(sprintf("Running Simulations for: %s\n", crop))
  cat(sprintf("=========================================\n"))
  
  crop_data <- nutrient_full %>% filter(Name == crop)
  if(nrow(crop_data) == 0) {
    cat("No data found for", crop, "skipping...\n")
    next
  }
  
  crop_data <- crop_data[order(crop_data$Yield), ]
  crop_data$Yield <- minMax(crop_data$Yield)
  crop_data$Nitrogen <- minMax(crop_data$Nitrogen)
  crop_data$Phosphorus <- minMax(crop_data$Phosphorus)
  crop_data$Potassium <- minMax(crop_data$Potassium)
  
  Y_bo <- crop_data$Yield
  X_bo <- crop_data %>% select(Nitrogen, Phosphorus, Potassium)
  C_bo <- crop_data %>% select(Temperature, Rainfall, pH, Light_Hours, Light_Intensity, Rh)
  C_bo[] <- lapply(C_bo, minMax)
  actual_max <- max(Y_bo)
  
  results_list <- mclapply(seeds, function(s) run_seed_for_crop(s, X_bo, C_bo, Y_bo, actual_max), mc.cores = detectCores() - 1)
  
  # Aggregate
  agg_ivan <- do.call(rbind, lapply(results_list, function(r) r$IVAN))
  agg_cbo <- do.call(rbind, lapply(results_list, function(r) r$CBO))
  agg_sbo <- do.call(rbind, lapply(results_list, function(r) r$SBO))
  agg_rand <- do.call(rbind, lapply(results_list, function(r) r$RANDOM))
  
  write.table(agg_ivan, sprintf("./results/%s_ivan.csv", crop), row.names=FALSE, col.names=FALSE, sep=",")
  write.table(agg_cbo, sprintf("./results/%s_cbo.csv", crop), row.names=FALSE, col.names=FALSE, sep=",")
  write.table(agg_sbo, sprintf("./results/%s_sbo.csv", crop), row.names=FALSE, col.names=FALSE, sep=",")
  write.table(agg_rand, sprintf("./results/%s_rand.csv", crop), row.names=FALSE, col.names=FALSE, sep=",")
  
  cat(sprintf("Completed %s! Results saved.\n", crop))
}

cat("\nALL CROPS COMPLETED SUCCESSFULLY.\n")
