MissingFCDComponentsG <- function(B, M, Sigma, mean_0, sd_0) {
  var_names <- rownames(M)
  covariate_names <- colnames(M)
  G <- dim(B)[3]
  
  FCDComponents <- vector("list", length(var_names))
  names(FCDComponents) <- var_names
  
  for (i in seq_along(var_names)) {
    target_var <- var_names[i]
    
    # Parents in fixed model
    parents_i <- covariate_names[which(M[i, ] == 1)]
    
    # Initialize FCD list: time 0 to T
    FCDComponents[[target_var]] <- vector("list", G)
    names(FCDComponents[[target_var]]) <- as.character(1:G)
    
    ## ---- Time 1 to T ----
    for (g in 1:G) {
      B_g <- B[, , g]  #B[, , G]
      
      if (target_var %in% covariate_names) {
        children_vars <- var_names[which(M[, target_var] == 1)]
        
        model_children_g <- lapply(children_vars, function(child_var) {
          child_parents <- covariate_names[which(M[child_var, ] == 1)]
          other_parents <- setdiff(child_parents, target_var)
          
          list(
            child = child_var,
            beta_i = B[child_var, target_var, g],
            beta_other = B[child_var, other_parents, g],
            sigma2_j = Sigma[child_var],
            other_covariate_names = other_parents
          )
        })
      } else {
        model_children_g <- list()
      }
      
      FCDComponents[[target_var]][[as.character(g)]] <- list(
        sigma2_i = Sigma[target_var],
        parents_i = parents_i,
        beta_parents_i = B_g[target_var, parents_i],
        model_children = model_children_g
      )
    }
  }
  
  return(FCDComponents)
}

MissingValuesUpdateG <- function(data, missing_mask, structure, group_segmentation_list, sample_id = 'Sample_Id', time_id = 'Time', mean_0, sd_0) {
  samples <- unique(data[[sample_id]])
  times <- unique(data[[time_id]])
  T_max <- max(times)
  #print(samples)
  #print(times)
  #print(T_max)
  # Get variable names for columns excluding sample_id and time_id
  var_names <- setdiff(colnames(data), c(sample_id, time_id))
  
  for (s in samples) {
    group_idx <- as.character(group_segmentation_list[[as.character(s)]])
    #print(group_idx)
    for (t in 0:T_max) {
      time_window <- c(t - 1, t, t + 1)
      valid_times <- time_window[time_window >= 0 & time_window <= T_max]
      
      sub_data <- data[data[[sample_id]] == s & data[[time_id]] %in% valid_times, ]
      if (nrow(sub_data) < length(valid_times)) next
      
      X <- as.matrix(sub_data[, var_names, drop = FALSE])
      X <- X[match(valid_times, sub_data[[time_id]]), , drop = FALSE]
      if (t == 0) {
        X <- rbind(rep(1, ncol(X)), X)
      }
      mask_row <- which(data[[sample_id]] == s & data[[time_id]] == t)
      missing_row <- missing_mask[mask_row, var_names, drop = FALSE]
      missing_vars <- var_names[as.logical(missing_row)]
      if (length(missing_vars) == 0) next
      
      for (var in missing_vars) {
        if (t == T_max){
          struc <- structure[[var]][[group_idx]]
          struc$model_children <- list()
          sampled_val <- MissingFCDSample(X, struc)
        } 
        else if (t == 0){
          struc <- structure[[var]][[group_idx]]
          struc$sigma2_i <- sd_0[var]
          struc$parents_i <- "1"
          struc$beta_parents_i <- mean_0[var]
          sampled_val <- MissingFCDSample(X, struc)
        } 
        else{
          sampled_val <- MissingFCDSample(X, structure[[var]][[group_idx]])
        }
        X[which(valid_times == t), var] <- sampled_val
        data[mask_row, var] <- sampled_val
      }
    }
  }
  
  data
}


MissingFCDSample <- function(data, structure_i, prior_precision = 1) {
  # data: matrix with rows = times (t-1, t, t+1), cols = variable names (no intercept col)
  sigma2_i <- structure_i$sigma2_i
  
  sum_precisions <- 1 / sigma2_i + prior_precision #prior precision
  weighted_sum <- sum(data[1, structure_i$parents_i] * structure_i$beta_parents_i) / sigma2_i
  
  for (child in structure_i$model_children) {
    j <- child$child
    mu_minus_i_j <- sum(data[2, child$other_covariate_names] * child$beta_other)  # time t
    x_j_tp1 <- data[3, j]  # time t+1
    
    beta_i <- child$beta_i
    sigma2_j <- child$sigma2_j
    
    sum_precisions <- sum_precisions + (beta_i^2) / sigma2_j
    weighted_sum <- weighted_sum + (beta_i / sigma2_j) * (x_j_tp1 - mu_minus_i_j)
  }
  sigma2_post <- 1 / sum_precisions
  rnorm(1, mean = sigma2_post * weighted_sum, sd = sqrt(sigma2_post))
}

InitialMissingImputation <- function(df, s = 5, sample_id = "Sample_Id", time_id = "Time"){
  df_imputed <- df
  
  # Get variable columns (all except ID columns)
  variable_cols <- setdiff(names(df), c(sample_id, time_id))
  
  # Extract time vector
  time_points <- sort(unique(df[[time_id]]))
  
  for (t in time_points) {
    for (var in variable_cols) {
      # Identify rows at time t where var is missing
      idx_missing <- which(df[[time_id]] == t & is.na(df[[var]]))
      if (length(idx_missing) == 0) next
      
      # 1. Mean at time t
      values_t <- df[df[[time_id]] == t, var, drop = TRUE]
      if (sum(!is.na(values_t)) >= s) {
        #impute_val <- mean(values_t, na.rm = TRUE)
        impute_val <- rnorm(1,mean(values_t, na.rm = TRUE),0.1)
        
      } else {
        # 2. Neighbors [t-1, t+1]
        neighbor_times <- c()
        if (t > min(time_points)) neighbor_times <- c(neighbor_times, t - 1)
        if (t < max(time_points)) neighbor_times <- c(neighbor_times, t + 1)
        
        values_neighbors <- df[df[[time_id]] %in% neighbor_times, var, drop = TRUE]
        if (sum(!is.na(values_neighbors)) >= s) {
          #impute_val <- mean(values_neighbors, na.rm = TRUE)
          impute_val <- rnorm(1,mean(values_neighbors, na.rm = TRUE),0.1)
          
        } else {
          # 3. Fallback: mean over full data
          #impute_val <- mean(df[[var]], na.rm = TRUE)
          impute_val <- rnorm(1,mean(df[[var]], na.rm = TRUE),0.1)
        }
      }
      
      # Impute missing values at time t
      df_imputed[idx_missing, var] <- impute_val
    }
  }
  
  df_imputed
}

# Introduce missingness
create_missing_mask <- function(data, k, missing_prob = 0.1) {
  variable_matrix <- as.matrix(data[, paste0("X", 1:k)])
  mask <- matrix(runif(nrow(variable_matrix) * k) < missing_prob,
                 nrow = nrow(variable_matrix), ncol = k)
  return(mask)
}

# Apply missingness
apply_missingness <- function(data, missing_mask, k) {
  data_missing <- data
  for (j in 1:k) {
    data_missing[missing_mask[, j], paste0("X", j)] <- NA
  }
  return(data_missing)
}
