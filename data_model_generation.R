CompleteLinearRegressionModelGeneration = function(n_strong_predictors = 3,
                                                   n_week_predictors = 0,
                                                   n_non_predictors = 0,
                                                   a_1 = 3,
                                                   a_2 = 5,
                                                   b_1 = 0.2,
                                                   b_2 = 0.5,
                                                   M = rep(0, n_strong_predictors + n_week_predictors + n_non_predictors),
                                                   V = rep(1, n_strong_predictors + n_week_predictors + n_non_predictors)) {
  # Randomly generating Betas for strong predictors from the Uniform distribution in the space [-a_2,-a_1] U [a_1,a_2]
  betas <-
    ifelse(rbinom(n_strong_predictors, 1, 0.5) == 1, 1,-1) * runif(n_strong_predictors, min = a_1, max = a_2)
  
  # Randomly generating Betas for week predictors from the Uniform distribution in the space [-b_2,-b_1] U [b_1,b_2]
  betas <- c(betas, ifelse(rbinom(n_week_predictors, 1, 0.5) == 1, 1,-1) * runif(n_week_predictors, min = b_1, max = b_2))
  
  # Setting betas for non predictors equal to 0  
  betas <- c(betas, rep(0, n_non_predictors))
  
  # Shuffle the parameters
  betas <- sample(betas)
  
  # Saving all the parameters of the model in a list
  list(
    # Intercept + Linear Coefficients
    betas = betas,
    # Means
    means = M,
    # Standard Deviations
    variances = V
  )
}

DataGeneration = function(N = 100, parameters) {
  # Number of Variables
  p <- length(parameters$means)
  if (length(parameters$variances) != p) {
    stop("Means and variances must have the same length.")
  }
  
  # Generating Covariates Values given their means and variances
  data <-
    mapply(function(mean, variance)
      rnorm(N, mean = mean, sd = sqrt(variance)),
      parameters$means[1:(p - 1)],
      parameters$variances[1:(p - 1)])
  
  df <- as.data.frame(cbind(1,data))
  colnames(df) <- c(1,paste0("X", 1:(p - 1)))
  
  # Generating Outcome Values by combining normals: linear combination of the covariates plus a normal error
  Y <-
    as.matrix(df) %*% parameters$betas + rnorm(N, parameters$means[p], sqrt(parameters$variances[p]))
  df$Y <- Y
  return(df)
}

TransitionNetwork <- function(k,
                              max_parent_set = 3,
                              a1 = 0.1, 
                              a2 = 1,
                              V = 1,
                              var_set = sapply(1:k, function(i) {paste0("X_", i)})){
  topological_order <- sample(var_set)
  transition_network <- list()
  for (i in 1:k){
    if (i == 1){
      transition_network[[topological_order[i]]] <- list(var_name = topological_order[i],
                                                         parents = c('1'),
                                                         params = sample(c(-1, 1), 1, replace = TRUE) * runif(1, a1, a2),
                                                         variance = V)
    }
    else{
      n_parents <- sample(0:min(i-1,max_parent_set),1)
      parent_set <- topological_order[sample(1:(i-1), n_parents)]
      transition_network[[topological_order[i]]] <- list(var_name = topological_order[i],
                                                         parents = c('1',parent_set),
                                                         params = sample(c(-1, 1), n_parents + 1, replace = TRUE) * runif(n_parents + 1, a1, a2),
                                                         variance = V)
    }
  }
  transition_network[order(names(transition_network))]
}

Network0 <- function(k,
                     a1 = 0.1,
                     a2 = 1,
                     V = 1,
                     var_set = sapply(1:k, function(i) {paste0("X_", i)})){
  setNames(
    lapply(var_set, function(i) list(
      variable = i,
      params = sample(c(-1, 1), 1) * runif(1, a1, a2),
      variance = V
    )),
    var_set
  )
}

DBNGeneration <- function(k,
                          max_parent_set = 3,
                          a1 = 0.1, 
                          a2 = 1,
                          V = 1,
                          var_set = sapply(1:k, function(i) {paste0("X_", i)})){
  list(Network0 = Network0(k = k, a1 = a1, a2 = a2, V = V, var_set = var_set),
       TransitionNetwork = TransitionNetwork(k = k, max_parent_set = max_parent_set, a1 = a1, a2 = a2, V = V, var_set = var_set))
}

DataGenerationFromDBN <- function(N, T, dbn, scaling = FALSE) {
  initial_net <- dbn$Network0
  transition_net <- dbn$TransitionNetwork
  
  k <- length(initial_net)
  variable_names <- names(initial_net)
  
  # --- Time 0: Generate from Network0 ---
  initial_data <- replicate(N, {
    sapply(initial_net, function(model) {
      rnorm(1, mean = model$params, sd = sqrt(model$variance))
    })
  })
  initial_df <- data.frame(Sample_Id = 1:N, Time = 0, t(initial_data))
  colnames(initial_df)[3:ncol(initial_df)] <- variable_names
  
  df <- initial_df
  
  # --- Time 1 to T: Generate from TransitionNetwork ---
  for (t in 1:T) {
    past_df <- df[df$Time == (t - 1), ]
    new_df <- data.frame(Sample_Id = 1:N, Time = t)
    
    for (var in variable_names) {
      model <- transition_net[[var]]
      parents <- model$parents
      params <- model$params
      
      intercept <- params[1]
      parent_vars <- parents[-1]
      parent_params <- if (length(parents) > 1) params[-1] else numeric(0)
      
      if (length(parent_vars) > 0) {
        predictors <- past_df[, parent_vars, drop = FALSE]
        if (scaling) {
          predictors <- scale(predictors)
        }
        lin_comb <- intercept + as.numeric(as.matrix(predictors) %*% parent_params)
      } else {
        lin_comb <- rep(intercept, N)
      }
      
      new_df[[var]] <- rnorm(N, mean = lin_comb, sd = sqrt(model$variance))
    }
    
    df <- rbind(df, new_df)
  }
  
  return(df)
}

# --- Helper: sample a perturbation from a spherical shell
sample_shell_signfree <- function(p, c_min, c_max) {
  u <- rnorm(p)
  v <- u / sqrt(sum(u^2))
  r <- runif(1, c_min, c_max)
  return(r * v)  # all positive magnitude, sign applied later
}


NH_DBNetwork <- function(k,
                         T,
                         changepoints = NULL,
                         max_parent_set = 3,
                         a1 = 0.1, 
                         a2 = 1,
                         V = 1,
                         var_set = sapply(1:k, function(i) paste0("X_", i)),
                         param_type = c("independent",
                                        "same_sign",
                                        "flip_sign",
                                        "similar",
                                        "global_coupling"),
                         diff_param = 0.1,
                         global_var = 1,
                         epsilon = 0.1) {
  
  param_type <- match.arg(param_type)
  
  # --------------------------------------------------
  # Handle changepoints
  # --------------------------------------------------
  if (is.null(changepoints)) {
    n_segments <- 3
    changepoints <- round(seq(1, T, length.out = n_segments + 1))[-1]
  }
  changepoints <- sort(unique(changepoints))
  segment_starts <- c(1, changepoints)
  segment_ends   <- c(changepoints - 1, T)
  n_segments     <- length(segment_starts)
  
  # --------------------------------------------------
  # Fixed DAG structure
  # --------------------------------------------------
  topological_order <- sample(var_set)
  structure <- list()
  
  for (i in 1:k) {
    if (i == 1) {
      structure[[topological_order[i]]] <- c("1")
    } else {
      n_parents <- sample(0:min(i - 1, max_parent_set), 1)
      structure[[topological_order[i]]] <-
        c("1", topological_order[sample(1:(i - 1), n_parents)])
    }
  }
  
  # --------------------------------------------------
  # Storage
  # --------------------------------------------------
  parameter_sets <- vector("list", n_segments)
  
  # For similar
  if (param_type == "similar") {
    global_means <- list()
    start_low    <- list()
  }
  
  # For global coupling
  if (param_type == "global_coupling") {
    global_mu <- list()
    raw_params <- vector("list", n_segments)
    for (s in 1:n_segments) raw_params[[s]] <- list()
  }
  
  # --------------------------------------------------
  # Generate parameters per segment
  # --------------------------------------------------
  for (s in 1:n_segments) {
    
    parameter_sets[[s]] <- list()
    
    for (var in topological_order) {
      
      parents    <- structure[[var]]
      n_parents  <- length(parents) - 1  # exclude intercept
      
      # -------------------------------
      # SIMILAR
      # -------------------------------
      if (param_type == "similar") {
        
        if (s == 1) {
          mag  <- runif(n_parents + 1, a1, a2)
          sign <- sample(c(-1, 1), n_parents + 1, replace = TRUE)
          global_means[[var]] <- mag * sign
          start_low[[var]]    <- rbinom(n_parents + 1, 1, 0.5)
        }
        
        params <- numeric(n_parents + 1)
        for (j in 1:(n_parents + 1)) {
          if (s == 1) {
            offset <- if (start_low[[var]][j] == 1)
              -diff_param / 2 else diff_param / 2
          } else {
            flip <- (s %% 2 == 0)
            offset <- if (xor(start_low[[var]][j] == 1, flip))
              -diff_param / 2 else diff_param / 2
          }
          params[j] <- global_means[[var]][j] + offset
        }
        
        # -------------------------------
        # GLOBAL COUPLING
        # -------------------------------
      } else if (param_type == "global_coupling") {
        
        if (s == 1) {
          global_mu[[var]] <-
            rnorm(n_parents + 1, mean = 0, sd = sqrt(global_var))
        }
        
        raw_params[[s]][[var]] <-
          rnorm(n_parents + 1,
                mean = global_mu[[var]],
                sd   = sqrt(epsilon))
        
        # temporary placeholder (overwritten after normalization)
        params <- raw_params[[s]][[var]]
        
        # -------------------------------
        # INDEPENDENT
        # -------------------------------
      } else if (s == 1 || param_type == "independent") {
        
        params <- sample(c(-1, 1), n_parents + 1, replace = TRUE) *
          runif(n_parents + 1, a1, a2)
        
        # -------------------------------
        # SAME SIGN
        # -------------------------------
      } else if (param_type == "same_sign") {
        
        prev_params <- parameter_sets[[s - 1]][[var]]$params
        params <- sign(prev_params) *
          runif(n_parents + 1, a1, a2)
        
        # -------------------------------
        # FLIP SIGN
        # -------------------------------
      } else if (param_type == "flip_sign") {
        
        prev_params <- parameter_sets[[s - 1]][[var]]$params
        params <- -prev_params
      }
      
      parameter_sets[[s]][[var]] <- list(
        var_name = var,
        parents  = parents,
        params   = params,
        variance = V
      )
    }
  }
  
  # --------------------------------------------------
  # Normalize global coupling across segments
  # --------------------------------------------------
  if (param_type == "global_coupling") {
    
    for (var in topological_order) {
      p <- length(raw_params[[1]][[var]])
      
      for (j in 1:p) {
        vec <- sapply(1:n_segments,
                      function(s) raw_params[[s]][[var]][j])
        vec <- vec / sqrt(sum(vec^2))
        
        for (s in 1:n_segments) {
          parameter_sets[[s]][[var]]$params[j] <- vec[s]
        }
      }
    }
  }
  
  # --------------------------------------------------
  # Output
  # --------------------------------------------------
  list(
    Network0 = Network0(k = k,
                        a1 = a1,
                        a2 = a2,
                        V  = V,
                        var_set = var_set),
    changepoints      = changepoints,
    segment_starts    = segment_starts,
    segment_ends      = segment_ends,
    TransitionNetworks = parameter_sets
  )
}




DataGenerationFromNHDBN <- function(N, T, nh_dbn, scaling = FALSE) {
  initial_net <- nh_dbn$Network0
  variable_names <- names(initial_net)
  segment_starts <- nh_dbn$segment_starts
  segments <- nh_dbn$TransitionNetworks
  n_segments <- length(segments)
  
  # Time 0
  initial_data <- replicate(N, {
    sapply(initial_net, function(model) {
      rnorm(1, mean = model$params, sd = sqrt(model$variance))
    })
  })
  df <- data.frame(Sample_Id = 1:N, Time = 0, t(initial_data))
  colnames(df)[3:ncol(df)] <- variable_names
  
  # Time 1 to T
  for (t in 1:T) {
    past_df <- df[df$Time == (t - 1), ]
    new_df <- data.frame(Sample_Id = 1:N, Time = t)
    
    # Identify segment
    s_idx <- max(which(segment_starts <= t))
    current_params <- segments[[s_idx]]
    
    for (var in variable_names) {
      model <- current_params[[var]]
      parents <- model$parents
      params <- model$params
      intercept <- params[1]
      parent_vars <- parents[-1]
      parent_params <- if (length(parents) > 1) params[-1] else numeric(0)
      
      if (length(parent_vars) > 0) {
        predictors <- past_df[, parent_vars, drop = FALSE]
        if (scaling) {
          predictors <- scale(predictors)
        }
        lin_comb <- intercept + as.numeric(as.matrix(predictors) %*% parent_params)
      } else {
        lin_comb <- rep(intercept, N)
      }
      
      new_df[[var]] <- rnorm(N, mean = lin_comb, sd = sqrt(model$variance))
    }
    
    df <- rbind(df, new_df)
  }
  
  return(df)
}
