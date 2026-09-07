ModelPlot <- function(strengths,
                      covariates = c('1',paste0("X", 1:(length(strengths)-1))),
                      outcome = 'Y',
                      threshold = 0.2) {
  k = length(strengths)
  # Defining covariates names
  covariates_df <- data.frame(name = covariates,
                              x = 0,
                              y = seq(k, 1))
  
  # Defining outcome node
  y_node <- data.frame(name = outcome,
                       x = 2 * k / 3,
                       y = mean(covariates_df$y))
  
  nodes <- rbind(covariates_df, y_node)
  
  # Create edges only where strength != 0
  edges <- data.frame(
    x = covariates_df$x,
    y = covariates_df$y,
    xend = y_node$x,
    yend = y_node$y,
    label = round(strengths, 2),
    strength = strengths
  )
  
  # Filter out zero-strength edges
  edges <- subset(edges, strength != 0)
  
  # Determine line type
  edges$lty <- ifelse(abs(edges$strength) >= threshold, "dotted", "solid")
  
  # Plot
  ggplot() +
    # Curved arrows with labels
    geom_textcurve(
      data = edges,
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend,
        label = label,
        linetype = lty
      ),
      curvature = 0,
      arrow = arrow(length = unit(0.2, "cm")),
      size = 4,
      color = "gray20",
      text_smoothing = 30,
      show.legend = FALSE
    ) +
    # Nodes
    geom_point(
      data = nodes,
      aes(x = x, y = y, color = name),
      size = 16,
      show.legend = FALSE
    ) +
    # Node labels
    geom_text(
      data = nodes,
      aes(x = x, y = y, label = name),
      color = "white",
      size = 5.5
    ) +
    scale_color_manual(values = c(rep("darkgreen", k), "tomato")) +
    theme_void() +
    coord_fixed(xlim = c(-0.5, 0.5 + 2 * k / 3),
                ylim = c(0.5, k + 0.5))
}

InclusionProbsToTransitionNetwork <- function(inclusion_probs_list, threshold = 0.05, default_variance = 1) {
  lapply(names(inclusion_probs_list), function(target_var) {
    probs <- inclusion_probs_list[[target_var]]
    
    # Filter based on threshold
    selected <- probs[abs(probs) > threshold]
    
    # Clean parent names (remove '_t-1'), keep '1' untouched
    cleaned_parents <- gsub("_t-1$", "", names(selected))
    
    # Clean target variable name (remove '_t')
    cleaned_target <- gsub("_t$", "", target_var)
    
    list(
      var_name = cleaned_target,
      parents = cleaned_parents,
      params = as.numeric(selected),
      variance = default_variance
    )
  }) |>
    setNames(gsub("_t$", "", names(inclusion_probs_list)))
}

Summarize_DBNSamples <- function(DBN_results, 
                                 true_values_list = NULL,
                                 varnames = NULL,
                                 burn_in_rate = 0.5,
                                 thin_out = 5,
                                 ncols = 3) {
  
  summaries <- list()
  
  for (j in seq_along(DBN_results)) {
    result <- DBN_results[[j]]
    varname <- if (!is.null(varnames)) varnames[j] else paste0("X_", j)
    
    # Apply MCMC selection (burn-in and thinning)
    selected <- RJ_MC_selection(
      betas = result$betas,
      sigma2s = result$sigma2s,
      delta2s = result$delta2s,
      models = result$models,
      allocation_vectors = result$allocation_vectors,
      burn_in = round(length(result$betas) * burn_in_rate),
      thin_out = thin_out
    )
    # Remove empty model if present in the thinned sample
    beta_samples <- selected$betas
    model_samples <- selected$betas
    # Posterior plot for betas
    cat("Variable:", varname, "\n")
    beta_summary <- PosteriorDistribution(
      posterior_samples = beta_samples,
      true_values = if (!is.null(true_values_list)) true_values_list[[j]] else NULL,
      ncols = ncols,
      covariates = result$covariates
    )
    
    # Posterior plot for sigma² and delta²
    sigma_delta_summary <- PosteriorDistributionSigmadelta(
      mcmc_output = selected,
      true_sigma2 = NULL,
      true_delta2 = NULL
    )
    
    summaries[[varname]] <- list(
      beta_summary = beta_summary,
      sigma_delta_summary = sigma_delta_summary
    )
  }
  
  return(summaries)
}

TransitionNetworkPlot <- function(transition_net,
                                  threshold = 0.05,
                                  node_size = 10,
                                  curvature = 0,
                                  text_size = 4.5) {
  library(ggplot2)
  library(geomtextpath)
  
  vars <- names(transition_net)
  k <- length(vars)
  
  # Layouts for t-1 and t
  layout_tminus1 <- data.frame(
    var = vars,
    label = vars,
    time = "t-1",
    x = 0,
    y = seq(k, 1)
  )
  layout_t <- data.frame(
    var = vars,
    label = vars,
    time = "t",
    x = 3,  # wider spacing between time steps
    y = seq(k, 1)
  )
  nodes <- rbind(layout_tminus1, layout_t)
  
  # Create edges (no labels)
  edges <- do.call(rbind, lapply(transition_net, function(entry) {
    target <- entry$var_name
    parents <- entry$parents
    params <- entry$params
    data.frame(
      from = parents,
      to = target,
      x = layout_tminus1$x[match(parents, layout_tminus1$var)],
      y = layout_tminus1$y[match(parents, layout_tminus1$var)],
      xend = layout_t$x[match(target, layout_t$var)],
      yend = layout_t$y[match(target, layout_t$var)],
      strength = params
    )
  }))
  
  # Filter by threshold
  edges <- subset(edges, abs(strength) > 0)
  edges$lty <- ifelse(abs(edges$strength) >= threshold, "solid", "dotted")
  
  # Plot
  ggplot() +
    # Edges
    geom_curve(
      data = edges,
      aes(x = x, y = y, xend = xend, yend = yend, linetype = lty),
      curvature = curvature,
      arrow = arrow(length = unit(0.2, "cm")),
      color = "gray30",
      show.legend = FALSE
    ) +
    # Time slice borders
    annotate("rect", xmin = -0.6, xmax = 0.6, ymin = 0.5, ymax = k + 0.5,
             fill = NA, color = "black", size = 0.6) +
    annotate("rect", xmin = 2.4, xmax = 3.6, ymin = 0.5, ymax = k + 0.5,
             fill = NA, color = "black", size = 0.6) +
    # Time labels above boxes
    annotate("text", x = 0, y = k + 1.2, label = expression(t - 1), size = 6) +
    annotate("text", x = 3, y = k + 1.2, label = expression(t), size = 6) +
    # Nodes
    geom_point(data = nodes, aes(x = x, y = y, fill = time),
               size = node_size, shape = 21, color = "black") +
    geom_text(data = nodes, aes(x = x, y = y, label = label),
              size = text_size) +
    scale_fill_manual(values = c("t-1" = "skyblue", "t" = "tomato")) +
    coord_fixed(xlim = c(-1, 4), ylim = c(0.5, k + 1.5)) +
    theme_void() +
    ggtitle("Transition Network") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.margin = margin(10, 10, 10, 10)
    )
}

PosteriorDistribution <- function(posterior_samples,
                                  true_values = NULL,
                                  ncols = 3,
                                  covariates = NULL,
                                  inclusion_probs = NULL,
                                  min_inclusion = 0.5,
                                  outcome = "Y") {
  
  nparams <- length(posterior_samples[[1]])
  
  # default covariates
  if (is.null(covariates))
    covariates <- paste0("V", seq_len(nparams))
  covariates[1] <- "Intercept"
  
  # -------------------------------
  # Combine all MCMC chains into one matrix
  # -------------------------------
  mcmc_list <- lapply(
    posterior_samples,
    function(x) matrix(x, ncol = nparams, byrow = TRUE)
  )
  mcmc_matrix <- do.call(rbind, mcmc_list)
  
  # -------------------------------
  # Parent selection based on inclusion probabilities
  # -------------------------------
  if (!is.null(inclusion_probs)) {
    
    if (is.null(names(inclusion_probs)))
      stop("inclusion_probs must be a named vector matching covariates")
    
    probs <- rep(0, nparams)
    names(probs) <- covariates
    
    common <- intersect(names(inclusion_probs), covariates)
    probs[common] <- inclusion_probs[common]
    
    # Always include intercept
    probs["Intercept"] <- 1
    keep_idx <- which(probs > min_inclusion)
    
  } else {
    keep_idx <- seq_len(nparams)
  }
  
  if (length(keep_idx) == 0)
    stop("No parameters selected for plotting.")
  
  covariates_plot <- covariates[keep_idx]
  mcmc_matrix     <- mcmc_matrix[, keep_idx, drop = FALSE]
  
  if (!is.null(true_values))
    true_values <- true_values[keep_idx]
  
  n_plot <- length(keep_idx)
  nrows  <- ceiling(n_plot / ncols)
  
  # -------------------------------
  # Posterior summaries
  # -------------------------------
  posterior_means <- colMeans(mcmc_matrix)
  
  credible_intervals <- apply(
    mcmc_matrix, 2,
    function(x) quantile(x, probs = c(0.025, 0.975))
  )
  colnames(credible_intervals) <- covariates_plot
  
  # -------------------------------
  # Plot layout
  # -------------------------------
  oldpar <- par(no.readonly = TRUE)
  par(
    mfrow = c(nrows, min(ncols, n_plot)),
    oma = c(4, 4, 4, 1),  # reserve outer margins (bottom, left, top, right)
    mar = c(4, 4, 2, 1)   # inner margins for each panel
  )
  
  # Global title (correctly drawn above all panels)
  mtext(
    paste0("Outcome: ", outcome),
    side = 3, outer = TRUE,
    line = 1.5, cex = 1.2, font = 2
  )
  
  # -------------------------------
  # Plot histograms
  # -------------------------------
  for (j in seq_len(n_plot)) {
    
    samples <- mcmc_matrix[, j]
    
    if (!is.null(true_values)) {
      
      x_min <- min(samples, true_values[j], na.rm = TRUE) - 0.05
      x_max <- max(samples, true_values[j], na.rm = TRUE) + 0.05
      
      hist_info <- hist(
        samples,
        main = bquote("Posterior of " ~ beta[.(keep_idx[j] - 1)] ~
                        "(" ~ .(covariates_plot[j]) ~ ")"),
        xlab = bquote(beta[.(keep_idx[j] - 1)]),
        col = "gray",
        border = "white",
        breaks = 20,
        freq = TRUE,
        xlim = c(x_min, x_max)
      )
      
      abline(v = true_values[j], col = "green", lwd = 2)
      y_max <- max(hist_info$counts)
      
    } else {
      
      hist_info <- hist(
        samples,
        main = bquote("Posterior of " ~ beta[.(keep_idx[j] - 1)] ~
                        "(" ~ .(covariates_plot[j]) ~ ")"),
        xlab = bquote(beta[.(keep_idx[j] - 1)]),
        col = "gray",
        border = "white",
        breaks = 20,
        freq = TRUE
      )
      
      y_max <- max(hist_info$counts)
    }
    
    # Posterior mean
    abline(v = posterior_means[j], col = "red", lwd = 2, lty = 2)
    
    # Credible interval
    rect(
      xleft  = credible_intervals[1, j],
      xright = credible_intervals[2, j],
      ybottom = 0,
      ytop    = y_max,
      col = rgb(0, 0, 0.2, alpha = 0.2),
      border = NA
    )
  }
  
  par(oldpar)
  credible_intervals
}




PosteriorDistributionSegmented <- function(posterior_samples,
                                           true_values = NULL,
                                           ncols = 3,
                                           covariates = NULL,
                                           segment_labels = NULL,
                                           colors = NULL,
                                           show_mean = FALSE,
                                           n_changepoints = NULL,
                                           inclusion_probs = NULL,
                                           min_inclusion = 0.5,
                                           outcome = "Y") {
  # Group posterior_samples by number of changepoints (segments - 1)
  segment_counts <- sapply(posterior_samples, length)
  changepoint_counts <- segment_counts - 1
  
  if (!is.null(n_changepoints)) {
    keep <- which(changepoint_counts == n_changepoints)
    if (length(keep) == 0)
      stop("No samples with the specified number of changepoints.")
    posterior_samples <- posterior_samples[keep]
    changepoint_counts <- changepoint_counts[keep]
    n_changepoints_list <- n_changepoints
  } else {
    n_changepoints_list <- sort(unique(changepoint_counts))
  }
  
  all_summaries <- list()
  
  for (k in n_changepoints_list) {
    samples_k <- posterior_samples[changepoint_counts == k]
    
    summary_k <- .plot_segmented_distributions(
      posterior_samples = samples_k,
      true_values = true_values,
      ncols = ncols,
      covariates = covariates,
      segment_labels = paste("Segment", 1:(k + 1)),
      colors = rainbow(k + 1, alpha = 0.7),
      show_mean = show_mean,
      inclusion_probs = inclusion_probs,
      min_inclusion = min_inclusion,
      outcome = outcome,
      n_changepoints = k
    )
    
    all_summaries[[paste0(k, "_changepoints")]] <- summary_k
  }
  
  return(all_summaries)
}



.plot_segmented_distributions <- function(posterior_samples,
                                          true_values,
                                          ncols,
                                          covariates,
                                          segment_labels,
                                          colors,
                                          show_mean,
                                          inclusion_probs = NULL,
                                          min_inclusion = 0.5,
                                          outcome,
                                          n_changepoints) {
  
  n_segments <- length(posterior_samples[[1]])
  n_iters <- length(posterior_samples)
  n_params <- length(posterior_samples[[1]][[1]])
  
  if (is.null(covariates))
    covariates <- paste0("V", 1:n_params)
  covariates[1] <- "Intercept"
  
  # Convert to matrix per segment
  segment_samples <- vector("list", n_segments)
  for (s in 1:n_segments) {
    segment_samples[[s]] <- matrix(NA, nrow = n_iters, ncol = n_params)
  }
  for (iter in seq_along(posterior_samples)) {
    for (s in 1:n_segments) {
      segment_samples[[s]][iter, ] <- as.numeric(posterior_samples[[iter]][[s]])
    }
  }
  
  # ---- Parent selection ----
  if (!is.null(inclusion_probs)) {
    if (is.null(names(inclusion_probs)))
      stop("inclusion_probs must be a named vector matching covariates")
    
    probs_full <- rep(0, n_params)
    names(probs_full) <- covariates
    
    common <- intersect(names(inclusion_probs), covariates)
    probs_full[common] <- inclusion_probs[common]
    probs_full["Intercept"] <- 1
    
    keep_idx <- which(probs_full > min_inclusion)
  } else {
    keep_idx <- 1:n_params
  }
  
  covariates_plot <- covariates[keep_idx]
  if (!is.null(true_values))
    true_values <- true_values[keep_idx]
  
  n_params_plot <- length(keep_idx)
  nrows <- ceiling(n_params_plot / ncols)
  par(mfrow = c(nrows, min(ncols, n_params_plot)))
  
  summary_list <- list()
  
  for (jj in seq_along(keep_idx)) {
    j <- keep_idx[jj]
    
    densities <- lapply(segment_samples, function(mat)
      density(mat[, j], na.rm = TRUE))
    
    x_min <- min(sapply(densities, function(d) min(d$x)))
    x_max <- max(sapply(densities, function(d) max(d$x)))
    y_max <- max(sapply(densities, function(d) max(d$y)))
    
    box_gap <- 0.06 * y_max
    y_box_base <- - (n_segments + 1) * box_gap
    y_lim <- c(y_box_base, 1.1 * y_max)
    
    plot(0, 0,
         type = "n",
         xlim = c(x_min, x_max),
         ylim = y_lim,
         main = paste0(
           "Outcome: ", outcome,
           " | β", j - 1, " (", covariates_plot[jj], ")",
           " | # changepoints = ", n_changepoints
         ),
         xlab = bquote(beta[.(j - 1)]),
         ylab = "Density")
    
    for (s in 1:n_segments)
      lines(densities[[s]], col = colors[s], lwd = 2)
    
    box_height <- 0.025 * y_max
    
    param_summary <- data.frame(
      Segment = segment_labels,
      Mean = NA,
      Median = NA,
      CI_2.5 = NA,
      CI_97.5 = NA
    )
    
    for (s in 1:n_segments) {
      samples_j <- segment_samples[[s]][, j]
      q <- quantile(samples_j, c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm = TRUE)
      stat_value <- if (show_mean) mean(samples_j, na.rm = TRUE) else q[3]
      y_box <- y_box_base + s * box_gap
      
      segments(q[1], y_box, q[5], y_box, col = colors[s], lwd = 2)
      rect(q[2], y_box - box_height / 2,
           q[4], y_box + box_height / 2,
           col = adjustcolor(colors[s], alpha.f = 0.5),
           border = colors[s])
      points(stat_value, y_box, pch = 19, col = colors[s], cex = 1.2)
      text(x_min, y_box, segment_labels[s], pos = 2, cex = 0.8)
      
      param_summary[s, -1] <- c(mean(samples_j, na.rm = TRUE), q[3], q[1], q[5])
    }
    
    if (!is.null(true_values))
      abline(v = true_values[jj], col = "green", lwd = 2)
    
    summary_list[[covariates_plot[jj]]] <- param_summary
  }
  
  par(mfrow = c(1, 1))
  summary_list
}



PosteriorDistributionSigmadelta <- function(mcmc_output, 
                                            true_sigma2 = NULL, 
                                            true_delta2 = NULL,
                                            trim_quantiles = c(0.01, 0.99)) {
  
  sigma2s <- as.numeric(unlist(mcmc_output$sigma2s))
  delta2s <- as.numeric(unlist(mcmc_output$delta2s))
  
  param_samples <- list(sigma2s = sigma2s, delta2s = delta2s)
  param_names <- c(expression(sigma^2), expression(delta^2))
  
  posterior_medians <- sapply(param_samples, median)
  credible_intervals <- lapply(param_samples, function(x)
    quantile(x, probs = c(0.025, 0.975)))
  
  par(mfrow = c(1, 2))
  
  for (j in 1:2) {
    samples <- param_samples[[j]]
    ci <- credible_intervals[[j]]
    true_val <- if (j == 1) true_sigma2 else true_delta2
    
    # Trim extreme outliers for visualization
    bounds <- quantile(samples, probs = trim_quantiles)
    samples_trimmed <- samples[samples >= bounds[1] & samples <= bounds[2]]
    
    # Define plotting range
    x_min <- min(samples_trimmed)* 0.99
    x_max <- max(samples_trimmed)* 1.01
    
    hist_info <- hist(samples_trimmed,
                      main = bquote("Posterior of " ~ .(param_names[[j]])),
                      xlab = as.expression(param_names[[j]]),
                      col = "gray",
                      border = "white",
                      breaks = 20,
                      freq = TRUE,
                      xlim = c(x_min, x_max))
    
    y_max <- max(hist_info$counts)
    
    # Credible interval
    rect(xleft = ci[1],
         xright = ci[2],
         ybottom = 0,
         ytop = y_max,
         col = rgb(0, 0, 0.2, alpha = 0.2),
         border = NA)
    
    # Posterior mean
    abline(v = posterior_medians[[j]], col = "blue", lwd = 2, lty = 2)
    
    # True value (if within trimmed range)
    if (!is.null(true_val) && true_val >= x_min && true_val <= x_max) {
      abline(v = true_val, col = "green", lwd = 2)
    }
  }
  
  par(mfrow = c(1, 1))
  
  # Return credible intervals nicely formatted
  credible_intervals <- do.call(rbind, credible_intervals)
  rownames(credible_intervals) <- c("sigma^2", "delta^2")
  print(credible_intervals)
  credible_intervals
}

CovariateConfigurations <-
  function(covariates, include_empty = FALSE) {
    k <- length(covariates)
    
    # Generate each combination of the covariates in the model
    combos <- lapply(0:k, function(i) {
      combn(covariates, i, simplify = FALSE)
    })
    all_configs <- unlist(combos, recursive = FALSE)
    
    # If include_empty also the model with just the Intercept is considered
    if (!include_empty) {
      all_configs <- Filter(length, all_configs)
    }
    all_configs
  }

DIC <-
  function(mcmc_output,
           data,
           outcome = "Y",
           covariates = colnames(data)[!colnames(data) %in% outcome]) {
    # Extract posteriors
    beta_samples <- mcmc_output$betas
    sigma2_samples <- mcmc_output$sigma2s
    
    # Prepare design matrix
    X <- as.matrix(data[, covariates])
    Y <- data[[outcome]]
    N <- length(Y)
    S <- length(beta_samples)  # number of posterior samples
    
    # Compute deviances for each posterior sample
    deviances <- numeric(S)
    for (s in 1:S) {
      beta <- beta_samples[[s]]
      sigma2 <- sigma2_samples[[s]]
      resid <- Y - X %*% beta
      deviances[s] <- N * log(2 * pi * sigma2) + sum(resid ^ 2) / sigma2
    }
    
    # Posterior means
    beta_mean <- Reduce("+", beta_samples) / S
    sigma2_mean <- mean(unlist(sigma2_samples))
    resid_mean <- Y - X %*% beta_mean
    dev_at_mean <-
      N * log(2 * pi * sigma2_mean) + sum(resid_mean ^ 2) / sigma2_mean
    
    # DIC Computation
    DIC <- 2 * mean(deviances) - dev_at_mean
    cat(DIC)
    cat('\n\n')
    DIC
  }

TracePlots = function(sigma2_samples, delta2_samples, 
                      true_sigma2 = NULL, true_delta2 = NULL) {
  # Convert lists to numeric vectors
  sigma2_vec <- unlist(sigma2_samples)
  delta2_vec <- unlist(delta2_samples)
  
  par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))  # Two plots stacked
  
  # Trace plot for log(sigma²)
  plot(sigma2_vec, type = "l", col = "blue",
       ylab = expression(log(sigma^2)), xlab = "Iteration",
       main = expression("Traceplot of log(" ~ sigma^2 ~ ")"),
       log = "y")
  if (!is.null(true_sigma2)) {
    abline(h = true_sigma2, col = "darkgreen", lwd = 2, lty = 2)
  }
  
  # Trace plot for log(delta²)
  plot(delta2_vec, type = "l", col = "purple",
       ylab = expression(log(delta^2)), xlab = "Iteration",
       main = expression("Traceplot of log(" ~ delta^2 ~ ")"),
       log = "y")
  if (!is.null(true_delta2)) {
    abline(h = true_delta2, col = "darkgreen", lwd = 2, lty = 2)
  }
  
  par(mfrow = c(1, 1))  # Reset layout
}

IsListOfVectors <- function(x, length) {
  is.list(x) &&
    all(sapply(x, is.vector)) &&
    length(x) == length
}


AverageM <- function(transition_network){
  all_nodes <- names(transition_network)
  # Initialize binary matrix
  parent_matrix <- matrix(0, nrow = length(all_nodes), ncol = length(all_nodes),
                          dimnames = list(all_nodes, all_nodes))
  # Fill matrix: row = target, col = parent
  for (target in names(transition_network)) {
    parents <- setdiff(transition_network[[target]]$parents, "1")
    parent_matrix[target, parents] <- 1
  }
  parent_matrix
}

PrecisionRecall <- function(real_transition_network, inclusion_prob_matrix){
  # Flatten both matrices to vectors for easy comparison
  true_vec <- as.vector(real_transition_network)
  pred_vec <- as.vector(inclusion_prob_matrix)
  
  # Compute confusion matrix components
  TP <- sum(true_vec == 1 & pred_vec == 1)  # True Positives
  FP <- sum(true_vec == 0 & pred_vec == 1)  # False Positives
  FN <- sum(true_vec == 1 & pred_vec == 0)  # False Negatives
  
  # Compute precision and recall
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  
  cat("Precision:", ifelse(is.nan(precision), NA, precision), "\n")
  cat("Recall:", ifelse(is.nan(recall), NA, recall), "\n")
  
  list(precision, recall)
}

InclusionProbMatrix <- function(inclusion_probs) {
  # Get all unique covariates, excluding "1"
  all_covariates <- sort(setdiff(unique(unlist(lapply(inclusion_probs, names))), "1"))
  all_targets <- names(inclusion_probs)
  
  # Initialize matrix with 0s
  prob_matrix <- matrix(0, nrow = length(all_targets), ncol = length(all_covariates),
                        dimnames = list(all_targets, all_covariates))
  
  # Fill in the inclusion probabilities
  for (target in all_targets) {
    probs <- inclusion_probs[[target]]
    probs <- probs[setdiff(names(probs), "1")]  # remove "1" if present
    covars <- names(probs)
    prob_matrix[target, covars] <- unlist(probs)
  }
  
  return(prob_matrix)
}

list_to_array <- function(model_list) {
  P <- nrow(model_list[[1]])
  Q <- ncol(model_list[[1]])
  T <- length(model_list)
  arr <- array(NA, dim = c(P, Q, T),
               dimnames = list(rownames(model_list[[1]]),
                               colnames(model_list[[1]]),
                               NULL))
  for (t in 1:T) {
    arr[,,t] <- model_list[[t]]
  }
  return(arr)
}


ComputePSRF <- function(adj_array, threshold = 1.1, every = 100, stop_after = NULL) {
  P <- dim(adj_array)[1]
  Q <- dim(adj_array)[2]
  T <- dim(adj_array)[3]
  M <- dim(adj_array)[4]
  
  # Set maximum iteration for evaluation
  max_iter <- if (!is.null(stop_after)) min(T, stop_after) else T
  selected_iters <- seq(every, max_iter, by = every)
  rates <- numeric(length(selected_iters))
  
  # Precompute observed edges (at least one nonzero over all chains & time)
  observed_edges <- apply(adj_array, c(1,2), function(x) any(!is.na(x) & x != 0))
  
  for (idx in seq_along(selected_iters)) {
    t <- selected_iters[idx]
    converged_count <- 0
    observed_count <- 0
    
    for (i in 1:P) {
      for (j in 1:Q) {
        if (!observed_edges[i, j]) next
        
        # Samples up to iteration t, across all chains
        samples <- matrix(NA, nrow = t, ncol = M)
        for (m in 1:M) {
          samples[, m] <- adj_array[i, j, 1:t, m]
        }
        
        chain_means <- colMeans(samples, na.rm = TRUE)
        within_vars <- apply(samples, 2, function(x) if (all(is.na(x))) NA_real_ else var(x, na.rm = TRUE))
        W <- mean(within_vars, na.rm = TRUE)
        B <- if (M > 1) t * var(chain_means, na.rm = TRUE) else 0
        
        if (is.na(W) || W == 0) {
          R_hat <- 1
        } else {
          V_hat <- ((t - 1) / t) * W + (1 / t) * B
          R_hat <- sqrt(V_hat / W)
        }
        
        observed_count <- observed_count + 1
        if (R_hat < threshold) {
          converged_count <- converged_count + 1
        }
      }
    }
    
    rates[idx] <- if (observed_count > 0) converged_count / observed_count else NA
  }
  
  return(data.frame(iteration = selected_iters, rate = rates))
}

ComputePSRF_incremental <- function(adj_array,
                                    threshold = 1.1,
                                    every = 100,
                                    stop_after = NULL,
                                    min_iter = 20,
                                    verbose = FALSE) {
  
  dims <- dim(adj_array)
  P <- dims[1]
  Q <- dims[2]
  T <- dims[3]
  M <- dims[4]
  
  max_iter <- if (!is.null(stop_after)) min(T, stop_after) else T
  
  selected_iters <- seq(every, max_iter, by = every)
  selected_iters <- selected_iters[selected_iters >= min_iter]
  
  E <- P * Q
  arr <- array(adj_array, dim = c(E, T, M))
  
  observed_edges <- apply(arr, 1, function(x) any(!is.na(x) & x != 0))
  arr <- arr[observed_edges, , , drop = FALSE]
  
  n_edges <- dim(arr)[1]
  
  # -------------------------
  # running stats
  # -------------------------
  sum_x  <- array(0, c(n_edges, M))
  sum_x2 <- array(0, c(n_edges, M))
  count  <- matrix(0, n_edges, M)
  
  rates <- numeric(length(selected_iters))
  eval_idx <- 1
  
  for (t in 1:max_iter) {
    
    xt <- matrix(arr[, t, ], nrow = n_edges, ncol = M)
    
    valid <- !is.na(xt)
    
    sum_x[valid]  <- sum_x[valid] + xt[valid]
    sum_x2[valid] <- sum_x2[valid] + xt[valid]^2
    count[valid]  <- count[valid] + 1
    
    if (t %in% selected_iters) {
      
      # per-chain means
      chain_mean <- sum_x / pmax(count, 1)
      
      # within-chain variance
      chain_var <- (sum_x2 / pmax(count, 1)) - chain_mean^2
      W <- rowMeans(chain_var, na.rm = TRUE)
      
      # between-chain variance
      n <- rowMeans(count)
      
      if (M > 1) {
        B <- n * apply(chain_mean, 1, var, na.rm = TRUE)
      } else {
        B <- rep(0, n_edges)
      }
      
      R_hat <- rep(NA_real_, n_edges)
      
      valid_edges <- is.finite(W) & W > 0 & is.finite(B)
      
      if (any(valid_edges)) {
        
        V_hat <- ((n[valid_edges] - 1) / n[valid_edges]) * W[valid_edges] +
          (B[valid_edges] / n[valid_edges]) +
          (B[valid_edges] / (M * n[valid_edges]))
        
        R_hat[valid_edges] <- sqrt(V_hat / W[valid_edges])
      }
      
      # enforce theoretical lower bound
      R_hat[R_hat < 1] <- 1
      
      # CRITICAL: correct threshold usage
      rates[eval_idx] <- mean(R_hat < threshold, na.rm = TRUE)
      
      if (verbose) {
        cat(
          "Iter:", t,
          "| threshold:", threshold,
          "| mean Rhat:", mean(R_hat, na.rm = TRUE),
          "| max Rhat:", max(R_hat, na.rm = TRUE),
          "\n"
        )
      }
      
      eval_idx <- eval_idx + 1
    }
  }
  
  data.frame(
    iteration = selected_iters,
    rate = rates
  )
}

EvaluateMissingImputation <- function(data = NULL,
                                      df_with_missing,
                                      missing_posterior_samples,
                                      sample_id = "Sample_Id",
                                      time_id = "Time",
                                      burn_in = 500,
                                      credible_level = 0.95,
                                      variable_order = NULL) {
  
  library(dplyr)
  library(ggplot2)
  
  # Order rows
  df_with_missing <- df_with_missing %>%
    arrange(.data[[sample_id]], .data[[time_id]])
  
  if (!is.null(data)) {
    data <- data %>%
      arrange(.data[[sample_id]], .data[[time_id]])
  }
  
  # Variable names
  df_with_missing <- df_with_missing[, sort(names(df_with_missing))]
  var_names <- colnames(df_with_missing)[!colnames(df_with_missing) %in% c(sample_id, time_id)]
  
  # Missing mask
  missing_mask <- df_with_missing %>%
    mutate(across(-c(sample_id, time_id), ~ is.na(.)))
  
  miss_idx <- which(
    as.matrix(missing_mask[var_names]),
    arr.ind = TRUE
  )
  
  node_per_missing <- miss_idx[, "col"]
  
  # Posterior samples
  posterior_matrix <- do.call(rbind, missing_posterior_samples)
  posterior_matrix <- posterior_matrix[(burn_in + 1):nrow(posterior_matrix), , drop = FALSE]
  
  n_missings <- ncol(posterior_matrix)
  
  # Credible intervals
  alpha <- (1 - credible_level) / 2
  credible_intervals <- t(
    apply(posterior_matrix, 2, quantile, probs = c(alpha, 1 - alpha))
  )
  
  # Variable names for missing entries
  variable_names_raw <- var_names[node_per_missing]
  
  # Apply custom ordering IF provided
  if (!is.null(variable_order)) {
    
    if (!all(variable_order %in% var_names)) {
      stop("All variable_order entries must be valid variable names.")
    }
    
    remaining_vars <- setdiff(var_names, variable_order)
    final_order <- c(variable_order, remaining_vars)
    
    variable_names <- factor(variable_names_raw, levels = final_order)
    
  } else {
    variable_names <- factor(variable_names_raw, levels = var_names)
  }
  
  # Local index
  local_indices <- ave(seq_along(variable_names), variable_names, FUN = seq_along)
  
  distributions_df <- data.frame(
    Variable = variable_names,
    LocalIndex = local_indices,
    Lower = credible_intervals[, 1],
    Upper = credible_intervals[, 2]
  )
  
  # Truth-based evaluation
  if (!is.null(data)) {
    
    true_values <- as.matrix(data[var_names])[miss_idx]
    
    quantiles <- sapply(seq_len(n_missings), function(i) {
      ecdf(posterior_matrix[, i])(true_values[i])
    })
    
    distributions_df$true_values <- true_values
    coverage <- mean(quantiles >= alpha & quantiles <= (1 - alpha))
    
  } else {
    quantiles <- NULL
    coverage <- NULL
  }
  
  # Plot
  p <- ggplot(
    distributions_df,
    aes(x = Variable, group = interaction(Variable, LocalIndex))
  ) +
    
    geom_linerange(
      aes(ymin = Lower, ymax = Upper, linetype = "Credible Interval"),
      position = position_dodge(width = 0.6),
      color = "gray40"
    ) +
    
    labs(
      title = paste0(
        "Posterior ",
        round(100 * credible_level),
        "% Credible Intervals for Missing Values"
      ),
      x = "Variable",
      y = "Value",
      linetype = NULL,
      shape = NULL
    ) +
    
    scale_linetype_manual(
      values = "solid",
      labels = paste0(round(100 * credible_level), "% Credible Interval")
    ) +
    
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(size = 22, hjust = 0.5),
      
      axis.title.x = element_text(size = 18, margin = margin(t = 12)),
      axis.title.y = element_text(size = 18, margin = margin(r = 12)),
      
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      
      legend.position = "bottom"
    ) +
    
    guides(
      linetype = guide_legend(order = 1),
      shape = guide_legend(order = 2),
      color = "none"
    )
  
  # Add true values (smaller black dots)
  if (!is.null(data)) {
    p <- p +
      geom_point(
        aes(y = true_values, shape = "True Value"),
        position = position_dodge(width = 0.6),
        size = 2,
        color = "black"
      ) +
      scale_shape_manual(
        values = 16,
        labels = "True Value"
      )
  }
  
  list(
    plot = p,
    distributions_df = distributions_df,
    coverage = coverage,
    quantiles = quantiles
  )
}